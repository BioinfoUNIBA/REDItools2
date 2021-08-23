#!/usr/bin/env python

import os
import glob
import sys
import re
import time
from mpi4py import MPI
import datetime
from collections import OrderedDict
import reditools
import argparse
import gc
import socket
import netifaces
from random import shuffle

ALIGN_CHUNK = 0
STOP_WORKING = 1
IM_FREE = 2
CALCULATE_COVERAGE = 3

STEP = 10000000

def weight_function(x):
    # x = math.log(1+x)
    # return 2.748*10**(-3)*x**3 -0.056*x**2 + 0.376*x + 2.093
    return x**3

def get_coverage(coverage_file, region = None):

    # Open the file and read i-th section (jump to the next '\n' character)
    n = float(sample_size)
    file_size = os.path.getsize(coverage_file)
    print("[{}] SIZE OF FILE {}: {} bytes".format(sample_rank, coverage_file, file_size))
    start = int(sample_rank*(file_size/n))
    end = int((sample_rank+1)*(file_size/n))
    print("[{}] [DEBUG] START={} END={}".format(sample_rank, start, end))

    f = open(coverage_file, "r")
    f.seek(start)
    loaded = start
    coverage_partial = 0
    with f as lines:
        line_no = 0
        for line in lines:
            if loaded >= end: continue
            loaded += len(line)

            line_no += 1
            if line_no == 1:
                if not line.startswith("chr"):
                    continue
            
            triple = line.rstrip().split("\t")
            
            if region is not None:
                if triple[0] != region[0]: continue
                if len(region) >= 2 and int(triple[1]) < region[1]: continue
                if len(region) >= 2 and int(triple[1]) > region[2]: continue
                
            #if line_no % 10000000 == 0:
            #    print("[{}] [DEBUG] Read {} lines so far".format(rank, line_no))
            cov = int(triple[2])
            coverage_partial += weight_function(cov)

    print("[{}] START={} END={} PARTIAL_COVERAGE={}".format(sample_rank, start, end, coverage_partial))    

    # Reduce
    coverage = None

    coverages = sample_comm.gather(coverage_partial)
    if sample_rank == 0:
        print("COVERAGES:", str(coverages))
        coverage = reduce(lambda x,y: x+y, coverages)
        
    coverage = sample_comm.bcast(coverage)
        
    # Return the total
    return coverage

def calculate_intervals(total_coverage, coverage_file, region):
    print("[SYSTEM] [{}] Opening coverage file={}".format(sample_rank, coverage_file))
    f = open(coverage_file, "r")

    chr = None
    start = None
    end = None
    C = 0
    max_interval_width = min(3000000, 3000000000 / sample_size)

    subintervals = []
    subtotal = total_coverage / sample_size
    print("[SYSTEM] TOTAL={} SUBTOTAL={} MAX_INTERVAL_WIDTH={}".format(total_coverage, subtotal, max_interval_width))

    line_no = 0
    with f as lines:
        for line in lines:
            line_no += 1
            if line_no % 1000000 == 0:
                print("[SYSTEM] [{}] Time: {} - {} lines loaded.".format(sample_rank, time.time(), line_no))

            fields = line.rstrip().split("\t")
            
            if region is not None:
                if fields[0] != region[0]: continue
                if len(region) >= 2 and int(fields[1]) < region[1]: continue
                if len(region) >= 3 and int(fields[1]) > region[2]: continue

            # If the interval has become either i) too large or ii) too heavy or iii) spans across two different chromosomes
            if C >= subtotal or (chr is not None and fields[0] != chr) or (end is not None and start is not None and (end-start) > max_interval_width):
                reason = None
                if C >= subtotal: reason = "WEIGHT"
                elif chr is not None and fields[0] != chr: reason = "END_OF_CHROMOSOME"
                elif end is not None and start is not None and (end-start) > max_interval_width: reason = "MAX_WIDTH"

                interval = (chr, start, end, C, end-start, reason)
                print("[SYSTEM] [{}] Time: {} - Discovered new interval={}".format(sample_rank, time.time(), interval))
                subintervals.append(interval)
                chr = None
                start = None
                end = None
                C = 0
            if len(fields) < 3: continue

            if chr is None: chr = fields[0]
            if start is None: start = int(fields[1])
            end = int(fields[1])
            C += weight_function(int(fields[2]))

        if C > 0:
            reason = "END_OF_CHROMOSOME"
            interval = (chr, start, end, C, end-start, reason)
            print("[SYSTEM] [{}] Time: {} - Discovered new interval={}".format(sample_rank, time.time(), interval))
            subintervals.append(interval)

    return subintervals

if __name__ == '__main__':

    # MPI init
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    options = reditools.parse_options()
    # Until tabix can properly index the output, the header is disruptive.
    options["remove_header"] = True
    
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('-D', '--coverage-dir', help='The coverage directory containing the coverage file of the sample to analyze divided by chromosome')
    parser.add_argument('-t', '--temp-dir', help='The temp directory where to store temporary data for this sample')
    parser.add_argument('-Z', '--chromosome-sizes', help='The file with the chromosome sizes')
    parser.add_argument('-g', '--region', help='The region of the bam file to be analyzed')
    parser.add_argument('-F', '--samples-file', help='The file listing each bam file to be analyzed on a separate line')
    args = parser.parse_known_args()[0]
    
    coverage_dir = args.coverage_dir
    temp_dir = args.temp_dir
    size_file = args.chromosome_sizes
    samples_filepath = args.samples_file
    
    samples = None
    if rank == 0:
        samples = []
        for line in open(samples_filepath, "r"):
            line = line.strip()
            samples.append(line)
        
    # Chronometer data structure
    if rank == 0:   
        chronometer = {}
        for sample in samples:
            sample = os.path.basename(sample)
            sample = ".".join(sample.split(".")[0:-1])
            chronometer[sample] = {
#                 "coverage": 0,
#                 "intervals": 0,
                "parallel": 0
            }
        
        print("CHRONOMETER", chronometer)
    
    samples = comm.bcast(samples, root=0)
    
    PROCS_PER_SAMPLE = size / len(samples)
    if rank == 0:
        print("[{}] PROCESSES_PER_SAMPLE={}".format(rank, PROCS_PER_SAMPLE))
    
    interface = 'ib0' if 'ib0' in netifaces.interfaces() else netifaces.interfaces()[0]
    hostname = socket.gethostbyaddr(netifaces.ifaddresses(interface)[netifaces.AF_INET][0]['addr'])
    pid = os.getpid()
#     print("[SYSTEM] [TECH] [NODE] RANK:{} HOSTNAME:{} PID:{}".format(rank, hostname, pid))
    
#     if rank == 0:
#         print("[SYSTEM] LAUNCHED PARALLEL REDITOOLS WITH THE FOLLOWING OPTIONS:", options, args)
        
    region = None
    if args.region:
        region = re.split("[:-]", args.region)
        if not region or len(region) == 2 or (len(region) == 3 and region[1] == region[2]):
            sys.stderr.write("[ERROR] Please provide a region of the form chrom:start-end (with end > start). Region provided: {}".format(region))
            exit(1)
        if len(region) >= 2:
            region[1] = int(region[1])
            region[2] = int(region[2])
    
    t1 = time.time()

#     print("I am rank #"+str(rank))
    
    # COVERAGE SECTION
    sample_index = rank/PROCS_PER_SAMPLE
    sample_filepath = samples[sample_index]
    sample = os.path.basename(sample_filepath)
    sample = ".".join(sample.split(".")[0:-1])
    
    if rank % PROCS_PER_SAMPLE == 0:
        print("[{}] SAMPLE_INDEX={} SAMPLE_FILEPATH={} SAMPLE={}".format(rank, sample_index, sample_filepath, sample))
    
    sample_comm = comm.Split(sample_index)
    sample_rank = sample_comm.Get_rank()
    sample_size = sample_comm.Get_size()
    
    coverage_dir += sample + "/"
    coverage_file = coverage_dir + sample + ".cov"
    temp_dir += sample + "/"
    
    if not os.path.isfile(coverage_file):
        print("[ERROR] Coverage file {} not existing!".format(coverage_file))
        exit(1)
    
    try:
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
    except Exception as e:
        print("[WARN] {}".format(e))

    interval_file = temp_dir + "/intervals.txt"
    homeworks = []
    
    if os.path.isfile(interval_file) and os.stat(interval_file).st_size > 0:
        if sample_rank == 0:
            print("[{}] [{}] [S{}] [RESTART] FOUND INTERVAL FILE {} ".format(rank, sample_rank, sample_index, interval_file))
            expected_total = 0
            for line in open(interval_file, "r"):
                line = line.strip()
                
                if expected_total == 0:
                    expected_total = int(line)
                    continue
                
                # Interval format: (chr, start, end, C, end-start, reason)
                fields = line.split("\t")
                for i in range(1, 5):
                    fields[i] = int(fields[i])
                homeworks.append([sample_index] + fields)
            print("[{}] [{}] [S{}] [RESTART] INTERVAL FILE #INTERVALS {} ".format(rank, sample_rank, sample_index, len(homeworks)))
            
    else:
        if sample_rank == 0:
            print("["+str(rank)+"] [S"+str(sample_index)+"] PRE-COVERAGE TIME " + str(datetime.datetime.now().time()))
        
        start_cov = time.time()        
        total_coverage = get_coverage(coverage_file, region)
        end_cov = time.time()
        elapsed = end_cov - start_cov
        print("[{}] [{}] [S{}] [{}] [TOTAL_COVERAGE] {}".format(rank, sample_rank, sample_index, sample, total_coverage))
        
#         if sample_rank == 0:
#             interval_time = [sample, elapsed]
#         else:
#             interval_time = []
            
#         interval_times = comm.gather(interval_time)
        
#         if rank == 0:
#             interval_times = list(filter(lambda x: x is not None, interval_times))
#             for interval_time in interval_times:
#                 if len(interval_time) > 0:
#                     print("INTERVAL_TIME[0]", interval_time[0])
#                     print("INTERVAL_TIME", interval_time)
#                     
#                     chronometer[interval_time[0]]["coverage"] = interval_time[1]

        if sample_rank == 0:
            now = datetime.datetime.now().time()
            elapsed = time.time() - t1
            print("[SYSTEM] [TIME] [MPI] [0] MIDDLE-COVERAGE [now:{}] [elapsed: {}]".format(now, elapsed))
        
        # Collect all the files with the coverage
        files = []
        for file in os.listdir(coverage_dir):
            if region is not None and file != region[0]: continue
            if file.startswith("."): continue
            if file.endswith(".cov"): continue
            if file.endswith(".txt"): continue
            if file == "chrM": continue
            files.append(file)
        files.sort()
    
        if sample_rank == 0:
            print("[0] [S"+str(sample_index)+"] " + str(len(files)) + " FILES => " + str(files))
        
        # Master: dispatches the work to the other slaves
        if sample_rank == 0:
            start_intervals = t1
            print("[0] Start time: {}".format(start_intervals))
            
            done = 0
            total = len(files)
              
            queue = set()
            for i in range(1, min(sample_size, total+1)):
                file = files.pop()
                print("[SYSTEM] [MPI] [0] Sending coverage data "+ str(file) +" to rank " + str(i))
                sample_comm.send(file, dest=i, tag=CALCULATE_COVERAGE)
                queue.add(i)
      
            while len(files) > 0:
                status = MPI.Status()
                subintervals = sample_comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
                for subinterval in subintervals:
                    homeworks.append([sample_index] + list(subinterval))
      
                done += 1
                who = status.Get_source()
                queue.remove(who)
                now = datetime.datetime.now().time()
                elapsed = time.time() - start_intervals
                print("[SYSTEM] [TIME] [MPI] [0] [S{}] COVERAGE RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [#intervals: {}] [{}/{}][{:.2f}%] [Queue:{}]".format(sample_index, str(who), now, elapsed, len(homeworks), done, total, 100 * float(done)/total, queue))
      
                file = files.pop()
                print("[SYSTEM] [MPI] [0] [S"+str(sample_index)+"] Sending coverage data "+ str(file) +" to rank " + str(who))
                sample_comm.send(file, dest=who, tag=CALCULATE_COVERAGE)
                queue.add(who)
      
            while len(queue) > 0:
                status = MPI.Status()
                print("[SYSTEM] [MPI] [0] [S{}] Going to receive data from slaves.".format(sample_index))
                subintervals = sample_comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
                for subinterval in subintervals:
                    homeworks.append([sample_index] + list(subinterval))
      
                done += 1
                who = status.Get_source()
                queue.remove(who)
                now = datetime.datetime.now().time()
                elapsed = time.time() - start_intervals
                print("[SYSTEM] [TIME] [MPI] [0] [S{}] COVERAGE RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [#intervals: {}] [{}/{}][{:.2f}%] [Queue:{}]".format(sample_index, str(who), now, elapsed, len(homeworks), done, total, 100 * float(done)/total, queue))
            
            # Let them know we finished calculating the coverage
            for i in range(1, sample_size):
                sample_comm.send(None, dest=i, tag=STOP_WORKING)
            
            now = datetime.datetime.now().time()
            elapsed = time.time() - start_intervals
            
            interval_file = temp_dir + "/intervals.txt"
            print("[SYSTEM] [TIME] [MPI] [0] [S{}] SAVING INTERVALS TO {} [now:{}] [elapsed: {}]".format(sample_index, interval_file, now, elapsed))
            writer = open(interval_file, "w")
            writer.write(str(len(homeworks)) + "\n")
            for homework in homeworks:
                writer.write("\t".join([str(x) for x in homework[1:]]) + "\n")
            writer.close()
            
            now = datetime.datetime.now().time()
            elapsed = time.time() - start_intervals
#             interval_time = [sample, elapsed]
            print("[SYSTEM] [TIME] [MPI] [0] [S{}] INTERVALS SAVED TO {} [now:{}] [elapsed: {}]".format(sample_index, interval_file, now, elapsed))
            print("[SYSTEM] [TIME] [MPI] [0] [S{}] FINISHED CALCULATING INTERVALS [now:{}] [elapsed: {}]".format(sample_index, now, elapsed))
        else:
            
#             interval_time = []
            
            while(True):
                status = MPI.Status()
                # Here data is the name of a chromosome.
                data = sample_comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                if tag == CALCULATE_COVERAGE:
                    intervals = calculate_intervals(total_coverage, coverage_dir + data, region)
                    sample_comm.send(intervals, dest=0, tag=IM_FREE)
                elif tag == STOP_WORKING:
                    print("[SYSTEM] [TIME] [MPI] [SEND/RECV] [RECV] [{}] received STOP calculating intervals SIGNAL FROM RANK 0 [{}]".format(str(rank), datetime.datetime.now().time()))
                    break
        
#         interval_times = comm.gather(interval_time)
#         if rank == 0:
#             for interval_time in interval_times:
#                 if len(interval_time) > 0:
#                     chronometer[interval_time[0]]["intervals"] = interval_time[1]
        
    print("[{}] [{}] [S{}] [{}] BEFORE GATHER HOMEWORKS #TOTAL={} intervals".format(rank, sample_rank, sample_index, sample, len(homeworks)))
    
    # Wait for all intervals to be collected
    homeworks = comm.gather(homeworks)
    homeworks_by_sample = {}
    homeworks_done = {}
    if rank == 0:
        homeworks = reduce(lambda x,y: x+y, homeworks)
        shuffle(homeworks)
        
        # Divide intervals by samples
        for homework in homeworks:
            sample_id = homework[0]
            if sample_id not in homeworks_by_sample: homeworks_by_sample[sample_id] = []
            homeworks_by_sample[sample_id].append(homework)
        
        for sample_id in homeworks_by_sample:
            homeworks_done[sample_id] = 0
                            
        print("[{}] [{}] [S{}] #TOTAL={} (all intervals)".format(rank, sample_rank, sample_index, len(homeworks)))

    ###########################################################
    ######### COMPUTATION SECTION #############################
    ###########################################################
    
    if rank == 0:
        done = 0
        print("[SYSTEM] [TIME] [MPI] [0] REDItools STARTED. MPI SIZE (PROCS): {} [now: {}]".format(size, datetime.datetime.now().time()))
        
        t1 = time.time()
        
        print("Loading chromosomes' sizes!")
        chromosomes = OrderedDict()
        for line in open(size_file):
            (key, val) = line.split()[0:2]
            chromosomes[key] = int(val)
        print("Sizes:")
        print(chromosomes)
  
        total = len(homeworks)
        #print("[SYSTEM] [MPI] [0] HOMEWORKS", total, homeworks)
  
        start = time.time()
        
        print("[SYSTEM] [TIME] [MPI] [0] REDItools PILEUP START: [now: {}]".format(datetime.datetime.now().time()))
        
        queue = set()
        for i in range(1, min(size, total)):
            interval = homeworks.pop()
            
            print("[SYSTEM] [MPI] [SEND/RECV] [SEND] [0] Sending data "+ str(interval) +" to rank " + str(i))
            comm.send(interval, dest=i, tag=ALIGN_CHUNK)
            queue.add(i)
  
        while len(homeworks) > 0:
            status = MPI.Status()
            response = comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.datetime.now().time()
            elapsed = time.time() - start
            print("[SYSTEM] [TIME] [MPI] [SEND/RECV] [RECV] [0] RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, done, total, 100 * float(done)/total, queue))
            
            sample_done = samples[response[0]]
            sample_done = os.path.basename(sample_done)
            sample_done = ".".join(sample_done.split(".")[0:-1])
            print(response)
            duration = response[-1] - response[-2]
            chronometer[sample_done]["parallel"] += duration
            homeworks_done[response[0]] += 1
            if homeworks_done[response[0]] == len(homeworks_by_sample[response[0]]):
                print("[SYSTEM] [MPI] [COMPLETE] [{}] [{}] [{}] now:{}".format(sample_done, chronometer[sample_done]["parallel"], str(datetime.timedelta(seconds=chronometer[sample_done]["parallel"])), now))
            
            interval = homeworks.pop()
            print("[SYSTEM] [MPI] [SEND/RECV] [SEND] [0] Sending data "+ str(interval) +" to rank " + str(who))
            comm.send(interval, dest=who, tag=ALIGN_CHUNK)
            queue.add(who)
  
        while len(queue) > 0:
            status = MPI.Status()
            response = comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.datetime.now().time()
            elapsed = time.time() - start
            print("[SYSTEM] [TIME] [MPI] [SEND/RECV] [RECV] [0] RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, done, total, 100 * float(done)/total, queue))
            
            sample_done = samples[response[0]]
            sample_done = os.path.basename(sample_done)
            sample_done = ".".join(sample_done.split(".")[0:-1])
            duration = response[-1] - response[-2]
            chronometer[sample_done]["parallel"] += duration
            homeworks_done[response[0]] += 1
            if homeworks_done[response[0]] == len(homeworks_by_sample[response[0]]):
                print("[SYSTEM] [MPI] [COMPLETE] [{}] [{}] [{}] now:{}".format(sample_done, chronometer[sample_done]["parallel"], str(datetime.timedelta(seconds=chronometer[sample_done]["parallel"])), now))
            
            print("[SYSTEM] [MPI] [SEND/RECV] [SEND] [0] Sending DIE SIGNAL TO RANK " + str(who))
            comm.send(None, dest=who, tag=STOP_WORKING)

        t2 = time.time()
        elapsed = t2-t1
        print("[SYSTEM] [TIME] [MPI] [0] WHOLE PARALLEL ANALYSIS FINISHED. CREATING SETUP FOR MERGING PARTIAL FILES - Total elapsed time [{:5.5f}] [{}] [now: {}]".format(elapsed, t2, datetime.datetime.now().time()))

        #####################################################################
        ######### RECOMBINATION OF SINGLE FILES #############################
        #####################################################################
        for s in samples:
            
            s = os.path.basename(s)
            s = ".".join(s.split(".")[0:-1])
            
            little_files = []
            print("Scanning all files in "+args.temp_dir + s +" matching " + ".*")
            for little_file in glob.glob(args.temp_dir + s + "/*"):
                if little_file.endswith("chronometer.txt"): continue
                if little_file.endswith("files.txt"): continue
                if little_file.endswith("intervals.txt"): continue
                if little_file.endswith("status.txt"): continue
                if little_file.endswith("progress.txt"): continue
                if little_file.endswith("times.txt"): continue
                if little_file.endswith("groups.txt"): continue
                
                print(little_file)
                pieces = re.sub(r"\.[^\.]*$", "", os.path.basename(little_file)).split("#")
                pieces.insert(0, little_file)
                little_files.append(pieces)
    
            # Sort the output files
            keys = chromosomes.keys()
            print("[SYSTEM] "+str(len(little_files))+" FILES TO MERGE: ", little_files)
            little_files = sorted(little_files, key = lambda x: (keys.index(x[1]), int(x[2])))
            print("[SYSTEM] "+str(len(little_files))+" FILES TO MERGE (SORTED): ", little_files)
            
            smallfiles_list_filename = args.temp_dir + s + "/" + "files.txt"
            f = open(smallfiles_list_filename, "w")
            for little_file in little_files:
                f.write(little_file[0] + "\n")
            f.close()
            
            # Chronometer data
            chronometer_filename = args.temp_dir + "/" + "chronometer.txt"
            f = open(chronometer_filename, "w")
            #f.write("\t".join(["SampleID", "Coverage", "Intervals", "Editing", "Coverage (human)", "Intervals (human)", "Coverage (human)"]))
            for s in chronometer:
#                 coverage_duration = str(datetime.timedelta(seconds=chronometer[s]["coverage"]))
#                 interval_duration = str(datetime.timedelta(seconds=chronometer[s]["intervals"]))
                parallel_duration = str(datetime.timedelta(seconds=chronometer[s]["parallel"]))
                f.write("\t".join([
                    s,
#                     str(chronometer[s]["coverage"]),
#                     str(chronometer[s]["intervals"]),
                    str(chronometer[s]["parallel"]),
#                     coverage_duration,
#                     interval_duration,
                    parallel_duration]) + "\n")
            f.close()
        
        t2 = time.time()
        print("[SYSTEM] [TIME] [MPI] [0] [END] - WHOLE ANALYSIS FINISHED - Total elapsed time [{:5.5f}] [{}] [now: {}]".format(t2-t1, t2, datetime.datetime.now().time()))
        
    # Slave processes
    else:
        
        while(True):
            status = MPI.Status()
            data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)

            tag = status.Get_tag()
            if tag == ALIGN_CHUNK:

                # Process it
                time_start = time.time()
                time_s = datetime.datetime.now().time()
                print("[SYSTEM] [TIME] [MPI] [SEND/RECV] [RECV] [{}] REDItools: STARTED {} from rank 0 [{}] Interval: {}".format(str(rank), str(data), time_s, data))
                
                local_sample_index = data[0]
                local_sample_filepath = samples[local_sample_index]
                local_sample = os.path.basename(local_sample_filepath)
                local_sample = ".".join(local_sample.split(".")[0:-1])
                print("[LAUNCH REDITOOLS] {} {} {}".format(local_sample_index, local_sample_filepath, local_sample))
                
                id = data[1] + "#" + str(data[2]) + "#" + str(data[3])
                
                options["bamfile"] = local_sample_filepath
                options["region"] = [data[1], data[2], data[3]]
                options["output"] = args.temp_dir + local_sample + "/" + id + ".gz"
                
                print("[MPI] [" + str(rank) + "] COMMAND-LINE:", options)
                
                if not os.path.exists(args.temp_dir + local_sample + "/" + id + ".gz"):
                    gc.collect()
                    reditools.analyze(options)

                time_end = time.time()
                time_e = datetime.datetime.now().time()
                print("[SYSTEM] [TIME] [MPI] [{}] REDItools: FINISHED {} [{}][{}] [TOTAL:{:5.2f}]".format(str(rank), str(data), time_s, datetime.datetime.now().time(), time_end - time_start))

                print("[SYSTEM] [TIME] [MPI] [SEND/RECV] [SEND] [{}] SENDING IM_FREE tag TO RANK 0 [{}]".format(str(rank), datetime.datetime.now().time()))
                comm.send(data + [time_s, time_e, time_start, time_end], dest=0, tag=IM_FREE)
            elif tag == STOP_WORKING:
                print("[SYSTEM] [TIME] [MPI] [SEND/RECV] [RECV] [{}] received DIE SIGNAL FROM RANK 0 [{}]".format(str(rank), datetime.datetime.now().time()))
                break
            
    print("[{}] EXITING [now:{}]".format(rank, time.time()))
