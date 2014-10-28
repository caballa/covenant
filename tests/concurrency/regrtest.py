#!/usr/bin/env python

# Usage: ./regrtest.py

import sys
import os, os.path
import subprocess
from subprocess import Popen, PIPE, call
import csv

brunch_cmmd = "../../scripts/brunch.py"
covenant = "../../build/tools/covenant"
timeout="5"

def printf(format, *args):
    sys.stdout.write(format % args)

def run_brunch(files, timeout, outdir, tool_cmmd, tool_args):
    args = []
    args.append(brunch_cmmd)
    args.append(" ")
    for f in files:
        args.append(f)
        args.append(" ")
    args.append(" --mem 4096")
    args.append(" --cpu " + timeout)
    args.append(" --out " + outdir)
    args.append(" --format base:Cpu:Result:Status")
    args.append(" -- " + tool_cmmd + " " + tool_args)

    strcmmd = ''.join(str(e) for e in args)

    ret = call(strcmmd, stdout=None, stderr=None, shell=True)
                          
    if ret != 0:
        printf("Error during %s\n",args)
        sys.exit(2)        

    return


def show_results(_f, print_ko, print_err):
    bench_data = {}
    with open(_f, 'rb') as csvfile:
        reader = csv.DictReader(csvfile)
        err, unsafe, safe, total_bench, imprecise, incorrect = 0, 0, 0, 0, 0, 0
        timeout, segfault, various, memory_out = 0, 0, 0, 0
        solved_safe = []
        solved_unsafe = []
        unsolved = []
        for row in reader:
            total_bench +=1
            if row['Result'] == "UNSAT":
                if 'safe' in row['base']:
                    safe +=1
                    solved_safe.append(row["base"])
                    bench_data.update({row['base']: row['Cpu']})
                else:
                    incorrect +=1
                    if print_ko: print row['File'] + " " + row['Result'] + ' --> KO'
            if row['Result'] == "SAT":
                if 'unsafe' in row['base']:
                    unsafe +=1
                    solved_unsafe.append(row["base"])
                    bench_data.update({row['base']: row['Cpu']})
                else:
                    imprecise +=1
                    if print_ko: print row['File'] + " " + row['Result'] + ' --> KO'
            if row['Result'] == "ERR":
                err +=1
                bench_data.update({row['base']: '50'})
            if row["Status"] == "-9":
                if print_err: print row['File']  + ' --> TIMEOUT'
                timeout += 1
            if row["Status"] == "-6":
                if print_err: print row['File']  + ' --> SEGFAULT'
                segfault += 1
            if row["Status"] == "-11":
                if print_err: print row['File']  + ' --> MEMORY-OUT'
                memory_out += 1
            if row["Result"] == "ERR" and (row["Status"] == "0" or row["Status"] == "1"):
                if print_err: print row['File']  + ' --> VARIOUS'
                various += 1
    print "\n\n==========SUMMARY of " +  _f  + " =========="
    print "TOTAL N. BENCHMARK: " + str(total_bench)
    print "SOLVED SAFE  : " + str(safe)
    print "SOLVED UNSAFE: " + str(unsafe)
    print "NOT SOLVED   : " + str(err)
    print "IMPRECISE (SAFE -> UNSAFE): " + str(imprecise)
    print "UNSOUND (UNSAFE -> SAFE)  : " + str(incorrect)
    print "ERROR SEGFAULT  : " + str(segfault)
    print "ERROR TIMEOUT   : " + str(timeout)
    print "ERROR MEMORY-OUT: " + str(memory_out)
    print "ERROR VARIOUS   : " + str(various)
    print "==============================="
    return bench_data, solved_safe, solved_unsafe

    
def main (argv):

    infiles= []
    for root, dirs, files in os.walk("./"):
        for file in files:
            if file.endswith(".cfg"):
                fullname = os.path.join(root, file)
                infiles.append(os.path.abspath(fullname))

    outdir =  "cycle-breaking--greedy"         
    tool_args = "-a cycle-breaking -g greedy"
    run_brunch(infiles, timeout, outdir , covenant, tool_args)
    show_results(os.path.join(outdir,'stats'), False, False)

if __name__ == "__main__":
    try:
        main (sys.argv)
    except: 
        print "Unexpected error:", sys.exc_info()[0]
        raise
