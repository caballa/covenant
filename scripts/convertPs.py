#!/usr/bin/env python

import sys
import os, os.path
import subprocess
from subprocess import Popen, PIPE

def printf(format, *args):
    sys.stdout.write(format % args)

def convert_ps(infile):
    outfile = infile + ".ps"
    ret = subprocess.call(["dot", "-Tps", infile,"-o", outfile], stdout=None, stderr=None)
    if ret != 0:
        raise IOError("command dot failed")
    return
    
def main (argv):
    for root, dirs, files in os.walk("./"):
        for file in files:
            if file.endswith(".dot"):
                fullname = os.path.join(root, file)
                convert_ps(fullname)


if __name__ == "__main__":
    try:
        main (sys.argv)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
    except: 
        print "Unexpected error:", sys.exc_info()[0]
        raise
