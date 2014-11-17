#!/usr/bin/env python

# Script to run in parallel different configurations (profiles) of
# covenant

import sys
import os
import os.path
import atexit
import tempfile
import shutil
import subprocess as sub
import threading
import signal
import resource

root = os.path.dirname (os.path.dirname (os.path.realpath (__file__)))
verbose = True

def initProfiles():
    base = []
    profiles = dict()
    profiles ['greedy'] = base + [ '--gen=greedy']
    profiles ['max_gen'] = base + ['--gen=max-gen']
    return profiles

profiles = initProfiles ()

def listProfiles ():
    for (k, v) in profiles.iteritems ():
        print k, ':', ' '.join (v)

def isexec (fpath):
    if fpath == None: return False
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def parseOpt (argv):
    from optparse import OptionParser

    parser = OptionParser ()
    parser.add_option ("--save-temps", dest="save_temps",
                       help="Do not delete temporary files",
                       action="store_true",
                       default=False)
    parser.add_option ("--temp-dir", dest="temp_dir",
                       help="Temporary directory",
                       default=None)
    parser.add_option ('--cpu', type='int', dest='cpu',
                       help='CPU time limit (seconds)', default=-1)
    parser.add_option ('--mem', type='int', dest='mem',
                       help='MEM limit (MB)', default=-1)
    parser.add_option ('--profiles', '-p', dest='profiles',
                       default='max_gen:greedy',
                       help='Colon separated list of profiles')
    parser.add_option ('--list-profiles', dest='list_profiles',
                       action='store_true', default=False)

    (options, args) = parser.parse_args (argv)

    if options.list_profiles:
        listProfiles ()
        sys.exit (0)

    return (options, args)

def createWorkDir (dname = None, save = False):
    if dname == None:
        workdir = tempfile.mkdtemp (prefix='covenantpar-')
    else:
        workdir = dname

    if verbose: print "Working directory", workdir

    if not save: atexit.register (shutil.rmtree, path=workdir)
    return workdir

def getCovenant (): 
    covenant = os.path.join (root, "covenant")
    if not isexec (covenant):
        covenant = os.path.join (root, "build/tools/covenant")
    if not isexec (covenant):
        raise IOError ("Cannot find covenant")
    return covenant

def cat (in_file, out_file): out_file.write (in_file.read ())

running = list()

def printf(format, *args):
    sys.stdout.write(format % args)

def kill (proc):
    try:
        printf("TIMEOUT\n")
        proc.terminate ()
        proc.kill ()
        proc.wait ()
        global running_process
        running_process = None
    except OSError:
        pass

def runCovenant (args, fname, stdout, stderr, cpu=-1, mem=-1):
    def set_limits ():
        if mem > 0:
            mem_bytes = mem * 1024 * 1024
            resource.setrlimit (resource.RLIMIT_AS, [mem_bytes, mem_bytes])

    in_name = fname
    if not os.path.isfile (in_name):
        raise IOError ("Cannot find input file")

    args = args + [in_name]
    if verbose: print ' '.join (args)
    p =  sub.Popen (args, preexec_fn=set_limits,
                    stdout=open (stdout, 'w'),
                    stderr=open (stderr, 'w'))

    # global running_process
    # running_process = p
    # timer = threading.Timer (cpu, kill, [p])
    # if cpu > 0: timer.start ()
    # try:
    #     (pid, returnvalue, ru_child) = os.wait4 (p.pid, 0)
    #     running_process = None
    # finally:
    #     ## kill the timer if the process has terminated already
    #     if timer.isAlive (): timer.cancel ()
    # ## if it did not terminate properly, propagate this error code
    # if returnvalue != 0: sys.exit (returnvalue)

    return p

def run (workdir, fname, covenant_args = [], profs = [],
         cpu=-1, mem=-1):

    sys.stdout.flush ()

    base_args = [getCovenant ()]
    base_args.extend (covenant_args)

    conf_name = list ()
    ufo = list ()

    for prof in profs:
        conf_name.append (prof)
        p_args = base_args + profiles [prof]
        ufo.append (p_args)

    name = os.path.splitext (os.path.basename (fname))[0]
    stdout = [os.path.join (workdir, name + '_covenant{}.stdout'.format (i))
              for i in range(len (ufo))]
    stderr = [os.path.join (workdir, name + '_covenant{}.stderr'.format (i))
              for i in range (len (ufo))]

    global running
    running.extend ([runCovenant (ufo [i], fname, stdout[i], stderr [i], cpu, mem)
                     for i in range (len (ufo))])


    orig_pids = [p.pid for p in running]
    pids = [p.pid for p in running ]
    pid = -1
    returnvalue = -1
    while len (pids) != 0:
        print 'Running: ', pids

        (pid, returnvalue, ru_child) = os.wait4 (-1, 0)

        print 'Finished pid {0} with'.format (pid),
        print ' code {0} and signal {1}'.format((returnvalue // 256),
                                                (returnvalue % 256))
        pids.remove (pid)

        if returnvalue == 0:
            for p in pids:
                try:
                    os.kill (p, signal.SIGTERM)
                except OSError: pass
                finally:
                    try:
                        os.waitpid (p, 0)
                    except OSError: pass
            break

    if returnvalue == 0:
        idx = orig_pids.index (pid)
        cat (open (stdout [idx]), sys.stdout)
        cat (open (stderr [idx]), sys.stderr)

        print 'WINNER: ', ' '.join (ufo [idx])
        print 'BRUNCH_STAT config {0}'.format (idx)
        print 'BRUNCH_STAT config_name {0}'.format (conf_name [idx])

    else:
        print "ALL INSTANCES FAILED"
        print 'Calling sys.exit with {}'.format (returnvalue // 256)
        sys.exit (returnvalue // 256)

    running[:] = []
    return returnvalue

def covenant_opt (x):
    if x.startswith ('-'):
        y = x.strip ('-')
        return y.startswith ('solutions') or  y.startswith ('abs')
    return False

def non_covenant_opt (x): return not covenant_opt (x)


def main (argv):
    covenant_args = filter (covenant_opt, argv [1:])
    argv = filter (non_covenant_opt, argv [1:])

    (opt, args) = parseOpt (argv)

    workdir = createWorkDir (opt.temp_dir, opt.save_temps)
    returnvalue = 0
    for fname in args:
        returnvalue = run (workdir, fname, covenant_args, opt.profiles.split (':'),
                           opt.cpu, opt.mem)
    return returnvalue

def killall ():
    global running
    for p in running:
        try:
            if p.poll () == None:
                p.terminate ()
                p.kill ()
                p.wait ()
                # no need to kill pg since it kills its children
        except OSError:   pass
    running[:] = []

if __name__ == '__main__':
    # unbuffered output
    sys.stdout = os.fdopen (sys.stdout.fileno (), 'w', 0)
    try:
        signal.signal (signal.SIGTERM, lambda x, y: killall())
        sys.exit (main (sys.argv))
    except KeyboardInterrupt: pass
    finally: killall ()
