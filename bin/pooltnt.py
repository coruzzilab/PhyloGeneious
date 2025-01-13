#!/usr/bin/env python
#
import sys
import time
import multiprocessing 
import subprocess
import shlex
import os
import glob
mydir = ""
runpth = ""
def do_align(subdir):
    print "starting align for %d\n" % (subdir)
    tstdir = "%s%d" % (mydir,subdir)
    if not os.access(tstdir,os.F_OK):
        print tstdir, "does not exist\n"
        return 1
    tstdir += "/oid.tre"
    if os.access(tstdir,os.F_OK):
        print tstdir, "already done, skipping\n"
        return 0
    atstr = "%sorthologid.pl -at \"^%d$\"" % (runpth,subdir)
    
    args = []
    args = shlex.split(atstr)
    print "str: %s \n" % atstr
    print "args",args
    p = subprocess.Popen(args)
    sys.stdout.flush()
    rc = p.wait()
    print "end align %s rc %d\n" % (subdir,rc)
    sys.stdout.flush()
## end do_align
mindir = 0
maxdir = 0
if __name__ == '__main__':
    
    print  "time is %s cpu %d \n" % (time.ctime(),time.clock())
    nargs = len(sys.argv)
    if nargs>2:
        maxdir = int(sys.argv[2])
        mindir = int(sys.argv[1])
    if nargs<2 or maxdir<mindir:
        print "pooltnt.py min max\n"
        print "must have range of dir to run\n"
        exit(1)
    mydir = os.getcwd() + '/data/'
    ncpu = multiprocessing.cpu_count()
    runpth = os.environ["OID_HOME"] + '/bin/'
    print "$OID_BIN ",runpth
    print "we were called with %d cpus\n" % (ncpu)
    if nargs>3:
        ncpu = int(sys.argv[3])
        print "we were explicitly given %d cpus to use\n" % (ncpu)
#    dirlst = list(range(mindir,maxdir+1)) # set up list of processes to run
    dirlst = list(range(maxdir,mindir-1,-1)) # set up list of processes to run
# highest is smallest and fastest
    p = multiprocessing.Pool(ncpu)
    poolrc = p.map(do_align,dirlst)
    p.close()
    p.join()
    prstr = "time is %s cpu %d end \n" % (time.ctime(),time.clock())
    print "pool: ",prstr
    ou = open ("log/job/pooldone",'w')
    print >>ou, "min %d, max %d\n" % (mindir,maxdir)
    print >>ou, prstr
    ou.close()
    sys.stdout.flush()
    


