import os
from subprocess import call

def convert_e2p(fdir='.'):
    """
    """
    if not fdir.endswith('/'): fdir += '/'
    dlist = os.listdir(fdir)
    for f in dlist:
        if f.endswith('.eps'):
            print "converting ", f, " to pdf"
            call(['epstopdf',fdir + f])
        else: continue

def rm_eps(fdir='.'):
    dlist = os.listdir(fdir)
    if not fdir.endswith('/'): fdir += '/'
    for f in dlist:
        if f.endswith('.eps'):
            print "removing ", f
            call(['rm', fdir + f])
        else: continue
