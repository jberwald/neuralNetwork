######################################################
# fileio.py
#
# Author: Jessse Berwald
#
# Wrapper around some simple file IO functions
#
######################################################

import os

slash = "/"

def savedir():
    uname = os.uname()[1]

    if uname == 'dell':
        prefix = '/home/jb-book2/'
        #savedir = '/home/jb-book2/Projects/neuralNet/temp/'
        savedir = '/home/data/'
    #'/var/tmp/search/'
    elif uname == 'simplex':
        prefix = '/home/jberwald/'
        savedir = '/data1/jberwald/'
    elif uname.endswith('rutgers.edu'):
        prefix = '/home/jberwald/'
        savedir = '/data/jberwald/neuralNet/'
    #'/home/jberwald/Projects/neuralNet/results/'

    return prefix, savedir

def mkdir(dir):
    """
    Make a directory if it does not exist. Do nothing otherwise.
    """
    try:
        os.mkdir(dir)
    except:
        print "could not make dir!" 

def make_sim_folder(fdir, **args):
    """
    Use str.join to create file name from params. Check if file
    exists. If it does, either prompt for continution or increment run
    number by 1.
    """
    nodes = str(args['num_ex']+args['num_in'])
    graph = str(args['graphnum'])
    eps = str(args['epsvec'])
    conn = str(args['conn'])
    simnum = str(args['sim_num'])

    plist = ['n', nodes, 'g', graph, 'e', eps, 'c', conn, 's']
    fname = ''.join(plist)
#     dlist = os.listdir(fdir)
    
#     if fname not in dlist:
#         os.mkdir(fdir+fname)

    return slash.join( [fdir, fname] )
    
def check_dir(fdir):
    """
    Check whether a fname exists in fdir. If not, call makedir().
    """
    return os.path.lexists(fdir)
    
