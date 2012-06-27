#!/usr/bin/python
"""
Simple module that reads text file into a numpy array and subsamples
that data at a given rate. Additional options are to compute the dfa
or mse on the data.
"""
import numpy
import os, sys, optparse
import subprocess as sp
from jb.neural_network.netutils import utils

D = utils.Data_Handlers()

dot = "."
slash = "/"

def subsample_data( dfile, srate, ndim=1 ):
    """
    Subsample the data in array stored in dfile.

    dfile :      Full path to data

    srate :      Take every <srate> value

    Writes subsampled data to <dfile>.sr<srate>.out

    Returns name of subsampled data file
    """
    data = numpy.fromfile( dfile, sep='\n' )

    if ndim > 1:
        data.resize( (len(data)/ndim, ndim) )

    sd = data[::srate]

    if dfile.endswith('.out') or dfile.endswith('.dat'):
        savename = dfile[:-3] + "sr" + str(srate) +\
                   dot + "out"
    else:
        savename = dfile + "sr" + str(srate) +\
                   dot + "out"
    numpy.savetxt( savename, sd )
    del data

    return savename

def run_dfa( dfile ):
    """
    Call dfa from utils.py
    """
    D.run_dfa( dfile )

def run_mse( dfile, **args ):
    """
    mse -m 2 -M 4 -b 1 -r 0.15 -R 0.2 -c 0.01 <nsr040.rr >nsr040.mse
    """
    fargs = {'-m': 2,
             '-r': 0.15 }
    for k in args:
        farg[k] = args[k]

    cmd = ["mse"]
    for k in fargs:
        cmd.append( k )
        cmd.append( str(fargs[k]) )

    if dfile.endswith('out') or dfile.endswith('dat'):
        outfile = dfile[:-3] + 'mse'
    else:
        outfile = dfile + '.mse'
    fin = open( dfile, 'r' )
    fout = open( outfile, 'wb' )
    sp.call( cmd, stdin=fin, stdout=fout )

if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.usage = "python subsample.py [options]"
        
    data_help = "Path to data file containing array data."
    srate_help = "Subsampling rate. [4]"
    dfa_help = "Call dfa.c on subsampled data. [False]"
    mse_help = "Call mse.c on subsampled data. All options set by "\
               "hand in run_mse()[False]"

    parser.add_option("--data", "-d",
                      help=data_help,
                      type="string",
                      action="store",
                      dest="data",
                      default=None)
    parser.add_option("--srate", "-s",
                      help=srate_help,
                      type="int",
                      action="store",
                      dest="srate",
                      default=4)
    parser.add_option("--dfa",
                      help=dfa_help,
                      action="store_true",
                      dest="dfa",
                      default=False)
    parser.add_option("--mse",
                      help=mse_help,
                      action="store_true",
                      dest="mse",
                      default=False)

    global options

    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        print parser.print_help()
    else:
       sdata = subsample_data( options.data, options.srate )
       if options.dfa:
           print "running dfa..."
           run_dfa( sdata )
       if options.mse:
           print "running mse..."
           run_mse( sdata )
