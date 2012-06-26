#----------------------------------------
# 
# setup.py 
#
# For use with neurons.c. Two source files: neurons.c and
# fitzhugh_nagumo.c
#
#----------------------------------------
'''python setup.py build will build the extension module 
fitzhugh_nagumo.so in ./build/lib.arch-id/

NOTE: On Mac OSX, the compile flags are incorrect. In particular, distutils automatically passes

'-fno-strict-aliasing -fno-common -dynamic -arch i386  -DNDEBUG -g -O3  -arch i386  '

to gcc. This forces the wrong arch for many dynamically linked
executables, which are built with 64bit arch. This is fixed below with
the hackish code forcing the sysconfig variables to their correct
values.
                                                                              
'''
from distutils.core import setup, Extension
from distutils.util import get_platform
from distutils import sysconfig
import numpy
import sys

space = ' '

# handles macport's default location
if sys.platform == 'darwin':
       GSL_INCLUDE = '/opt/local/include/'
       GSL_LIB = '/opt/local/lib/'
       ARCH = 'x86_64'
       compile_args = ['-O3', '-DHAVE_INLINE', '-msse3' ]#, '-march=intel']
       
       # vars = sysconfig.get_config_vars()

       # # change OPT
       # opt = vars['OPT']
       # optsplit = opt.split( ' ' )
       # idx = optsplit.index( '-arch' )
       # optsplit[ idx+1 ] = 'x86_64'
       # new_opt = space.join( optsplit )
       # vars['OPT'] = new_opt

       # # change PY_CFLAGS
       # for key in ['PY_CFLAGS', 'CFLAGS']:
       #        opt = vars[ key ]
       # optsplit = opt.split( ' ' )
       # for arg in optsplit:
       #        if arg == '-arch':
       #               idx = optsplit.index( '-arch' )
       #               optsplit[ idx+1 ] = 'x86_64'
       #               new_opt = space.join( optsplit )
       #               vars[ key ] = new_opt
       
# sets paths to typical linux locations
elif sys.platform == 'linux2':
       GSL_INCLUDE = '${HOME}/local/include/'
       GSL_LIB = '${HOME}/local/lib/'
       ARCH = ''
       compile_args = ['-O3', '-DHAVE_INLINE', '-msse3'] #'-march=native']

       print 'home', GSL_INCLUDE

fname = "neurons_file_cells"

includen=[numpy.get_include(), GSL_INCLUDE] # or /usr/include/gsl

print "COMPILE ARGS", compile_args

module1 = Extension( fname, include_dirs=includen,
                    library_dirs=[GSL_LIB],
                    libraries=['gsl', 'gslcblas'],
                    extra_compile_args = compile_args, # ['-O3', '-DHAVE_INLINE'], ##, ARCH],  ##'-msse3'], #, '-march=native'],
                    sources = [ 'fitzhugh_nagumo.c', fname+'.c' ] )

setup (name = fname,
       version = 'epsilon',
       description = 'C api module for neurons*.py',
       ext_modules = [module1])


