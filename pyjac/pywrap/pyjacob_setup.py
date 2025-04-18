from distutils.core import setup, Extension
import distutils.ccompiler

from Cython.Distutils import build_ext
import parallel_compiler as pcc
import numpy
import os

sources = ['/home/syoukera/github/pyJac_syoukera/pyjac/pywrap/pyjacob_wrapper.pyx']
includes = ['out_nh3_edit/']

distutils.ccompiler.CCompiler.compile = pcc.parallel_compile

ext_modules=[Extension("pyjacob",
     sources=sources,
     include_dirs=includes + [numpy.get_include()],
     extra_compile_args=['-frounding-math', '-fsignaling-nans'],
     language='c',
     extra_objects=[os.path.join('build/temp.linux-x86_64-3.6', '/home/syoukera/github/pyJac_syoukera/build/temp.linux-x86_64-3.6/libc_pyjac.a')]
     )]

setup(
    name='pyjacob',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext}
)
