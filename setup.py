from ctypes.util import find_library
import glob
import os
import pkgconfig
import setuptools
import sys

from setuptools.command.build_ext import build_ext

__version__ = "0.0.1"

BOOST_LIBRARIES = ["boost_system", "boost_filesystem", "boost_chrono",
                   "boost_regex"]

IS_DEV_BUILD = False

def check_build_dependencies():
    """Check that required libraries are present before trying to build."""
    for library in BOOST_LIBRARIES:
        if not find_library(library):
            raise Exception("Could not find {} library".format(library))

    if (not pkgconfig.exists("hdf5")) and (not find_library("hdf5")):
        raise Exception("Could not find hdf5")

    if not pkgconfig.exists("eigen3"):
        raise Exception("Could not find hdf5")

def get_libraries():
    """Get libraries for building EM2.

    This is mostly just finding hdf5. Boost doesn't seem to work with
    pkg-config.
    """

    all_libraries = list(BOOST_LIBRARIES)
    if pkgconfig.exists("hdf5"):
        all_libraries.extend(pkgconfig.parse("hdf5")['libraries'])
    else:
        all_libraries.append("hdf5")
    all_libraries.append("hdf5_cpp")

    return all_libraries

def get_library_dirs():
    """Get any library dirs to pass to the compiler."""

    library_dirs = pkgconfig.parse("hdf5")['library_dirs']
    return library_dirs

def get_include_dirs():
    """Get includes to pass to the compiler."""
    include_dirs = []
    include_dirs.extend(pkgconfig.parse("eigen3")["include_dirs"])
    include_dirs.extend(pkgconfig.parse("hdf5")["include_dirs"])
    include_dirs.append(get_pybind_include())
    include_dirs.append(get_pybind_include(user=True))
    return include_dirs


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']


    def build_extensions(self):
        # Avoid lots of warnings
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except ValueError:
            pass


        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())

        # Opts for building EM2 from Paolo
        try:
            self.compiler.compiler_so.remove("-O2")
        except ValueError:
            pass
        if has_flag(self.compiler, '-msse4.2'):
            opts.append("-msse4.2")
        opts.append("-O3")

        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

check_build_dependencies()

ext_modules = [
    setuptools.Extension(
        'pyEM2',
        sources=glob.glob(os.path.join("src", "pyEM2*.cpp")) +
        (glob.glob("ExpressionMatrix2/src/*.cpp") if not IS_DEV_BUILD else []),
        include_dirs=get_include_dirs() + [os.path.join("ExpressionMatrix2", "src")],
        library_dirs=(["ExpressionMatrix2"] if IS_DEV_BUILD else []) + get_library_dirs(),
        libraries=(["em2"] if IS_DEV_BUILD else []) + get_libraries(),
        language="c++")
]

setuptools.setup(
    name='pyEM2',
    version=__version__,
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)
