# Note
This is the `README` of **lossyWAV for POSIX**.
The README of lossyWAV (proper) was moved to  `README.lossyWAV.md`.

# lossyWAV for POSIX

This is an attempt to build lossyWAV natively in POSIX systems.

The Windows API is used in the official lossyWAV sources. But not extensively.
lossyWAV is only offered as Windows binaries or source files which can't be
compiled and linked to run in any non-Windows environment natively.

Luckily, I was able to replace the instances where the Windows API is used
with C++11-compliant or POSIX-compliant code that provides similar
functionality.

None of the changes would be enabled at the preprocessing stage if `_WIN32` is
defined.

**Tested platforms**: GNU/Linux, OSX, FreeBSD, NetBSD, OpenBSD(needs latest eg++).

# How to Build and Install

## Mac OSX

Homebrew users can install lossywav with a simple command:

     brew install --HEAD https://raw.githubusercontent.com/MoSal/lossywav-for-posix/master/lossywav.rb

## General

Building lossywav should be as simple as:

```
./waf configure [OPTIONS]
./waf build [OPTIONS]
./waf install [OPTIONS]
```

A typical example would be:
```
./waf configure --prefix=/usr --enable-fftw3
./waf build
# As root
./waf install --destdir=/
```

**Note:** If there is no `python` in `PATH`. You can invoke `waf` with whatever python
executable you have. Both *Python 2* and *Python 3* are supported. For example:

     python2.7 waf configure

---

The output of `./waf -h` including all OPTIONS:

```
waf [commands] [options]

Main commands (example: ./waf build -j4)
  build    : executes the build
  clean    : cleans the project
  configure: configures the project
  dist     : makes a tarball for redistributing the sources
  distcheck: checks if the project compiles (tarball from 'dist')
  distclean: removes the build directory
  install  : installs the targets on the system
  list     : lists the targets to execute
  step     : executes tasks in a step-by-step fashion, for debugging
  uninstall: removes the targets installed

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -c COLORS, --color=COLORS
                        whether to use colors (yes/no/auto) [default: auto]
  -j JOBS, --jobs=JOBS  amount of parallel jobs (2)
  -k, --keep            continue despite errors (-kk to try harder)
  -v, --verbose         verbosity level -v -vv or -vvv [default: 0]
  --zones=ZONES         debugging zones (task_gen, deps, tasks, etc)

  Configuration options:
    -o OUT, --out=OUT   build dir for the project
    -t TOP, --top=TOP   src dir for the project
    --prefix=PREFIX     installation prefix [default: '/usr/local/']
    --bindir=BINDIR     bindir
    --libdir=LIBDIR     libdir
    --check-cxx-compiler=CHECK_CXX_COMPILER
                        list of C++ compilers to try [g++ clang++ icpc]
    --enable-compiler-warnings
                        Enable compiler warnings. (default: False)
    --werror            Consider warnings fatal. (default: False)
    --disable-compile-optimizations
                        Don't check/set compile optimization flags. (default: False)
    --disable-link-optimizations
                        Don't check/set link optimization flags. (default: False)
    --disable-lto       Don't check/set lto flags. (default: False)
    --enable-debug      Set debug flags. (default: False)
    --enable-fftw3      Compile and link against libfftw3. (default: False)
    --fftw3-cxxflags=FFTW3_CXXFLAGS
                        Skip pkg-config and set fftw3 cxxflags explicitly (default: None)
    --fftw3-libs=FFTW3_LIBS
                        Skip pkg-config and set fftw3 libs explicitly (default: None)

  Build and installation options:
    -p, --progress      -p: progress bar; -pp: ide output
    --targets=TARGETS   task generators, e.g. "target1,target2"

  Step options:
    --files=FILES       files to process, by regexp, e.g. "*/main.c,*/test/main.o"

  Installation and uninstallation options:
    --destdir=DESTDIR   installation root [default: '']
    -f, --force         force file installation
    --distcheck-args=ARGS
                        arguments to pass to distcheck
```

A simple `Makefile.unix` is also available as a last resort alternative.

# Credits
* All lossyWAV authors and contributors.
* HydrogenAudio community.
