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

None of the changes would be enabled at the preprocessing stage if _WIN32 is
defined. Also, using libfftw3 is enabled with a simple define at build time.

Build-time checks and improvements are left for whoever is interested and
experienced in that area.

# How to Build and Install
A simple `Makefile` is offered. Besides the default(*all*) target,
The following targets are offered:

* **optimized**: Add optimization flags as follows:
  * *CXXFLAGS*:  -march=native -Ofast -flto
  * *LDFLAGS*:  -flto

* **fftw**: Enable libfftw3. Needs libfftw3 header and library.

* **fftw-optimized**: *fftw* + *optimized*

Passed **CXXFLAGS** and **LDFLAGS** (from shell or make) will be appended.

No **install** target is provided. You would only need one file(`lossywav`)
anyway.

A **clean** target is available to remove all object files and the `lossywav` 
executable.

# Warning
I take no responsibility for the patches offered in this repository.
I have near-zero experience with C++ and the Windows API. As I usually
only code in C/C99 in POSIX environments.

# Credits
* All lossyWAV authors and contributors.
* HydrogenAudio community.
