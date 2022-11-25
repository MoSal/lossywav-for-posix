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

## Building

A simple `Makefile` is available. Run `make` to build the `lossywav` executable,
`make clean` to clean up. `lossywav` can be moved to e.g. `~/.local/bin`.

Note: the Python build scripts from the original repo are removed, due to not
functioning on Debian 11.

# Usage

See output of `lossywav -h`, reproduced here in full:

```
lossyWAV 1.4.2, Copyright (C) 2007-2016 Nick Currie. Copyleft.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful,but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

Usage   : lossyWAV <input wav file> <options>

Example : lossyWAV musicfile.wav

Quality Options:

-q, --quality <t>    where t is one of the following (default = standard):
    I, insane        highest quality output, suitable for transcoding;
    E, extreme       higher quality output, suitable for transcoding;
    H, high          high quality output, suitable for transcoding;
    S, standard      default quality output, considered to be transparent;
    C, economic      intermediate quality output, likely to be transparent;
    P, portable      good quality output for DAP use, may not be transparent;
    X, extraportable lowest quality output, probably not transparent.

Standard Options:

-C, --correction     write correction file for processed WAV file; default=off.
-f, --force          forcibly over-write output file if it exists; default=off.
-h, --help           display help.
-L, --longhelp       display extended help.
-M, --merge          merge existing lossy.wav and lwcdf.wav files.
-o, --outdir <t>     destination directory for the output file(s).
-v, --version        display the lossyWAV version number.
-w, --writetolog     create (or add to) lossyWAV.log in the output directory.

Special thanks go to:

David Robinson       for the publication of his lossyFLAC method, guidance, and
                     the motivation to implement his method as lossyWAV.

Horst Albrecht       for ABX testing, valuable support in tuning the internal
                     presets, constructive criticism and all the feedback.

Sebastian Gesemann   for the adaptive noise shaping method and the amount of
                     help received in implementing it and also for the basis of
                     the fixed noise shaping method.

Tyge Lovset          for the C++ translation initiative.

Matteo Frigo and     for libfftw3-3.dll contained in the FFTW distribution
Steven G Johnson     (v3.2.1 or v3.2.2).

Mark G Beckett       for the Delphi unit that provides an interface to the
(Univ. of Edinburgh) relevant fftw routines in libfftw3-3.dll.

Don Cross            for the Complex-FFT algorithm originally used.
```

## Example

If you have a file hierachy containing your music stored in the FLAC format, you can use the following script in the directory of the library to convert the collection to LossyFLAC, at standard LossyWAV processing, and standard FLAC compression:

```
for i in **/*.flac; do
  flac -d "$i" --stdout --silent|lossywav - --stdout --quality standard --stdinname -- |flac - -b 512 -o "${i%.*}.lossy.flac" --silent &&  metaflac --export-tags-to=- "$i" | metaflac --import-tags-from=- "${i%.*}.lossy.flac"
done
```

# Credits

* All lossyWAV authors and contributors.
* HydrogenAudio community.
