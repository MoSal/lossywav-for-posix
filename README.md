# Description
lossyWAV is a near-lossless audio processor which dynamically reduces the
bitdepth of the signal on a block-by-block basis. Bitdepth reduction adds noise
to the processed output. The added noise is adaptively shaped by default and can
alternatively be fixed noise shaped or white noise depending on command line
parameters.

When lossyWAV-processed output is compressed with supported lossless codecs, the
bitrate of the output file is significantly reduced compared to the lossless
original.

# Supported lossless codecs
 * FLAC
 * Wavpack
 * TAK
 * LPAC
 * MPEG-4 ALS
 * WMA Lossless

# Authors
Nick Currie: principal developer.

Special thanks go to:
* David Robinson for the publication of his lossyFLAC method, guidance, and the motivation to implement his method as lossyWAV.
* Horst Albrecht for ABX testing, valuable support in tuning the internal presets, constructive criticism and all the feedback.
* Sebastian Gesemann for the adaptive noise shaping method and the amount of help received in implementing it and also for the basis of the fixed noise shaping method.
* Tyge Lovset for the C++ translation initiative.
* Matteo Frigo and Steven G Johnson for libfftw3-3.dll contained in the FFTW distribution (v3.2.1 or v3.2.2).
* Mark G Beckett (Univ. of Edinburgh) for the Delphi unit that provides an interface to the relevant fftw routines in libfftw3-3.dll.
* Don Cross for the Complex-FFT algorithm originally used.

# Discussion
 * v1.4.0 discussion: http://www.hydrogenaud.io/forums/index.php?showtopic=107081
