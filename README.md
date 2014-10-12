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

# Discussion
 * v1.4.0 discussion: http://www.hydrogenaud.io/forums/index.php?showtopic=107081
