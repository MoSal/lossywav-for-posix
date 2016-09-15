/**===========================================================================

    lossyWAV: Added noise WAV bit reduction method by David Robinson;
              Noise shaping coefficients by Sebastian Gesemann;

    Copyright (C) 2007-2016 Nick Currie, Copyleft.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: lossywav <at> hotmail <dot> co <dot> uk

==============================================================================
    Initial translation to C++ from Delphi
    Copyright (C) Tyge Løvset (tycho), Aug. 2012
===========================================================================**/

#ifndef nWav_H_
#define nWav_H_

void nWAV_Init();

void nWAV_Cleanup();

void MergeFiles();                 // Merge lossy.wav and lwcdf.wav files.

bool openWavIO();                  // Zero if operation not succesful.
                                   // Reads the header chunks of inWav and memorizes the
                                   // Headerrmation from the fmt chunk.
                                   // Creates outWav and fills its header chunks with an
                                   // identical copy of that of inWav.
                                   // Memorizes blocksize (nr of samples to read/write with subsequent calls to
                                   // readNextSampleBlock resp. writeNextSampleBlock).

bool readNextNextCodecBlock();     // False if not possible (usually because all the samples have been read).
                                   // sampleBlock: next block of samples from input wav file for all the channels.
                                   //              first sample is in element with index 0.
                                   // readsamplesCount: number of samples read (in each of the channels).

bool writeNextBTRDcodecblock();    // False if not possible.
                                   // Writes sampleBlock as the next sample block of all the channels to the output file.
                                   // The number of samples to be written is given by samplesCountToBeWritten.
                                   // First sample written is sampleBlock element with index overlapsize
                                   // - overlapsize from openWavIO call.
                                   // Last sample written is in element with index samplesCountToBeWritten - 1.
                                   // samplesCountToBeWritten should be the readsamplesCount value from the
                                   // corresponding readNextSampleBlock function call.

bool writeNextCORRcodecblock();    // False if not possible.
                                   // Writes sampleBlock as the next sample block of all the channels to the output file.
                                   // The number of samples to be written is given by samplesCountToBeWritten.
                                   // First sample written is sampleBlock element with index overlapsize
                                   // - overlapsize from openWavIO call.
                                   // Last sample written is in element with index samplesCountToBeWritten - 1.
                                   // samplesCountToBeWritten should be the readsamplesCount value from the
                                   // corresponding readNextSampleBlock function call.

bool closeWavIO();                 // False if not possible.
                                   // Finalizes the wavIO operations.

#endif // nWav_h_
