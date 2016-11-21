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

#include "nShiftBlocks.h"
#include "nCore.h"
#include "nMaths.h"

inline double Channel_RMS(int32_t this_channel)
{
    double Channel_Temp = 0.0;

    for (int32_t sc_i = 0; sc_i < AudioData.Size.This; sc_i++)
        Channel_Temp += nsqrd(AudioData.WAVEPTR[THIS_CODEC_BLOCK][this_channel][sc_i].Integers[0]);

    return nlog2(Channel_Temp * Global.Codec_Block.Size_recip) * 0.50;
}

void Shift_Codec_Blocks()
{
    int32_t sc_i = AudioData.Rev_LUT[PREV_CODEC_BLOCK];
    AudioData.Rev_LUT[PREV_CODEC_BLOCK] = AudioData.Rev_LUT[LAST_CODEC_BLOCK];
    AudioData.Rev_LUT[LAST_CODEC_BLOCK] = AudioData.Rev_LUT[THIS_CODEC_BLOCK];
    AudioData.Rev_LUT[THIS_CODEC_BLOCK] = AudioData.Rev_LUT[NEXT_CODEC_BLOCK];
    AudioData.Rev_LUT[NEXT_CODEC_BLOCK] = sc_i;

    MultiChannelCodecBlockPtr sc_p = AudioData.WAVEPTR[PREV_CODEC_BLOCK];
    AudioData.WAVEPTR[PREV_CODEC_BLOCK] = AudioData.WAVEPTR[LAST_CODEC_BLOCK];
    AudioData.WAVEPTR[LAST_CODEC_BLOCK] = AudioData.WAVEPTR[THIS_CODEC_BLOCK];
    AudioData.WAVEPTR[THIS_CODEC_BLOCK] = AudioData.WAVEPTR[NEXT_CODEC_BLOCK];
    AudioData.WAVEPTR[NEXT_CODEC_BLOCK] = sc_p;

    sc_p = AudioData.BTRDPTR[PREV_CODEC_BLOCK];
    AudioData.BTRDPTR[PREV_CODEC_BLOCK] = AudioData.BTRDPTR[LAST_CODEC_BLOCK];
    AudioData.BTRDPTR[LAST_CODEC_BLOCK] = AudioData.BTRDPTR[THIS_CODEC_BLOCK];
    AudioData.BTRDPTR[THIS_CODEC_BLOCK] = AudioData.BTRDPTR[NEXT_CODEC_BLOCK];
    AudioData.BTRDPTR[NEXT_CODEC_BLOCK] = sc_p;

    sc_p = AudioData.CORRPTR[PREV_CODEC_BLOCK];
    AudioData.CORRPTR[PREV_CODEC_BLOCK] = AudioData.CORRPTR[LAST_CODEC_BLOCK];
    AudioData.CORRPTR[LAST_CODEC_BLOCK] = AudioData.CORRPTR[THIS_CODEC_BLOCK];
    AudioData.CORRPTR[THIS_CODEC_BLOCK] = AudioData.CORRPTR[NEXT_CODEC_BLOCK];
    AudioData.CORRPTR[NEXT_CODEC_BLOCK] = sc_p;

    AudioData.Size.Prev = AudioData.Size.Last;
    AudioData.Size.Last = AudioData.Size.This;
    AudioData.Size.This = AudioData.Size.Next;
    AudioData.Size.Next = 0;

    if (AudioData.Size.This > 0)
    {
        Global.Codec_Block.Size_recip = OneOver[AudioData.Size.This];
        Global.blocks_processed++;
        Global.samples_processed += AudioData.Size.This;

        Global.blocks_processed_recip = 1.0 / Global.blocks_processed;

        for (int32_t channel = 0; channel < Global.Channels; ++channel)
            AudioData.Channel_Log2_RMS[channel] = Channel_RMS(channel);

        if (parameters.midside && (Global.Channels == 2))
        {
            for (sc_i = 0; sc_i < AudioData.Size.This; ++sc_i)
            {
                double left_sample = AudioData.WAVEPTR[THIS_CODEC_BLOCK][0][sc_i].Integers[0];
                double right_sample = AudioData.WAVEPTR[THIS_CODEC_BLOCK][1][sc_i].Integers[0];

                AudioData.WAVEPTR[THIS_CODEC_BLOCK][2][sc_i].Integers[0] = nRoundEvenInt32(0.50f * (left_sample + right_sample));
                AudioData.WAVEPTR[THIS_CODEC_BLOCK][3][sc_i].Integers[0] = nRoundEvenInt32(0.50f * (left_sample - right_sample));
            }

            AudioData.Channel_Log2_RMS[2] = Channel_RMS(2);
            AudioData.Channel_Log2_RMS[3] = Channel_RMS(3);
        }
    }
    else
    {
        for (int32_t channel = 0; channel < Global.Channels; ++channel)
            AudioData.Channel_Log2_RMS[channel] = 0;
    }
}
