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
    along with this program.  If not, see <http:www.gnu.org/licenses/>.

    Contact: lossywav <at> hotmail <dot> co <dot> uk

==============================================================================
    Initial translation to C++ from Delphi
    Copyright (C) Tyge Løvset (tycho), Aug. 2012
===========================================================================**/

#include "units/nCore.h"
#include "units/fftw_interface.h"
#include "units/nFFT.h"
#include "units/nFillFFT.h"
#include "units/nInitialise.h"
#include "units/nMaths.h"
#include "units/nOutput.h"
#include "units/nProcess.h"
#include "units/nParameter.h"
#include "units/nRemoveBits.h"
#include "units/nSGNS.h"
#include "units/nShiftBlocks.h"
#include "units/nSpreading.h"
#include "units/nWav.h"

class Init
{
public:
    Init()
    {
        nCore_Init();

        nWAV_Init();

        FFTW_Initialise();

        if (!FFTW_Initialised())
        {
            nFFT_Init(MAX_FFT_BIT_LENGTH);
        }
    }

    ~Init()
    {
        nWAV_Cleanup();

        nFillFFT_Cleanup();

        nSpreading_Cleanup();

        nSGNS_Cleanup();

        nParameter_Cleanup();

        if (FFTW_Initialised())
        {
            FFTW_Cleanup();
        }
        else
        {
            nFFT_Cleanup();
        }

        nProcess_Cleanup();
    }
};

int main(int32_t argc, char* argv[])
{
    Init init;

    try
    {
        nParameter_Init(argc, argv);

        nCheck_Switches();

        if (parameters.merging)
        {
            MergeFiles();
        }
        else
        {
            if (!openWavIO())
            {
                lossyWAVError("Error initialising wavIO unit.", 0x11);
            }

            if (Global.Codec_Block.Size == 0)
            {
                lossyWAVError("Error initialising wavIO unit.", 0x11);
            }

            nInitial_Setup();

            nSpreading_Init();

            nProcess_Init();

            nFillFFT_Init();          // dependent on Codec_Block_Size.

            nRemoveBits_Init();       // bitdepth and samplerate dependent.

            nOutput_Init();

            if (!readNextNextCodecBlock())
            {
                lossyWAVError("Error reading from input file.", 0x21);
            }

            Global.blocks_processed = 0;

            //==========================================================================
            // Main processing loop.
            //==========================================================================
            while (AudioData.Size.Next > 0)
            {
                Global.last_codec_block = (AudioData.Size.Next == 0);

                Global.first_codec_block = (AudioData.Size.Last == 0);

                Shift_Codec_Blocks();

                readNextNextCodecBlock();

                Process_This_Codec_Block();

                if (!writeNextBTRDcodecblock())
                {
                    lossyWAVError("Error writing to output file.", 0x21);
                }

                if (parameters.correction)
                {
                     if (!writeNextCORRcodecblock())
                    {
                        lossyWAVError("Error writing to correction file.", 0x22);
                    }
                }
            }

            if (!closeWavIO())
            {
                lossyWAVError("Error closing wavIO unit.", 0x11);
            }

            write_cleanup();
        }
    }

    catch (int32_t ret)
    {
        switch (ret)
        {
        case 0:
            return 0;
        case 1:
            return 1;
        default:
            return -1;
        }
        return -1;
    }
}
