//==============================================================================
//
//    lossyWAV: Added noise WAV bit reduction method by David Robinson;
//              Noise shaping coefficients by Sebastian Gesemann;
//
//    Copyright (C) 2007-2013 Nick Currie, Copyleft.
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//    Contact: lossywav <at> hotmail <dot> co <dot> uk
//
//==============================================================================
//    Initial translation to C++ from Delphi
//    by Tyge Løvset (tycho), Aug. 2012
//==============================================================================

#include "math.h"
#include "nCore.h"
#include "nSpreading.h"
#include "nParameter.h"
#include "nInitialise.h"
#include "nSGNS.h"
#include "nFillFFT.h"

static const char UnitName [] = "nInitialise";

//============================================================================
// Quality Preset Data.
//============================================================================
static const int32_t QUALITY_PRESET_MIN = -5;
static const int32_t QUALITY_PRESET_MAX = 10;

//============================================================================
// constants directly related to quality_integer settings.
//============================================================================                        //        -5,        -4,        -3,        -2,        -1,         0,         1,         2,         3,         4,         5,         6,         7,         8,         9,        10 //

const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_MINIMUM_BITS_TO_KEEP       = {{    2.710f,    2.848f,    2.986f,    3.124f,    3.262f,    3.400f,    3.470f,    3.540f,    3.610f,    3.680f,    3.750f,    4.000f,    4.250f,    4.500f,    4.750f,    5.000f}};
const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_NOISE_THRESHOLD_SHIFTS     = {{   16.667f,   14.444f,   12.222f,   10.000f,    7.778f,    5.556f,    4.444f,    3.333f,    2.222f,    1.111f,    0.000f,   -2.222f,   -4.444f,   -6.667f,   -8.889f,  -11.111f}};
const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_NOISE_THRESHOLD_SHIFTS_ALT = {{  0.89195f, -0.96843f, -2.82881f, -4.70120f, -6.58560f, -8.47000f, -9.54245f,-10.61490f,-11.66250f,-12.68525f,-13.70800f,-15.72504f,-17.74209f,-19.80792f,-21.92253f,-24.03715f}};
const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_SIGNAL_TO_NOISE_RATIOS     = {{  -19.500f,  -20.500f,  -21.500f,  -22.500f,  -23.500f,  -24.500f,  -25.000f,  -25.500f,  -26.000f,  -26.500f,  -27.000f,  -28.000f,  -29.000f,  -30.000f,  -31.000f,  -32.000f}};
const range_array<int,   QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_CLIPS_PER_CHANNEL          = {{         3,         3,         3,         3,         3,         2,         2,         2,         2,         2,         1,         1,         1,         0,         0,         0}};
const range_array<int,   QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_UPPER_CALC_FREQ_LIMIT      = {{     15159,     15332,     15504,     15676,     15848,     16021,     16107,     16193,     16279,     16365,     16451,     16538,     16624,     16710,     16796,     16882}};
const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_AUTO_SHAPING_FACTOR        = {{ 0.077903f, 0.104718f, 0.133335f, 0.162936f, 0.195383f, 0.239690f, 0.265891f, 0.304356f, 0.340499f, 0.380561f, 0.422840f, 0.517277f, 0.621117f, 0.734487f, 0.863206f, 1.000000f}};


//const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_NOISE_THRESHOLD_SHIFTS_MUL = {{  -2.000f,  -8.538f, -11.565f, -13.903f, -15.551f, -17.200f, -18.335f, -19.470f, -20.530f, -21.515f, -22.500f, -24.540f, -26.580f, -28.620f, -30.660f, -31.600f }};
//const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_NOISE_THRESHOLD_SHIFTS_ALT = {{  -5.511f,  -8.401f, -11.290f, -14.178f, -17.067f, -19.956f, -21.734f, -23.512f, -25.123f, -26.567f, -28.012f, -30.567f, -33.122f, -35.679f, -38.234f, -40.789f }};

//============================================================================
// Quality_Noise_Threshold_Delta is the correction to the calculated noise
// average when using DC offset corrected FFT input data (default).
//============================================================================
const range_array<float, QUALITY_PRESET_MIN, QUALITY_PRESET_MAX> QUALITY_NOISE_THRESHOLD_DELTA      = {{ 2.420359f, 2.296165f, 2.141537f, 1.964923f, 1.753599f, 1.424952f, 1.337800f, 1.191097f, 1.096227f, 0.923858f, 0.765617f, 0.350196f, 0.000000f, 0.000000f, 0.000000f,  0.000000f}};

//============================================================================
// Check command line parameters from nParameter;
//============================================================================
void nCheck_Switches()
{
    settings.quality_double = parameters.quality;
    settings.quality_integer = floor(settings.quality_double);
    settings.quality_fraction = settings.quality_double - settings.quality_integer;

    if (parameters.fft.analyses != 3)
    {
        strings.parameter += " --analyses ";
        strings.parameter += IntToStr(parameters.fft.analyses);
    }

    if (parameters.scaling != -1)
    {
        settings.scaling_factor = parameters.scaling;
        strings.parameter += " --scale ";
        std::string temp_STR = floattostrf(settings.scaling_factor, 6);
        strings.parameter += temp_STR;
        settings.scaling_factor = atof(temp_STR.c_str());
    }
    else
    {
        settings.scaling_factor = 1.0d;
    }

    settings.scaling_factor_inv = 1.0d / settings.scaling_factor;

    if (parameters.linkchannels || parameters.midside)
    {
        strings.parameter += " --linkchannels";
    }

    if (!parameters.shaping.active)
    {
        strings.parameter += " --shaping off";
    }
    else
    {
        std::string this_parameter = " --shaping";

       if (parameters.shaping.altfilter)
        {
            this_parameter += " altfilter";
        }

        if (parameters.shaping.interp_cubic)
        {
            this_parameter += " cubic";
        }

        if (parameters.shaping.extra != -1)
        {
            this_parameter += " extra " + floattostrf(parameters.shaping.extra,0);
            parameters.shaping.extra = npower2(parameters.shaping.extra);
        }
        else
            parameters.shaping.extra = 0;

        if (parameters.shaping.fixed)
        {
            this_parameter += " fixed";
        }

        if (!parameters.shaping.warp)
        {
            this_parameter += " nowarp";
        }

        if (parameters.shaping.scale != -99)
        {
            this_parameter += " scale " + floattostrf(parameters.shaping.scale,4);
        }

        if (parameters.shaping.taps != -1)
        {
            this_parameter += " taps " + IntToStr(parameters.shaping.taps);
        }

        if (parameters.shaping.interp_warp)
        {
            this_parameter += " warp";
        }


        if (this_parameter != " --shaping")
            strings.parameter += this_parameter;
    }


    if (parameters.altspread)
    {
        strings.parameter += " --altspread";
    }


    if (parameters.feedback.active)
    {
        strings.parameter += " --feedback";

        if (parameters.feedback.numeric != -99)
        {
            strings.parameter += " " + floattostrf(parameters.feedback.numeric,2);
        }
        else
        {
            parameters.feedback.numeric = 0;
        }

        double numeric_factor = parameters.feedback.numeric * OneOver[10];
        double numeric_factor_sqrt = sqrt(numeric_factor);

        if (parameters.feedback.round != -99)
        {
            strings.parameter += " round " + floattostrf(parameters.feedback.round,3);
        }
        else
        {
            parameters.feedback.round = numeric_factor_sqrt * (-1.0d);
        }

        if (parameters.feedback.noise != -99)
        {
            strings.parameter += " noise " + floattostrf(parameters.feedback.noise,3);
        }
        else
        {
            parameters.feedback.noise = numeric_factor_sqrt * (-2.0d) + 1.625d * (parameters.shaping.altfilter);
        }

        if (parameters.feedback.aclips != -99)
        {
            strings.parameter += " aclips " + IntToStr(parameters.feedback.aclips);
        }
        else
        {
            parameters.feedback.aclips = 32 - (numeric_factor * 20) + 8 * (parameters.shaping.altfilter);
        }

        if (parameters.feedback.alevel != -99)
        {
            strings.parameter += " alevel " + floattostrf(parameters.feedback.alevel,3);
        }
        else
        {
            parameters.feedback.alevel = numeric_factor_sqrt * (-1.025d)  + 1.375d * (parameters.shaping.altfilter);
        }
    }
    else
    {
        parameters.feedback.numeric = 0;
        parameters.feedback.round = 0.0d;
        parameters.feedback.noise = 0.0d;
        parameters.feedback.aclips = 32;
        parameters.feedback.alevel = 0.0d;
    }

    parameters.feedback.round += 2.5d;
    parameters.feedback.round *= OneOver[50];
    parameters.feedback.noise += -3.50d + 0.625d;
    parameters.feedback.alevel += -3.25d;

    if (parameters.midside)
    {
        strings.parameter += " --midside";
    }

    if (parameters.fft.dccorrect)
    {
        settings.dccorrect_multiplier = 1;
    }
    else
    {
        settings.dccorrect_multiplier = 0;
        strings.parameter += " --nodccorrect";
    }

    if (parameters.fft.underlap > 0)
    {
        strings.parameter += " --underlap ";
        strings.parameter += IntToStr(parameters.fft.underlap);
    }

    if (parameters.feedback.rclips > -1)
    {
        strings.parameter += " --maxclips ";
        strings.parameter += IntToStr(parameters.feedback.rclips);
    }
    else
    {
        parameters.feedback.rclips = QUALITY_CLIPS_PER_CHANNEL[std::min(QUALITY_PRESET_MAX, settings.quality_integer + int32_t(settings.quality_fraction != 0.0d))];
    }

    if (!parameters.skewing)
    {
        strings.parameter += " --noskew";
    }
}


//============================================================================
// Setup most of the internal arrays of values to be used in the analyses;
//============================================================================

void nInitial_Setup()
{
    double this_threshold_shift;
    double this_dccorrect_shift;

    if (settings.quality_double < -5.0d)
    {
        settings.static_minimum_bits_to_keep = STATIC_MINIMUM_BITS_TO_KEEP - 1;
    }
    else
    {
        settings.static_minimum_bits_to_keep = STATIC_MINIMUM_BITS_TO_KEEP;
    }

    if (parameters.limit != -1)
    {
        Global.upper_freq_limit = std::min(floor(Global.sample_rate * 0.453515), parameters.limit * 1.0d);
        strings.parameter += " --limit ";
        strings.parameter += floattostrf(Global.upper_freq_limit, 0);
    }
    else
    {
        Global.upper_freq_limit = std::min(floor(Global.sample_rate * 0.453515), QUALITY_UPPER_CALC_FREQ_LIMIT[settings.quality_integer] * 1.0d);
    }

    Frequency_Limits[SPREAD_ZONES + 1] = Global.upper_freq_limit;
    Global.lower_freq_limit = LOWER_FREQ_LIMIT;

    if (Global.Channels <= 2)
    {
        parameters.output.perchannel = true;
    }

    history.Histogram_Multiplier = PowersOf.Two[6 - Global.bits_per_sample];
    history.Histogram_Length = 64;
    history.Histogram_Offset = history.Histogram_Length >> 1;

    //==========================================================================
    // Make an sample-rate adjustment to noise_threshold_shift.
    //==========================================================================
    this_threshold_shift = -log10_2x20 * Global.Codec_Block.bit_shift * 0.50f;

    //==========================================================================
    // Set noise-thresholds, adjusting  for samplerate - NB: most tuning was
    // performed using 44.1kHz samples!! (second tweak)
    //==========================================================================
    if (settings.quality_integer == QUALITY_PRESET_MAX)
    {
        this_dccorrect_shift = QUALITY_NOISE_THRESHOLD_DELTA[QUALITY_PRESET_MAX];
        settings.noise_threshold_shift_minimum = QUALITY_NOISE_THRESHOLD_SHIFTS[QUALITY_PRESET_MAX];
        settings.noise_threshold_shift_minimum_alt = QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[QUALITY_PRESET_MAX];
        settings.noise_threshold_shift_average = QUALITY_SIGNAL_TO_NOISE_RATIOS[QUALITY_PRESET_MAX];
        settings.dynamic_minimum_bits_to_keep = QUALITY_MINIMUM_BITS_TO_KEEP[QUALITY_PRESET_MAX];
        settings.fixed_noise_shaping_factor = QUALITY_AUTO_SHAPING_FACTOR[QUALITY_PRESET_MAX];
    }
    else
    {
        this_dccorrect_shift = (settings.quality_fraction * (QUALITY_NOISE_THRESHOLD_DELTA[settings.quality_integer + 1] - QUALITY_NOISE_THRESHOLD_DELTA[settings.quality_integer]) + QUALITY_NOISE_THRESHOLD_DELTA[settings.quality_integer]);
        settings.noise_threshold_shift_minimum = settings.quality_fraction * (QUALITY_NOISE_THRESHOLD_SHIFTS[settings.quality_integer + 1] - QUALITY_NOISE_THRESHOLD_SHIFTS[settings.quality_integer]) + QUALITY_NOISE_THRESHOLD_SHIFTS[settings.quality_integer];
        settings.noise_threshold_shift_minimum_alt = settings.quality_fraction * (QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[settings.quality_integer + 1] - QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[settings.quality_integer]) + QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[settings.quality_integer];
        settings.noise_threshold_shift_average = settings.quality_fraction * (QUALITY_SIGNAL_TO_NOISE_RATIOS[settings.quality_integer + 1] - QUALITY_SIGNAL_TO_NOISE_RATIOS[settings.quality_integer]) + QUALITY_SIGNAL_TO_NOISE_RATIOS[settings.quality_integer];
        settings.dynamic_minimum_bits_to_keep = settings.quality_fraction * (QUALITY_MINIMUM_BITS_TO_KEEP[settings.quality_integer + 1] - QUALITY_MINIMUM_BITS_TO_KEEP[settings.quality_integer]) + QUALITY_MINIMUM_BITS_TO_KEEP[settings.quality_integer];
        settings.fixed_noise_shaping_factor = settings.quality_fraction * (QUALITY_AUTO_SHAPING_FACTOR[settings.quality_integer + 1] - QUALITY_AUTO_SHAPING_FACTOR[settings.quality_integer]) + QUALITY_AUTO_SHAPING_FACTOR[settings.quality_integer];
    }


    this_dccorrect_shift *= settings.dccorrect_multiplier;
    settings.noise_threshold_shift_average += this_dccorrect_shift + this_threshold_shift;
    settings.noise_threshold_shift_minimum += this_threshold_shift;
    settings.noise_threshold_shift_minimum_alt += this_threshold_shift;


    //==========================================================================
    // Select which FFT analyses to perform.
    //==========================================================================
    settings.FFT_analysis_switches = (SHORT_FFT_LENGTH | LONG_FFT_LENGTH);

    if (parameters.fft.analyses > 2)
    {
        settings.FFT_analysis_switches |= IMPULSE_FFT_LENGTH;
    }

    if (parameters.fft.analyses > 3)
    {
        settings.FFT_analysis_switches |= ADDED_FFT_LENGTH_16;
    }

    if (parameters.fft.analyses > 4)
    {
        settings.FFT_analysis_switches |= ADDED_FFT_LENGTH_128;
    }

    if (parameters.fft.analyses > 5)
    {
        settings.FFT_analysis_switches |= ADDED_FFT_LENGTH_256;
    }

    if (parameters.fft.analyses > 6)
    {
        settings.FFT_analysis_switches |= ADDED_FFT_LENGTH_512;
    }


    //==========================================================================
    // Shift analysis flags to reflect actual codec-block-size.
    //==========================================================================
    settings.FFT_analysis_switches = floor(settings.FFT_analysis_switches * PowersOf.Two[Global.Codec_Block.bit_shift]);

    //==========================================================================
    history.FFT = FFT_PreCalc_Data_Rec[SHORT_FFT_BIT_LENGTH + Global.Codec_Block.bit_shift];
    LongDist.FFT = FFT_PreCalc_Data_Rec[Global.Codec_Block.bits];



    for (int32_t sa_i = 1; sa_i < (PRECALC_ANALYSES + 1); ++sa_i)
    {
        Global.analysis.bits[sa_i] = history.FFT.bit_length - 3 + sa_i;
        Global.analysis.length[sa_i] = PowersOf.TwoInt64[Global.analysis.bits[sa_i]];
    }

    if (parameters.output.freqdist || parameters.shaping.active)
    {
        settings.raw_result_short = history.FFT.length;
    }
    else
    {
        settings.raw_result_short = 0;
    }

    if (parameters.shaping.active)
    {
        SGNSFFT.FFT = FFT_PreCalc_Data_Rec[LONG_FFT_BIT_LENGTH + Global.Codec_Block.bit_shift];
        settings.raw_result_sgns = SGNSFFT.FFT.length;
        nSGNS_Initialise();
    }
    else
    {
        settings.raw_result_sgns = 0;
    }

    if (Global.upper_freq_limit > Global.sample_rate * 0.475)
    {
        Global.upper_freq_limit = floor(Global.sample_rate * 0.475);
    }

    if (!parameters.output.silent)
    {
        if (parameters.output.verbosity)
        {
            std::cerr   << "Filename  : " << WAVFilePrintName() << std::endl
                        << "File Info : " << std::setprecision(2) << std::fixed << (Global.sample_rate * OneOver[1000]) << "kHz; "
                        << Global.Channels << " channel; "
                        << Global.bits_per_sample << " bit; ";

            time_string_make(strings.time, double(Global.Total_Samples) / Global.sample_rate);
            size_string_make(strings.Size, double(Global.Total_Samples) * Global.Channels * Global.bytes_per_sample);

            std::cerr << strings.time << ", " << strings.Size;

            std::cerr << std::endl << "Progress  : ";
        }
        else
        {
            std::cerr << WAVFilePrintName() << ';';
        }
    }

    for (int32_t this_channel = 0; this_channel < Global.Channels; ++this_channel)
    {
        Global.Channel = this_channel;

        for (int32_t this_analysis_number = 1; this_analysis_number <= PRECALC_ANALYSES; ++this_analysis_number)
        {
            results.saved_FFT_results[Global.Channel][this_analysis_number].start = -1;
        }
    }

    Global.output_blocks = floor(PowersOf.Two[20 - Global.Codec_Block.bits] * OneOver[Global.Channels]);
    Global.last_print = 0;

    Global.blocks_processed = 1;
    Global.blocks_processed_recip = 1.0;

    for (int32_t sa_i = 0; sa_i < Global.Channels; sa_i++)
    {
        process.Channel_Data[sa_i].calc_bits_to_remove = 0;
        process.Channel_Data[sa_i].bits_removed = 0;

        process.Channel_Data[sa_i].Total.eclip = 0;
        process.Channel_Data[sa_i].Total.sclip = 0;
        process.Channel_Data[sa_i].Total.rclip = 0;
        process.Channel_Data[sa_i].Total.aclip = 0;
        process.Channel_Data[sa_i].Total.round = 0;
        process.Channel_Data[sa_i].Total.noise = 0;

        process.Channel_Data[sa_i].Count.eclips = 0;
        process.Channel_Data[sa_i].Count.sclips = 0;
        process.Channel_Data[sa_i].Count.rclips = 0;
        process.Channel_Data[sa_i].Count.aclips = 0;
    }
}

