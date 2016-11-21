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

static const int32_t QUALITY_OFFSET = -QUALITY_PRESET_MIN;
static const int32_t QUALITY_RANGE = QUALITY_PRESET_MAX - QUALITY_PRESET_MIN;

//============================================================================
// constants directly related to quality_integer settings.
//============================================================================                        //        -5,        -4,        -3,        -2,        -1,         0,         1,         2,         3,         4,         5,         6,         7,         8,         9,        10 //
const double  QUALITY_MINIMUM_BITS_TO_KEEP       [QUALITY_RANGE + 1] = {    2.710f,    2.848f,    2.986f,    3.124f,    3.262f,    3.400f,    3.470f,    3.540f,    3.610f,    3.680f,    3.750f,    4.000f,    4.250f,    4.500f,    4.750f,    5.000f};
const double  QUALITY_NOISE_THRESHOLD_SHIFTS     [QUALITY_RANGE + 1] = {   16.667f,   14.444f,   12.222f,   10.000f,    7.778f,    5.556f,    4.444f,    3.333f,    2.222f,    1.111f,    0.000f,   -2.222f,   -4.444f,   -6.667f,   -8.889f,  -11.111f};
const double  QUALITY_NOISE_THRESHOLD_SHIFTS_ALT [QUALITY_RANGE + 1] = {  0.89195f, -0.96843f, -2.82881f, -4.70120f, -6.58560f, -8.47000f, -9.54245f,-10.61490f,-11.66250f,-12.68525f,-13.70800f,-15.72504f,-17.74209f,-19.80792f,-21.92253f,-24.03715f};
const double  QUALITY_SIGNAL_TO_NOISE_RATIOS     [QUALITY_RANGE + 1] = {  -19.500f,  -20.500f,  -21.500f,  -22.500f,  -23.500f,  -24.500f,  -25.000f,  -25.500f,  -26.000f,  -26.500f,  -27.000f,  -28.000f,  -29.000f,  -30.000f,  -31.000f,  -32.000f};
const int32_t QUALITY_CLIPS_PER_CHANNEL          [QUALITY_RANGE + 1] = {         3,         3,         3,         3,         3,         2,         2,         2,         2,         2,         1,         1,         1,         0,         0,         0};
const int32_t QUALITY_UPPER_CALC_FREQ_LIMIT      [QUALITY_RANGE + 1] = {     15159,     15332,     15504,     15676,     15848,     16021,     16107,     16193,     16279,     16365,     16451,     16538,     16624,     16710,     16796,     16882};
const double  QUALITY_AUTO_SHAPING_FACTOR        [QUALITY_RANGE + 1] = { 0.077903f, 0.104718f, 0.133335f, 0.162936f, 0.195383f, 0.239690f, 0.265891f, 0.304356f, 0.340499f, 0.380561f, 0.422840f, 0.517277f, 0.621117f, 0.734487f, 0.863206f, 1.000000f};
//============================================================================
// Quality_Noise_Threshold_Delta is the correction to the calculated noise
// average when using DC offset corrected FFT input data (default).
//============================================================================
const double  QUALITY_NOISE_THRESHOLD_DELTA      [QUALITY_RANGE + 1] = { 2.420359f, 2.296165f, 2.141537f, 1.964923f, 1.753599f, 1.424952f, 1.337800f, 1.191097f, 1.096227f, 0.923858f, 0.765617f, 0.350196f, 0.000000f, 0.000000f, 0.000000f,  0.000000f};


//============================================================================
// Check command line parameters from nParameter;
//============================================================================
void nCheck_Switches()
{
    settings.quality_double = parameters.quality;
    settings.quality_integer = int32_t(settings.quality_double);
    settings.quality_fraction = settings.quality_double - settings.quality_integer;

    if (parameters.fft.analyses != 3)
    {
        strings.parameter += " --analyses ";
        strings.parameter += NumToStr(parameters.fft.analyses);
    }

    if (parameters.scaling != -1)
    {
        settings.scaling_factor = parameters.scaling;
        strings.parameter += " --scale ";
        std::string temp_STR = NumToStr(settings.scaling_factor, 6);
        strings.parameter += temp_STR;
        settings.scaling_factor = atof(temp_STR.c_str());
    }
    else
    {
        settings.scaling_factor = 1.0;
    }

    settings.scaling_factor_inv = 1.0 / settings.scaling_factor;

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
            this_parameter += " extra " + NumToStr(parameters.shaping.extra,0);
            parameters.shaping.extra = npower2(parameters.shaping.extra);
        }
        else
            parameters.shaping.extra = 0;

        if (parameters.shaping.fixed)
        {
            this_parameter += " fixed";
        }

        if (parameters.shaping.hybrid)
        {
            this_parameter += " hybrid";
        }

        if (!parameters.shaping.warp)
        {
            this_parameter += " nowarp";
        }

        if (parameters.shaping.scale != -99)
        {
            this_parameter += " scale " + NumToStr(parameters.shaping.scale,4);
        }

        if (parameters.shaping.average != -99)
        {
            this_parameter += " average " + NumToStr(parameters.shaping.average,6);
        }

        if (parameters.shaping.taps != -1)
        {
            this_parameter += " taps " + NumToStr(parameters.shaping.taps);
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
        if (parameters.altspread_value == -1)
        {
            strings.parameter += " --altspread";
            parameters.altspread_value = 1.0;
        }
        else
        {
            strings.parameter += " --altspread " + NumToStr(parameters.altspread_value,6);
        }

    }


    if (parameters.feedback.active)
    {
        strings.parameter += " --feedback";

        if (parameters.feedback.numeric != -99)
        {
            strings.parameter += " " + NumToStr(parameters.feedback.numeric,2);
        }
        else
        {
            parameters.feedback.numeric = 0;
        }

        double numeric_factor = parameters.feedback.numeric * OneOver[10];
        double numeric_factor_sqrt = sqrt(numeric_factor);

        if (parameters.feedback.round != -99)
        {
            strings.parameter += " round " + NumToStr(parameters.feedback.round,3);
        }
        else
        {
            parameters.feedback.round = numeric_factor_sqrt * (-1.0);
        }

        if (parameters.feedback.noise != -99)
        {
            strings.parameter += " noise " + NumToStr(parameters.feedback.noise,3);
        }
        else
        {
            parameters.feedback.noise = numeric_factor_sqrt * (-2.0) + 1.625 * (parameters.shaping.altfilter);
        }

        if (parameters.feedback.aclips != -99)
        {
            strings.parameter += " aclips " + NumToStr(parameters.feedback.aclips);
        }
        else
        {
            parameters.feedback.aclips = 32 - nRoundEvenInt32(numeric_factor * 20) + 8 * (parameters.shaping.altfilter) + 32 * (parameters.shaping.hybrid);
        }

        if (parameters.feedback.alevel != -99)
        {
            strings.parameter += " alevel " + NumToStr(parameters.feedback.alevel,3);
        }
        else
        {
            parameters.feedback.alevel = numeric_factor_sqrt * (-1.025)  + 1.375 * (parameters.shaping.altfilter) + 1.50 * (parameters.shaping.hybrid);
        }
    }
    else
    {
        parameters.feedback.numeric = 0;
        parameters.feedback.round = 0.0;
        parameters.feedback.noise = 0.0;
        parameters.feedback.aclips = 32;
        parameters.feedback.alevel = 0.0;
    }

    parameters.feedback.round += 2.5;
    parameters.feedback.round *= OneOver[50];
    parameters.feedback.noise += -3.50 + 0.625;
    parameters.feedback.alevel += -3.25;

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
        strings.parameter += NumToStr(parameters.fft.underlap);
    }

    if (parameters.feedback.rclips > -1)
    {
        strings.parameter += " --maxclips ";
        strings.parameter += NumToStr(parameters.feedback.rclips);
    }
    else
    {
        parameters.feedback.rclips = QUALITY_CLIPS_PER_CHANNEL[QUALITY_OFFSET + std::min(QUALITY_PRESET_MAX, settings.quality_integer + int32_t(settings.quality_fraction != 0.0))];
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

    if (parameters.shaping.hybrid)
    {
        settings.static_minimum_bits_to_keep = _STATIC_MINIMUM_BITS_TO_KEEP - 1;
        settings.dynamic_maximum_bits_to_remove -= 1;
    }
    else
    {
        settings.static_minimum_bits_to_keep = _STATIC_MINIMUM_BITS_TO_KEEP;
    }

    if (parameters.Static > -1)
    {
        parameters.Static = std::min(parameters.Static,Global.bits_per_sample - 4);
        settings.static_minimum_bits_to_keep = parameters.Static;

        strings.parameter += " --static ";
        strings.parameter += NumToStr(parameters.Static);
    }

    if (parameters.limit != -1)
    {
        Global.upper_freq_limit = std::min(floor(Global.sample_rate * 0.453515), parameters.limit * 1.0);
        strings.parameter += " --limit ";
        strings.parameter += NumToStr(Global.upper_freq_limit, 0);
    }
    else
    {
        Global.upper_freq_limit = std::min(floor(Global.sample_rate * 0.453515), QUALITY_UPPER_CALC_FREQ_LIMIT[QUALITY_OFFSET + settings.quality_integer] * 1.0);
    }

    Frequency_Limits[SPREAD_ZONES + 1] = Global.upper_freq_limit;
    Global.lower_freq_limit = LOWER_FREQ_LIMIT;

    if (Global.Channels <= 2)
    {
        parameters.output.perchannel = true;
    }

    history.Histogram_Multiplier = PowersOf.TwoX[TWO_OFFSET + 6 - Global.bits_per_sample];
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
        this_dccorrect_shift = QUALITY_NOISE_THRESHOLD_DELTA[QUALITY_OFFSET + QUALITY_PRESET_MAX];
        settings.noise_threshold_shift_minimum = QUALITY_NOISE_THRESHOLD_SHIFTS[QUALITY_OFFSET + QUALITY_PRESET_MAX];
        settings.noise_threshold_shift_minimum_alt = QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[QUALITY_OFFSET + QUALITY_OFFSET + QUALITY_PRESET_MAX];
        settings.noise_threshold_shift_average = QUALITY_SIGNAL_TO_NOISE_RATIOS[QUALITY_OFFSET + QUALITY_PRESET_MAX];
        settings.dynamic_minimum_bits_to_keep = QUALITY_MINIMUM_BITS_TO_KEEP[QUALITY_OFFSET + QUALITY_PRESET_MAX];
        settings.fixed_noise_shaping_factor = QUALITY_AUTO_SHAPING_FACTOR[QUALITY_OFFSET + QUALITY_PRESET_MAX];
    }
    else
    {
        this_dccorrect_shift = (settings.quality_fraction * (QUALITY_NOISE_THRESHOLD_DELTA[QUALITY_OFFSET + settings.quality_integer + 1] - QUALITY_NOISE_THRESHOLD_DELTA[QUALITY_OFFSET + settings.quality_integer]) + QUALITY_NOISE_THRESHOLD_DELTA[QUALITY_OFFSET + settings.quality_integer]);
        settings.noise_threshold_shift_minimum = settings.quality_fraction * (QUALITY_NOISE_THRESHOLD_SHIFTS[QUALITY_OFFSET + settings.quality_integer + 1] - QUALITY_NOISE_THRESHOLD_SHIFTS[QUALITY_OFFSET + settings.quality_integer]) + QUALITY_NOISE_THRESHOLD_SHIFTS[QUALITY_OFFSET + settings.quality_integer];
        settings.noise_threshold_shift_minimum_alt = settings.quality_fraction * (QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[QUALITY_OFFSET + QUALITY_OFFSET + settings.quality_integer + 1] - QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[QUALITY_OFFSET + QUALITY_OFFSET + settings.quality_integer]) + QUALITY_NOISE_THRESHOLD_SHIFTS_ALT[QUALITY_OFFSET + QUALITY_OFFSET + settings.quality_integer];
        settings.noise_threshold_shift_average = settings.quality_fraction * (QUALITY_SIGNAL_TO_NOISE_RATIOS[QUALITY_OFFSET + settings.quality_integer + 1] - QUALITY_SIGNAL_TO_NOISE_RATIOS[QUALITY_OFFSET + settings.quality_integer]) + QUALITY_SIGNAL_TO_NOISE_RATIOS[QUALITY_OFFSET + settings.quality_integer];
        settings.dynamic_minimum_bits_to_keep = settings.quality_fraction * (QUALITY_MINIMUM_BITS_TO_KEEP[QUALITY_OFFSET + settings.quality_integer + 1] - QUALITY_MINIMUM_BITS_TO_KEEP[QUALITY_OFFSET + settings.quality_integer]) + QUALITY_MINIMUM_BITS_TO_KEEP[QUALITY_OFFSET + settings.quality_integer];
        settings.fixed_noise_shaping_factor = settings.quality_fraction * (QUALITY_AUTO_SHAPING_FACTOR[QUALITY_OFFSET + settings.quality_integer + 1] - QUALITY_AUTO_SHAPING_FACTOR[QUALITY_OFFSET + settings.quality_integer]) + QUALITY_AUTO_SHAPING_FACTOR[QUALITY_OFFSET + settings.quality_integer];
    }


    this_dccorrect_shift *= settings.dccorrect_multiplier;
    settings.noise_threshold_shift_average += this_dccorrect_shift + this_threshold_shift;
    settings.noise_threshold_shift_minimum += this_threshold_shift;
    settings.noise_threshold_shift_minimum_alt += this_threshold_shift;


    if (parameters.dynamic != -1)
    {
        strings.parameter += " --dynamic "+ NumToStr(parameters.dynamic, 4);
        settings.dynamic_minimum_bits_to_keep = parameters.dynamic;
    }

    for (int32_t sa_i =0; sa_i <= PRECALC_ANALYSES; sa_i++)
    {
        settings.analysis[sa_i].active = false;
        settings.analysis[sa_i].FFT = FFT_PreCalc_Data_Rec[PRECALC_ANALYSES_BITLENGTHS[sa_i] + Global.Codec_Block.bit_shift];
    }

    //==========================================================================
    // Select which FFT analyses to perform.
    //==========================================================================

    settings.analysis[SHORT_ANALYSIS].active = true;
    settings.analysis[LONG_ANALYSIS].active = true;

    if ((!parameters.shaping.hybrid) && (!parameters.altspread))
    {
        settings.analysis[2].active = (parameters.fft.analyses > 2);
        settings.analysis[1].active = (parameters.fft.analyses > 3);
        settings.analysis[4].active = (parameters.fft.analyses > 4);
        settings.analysis[5].active = (parameters.fft.analyses > 5);
        settings.analysis[6].active = (parameters.fft.analyses > 6);
    }
    else
    {
        settings.analysis[2].active = (parameters.fft.analyses > 2);
        settings.analysis[6].active = (parameters.fft.analyses > 3);
        settings.analysis[5].active = (parameters.fft.analyses > 4);
        settings.analysis[4].active = (parameters.fft.analyses > 5);
        settings.analysis[1].active = (parameters.fft.analyses > 6);
    }

    //==========================================================================

    for (int32_t sa_i = 1; sa_i < (PRECALC_ANALYSES + 1); ++sa_i)
    {
        Current.Analysis.bits[sa_i] = settings.analysis[SHORT_ANALYSIS].FFT.bit_length - 3 + sa_i;
        Current.Analysis.length[sa_i] = PowersOf.TwoInt64[Current.Analysis.bits[sa_i]];
    }

    if (parameters.shaping.active)
    {
        nSGNS_Initialise();
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
                        << "Settings  : " << strings.parameter << std::endl
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
        Current.Channel = this_channel;

        for (int32_t this_analysis_number = 1; this_analysis_number <= PRECALC_ANALYSES; ++this_analysis_number)
        {
            results.saved_FFT_results[Current.Channel][this_analysis_number].start = -1;
        }
    }

    Global.output_blocks = floor(PowersOf.TwoX[TWO_OFFSET + 20 - Global.Codec_Block.bits] * OneOver[Global.Channels]);
    Global.last_print = 0;

    Global.blocks_processed = 1;
    Global.blocks_processed_recip = 1.0;

    for (int32_t this_channel = 0; this_channel < Global.Channels; this_channel++)
    {
        process.Channel_Data[this_channel].calc_bits_to_remove = 0;
        process.Channel_Data[this_channel].bits_removed = 0;

        process.Channel_Data[this_channel].Total.eclip = 0;
        process.Channel_Data[this_channel].Total.sclip = 0;
        process.Channel_Data[this_channel].Total.rclip = 0;
        process.Channel_Data[this_channel].Total.aclip = 0;
        process.Channel_Data[this_channel].Total.round = 0;
        process.Channel_Data[this_channel].Total.noise = 0;

        process.Channel_Data[this_channel].Count.eclips = 0;
        process.Channel_Data[this_channel].Count.sclips = 0;
        process.Channel_Data[this_channel].Count.rclips = 0;
        process.Channel_Data[this_channel].Count.aclips = 0;
    }
}
