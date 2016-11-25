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

#if !defined(_WIN32) && defined(HAVE_SETPRIORITY)
#include <sys/resource.h> // setpriority()
#elif !defined(_WIN32)
#error Neither Windows API nor POSIX setpriority() seems to be available.
#endif

#include <cstdlib>

#include "nCore.h"

#include <string>

#include "nParameter.h"
#include "nMaths.h"
#include "nOutput.h"

int    main_argc = 0;
char** main_argv = nullptr;

namespace {

const char lossyWAV_GPL [] =
    "This program is free software: you can redistribute it and/or modify it under\n"
    "the terms of the GNU General Public License as published by the Free Software\n"
    "Foundation, either version 3 of the License, or (at your option) any later\n"
    "version.\n\n"
    "This program is distributed in the hope that it will be useful,but WITHOUT ANY\n"
    "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n"
    "PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n\n"
    "You should have received a copy of the GNU General Public License along with\n"
    "this program.  If not, see <http://www.gnu.org/licenses/>.\n";

const char lossyWAV_process_description [] =
    "Process Description:\n\n"
    "lossyWAV is a near lossless audio processor which dynamically reduces the\n"
    "bitdepth of the signal on a block-by-block basis. Bitdepth reduction adds noise\n"
    "to the processed output. The amount of permissible added noise is based on\n"
    "analysis of the signal levels in the default frequency range 20Hz to 16kHz.\n\n"
    "If signals above the upper limiting frequency are at an even lower level, they\n"
    "can be swamped by the added noise. This is usually inaudible, but the behaviour\n"
    "can be changed by specifying a different --limit (in the range 10kHz to 20kHz).\n\n"
    "For many audio signals there is little content at very high frequencies and\n"
    "forcing lossyWAV to keep the added noise level lower than the content at these\n"
    "frequencies can increase the bitrate of the losslessly compressed output\n"
    "dramatically for no perceptible benefit.\n\n"
    "The noise added by the process is shaped using an adaptive method provided by\n"
    "Sebastian Gesemann. This method, as implemented in lossyWAV, aims to use the\n"
    "signal itself as the basis of the filter used for noise shaping. Adaptive noise\n"
    "shaping is enabled by default.\n";

const char lossyWAV_standard_help [] =
    "Usage   : lossyWAV <input wav file> <options>\n"
    "\nExample : lossyWAV musicfile.wav\n"
    "\nQuality Options:\n\n"
    "-q, --quality <t>    where t is one of the following (default = standard):\n"
    "    I, insane        highest quality output, suitable for transcoding;\n"
    "    E, extreme       higher quality output, suitable for transcoding;\n"
    "    H, high          high quality output, suitable for transcoding;\n"
    "    S, standard      default quality output, considered to be transparent;\n"
    "    C, economic      intermediate quality output, likely to be transparent;\n"
    "    P, portable      good quality output for DAP use, may not be transparent;\n"
    "    X, extraportable lowest quality output, probably not transparent.\n"
    "\nStandard Options:\n\n"
    "-C, --correction     write correction file for processed WAV file; default=off.\n"
    "-f, --force          forcibly over-write output file if it exists; default=off.\n"
    "-h, --help           display help.\n"
    "-L, --longhelp       display extended help.\n"
    "-M, --merge          merge existing lossy.wav and lwcdf.wav files.\n"
    "-o, --outdir <t>     destination directory for the output file(s).\n"
    "-v, --version        display the lossyWAV version number.\n"
    "-w, --writetolog     create (or add to) lossyWAV.log in the output directory.\n";

const char lossyWAV_advanced_help [] =
    "\n"
    "Advanced Options:\n"
    "\n"
    "-                    take WAV input from STDIN.\n"
    "-c, --check          check if WAV file has already been processed; default=off.\n"
    "                     errorlevel=16 if already processed, 0 if not.\n"
    "-q, --quality <n>    quality preset (-5.0<=n<=10.0); (-5=lowest, 10=highest;\n"
    "                     default=2.5; I=10.0; E=7.5; H=5.0; S=2.5; C=0.0; P=-2.5;\n"
    "                     X=-5.0.\n"//; T=-7.5; D=-10; G=-12.5; A=-15).\n"
    "--, --stdout         write WAV output to STDOUT.\n"
    "    --stdinname <t>  pseudo filename to use when input from STDIN.\n"
    "\n"
    "Advanced Quality Options:\n"
    "\n"
    "-A, --altspread [n]  disables 'old' sperading mechanism in favour of 'new'\n"
    "                     mechanism (default spreading uses both 'old' and 'new'\n"
    "                     mechanisms). Takes an optional parameter, n, which relates\n"
    "                     to the proportion of adjacent bins taken into account when\n"
    "                     calculating spread average for a particular bin (0<=n<=1;\n"
    "                     default = 0.768544).\n"
    "-a, --analyses <n>   set number of FFT analysis lengths, (2<=n<=7; default=3,\n"
    "                     i.e. 32, 64 & 1024 samples. n = 2, remove 32 sample FFT;\n"
    "                     n > 3 add 16; n > 4, add 128; n > 5, add 256, n > 6, add\n"
    "                     512) n.b. FFT lengths stated are for 44.1/48kHz audio,\n"
    "                     higher sample rates will automatically increase all FFT\n"
    "                     lengths as required.\n"
    "-D, --dynamic <n>    select minimum_bits_to_keep_dynamic to n bits (default\n"
    "                     2.71 at -q X and 5.00 at -q I, 1.0 <= n <= 7.0.\n"
    "    --feedback [n]   enable experimental bit removal / adaptive noise shaping\n"
    "                     noise limiter. Tuning has been carried out at -q X and\n"
    "                     should have a negligible effect at -q S. Optional setting\n"
    "                     (0.0 <= n <= 10.0, default = 0.0) automatically selects\n"
    "                     the following parameters (0 = least effect, 10 = most):\n"
    "       r, round <n>  limit deviation from expected added noise due to rounding\n"
    "                     (-2.0 <= n <= 2.0, default = 0.0).\n"
    "       n, noise <n>  limit added noise due to adaptive noise shaping\n"
    "                     (-2.5 <= n <= 7.5, default = 0.0).\n"
    "       a, aclips <n> number of permissible exceedences of adaptive noise\n"
    "                     shaping level limit (0 <= n <= 64, default = 32).\n"
    "       A, alevel <n> adaptive noise shaping level limit (-2.0 <= n <= 2.5,\n"
    "                     default = 0.0).\n"
    "       V, verbose    enable more detailed feedback information in output.\n"
    "-I, --ignore-chunk-sizes.\n"
    "                     ignore 'RIFF' and 'data' chunk sizes in input.\n"
    "-l, --limit <n>      set upper frequency limit to be used in analyses to n Hz;\n"
    "                     (12500 <= n <= 20000*; default=16000).\n"
    "                     *: for 44.1/48 kHz audio. Upper limit for audio of\n"
    "                     other sampling rates is limited to sample-rate x 45.35%\n"
    "    --linkchannels   revert to original single bits-to-remove value for all\n"
    "                     channels rather than channel dependent bits-to-remove.\n"
    "    --maxclips <n>   set max. number of acceptable clips per channel per block;\n"
    "                     (0 <= n <= 16; default = 3,3,3,3,3,2,2,2,2,2,1,1,1,0,0,0).\n"
    "-m, --midside        analyse 2 channel audio for mid/side content.\n"
    "    --nodccorrect    disable DC correction of audio data prior to FFT analysis,\n"
    "                     default=on; (DC offset calculated per FFT data set).\n"
    "-n, --noskew         disable application of low frequency level reduction prior\n"
    "                     to determination of bits-to-remove.\n"
    "    --scale <n>      factor to scale audio by; (0.03125 < n <= 8.0; default=1).\n"
    "-s, --shaping        modify settings for noise shaping used in bit-removal:\n"
    "       a, altfilter  enable alternative adaptive shaping filter method.\n"
    "       A, average    set factor of shape modification above upper calculation\n"
    "                     frequency limit (0.00000 <= n <= 1.00000)\n"
    "       c, cubic      enable cubic interpolation when defining filter shape\n"
    "       e, extra      additional white noise to add during creation of filter\n"
    "       f, fixed      disable adaptive noise shaping (use fixed shaping)\n"
    "       h, hybrid     enable hybrid alternative to default adaptive noise shaping\n"
    "                     method. Uses all available calculated analyses to create\n"
    "                     the desired noise filter shape rather than only those for\n"
    "                     1.5ms and 20ms FFT analyses.\n"
    "       n, nowarp     disable warped noise shaping (use linear frequency shaping)\n"
    "       o, off        disable noise shaping altogether (use simple rounding)\n"
    "       s, scale <n>  change effectiveness of noise shaping (0 < n <= 2; default\n"
    "                     = 1.0)\n"
    "       t, taps <n>   select number of taps to use in FIR filter (8 <= n <= 256;\n"
    "                     default = 64)\n"
    "       w, warp       enable cubic interpolation when creating warped filter\n"
    "    --static <n>     set minimum-bits-to-keep-static to n bits (default=6;\n"
    "                     3<=n<=28, limited to bits-per-sample - 3).\n"
    "-U, --underlap <n>   enable underlap mode to increase number of FFT analyses\n"
    "                     performed at each FFT length, (n = 2, 4 or 8, default=2).\n"
    "\n"
    "Output Options:\n"
    "\n"
    "    --bitdist        show distrubution of bits to remove.\n"
    "    --blockdist      show distribution of lowest / highest significant bit of\n"
    "                     input codec-blocks and bit-removed codec-blocks.\n"
    "-d, --detail         enable per block per channel bits-to-remove data display.\n"
    "-F, --freqdist [all] enable frequency analysis display of input data. Use of \n"
    "                     'all' parameter displays all calculated analyses.\n"
    "-H, --histogram      show sample value histogram (input, lossy and correction).\n"
    "    --perchannel     show selected distribution data per channel.\n"
    "-p, --postanalyse    enable frequency analysis display of output and\n"
    "                     correction data in addition to input data.\n"
    "    --sampledist     show distribution of lowest / highest significant bit of\n"
    "                     input samples and bit-removed samples.\n"
    "    --spread [full]  show detailed [more detailed] results from the spreading/\n"
    "                     averaging algorithm.\n"
    "-W, --width <n>      select width of output options (79<=n<=255).\n"
    "\n"
    "System Options:\n"
    "\n"
    "-B, --below          set process priority to below normal.\n"
    "    --low            set process priority to low.\n"
    "-N, --nowarnings     suppress lossyWAV warnings.\n"
    "-Q, --quiet          significantly reduce screen output.\n"
    "-S, --silent         no screen output.\n";

const char lossyWAV_special_thanks [] =
    "\n"
    "Special thanks go to:\n"
    "\n"
    "David Robinson       for the publication of his lossyFLAC method, guidance, and\n"
    "                     the motivation to implement his method as lossyWAV.\n"
    "\n"
    "Horst Albrecht       for ABX testing, valuable support in tuning the internal\n"
    "                     presets, constructive criticism and all the feedback.\n"
    "\n"
    "Sebastian Gesemann   for the adaptive noise shaping method and the amount of\n"
    "                     help received in implementing it and also for the basis of\n"
    "                     the fixed noise shaping method.\n"
    "\n"
    "Tyge Lovset          for the C++ translation initiative.\n"
    "\n"
    "Matteo Frigo and     for libfftw3-3.dll contained in the FFTW distribution\n"
    "Steven G Johnson     (v3.2.1 or v3.2.2).\n"
    "\n"
    "Mark G Beckett       for the Delphi unit that provides an interface to the\n"
    "(Univ. of Edinburgh) relevant fftw routines in libfftw3-3.dll.\n"
    "\n"
    "Don Cross            for the Complex-FFT algorithm originally used.\n";

const int32_t num_quality_synonyms = 6;

const double quality_synonyms_vals [num_quality_synonyms + 1] = { -5.0, -2.5, 0.0, 2.5, 5.0, 7.5, 10.0 };
const char* quality_synonyms_short [num_quality_synonyms + 1] = { "X", "P", "C", "S", "H", "E", "I" };
const char* quality_synonyms_long  [num_quality_synonyms + 1] = { "extraportable", "portable", "economic", "standard", "high", "extreme", "insane" };

std::string this_quality_synonym_long;
std::string current_parameter;
std::string parmError;
int ThisParameterNumber;
//int max_open_retries = 0;

} // namespace


std::string ParamStr(int32_t index)
{
    return main_argv[index];
}

char parmchar()
{
    char result;
    result = '@';
    // remembering Delphi strings first char @ ordinate 1; C++ @ ordinate 0.
    if (current_parameter.length() == 2)
    {
        if (current_parameter[0] == '-')
        {
            result = current_parameter[1];
        }
    }
    return result;
}


char parmchar_II()
{
    char result;
    result = '@';

    if (current_parameter.length() == 1)
    {
        result = current_parameter[0];
    }

    return result;
}


void parmsError()
{
    if (((main_argc - 1) > 0) && (parameters.help == 0))
    {
        std::cerr << parmError << std::endl;
    }
    else
    {
        if (parameters.help > 1)
        {
            std::cerr << lossyWAV_process_description << std::endl;
        }

        std::cerr << lossyWAV_standard_help;

        if (parameters.help > 1)
        {
            std::cerr << lossyWAV_advanced_help;
        }

        std::cerr << lossyWAV_special_thanks;
    }
}


bool GetNextParamStr()
{
    if (ThisParameterNumber < (main_argc - 1))
    {
        ++ ThisParameterNumber;
        current_parameter = ParamStr(ThisParameterNumber);
        return true;
    }

    return false;
}


bool StringIsANumber(std::string this_str)
{
    int32_t num_minus = 0;
    int32_t num_point = 0;
    int32_t num_digit = 0;
    int32_t this_length = this_str.length();

    if (this_length == 0)
    {
        return false;
    }

    num_minus += (this_str[0] == '-');

    if (this_length == num_minus)
    {
        return false;
    }

    for (int32_t sa_i = num_minus; sa_i < this_length; sa_i++)
    {
        num_point += (this_str[sa_i] == '.');
        num_digit += ((this_str[sa_i] >='0') * (this_str[sa_i] <= '9'));
    }

    int32_t num_other = this_length - num_minus - num_point - num_digit;

    return ((num_digit>0) & (num_point<=1) & (num_minus<=1) & (num_other == 0));
}


bool NextParameterIsParameterOrEnd()
{
    std::string NextParamStr;

    if (ThisParameterNumber == (main_argc - 1))
    {
        return true;
    }
    else
    {
        NextParamStr = ParamStr(ThisParameterNumber + 1);
        int32_t this_length = NextParamStr.length();

        if (this_length==0)
        {
            return true;
        }

        if ((this_length==1) & (NextParamStr[0] == '-'))
        {
            return true;
        }

        if ((this_length==2) & (NextParamStr=="--"))
        {
            return true;
        }

        if (NextParamStr[0] == '-')
        {
            return !StringIsANumber(NextParamStr);
        }
    }
    return false;
}


void parmerror_multiple_selection()
{
    lossyWAVError(std::string("Multiple ") + parmError + " switches given.", 0x31);
}


void parmerror_mutually_incompatible_selection()
{
    lossyWAVError(std::string("Switches ") + parmError + " are incompatible.", 0x31);
}


void parmerror_no_value_given()
{
    lossyWAVError(std::string("No ") + parmError + " value given.", 0x31);
}


void parmerror_val_error()
{
    lossyWAVError(std::string("Error evaluating ") + parmError + " value given.", 0x31);
}


void check_permitted_values(int32_t cpv_var, int32_t cpv_low, int32_t cpv_high)    /* overload */
{
    if ((cpv_var < cpv_low) || (cpv_var > cpv_high))
    {
        lossyWAVError(std::string("Permitted ") + parmError + " values : " + NumToStr(cpv_low) + "<=n<=" + NumToStr(cpv_high), 0x31);
    }
}


void check_permitted_values(double cpv_var, double cpv_low, double cpv_high, int32_t Digits)    /* overload */
{
    if ((cpv_low > cpv_var) || (cpv_var > cpv_high))
    {
        lossyWAVError(std::string("Permitted ") + parmError + " values : " + NumToStr(cpv_low, Digits) + "<=n<=" + NumToStr(cpv_high, Digits), 0x31);
    }
}


bool check_parameter()
{
    if ((parmchar() == 'v') || (current_parameter == "--version"))
    {
        std::cout << version_string << strings.version_short << std::endl;
        throw (0);
    }

    if ((parmchar() == 'h') || (current_parameter == "--help"))
    {
        std::cerr << version_string << strings.version_short << lossyWAVHead1 << std::endl;
        std::cerr << lossyWAV_GPL << std::endl;
        parameters.help = 1;
        parmsError();
        throw (0);
    }

    if ((parmchar() == 'L') || (current_parameter == "--longhelp"))
    {
        std::cerr << version_string << strings.version_short << lossyWAVHead1 << std::endl;
        std::cerr << lossyWAV_GPL << std::endl;
        parameters.help = 2;
        parmsError();
        throw (0);
    }

    if ((parmchar() == 'q') || (current_parameter == "--quality"))
    {
        parmError = "quality preset";

        if (parameters.quality != -99)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        if (StringIsANumber(current_parameter))
        {
            if (parameters.quality == -99)
            {
                parameters.quality = std::atof(current_parameter.c_str());

                parameters.quality = nround10(parameters.quality, 4);

                check_permitted_values(parameters.quality, -5, 10, 4);

                if (parameters.quality == int32_t (parameters.quality))
                {
                    strings.parameter = std::string("--quality ") + NumToStr(int(parameters.quality));
                }
                else
                {
                    strings.parameter = std::string("--quality ") + NumToStr(parameters.quality, 4);
                }

                return true;
            }
        }

        for (int32_t qs_i = 0; qs_i <= num_quality_synonyms;  ++qs_i)
        {
            this_quality_synonym_long = quality_synonyms_long[qs_i];

            if ((current_parameter == quality_synonyms_short[qs_i]) || (current_parameter == this_quality_synonym_long))
            {
                parmError = "quality preset";
                parameters.quality = quality_synonyms_vals[qs_i];
                strings.parameter = std::string("--quality ") + this_quality_synonym_long;

                return true;
            }
        }

        if (parameters.quality == -99)
        {
            parmError = "quality preset";
            parmerror_no_value_given();
        }
    }

    if (current_parameter == "--maxclips")
    {
        parmError = "maximum clips";

        if (parameters.feedback.rclips
            != -1)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        if (!StringIsANumber(current_parameter))
        {
            parmerror_val_error();
        }

        parameters.feedback.rclips = std::atoi(current_parameter.c_str());

        check_permitted_values(parameters.feedback.rclips, 0, 16);

        return true;
    }

    if ((parmchar() == '-') || (current_parameter == "--stdout"))
    {
        parmError = "output to STDOUT";

        if (parameters.STDOUTPUT == true)
        {
            parmerror_multiple_selection();
        }

        parameters.STDOUTPUT = true;

        return true;
    }

    if ((parmchar() == 'I') || (current_parameter == "--ignore-chunk-sizes"))
    {
        parmError = "ignore RIFF chunk sizes";

        if (parameters.ignorechunksizes == true)
        {
            parmerror_multiple_selection();
        }

        parameters.ignorechunksizes = true;

        return true;
    }

    if ((parmchar() == 'B') || (current_parameter == "--below"))
    {
        parmError = "process priority";

        if (parameters.priority > 0)
        {
            parmerror_multiple_selection();
        }
#ifdef _WIN32
        setpriority(BELOW_NORMAL_PRIORITY_CLASS);
#elif defined(HAVE_SETPRIORITY)
        setpriority(PRIO_PROCESS, 0, 10);
#else
#error Neither Windows API nor POSIX setpriority() seems to be available.
#endif
        parameters.priority = 1;

        return true;
    }

    if (current_parameter == "--low")
    {
        parmError = "process priority";

        if (parameters.priority > 0)
        {
            parmerror_multiple_selection();
        }

        parameters.priority = 2;
#ifdef _WIN32
        setpriority(IDLE_PRIORITY_CLASS);
#elif defined(HAVE_SETPRIORITY)
        setpriority(PRIO_PROCESS, 0, 19); // Not IDLE, but avoids Linuxisms.
#else
#error Neither Windows API nor POSIX setpriority() seems to be available.
#endif

        return true;
    }

    if ((parmchar() == 'd') || (current_parameter == "--detail"))
    {
        parmError = "detail";

        if (parameters.output.detail)
        {
            parmerror_multiple_selection();
        }

        parameters.output.detail = true;

        return true;
    }

    if (current_parameter == "--perchannel")
    {
        parmError = "per-channel data output";

        if (parameters.output.perchannel)
        {
            parmerror_multiple_selection();
        }

        parameters.output.perchannel = true;

        return true;
    }

    if (current_parameter == "--bitdist")
    {
        parmError = "bits-to-remove distribution";

        if (parameters.output.bitdist)
        {
            parmerror_multiple_selection();
        }

        parameters.output.bitdist = true;

        return true;
    }

    if (current_parameter == "--spread")
    {
        parmError = "spreading results";

        if (parameters.output.spread != -1)
        {
            parmerror_multiple_selection();
        }

        if (NextParameterIsParameterOrEnd())
        {
            parameters.output.spread = 1;
        }
        else
        {
            GetNextParamStr();

            if (!(current_parameter == "full"))
            {
                parmError = "spreading results";
                parmerror_val_error();
            }
            else
            {
                parameters.output.spread = 2;
            }
        }

        return true;
    }

    if (current_parameter == "--blockdist")
    {
        parmError = "block lsb distribution";

        if (parameters.output.blockdist)
        {
            parmerror_multiple_selection();
        }

        parameters.output.blockdist = true;

        return true;
    }

    if (current_parameter == "--sampledist")
    {
        parmError = "sample lsb distribution";

        if (parameters.output.sampledist)
        {
            parmerror_multiple_selection();
        }

        parameters.output.sampledist = true;

        return true;
    }

    if ((parmchar() == 'F') || (current_parameter == "--freqdist"))
    {
        parmError = "input frequency distribution";

        if (parameters.output.freqdist)
        {
            parmerror_multiple_selection();
        }

        parameters.output.freqdist = true;

        if (!NextParameterIsParameterOrEnd())
        {
            if (!GetNextParamStr())
            {
                parmerror_no_value_given();
            }

            if (current_parameter == "all")
            {
                parameters.output.longdist = true;
            }
            else
            {
                parmerror_val_error();
            }
        }

        return true;
    }

    if ((parmchar() == 'W') || (current_parameter == "--width"))
    {
        if (parameters.output.width != -1)
        {
            parmError = "wide output";
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmError = "wide output";
            parmerror_no_value_given();
        }

        if (!StringIsANumber(current_parameter))
        {
            parmerror_val_error();
        }

        parameters.output.width = std::atoi(current_parameter.c_str());

        check_permitted_values(parameters.output.width, 79, 255);

        return true;
    }

    if ((parmchar() == 'H') || (current_parameter == "--histogram"))
    {
        if (parameters.output.histogram)
        {
            parmError = "sample value histogram";
            parmerror_multiple_selection();
        }

        parameters.output.histogram = true;

        return true;
    }

    if ((parmchar() == 'p') || (current_parameter == "--postanalyse"))
    {
        parmError = "output frequency distribution";

        if (parameters.output.postanalyse)
        {
            parmerror_multiple_selection();
        }

        if (!parameters.output.freqdist)
        {
            lossyWAVWarning("freqdist not selected before postanalyse");
            parameters.output.freqdist = true;
        }

        parameters.output.postanalyse = true;

        return true;
    }

    if  ((parmchar() == 'A') || (current_parameter == "--altspread"))
    {
        if (parameters.altspread)
        {
            parmError = "alternative spreading function";
            parmerror_multiple_selection();
        }

        parameters.altspread = true;

        if (!NextParameterIsParameterOrEnd())
        {
            GetNextParamStr();

            if (!StringIsANumber(current_parameter))
            {
                parmerror_val_error();
            }

            parameters.altspread_value = std::atof(current_parameter.c_str());

            parameters.altspread_value = nround10(parameters.altspread_value, 6);
            check_permitted_values(parameters.altspread_value, 0.000000, 1.000000, 6);
        }

        return true;
    }

    if ((parmchar() == 'Q') || (current_parameter == "--quiet"))
    {
        parmError = "reduced output";

        if (!parameters.output.verbosity)
        {
            parmerror_multiple_selection();
        }

        parameters.output.verbosity = false;

        return true;
    }

    if ((parmchar() == 'S') || (current_parameter == "--silent"))
    {

        if (parameters.output.silent)
        {
            parmError = "silent running";
            parmerror_multiple_selection();
        }

        parameters.output.silent = true;

        return true;
    }


    if ((parmchar() == 'n') || (current_parameter == "--noskew"))
    {
        if (!parameters.skewing)
        {
            parmError = "disable skewing";
            parmerror_multiple_selection();
        }

        parameters.skewing = false;

        return true;
    }


    if ((parmchar() == 'N') || (current_parameter == "--nowarnings"))
    {
        if (!parameters.output.warnings)
        {
            parmError = "no warnings";
            parmerror_multiple_selection();
        }

        parameters.output.warnings = false;

        return true;
    }

    if ((parmchar() == 'c') || (current_parameter == "--check"))
    {
        if (parameters.checking)
        {
            parmError = "check";
            parmerror_multiple_selection();
        }

        parameters.checking = true;

        return true;
    }

    if ((parmchar() == 'M') || (current_parameter == "--merge"))
    {
        if (parameters.merging)
        {
            parmError = "merge";
            parmerror_multiple_selection();
        }

        parameters.merging = true;

        return true;
    }

    if ((parmchar() == 'w') || (current_parameter == "--writetolog"))
    {
        parmError = "write log output";

        if (parameters.output.writetolog)
        {
            parmerror_multiple_selection();
        }

        parameters.output.writetolog = true;

        if (!NextParameterIsParameterOrEnd())
        {
            GetNextParamStr();

            parameters.output.logfilename = current_parameter;

            return true;
        }

        parameters.output.logfilename = "lossyWAV.log";

        return true;
    }

    if (current_parameter == "--stdinname")
    {
        parmError = "pseudo filename";

        if (parameters.stdinname != "")
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        parameters.stdinname = current_parameter;

        return true;
    }

    if ((parmchar() == 'f') || (current_parameter == "--force"))
    {
        parmError = "forced over-write";

        if (parameters.forcing)
        {
            parmerror_multiple_selection();
        }

        parameters.forcing = true;

        return true;
    }

    if ((parmchar() == 'o') || (current_parameter == "--outdir"))
    {
        parmError = "destination folder";

        if (parameters.WavOutDir != "")
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        parameters.WavOutDir = current_parameter;
#ifdef _WIN32
        if (!((parameters.WavOutDir[parameters.WavOutDir.length()] == '\\') || (parameters.WavOutDir[parameters.WavOutDir.length()] == ':')))
        {
            parameters.WavOutDir = parameters.WavOutDir + '\\';
        }
#else
        if (parameters.WavOutDir[parameters.WavOutDir.length()] != '/')
        {
            parameters.WavOutDir = parameters.WavOutDir + '/';
        }
#endif

        return true;
    }

    if ((parmchar() == 'C') || (current_parameter == "--correction"))
    {
        if (parameters.correction)
        {
            parmError = "correction file";
            parmerror_multiple_selection();
        }

        parameters.correction = true;

        return true;
    }

    if ((current_parameter == "--feedback"))
    {
        if (parameters.feedback.active)
        {
            parmError = "feedback";
            parmerror_multiple_selection();
        }

        parameters.feedback.active = true;

        if (NextParameterIsParameterOrEnd())
        {
            return true;
        }

        do
        {
            GetNextParamStr();

            if (StringIsANumber(current_parameter))
            {
                parmError = "numeric feedback";

                if (parameters.feedback.numeric != -99)
                {
                    parmerror_multiple_selection();
                }

                parameters.feedback.numeric = nround10(std::atof(current_parameter.c_str()),2);

                check_permitted_values(parameters.feedback.numeric, 0.0, 10.0, 2);

                continue;
            }
            else
            {
                if ((parmchar_II() == 'r') || (current_parameter == "round"))
                {
                    parmError = "rounding feedback";

                    if (parameters.feedback.round != -99)
                    {
                        parmerror_multiple_selection();
                    }

                    if (NextParameterIsParameterOrEnd())
                    {
                        parmerror_no_value_given();
                    }

                    GetNextParamStr();

                    if (!StringIsANumber(current_parameter))
                    {
                        parmerror_val_error();
                    }

                    parameters.feedback.round = nround10(std::atof(current_parameter.c_str()),3);

                    check_permitted_values(parameters.feedback.round, -2.000, 2.000, 3);

                    continue;
                }

                if ((parmchar_II() == 'n') || (current_parameter == "noise"))
                {
                    parmError = "noise shaping feedback";

                    if (!parameters.shaping.active)
                    {
                        parmError = "adaptive shaping off and noise shaping feedback";
                        parmerror_mutually_incompatible_selection();
                    }

                    if (parameters.feedback.noise != -99)
                    {
                        parmerror_multiple_selection();
                    }

                    if (NextParameterIsParameterOrEnd())
                    {
                        parmerror_no_value_given();
                    }

                    GetNextParamStr();

                    if (!StringIsANumber(current_parameter))
                    {
                        parmerror_val_error();
                    }

                    parameters.feedback.noise = nround10(std::atof(current_parameter.c_str()),3);

                    check_permitted_values(parameters.feedback.noise, -2.501, 7.501, 3);

                    continue;
                }

                if ((parmchar_II() == 'a') || (current_parameter == "aclips"))
                {
                    parmError = "noise shaping clipping";

                    if (!parameters.shaping.active)
                    {
                        parmError = "adaptive shaping off and noise shaping clipping";
                        parmerror_mutually_incompatible_selection();
                    }

                    if (parameters.feedback.aclips != -99)
                    {
                        parmerror_multiple_selection();
                    }

                    if (NextParameterIsParameterOrEnd())
                    {
                        parmerror_no_value_given();
                    }

                    GetNextParamStr();

                    if (!StringIsANumber(current_parameter))
                    {
                        parmerror_val_error();
                    }

                    parameters.feedback.aclips = std::atoi(current_parameter.c_str());

                    check_permitted_values(parameters.feedback.aclips, 0, 64);

                    continue;
                }

                if ((parmchar_II() == 'A') || (current_parameter == "alevel"))
                {
                    parmError = "noise shaping clipping threshold";

                    if (!parameters.shaping.active)
                    {
                        parmError = "adaptive shaping off and noise shaping clipping threshold";
                        parmerror_mutually_incompatible_selection();
                    }

                    if (parameters.feedback.alevel != -99)
                    {
                        parmerror_multiple_selection();
                    }

                    if (NextParameterIsParameterOrEnd())
                    {
                        parmerror_no_value_given();
                    }

                    GetNextParamStr();

                    if (!StringIsANumber(current_parameter))
                    {
                        parmerror_val_error();
                    }

                    parameters.feedback.alevel = nround10(std::atof(current_parameter.c_str()),3);

                    check_permitted_values(parameters.feedback.alevel, -2.000, 2.000, 3);

                    continue;
                }

                if ((parmchar_II() == 'V') || (current_parameter == "verbose"))
                {
                    parmError = "feedback verbose information";

                    if (parameters.feedback.verbose)
                    {
                        parmerror_multiple_selection();
                    }

                    parameters.feedback.verbose = true;

                    continue;
                }

                parmError = "feedback";
                parmerror_val_error();
            }
        }
        while (!NextParameterIsParameterOrEnd());

        return true;
    }


    if ((parmchar() == 's') || (current_parameter == "--shaping"))
    {
        if (NextParameterIsParameterOrEnd())
        {
            parmError = "noise shaping";
            parmerror_no_value_given();
        }

        do
        {
            GetNextParamStr();

            if ((parmchar_II() == 'h') || (current_parameter == "hybrid"))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "adaptive shaping off and hybrid noise shaping";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.altfilter)
                {
                    parmError = "altfilter and hybrid noise shaping";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.fixed)
                {
                    parmError = "hybrid noise shaping and fixed noise shaping";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.hybrid)
                {
                    parmError = "hybrid noise shaping";
                    parmerror_multiple_selection();
                }

                parameters.shaping.hybrid = true;

                continue;
            }


            if ((parmchar_II() == 'c') || (current_parameter == "cubic"))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "adaptive shaping off and cubic interpolation";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.interp_cubic)
                {
                    parmError = "cubic interpolation";
                    parmerror_multiple_selection();
                }

                parameters.shaping.interp_cubic = true;

                continue;
            }


            if ((parmchar_II() == 'w') || (current_parameter == "warp"))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "adaptive shaping off and curvilinear warp interpolation";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.interp_warp)
                {
                    parmError = "curvilinear warp interpolation";
                    parmerror_multiple_selection();
                }

                parameters.shaping.interp_warp = true;

                continue;
            }


            if ((parmchar_II() == 'a') || ((current_parameter == "altfilter")))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "altfilter and noise shaping off";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.hybrid)
                {
                    parmError = "altfilter and hybrid noise shaping";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.fixed)
                {
                    parmError = "altfilter and fixed noise shaping";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.altfilter)
                {
                    parmError = "altfilter";
                    parmerror_multiple_selection();
                }

                parameters.shaping.altfilter = true;

                continue;
            }


            if ((parmchar_II() == 'f') || (current_parameter == "fixed"))
            {
                parmError = "fixed noise shaping";

                if (parameters.shaping.fixed)
                {
                    parmerror_multiple_selection();
                }

                parameters.shaping.fixed = true;

                continue;
            }

            if ((parmchar_II() == 'n') || (current_parameter == "nowarp"))
            {
                parmError = "noise shaping off and nowarp";

                if (!parameters.shaping.active)
                {
                    parmerror_mutually_incompatible_selection();
                }

                if (!parameters.shaping.warp)
                {
                    parmerror_multiple_selection();
                }

                parameters.shaping.warp = false;

                continue;
            }


            if ((parmchar_II() == 'o') || (current_parameter == "off"))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "noise shaping off";
                    parmerror_multiple_selection();
                }

                if (parameters.shaping.altfilter)
                {
                    parmError = "noise shaping off and altfilter";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.interp_cubic)
                {
                    parmError = "noise shaping off and cubic interpolation";
                    parmerror_mutually_incompatible_selection();
                }

                if (!parameters.shaping.warp)
                {
                    parmError = "noise shaping off and nowarp";
                    parmerror_mutually_incompatible_selection();
                }

                parameters.shaping.active = false;

                continue;
            }


            if ((parmchar_II() == 's') || (current_parameter == "scale"))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "noise shaping off and shaping scale";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.scale != -99)
                {
                    parmError = "noise shaping scale";
                    parmerror_multiple_selection();
                }


                if (NextParameterIsParameterOrEnd())
                {
                    parmError = "noise shaping scale value";
                    parmerror_no_value_given();
                }

                GetNextParamStr();

                parameters.shaping.scale = std::atof(current_parameter.c_str());

                parameters.shaping.scale = nround10(parameters.shaping.scale,4);

                check_permitted_values(parameters.shaping.scale,-3,2,4);

                continue;
            }


            if ((parmchar_II() == 'A') || (current_parameter == "average"))
            {
                if (!parameters.shaping.active)
                {
                    parmError = "noise shaping off and average factor";
                    parmerror_mutually_incompatible_selection();
                }

                if (parameters.shaping.average != -99)
                {
                    parmError = "average factor";
                    parmerror_multiple_selection();
                }


                if (NextParameterIsParameterOrEnd())
                {
                    parmError = "average factor value";
                    parmerror_no_value_given();
                }

                GetNextParamStr();

                parameters.shaping.average = std::atof(current_parameter.c_str());

                parameters.shaping.average = nround10(parameters.shaping.average,6);

                check_permitted_values(parameters.shaping.average,0.0,1.0,6);

                continue;
            }


            if ((parmchar_II() == 't') || (current_parameter == "taps"))
            {
                if (parameters.shaping.taps != -1)
                {
                    parmerror_multiple_selection();
                }

                if (NextParameterIsParameterOrEnd())
                {
                    parmError = "noise shaping taps value";
                    parmerror_no_value_given();
                }

                GetNextParamStr();

                parameters.shaping.taps = std::atoi(current_parameter.c_str());

                check_permitted_values(parameters.shaping.taps, 3, 256);

                continue;
            }


            if ((parmchar_II() == 'e') || (current_parameter == "extra"))
            {
                if (parameters.shaping.extra != -1)
                {
                    parmError = "noise shaping curve modifier";
                    parmerror_multiple_selection();
                }

                if (NextParameterIsParameterOrEnd())
                {
                    parmError = "noise shaping curve modifier value";
                    parmerror_no_value_given();
                }

                GetNextParamStr();

                if (!StringIsANumber(current_parameter))
                {
                    parmerror_val_error();
                }

                parameters.shaping.extra = std::atof(current_parameter.c_str());

                check_permitted_values(parameters.shaping.extra, 0.0, 25.0,0);

                continue;
            }

            parmError = "invalid shaping option";
            parmerror_no_value_given();
        }
        while (!NextParameterIsParameterOrEnd());

        return true;
    }


    if ((parmchar() == 'D') || (current_parameter == "--dynamic"))
    {
        parmError = "minimum bits to keep (dynamic)";

        if (parameters.dynamic != -1)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        parameters.dynamic = nround10(std::atof(current_parameter.c_str()),4);

        check_permitted_values(parameters.dynamic, 1.0, 7.0, 4);

        return true;
    }

    if (current_parameter == "--static")
    {
        parmError = "minimum bits to keep (static)";

        if (parameters.Static != -1)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        parameters.Static = std::atoi(current_parameter.c_str());

        check_permitted_values(parameters.Static, 3, 28);

        return true;
    }


    if (current_parameter == "--scale")
    {
        parmError = "scaling factor";

        if (parameters.scaling != -1)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        if (!StringIsANumber(current_parameter))
        {
            parmerror_val_error();
        }

        parameters.scaling = std::atof(current_parameter.c_str());

        parameters.scaling = nround10(parameters.scaling, 6);
        check_permitted_values(parameters.scaling, 0.03125, 8.0, 6);

        return true;
    }

    if ((parmchar() == 'l') || (current_parameter == "--limit"))
    {
        parmError = "limit frequency";

        if (parameters.limit != -1)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        if (!StringIsANumber(current_parameter))
        {
            parmerror_val_error();
        }

        parameters.limit = std::atoi(current_parameter.c_str());

        check_permitted_values(parameters.limit, 12500, 999999);

        return true;
    }

    if ((parmchar() == 'a') || (current_parameter == "--analyses"))
    {
        parmError = "number of FFT analyses";

        if (parameters.fft.analyses != 3)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        if (!StringIsANumber(current_parameter))
        {
            parmerror_val_error();
        }

        parameters.fft.analyses = std::atoi(current_parameter.c_str());

        check_permitted_values(parameters.fft.analyses, 2, 7);

        return true;
    }

    if ((parmchar() == 'm') || (current_parameter == "--midside"))
    {
        parmError = "midside";

        if (parameters.midside)
        {
            parmerror_multiple_selection();
        }

        parameters.midside = true;

        return true;
    }

    if (current_parameter == "--nodccorrect")
    {
        parmError = "DC correction";

        if (!parameters.fft.dccorrect)
        {
            parmerror_multiple_selection();
        }

        parameters.fft.dccorrect = false;

        return true;
    }

    if (current_parameter == "--linkchannels")
    {
        parmError = "link channels";

        if (parameters.linkchannels)
        {
            parmerror_multiple_selection();
        }

        parameters.linkchannels = true;

        return true;
    }

    if ((parmchar() == 'U') || (current_parameter == "--underlap"))
    {
        parmError = "underlap";

        if (parameters.fft.underlap != 0)
        {
            parmerror_multiple_selection();
        }

        if (!GetNextParamStr())
        {
            parmerror_no_value_given();
        }

        if (!StringIsANumber(current_parameter))
        {
            parmerror_val_error();
        }

        parameters.fft.underlap = std::atoi(current_parameter.c_str());

        check_permitted_values(parameters.fft.underlap, 2, 8);

        if ((parameters.fft.underlap & (parameters.fft.underlap - 1)) > 0)
        {
            lossyWAVError("Permitted underlap values : 2, 4, 8", 0x31);
        }

        return true;
    }

    return false;
}


bool getParms()
{
    bool result;
    result = false;
    ThisParameterNumber = 0;
    parameters.WavInpDir = "";
    parameters.WavOutDir = "";
    parameters.wavName = "";
    parameters.stdinname = "";
    parameters.priority = 0;
    parameters.quality = -99;
    parameters.scaling = -1;
    parameters.skewing = true;

    parameters.altspread = false;
    parameters.altspread_value = -1;

    parameters.fft.analyses = 3;
    parameters.fft.underlap = 0;
    parameters.fft.dccorrect = true;

    parameters.feedback.active = false;
    parameters.feedback.numeric = -99;
    parameters.feedback.round = -99;
    parameters.feedback.noise = -99;
    parameters.feedback.aclips = -99;
    parameters.feedback.alevel = -99;
    parameters.feedback.verbose = false;

    parameters.output.logfilename = "";
    parameters.shaping.active = true;
    parameters.shaping.warp = true;
    parameters.shaping.hybrid = false;
    parameters.shaping.scale = -99;
    parameters.shaping.average = -99;
    parameters.shaping.fixed = false;
    parameters.shaping.taps = -1;
    parameters.shaping.extra = -1;
    parameters.shaping.interp_cubic = false;
    parameters.shaping.interp_warp = false;
    parameters.shaping.altfilter = false;

    parameters.feedback.rclips = -1;
    parameters.help = 0;
    parameters.limit = -1;
    parameters.Static = -1;
    parameters.dynamic = -1;

    parameters.output.perchannel = false;
    parameters.output.bitdist = false;
    parameters.output.spread = -1;
    parameters.output.sampledist = false;
    parameters.output.freqdist = false;
    parameters.output.longdist = false;
    parameters.output.histogram = false;
    parameters.output.blockdist = false;
    parameters.output.width = -1;
    parameters.output.detail = false;
    parameters.output.sampledist = false;
    parameters.output.silent = false;
    parameters.output.warnings = true;
    parameters.output.writetolog = false;
    parameters.output.verbosity = true;

    parameters.checking = false;
    parameters.correction = false;
    parameters.forcing = false;
    parameters.linkchannels = false;
    parameters.merging = false;
    parameters.midside = false;
    parameters.STDINPUT = false;
    parameters.STDOUTPUT = false;

    if (!GetNextParamStr())
    {
        return result;
    }

    if (current_parameter == "-")
    {
        parameters.STDINPUT = true;
    }
    else if (current_parameter[0] == '-')
    {
        if ((parmchar() == 'v') || (current_parameter == "--version"))
        {
            // version number should be output to stdout, NOT stderr!
            std::cout << version_string << strings.version_short << std::endl;
            throw (0);
        }
        else
        {
            std::cerr << version_string << strings.version_short << lossyWAVHead1 << std::endl;
            std::cerr << lossyWAV_GPL << std::endl;

            if ((parmchar() == 'h') || (current_parameter == "--help") || (parmchar() == 'L') || (current_parameter == "--longhelp"))
            {
                parameters.help = 1 + int32_t ((current_parameter == "--longhelp") || (parmchar() == 'L'));
                parmsError();
                throw (0);
            }
            else
            {
                lossyWAVError("No input file given.", 0x01);
            }
        }
    }
    else
    {
        parameters.wavName = current_parameter;
    }

    if (parameters.wavName != "")
    {
        size_t found = parameters.wavName.find_last_of("/\\");
        if (found != std::string::npos)
        {
            parameters.WavInpDir = parameters.wavName.substr(0, found + 1);
            parameters.wavName = parameters.wavName.substr(found + 1);
        }
/*
        if (parameters.wavName[0] == '"')
        {
            parameters.wavName = parameters.wavName.substr(1, parameters.wavName.length() - 1);
        }

        if (parameters.wavName[parameters.wavName.length() - 1] == '"')
        {
            parameters.wavName = parameters.wavName.substr(0, parameters.wavName.length() - 1);
        }

        colonpos = parameters.wavName.find(':');

        if (colonpos != std::string::npos)
        {
            parameters.WavInpDir = parameters.wavName.substr(0, colonpos + 1);
            parameters.wavName = parameters.wavName.substr(colonpos, parameters.wavName.length() - (colonpos + 1));
        }

        slashpos = parameters.wavName.find('\\');

        while (slashpos != std::string::npos)
        {
            parameters.WavInpDir += parameters.wavName.substr(0, slashpos + 1);
            parameters.wavName = parameters.wavName.substr(slashpos, parameters.wavName.length() - (slashpos + 1));
            slashpos = parameters.wavName.find('\\');
        }
*/
    }

    while (GetNextParamStr())
    {
        if (!check_parameter())
        {
            parmError = "%lossyWAV Error% : Incorrect option: " + current_parameter + "\n";
            return result;
        }
    }

    result = true;
    return result;
}


std::string WAVFilePrintName()
{
    std::string result;

    if (parameters.STDINPUT)
    {
        if (parameters.stdinname == "")
        {
            result = "from STDIN";
        }
        else
        {
            result = parameters.stdinname;
        }
    }
    else
    {
        result = parameters.wavName;
    }

    return result;
}


std::string wavInp()
{
    return parameters.WavInpDir + parameters.wavName;
}


std::string WAVOut()
{
    return parameters.WavOutDir + parameters.wavName;
}


std::string lossyInp()
{
    return parameters.WavInpDir + parameters.lossyName;
}


std::string lossyOut()
{
    return parameters.WavOutDir + parameters.lossyName;
}


std::string lwcdfInp()
{
    return parameters.WavInpDir + parameters.lwcdfName;
}


std::string lwcdfOut()
{
    return parameters.WavOutDir + parameters.lwcdfName;
}


void open_log_file()
{
    if (!parameters.output.logfileopened)
    {
        LogOutput.open(parameters.output.logfilename.c_str(), std::ios::out | std::ios::app);
        if (LogOutput.good())
        {
            parameters.output.logfileopened = true;
        }
        else
        {
            parameters.output.logfileopened = false;
            lossyWAVError("Cannot gain write access to " + parameters.output.logfilename, 0x31);
        }
    }
}


void close_log_file()
{
    if (parameters.output.logfileopened)
    {
        LogOutput.close();
        parameters.output.logfileopened = false;
    }
}


void nParameter_Init(int32_t argc, char* argv[])
{
    main_argc = argc;
    main_argv = argv;

    parameters.parameters_checked = false;
    parameters.output.logfileopened = false;

    if (!getParms())
    {
        std::cerr << version_string << strings.version_short << lossyWAVHead1 << std::endl;
        std::cerr << lossyWAV_GPL << std::endl;

        if (parmError != "")
        {
            std::cerr << parmError;
        }

        std::cerr << "Syntax: lossyWAV <wavfilename>\n";
        std::cerr << std::endl << "Type 'lossyWAV --help' for more detail.\n";

        throw (0);
    }

    if (parameters.quality == -99)
    {
        parameters.quality = 2.5;
        strings.parameter = "--quality standard";
    }

    if (parameters.help >= 1)
    {
        std::cerr << version_string << strings.version_short << lossyWAVHead1 << std::endl;
        std::cerr << lossyWAV_GPL << std::endl;
        parmsError();
        throw (0);
    }

    if ((parameters.output.verbosity) && (!parameters.output.silent))
    {
        std::cerr << version_string << strings.version_short << lossyWAVHead1 << lossyWAVHead2;
    }

    if ((!parameters.STDINPUT) && (parameters.wavName.length() == 0))
    {
        lossyWAVError("Name of input or output wav file missing.", 0x31);
    }

    if (parameters.STDINPUT)
    {
        parameters.wavName = "";
        parameters.lossyName = "lossywav.lossy.wav";
        parameters.lwcdfName = "lossywav.lwcdf.wav";

        if (parameters.merging)
        {
            lossyWAVError("Merge parameter is incompatible\n"
                          "                   with STDIN file input mode.", 0x31);
        }
    }
    else
    {
        if (!FileExists(wavInp()))
        {
            if (!FileExists(wavInp() + ".wav"))
            {
                if (!FileExists(wavInp() + ".rf64"))
                {
                    if (!FileExists(wavInp() + ".w64"))
                    {
                        lossyWAVError(std::string("Input file: ") + parameters.wavName + " does not exist.", 0x31);
                    }
                    else
                    {
                        parameters.wavName += ".w64";
                    }
                }
                else
                {
                    parameters.wavName += ".rf64";
                }
            }
            else
            {
                parameters.wavName += ".wav";
            }
        }
    }

    if (parameters.STDOUTPUT == true)
    {
        parameters.lossyName = "";
        parameters.lwcdfName = "";
        parameters.WavOutDir = "";

        if (parameters.correction)
        {
            lossyWAVError("Correct parameter is incompatible\n"
                          "                   with STDOUT file output mode.", 0x31);
        }
    }
    else
    {
        if (parameters.STDINPUT)
        {
            if (parameters.stdinname != "")
            {
                if (parameters.stdinname.length() < 5)
                {
                    if (parameters.stdinname[parameters.stdinname.length() - 1] == '.') // last char
                    {
                        parameters.stdinname = parameters.stdinname + "wav";
                    }
                    else
                    {
                        parameters.stdinname = parameters.stdinname + ".wav";
                    }
                }

                parameters.lossyName = parameters.stdinname.substr(0, parameters.stdinname.length() - 3) + "lossy." + parameters.stdinname.substr(parameters.stdinname.length() - 3, 3);
                parameters.lwcdfName = parameters.stdinname.substr(0, parameters.stdinname.length() - 3) + "lwcdf." + parameters.stdinname.substr(parameters.stdinname.length() - 3, 3);
            }
            else
            {
                parameters.lossyName = "stdin.lossy.wav";
                parameters.lwcdfName = "stdin.lwcdf.wav";
            }
        }
        else
        {
            int32_t sa_i = parameters.wavName.length();
            while ((parameters.wavName[sa_i-1] != '.') && (sa_i > 0))
                sa_i--;

            parameters.lossyName = parameters.wavName.substr(0, sa_i) + "lossy." + parameters.wavName.substr(sa_i, parameters.wavName.length() - sa_i);
            parameters.lwcdfName = parameters.wavName.substr(0, sa_i) + "lwcdf." + parameters.wavName.substr(sa_i, parameters.wavName.length() - sa_i);
        }

        if ((parameters.WavOutDir != "") && (DirectoryExists(parameters.WavOutDir) == false))
        {
            lossyWAVError("Directory : " + parameters.WavOutDir + " cannot be accessed.", 0x31);
        }
    }

    if (!parameters.merging == true)
    {
        if ((parameters.STDOUTPUT == false) && (FileExists(lossyOut()) == true))
        {
            if ((parameters.forcing == true) || (parameters.checking == true))
            {
                if (FileIsReadOnly(lossyOut()))
                {
                    lossyWAVWarning("Output lossy file is read-only.");

                    if (!SetFileReadWrite(lossyOut()))
                    {
                        lossyWAVError(std::string("Cannot gain write access to ") + lossyOut(), 0x31);
                    }
                }
            }
            else
            {
                lossyWAVError(std::string("Output file: ") + lossyOut() + " exists.", 0x31);
            }
        }

        if ((parameters.correction) && (FileExists(lwcdfOut()) == true))
        {
            if ((parameters.forcing == true) || (parameters.checking == true))
            {
                if (FileIsReadOnly(lwcdfOut()))
                {
                    lossyWAVWarning("Output correction file is read-only.");

                    if (!SetFileReadWrite(lwcdfOut()))
                    {
                        lossyWAVError(std::string("Cannot gain write access to ") + lossyOut(), 0x31);
                    }
                }
            }
            else
            {
                lossyWAVError(std::string("Output file: ") + lwcdfOut() + " exists.", 0x31);
            }
        }
    }
    else
    {
        parameters.lossyName = parameters.wavName;

        if (parameters.lossyName.substr(parameters.lossyName.length() - 10, 10) == ".lossy.wav")
        {
            parameters.lwcdfName = parameters.wavName.substr(0, parameters.wavName.length() - 10) + ".lwcdf.wav";
            parameters.wavName = parameters.wavName.substr(0, parameters.wavName.length() - 10) + ".wav";
        }

        if (parameters.lossyName.substr(parameters.lossyName.length() - 10, 10) == ".lossy.w64")
        {
            parameters.lwcdfName = parameters.wavName.substr(0, parameters.wavName.length() - 10) + ".lwcdf.w64";
            parameters.wavName = parameters.wavName.substr(0, parameters.wavName.length() - 10) + ".w64";
        }

        if (parameters.lossyName.substr(parameters.lossyName.length() - 11, 11) == ".lossy.rf64")
        {
            parameters.lwcdfName = parameters.wavName.substr(0, parameters.wavName.length() - 11) + ".lwcdf.rf64";
            parameters.wavName = parameters.wavName.substr(0, parameters.wavName.length() - 11) + ".rf64";
        }

        if (FileExists(WAVOut()) && (!parameters.forcing))
        {
            lossyWAVError("Output merged file already exists", 0x31);
        }

        if (!FileExists(lossyInp()))
        {
            lossyWAVError("Input lossy file not found.", 0x31);
        }

        if (!FileExists(lwcdfInp()))
        {
            lossyWAVError("Input correction file not found.", 0x31);
        }
    }

    if (parameters.output.logfilename == "")
    {
        parameters.output.logfilename = "lossyWAV.log";
        parameters.output.logisunique = false;
    }
    else
        parameters.output.logisunique = true;

    parameters.output.logfilename = parameters.WavOutDir + parameters.output.logfilename;
}

void nParameter_Cleanup()
{
    if (bit_removal_history != nullptr)
        delete[] bit_removal_history;
}
