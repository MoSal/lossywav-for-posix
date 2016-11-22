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

#include <sstream>
#include <iomanip> // std::hex etc.
#include <iostream>
#include <cmath>
#include <ctime>

#include "nCore.h"
#include "nMaths.h"

#if !defined(_WIN32) && defined(HAVE_STD_CHRONO_STEADY_CLOCK_NOW)

void QueryPerformanceCounter(std::chrono::steady_clock::time_point *t)
{
    *t = std::chrono::steady_clock::now();
}

double GetTimeElapsed(std::chrono::steady_clock::time_point *Stop, std::chrono::steady_clock::time_point *Start)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(*Stop - *Start).count();
}

#endif

//============================================================================
// Define global vars
//============================================================================

std::ofstream LogOutput;

int ExitCode = 0;

FFT_Data_Rec FFT_PreCalc_Data_Rec[MAX_FFT_BIT_LENGTH + 2] __attribute__ ((aligned(16)));

FFT_Array_type FFT_Array __attribute__ ((aligned (16)));

double FFT_unity_result[MAX_FFT_LENGTH]      __attribute__ ((aligned(16)));
double FFT_root_result[MAX_FFT_LENGTH]       __attribute__ ((aligned(16)));

double FFT_last_unity_result[MAX_FFT_LENGTH] __attribute__ ((aligned(16)));
double FFT_last_root_result[MAX_FFT_LENGTH]  __attribute__ ((aligned(16)));

double FFT_skew_result[MAX_FFT_LENGTH]       __attribute__ ((aligned(16)));

double OneOver[MAX_FFT_LENGTH + 2]    __attribute__ ((aligned(16)));

AudioData_type  AudioData   __attribute__ ((aligned(16)));
timer_type      timer       __attribute__ ((aligned(16)));
Current_type    Current     __attribute__ ((aligned(16)));
Global_type     Global      __attribute__ ((aligned(16)));
parameters_type parameters  __attribute__ ((aligned(16)));
settings_type   settings    __attribute__ ((aligned(16)));
process_type    process     __attribute__ ((aligned(16)));
results_type    results     __attribute__ ((aligned(16)));
history_type    history     __attribute__ ((aligned(16)));
strings_type    strings     __attribute__ ((aligned(16)));
version_type    version     __attribute__ ((aligned(16)));
PowersOf_type   PowersOf    __attribute__ ((aligned(16)));
Stats_type      Stats       __attribute__ ((aligned(16)));

void lossyWAVError(std::string lwe_string, int32_t lwe_value)
{
    if ((!parameters.output.silent) && (lwe_string != ""))
    {
        std::cerr << "%lossyWAV Error%: " << lwe_string << std::endl;
    }
    throw (lwe_value);
}


void lossyWAVWarning(std::string lww_string)
{
    if ((parameters.output.warnings) && (!parameters.output.silent))
    {
        std::cerr << "%lossyWAV Warning%: " << lww_string << std::endl;
    }
}


void GetVersionInfo()
{
    strings.version = AutoVersion::FULLVERSION_STRING;
    version.Major   = AutoVersion::MAJOR;
    version.Minor   = AutoVersion::MINOR;
    version.Release = AutoVersion::BUILD;
    version.Build   = AutoVersion::REVISION;
}


void gettimer()
{
    QueryPerformanceCounter(&timer.Stop);
#ifdef _WIN32
    timer.Elapsed = (timer.Stop64 - timer.Start64) * timer.Period;
#elif defined (HAVE_STD_CHRONO_STEADY_CLOCK_NOW)
    timer.Elapsed = GetTimeElapsed(&timer.Stop, &timer.Start);
#else
#error Neither Windows API nor std::chrono seems to be available.
#endif

    if (timer.Elapsed != 0)
        Global.processing_rate = Global.samples_processed * Global.sample_rate_recip / timer.Elapsed;
}


#ifdef _WIN32
void setpriority(int32_t spv)
{
    SetPriorityClass(GetCurrentProcess(), spv);
}
#endif


std::string CardinalToHex(uint32_t cth_c)
{
    // output an 8 char uppercase hex string
    std::ostringstream ss;
    ss << std::hex << std::uppercase << std::setfill('0') << std::setw(8) << cth_c;
    return ss.str();
}


std::string WordToHex(uint16_t wth_w)
{
    // output an 8 char uppercase hex string
    std::ostringstream ss;
    ss << std::hex << std::uppercase << std::setfill('0') << std::setw(4) << wth_w;
    return ss.str();
}


std::string ByteToHex(uint8_t bth_b)
{
    // output an 8 char uppercase hex string
    std::ostringstream ss;
    ss << std::hex << std::uppercase << std::setfill('0') << std::setw(2) << bth_b;
    return ss.str();
}


void date_time_string_make(std::string& tsms, double tsmt)
{
    tsms.resize(20);
    std::time_t t = std::time_t(tsmt);
    struct tm* tmdata = std::localtime(&t);
    std::strftime(&tsms[0], 20, "%d/%m/%Y %H:%M:%S", tmdata); // fails with GNU g++ when giving 8.
    tsms = tsms.substr(0,tsms.length()-1);
}


void time_string_make(std::string &tsms, double tsmt)
{
    double factor = 1.0;

    if (tsmt < 3600)
    {
        factor = 100.0;
    }

    if (tsmt >= 86400)
    {
        factor = OneOver[60];
    }

    double this_time = nRoundEvenInt64(tsmt * factor) / factor;

    int32_t secs = int(this_time);
    int32_t hunds = (this_time - secs) * 100;

    int32_t mins = int(secs * OneOver[60]);
    secs -= (mins * 60);

    int32_t hours = int(mins * OneOver[60]);
    mins -= (hours * 60);

    int32_t days = int(hours * OneOver[24]);
    hours -= days * 24;

    std::ostringstream ss;

    if (days > 0)
    {
        ss << std::setfill('0') << std::setw(2) << days << "d" << std::setw(2) << hours << ":" << std::setw(2) << mins;
    }
    else
    {
        if (hours > 0)
        {
            ss << std::setfill('0') << std::setw(2) << hours << ":" << std::setw(2) << mins << ":" << std::setw(2) << secs;
        }
        else
        {
            ss << std::setfill('0') << std::setw(2) << mins << ":" << std::setw(2) << secs << "." << std::setw(2) << hunds;
        }
    }
    tsms = ss.str();
}


void size_string_make(std::string& ssms, double ssmt)
{
    const char size_chars[] = "kMGTPEZ";
    int32_t ssmi, ssmj;

    if (ssmt > 99999)
    {
        ssmj = int(nlog2(ssmt) * OneOver[10]);
        ssmi = int(log10(ssmt * PowersOf.TwoX[TWO_OFFSET + -ssmj * 10]));

        if (ssmi > 2)
        {
            ssmj ++;
            ssmi = int(log10(ssmt * PowersOf.TwoX[TWO_OFFSET + -ssmj * 10]));
        }

        ssmi = std::max(0, 3 - ssmi);
        ssmt = std::floor(ssmt * PowersOf.TwoX[TWO_OFFSET + -ssmj * 10] * PowersOf.TenX[TEN_OFFSET + ssmi]) * PowersOf.TenX[TEN_OFFSET + - ssmi];
        ssms = NumToStr(ssmt, ssmi) + size_chars[ssmj-1]+"iB";
    }
    else
    {
        ssmi = 0;
        ssms = NumToStr(ssmt, ssmi) + "B";
    }
}


void nCore_Init()
{
    int32_t nt_i, nt_j;

    for (nt_i = 0; nt_i <= 3; ++nt_i)
    {
        AudioData.WAVEPTR[nt_i] = AudioData.WAVEDATA[nt_i];
        AudioData.BTRDPTR[nt_i] = AudioData.BTRDDATA[nt_i];
        AudioData.CORRPTR[nt_i] = AudioData.CORRDATA[nt_i];
        AudioData.Rev_LUT[nt_i] = nt_i;
    }

    AudioData.Size.Prev = 0;
    AudioData.Size.Last = 0;
    AudioData.Size.This = 0;
    AudioData.Size.Next = 0;

    settings.scaling_factor = 1;
    settings.scaling_factor_inv = 1;
    GetVersionInfo();

    /*with version do*/
    if (version.Build < 27)
    {
        strings.version_short = NumToStr(version.Major) + '.' +  NumToStr(version.Minor) + '.' +  NumToStr(version.Release);

        if (version.Build > 0)
        {
            strings.version_short = strings.version_short + char(version.Build + 96);
        }
    }
    else
    {
        strings.version_short = strings.version;
    }

    timer.StartTime = std::time(nullptr);
    QueryPerformanceCounter(&timer.Start);
#ifdef _WIN32
    QueryPerformanceFrequency(&timer.Frequency);
    timer.Period = 1.0 / timer.Frequency64;
#endif

    for (nt_i = 0; nt_i < int(MAX_FFT_LENGTH); ++nt_i)
    {
        FFT_unity_result[nt_i] = 0;
        FFT_root_result[nt_i] = 0;
    }

    for (nt_j = 0; nt_j < MAX_CHANNELS; nt_j ++)
    {
        for (nt_i = 0; nt_i <= 33; ++nt_i)
        {
            Stats.bits_removed[nt_j][nt_i] = 0;
            Stats.bits_lost[nt_j][nt_i] = 0;
        }

        for (nt_i = 1; nt_i <= PRECALC_ANALYSES + 2; ++nt_i)
        {
            results.minima[nt_j][nt_i] = 0;
        }
    }

    Stats.Incidence.eclip = 0;
    Stats.Incidence.sclip = 0;
    Stats.Incidence.rclip = 0;
    Stats.Incidence.aclip = 0;
    Stats.Incidence.round = 0;
    Stats.Incidence.noise = 0;

    Stats.Count.eclips = 0;
    Stats.Count.sclips = 0;
    Stats.Count.rclips = 0;
    Stats.Count.aclips = 0;

    Stats.total_bits_removed = 0;
    Stats.total_bits_lost = 0;

    for (nt_i = 0; nt_i <= 1025; ++nt_i)
    {
        history.Histogram_DATA[nt_i] = 0;
        history.Histogram_BTRD[nt_i] = 0;
        history.Histogram_CORR[nt_i] = 0;
    }

    PowersOf.TwoInt64[0] = 1;
    PowersOf.TwoM1[0] = 0;

    for (nt_i = 1; nt_i <= 62; ++nt_i)
    {
        PowersOf.TwoInt64[nt_i] = PowersOf.TwoInt64[nt_i - 1] + PowersOf.TwoInt64[nt_i - 1];
        PowersOf.TwoM1[nt_i] = PowersOf.TwoInt64[nt_i] - 1;
    }

    for (nt_i = 0; nt_i < 32; ++nt_i)
        PowersOf.TwoInt32[nt_i] = PowersOf.TwoInt64[nt_i];

    for (nt_i = -1024; nt_i <= 1023; ++nt_i)
    {
        PowersOf.TwoX[TWO_OFFSET + nt_i] = ldexp(1, nt_i);
    }

    for (nt_i = -308; nt_i <= 307; ++nt_i)
    {
         PowersOf.TenX[TEN_OFFSET + nt_i] = std::pow((long double) 10.0, (long double) nt_i);
    }

    OneOver[0] = 1;

    for (nt_i = 1; nt_i <= MAX_FFT_LENGTH + 1; ++nt_i)
    {
        OneOver[nt_i] = 1.0 / nt_i;
    }

}
