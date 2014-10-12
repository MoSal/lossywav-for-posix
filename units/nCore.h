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

#ifndef nCore_h_
#define nCore_h_
#define NOMINMAX 1

#include <windows.h> // For now...
#include <cmath>
#include <string>
#include "..\version.h"
#include "nSupport.h"
#include "nComplex.h"


#if defined(_MSC_VER) && (_MSC_VER <= 1500)
typedef   signed    char                     int8_t;
typedef unsigned    char                    uint8_t;
typedef   signed    short                   int16_t;
typedef unsigned    short                  uint16_t;
typedef   signed    long    int             int32_t;
typedef unsigned    long    int            uint32_t;
typedef   signed    long    long    int     int64_t;
typedef unsigned    long    long    int    uint64_t;
#else
#include <stdint.h>
#endif

static const char version_string[] = "lossyWAV ";

enum
{
    PRECALC_ANALYSES = 7,
    BITS_TO_CALCULATE = 32,

    WORD_SIZE = 2,
    WORD_SHIFT = 1,

    SINGLE_SIZE = 4,
    SINGLE_SHIFT = 2,

    DOUBLE_SIZE = 8,
    DOUBLE_SHIFT = 3,

    COMPLEX_SIZE = 16,
    COMPLEX_SHIFT = 4,

    INTEGER_SIZE = 4,
    INTEGER_SHIFT = 2,

    CARDINAL_SIZE = 4,
    CARDINAL_SHIFT = 2,

    INT64_SIZE = 8,
    INT64_SHIFT = 3,

    POINTER_SIZE = 4,
    POINTER_SHIFT = 2
};

//==========================================================================================================
// Many thanks to "http://keisan.casio.com/has10/Free.cgi" for the online high precision calculator.
//==========================================================================================================
static const double e              =  2.71828182845904523536028747135266249775724709369995L;    // e
static const double OneOverTwoPi   =  0.15915494309189533576888376337251436203445964574046L;    // 1/(2*Pi)
static const double OneOverPi      =  0.31830988618379067153776752674502872406891929148091L;    // 1/Pi
static const double HalfPi         =  1.57079632679489661923132169163975144209858469968760L;    // Pi/2
static const double Pi             =  3.14159265358979323846264338327950288419716939937510L;    // Pi
static const double TwoPi          =  6.28318530717958647692528676655900576839433879875020L;    // Pi*2

static const double OneOverRoot2   =  0.70710678118654752440084436210484903928483593768847L;    // 2^-0.5
static const double RootTwo        =  1.41421356237309504880168872420969807856967187537690L;    // 2^0.5
static const double log2_10        =  3.32192809488736234787031942948939017586483139302460L;    // 1/log10(2)
static const double log10_2        =  0.30102999566398119521373889472449302676818988146211L;    // log10(2)
static const double log10_2x20     =  6.02059991327962390427477789448986053536379762924220L;    // 20 * log10(2)
static const double Max_dB         =  5623413251903.49080394951039776481231468251043098690L;    // 10^(255/20)
static const double GoldenRatio    =  1.61803398874989484820458683436563811772030917980580L;    // Golden Ratio
static const double HannWindowRMS  = -4.2596873227228114834618878091836377110444078019563L;  // == RMS of Hann Window x 20 x log10(2).

//==========================================================================================================

//============================================================================
// Core settings relating to WAV limits and processing limit.
//============================================================================

static const uint64_t MAX_WAVE_SIZE = 0xFFFFFFFFuLL;

static const int32_t MAX_CHANNELS = 8;
static const int32_t MAX_BLOCK_BITS = 12;

static const int32_t PREV_CODEC_BLOCK = 0;
static const int32_t LAST_CODEC_BLOCK = 1;
static const int32_t THIS_CODEC_BLOCK = 2;
static const int32_t NEXT_CODEC_BLOCK = 3;

static const int32_t MAX_BLOCK_SIZE = 1 << MAX_BLOCK_BITS;
static const int32_t CHANNEL_BYTE_SIZE = MAX_BLOCK_SIZE * INT64_SIZE;
static const int32_t BUFFER_SIZE = MAX_CHANNELS * CHANNEL_BYTE_SIZE;

//============================================================================
// MAX_FFT_LENGTH **MUST** be double (or greater) MAX_BLOCK_SIZE (!!)
//============================================================================
static const int32_t MAX_FFT_BIT_LENGTH = MAX_BLOCK_BITS + 1;
static const int32_t MAX_FFT_LENGTH = 1 << MAX_FFT_BIT_LENGTH;
static const int32_t MAX_FFT_LENGTH_M1 = MAX_FFT_LENGTH - 1;
static const int32_t MAX_FFT_LENGTH_HALF = MAX_FFT_LENGTH / 2;

static const double MAX_FFT_LENGTH_RECIP = 1.0 / double(MAX_FFT_LENGTH);

static const int32_t LONG_FFT_LENGTH          = 1024; // for 44.1/48kHz content
static const int32_t ADDED_FFT_LENGTH_512     =  512; // for 44.1/48kHz content
static const int32_t ADDED_FFT_LENGTH_256     =  256; // for 44.1/48kHz content
static const int32_t ADDED_FFT_LENGTH_128     =  128; // for 44.1/48kHz content
static const int32_t SHORT_FFT_LENGTH         =   64; // for 44.1/48kHz content
static const int32_t IMPULSE_FFT_LENGTH       =   32; // for 44.1/48kHz content
static const int32_t ADDED_FFT_LENGTH_16      =   16; // for 44.1/48kHz content

static const int32_t LONG_FFT_BIT_LENGTH      =   10;
static const int32_t SHORT_FFT_BIT_LENGTH     =    6;

static const int32_t MAX_MUL = 8;
static const int32_t MAX_SHIFT = 3;
static const int32_t MAX_HIST_FFT_SHORT_LENGTH = SHORT_FFT_LENGTH << MAX_SHIFT;

//============================================================================
// Core Settings.
//============================================================================
static const int32_t STATIC_MINIMUM_BITS_TO_KEEP = 6;

//============================================================================
// Frequency range over which to check for lowest amplitude signal
//============================================================================
static const int32_t LOWER_FREQ_LIMIT    = 20;
static const int32_t MID_FREQ_LIMIT      = 3700;
static const int32_t UPPER_FREQ_LIMIT    = 16000;
static const int32_t HUMAN_FREQ_LIMIT    = 20000;

//============================================================================
// Types
//============================================================================

union DATA128
{
    uint64_t UInt64s[2];
    int64_t Int64s[2];
    int32_t Integers[4];
    uint32_t Cardinals[4];
    uint8_t Bytes[16];
    char Chars[16];
    float Singles[4];
    uint16_t Words[8];
    double Doubles[2];
    tDComplex Complex;
};

union DATA64
{
    uint64_t UInt64 ;
    int64_t Int64 ;
    int32_t Integers[2] ;
    uint32_t Cardinals[2] ;
    uint8_t Bytes[8] ;
    char Chars[8] ;
    float Singles[2] ;
    uint16_t Words[4] ;
    double Double ;
};

union DATA32
{
    int32_t Integer ;
    uint32_t Cardinal ;
    uint8_t Bytes[4] ;
    int8_t SmallInts[4];
    char Chars[4] ;
    float Single ;
    uint16_t Words[2] ;
    int16_t ShortInts[2];
};

union DATA16
{
    uint8_t Bytes[2];
    int8_t SmallInts[2];
    char Chars[2];
    uint16_t Word;
    int16_t ShortInt;
};

struct FFT_results_rec
{
    int16_t analysis;
    int16_t btr;
    int32_t start;
    double spreading;
};

struct FFT_Data_Rec
{
    int32_t bit_length;
    int32_t bit_shift_from_max;
    int32_t bit_shift_from_32;
    int32_t length;
    int32_t length_m1;
    int32_t length_half;
    int32_t length_half_m1;
    float length_recip;
    float length_half_recip;
    float threshold_shift;
};

typedef DATA64 SingleChannelCodecBlock[MAX_BLOCK_SIZE]                  __attribute__ ((aligned(16)));
typedef SingleChannelCodecBlock MultiChannelCodecBlock[MAX_CHANNELS]    __attribute__ ((aligned(16)));
typedef SingleChannelCodecBlock* MultiChannelCodecBlockPtr;

extern FFT_Data_Rec FFT_PreCalc_Data_Rec[MAX_FFT_BIT_LENGTH + 1]        __attribute__ ((aligned(16)));

//==============================================================================
// FFT control record - used to define and execute FFT analyses.
//==============================================================================

struct FFT_Proc_Rec
{
    int32_t NumberOfBitsNeeded;        // bit-length of array - set by programmer's code.
    int32_t BlockBitLen;            // progress counter - set and used in FFT code.

    union
    {
        void* FFT_Array;         //
        tDComplex* DComplex;     // dynamic array to contain FFT data - set by programmer's code.
        double* DReal;           //
    };

    FFT_Data_Rec* FFT;            // pointer to FFT bit-length related constants - set and used in FFT code.

    struct
    {
        int32_t block_start;
        int32_t analyses_performed;
        double (* Fill_FFT_Proc)(FFT_Proc_Rec*) = nullptr;
        double* Dist_Array;
    } Task;
};

//==============================================================================

//============================================================================
// Array for in-place FFT analysis.
//============================================================================
extern struct FFT_Array_type
{
    union
    {
        tDComplex DComplex[MAX_FFT_LENGTH];
        double DReal[MAX_FFT_LENGTH * 2];
    };
} FFT_Array  __attribute__ ((aligned(16)));


//============================================================================
// Processed FFT Result storage (FFT procedure is in-place).
//============================================================================
extern double FFT_unity_result[MAX_FFT_LENGTH]          __attribute__ ((aligned(16)));
extern double FFT_root_result[MAX_FFT_LENGTH]           __attribute__ ((aligned(16)));

extern double FFT_last_unity_result[MAX_FFT_LENGTH]     __attribute__ ((aligned(16)));
extern double FFT_last_root_result[MAX_FFT_LENGTH]      __attribute__ ((aligned(16)));

extern double FFT_skew_result[MAX_FFT_LENGTH]           __attribute__ ((aligned(16)));


//============================================================================
// All audio data
//===========   =================================================================
extern struct AudioData_type
{
    MultiChannelCodecBlock WAVEDATA[4], BTRDDATA[4], CORRDATA[4]    __attribute__ ((aligned(16)));

    double Channel_Log2_RMS[MAX_CHANNELS];
    struct
    {
        int32_t Prev, Last, This, Next;
    }
    Size;

    MultiChannelCodecBlockPtr WAVEPTR[4], BTRDPTR[4], CORRPTR[4];

    int32_t Rev_LUT[4];
} AudioData     __attribute__ ((aligned(16)));


//============================================================================
// Lookup table for integer reciprocal.
//============================================================================
extern double OneOver[MAX_FFT_LENGTH + 2]     __attribute__ ((aligned(16)));

//============================================================================
// Global variables used by the Units.
//============================================================================
extern struct timer_type
{
    union
    {
        LARGE_INTEGER Frequency;
        int64_t Frequency64;
    };
    double Period;
    union
    {
        LARGE_INTEGER Start;
        int64_t Start64;
    };
    union
    {
        LARGE_INTEGER Stop;
        int64_t Stop64;
    };
    double Elapsed;
    time_t StartTime;
} timer     __attribute__ ((aligned(16)));


extern struct Global_type
{
    struct
    {
        int32_t number;
        int32_t bits[PRECALC_ANALYSES + 1];
        int32_t length[PRECALC_ANALYSES + 1];
    }
    analysis;

    struct
    {
        int32_t aclip;
        int32_t eclip;
        int32_t noise;
        int32_t rclip;
        int32_t round;
        int32_t sclip;
    } feedback;

    struct
    {
        int32_t Size;
        double  Size_recip;
        double  duration;
        int32_t bits;
        int32_t bit_shift;
        int64_t Total;
    } Codec_Block;

    FFT_Data_Rec FFT;

    int32_t Channel;
    int32_t Channels;
    int32_t sample_rate;
    int32_t bits_per_sample;
    int32_t bytes_per_sample;
    uint64_t WAVE_size;
    uint64_t Total_Samples;
    int32_t lower_freq_limit;
    int32_t upper_freq_limit;
    uint64_t blocks_processed;

    uint64_t last_print;
    uint64_t output_blocks;
    double sample_rate_recip;
    double blocks_processed_recip;
    double processing_rate;
    uint64_t samples_processed;
} Global     __attribute__ ((aligned(16)));


extern struct parameters_type
{
    std::string wavName;
    std::string stdinname;
    std::string lossyName;
    std::string lwcdfName;
    std::string WavInpDir;
    std::string WavOutDir;
    bool forcing;
    bool correction;
    bool checking;
    bool merging;
    bool linkchannels;
    bool ignorechunksizes;
    bool STDINPUT;
    bool STDOUTPUT;
    bool skewing;

    bool midside;

    struct
    {
        bool dccorrect;
        int32_t analyses;
        int32_t underlap;
    } fft;

    bool parameters_checked;

    struct
    {
        std::string logfilename;
        bool silent;
        bool detail;
        bool verbosity;
        bool warnings;
        bool writetolog;
        bool logisunique;
        int32_t spread;
        bool blockdist;
        bool sampledist;
        bool logfileopened;
        bool perchannel;
        bool bitdist;
        bool freqdist;
        bool postanalyse;
        bool histogram;
        int32_t width;
        bool longdist;
    } output;

    bool altspread;

    int32_t noisecalc;

    struct
    {
        bool    active;
        double  numeric;
        double  round;
        double  noise;
        int32_t rclips;
        double  alevel;
        int32_t aclips;
        bool    verbose;
    } feedback;

    struct
    {
        bool    active;
        bool    warp;
        double  scale;
        bool    fixed;
        int32_t taps;
        double  extra;
        bool    altfilter;
        bool    interp_cubic;
        bool    interp_warp;
    } shaping;

    double quality;
    double scaling;
    int32_t priority;
    int32_t help;
    int32_t limit;
} parameters     __attribute__ ((aligned(16)));


extern struct settings_type
{
    double noise_threshold_shift_minimum;
    double noise_threshold_shift_minimum_alt;
    double noise_threshold_shift_average;

    double dccorrect_multiplier;

    double dynamic_minimum_bits_to_keep;
    double dynamic_maximum_bits_to_remove;

    int32_t static_minimum_bits_to_keep;
    int32_t static_maximum_bits_to_remove;

    double scaling_factor;
    double scaling_factor_inv;

    int32_t quality_integer;
    double quality_double;
    double quality_fraction;

    double fixed_noise_shaping_factor;
    int32_t FFT_analysis_switches;
    int32_t raw_result_short;
    int32_t raw_result_sgns;
    uint64_t program_expiry_is_checked;
} settings     __attribute__ ((aligned(16)));


extern struct process_type
{
    void (* Post_Analysis[PRECALC_ANALYSES + 1])(FFT_Proc_Rec*) ;

    struct
    {
        int32_t maximum_bits_to_remove;
        int32_t calc_bits_to_remove;
        int32_t bits_to_remove;
        int32_t bits_removed;
        int32_t bits_lost;

        struct
        {
            bool eclip;
            bool sclip;
            bool rclip;
            bool aclip;
            bool noise;
            bool round;
            bool retry;
        } Incidence;

        struct
        {
            int64_t eclip;
            int64_t sclip;
            int64_t rclip;
            int64_t aclip;
            int64_t noise;
            int64_t round;
        } Total;

        struct
        {
            int64_t eclips;
            int64_t sclips;
            int64_t rclips;
            int64_t aclips;
        } Count;

        FFT_results_rec min_FFT_result;

    } Channel_Data[MAX_CHANNELS]    __attribute__ ((aligned(16)));

    struct
    {
        double old_minimum;
        double new_minimum;
        double alt_average;
    } FFT_spreading[PRECALC_ANALYSES + 1][MAX_CHANNELS];

    struct
    {
        int32_t minstart, maxend;
    } limits;

    double   FFT_underlap_length[MAX_FFT_BIT_LENGTH + 1];
    int32_t  analysis_blocks[MAX_FFT_BIT_LENGTH + 1];
    int32_t  end_overlap_length[MAX_FFT_BIT_LENGTH + 1];
    int32_t  actual_analysis_blocks_start[MAX_FFT_BIT_LENGTH + 1];
    uint64_t Analyses_Completed[PRECALC_ANALYSES + 1];

    uint64_t Old_Min_Used_History[PRECALC_ANALYSES + 1][MAX_FFT_LENGTH_HALF]    __attribute__ ((aligned(16)));
    uint64_t New_Min_Used_History[PRECALC_ANALYSES + 1][MAX_FFT_LENGTH_HALF]    __attribute__ ((aligned(16)));
    uint64_t Alt_Ave_Used[PRECALC_ANALYSES + 1];
    uint64_t Old_Min_Used[PRECALC_ANALYSES + 1];
    uint64_t New_Min_Used[PRECALC_ANALYSES + 1];
    uint64_t Over_Static[PRECALC_ANALYSES + 1];
    uint64_t Over_Dynamic[PRECALC_ANALYSES + 1];
    int32_t  old_min_bin;
    int32_t  new_min_bin;
    int32_t  dynamic_maximum_bits_to_remove;
} process     __attribute__ ((aligned(16)));


extern struct results_type
{
    FFT_results_rec saved_FFT_results[MAX_CHANNELS][PRECALC_ANALYSES + 1];
    int32_t minima[MAX_CHANNELS][PRECALC_ANALYSES + 3];
    FFT_results_rec this_FFT_result;
} results     __attribute__ ((aligned(16)));


extern struct history_type
{
    double WAVE_results[MAX_CHANNELS][MAX_HIST_FFT_SHORT_LENGTH / 2 + 2]     __attribute__ ((aligned(16)));
    double BTRD_results[MAX_CHANNELS][MAX_HIST_FFT_SHORT_LENGTH / 2 + 2]     __attribute__ ((aligned(16)));
    double CORR_results[MAX_CHANNELS][MAX_HIST_FFT_SHORT_LENGTH / 2 + 2]     __attribute__ ((aligned(16)));
    double Last_Results[MAX_CHANNELS][MAX_HIST_FFT_SHORT_LENGTH / 2 + 2]     __attribute__ ((aligned(16)));

    uint64_t Histogram_DATA[1026]     __attribute__ ((aligned(16)));
    uint64_t Histogram_BTRD[1026]     __attribute__ ((aligned(16)));
    uint64_t Histogram_CORR[1026]     __attribute__ ((aligned(16)));

    double Histogram_Multiplier;
    int32_t Histogram_Length;
    int32_t Histogram_Offset;

    FFT_Data_Rec FFT;
}
history     __attribute__ ((aligned(16)));

extern struct LongDist_type
{
    double WAVE_results[MAX_CHANNELS][MAX_FFT_LENGTH_HALF + 2]  __attribute__ ((aligned(16)));
    double BTRD_results[MAX_CHANNELS][MAX_FFT_LENGTH_HALF + 2]  __attribute__ ((aligned(16)));
    double CORR_results[MAX_CHANNELS][MAX_FFT_LENGTH_HALF + 2]  __attribute__ ((aligned(16)));
    FFT_Data_Rec FFT;
} LongDist    __attribute__ ((aligned(16)));

extern struct SGNSFFT_type
{
    double FFT_results_short[MAX_FFT_LENGTH_HALF + 2]           __attribute__ ((aligned(16)));
    double FFT_results_short_root[MAX_FFT_LENGTH_HALF + 2]      __attribute__ ((aligned(16)));
    double FFT_results_short_log[MAX_FFT_LENGTH_HALF + 2]       __attribute__ ((aligned(16)));
    double FFT_results_long[MAX_FFT_LENGTH_HALF + 2]            __attribute__ ((aligned(16)));
    double FFT_results_long_root[MAX_FFT_LENGTH_HALF + 2]       __attribute__ ((aligned(16)));
    double Last_Results[MAX_CHANNELS][MAX_FFT_LENGTH_HALF + 2]  __attribute__ ((aligned(16)));
    FFT_Data_Rec FFT;
} SGNSFFT    __attribute__ ((aligned(16)));

extern struct strings_type
{
    std::string parameter;
    std::string datestamp;
    std::string Size;
    std::string time;
    std::string Elapsed;
    std::string estimate;
    std::string version;
    std::string version_short;
} strings;

extern struct version_type
{
    union
    {
        uint64_t MajorMinorReleaseBuild;
        struct { uint32_t ReleaseBuild, MajorMinor; };
        struct { uint16_t Build, Release, Minor, Major; };
    };
} version;

extern struct PowersOf_type
{
    int32_t TwoInt32[32]                     __attribute__ ((aligned(16)));
    int64_t TwoInt64[64]                 __attribute__ ((aligned(16)));
    int64_t TwoM1[64]                    __attribute__ ((aligned(16)));
    range_array< double, -1024, 1023 > Two __attribute__ ((aligned(16)));
    range_array< double, -308, 307 > Ten   __attribute__ ((aligned(16)));
} PowersOf;

extern struct Stats_type
{
    int64_t bits_removed[MAX_CHANNELS][34]    __attribute__ ((aligned(16)));
    int64_t bits_lost[MAX_CHANNELS][34]    __attribute__ ((aligned(16)));
    int64_t Skipped_Filters;

    struct
    {
        int64_t eclip;
        int64_t sclip;
        int64_t rclip;
        int64_t aclip;
        int64_t xclip;
        int64_t noise;
        int64_t round;
    } Incidence;

    struct
    {
        int64_t eclips;
        int64_t sclips;
        int64_t rclips;
        int64_t aclips;
        int64_t xclips;
    } Count;

    int64_t total_bits_removed;
    int64_t total_bits_lost;
}
Stats;


void lossyWAVError(std::string lwe_string, int32_t lwe_value);
void lossyWAVWarning(std::string lww_string);

void gettimer();

void setpriority(int32_t spv);

std::string CardinalToHex(uint32_t cth_c);
std::string WordToHex(uint16_t wth_w);
std::string ByteToHex(uint8_t bth_b);

void date_time_string_make(std::string& tsms, double tsmt);
void time_string_make(std::string& tsms, double tsmt);
void size_string_make(std::string& ssms, double ssmt);

void nCore_Init();

#endif // nCore_h_
