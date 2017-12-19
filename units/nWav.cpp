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

#ifdef _WIN32
    #include <fcntl.h>
#endif

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include "nWav.h"
#include "nCore.h"
#include "nParameter.h" // filemode
#include "nMaths.h"
#include "nOutput.h"

#ifndef _WIN32

#define GetLastError() errno
#define ERROR_BROKEN_PIPE EPIPE

#ifdef HAVE_NANOSLEEP
#include <ctime>

int Sleep(double interval)
{
    struct timespec tm;
    interval /= 1000.0;
    tm.tv_sec = (time_t)interval;
    tm.tv_nsec = (interval - tm.tv_sec)*1000*1000*1000;
    return nanosleep(&tm, NULL);
}
#else
#error Neither Windows API nor nanosleep() seems to be available.
#endif // ifdef HAVE_NANOSLEEP

struct GUID
{
    uint32_t Data1;
    uint16_t Data2;
    uint16_t Data3;
    uint8_t Data4[8];
};


bool GUID_compare(GUID g1, GUID g2)
{
    // true if not equal, false if equal

    bool b1 = g1.Data1 != g2.Data1;
    bool b2 = g1.Data2 != g2.Data2;
    bool b3 = g1.Data3 != g2.Data3;

    if (b1 || b2 || b3) {
        return true;
    }

    for (int idx = 0; idx < 8; idx++) {
        if (g1.Data4[idx] != g2.Data4[idx]) {
            return true;
        }
    }
    return false;
}

#endif

namespace {

static const uint64_t MAX_uint64_t = 0xFFFFFFFFFFFFFFFFull;
static const uint32_t MAX_uint32_t = 0xFFFFFFFFl;

static const uint32_t ChunkDATA_array_size = 33554432;


static const GUID GuidData[] = {
{0x66666972, 0x912E, 0x11CF, 0xA5, 0xD6, 0x28, 0xDB, 0x04, 0xC1, 0x00, 0x00},  // 'riff'
{0x7473696C, 0x912F, 0x11CF, 0xA5, 0xD6, 0x28, 0xDB, 0x04, 0xC1, 0x00, 0x00},  // 'list'
{0x65766177, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'wave'
{0x20746D66, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'fmt '
{0x74636166, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'fact'
{0x61746164, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'data'
{0x6C76656C, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'levl'
{0x6b6E756A, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'junk'
{0x74786562, 0xACF3, 0x11D3, 0x8C, 0xD1, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // 'bext'
{0xABF76256, 0x392D, 0x11D2, 0x86, 0xC7, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // MARKER
{0x925F94BC, 0x525A, 0x11D2, 0x86, 0xDC, 0x00, 0xC0, 0x4F, 0x8E, 0xDB, 0x8A},  // SUMMARYLIST
{0x00000001, 0x0000, 0x0010, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71},  // WAVEFORMATEXTENSIBLE - PCM;
{0x00000003, 0x0000, 0x0010, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71},  // WAVEFORMATEXTENSIBLE - IEEE FLOAT;
{0x00000000, 0x0000, 0x0000, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}}; // VOID.


static const int32_t GUID_riff = 0;
static const int32_t GUID_list = 1;
static const int32_t GUID_wave = 2;
static const int32_t GUID_fmt  = 3;
static const int32_t GUID_fact = 4;
static const int32_t GUID_data = 5;
static const int32_t GUID_levl = 6;
static const int32_t GUID_junk = 7;
static const int32_t GUID_bext = 8;
static const int32_t GUID_PCM  = 11;
static const int32_t GUID_IEEE = 12;
static const int32_t GUID_VOID = 13;

static const int32_t CUEChunkMaxSize = 65536;
static const int32_t FACTChunkMaxSize = 65536;
static const int32_t FMTChunkMaxSize = 18 + 22;

typedef char tChunkID[4];

struct tChunk32Header
{
    union
    {
      tChunkID ID;
      uint32_t CNum;
    };

    uint32_t Size;

    union
    {
        tChunkID RIFFType;
        uint32_t RNum;
    };
};

struct tChunkHeader
{
    union
    {
        tChunkID ID;
        uint32_t CNum;
        GUID Guid;
    };

    union
    {
        uint64_t Size;
        uint32_t Size_Lo, Size_Hi;
    };
};


struct tChunkMapRecord
{
    tChunkHeader Header;
    uint8_t* DATA;
};

struct tWAVHeader
{
    tChunkHeader Header; // "RIFF" + Size
    union
    {
        tChunkID RIFFType;   // "WAVE"

        uint32_t RNum;
        GUID RIFFGuid;
    };
};


struct tChunkSize64 // declare ChunkSize64 structure
{
    char chunkId[4]; // chunk ID (i.e. “big1” – this chunk is a big one)
    union
    {
        uint64_t ChunkSize;
        uint32_t ChunkSize_Lo, ChunkSize_Hi; // high 4 byte chunk size
    };
};


struct tDS64Header
{
    tChunkHeader Header; // "ds64" + Size
    union
    {
        uint64_t RIFFSize;
        uint32_t RIFFSize_Lo, RIFFSize_Hi;
    };

    union
    {
        uint64_t DATASize;
        uint32_t DATASize_Lo, DATASize_Hi;
    };

    union
    {
        uint64_t SampleCount;
        uint32_t SampleCount_Lo, SampleCount_Hi;
    };

    uint32_t TableLength;

    tChunkSize64 Table[0];
};


struct tFMTChunk
{
    tChunkHeader Header; // "fmt " + Size
    uint16_t wFormatTag;
    uint16_t wChannels;
    uint32_t nSamplesPerSec;
    uint32_t nAvgBytesPerSec;
    uint16_t nBlockAlign;
    uint16_t wBitsPerSample;
    uint16_t cbSize;

    struct
    {
        union
        {
            uint16_t wValidBitsPerSample;
            uint16_t wSamplesPerBlock;
            uint16_t wReserved;
        };
    } Samples;

    uint32_t dwChannelMask;

    union
    {
        GUID Guid;

        struct
        {
            uint16_t SubFormat_Word[4];
            uint8_t SubFormat_Byte[8];
        };
    };
};


struct CuePoint // declare CuePoint structure
{
    uint32_t identifier;    // unique identifier for the cue point
    uint32_t position;      // position of the cue point in the play order
    char dataChunkId[4];    // normally ‘data’
    uint32_t chunkStart;    // used for wave lists
    uint32_t blockStart;    // Start of compressed data block containing the cue point
                            // (not used for PCM)
    uint32_t sampleOffset;  // sample offset of cue point (absolute for PCM,
                            // relative to block start for compressed data)
};


struct CueChunk // declare CueChunk structure
{
    tChunkHeader Header; // "cue " + Size
    uint32_t cuePointCount; // number of cue points (markers)
    CuePoint cuePoints[]; // cue points
};


struct JunkChunk // declare JunkChunk structure
{
    tChunkHeader Header; // "JUNK" + Size
    // least 28 if the chunk is intended as a place-holder for a ‘ds64’ chunk.
    char ChunkDATA[]; // dummy bytes
};


struct ListChunk // declare ListChunk structure
{
    tChunkHeader Header; // "list" + Size
    char typeId[4]; // ‘adtl’ associated data list
};


struct LabelChunk // declare LabelChunk structure
{
    tChunkHeader Header; // "labl" + Size
    uint32_t identifier; // unique identifier for the cue point
    char text[]; // label text: nullptr terminated string (ANSI)
};


struct MarkerEntry // declare MarkerEntry structure
{
    uint32_t flags; // flags field
    uint32_t sampleOffsetLow; // low 4 byte marker’s offset in samples in data chunk
    uint32_t sampleOffsetHigh; // high 4 byte marker’s offset
    uint32_t byteOffsetLow; // low and high 4 byte of the beginning of the nearest
    uint32_t byteOffsetHigh; // compressed frame next to marker (timely before)
    uint32_t intraSmplOffsetHigh; // low and high 4 byte of marker’s offset in samples
    uint32_t intraSmplOffsetLow; // relative to the position of the first sample in frame
    char     labelText[256]; // nullptr terminated label string
    uint32_t lablChunkIdentifier; // link to ‘labl’ subchunk of ‘list’ chunk8
    GUID     vendorAndProduct; // GUID identifying specific vendor application
    uint32_t userData1; // 4 byte application specific user data
    uint32_t userData2; // 4 byte application specific user data
    uint32_t userData3; // 4 byte application specific user data
    uint32_t userData4; // 4 byte application specific user data
};


struct MarkerChunk // declare MarkerChunk structure
{
    tChunkHeader Header; // "r64m" + Size
    MarkerEntry markers[]; // marker entries
};


struct tFACTChunk
{
    tChunkHeader Header; // "fact" + Size
    char DATA[65536];
};


struct tDATAChunk
{
    tChunkHeader Header; // "data" + Size
};


struct tWAVEChunks
{
    tWAVHeader WAV;
    tDS64Header DS64;
    tFMTChunk FMT;
    tFACTChunk FACT;
    tDATAChunk DATA;

    struct
    {
        struct
        {
            uint32_t Write;
            uint32_t Read;
            uint32_t Free;
        } Current;

        struct
        {
            uint32_t RIFF;
            uint32_t RF64;
            uint32_t riff;
            uint32_t ds64;
            uint32_t fmt;
            uint32_t fact;
            uint32_t DATA;
            uint32_t LAST;
        } Pos;

        struct
        {
            uint32_t Reduction;
            uint32_t Header;
            uint32_t Type;
            uint32_t Size;
            uint32_t RIFF;
        } Size;

        tChunkMapRecord Map[384];
    };
};


struct CUE64_HEADER
{
    GUID Guid;
    uint64_t llSize;
    uint32_t dwMarkerCount;
};


struct tWAVEBuffer
{
    union
    {
         uint8_t Bytes[BUFFER_SIZE];
          int8_t SmallInts[BUFFER_SIZE];
         int16_t ShortInts[BUFFER_SIZE >> 1];
        uint16_t Words[BUFFER_SIZE >> 1];
         int32_t Integers[BUFFER_SIZE >> 2];
        uint32_t Cardinals[BUFFER_SIZE >> 2];
         int64_t Int64s[BUFFER_SIZE >> 3];
    };
};

struct tRIFF_Rec
{
    tWAVEBuffer Buffer;

    tWAVEChunks Chunks;

    uint64_t BytesInBuffer;

    void (* ReadTransfer)(unsigned char* pB, unsigned char* pEndB);
    void (* WriteTransfer)(MultiChannelCodecBlock& outputcodecblock, unsigned char* pB);

    uint64_t samplebytesread;
    uint64_t firstsamplesread;
    uint64_t sampleByteLeftToRead;
    uint64_t samplesLeftToRead;

    int32_t  wBytesPerSample;

    std::fstream RIFF_File;

    int32_t ID = -1;

    struct
    {
        std::string Name;

        uint64_t     (*Read)(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestoread);

        bool         (*Write)(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestowrite);

        struct
        {
            bool Open;
            bool File;
            bool Pipe;
        } Is;

        struct
        {
            bool RIFF64;
            bool WAVE64;
        } Type;

        struct
        {
            uint64_t Read;
            uint64_t Write;
        } Last;

        struct
        {
            uint64_t DS64;
            uint64_t FMT;
            uint64_t FACT;
            uint64_t DATA;
            uint64_t CURRENT;
        } Pos;

        uint32_t Retries;

        struct
        {
            bool Read;
            bool Write;
        } Cant;

        struct
        {
            uint64_t Read;
            uint64_t Written;
        } Total_Bytes;
    } File;
};

struct RIFF
{
    tRIFF_Rec WAVE;
    tRIFF_Rec BTRD;
    tRIFF_Rec CORR;
    uint32_t  EndOfChunkMap = 0;
    uint8_t*  ChunkDATA;
} RIFF;

std::string DateTimeString;

uint64_t BUFFER_SIZEread;                    // number of byte that make up for inbuff
uint64_t BUFFER_SIZEwrite;                    // number of byte that make up for outbuff
uint64_t BUFFER_SIZEtouse;

uint64_t nrOfByteInOneBlockInBuff;            // number of byte that make up a complete block in inBuff and outBuff.
uint64_t nrOfBlocksInReadBuff;
uint64_t nrOfBlocksInWriteBuff;

uint64_t nrOfFullBlocksReadIntoInBuff;        // with the last BlockRead: how many full blocks managed it into inBuff ?
uint64_t nrOfByteOfTheLastIncompleteBlock;    // when EOF was encountered with the last BlockRead: the number of byte that didn't make it
                                              // into the last full block
uint64_t nrOfBlockInBuffNotYetFetched;        // the next number of the block in inBuff that was not fetched yet by readNextSampleBlock.
                                              // 0: no data in inBuff;

} // namespace

uint64_t readfrom_stdin(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestoread)
{
    uint64_t bytesthisread, bytesread, byteslefttoread;
    uint64_t lastwinerror;

    bytesread = 0;
    thisRIFF.File.Retries = 0;

    while ((bytesread < bytestoread) && (!thisRIFF.File.Cant.Read))
    {
        byteslefttoread = bytestoread - bytesread;

        std::cin.read(static_cast<char *>(buffpointer) + bytesread, byteslefttoread);

        bytesthisread = std::cin.gcount();

        if (std::cin.good())
        {
            lastwinerror = 0;
        }
        else
        {
            lastwinerror = GetLastError();
        }

        if (bytesthisread == 0)
        {
            ++ thisRIFF.File.Retries;
            Sleep(PowersOf.TwoX[TWO_OFFSET + std::min(7u, thisRIFF.File.Retries) -11]*1000);
        }
        else
        {
            thisRIFF.File.Retries = 0;
        }

        if ((lastwinerror == ERROR_BROKEN_PIPE) || (thisRIFF.File.Retries == 32))
        {
            thisRIFF.File.Cant.Read = true;
        }

        bytesread += bytesthisread;
    }

    thisRIFF.File.Last.Read=bytesread;

    thisRIFF.File.Total_Bytes.Read+=bytesread;

    return bytesread;
}


uint64_t readfrom_file(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestoread)
{
    uint64_t bytesread = 0;
    thisRIFF.File.Retries = 0;
    thisRIFF.RIFF_File.read((char *) buffpointer, bytestoread);
    bytesread = thisRIFF.RIFF_File.gcount();

    if ((thisRIFF.RIFF_File.bad()) || (bytesread < bytestoread))
    {
        thisRIFF.File.Cant.Read = true;
    }

    thisRIFF.File.Last.Read=bytesread;

    thisRIFF.File.Total_Bytes.Read+=bytesread;

    return bytesread;
}


uint64_t readfrom_nullptr(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestoread)
{
    thisRIFF.File.Retries = 99;
    thisRIFF.File.Cant.Read = true;
    thisRIFF.File.Last.Read = 0;
    thisRIFF.File.Total_Bytes.Read = 0;
    buffpointer = buffpointer;
    bytestoread = bytestoread;

    return 0;
}


bool writeto_stdout(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestowrite)
{
    std::cout.write(static_cast<char *>(buffpointer), bytestowrite);
    if (std::cout.bad())
    {
        thisRIFF.File.Cant.Write = true;
        return false;
    }
    thisRIFF.File.Total_Bytes.Written+=bytestowrite;
    thisRIFF.File.Last.Write=bytestowrite;

    return true;
}


bool writeto_file(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestowrite)
{
    thisRIFF.RIFF_File.write((char *) buffpointer, bytestowrite);

    if (thisRIFF.RIFF_File.bad())
    {
        thisRIFF.File.Cant.Write = false;
        return false;
    }

    thisRIFF.File.Total_Bytes.Written+=bytestowrite;
    thisRIFF.File.Last.Write=bytestowrite;

    return true;
}


bool writeto_nullptr(tRIFF_Rec &thisRIFF, void* buffpointer, uint64_t bytestowrite)
{
    thisRIFF.File.Last.Write = 0;
    thisRIFF.File.Total_Bytes.Written = 0;
    thisRIFF.File.Cant.Write = true;
    thisRIFF.File.Retries = 99;
    buffpointer = buffpointer;
    bytestowrite = bytestowrite;
    return false;
}


uint32_t WordAlign(uint32_t WAvalue)
{
    return (WAvalue + 1) & 0xFFFFFFFE;
}


uint64_t QWordAlign(uint64_t QWAvalue)
{
    return (QWAvalue + 7) & 0xFFFFFFFFFFFFFFF8;
}


bool WritePaddingToFile(tRIFF_Rec &thisRIFF)
{
    uint64_t paddingbytes = 0;

    if (!thisRIFF.File.Type.WAVE64)
    {
        paddingbytes = WordAlign(thisRIFF.File.Total_Bytes.Written) - thisRIFF.File.Total_Bytes.Written;
    }
    else
    {
        paddingbytes = QWordAlign(thisRIFF.File.Total_Bytes.Written) - thisRIFF.File.Total_Bytes.Written;
    }

    return ((paddingbytes == 0) || ((paddingbytes != 0) && (thisRIFF.File.Write(thisRIFF, (char*) &GuidData[GUID_VOID], paddingbytes))));
}


bool ReadPaddingFromFile(tRIFF_Rec &thisRIFF)
{
    uint64_t paddingbytes = 0;

    uint8_t paddingtempdata[8];

    if (!thisRIFF.File.Type.WAVE64)
    {
        paddingbytes = WordAlign(thisRIFF.File.Total_Bytes.Read) - thisRIFF.File.Total_Bytes.Read;
    }
    else
    {
        paddingbytes = QWordAlign(thisRIFF.File.Total_Bytes.Read) - thisRIFF.File.Total_Bytes.Read;
    }

    if (paddingbytes == 0)
    {
        return true;
    }
    else
    {
        return (thisRIFF.File.Read(thisRIFF, (char*) &paddingtempdata, paddingbytes) == paddingbytes);
    }
}


std::string GuidToString(GUID thisGuid)
{
    std::ostringstream ss;

    ss << CardinalToHex(thisGuid.Data1) << "-";
    ss << WordToHex(thisGuid.Data2) << "-";
    ss << WordToHex(thisGuid.Data3) << "-";
    ss << WordToHex((uint32_t)thisGuid.Data4[1] + ((uint32_t)thisGuid.Data4[0] << 8)) << "-";
    ss << WordToHex((uint32_t)thisGuid.Data4[3] + ((uint32_t)thisGuid.Data4[2] << 8));
    ss << CardinalToHex(((uint32_t)thisGuid.Data4[7]) + ((uint32_t)thisGuid.Data4[6] << 8) + ((uint32_t)thisGuid.Data4[5] << 16) + ((uint32_t)thisGuid.Data4[4] << 24));
    return ss.str();
}


void wavIOExitProc(std::string wavIOString, int32_t wavIOCode)
{
    for (uint32_t wi = 0; wi < RIFF.WAVE.Chunks.Current.Free; ++wi)
    {
        tChunkHeader* thisHeader = &RIFF.WAVE.Chunks.Map[wi].Header;
        std::cerr << std::string(thisHeader->ID,4) << "; " << GuidToString(thisHeader->Guid) << "; " << CardinalToHex(thisHeader->Size) << std::endl;
    }

    std::cerr << std::endl << wavIOString << std::endl;

    closeWavIO();

    throw (wavIOCode);
}


void CombineLossyBuffers(int32_t bytes_in_buffer)
{
    int32_t cb_i, cb_j;
    cb_i = bytes_in_buffer / ((int32_t) Global.bytes_per_sample);
    int32_t btrd_sample, corr_sample;
    DATA32 this_sample;

    switch (Global.bytes_per_sample)
    {
        case 1:
            for (cb_j = 0; cb_j < cb_i; ++cb_j)
            {
                RIFF.WAVE.Buffer.Bytes[cb_j] = (int32_t)(nRoundEvenInt64((RIFF.BTRD.Buffer.Bytes[cb_j] - 128) * settings.scaling_factor_inv) + RIFF.CORR.Buffer.Bytes[cb_j]);
            }

            break;


        case 2:
            for (cb_j = 0; cb_j < cb_i; ++cb_j)
            {
                RIFF.WAVE.Buffer.ShortInts[cb_j] = (int32_t)(nRoundEvenInt64(RIFF.BTRD.Buffer.ShortInts[cb_j] * settings.scaling_factor_inv) + RIFF.CORR.Buffer.ShortInts[cb_j]);
            }

            break;


        case 3:
            for (int32_t cb_x = 0; cb_x < bytes_in_buffer; cb_x+=3)
            {
                btrd_sample = RIFF.BTRD.Buffer.Bytes[cb_x] | (RIFF.BTRD.Buffer.Bytes[cb_x+1] << 8) | (RIFF.BTRD.Buffer.SmallInts[cb_x+2] << 16);

                corr_sample = RIFF.CORR.Buffer.Bytes[cb_x] | (RIFF.CORR.Buffer.Bytes[cb_x+1] << 8) | (RIFF.CORR.Buffer.SmallInts[cb_x+2] << 16);

                this_sample.Integer = nRoundEvenInt32(btrd_sample * settings.scaling_factor_inv) + corr_sample;

                RIFF.WAVE.Buffer.Bytes[cb_x]   = this_sample.Bytes[0];
                RIFF.WAVE.Buffer.Bytes[cb_x+1] = this_sample.Bytes[1];
                RIFF.WAVE.Buffer.Bytes[cb_x+2] = this_sample.Bytes[2];
            }

            break;


        case 4:
            for (cb_j = 0; cb_j < cb_i; ++cb_j)
            {
                RIFF.WAVE.Buffer.Integers[cb_j] = (int32_t)(nRoundEvenInt64(RIFF.BTRD.Buffer.Integers[cb_j] * settings.scaling_factor_inv) + RIFF.CORR.Buffer.Integers[cb_j]);
            }

            break;


        case 8:
            for (cb_j = 0; cb_j < cb_i; ++cb_j)
                RIFF.WAVE.Buffer.Int64s[cb_j] = (int64_t)(nRoundEvenInt64(RIFF.BTRD.Buffer.Int64s[cb_j] * settings.scaling_factor_inv) + RIFF.CORR.Buffer.Int64s[cb_j]);

            break;

        default:
            wavIOExitProc("Error: Invalid bytes per sample.",0x12);
    }
}


bool nOpenFile(tRIFF_Rec &thisRIFF, std::string thisname, uint32_t thisfilemode)
{
    switch (thisfilemode)
    {
        case 0:
            //Reset(thisfile, 1);
            thisRIFF.File.Name = thisname.c_str();
            thisRIFF.RIFF_File.open(thisname.c_str(), std::ios::in | std::ios::binary);
            thisRIFF.File.Is.File = true;
            thisRIFF.File.Is.Pipe = false;
            thisRIFF.File.Read = readfrom_file;
            thisRIFF.File.Write = writeto_nullptr;
            thisRIFF.File.Cant.Read = false;
            thisRIFF.File.Cant.Write = true;
            thisRIFF.File.Type.RIFF64 = false;
            thisRIFF.File.Type.WAVE64 = false;
            break;

        case 2:
            //Rewrite(thisfile, 1);
            thisRIFF.File.Name = "STDIN";
            thisRIFF.RIFF_File.open(thisname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            thisRIFF.File.Is.File = true;
            thisRIFF.File.Is.Pipe = false;
            thisRIFF.File.Read = readfrom_nullptr;
            thisRIFF.File.Write = writeto_file;
            thisRIFF.File.Cant.Read = true;
            thisRIFF.File.Cant.Write = false;
            thisRIFF.File.Type.RIFF64 = false;
            thisRIFF.File.Type.WAVE64 = false;
            break;

        default:
            wavIOExitProc("Error: Invalid FileMode!", 0x12);
    }

    thisRIFF.File.Is.Open = thisRIFF.RIFF_File.good();

    return thisRIFF.File.Is.Open;
}


bool WriteChunkHeader(tRIFF_Rec &thisRIFF)
{

    if (!thisRIFF.File.Type.WAVE64)
    {
        tChunk32Header tempChunk32;
        tempChunk32.CNum = thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.CNum;
        tempChunk32.Size = thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.Size;

        if (!thisRIFF.File.Write(thisRIFF, (char*) &tempChunk32, thisRIFF.Chunks.Size.Header))
        {
            if (!thisRIFF.File.Type.RIFF64)
            {
                wavIOExitProc("Error writing 'RIFF' chunk.", 0x12);
            }
            else
            {
                wavIOExitProc("Error writing 'RF64' chunk.", 0x12);
            }
        }
    }
    else
    {
        if (!thisRIFF.File.Write(thisRIFF, (char*) &thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write], thisRIFF.Chunks.Size.Header))
        {
            wavIOExitProc("Error writing 'riff' chunk.", 0x12);
        }
    }

    return true;
}


void WriteRIFFHeader(tRIFF_Rec &thisRIFF)
{
    if (!thisRIFF.File.Type.WAVE64)
    {
        tChunk32Header thisChunk32;
        thisChunk32.CNum = thisRIFF.Chunks.WAV.Header.CNum;
        thisChunk32.Size = thisRIFF.Chunks.WAV.Header.Size;
        thisChunk32.RNum = thisRIFF.Chunks.WAV.RNum;

        if (!thisRIFF.File.Write(thisRIFF, (char*) &thisChunk32, thisRIFF.Chunks.Size.RIFF))
        {
            if (!thisRIFF.File.Type.RIFF64)
            {
                wavIOExitProc("Error writing 'RIFF' chunk.", 0x12);
            }
            else
            {
                wavIOExitProc("Error writing 'RF64' chunk.", 0x12);
            }
        }
    }
    else
    {
        if (!thisRIFF.File.Write(thisRIFF, (char*) &thisRIFF.Chunks.WAV, thisRIFF.Chunks.Size.RIFF))
            wavIOExitProc("Error writing 'riff' chunk.", 0x12);
    }
}


bool WriteChunkHeaderAndPayload(tRIFF_Rec &thisRIFF, void* thispayload)
{
    if (thispayload == nullptr)
    {
        return false;
    }

    if (!WriteChunkHeader(thisRIFF))
    {
        return false;
    }

    if (!thisRIFF.File.Write(thisRIFF, thispayload, thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.Size - thisRIFF.Chunks.Size.Reduction))
    {
        return false;
    }

    return WritePaddingToFile(thisRIFF);
}


bool WriteChunk(tRIFF_Rec &thisRIFF)
{
    if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].DATA == nullptr)
    {
        if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.CNum == 0x20746D66) //"fmt "
        {
            if (!parameters.STDOUTPUT)
                thisRIFF.File.Pos.FMT = thisRIFF.File.Total_Bytes.Written;

            WriteChunkHeaderAndPayload(thisRIFF, &thisRIFF.Chunks.FMT.wFormatTag);

            return true;
        }

        if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.CNum == 0x61746164) //'data'
        {
            if (!parameters.STDOUTPUT)
            {
                thisRIFF.File.Pos.DATA = thisRIFF.File.Total_Bytes.Written;
            }

            return WriteChunkHeader(thisRIFF);
        }

        if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.CNum == 0x34367364) //'ds64'
        {
            if (thisRIFF.File.Type.WAVE64)
                return false;

            if (!parameters.STDOUTPUT)
                thisRIFF.File.Pos.DS64 = thisRIFF.File.Total_Bytes.Written;

            return WriteChunkHeaderAndPayload(thisRIFF, &thisRIFF.Chunks.DS64.RIFFSize);
        }

        if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.CNum == 0x74636166) //"fact"
        {
            if ((!parameters.merging) && (thisRIFF.Chunks.FACT.Header.Size != 0) && (!parameters.STDOUTPUT))
            {
                if (!parameters.STDOUTPUT)
                    thisRIFF.File.Pos.FACT = thisRIFF.File.Total_Bytes.Written;

                return WriteChunkHeaderAndPayload(thisRIFF, &thisRIFF.Chunks.FACT.DATA);
            }

            return true;
        }

        return true;
    }
    else
    {
        return WriteChunkHeaderAndPayload(thisRIFF, thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].DATA);
    }
}


bool WriteChunksUpToData(tRIFF_Rec &thisRIFF)
{
    WriteRIFFHeader(thisRIFF);

    for (thisRIFF.Chunks.Current.Write = 0; thisRIFF.Chunks.Current.Write <= thisRIFF.Chunks.Pos.DATA; thisRIFF.Chunks.Current.Write++)
    {
        if (!WriteChunk(thisRIFF))
        {
            wavIOExitProc("Error writing '" + std::string(thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.ID,4) + "' chunk.", 0x12);
        }
    }

    return true;
}


bool WriteChunksAfterData(tRIFF_Rec &thisRIFF)
{
    for (thisRIFF.Chunks.Current.Write = thisRIFF.Chunks.Pos.DATA + 1; thisRIFF.Chunks.Current.Write <= thisRIFF.Chunks.Pos.LAST; thisRIFF.Chunks.Current.Write++)
    {
        if (!WriteChunk(thisRIFF))
        {
            wavIOExitProc("Error writing '" + std::string(thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Write].Header.ID,4) + "' chunk.", 0x12);
        }
    }

    return true;
}


void RewriteDataChunk(tRIFF_Rec &thisRIFF)
{
    tChunk32Header thisChunk32;

    thisRIFF.File.Pos.CURRENT = thisRIFF.File.Total_Bytes.Written;

    if (!thisRIFF.File.Type.WAVE64)
    {
        if (!thisRIFF.File.Type.RIFF64)
        {
            thisRIFF.RIFF_File.seekp(0);

            thisRIFF.Chunks.WAV.Header.Size = uint32_t(thisRIFF.File.Pos.CURRENT - 8);
            thisChunk32.CNum = thisRIFF.Chunks.WAV.Header.CNum;
            thisChunk32.Size = (uint32_t) thisRIFF.Chunks.WAV.Header.Size;

            thisRIFF.RIFF_File.write((const char*) &thisChunk32, 0x08);

            thisRIFF.RIFF_File.seekp(thisRIFF.File.Pos.DATA);
            thisRIFF.Chunks.DATA.Header.Size = RIFF.WAVE.samplebytesread;

            thisChunk32.CNum = thisRIFF.Chunks.DATA.Header.CNum;
            thisChunk32.Size = (uint32_t) thisRIFF.Chunks.DATA.Header.Size;

            thisRIFF.RIFF_File.write((const char*) &thisChunk32, 0x08);
        }
        else
        {
            thisRIFF.RIFF_File.seekp(0);

            thisRIFF.Chunks.WAV.Header.Size = MAX_uint32_t; //uint32_t(thisRIFF.File.Pos.CURRENT - 8);
            thisChunk32.CNum = thisRIFF.Chunks.WAV.Header.CNum;
            thisChunk32.Size = (uint32_t) thisRIFF.Chunks.WAV.Header.Size;

            thisRIFF.RIFF_File.write((const char*) &thisChunk32, 0x08);

            thisRIFF.RIFF_File.seekp(thisRIFF.File.Pos.DS64);

            thisChunk32.CNum = thisRIFF.Chunks.DS64.Header.CNum;
            thisChunk32.Size = (uint32_t) thisRIFF.Chunks.DS64.Header.Size;

            thisRIFF.RIFF_File.write((const char*) &thisChunk32, 0x08);
            thisRIFF.Chunks.DS64.DATASize = RIFF.WAVE.samplebytesread;
            thisRIFF.Chunks.DS64.RIFFSize = thisRIFF.File.Pos.CURRENT;
            thisRIFF.Chunks.DS64.SampleCount = RIFF.WAVE.samplebytesread / Global.Channels / thisRIFF.wBytesPerSample;
            thisRIFF.RIFF_File.write((const char*) &thisRIFF.Chunks.DS64.RIFFSize,thisChunk32.Size);

            thisRIFF.RIFF_File.seekp(thisRIFF.File.Pos.DATA);

            thisChunk32.CNum = thisRIFF.Chunks.DATA.Header.CNum;
            thisChunk32.Size = MAX_uint32_t;

            thisRIFF.RIFF_File.write((const char*) &thisChunk32, 0x08);
        }
    }
    else
    {
        thisRIFF.RIFF_File.seekp(0);

        thisRIFF.RIFF_File.write((const char*) &thisRIFF.Chunks.WAV, 0x18);

        thisRIFF.RIFF_File.seekp(thisRIFF.File.Pos.DATA);

        thisRIFF.Chunks.DATA.Header.Size = RIFF.WAVE.samplebytesread;

        thisRIFF.RIFF_File.write((const char*) &thisRIFF.Chunks.DATA, 0x18);
    }

    thisRIFF.RIFF_File.seekp(thisRIFF.File.Pos.CURRENT);
}


void SetChunkSizeData(tRIFF_Rec &thisRIFF)
{
    if (!thisRIFF.File.Type.WAVE64)
    {
        thisRIFF.Chunks.Size.Reduction = 0x00;
        thisRIFF.Chunks.Size.Header = 0x08;
        thisRIFF.Chunks.Size.Type = 0x04;
        thisRIFF.Chunks.Size.Size = 0x04;
        thisRIFF.Chunks.Size.RIFF = 0x0C;
    }
    else
    {
        thisRIFF.Chunks.Size.Reduction = 0x18;
        thisRIFF.Chunks.Size.Header = 0x18;
        thisRIFF.Chunks.Size.Type = 0x10;
        thisRIFF.Chunks.Size.Size = 0x08;
        thisRIFF.Chunks.Size.RIFF = 0x28;
    }
}

void DetermineFileType(tRIFF_Rec &thisRIFF)
{
    tChunk32Header thisChunk32;

    std::memset((void*) &thisRIFF.Chunks, 0, sizeof(tWAVEChunks));

    if (thisRIFF.File.Read(thisRIFF, (void*) &thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header, 4) != 4)
    {
        wavIOExitProc("Error Reading WAV file Information.", 0x11);
    }

    //============================================================================================
    // Determine WAVE file type.
    //============================================================================================
    if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.CNum == 0x46464952) //"RIFF"
    {
        thisRIFF.File.Type.WAVE64 = false;
        thisRIFF.File.Type.RIFF64 = false;
    }
    else
        if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.CNum == 0x34364652) //"RF64"
        {
            thisRIFF.File.Type.WAVE64 = false;
            thisRIFF.File.Type.RIFF64 = true;
        }
        else
            if (thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.CNum == 0x66666972) //"riff"
            {
                thisRIFF.File.Type.WAVE64 = true;
                thisRIFF.File.Type.RIFF64 = false;
            }
            else
            {
                wavIOExitProc("File type is neither RIFF, RF64 nor riff.", 0x12);
            }
    //============================================================================================

    SetChunkSizeData(thisRIFF);

    if (!thisRIFF.File.Type.WAVE64)
    {
        thisRIFF.Chunks.WAV.Header.CNum = thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.CNum;

        if (thisRIFF.File.Read(thisRIFF, (void*) &thisChunk32.Size, 0x08) != 0x08)
        {
            wavIOExitProc("Error Reading WAV file Information.", 0x11);
        }

        thisRIFF.Chunks.WAV.Header.Size = (uint64_t) thisChunk32.Size;
        thisRIFF.Chunks.WAV.RNum = (uint32_t) thisChunk32.RNum;

        if (thisRIFF.Chunks.WAV.RNum != 0x45564157) //"WAVE"
        {
            wavIOExitProc("RIFF type is not 'WAVE'.", 0x12);
        }
    }
    else
        {
            thisRIFF.Chunks.WAV.Header.CNum = thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.CNum;

            if (thisRIFF.File.Read(thisRIFF, (void*) &thisRIFF.Chunks.WAV.Header.Guid.Data2,0x24) != 0x24)
            {
                wavIOExitProc("Error Reading WAV64 file Information.", 0x11);
            }

#ifndef _WIN32
            if (GUID_compare(thisRIFF.Chunks.WAV.Header.Guid, GuidData[GUID_riff]))
#else
            if (thisRIFF.Chunks.WAV.Header.Guid != GuidData[GUID_riff])
#endif
            {
                std::cerr << GuidToString(thisRIFF.Chunks.WAV.Header.Guid) << "; " << GuidToString(GuidData[GUID_riff]) << std::endl;
                wavIOExitProc("Error Reading WAV64 file Information.", 0x11);
            }

#ifndef _WIN32
            if (GUID_compare(thisRIFF.Chunks.WAV.RIFFGuid, GuidData[GUID_wave]))
#else
            if (thisRIFF.Chunks.WAV.RIFFGuid != GuidData[GUID_wave])
#endif
            {
                wavIOExitProc("WAVE64 'riff' type is not 'wave'.", 0x12);
            }
        }

    thisRIFF.Chunks.Pos.LAST = thisRIFF.Chunks.Current.Free++;
}


bool ReadChunkHeader(tRIFF_Rec &thisRIFF)
{
    if (!thisRIFF.File.Type.WAVE64)
    {
        tChunk32Header thisChunk32;

        if (thisRIFF.File.Read(thisRIFF, (char*) &thisChunk32, 0x08) != 0x08)
        {
            thisRIFF.File.Cant.Read = true;
            return false;
        }

        thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.CNum = thisChunk32.CNum;
        thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.Size = thisChunk32.Size;
    }
    else
    {
        if (thisRIFF.File.Read(thisRIFF, (char*) &thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header.Guid, 0x18) != 0x18)
        {
            thisRIFF.File.Cant.Read = true;
            return false;
        }
    }
    return true;
}

bool ReadChunkData(tRIFF_Rec& thisRIFF)
{
    tChunkMapRecord* thisMapRecord = &thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free];

    uint32_t thisChunkSize = thisMapRecord->Header.Size  - thisRIFF.Chunks.Size.Reduction;

    if (thisMapRecord->Header.CNum == 0x20746D66) // 'fmt '
    {
        thisRIFF.Chunks.Pos.fmt = thisRIFF.Chunks.Current.Free;

        thisMapRecord->DATA = nullptr;

        thisRIFF.Chunks.FMT.Header = thisMapRecord->Header;

        if ((thisChunkSize != 0x10) && (thisChunkSize != 0x12) && (thisChunkSize != 0x28))
        {
            wavIOExitProc("FMT Chunk incorrect size.", 0x12);
        }

        if (thisRIFF.File.Read(thisRIFF, (char*) &thisRIFF.Chunks.FMT.wFormatTag, thisChunkSize) != thisChunkSize)
        {
            wavIOExitProc("Error reading 'fmt ' chunk.", 0x12);
        }

        if (!parameters.merging)
        {
            thisRIFF.Chunks.Current.Free++;
            thisRIFF.Chunks.FACT.Header.Guid = GuidData[GUID_fact];
            thisRIFF.Chunks.FACT.Header.Size = 0;
            thisRIFF.Chunks.Pos.fact = thisRIFF.Chunks.Current.Free;
            thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].Header = thisRIFF.Chunks.FACT.Header;
            thisRIFF.Chunks.Map[thisRIFF.Chunks.Current.Free].DATA = nullptr;
        }

        return true;
    }


    if (thisMapRecord->Header.CNum == 0x34367364) // 'ds64'
    {
        thisRIFF.Chunks.Pos.ds64 = thisRIFF.Chunks.Current.Free;

        thisMapRecord->DATA = nullptr;

        thisRIFF.Chunks.DS64.Header = thisMapRecord->Header;

        if (thisMapRecord->Header.Size < 0x1C)
        {
            wavIOExitProc("Invalid 'ds64' chunk.", 0x12);
        }

        if (thisRIFF.Chunks.WAV.Header.Size != MAX_uint32_t)
        {
            lossyWAVWarning("File type is RF64 and 'RIFF' size <> 0xFFFFFFFF.");
        }

        if (thisRIFF.File.Read(thisRIFF, (char*) &thisRIFF.Chunks.DS64.RIFFSize, thisRIFF.Chunks.DS64.Header.Size) != thisRIFF.Chunks.DS64.Header.Size)
        {
            wavIOExitProc("Error reading 'ds64' chunk.", 0x12);
        }

        if (!ReadPaddingFromFile(thisRIFF))
        {
            wavIOExitProc("Error reading 'ds64' chunk.", 0x12);
        }

        return true;
    }


    if (thisMapRecord->Header.CNum == 0x61746164) // 'data'
    {
        thisRIFF.Chunks.Pos.DATA = thisRIFF.Chunks.Current.Free;

        thisRIFF.Chunks.DATA.Header = thisMapRecord->Header;

        thisMapRecord->DATA = nullptr;

        if ((thisRIFF.File.Type.RIFF64) && (thisMapRecord->Header.Size != MAX_uint32_t))
        {
            lossyWAVWarning("File type is RF64 and 'data' size <> 0xFFFFFFFF.");
        }

        return true;
    }

    if (thisMapRecord->Header.CNum == 0x74636166) // 'fact'
    {
        thisMapRecord->Header.Guid = GuidData[GUID_fact];

        if (thisRIFF.Chunks.Pos.fact != 0)
        {
            thisRIFF.Chunks.Pos.fact = thisRIFF.Chunks.Pos.fact;
        }
        else
        {
            thisRIFF.Chunks.Pos.fact = thisRIFF.Chunks.Current.Free;
        }

        thisMapRecord->DATA = &RIFF.ChunkDATA[RIFF.EndOfChunkMap];
        RIFF.EndOfChunkMap = RIFF.EndOfChunkMap + QWordAlign(thisChunkSize);

        if (RIFF.EndOfChunkMap > ChunkDATA_array_size)
        {
            wavIOExitProc("Too much WAV information before 'data' chunk!", 0x12);
        }

        thisRIFF.Chunks.FACT.Header = thisMapRecord->Header;

        if (!(thisRIFF.File.Read(thisRIFF, (char*) thisMapRecord->DATA, thisChunkSize) == thisChunkSize))
        {
            wavIOExitProc("Error reading 'fact' chunk", 0x12);
        }

        if (int32_t(thisRIFF.Chunks.FACT.Header.Size) <= FACTChunkMaxSize)
            std::memcpy(thisRIFF.Chunks.FACT.DATA, (char*) thisMapRecord->DATA, thisChunkSize);

        if (!ReadPaddingFromFile(thisRIFF))
        {
            wavIOExitProc("Error reading 'fact' chunk", 0x12);
        }

        if (std::string(thisRIFF.Chunks.FACT.DATA, 8) != "lossyWAV")
        {
            std::cerr << "!!!" << std::endl;
            thisRIFF.Chunks.FACT.Header.Size = 0;
        }
        else
        {
            thisRIFF.Chunks.Map[thisRIFF.Chunks.Pos.fact].DATA = nullptr;
        }

        return true;
    }

    thisMapRecord->DATA = &RIFF.ChunkDATA[RIFF.EndOfChunkMap];
    RIFF.EndOfChunkMap = RIFF.EndOfChunkMap + QWordAlign(thisChunkSize);

    if (RIFF.EndOfChunkMap > ChunkDATA_array_size)
    {
        wavIOExitProc("Too much WAV information before 'data' chunk!", 0x12);
    }
    else
    {
        if (!(thisRIFF.File.Read(thisRIFF, (char*) thisMapRecord->DATA, thisChunkSize) == thisChunkSize))
        {
            wavIOExitProc("Error reading '"+std::string(thisMapRecord->Header.ID,4)+"' chunk.", 0x12);
        }

        if (!ReadPaddingFromFile(thisRIFF))
        {
            wavIOExitProc("Error reading '"+std::string(thisMapRecord->Header.ID,4)+"' chunk.", 0x12);
        }
    }

    return true;
}


void ReadChunksUpToDATA(tRIFF_Rec &thisRIFF)
{
    DetermineFileType(thisRIFF);

    while (thisRIFF.Chunks.Pos.DATA == 0) // Have not yet found 'data' chunk.
    {
        if (!ReadChunkHeader(thisRIFF))
        {
            return;
        }

        if (!ReadChunkData(thisRIFF))
        {
            wavIOExitProc("Error Reading WAV file Information.", 0x11);
        }
        else
        {
            thisRIFF.Chunks.Pos.LAST = thisRIFF.Chunks.Current.Free++;
        }
    }

    if ((thisRIFF.Chunks.FMT.wFormatTag != 1) && ((thisRIFF.Chunks.FMT.wFormatTag == 0xFFFE) && (thisRIFF.Chunks.FMT.SubFormat_Word[0] != 1)))
    {
        wavIOExitProc("WAV Data is not uncompressed PCM.", 0x12);
    }

    if (false)
    {
        for (uint32_t thischunk = 0; thischunk < thisRIFF.Chunks.Current.Free; thischunk++)
        {
            if (!RIFF.WAVE.File.Type.WAVE64)
            {
                std::cerr << std::string(thisRIFF.Chunks.Map[thischunk].Header.ID,4) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].Header.Size) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].DATA[0]) << std::endl;
            }
            else
            {
                std::cerr << GuidToString(thisRIFF.Chunks.Map[thischunk].Header.Guid) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].Header.Size) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].DATA[0]) << std::endl;
            }
        }
    }
}

bool ReadChunksAfterDATA(tRIFF_Rec &thisRIFF)
{
    while (!thisRIFF.File.Cant.Read)
    {
        if (ReadChunkHeader(thisRIFF))
        {
            if (ReadChunkData(thisRIFF))
            {
                thisRIFF.Chunks.Pos.LAST = thisRIFF.Chunks.Current.Free++;
            }
        }
    }

    if (false)
    {
        for (uint32_t thischunk = thisRIFF.Chunks.Pos.DATA+1; thischunk < thisRIFF.Chunks.Current.Free; thischunk++)
        {
            if (!RIFF.WAVE.File.Type.WAVE64)
            {
                std::cerr << std::string(thisRIFF.Chunks.Map[thischunk].Header.ID,4) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].Header.Size) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].DATA[0]) << std::endl;
            }
            else
            {
                std::cerr << GuidToString(thisRIFF.Chunks.Map[thischunk].Header.Guid) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].Header.Size) << "; " << CardinalToHex(thisRIFF.Chunks.Map[thischunk].DATA[0]) << std::endl;
            }
        }
    }

    return true;
}


void WCH(tChunkHeader thisHeader, bool isWAVE64)
{
    if (!isWAVE64)
    {
        std::cerr << std::string(thisHeader.ID,4) << "; " << thisHeader.Size << "; " << std::endl;
    }
    else
    {
        std::cerr << std::string(thisHeader.ID,4) << "; " << GuidToString(thisHeader.Guid) << "; " << thisHeader.Size << "; " << std::endl;
    }
}


void MergeFiles()
{
    uint32_t thisread;
    uint64_t numreads;
    uint32_t lastread;
    uint64_t readcount;
    uint64_t readbyteCount = 0;
    uint64_t DATASize;
    uint64_t RIFFSize;

    std::string temp_str;

    RIFF.EndOfChunkMap = 0;

    nOpenFile(RIFF.BTRD, parameters.WavInpDir + parameters.lossyName, 0);

    ReadChunksUpToDATA(RIFF.BTRD);

    nOpenFile(RIFF.CORR, parameters.WavInpDir + parameters.lwcdfName, 0);

    ReadChunksUpToDATA(RIFF.CORR);

    std::cerr << std::endl << "Merging: " << parameters.lossyName << std::endl;
    std::cerr << "       & " << parameters.lwcdfName << std::endl;

    if (RIFF.BTRD.Chunks.WAV.Header.Size != RIFF.CORR.Chunks.WAV.Header.Size)
    {
        WCH(RIFF.BTRD.Chunks.WAV.Header, RIFF.BTRD.File.Type.WAVE64);
        WCH(RIFF.CORR.Chunks.WAV.Header, RIFF.CORR.File.Type.WAVE64);
        lossyWAVError("Files not same length.", 0x12);
    }

    if (RIFF.BTRD.Chunks.DATA.Header.Size != RIFF.CORR.Chunks.DATA.Header.Size)
    {
        WCH(RIFF.BTRD.Chunks.DATA.Header, RIFF.BTRD.File.Type.WAVE64);
        WCH(RIFF.CORR.Chunks.DATA.Header, RIFF.CORR.File.Type.WAVE64);
        lossyWAVError("Data not same length.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FACT.Header.Size == 0)
    {
        WCH(RIFF.BTRD.Chunks.FACT.Header, RIFF.BTRD.File.Type.WAVE64);
        WCH(RIFF.CORR.Chunks.FACT.Header, RIFF.CORR.File.Type.WAVE64);
        lossyWAVError("FACT Data missing (lossy.wav).", 0x12);
    }

    if (RIFF.CORR.Chunks.FACT.Header.Size == 0)
    {
        WCH(RIFF.BTRD.Chunks.FACT.Header, RIFF.BTRD.File.Type.WAVE64);
        WCH(RIFF.CORR.Chunks.FACT.Header, RIFF.CORR.File.Type.WAVE64);
        lossyWAVError("FACT Data missing (lwcdf.wav).", 0x12);
    }

    if (RIFF.BTRD.Chunks.FACT.Header.Size != RIFF.CORR.Chunks.FACT.Header.Size)
    {
        std::string BTRD_FACT_DATA = std::string(RIFF.BTRD.Chunks.FACT.DATA,RIFF.BTRD.Chunks.FACT.Header.Size);
        std::string CORR_FACT_DATA = std::string(RIFF.CORR.Chunks.FACT.DATA,RIFF.CORR.Chunks.FACT.Header.Size);

        if (BTRD_FACT_DATA != CORR_FACT_DATA)
        {
            lossyWAVError("FACT Data mismatch.", 0x12);
        }
    }

    if (std::string(RIFF.BTRD.Chunks.FACT.DATA, 8) != "lossyWAV")
    {
        std::cerr << std::string(RIFF.BTRD.Chunks.FACT.DATA, 8) << "; " << std::string(RIFF.CORR.Chunks.FACT.DATA, 8) << std::endl;

        lossyWAVError("Not a lossyWAV FACT chunk.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.Header.Size != RIFF.CORR.Chunks.FMT.Header.Size)
    {
        lossyWAVError("FMT DATA: Chunk size mismatch.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.wFormatTag != RIFF.CORR.Chunks.FMT.wFormatTag)
    {
        lossyWAVError("FMT Data: wFormatTag mismatch.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.wChannels != RIFF.CORR.Chunks.FMT.wChannels)
    {
        lossyWAVError("FMT Data: channel mismatch.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.nSamplesPerSec != RIFF.CORR.Chunks.FMT.nSamplesPerSec)
    {
        lossyWAVError("FMT Data: nSamplesPerSec mismatch.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.nAvgBytesPerSec != RIFF.CORR.Chunks.FMT.nAvgBytesPerSec)
    {
        lossyWAVError("FMT Data: nAvgBytesPerSec mismatch.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.nBlockAlign != RIFF.CORR.Chunks.FMT.nBlockAlign)
    {
        lossyWAVError("FMT Data: nBlockAlign mismatch.", 0x12);
    }

    if (RIFF.BTRD.Chunks.FMT.wBitsPerSample != RIFF.CORR.Chunks.FMT.wBitsPerSample)
    {
        lossyWAVError("FMT Data: wBitsPerSample mismatch.", 0x12);
    }

    Global.bytes_per_sample = (RIFF.BTRD.Chunks.FMT.wBitsPerSample + 7) >> 3;

    if ((Global.bytes_per_sample < 1) || (Global.bytes_per_sample > 4))
    {
        lossyWAVError(std::string("Invalid bitdepth: ") + NumToStr(RIFF.BTRD.Chunks.FMT.wBitsPerSample), 0x12);
    }

    if (!RIFF.BTRD.File.Type.WAVE64)
    {
        if (!RIFF.BTRD.File.Type.RIFF64)
        {
            DATASize = RIFF.BTRD.Chunks.DATA.Header.Size;

            if (RIFF.BTRD.Chunks.WAV.Header.Size != MAX_uint32_t)
            {
                RIFF.BTRD.Chunks.WAV.Header.Size -= WordAlign(RIFF.BTRD.Chunks.FACT.Header.Size + 0x08);
            }

            RIFFSize = RIFF.BTRD.Chunks.WAV.Header.Size + 0x08;
        }
        else
        {
            DATASize = RIFF.BTRD.Chunks.DS64.DATASize;

            RIFF.BTRD.Chunks.DS64.RIFFSize -= WordAlign(RIFF.BTRD.Chunks.FACT.Header.Size + 0x08);

            RIFFSize = RIFF.BTRD.Chunks.DS64.RIFFSize + 0x08;
        }
    }
    else
    {
        DATASize = RIFF.BTRD.Chunks.DATA.Header.Size - 0x18;

        RIFF.BTRD.Chunks.WAV.Header.Size -= QWordAlign(RIFF.BTRD.Chunks.FACT.Header.Size);

        RIFFSize = RIFF.BTRD.Chunks.WAV.Header.Size;
    }

    BUFFER_SIZEtouse = BUFFER_SIZE - (BUFFER_SIZE % (Global.bytes_per_sample * RIFF.BTRD.Chunks.FMT.wChannels));

    lastread = DATASize % BUFFER_SIZEtouse;
    numreads = (DATASize / BUFFER_SIZEtouse + ((uint32_t) lastread > 0)) - 1;

    nOpenFile(RIFF.WAVE, parameters.WavOutDir + parameters.wavName, 2);

    RIFF.WAVE.Chunks = RIFF.BTRD.Chunks;

    RIFF.WAVE.File.Type = RIFF.BTRD.File.Type; // Copy file type information AFTER opening file.

    if (!WriteChunksUpToData(RIFF.WAVE))
    {
        lossyWAVError("Writing to output file.", 0x12);
    }

    std::string source(&RIFF.WAVE.Chunks.FACT.DATA[0]), key("--scale ");
    size_t pos = source.find(key);

    if (pos != std::string::npos)
    {
        settings.scaling_factor = std::atof(&source[pos + key.length()]);
    }
    else
    {
        settings.scaling_factor = 1.0;
    }

    if (settings.scaling_factor <= 0.0)
    {
        lossyWAVError(std::string("Invalid scale value: ") + source.substr(pos + 8, 8), 0x12);
    }

    settings.scaling_factor_inv = 1.0 / settings.scaling_factor;

    for (readcount = 0; readcount <= numreads; ++readcount)
    {
        if (readcount == numreads)
        {
            thisread = lastread;
        }
        else
        {
            thisread = BUFFER_SIZEtouse;
        }

        RIFF.BTRD.File.Read(RIFF.BTRD, (char*) &RIFF.BTRD.Buffer, thisread);
        if (!RIFF.BTRD.RIFF_File.gcount())
        {
            lossyWAVError("Reading lossy audio data.",0x12);
        }

        RIFF.CORR.File.Read(RIFF.CORR, (char*) &RIFF.CORR.Buffer, thisread);

        if (!RIFF.CORR.RIFF_File.gcount())
        {
            lossyWAVError("Reading lwcdf audio data.",0x12);
        }

        readbyteCount = readbyteCount + thisread;

        CombineLossyBuffers(thisread);

        RIFF.WAVE.File.Write(RIFF.WAVE, (char*) &RIFF.WAVE.Buffer, thisread);

        std::cerr << "\rProgress: " << (double(RIFF.WAVE.File.Total_Bytes.Written) / RIFFSize * 100)/*: 0: 2*/ << "%     ";
    }

    RIFF.CORR.RIFF_File.close();

    if (RIFF.CORR.RIFF_File.bad())
        lossyWAVError("Closing correction file.", 0x12);

    if (!ReadPaddingFromFile(RIFF.BTRD))
        lossyWAVError("Reading chunk padding from lossy file.", 0x12);

    if (!ReadChunksAfterDATA(RIFF.BTRD))
        lossyWAVError("Reading chunks after data chunk from lossy file.", 0x12);

    RIFF.BTRD.RIFF_File.close();

    if (RIFF.BTRD.RIFF_File.bad())
        lossyWAVError("Closing lossy file.", 0x12);

    RIFF.WAVE.Chunks = RIFF.BTRD.Chunks;

    if (!WritePaddingToFile(RIFF.WAVE))
        lossyWAVError("Writing chunk padding to merged file.", 0x12);

    if (!WriteChunksAfterData(RIFF.WAVE))
        lossyWAVError("Writing chunk data after data chunk to merged file.", 0x12);

    if (!WritePaddingToFile(RIFF.WAVE))
        lossyWAVError("Writing chunk padding to merged file.", 0x12);

    RIFF.WAVE.RIFF_File.close();

    if (RIFF.WAVE.RIFF_File.bad())
        lossyWAVError("Closing merged file.", 0x12);

    std::cerr << "\rProgress: " << (double(RIFF.WAVE.File.Total_Bytes.Written) / RIFFSize * 100)/*: 0: 2*/ << "%     ";
}


void ReadTransfer_One(unsigned char* pB, unsigned char* pEndB)
{
    int32_t iSample = 0;

    while (pB < pEndB)
    {
        for (int32_t iChannel = 0; iChannel < Global.Channels; ++iChannel)
        {
            AudioData.WAVEPTR[NEXT_CODEC_BLOCK][iChannel][iSample].Int64 = int64_t(*pB) - 128;
            pB++;
        }

        ++ iSample;
    }
}


void ReadTransfer_Two(unsigned char* pB, unsigned char* pEndB)
{
    int32_t iSample = 0;

    while (pB < pEndB)
    {
        for (int32_t iChannel = 0; iChannel < Global.Channels; ++iChannel)
        {
            AudioData.WAVEPTR[NEXT_CODEC_BLOCK][iChannel][iSample].Int64 = int64_t(*(short*) pB);
            pB += 2;
        }

        ++ iSample;
    }
}


void ReadTransfer_Three(unsigned char* pB, unsigned char* pEndB)
{
    int32_t iSample = 0;
    DATA32 this32;

    while (pB < pEndB)
    {
        for (int32_t iChannel = 0; iChannel < Global.Channels; ++iChannel)
        {
            this32.Words[0] = (*(uint16_t*) pB);
            pB += 2;
            this32.ShortInts[1] = (*(int8_t*) pB);
            pB++;

            AudioData.WAVEPTR[NEXT_CODEC_BLOCK][iChannel][iSample].Int64 = (int64_t) this32.Integer;
        }

        ++ iSample;
    }
}


void ReadTransfer_Four(unsigned char* pB, unsigned char* pEndB)
{
    int32_t iSample = 0;

    while (pB < pEndB)
    {
        for (int32_t iChannel = 0; iChannel < Global.Channels; ++iChannel)
        {
            AudioData.WAVEPTR[NEXT_CODEC_BLOCK][iChannel][iSample].Int64 = (int64_t)(*(int32_t*) pB);
            pB += 4;
        }

        ++ iSample;
    }
}


bool readNextNextCodecBlock()
{
    unsigned char* pB;
    unsigned char* pStartB;
    uint64_t thisblockreadlength;

    AudioData.Size.Next = 0;

    if (RIFF.WAVE.File.Cant.Read)
    {
        return false;
    }

    if (nrOfBlockInBuffNotYetFetched == 0)  // samples to be fetched aren't in inbuff
    {
        if ((nrOfFullBlocksReadIntoInBuff < nrOfBlocksInReadBuff) || (RIFF.WAVE.sampleByteLeftToRead == 0))  // EOF encountered with last call to readNextSampleBlock
        {
            return false; // no more data
        }

        thisblockreadlength = std::min(BUFFER_SIZEread, RIFF.WAVE.sampleByteLeftToRead);
        RIFF.WAVE.File.Read(RIFF.WAVE, &RIFF.WAVE.Buffer, thisblockreadlength);

        if (!parameters.ignorechunksizes)
        {
            RIFF.WAVE.sampleByteLeftToRead -= RIFF.WAVE.File.Last.Read;
        }

        RIFF.WAVE.samplebytesread = RIFF.WAVE.samplebytesread + RIFF.WAVE.File.Last.Read;
        nrOfFullBlocksReadIntoInBuff = RIFF.WAVE.File.Last.Read / nrOfByteInOneBlockInBuff;
        nrOfByteOfTheLastIncompleteBlock = int32_t(RIFF.WAVE.File.Last.Read % nrOfByteInOneBlockInBuff);
        nrOfBlockInBuffNotYetFetched = 1;
    }

    pB = (unsigned char*)(&RIFF.WAVE.Buffer);
    pB += (nrOfBlockInBuffNotYetFetched - 1) * nrOfByteInOneBlockInBuff;

    pStartB = pB; // points to first sample in Block

    if (nrOfBlockInBuffNotYetFetched <= nrOfFullBlocksReadIntoInBuff)
    {
        pB += nrOfByteInOneBlockInBuff; // pB points to the spot after last sample in Block
        AudioData.Size.Next = uint32_t(nrOfByteInOneBlockInBuff / (((uint32_t) Global.bytes_per_sample * Global.Channels)));
        nrOfBlockInBuffNotYetFetched = nrOfBlockInBuffNotYetFetched + 1;

        if (nrOfBlockInBuffNotYetFetched > nrOfBlocksInReadBuff)
        {
            nrOfBlockInBuffNotYetFetched = 0;
        }
    }
    else
    {
        if (nrOfByteOfTheLastIncompleteBlock == 0)
        {
            return false;
        }

        pB += nrOfByteOfTheLastIncompleteBlock;  // pB points to the spot after last sample in Block
        AudioData.Size.Next = uint32_t(nrOfByteOfTheLastIncompleteBlock / (((uint32_t) Global.bytes_per_sample * Global.Channels)));
        nrOfBlockInBuffNotYetFetched = 0; // force a try for a new Blockread with next call which will lead to a return False
    }

    RIFF.WAVE.ReadTransfer(pStartB, pB);

    return true;
}


void WriteTransfer_One(MultiChannelCodecBlock& outputcodecblock, unsigned char* pB)
{
    for (int32_t wt_i = 0; wt_i < AudioData.Size.This; ++wt_i)
        for (int32_t wt_j = 0; wt_j < Global.Channels; ++wt_j)
        {
            (*pB) = outputcodecblock[wt_j][wt_i].Bytes[0] + 128;
            pB++;
        }
}


void WriteTransfer_Two(MultiChannelCodecBlock& outputcodecblock, unsigned char* pB)
{
    for (int32_t wt_i = 0; wt_i < AudioData.Size.This; ++wt_i)
        for (int32_t wt_j = 0; wt_j < Global.Channels; ++wt_j)
        {
            (*(short*) pB) = outputcodecblock[wt_j][wt_i].Words[0];
            pB+=2;
        }
}


void WriteTransfer_Three(MultiChannelCodecBlock& outputcodecblock, unsigned char* pB)
{
    DATA32 this32;

    for (int32_t wt_i = 0; wt_i < AudioData.Size.This; ++wt_i)
        for (int32_t wt_j = 0; wt_j < Global.Channels; ++wt_j)
        {
            this32.Integer = outputcodecblock[wt_j][wt_i].Integers[0];

            (*(short*) pB) = this32.Words[0];
            pB+=2;
            (*pB) = this32.Bytes[2];
            pB++;
        }
}


void WriteTransfer_Four(MultiChannelCodecBlock& outputcodecblock, unsigned char* pB)
{
    for (int32_t wt_i = 0; wt_i < AudioData.Size.This; ++wt_i)
        for (int32_t wt_j = 0; wt_j < Global.Channels; ++wt_j)
        {
            (*(int32_t*) pB) = outputcodecblock[wt_j][wt_i].Integers[0];
            pB+=4;
        }
}


static void (* WriteTransferProcs[5])(MultiChannelCodecBlock & outputcodecblock, unsigned char* pB) = {nullptr, WriteTransfer_One, WriteTransfer_Two, WriteTransfer_Three, WriteTransfer_Four};
static void (* ReadTransferprocs [5])(unsigned char* pB, unsigned char* pEndB)                               = {nullptr, ReadTransfer_One,  ReadTransfer_Two,  ReadTransfer_Three,  ReadTransfer_Four };


bool writeNextCodecBlock(tRIFF_Rec &thisRIFF, MultiChannelCodecBlock& outputcodecblock)
{
    uint32_t nrOfByteForOutBuff = ((int32_t) Global.bytes_per_sample) * thisRIFF.Chunks.FMT.wChannels * AudioData.Size.This;
    unsigned char* pB;

    if ((thisRIFF.BytesInBuffer + nrOfByteForOutBuff) > BUFFER_SIZEwrite)
    {
        if (!thisRIFF.File.Write(thisRIFF, (char*)&thisRIFF.Buffer, thisRIFF.BytesInBuffer))
            return false;

        thisRIFF.BytesInBuffer = 0;
    }

    pB = (unsigned char*) &thisRIFF.Buffer; // Check this // TY
    pB += thisRIFF.BytesInBuffer;

    thisRIFF.WriteTransfer(outputcodecblock, pB);
    thisRIFF.BytesInBuffer += nrOfByteForOutBuff;

    return true;
}


bool writeNextBTRDcodecblock()
{
    return writeNextCodecBlock(RIFF.BTRD, AudioData.BTRDDATA[AudioData.Rev_LUT[THIS_CODEC_BLOCK]]);
}


bool writeNextCORRcodecblock()
{
    return writeNextCodecBlock(RIFF.CORR, AudioData.CORRDATA[AudioData.Rev_LUT[THIS_CODEC_BLOCK]]);
}


bool FlushRIFFBuffer(tRIFF_Rec &thisRIFF)
{
    if (thisRIFF.BytesInBuffer > 0)
    {
        if (!thisRIFF.File.Write(thisRIFF, (char*) &thisRIFF.Buffer, thisRIFF.BytesInBuffer))
        {
            return false;
        }

        return WritePaddingToFile(thisRIFF);
    }
    else
        return true;
}

bool CopyChunksAndWrite(tRIFF_Rec &thisRIFF)
{
    thisRIFF.Chunks = RIFF.WAVE.Chunks;

    if (!WriteChunksAfterData(thisRIFF))
    {
        return false;
    }

    return WritePaddingToFile(thisRIFF);
}

bool FlushWriteChunksAndClose(tRIFF_Rec &thisRIFF)
{
    if (thisRIFF.File.Is.Open)
    {
        if (!FlushRIFFBuffer(thisRIFF))
        {
            return false;
        }

        if (!CopyChunksAndWrite(thisRIFF))
        {
            return false;
        }

        if ((parameters.STDINPUT) && (!parameters.STDOUTPUT) && (!parameters.ignorechunksizes))
            RewriteDataChunk(thisRIFF);

        if (thisRIFF.File.Is.File)
        {
            thisRIFF.RIFF_File.close();
            thisRIFF.File.Is.Open = false;
            thisRIFF.File.Is.File = false;

            if (thisRIFF.RIFF_File.bad())
                return false;
        }
    }

    return true;
}

bool closeWavIO()
{
    if (!RIFF.WAVE.File.Cant.Read)
    {
        ReadChunksAfterDATA(RIFF.WAVE);
    }

    if (RIFF.WAVE.File.Is.Open)
    {
        if (RIFF.WAVE.File.Is.File)
        {
            RIFF.WAVE.RIFF_File.close();
            RIFF.WAVE.File.Is.Open = false;
            RIFF.WAVE.File.Is.File = false;

            if (RIFF.WAVE.RIFF_File.bad())
                return false;
        }
    }

    if (!FlushWriteChunksAndClose(RIFF.BTRD))
    {
        return false;
    }

    if (!FlushWriteChunksAndClose(RIFF.CORR))
    {
        return false;
    }

//    for (uint32_t sa_i = 0; sa_i<RIFF.WAVE.Chunks.Current.Free; sa_i++)
//    {
//        std::cerr << std::string(RIFF.WAVE.Chunks.Map[sa_i].Header.ID,4) << "; " << std::endl;
//    }

    return true;
}


void Create_FACT_Chunk(tRIFF_Rec &thisRIFF)
{
    date_time_string_make(strings.datestamp,timer.StartTime);

    thisRIFF.Chunks.FACT.Header.Guid = GuidData[GUID_fact];

    DateTimeString = version_string + strings.version_short + " @ " +
                     strings.datestamp + ", " +
                     strings.parameter + char(13) + char(10) + char(0);

    if (DateTimeString.length() && 1)
        DateTimeString = DateTimeString + char(0);

    std::memcpy(thisRIFF.Chunks.FACT.DATA, &DateTimeString[0], DateTimeString.length()); // NB! index 0 is first char

    if (!thisRIFF.File.Type.WAVE64)
        thisRIFF.Chunks.FACT.Header.Size = DateTimeString.length();
    else
        thisRIFF.Chunks.FACT.Header.Size = DateTimeString.length() + 0x18;

    thisRIFF.Chunks.Map[thisRIFF.Chunks.Pos.fact].Header = thisRIFF.Chunks.FACT.Header;
}


bool openWavIO()
{
    Global.Codec_Block.Size = 0;

    if (parameters.STDINPUT)
    {
        #ifdef _WIN32
        _setmode(STDIN_FILENO, _O_BINARY);
        #endif

        RIFF.WAVE.File.Read = readfrom_stdin;
        RIFF.WAVE.File.Is.Open = true;
        RIFF.WAVE.File.Is.Pipe = true;
        RIFF.WAVE.File.Is.File = false;
    }
    else
    {
        if (!nOpenFile(RIFF.WAVE, parameters.WavInpDir + parameters.wavName, 0))
            return false;
    }

    RIFF.EndOfChunkMap = 0;

    ReadChunksUpToDATA(RIFF.WAVE);

    if ((RIFF.WAVE.Chunks.WAV.Header.Size == 0) && (!parameters.ignorechunksizes))
    {
        wavIOExitProc("Zero sample WAV file!?!", 0x12);
    }

    Global.bits_per_sample = RIFF.WAVE.Chunks.FMT.wBitsPerSample;

    RIFF.WAVE.wBytesPerSample = (RIFF.WAVE.Chunks.FMT.wBitsPerSample + 7) >> 3;

    Global.bytes_per_sample = RIFF.WAVE.wBytesPerSample;

    if ((Global.bytes_per_sample < 1) || (Global.bytes_per_sample > 4))
    {
        lossyWAVError(std::string("Invalid bitdepth: ") + NumToStr(RIFF.BTRD.Chunks.FMT.wBitsPerSample), 0x12);
    }

    RIFF.WAVE.ReadTransfer = ReadTransferprocs[Global.bytes_per_sample];
    RIFF.WAVE.WriteTransfer = WriteTransferProcs[Global.bytes_per_sample];

    Global.Channels = RIFF.WAVE.Chunks.FMT.wChannels;
    Global.sample_rate = RIFF.WAVE.Chunks.FMT.nSamplesPerSec;
    Global.sample_rate_recip = 1.0 / Global.sample_rate;

    if (!RIFF.WAVE.File.Type.RIFF64)
    {
        Global.Total_Samples = uint64_t((1.0 * RIFF.WAVE.Chunks.DATA.Header.Size) / Global.Channels / Global.bytes_per_sample);

        Global.WAVE_size = RIFF.WAVE.Chunks.WAV.Header.Size;

        RIFF.WAVE.Chunks.DS64.Header.Guid = GuidData[GUID_VOID];

        RIFF.WAVE.Chunks.DS64.Header.Size = 0;

        RIFF.WAVE.Chunks.DS64.RIFFSize = 0;
        RIFF.WAVE.Chunks.DS64.DATASize = 0;
        RIFF.WAVE.Chunks.DS64.SampleCount = 0;
        RIFF.WAVE.Chunks.DS64.TableLength = 0;
    }
    else
    {
        Global.Total_Samples = RIFF.WAVE.Chunks.DS64.SampleCount;
        Global.WAVE_size = RIFF.WAVE.Chunks.DS64.RIFFSize;
    }

    if (Global.sample_rate < 22050)
    {
        std::cerr << "Sample Rate too low : " << NumToStr(Global.sample_rate * OneOver[1000], 2) << "kHz (Min=22.05kHz).\n";
        wavIOExitProc("This is not supported.", 0x12);
    }

    if (Global.sample_rate > 409600)
    {
        std::cerr << "Sample Rate too high : " << NumToStr(Global.sample_rate * OneOver[1000], 2) << "kHz (Max=409.6kHz).\n";
        wavIOExitProc("This is not supported.", 0x12);
    }

    Global.Codec_Block.duration = OneOver[100];
    Global.Codec_Block.bits = std::max(8, nRoundEvenInt32(nlog2(Global.Codec_Block.duration * Global.sample_rate)));
    Global.Codec_Block.Size = PowersOf.TwoInt64[Global.Codec_Block.bits];
    Global.Codec_Block.Size_recip = PowersOf.TwoX[TWO_OFFSET + -Global.Codec_Block.bits];
    Global.Codec_Block.bit_shift = Global.Codec_Block.bits - 9;
    Global.Codec_Block.duration = double(Global.Codec_Block.Size) / Global.sample_rate;

    //========================================================================
    // Probably totally unnecessary double check on maximum FFT length
    //========================================================================
    if (Global.Codec_Block.bits > MAX_FFT_BIT_LENGTH - 1)
    {
        std::cerr << "Sample Rate too high : " << NumToStr(Global.sample_rate * OneOver[1000], 2) << "kHz (Max=409.6kHz).\n";
        wavIOExitProc("This is not supported.", 0x12);
    }
    //========================================================================

    //========================================================================
    // Pre-calculate variables used in FFT and others.
    //========================================================================
    for (int32_t this_bit = 1; this_bit <= (MAX_FFT_BIT_LENGTH + 1); this_bit++)
    {
        FFT_PreCalc_Data_Rec[this_bit].bit_length = this_bit;
        FFT_PreCalc_Data_Rec[this_bit].bit_shift_from_max = MAX_FFT_BIT_LENGTH - FFT_PreCalc_Data_Rec[this_bit].bit_length;
        FFT_PreCalc_Data_Rec[this_bit].bit_shift_from_32 = 32 - FFT_PreCalc_Data_Rec[this_bit].bit_length;
        FFT_PreCalc_Data_Rec[this_bit].threshold_shift = (Global.Codec_Block.bits - FFT_PreCalc_Data_Rec[this_bit].bit_length) * 0.5;
        FFT_PreCalc_Data_Rec[this_bit].length = PowersOf.TwoInt64[FFT_PreCalc_Data_Rec[this_bit].bit_length];
        FFT_PreCalc_Data_Rec[this_bit].length_m1 = FFT_PreCalc_Data_Rec[this_bit].length - 1;
        FFT_PreCalc_Data_Rec[this_bit].length_half = FFT_PreCalc_Data_Rec[this_bit].length >> 1;
        FFT_PreCalc_Data_Rec[this_bit].length_half_m1 = FFT_PreCalc_Data_Rec[this_bit].length_half - 1;
        FFT_PreCalc_Data_Rec[this_bit].length_recip = float(PowersOf.TwoX[TWO_OFFSET + - FFT_PreCalc_Data_Rec[this_bit].bit_length]);
        FFT_PreCalc_Data_Rec[this_bit].length_half_recip = FFT_PreCalc_Data_Rec[this_bit].length_recip * 2;
    }
    //========================================================================

    Global.Codec_Block.Total = (uint64_t)(1. * (Global.Total_Samples + Global.Codec_Block.Size - 1.) * Global.Channels * Global.Codec_Block.Size_recip);

    if (parameters.ignorechunksizes)
    {
        RIFF.WAVE.sampleByteLeftToRead = MAX_uint64_t;
        RIFF.WAVE.samplesLeftToRead = MAX_uint64_t;
    }
    else
    {
        if (RIFF.WAVE.File.Type.RIFF64)
        {
            RIFF.WAVE.sampleByteLeftToRead = RIFF.WAVE.Chunks.DS64.DATASize;
            RIFF.WAVE.samplesLeftToRead = RIFF.WAVE.Chunks.DS64.SampleCount;
        }
        else
        {
            if (RIFF.WAVE.File.Type.WAVE64)
            {
                RIFF.WAVE.sampleByteLeftToRead = RIFF.WAVE.Chunks.DATA.Header.Size - 0x18;
                RIFF.WAVE.samplesLeftToRead = RIFF.WAVE.sampleByteLeftToRead / (Global.Channels * Global.bytes_per_sample);
            }
            else
            {
                RIFF.WAVE.sampleByteLeftToRead = RIFF.WAVE.Chunks.DATA.Header.Size;
                RIFF.WAVE.samplesLeftToRead = RIFF.WAVE.sampleByteLeftToRead / (Global.Channels * Global.bytes_per_sample);

                if ((RIFF.WAVE.sampleByteLeftToRead == MAX_uint32_t) || (RIFF.WAVE.Chunks.WAV.Header.Size == MAX_uint32_t))
                {
                    RIFF.WAVE.sampleByteLeftToRead = MAX_uint64_t;
                    RIFF.WAVE.samplesLeftToRead = MAX_uint64_t;
                }

                if ((RIFF.WAVE.sampleByteLeftToRead == 0) || (RIFF.WAVE.Chunks.WAV.Header.Size == 0))
                {
                    RIFF.WAVE.sampleByteLeftToRead = MAX_uint64_t;
                    RIFF.WAVE.samplesLeftToRead = MAX_uint64_t;
                }
            }
        }
    }

    if (parameters.output.detail)
    {
        if (Global.Codec_Block.Total <= 0x10000000)
        {
            bit_removal_history = new (std::nothrow) uint8_t [(Global.Codec_Block.Total + 0x0F) & 0xFFFFFFF0];

            if (bit_removal_history == nullptr)
            {
                parameters.output.detail = false;
                lossyWAVWarning("Not enough memory to retain bit-removal history.");
            }
        }
        else
        {
            lossyWAVWarning("File too large to retain bit-removal history.");
            parameters.output.detail = false;
        }
    }

    nrOfByteInOneBlockInBuff = Global.bytes_per_sample * Global.Channels * Global.Codec_Block.Size;


    //========================================================================
    // Define how many codec-blocks of data to buffer.
    //========================================================================
    if (parameters.STDINPUT)
        nrOfBlocksInReadBuff = 1;
    else
        nrOfBlocksInReadBuff = (BUFFER_SIZE / nrOfByteInOneBlockInBuff);

    if (parameters.STDOUTPUT)
        nrOfBlocksInWriteBuff = 1;
    else
        nrOfBlocksInWriteBuff = (BUFFER_SIZE / nrOfByteInOneBlockInBuff);
    //========================================================================

    BUFFER_SIZEread = nrOfByteInOneBlockInBuff * nrOfBlocksInReadBuff;
    BUFFER_SIZEwrite = nrOfByteInOneBlockInBuff * nrOfBlocksInWriteBuff;

    //========================================================================
    // If checking then output relevant message and exit, else error if FACT.
    //========================================================================
    if (parameters.checking)
    {
        if ((RIFF.WAVE.Chunks.FACT.Header.Size > 0) && (std::string(&RIFF.WAVE.Chunks.FACT.DATA[0], 8) == "lossyWAV"))
        {
            if ((parameters.output.verbosity) && (!parameters.output.silent))
            {
                std::cerr << "lossyWAV FACT Chunk found. File already processed.\n";
                std::cerr << RIFF.WAVE.Chunks.FACT.DATA;
            }
            lossyWAVError("", 0x10);
        }
        else
        {
            if ((parameters.output.verbosity) && (!parameters.output.silent))
                std::cerr << "lossyWAV FACT Chunk not found. File not marked as processed.\n";

            lossyWAVError("", 0x00);
        }
    }
    else
    {
        if (std::string(&RIFF.WAVE.Chunks.FACT.DATA[0], 8) == "lossyWAV")
        {
            std::cerr << RIFF.WAVE.Chunks.FACT.DATA;
            lossyWAVError("lossyWAV FACT Chunk found. File already processed.", 0x10);
        }
    }
    //========================================================================

    //========================================================================
    // Create FACT chunk and add to size of RIFF file.
    //========================================================================
    Create_FACT_Chunk(RIFF.WAVE);

    if (!RIFF.WAVE.File.Type.WAVE64)
    {
        if (!RIFF.WAVE.File.Type.RIFF64)
        {
            if (RIFF.WAVE.Chunks.WAV.Header.Size < MAX_WAVE_SIZE)
            {
                RIFF.WAVE.Chunks.WAV.Header.Size += WordAlign(RIFF.WAVE.Chunks.FACT.Header.Size + 0x08);
            }
        }
        else
        {
            RIFF.WAVE.Chunks.DS64.RIFFSize += WordAlign(RIFF.WAVE.Chunks.FACT.Header.Size + 0x08);
        }
    }
    else
    {
        RIFF.WAVE.Chunks.WAV.Header.Size += QWordAlign(RIFF.WAVE.Chunks.FACT.Header.Size);
    }
    //========================================================================

    RIFF.WAVE.firstsamplesread = std::min(RIFF.WAVE.sampleByteLeftToRead, BUFFER_SIZEread);

    if (RIFF.WAVE.File.Read(RIFF.WAVE, &RIFF.WAVE.Buffer, RIFF.WAVE.firstsamplesread) < RIFF.WAVE.firstsamplesread)
    {
        return false;
    }

    RIFF.WAVE.samplebytesread = RIFF.WAVE.File.Last.Read;

    if (!parameters.ignorechunksizes)
    {
        RIFF.WAVE.sampleByteLeftToRead -= RIFF.WAVE.File.Last.Read;
    }

    nrOfFullBlocksReadIntoInBuff = RIFF.WAVE.File.Last.Read / nrOfByteInOneBlockInBuff;
    nrOfByteOfTheLastIncompleteBlock = RIFF.WAVE.File.Last.Read - nrOfFullBlocksReadIntoInBuff * nrOfByteInOneBlockInBuff;
    nrOfBlockInBuffNotYetFetched = 1;

    RIFF.BTRD.Chunks = RIFF.WAVE.Chunks;
    RIFF.BTRD.wBytesPerSample = (RIFF.BTRD.Chunks.FMT.wBitsPerSample + 7) >> 3;
    RIFF.BTRD.ReadTransfer  = ReadTransferprocs[RIFF.BTRD.wBytesPerSample];
    RIFF.BTRD.WriteTransfer = WriteTransferProcs[RIFF.BTRD.wBytesPerSample];
    RIFF.BTRD.BytesInBuffer = 0;

    RIFF.CORR.Chunks = RIFF.WAVE.Chunks;
    RIFF.CORR.wBytesPerSample = (RIFF.CORR.Chunks.FMT.wBitsPerSample + 7) >> 3;
    RIFF.CORR.ReadTransfer  = ReadTransferprocs[RIFF.CORR.wBytesPerSample];
    RIFF.CORR.WriteTransfer = WriteTransferProcs[RIFF.CORR.wBytesPerSample];
    RIFF.CORR.BytesInBuffer = 0;

    if (parameters.STDOUTPUT)
    {
        #ifdef _WIN32
        _setmode(STDOUT_FILENO, _O_BINARY);
        #endif

        RIFF.BTRD.File.Write = writeto_stdout;
        RIFF.BTRD.File.Is.Pipe = true;
        RIFF.BTRD.File.Is.File = false;
        RIFF.BTRD.File.Is.Open = true;
    }
    else
        if (!nOpenFile(RIFF.BTRD, parameters.WavOutDir + parameters.lossyName, 2))
            return false;

    RIFF.BTRD.File.Type = RIFF.WAVE.File.Type;
    if (!WriteChunksUpToData(RIFF.BTRD))
        lossyWAVError("Writing to output file.", 0x12);


    if ((!parameters.STDOUTPUT) && (parameters.correction))
    {
        if (!nOpenFile(RIFF.CORR, parameters.WavOutDir + parameters.lwcdfName, 2))
            return false;

        RIFF.CORR.File.Type = RIFF.WAVE.File.Type;

        if (!WriteChunksUpToData(RIFF.CORR))
            lossyWAVError("Writing to output file.", 0x12);
    }

    RIFF.WAVE.ID = 0;
    RIFF.BTRD.ID = 1;
    RIFF.CORR.ID = 2;

    return true;
}


void nWAV_Init()
{
    RIFF.ChunkDATA = new (std::nothrow) uint8_t[ChunkDATA_array_size];

    if (RIFF.ChunkDATA == nullptr)
    {
        lossyWAVError("Initialising WAV Chunk Storage", 0x12);
    }
}


void nWAV_Cleanup()
{
    if (RIFF.ChunkDATA != nullptr)
    {
        delete[] RIFF.ChunkDATA;
    }
}
