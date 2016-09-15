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

#ifndef nSupport_h_
#define nSupport_h_

#include <ctime>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

// DELPHI EMULATION

inline std::string NumToStr(int32_t number) // overload
{
    std::ostringstream ss;

    ss << number;

    return ss.str();
}


inline std::string NumToStr(double number, int32_t digits) // overload
{
    std::ostringstream ss;

    ss << std::fixed << std::setprecision(digits) << std::setfill('0') << std::left << number;

    return ss.str();
}


inline bool FileIsReadOnly(const std::string& fileName_in)
{
    DWORD fAttrib = GetFileAttributesA(fileName_in.c_str());

    return (!(fAttrib == INVALID_FILE_ATTRIBUTES) && (fAttrib && FILE_ATTRIBUTE_READONLY > 0));
}

inline bool SetFileReadWrite(const std::string& fileName_in)
{
    DWORD fAttrib = GetFileAttributesA(fileName_in.c_str());

    if (fAttrib == INVALID_FILE_ATTRIBUTES)
    {
        return false;  //something is wrong with your path!
    }

    return SetFileAttributesA(fileName_in.c_str(), fAttrib & (0xFFFFFFFF ^ FILE_ATTRIBUTE_READONLY));
}



inline bool DirectoryExists(const std::string& dirName_in)
{
    DWORD fAttrib = GetFileAttributesA(dirName_in.c_str());

    return (!(fAttrib == INVALID_FILE_ATTRIBUTES) && (fAttrib && FILE_ATTRIBUTE_DIRECTORY > 0));
}


inline bool FileExists(std::string file_name)
{
    std::fstream filetest(file_name.c_str(), std::ios_base::in | std::ios_base::binary);

    if (filetest.is_open())
    {
        filetest.close();

        return true;
    }

    return false;
}

#endif // nSupport_h
