/*
 * File: StringUtils.cpp
 * Yezheng (Jim) HU
 * History:
 * May 01 2015 created
 * May 01 2015 borrow convertIntToFixLenString from BSssfMig2D.cpp
 * June 11 2015 add convertFloatToString
 * June 16 2015 add removeExtension, removeNamedExtension and removeNamedExtensions
 * June 17 2015 add filesWithExtensionsExist
 * 1 April 2016 add add FilterFileLists
 */
#include <algorithm>
#include "StringUtils.h"
const std::string convertIntToFixLenString(const int number, const int len)
{
    std::stringstream ss;
    ss << std::setw(len) << std::setfill('0') << number;

    return ss.str();
}

const std::string convertFloatToString(const float & number)
{
   //std::stringstream ss;
   char ss[1024];
   sprintf(ss, "%g", number);
   return std::string(ss);
}

const std::string Trim(const std::string str)
{
    std::string midman(str.c_str());
    size_t first = midman.find_first_not_of(' ');
    size_t last = midman.find_last_not_of(' ');
    std::string dummy;
    if (std::string::npos != first && std::string::npos != last) {
       dummy = midman.substr(first, (last-first+1));
    } else {
       dummy = std::string("");
    }
    return dummy;
}

const std::string removeExtension(const std::string& filename)
{
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
}

bool StringEndWithExtension(const std::string target, const std::string extension)
{
    size_t targetLength = target.length();
    size_t extensionLength = extension.length();
    return (targetLength > extensionLength ? (target.substr(targetLength - extensionLength, extensionLength).compare(extension) == 0) : 0);
}

const std::string removeSlashIfAtEnd(const std::string& filename)
{
    size_t lastslash = filename.find_last_of("/");
    if (lastslash == std::string::npos) return filename;
    if (lastslash == filename.length()-1) {
       return filename.substr(0, lastslash); 
    } else {
       return filename;
    }
}

const std::string removeNamedExtension(const std::string& filename, const std::string &extension)
{
    if (StringEndWithExtension(filename, extension) ) {
       return removeExtension(filename);
    } else {
       return (filename);
    }
}

const std::string removeNamedExtensions(const std::string& filename, const std::vector<std::string> &extensions)
{
    for (std::vector<std::string>::const_iterator it = extensions.begin(); it!= extensions.end(); ++it) {
        if (StringEndWithExtension(filename, *it) ) {
            return removeExtension(filename);
        }
    }
    return (filename);
}

void getComplementLabels(const std::vector<std::string> & srcLabels, const std::vector<std::string> &targetlabels, std::vector<std::string> &complementLabels)
{
    complementLabels.resize(0);
    for (size_t di = 0; di<srcLabels.size(); di++) {
        if (std::find(targetlabels.begin(), targetlabels.end(), srcLabels[di]) == targetlabels.end()) {
           complementLabels.push_back(srcLabels[di]);
        }
    }
}

void getLabelDimensionIndex(const std::vector<std::string> & srcLabels, const std::vector<std::string> & targetLabels, std::vector<size_t> &targetIndex)
{
    targetIndex.resize(0);
    for (size_t it=0; it<targetLabels.size(); it++) {
        std::vector<std::string>::const_iterator ip = std::find(srcLabels.begin(), srcLabels.end(), targetLabels[it]);
        if (ip == srcLabels.end()) {
           targetIndex.push_back(std::string::npos);
        } else {
           size_t idx = ip - srcLabels.begin();
           targetIndex.push_back(idx);
        }
    }
}
