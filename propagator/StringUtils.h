/*
 * File: StringUtils.h
 * Yezheng (Jim) HU
 * History:
 * May 01 2015 created
 * May 01 2015 borrow convertIntToFixLenString from BSssfMig2D.cpp
 * June 11 2015 add convertFloatToString
 * June 16 2015 add removeExtension, removeNamedExtension and removeNamedExtensions
 * June 17 2015 add filesWithExtensionsExist
 * March 30 2016 add AllfilesWithExtensionsExist and RetrieveFileListsFromStringCards
 */

#ifndef STRING_UTILITY_H
#define STRING_UTILITY_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
const std::string convertIntToFixLenString(const int number, const int len);
const std::string convertFloatToString(const float &number);
const std::string Trim(const std::string str);
const std::string removeExtension(const std::string& filename);
bool StringEndWithExtension(const std::string target, const std::string extension);
const std::string removeSlashIfAtEnd(const std::string& filename);
const std::string removeNamedExtension(const std::string& filename, const std::string &extension);
const std::string removeNamedExtensions(const std::string& filename, const std::vector<std::string> &extensions);
void getComplementLabels(const std::vector<std::string> & srcLabels, const std::vector<std::string> &targetlabels, std::vector<std::string> &complementLabels);
void getLabelDimensionIndex(const std::vector<std::string> & srcLabels, const std::vector<std::string> & targetLabels, std::vector<size_t> &targetIndex);
#endif
