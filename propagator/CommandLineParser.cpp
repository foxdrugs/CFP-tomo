/* 
 * File:   CommandLineParser.cpp
 * Author: Yezheng (Jim) Hu
 * June 23 2017 rewritten using multimap and some boost utility
 */

#include "CommandLineParser.h"
#include "StringUtils.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;

#define PRINT_PARAMETERS 1

/**
 * Construct the Object with two parameters
 * @param num : number of parameter groups in command line
 * @param argv :parameter groups separated by blank space
 * 
 */
CommandLineParser::CommandLineParser(int num, char** argv) {
  #if PRINT_PARAMETERS
    for (int i = 0; i < num; i++) {
      std::cout << "argv[" << i << "]=" << argv[i] << std::endl;
    }
  #endif
    setParameterValuePair(num, argv);
}

/**
 * print parameters for QC
 */
void CommandLineParser::printParameters()
{
    std::cout << "The valid parameters after command line:" << "\n";
    for (std::multimap<std::string, std::string>::iterator it=m_parValue.begin(); it!= m_parValue.end(); ++it) {
        std::cout << it->first << "=" << it->second << std::endl;
    }
}

bool CommandLineParser::ParameterNameExists(const std::string& parName)
{
    return (m_parValue.find(parName) != m_parValue.end());
}

/*
 * Usage:
 * strPar=strValue
 */
bool CommandLineParser::getStringValue(const std::string& parName, std::string& parValue)
{
    parValue = "";
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       parValue = it->second;
       parExist = true;
    }
    return (parExist && (!parValue.empty()));
}

/*
 * Usage: strPar=strValue1,strValue2,strValue3,... strPar=strValuea,strValueb,strValuec,...
 * Example:
 * infile=file1,file2,file3 infile=file_a,file_b,file_c
 */
std::vector<std::string> CommandLineParser::getStringVectorValue(const std::string& parName, const std::string& delimeter, bool required)
{
    if (required && !ParameterNameExists(parName)) {
       std::cerr << "No value related to parameter name: " << parName << " exists!" << std::endl; 
       exit(EXIT_FAILURE);
    }
    std::pair<std::multimap<std::string, std::string>::iterator, std::multimap<std::string, std::string>::iterator> itRange;
    itRange = m_parValue.equal_range(parName);
    std::vector<std::string> stringList;
    for (std::multimap<std::string, std::string>::iterator it=itRange.first; it!=itRange.second; ++it) {
        std::string strValue = it->second;
        std::vector <std::string> parsVectorString = convertParStringToVector(strValue, delimeter);
        stringList.insert(stringList.end(), parsVectorString.begin(), parsVectorString.end());
    }
    size_t size = stringList.size();
    if (size == 0 && required) {
       std::cerr << std::endl;
       std::cerr << "Error: during parameter parse stage, an error happens:" << std::endl;
       std::cerr << "Error: input parameter number of " << parName << " is  0." << std::endl;
       exit(EXIT_FAILURE);
    }
    return stringList;
}

/* As above but equires at least num of string values and result left as a parameter*/
bool CommandLineParser::getStringVectorValue(const std::string& parName,
           const std::string& delimeter, const int& num, std::vector<std::string>& stringTagVector)
{
    stringTagVector = getStringVectorValue(parName, delimeter);
    int size = stringTagVector.size();

    bool stringTagValid = true;
    if (size < num) {
        std::cerr << std::endl;
        std::cerr << "Error: during parameter parse stage, an error happens:" << std::endl;
        std::cerr << "Error: input parameter number of " << parName << " is less " << num << std::endl;
        stringTagValid = false;
  #if 0
    } else if (size > num) {
        stringTagVector.resize(num);
  #endif
    }
    return stringTagValid;
}

/* Usage: intPar=intValue */
bool CommandLineParser::getIntValue(const std::string& parName, int& parValue)
{
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       parExist = num_valid<int>(it->second.c_str(), parValue);
    }
    return parExist;
}

/*
 * Usage: intPar=iv1,iv2,iv3
 */
bool CommandLineParser::getIntVectorValue(const std::string& parName, const std::string& delimeter, std::vector<int> &intVectorPar)
{
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       std::string strValue = it->second;
       std::vector<std::string> parsVector = convertParStringToVector(strValue, delimeter);
        int size = parsVector.size();
	if (size > 0) {
           intVectorPar.resize(0);
           parExist = true;
           for (int i = 0; i < size; i++) {
              int itemValue = 0.0;
              if (num_valid<int>(parsVector.at(i).c_str(), itemValue) ) {
                 intVectorPar.push_back(itemValue);
              } else {
                 intVectorPar.clear();
                 fprintf(stderr, "%s: Failed to convert item %s to a double value\n", parName.c_str(), parsVector.at(i).c_str());
                 parExist = false;
                 break;
              }
           }
        } else {
           std::cerr << "\n";
           std::cerr << "Error: during parameter parse stage, an error happens:" << "\n";
           std::cerr << "Error: input parameter number of " << parName << " is 0." << "\n";
        }
    }
    return parExist;
}

/*
 * Usage:
 * intPar=iv1,iv2,iv3
 * However require at least num of integer values
 */
bool CommandLineParser::getIntVectorValue(const std::string& parName, const std::string& delimeter, const int& num,
                                                        std::vector<int>& intVectorPar)
{
    bool vectorValid = false;
    bool parExist = getIntVectorValue(parName, delimeter, intVectorPar);
    if (parExist) {
       int size = intVectorPar.size();
       vectorValid = true;
       if (size < num) {
          std::cerr << "\n";
          std::cerr << "Error: during parameter parse stage, an error happens:" << "\n";
          std::cerr << "Error: input parameter number of " << parName << " is less " << num << "\n";
          vectorValid = false;
     #if 0
       } else if (size > num) {
          intVectorPar.resize(num);
     #endif
       }
    }
    return vectorValid;
}

/*
 * Usage: fltPar=fltValue
 */
bool CommandLineParser::getFloatValue(const std::string & parName, float& parValue)
{
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       parExist = num_valid<float>(it->second.c_str(), parValue);
    }
    return parExist;
}

/*
 * Usage:
 * fltPar=fv1,fv2,fv3
 */
bool CommandLineParser::getFloatVectorValue(const std::string& parName, const std::string& delimeter, std::vector<float> &fltVectorPar)
{
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       std::string strValue = it->second;
       std::vector<std::string> parsVector = convertParStringToVector(strValue, delimeter);
       int size = parsVector.size();
       if (size > 0) {
          fltVectorPar.resize(0);
          parExist = true;
          for (int i = 0; i < size; i++) {
              float itemValue = 0.0;
              if (num_valid<float>(parsVector.at(i).c_str(), itemValue) ) {
                 fltVectorPar.push_back(itemValue);
              } else {
                 fltVectorPar.clear();
                 fprintf(stderr, "%s: Failed to convert item %s to a float value\n", parName.c_str(), parsVector.at(i).c_str());
                 parExist = false;
                 break;
              }
          }
       } else {
          std::cerr << "\n";
          std::cerr << "Error: during parameter parse stage, an error happens:" << "\n";
          std::cerr << "Error: input parameter number of " << parName << " is 0." << "\n";
       }
 #if 0
    } else {
       fprintf(stderr, "Parameter %s does not exist or its value string is empty\n", parName.c_str());
 #endif
    }
    return parExist;
}

/*
 * Usage:
 * fltPar=fv1,fv2,fv3
 * However require at least num of integer values
 */
bool CommandLineParser::getFloatVectorValue(const std::string& parName, const std::string& delimeter, const int& num,
                                                    std::vector<float> &fltVectorPar)
{
    bool vectorValid = false;
    bool parExist = getFloatVectorValue(parName, delimeter, fltVectorPar);
    if (parExist) {
       int size = fltVectorPar.size();
       vectorValid = true;
       if (size < num) {
          std::cerr << "\n";
          std::cerr << "Error: during parameter parse stage, an error happens:" << "\n";
          std::cerr << "Error: input parameter number of " << parName << " is less " << num << "\n";
          vectorValid = false;
     #if 0
       } else if (size > num) {
          fltVectorPar.resize(num);
     #endif
       }
    }
    return vectorValid;
}

/*
 * Usage: dblPar=dblValue
 */
bool CommandLineParser::getDoubleValue(const std::string & parName, double& parValue)
{
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       parExist = num_valid<double>(it->second.c_str(), parValue);
    }
    return parExist;
}

/*
 * Usage:
 * dblPar=dv1,dv2,dv3
 */
bool CommandLineParser::getDoubleVectorValue(const std::string& parName, const std::string& delimeter, std::vector<double> &dblVectorPar)
{
    bool parExist = false;
    std::multimap<std::string, std::string>::iterator it = m_parValue.find(parName);
    if (it != m_parValue.end()) {
       std::string strValue = it->second;
       std::vector<std::string> parsVector = convertParStringToVector(strValue, delimeter);
       int size = parsVector.size();
       if (size > 0) {
          dblVectorPar.resize(0);
          parExist = true;
          for (int i = 0; i < size; i++) {
              double itemValue = 0.0;
              if (num_valid<double>(parsVector.at(i).c_str(), itemValue) ) {
                 dblVectorPar.push_back(itemValue);
              } else {
                 dblVectorPar.clear();
                 fprintf(stderr, "%s: Failed to convert item %s to a double value\n", parName.c_str(), parsVector.at(i).c_str());
                 parExist = false;
                 break;
              }
          }
       } else {
          std::cerr << "\n";
          std::cerr << "Error: during parameter parse stage, an error happens:" << "\n";
          std::cerr << "Error: input parameter number of " << parName << " is 0." << "\n";
       }
    }
    return parExist;
}

/*
 * Usage:
 * fltPar=fv1,fv2,fv3
 * However require at least num of integer values
 */
bool CommandLineParser::getDoubleVectorValue(const std::string& parName, const std::string& delimeter, const int& num,
                                                    std::vector<double> &dblVectorPar)
{
    bool vectorValid = false;
    bool parExist = getDoubleVectorValue(parName, delimeter, dblVectorPar);
    if (parExist) {
       int size = dblVectorPar.size();
       vectorValid = true;
       if (size < num) {
          std::cerr << "\n";
          std::cerr << "Error: during parameter parse stage, an error happens:" << "\n";
          std::cerr << "Error: input parameter number of " << parName << " is less " << num << "\n";
          vectorValid = false;
     #if 0
       } else if (size > num) {
          dblVectorPar.resize(num);
     #endif
       }
    }
    return vectorValid;
}

void CommandLineParser::setParameterValuePair(const int argc, char **argv)
{
    if (argc > 1) {
       m_parValue.clear();
       for  (int i = 1; i < argc; i++) {
            std::string parValueString = std::string(argv[i]);
            std::size_t firstEqPosition = parValueString.find("=");
            if (firstEqPosition != std::string::npos) {
               std::size_t secondEqPosition = parValueString.find("=", firstEqPosition+1);
               // Only one "=" is allowed
               if (secondEqPosition == std::string::npos) {
                  std::string parName = std::string(parValueString.begin(), parValueString.begin()+firstEqPosition);
                  std::string parValue = std::string(parValueString.begin()+firstEqPosition+1, parValueString.end());
                  if (!parName.empty() && !parValue.empty()) {
                     //std::cout << "parName: " << parName <<"; parValue: " << parValue << std::endl;
                     m_parValue.insert(std::pair<std::string, std::string>(parName, parValue));
                  }
               }
            }
       }
    }
}

std::vector<std::string> CommandLineParser::convertParStringToVector(const std::string& sourceString, const std::string& delimeter)
{
    std::vector <std::string> parsVector;
    boost::split(parsVector, sourceString, boost::is_any_of(delimeter));
    return parsVector;
}


