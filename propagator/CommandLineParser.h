/* 
 * File:   CommandLineParser.h
 * Author: Yezheng (Jim) HU
 * Date: June 23 2017
 */

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H
/**
 * This class receive command line parameters and sparse it as user defined 
 * parameter names.
 * 
 * If multiple values are expected from the command line. A vector <type> will 
 * be returned which stores the corresponding values.
 * 
 */
#include <string>
#include <vector>
#include <map>

#include <boost/lexical_cast.hpp>

//using namespace boost;

class CommandLineParser
{
  public:
    CommandLineParser(int num, char** argv);

    void printParameters();

    bool ParameterNameExists(const std::string& parName);

    /* Usage: strPar=strValue */
    bool getStringValue(const std::string& parName, std::string& par);

    /*
     * Usage: strPar=strValue1,strValue2,strValue3,... strPar=strValuea,strValueb,strValuec,...
     * Example:
     * infile=file1,file2,file3 infile=file_a,file_b,file_c
     */
    std::vector<std::string> getStringVectorValue(const std::string& parName, const std::string& delimeter=",", bool required=false);

    /* As above but equires at least num of string values */
    bool getStringVectorValue(const std::string& parName,
           const std::string& delimeter, const int& num, std::vector<std::string>& stringTagVector);

    /* Usage: intPar=intValue */
    bool getIntValue(const std::string& parName, int& parValue);

    /* Usage: intPar=iv1,iv2,iv3 */
    bool getIntVectorValue(const std::string& parName, const std::string& delimeter, std::vector<int> &intVectorPar);

    /*
     * Usage: intPar=iv1,iv2,iv3
     * However requires at least num of integer values
     */
    bool getIntVectorValue(const std::string& parName, const std::string& delimeter, const int& num,
                              std::vector<int>& intVectorPar);

    /* Usage: fltPar=fltValue */
    bool getFloatValue(const std::string& parName, float& parValue);

    /* Usage: fltPar=fv1,fv2,fv3 */
    bool getFloatVectorValue(const std::string& parName, const std::string& delimeter, const int& num, std::vector<float> &dblVectorPar);

    /*
     * Usage: fltPar=fv1,fv2,fv3
     * However requires at least num of integer values
     */
    bool getFloatVectorValue(const std::string& parName, const std::string& delimeter, std::vector<float> &dblVectorPar);

    /* Usage: dblPar=dblValue */
    bool getDoubleValue(const std::string& parName, double& parValue);

    /* Usage: dblPar=dv1,dv2,dv3 */
    bool getDoubleVectorValue(const std::string& parName, const std::string& delimeter, const int& num, std::vector<double> &dblVectorPar);
    /*
     * Usage: dblPar=dv1,dv2,dv3
     * However requires at least num of integer values
     */
    bool getDoubleVectorValue(const std::string& parName, const std::string& delimeter, std::vector<double> &dblVectorPar);

  private:
    std::multimap<std::string, std::string> m_parValue;

    void setParameterValuePair(const int argc, char **argv);

    /** convert parameter elements which are separated by a specific delimeter in command line to vectors.*/
    std::vector<std::string> convertParStringToVector(const std::string& source, const std::string& delimeter);
};

template<typename T>
bool num_valid(const char *str, T& strValue)
{
   try
   {
       strValue = boost::lexical_cast<T>(str);
       return true;
   }
   catch(boost::bad_lexical_cast&)
   {
      return false;
   };
}

#endif	/* COMMANDLINEPARSER_H */

