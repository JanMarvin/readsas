#ifndef SAS_H
#define SAS_H

#include <Rcpp.h>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/regex.hpp>

#include "swap_endian.h"


template <typename T>
T readbin( T t , std::istream& sas, bool swapit)
{
  if (!sas.read ((char*)&t, sizeof(t)))
    Rcpp::stop("readbin: a binary read error occurred");
  if (swapit==0)
    return(t);
  else
    return(swap_endian(t));
}

double readbinlen(double d, std::istream& sas, bool swapit, int len)
{

  unsigned char buffer[sizeof d] = {0};

  for (int i = 0; i < len; ++i) {

    unsigned char tmp;
    if(!sas.read ((char*)&tmp, sizeof(tmp)))
      Rcpp::stop("readbin: a binary read error occurred");

    int pos = 8 - (len-i);
    buffer[pos] = tmp;

  }

  memcpy(&d , buffer, sizeof(d));

  if (swapit==0)
    return(d);
  else
    return(swap_endian(d));
}


template <typename T>
inline std::string readstring(std::string &mystring, T& sas)
{

  if (!sas.read(&mystring[0], mystring.size()))
    Rcpp::stop("char: a binary read error occurred");

  return(mystring);
}


// std::string decode ( const std::string & to_decode ) {
//   boost::regex e ( "(\\d+)(\\w)" ) ;
//   boost::match_results<std::string::const_iterator> matches ;
//   std::ostringstream oss ;
//   std::string::const_iterator start = to_decode.begin( ) , end = to_decode.end( ) ;
//   while ( boost::regex_search ( start , end , matches , e ) ) {
//     std::string numberstring ( matches[ 1 ].first , matches[ 1 ].second ) ;
//     int number = atoi( numberstring.c_str( ) ) ;
//     std::string character ( matches[ 2 ].first , matches[ 2 ].second ) ;
//     for ( int i = 0 ; i < number ; i++ )
//       oss << character ;
//     start = matches[ 2 ].second ;
//   }
//   return oss.str( ) ;
// }

// PAGE_OFFSET_TABLE
struct PO_Tab {
  int64_t SH_OFF = 0;
  int64_t SH_LEN = 0;
  int8_t COMPRESSION = 0;
  int8_t SH_TYPE = 0;
};

// COLUMN_NAME_POINTER
struct CN_Poi {
  int16_t CN_IDX = 0;
  int16_t CN_OFF = 0;
  int16_t CN_LEN = 0;
  int16_t zeros  = 0;
};

struct CN_Att {
  int64_t CN_OFF  = 0;
  int32_t CN_WID  = 0;
  int16_t NM_FLAG = 0;
  int8_t  CN_TYP  = 0;
  int8_t  UNK8    = 0;
};


struct SCV {
  int64_t SIG = 0;
  int64_t FIRST = 0;
  int16_t F_POS = 0; // why is this pos not int64_?
  int64_t LAST = 0;
  int16_t L_POS = 0;
};


inline std::string int32_to_hex (int32_t val) {
  std::stringstream stream;
  stream << std::hex << val;
  std::string res( stream.str() );
  return (res);
}



#endif
