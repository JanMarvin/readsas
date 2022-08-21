#ifndef SAS_H
#define SAS_H

#include <Rcpp.h>
#include <fstream>
#include <string>
#include <sstream>

#include "swap_endian.h"

inline void writestr(std::string val_s, int32_t len, std::fstream& sas)
{

  std::stringstream val_stream;
  val_stream << std::left << std::setw(len) << std::setfill(' ') << val_s;
  std::string val_strl = val_stream.str();

  sas.write(val_strl.c_str(), val_strl.length());

}

// return only the matched positions. Either Rcpps in() can't handle Character-
// Vectors or I could not make it work. Wanted to select the selected varname
// position from the varnames vector.
inline Rcpp::IntegerVector choose(Rcpp::CharacterVector x,
                                  Rcpp::CharacterVector y)
{
  Rcpp::IntegerVector mm = Rcpp::match(x, y);

  if (Rcpp::any(Rcpp::is_na(mm))) {
    Rcpp::LogicalVector ll = !Rcpp::is_na(mm);

    Rcpp::CharacterVector ms = x[ll==0];

    // does not work if ms contains multiple names: Rcpp::as<std::string>(ms)
    Rcpp::Rcout << "Variable " << ms <<
      " was not found in sas-file." << std::endl;
  }

  // report position for found cases
  mm = Rcpp::match(y, x);

  return(mm);
}

inline Rcpp::IntegerVector order(Rcpp::IntegerVector x) {
  Rcpp::IntegerVector sorted = Rcpp::clone(x).sort();
  return match(sorted, x) -1;
}

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

// PAGE_OFFSET_TABLE
struct PO_Tab {
  uint64_t SH_OFF = 0;
  uint64_t SH_LEN = 0;
  uint8_t COMPRESSION = 0;
  uint8_t SH_TYPE = 0;
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

struct idxofflen {
  int16_t IDX = 0;
  int16_t OFF = 0;
  int16_t LEN = 0;
};


template <typename T>
inline std::string int_to_hex (T val) {
  std::stringstream stream;
  stream << std::hex << val;
  std::string res( stream.str() );
  return (res);
}

inline std::string SASEncoding(uint8_t encval) {

  std::string enc = "";

  switch(encval)
  {
  case 70:
    enc = "x-MacArabic";
    break;
  case 245:
    enc = "x-MacCroatian";
    break;
  case 246:
    enc = "x-MacCyrillic";
    break;
  case 72:
    enc = "x-MacGreek";
    break;
  case 71:
    enc = "x-MacHebrew";
    break;
  case 163:
    enc = "x-MacIceland";
    break;
  case 34:
    enc = "ISO-8859-6";
    break;
  case 69:
    enc = "x-MacRoman";
    break;
  case 247:
    enc = "x-MacRomania";
    break;
  case 73:
    enc = "x-MacThai";
    break;
  case 75:
    enc = "x-MacTurkish";
    break;
  case 76:
    enc = "x-MacUkraine";
    break;
  case 123:
    enc = "Big5";
    break;
  case 33:
    enc = "ISO-8859-5";
    break;
  case 78:
    enc = "IBM037";
    break;
  case 95:
    enc = "x-IBM1025";
    break;
  case 207:
    enc = "x-IBM1097";
    break;
  case 98:
    enc = "x-IBM1112";
    break;
  case 99:
    enc = "x-IBM1122";
    break;
  case 183:
    enc = "IBM01140";
    break;
  case 184:
    enc = "IBM01141";
    break;
  case 185:
    enc = "IBM01142";
    break;
  case 186:
    enc = "IBM01143";
    break;
  case 187:
    enc = "IBM01144";
    break;
  case 188:
    enc = "IBM01145";
    break;
  case 189:
    enc = "IBM01146";
    break;
  case 190:
    enc = "IBM01147";
    break;
  case 191:
    enc = "IBM01148";
    break;
  case 211:
    enc = "IBM01149";
    break;
  case 87:
    enc = "IBM424";
    break;
  case 88:
    enc = "IBM500";
    break;
  case 89:
    enc = "IBM-Thai";
    break;
  case 90:
    enc = "IBM870";
    break;
  case 91:
    enc = "x-IBM875";
    break;
  case 125:
    enc = "GBK";
    break;
  case 134:
    enc = "EUC-JP";
    break;
  case 140:
    enc = "EUC-KR";
    break;
  case 119:
    enc = "x-EUC-TW";
    break;
  case 205:
    enc = "GB18030";
    break;
  case 35:
    enc = "ISO-8859-7";
    break;
  case 36:
    enc = "ISO-8859-8";
    break;
  case 128:
    enc = "x-IBM1381";
    break;
  case 130:
    enc = "x-IBM930";
    break;
  case 139:
    enc = "x-IBM933";
    break;
  case 124:
    enc = "x-IBM935";
    break;
  case 117:
    enc = "x-IBM937";
    break;
  case 129:
    enc = "x-IBM939";
    break;
  case 137:
    enc = "x-IBM942";
    break;
  case 142:
    enc = "x-IBM949";
    break;
  case 172:
    enc = "x-ISO2022-CN-CNS";
    break;
  case 169:
    enc = "x-ISO2022-CN-GB";
    break;
  case 167:
    enc = "ISO-2022-JP";
    break;
  case 168:
    enc = "ISO-2022-KR";
    break;
  case 29:
    enc = "ISO-8859-1";
    break;
  case 30:
    enc = "ISO-8859-2";
    break;
  case 31:
    enc = "ISO-8859-3";
    break;
  case 32:
    enc = "ISO-8859-4";
    break;
  case 37:
    enc = "ISO-8859-9";
    break;
  case 242:
    enc = "ISO-8859-13";
    break;
  case 40:
    enc = "ISO-8859-15";
    break;
  case 136:
    enc = "x-windows-iso2022jp";
    break;
  case 126:
    enc = "x-mswin-936";
    break;
  case 141:
    enc = "x-windows-949";
    break;
  case 118:
    enc = "x-windows-950";
    break;
  case 173:
    enc = "IBM037";
    break;
  case 108:
    enc = "x-IBM1025";
    break;
  case 109:
    enc = "IBM1026";
    break;
  case 110:
    enc = "IBM1047";
    break;
  case 208:
    enc = "x-IBM1097";
    break;
  case 111:
    enc = "x-IBM1112";
    break;
  case 112:
    enc = "x-IBM1122";
    break;
  case 192:
    enc = "IBM01140";
    break;
  case 193:
    enc = "IBM01141";
    break;
  case 194:
    enc = "IBM01142";
    break;
  case 195:
    enc = "IBM01143";
    break;
  case 196:
    enc = "IBM01144";
    break;
  case 197:
    enc = "IBM01145";
    break;
  case 198:
    enc = "IBM01146";
    break;
  case 199:
    enc = "IBM01147";
    break;
  case 200:
    enc = "IBM01148";
    break;
  case 212:
    enc = "IBM01149";
    break;
  case 102:
    enc = "IBM424";
    break;
  case 103:
    enc = "IBM-Thai";
    break;
  case 104:
    enc = "IBM870";
    break;
  case 105:
    enc = "x-IBM875";
    break;
  case 234:
    enc = "x-IBM930";
    break;
  case 235:
    enc = "x-IBM933";
    break;
  case 236:
    enc = "x-IBM935";
    break;
  case 237:
    enc = "x-IBM937";
    break;
  case 238:
    enc = "x-IBM939";
    break;
  case 43:
    enc = "IBM437";
    break;
  case 44:
    enc = "IBM850";
    break;
  case 45:
    enc = "IBM852";
    break;
  case 58:
    enc = "IBM857";
    break;
  case 46:
    enc = "IBM00858";
    break;
  case 47:
    enc = "IBM862";
    break;
  case 51:
    enc = "IBM866";
    break;
  case 138:
    enc = "Shift_JIS";
    break;
  case 248:
    enc = "JIS_X0201";
    break;
  case 39:
    enc = "x-iso-8859-11";
    break;
  case 28:
    enc = "US-ASCII";
    break;
  case 20:
    enc = "UTF-8";
    break;
  case 66:
    enc = "windows-1256";
    break;
  case 67:
    enc = "windows-1257";
    break;
  case 61:
    enc = "windows-1251";
    break;
  case 63:
    enc = "windows-1253";
    break;
  case 65:
    enc = "windows-1255";
    break;
  case 62:
    enc = "windows-1252";
    break;
  case 60:
    enc = "windows-1250";
    break;
  case 64:
    enc = "windows-1254";
    break;
  case 68:
    enc = "windows-1258";
    break;
  }

  return enc;
}

// order only the valid options
Rcpp::IntegerVector order_(Rcpp::IntegerVector v) {
  Rcpp::IntegerVector z = clone(v);
  if (any(Rcpp::is_na(z))) {
    return z[z >= 0] = order(z[z >= 0]);
  } else{
    // Rcpp::Rcout << "no missings" << order(z) << std::endl;
    return order(z);
  }
}

// create new column vector. has the expected length, but counts new
Rcpp::IntegerVector cvec_(Rcpp::LogicalVector &l) {

  Rcpp::IntegerVector out(l.size());

  auto idx = 0;
  for (auto i = 0; i < l.size(); ++i) {
    if (l[i] == 0) {
      out[i] = idx;
      ++idx;
    } else
      out[i] = -1;
  }

  return out;
}


#endif
