/*
 * Both functions are c++ conversions of java functions from the parso library
 * Licensed under the Apache License, Version 2.0. Copyright 2015 EPAM
 *
 */

#include <Rcpp.h>
#include <math.h>
#include "uncompress.h"

std::string SASYZCRL(uint64_t rowlen, uint64_t reslen, std::string rowstr,
                     bool debug)
{

  // Minimal change for VLA: Use new[] for dynamic allocation
  size_t rowstr_size = rowstr.size();
  char* row = new char[rowstr_size];
  memcpy(row, rowstr.data(), rowstr_size); // Copy data to the allocated memory


  std::string res = "";
  uint32_t rowoff = 0;

  while (rowoff < rowlen) {

    // unsure
    int32_t cbyte = row[rowoff] & 0xF0;
    int32_t ebyte = row[rowoff] & 0x0F;

    int32_t len = 0;

    switch (cbyte) {

    case 0:
    case 16:
    case 32:
    case 48:
    {
      if (rowoff != rowlen - 1) {

      len = (row[rowoff + 1] & 0xFF) + 64 + row[rowoff] * 256;

      res +=  rowstr.substr(rowoff+2, len);

      rowoff += len + 1;
    }

      break;
    }

    case 64:
    {
      int32_t copyCounter = ebyte * 16 + (row[rowoff + 1] & 0xFF);

      for (auto i = 0; i < copyCounter + 18; ++i) {
        res += row[rowoff + 2];
      }

      rowoff += 2;

      break;
    }

    case 96:
    {
      uint8_t val = 32;

      for (auto i = 0; i < ebyte * 256 + (row[rowoff + 1] & 0xFF) + 17; ++i) {
        res += val;
      }

      ++rowoff;

      break;
    }

    case 112:
    {
      uint8_t val = 0;

      for (auto i = 0; i < ebyte * 256 + (row[rowoff + 1] & 0xFF) + 17; ++i) {
        res += val;
      }

      ++rowoff;

      break;
    }

    case 128:
    case 144:
    case 160:
    case 176:
    {
      int32_t p0 = ((ebyte + 1) + (cbyte - 128));
      int32_t p1 = (rowlen - (rowoff + 1));
      len = std::min(p0 , p1);

      res += rowstr.substr(rowoff+1, len);
      rowoff += len;

      break;
    }

    case 192:
    {
      for (auto i = 0; i < ebyte + 3; ++i) {
      res += row[rowoff + 1];
    }
      ++rowoff;

      break;
    }

    case 208:
    {
      uint8_t val = 64;
      for (auto i = 0; i < ebyte + 2; ++i) {
        res += val;
      }

      break;
    }

    case 224:
    {
      uint8_t val = 32;
      for (auto i = 0; i < ebyte + 2; ++i) {
        res += val;
      }

      break;
    }

    case 240:
    {
      uint8_t val = 0;
      for (auto i = 0; i < ebyte + 2; ++i) {
        res += val;
      }

      break;
    }

    default:
    {
      Rcpp::warning("SASYZCRL: Error control byte: %d", cbyte);
      break; // Exit loop on error
    }
    }
    ++rowoff;

  }

  // comp_deleted throws warning
  if (res.size() != reslen)
    Rcpp::warning("SASYZCRL: string size %zu != %zu\n",
                  res.size(), reslen);

  // Clean up allocated memory
  delete[] row;

  return res;
}

std::string SASYZCR2(uint64_t rowlen, uint64_t reslen, const std::string rowstr,
                     bool debug) {

  if (debug)
    Rcpp::Rcout << "row ----------------------------- " << std::endl;

  // Minimal change for VLA: Use new[] for dynamic allocation
  size_t rowstr_size = rowstr.size();
  char* row = new char[rowstr_size];
  memcpy(row, rowstr.data(), rowstr_size); // Copy data to the allocated memory


  std::string res = "";
  uint32_t rowoff = 0;
  int32_t resoff = 0, ofs = 0, cbit = 0, cmsk = 0;

  while (rowoff < rowlen) {

    cmsk >>= 1;
    if (cmsk == 0) {
      cbit = (((row[rowoff]) & 0xff) << 8) | (row[rowoff + 1] & 0xff);
      rowoff += 2;
      cmsk = 0x8000;
    }

    // copy char if control bit is zero
    if ((cbit & cmsk) == 0) {
      res += row[rowoff++];
      ++resoff;

      continue;
    }

    int32_t cmd = (row[rowoff] >> 4) & 0x0F;
    int32_t len = row[rowoff++] & 0x0F;

    switch (cmd)
    {
    case 0: // short RLE
    {
      len += 3;
      for (int32_t i = 0; i < len; ++i) {
        res += row[rowoff];
      }
      ++rowoff;
      resoff += len;

      // if (debug)
        // Rcpp::Rcout << "d1 " << std::endl;

      break;
    }

    case 1: // long RLE
    {
      len += ((row[rowoff++] & 0xff) << 4) + 19;
      for (int32_t i = 0; i < len; ++i) {
        res += row[rowoff];
      }
      ++rowoff;
      resoff += len;

      // if (debug)
        // Rcpp::Rcout << "d2 " << std::endl;

      break;
    }

    case 2: // long
    {
      ofs = len + 3;
      ofs += ((row[rowoff++] & 0xff) << 4);
      len = (row[rowoff++] & 0xff) + 16;

      auto pos = resoff - ofs;
      if (pos< 0)
        pos = std::abs(pos);

      // if (debug)
        // Rcpp::Rcout << "d3 " << resoff << " " << ofs << " " << pos << std::endl;

      res += res.substr(pos, len);
      resoff += len;

      break;
    }

    default: // short
    {
      ofs = len + 3;
      ofs += ((row[rowoff++] & 0xff) << 4);

      auto pos = resoff - ofs;
      if (pos< 0)
        pos = std::abs(pos);

      // if (debug)
        // Rcpp::Rcout << "d4 " << resoff << " " << ofs << " " << pos << std::endl;

      res += res.substr(pos, cmd);
      resoff += cmd;

      break;
    }
    }
  }

  if (res.size() != reslen)
    Rcpp::warning("SASYZCR2: string size %zu != %zu\n", // Use %zu for size_t
                  res.size(), reslen);

  // Clean up allocated memory
  delete[] row;

  return res;
}
