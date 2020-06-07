/*
 * Copyright (C) 2020 Jan Marvin Garbuszus
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <Rcpp.h>
#include <cstdio>
#include <stdint.h>
#include <string>
#include <fstream>
#include <streambuf>
#include <regex>

#include "sas.h"
#include "uncompress.h"
using namespace Rcpp;

//' Writes the binary SAS file
//'
//' @param filePath The full systempath to the sas7bdat file you want to export.
//' @param dat an R-Object of class data.frame.
//' @param compress the file
//' @param debug print debug information
//' @export
// [[Rcpp::export]]
void writesas(const char * filePath, Rcpp::DataFrame dat, uint8_t compress,
              bool debug) {

  uint32_t k = dat.size();
  uint64_t n = dat.nrows();


  std::fstream sas (filePath, std::ios::out | std::ios::binary);
  if (sas.is_open())
  {

    bool swapit = 0;
    int8_t ALIGN_2_VALUE = 0;
    int8_t u64 = 0;
    int8_t zero = 0;
    int16_t PAGE_TYPE = 0, BLOCK_COUNT = 0, SUBHEADER_COUNT = 0;

    // partially known: value is known meaning is unknown
    int8_t
      pkt1 = 0, pkt2 = 0, pkt3 = 0, pkt4 = 0;

    // unknown
    int8_t unk8 = 0;
    int16_t unk16 = 0;
    int32_t unk32 = 0;
    uint32_t uunk32 = 0;
    int64_t unk1 = 0, unk2 = 0, unk3 = 0;
    double unkdub = 0;

    double created = 0, created2 = 0;   // 8
    double modified = 0, modified2 = 0; // 16
    double thrdts = 0;

    // possibly make headersize and pagesize variable
    int32_t headersize = 65536;
    int32_t pagesize = 65536;
    int64_t pagecount = 0;
    uint32_t pageseqnum32 = 0;

    // Begin: Magic Number
    int32_t
      mn1 = 0, mn2 = 0, mn3 = 0, mn4 = 1619126978, mn5 = -820964173,
        mn6 = 561853, mn7 = -1942894839, mn8 = 286269208;

    writebin(mn1, sas, 0);
    writebin(mn2, sas, 0);
    writebin(mn3, sas, 0);
    writebin(mn4, sas, 0);
    writebin(mn5, sas, 0);
    writebin(mn6, sas, 0);
    writebin(mn7, sas, 0);
    writebin(mn8, sas, 0);
    // End: Magic Number

    // packet 1 --------------------------------------- //
    int8_t ALIGN_1_CHECKER_VALUE = 51;
    if (ALIGN_1_CHECKER_VALUE == 51) u64 = 4;
    int8_t U64_BYTE_CHECKER_VALUE = 51;
    if (U64_BYTE_CHECKER_VALUE == 51) ALIGN_2_VALUE = 4;

    writebin(ALIGN_1_CHECKER_VALUE, sas, 0);
    writebin(unk8, sas, 0); // 34
    writebin(unk8, sas, 0); // 0
    writebin(U64_BYTE_CHECKER_VALUE, sas, 0);

    // packet 2 --------------------------------------- //
    int8_t ENDIANNESS = 1; // 0 is swapit = 1
    uint8_t PLATFORM = 49; // (1) 49 Unix (2) 50 Win

    writebin(unk8, sas, 0); // 51
    writebin(ENDIANNESS, sas, 0);
    writebin(unk8, sas, 0); // 2
    writebin(PLATFORM, sas, 0);

    // packet 3 --------------------------------------- //
    writebin(unk8, sas, 0); // 1
    writebin(unk8, sas, 0);
    writebin(unk8, sas, 0);
    writebin(unk8, sas, 0);

    // packet 5 --------------------------------------- //
    writebin(unk8, sas, 0);
    writebin(unk8, sas, 0);
    writebin(unk8, sas, 0);
    writebin(unk8, sas, 0); // 20

    // packet 5 --------------------------------------- //
    pkt1 = 0, pkt2 = 0, pkt3 = 3, pkt4 = 1;
    writebin(pkt1, sas, 0);
    writebin(pkt2, sas, 0);
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0);

    // packet 6 --------------------------------------- //
    writebin(unk8, sas, 0); // 24
    writebin(unk8, sas, 0); // 31
    writebin(unk8, sas, 0); // 16
    writebin(unk8, sas, 0); // 17

    // packet 7 --------------------------------------- //
    int8_t U64_BYTE_CHECKER_VALUE2 = 51;

    writebin(ALIGN_1_CHECKER_VALUE, sas, 0);
    writebin(unk8, sas, 0); // 34
    writebin(unk8, sas, 0); // 0
    writebin(U64_BYTE_CHECKER_VALUE2, sas, 0);

    // packet 8 --------------------------------------- //
    writebin(unk8, sas, 0); // 51
    writebin(ENDIANNESS, sas, 0);
    writebin(unk8, sas, 0); // 2
    writebin(PLATFORM, sas, 0);

    // packet 9 --------------------------------------- //
    pkt1 = 1, pkt2 = 51, pkt3 = 1, pkt4 = 35;
    writebin(pkt1, sas, 0); // 1 | 4
    writebin(pkt2, sas, 0);
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0);

    // packet 10 -------------------------------------- //
    int8_t encoding = 20; // utf8
    writebin(unk8, sas, 0); // 51
    writebin(unk8, sas, 0); // 0
    writebin(encoding, sas, 0);
    writebin(unk8, sas, 0); // encoding again?

    // packet 11 --------------------------------------- //
    pkt1 = 0, pkt2 = 32, pkt3 = 3, pkt4 = 1;
    writebin(pkt1, sas, 0);
    writebin(pkt2, sas, 0); // 16 | 32 (int16 or 32?)
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0);

    writebin(unk32, sas, 0);
    writebin(unk32, sas, 0);

    std::string sasfile = "SAS FILE";
    writestr(sasfile, sasfile.size(), sas);

    std::string dataset = "test";
    writestr(dataset, 64, sas);

    std::string filetype = "DATA";
    writestr(filetype, 8, sas);

    if (ALIGN_2_VALUE == 4) {
      writebin(unk32, sas, 0);
    }

    writebin(created, sas, 0);
    writebin(created2, sas, 0);
    writebin(modified, sas, 0);
    writebin(modified2, sas, 0);

    writebin(headersize, sas, 0);
    writebin(pagesize, sas, 0);
    if (u64 == 4) {
      writebin(pagecount, sas, swapit);
    } else {
      writebin((int32_t)pagecount, sas, swapit);
    }

    writebin(unkdub, sas, 0);

    std::string sasrel = "9.0401M6";
    std::string sasserv = "Linux";
    std::string osver = "5.6.15-arch1-1";
    std::string osmaker = "";
    std::string osname = "x86_64";

    writestr(sasrel, 8, sas);
    writestr(sasserv, 16, sas);
    writestr(osver, 16, sas);
    writestr(osmaker, 16, sas);
    writestr(osname, 16, sas);

    // large unk
    writebin(uunk32, sas, 0);

    // three identical smaller unks
    writebin(uunk32, sas, 0);
    writebin(uunk32, sas, 0);
    writebin(uunk32, sas, 0);

    // 2 * zero
    writebin(unkdub, sas, 0);
    writebin(unkdub, sas, 0);

    writebin(pageseqnum32, sas, 0);
    writebin(unk32, sas, 0);

    writebin(thrdts, sas, 0);

    uint64_t num_zeros = headersize - sas.tellg();

    for (uint64_t i = 0; i < num_zeros; ++i) {
      writebin(zero, sas, swapit);
    }
    // end of header -------------------------------------------------------- //

    /*** write pages **********************************************************/
    for (auto pg = 0; pg < pagecount; ++pg) {
      checkUserInterrupt();

      // Page Offset Table
      if (u64 == 4) {
        writebin(pageseqnum32, sas, swapit);
        writebin(unk32, sas, swapit);
        writebin(unk1, sas, swapit);
        writebin(unk2, sas, swapit);
        writebin(unk3, sas, swapit);
      } else {
        writebin(pageseqnum32, sas, swapit);
        writebin((int32_t)unk1, sas, swapit);
        writebin((int32_t)unk2, sas, swapit);
        writebin((int32_t)unk3, sas, swapit);
      }


      // all 3 values need to be defined depending on the input
      writebin(PAGE_TYPE, sas, swapit);
      writebin(BLOCK_COUNT, sas, swapit);
      writebin(SUBHEADER_COUNT, sas, swapit);
      writebin(unk16, sas, swapit);

      // not really required, but maybe usefull to track positions?
      std::vector<PO_Tab> potabs(SUBHEADER_COUNT);

      // write entries for all subheaders
      // write them as dummies first, once their position is known fill sh_off.
      //  requires to track the position of each
      for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
        if (u64 == 4) {

          writebin(potabs[i].SH_OFF, sas, swapit);
          writebin(potabs[i].SH_LEN, sas, swapit);
          writebin(potabs[i].COMPRESSION, sas, swapit);
          writebin(potabs[i].SH_TYPE, sas, swapit);

          writebin(zero, sas, swapit);
          writebin(zero, sas, swapit);
          writebin(zero, sas, swapit);

        } else {

          writebin((uint32_t)potabs[i].SH_OFF, sas, swapit);
          writebin((uint32_t)potabs[i].SH_LEN, sas, swapit);
          writebin(potabs[i].COMPRESSION, sas, swapit);
          writebin(potabs[i].SH_TYPE, sas, swapit);

          writebin(zero, sas, swapit);
        }

      }

    }
  }
}
