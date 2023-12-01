/*
 * Copyright (C) 2020-2021, 2023 Jan Marvin Garbuszus
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
//' @param bit32 write smaller 32 bit file
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
void writesas(const char * filePath, Rcpp::DataFrame dat, uint8_t compress,
              bool debug, bool bit32, int32_t headersize, int32_t pagesize) {

  uint32_t k = dat.size();
  uint64_t n = dat.nrows();
  CharacterVector nvarnames = dat.attr("names");
  CharacterVector nlabels = dat.attr("labels");
  CharacterVector nformats = dat.attr("formats");
  IntegerVector width = dat.attr("width");
  IntegerVector decim = dat.attr("decim");
  IntegerVector vartypes = dat.attr("vartypes");
  IntegerVector colwidth = dat.attr("colwidth");
  CharacterVector nvarlabels = dat.attr("varlabels");

  auto pagecount_pos1 = 0, pagecount_pos2 = 0, rcmix_pos = 0,
    block_count_pos1 = 0, rowcount_pos = 0;
  int64_t pagecount = 1;

  // trinity of varnames varlabels and varformats
  std::vector<std::string> varnames (k);
  std::vector<std::string> varlabels (k);
  std::vector<std::string> varformats (k);

  auto totalvarnamesize = 0;
  auto totalvarlabelssize = 0;
  auto totalvarformatssize = 0;

  for (auto i = 0; i < k; ++i) {
    // varnames
    std::string varname = as<std::string>(nvarnames[i]);
    if (varname.size() % 4 > 0 && varname.size() <= 32) { // for < 4, why not < 32?
      int8_t div4len = 4 * ceil((double)varname.size()/4);
      varname.resize(div4len, '\0');
      // for now resize everything to 8 char
    } else if (varname.size() > 32) {
      warning("varname was shorten to 32 characters");
      varname.resize(32, '\0');
    }

    varnames[i] = varname;
    totalvarnamesize += varname.size();

    // varlabels
    std::string varlabel = as<std::string>(nvarlabels[i]);
    if (varlabel.size() % 4 > 0 && varlabel.size() <= 32) { // for < 4, why not < 32?
      int8_t div4len = 4 * ceil((double)varlabel.size()/4);
      varlabel.resize(div4len, '\0');
      // for now resize everything to 8 char
    } else if (varlabel.size() > 32) {
      warning("varlabel was shorten to 32 characters");
      varlabel.resize(32, '\0');
    }

    varlabels[i] = varlabel;
    totalvarlabelssize += varlabel.size();

    // varformats
    std::string varformat = as<std::string>(nformats[i]);
    if (debug) Rcout << varformat << " : " << varformat.size() << std::endl;

    if (varformat.size() == 0) {
      // Rcout << "size 0" << std::endl;
    }
    if (varformat.size() > 0 && varformat.size() <= 4) {
      varformat.resize(4, '\0');
    } else if (varformat.size() > 4) { // resize everything to 8 characters
      varformat.resize(8, '\0');
    }

    varformats[i] = varformat;
    totalvarformatssize += varformat.size();

    // Rcout << "varformat size: " << varformat.size() << std::endl;

  }



  std::fstream sas (filePath, std::ios::out | std::ios::binary);
  if (sas.is_open())
  {

    auto pos_SH_C = 0;

    bool swapit = 0;
    int8_t ALIGN_2_VALUE = 0;
    int8_t u64 = 0;
    int8_t zero8 = 0;
    int16_t zero16 = 0;
    int16_t PAGE_TYPE = 512, SUBHEADER_COUNT = (7 + k);
    if (k > 1) SUBHEADER_COUNT++;
    int16_t BLOCK_COUNT = n + SUBHEADER_COUNT;

    // partially known: value is known meaning is unknown
    int8_t
    pkt1 = 0, pkt2 = 0, pkt3 = 0, pkt4 = 0;

    // unknown
    int8_t unk8 = 0;
    int16_t unk16 = 0;
    int32_t unk32 = 0;
    uint32_t uunk32 = 0;
    int64_t unk64 = 0;
    int64_t unk1 = 0, unk2 = 0, unk3 = 0;
    double unkdub = 0;
    auto inioff = 0;

    double created = 0, created2 = 0;   // 8
    double modified = 0, modified2 = 0; // 16
    double thrdts = 0;

    // possibly make headersize and pagesize variable
    // int32_t headersize = 65536;
    // if (bit32 == 1) headersize = 1024;
    // int32_t pagesize = 65536;
    // if (bit32 == 1) pagesize = 8192;
    // int64_t pagecount = 1;
    int32_t pgidx = 0;
    int32_t pageseqnum32 = pgidx;

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
    if (bit32 == 1) ALIGN_1_CHECKER_VALUE = 34;
    if (ALIGN_1_CHECKER_VALUE == 51) u64 = 4;
    int8_t U64_BYTE_CHECKER_VALUE = 51;
    if (bit32 == 1) U64_BYTE_CHECKER_VALUE = 50;
    if (U64_BYTE_CHECKER_VALUE == 51) ALIGN_2_VALUE = 4;

    pkt2 = 34, pkt3 = 0;

    writebin(ALIGN_1_CHECKER_VALUE, sas, 0);
    writebin(pkt2, sas, 0); // 34
    writebin(pkt3, sas, 0); // 0
    writebin(U64_BYTE_CHECKER_VALUE, sas, 0);

    // packet 2 --------------------------------------- //
    int8_t ENDIANNESS = 1; // 0 is swapit = 1
    uint8_t PLATFORM = 49; // (1) 49 Unix (2) 50 Win

    pkt1 = 51, pkt3 = 2;
    if (bit32 == 1) pkt1 = 50;
    writebin(pkt1, sas, 0); // 51
    writebin(ENDIANNESS, sas, 0);
    writebin(pkt3, sas, 0); // 2
    writebin(PLATFORM, sas, 0);

    // packet 3 --------------------------------------- //
    pkt1 = 1, pkt2 = 0, pkt3 = 0, pkt4 = 0;
    writebin(pkt1, sas, 0); // 1
    writebin(pkt2, sas, 0);
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0);

    // packet 5 --------------------------------------- //
    pkt1 = 0, pkt2 = 0, pkt3 = 0, pkt4 = 20;
    writebin(pkt1, sas, 0);
    writebin(pkt2, sas, 0);
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0); // 20

    // packet 5 --------------------------------------- //
    pkt1 = 0, pkt2 = 0, pkt3 = 3, pkt4 = 1;
    writebin(pkt1, sas, 0);
    writebin(pkt2, sas, 0);
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0);

    // packet 6 --------------------------------------- //
    pkt1 = 24, pkt2 = 31, pkt3 = 16, pkt4 = 17;
    writebin(pkt1, sas, 0); // 24
    writebin(pkt2, sas, 0); // 31
    writebin(pkt3, sas, 0); // 16
    writebin(pkt4, sas, 0); // 17

    // packet 7 --------------------------------------- //
    int8_t U64_BYTE_CHECKER_VALUE2 = 51;
    if (bit32 == 1) U64_BYTE_CHECKER_VALUE2 = 50;

    pkt1 = 0, pkt2 = 34, pkt3 = 0, pkt4 = 0;

    writebin(ALIGN_1_CHECKER_VALUE, sas, 0);
    writebin(pkt2, sas, 0); // 34
    writebin(pkt3, sas, 0); // 0
    writebin(U64_BYTE_CHECKER_VALUE2, sas, 0);

    // packet 8 --------------------------------------- //

    pkt1 = 51, pkt2 = 0, pkt3 = 2, pkt4 = 0;
    if (bit32 == 1) pkt1 = 50;
    writebin(pkt1, sas, 0); // 51
    writebin(ENDIANNESS, sas, 0);
    writebin(pkt3, sas, 0); // 2
    writebin(PLATFORM, sas, 0);

    // packet 9 --------------------------------------- //
    pkt1 = 1, pkt2 = 51, pkt3 = 1, pkt4 = 35;
    if (bit32 == 1) pkt2 = 50;
    writebin(pkt1, sas, 0); // 1 | 4
    writebin(pkt2, sas, 0);
    writebin(pkt3, sas, 0);
    writebin(pkt4, sas, 0);

    // packet 10 -------------------------------------- //
    int8_t encoding = 20; // utf8
    pkt1 = 51, pkt2 = 0, pkt3 = 0, pkt4 = 20;
    if (bit32 == 1) pkt1 = 50;
    writebin(pkt1, sas, 0); // 51
    writebin(unk8, sas, 0); // 0
    writebin(encoding, sas, 0);
    writebin(pkt4, sas, 0); // encoding again?

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

    std::string dataset = "TEST";
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

    if (debug) {
      Rcout << "headersize: "  << headersize <<
        " pagesize: " << pagesize << std::endl;
    }

    pagecount_pos1 = sas.tellg();
    if (u64 == 4) {
      writebin(pagecount, sas, swapit);
    } else {
      writebin((int32_t)pagecount, sas, swapit);
    }

    writebin(unkdub, sas, 0);

    std::string sasrel = "9.0401M7";
    sasrel.resize(8);
    sasrel[sasrel.size()] = '\0';

    std::string sasserv = "Linux";
    sasserv.resize(16);
    sasserv[sasserv.size()] = '\0';

    std::string osver = "5.6.15-arch1-1";
    osver.resize(16);
    osver[osver.size()] = '\0';

    std::string osmaker = "";
    osmaker.resize(16);
    osmaker[osmaker.size()] = '\0';

    std::string osname = "x86_64";
    osname.resize(16);
    osname[osname.size()] = '\0';

    writestr(sasrel, 8, sas);
    writestr(sasserv, 16, sas);
    writestr(osver, 16, sas);
    writestr(osmaker, 16, sas);
    writestr(osname, 16, sas);

    // large unk
    uint32_t pktu32 = 1157289805;
    writebin(pktu32, sas, 0);

    pktu32 = 563452161;
    // three identical smaller unks
    writebin(pktu32, sas, 0);
    writebin(pktu32, sas, 0);
    writebin(pktu32, sas, 0);

    // 2 * zero
    writebin(unkdub, sas, 0);
    writebin(unkdub, sas, 0);


    if (debug) Rcout << "pageseqnum32: " << pageseqnum32 << std::endl;
    writebin(pageseqnum32, sas, 0);
    writebin(unk32, sas, 0);

    writebin(thrdts, sas, 0);

    uint64_t num_zeros = headersize - sas.tellg();

    for (uint64_t i = 0; i < num_zeros; ++i) {
      writebin(zero8, sas, swapit);
    }
    // end of header -------------------------------------------------------- //

    auto pos = sas.tellg();

    // sas files consist of pages. pages are written completely. even if data
    // consists of only a few bytes, a full page is written. since SAS requires
    // this, it seems plausible to always write a full page at first and seek to
    // the positions SAS expects afterwards.
    for (auto i = 0; i < (pagesize/sizeof(double)); ++i) {
      writebin(unkdub, sas, 0);
    }

    // rewind
    sas.seekg(pos, sas.beg);

    auto pg = 0;
    uint64_t rowlength = 0;
    rowlength = sum(colwidth);
    auto rows_written = 0;
    uint64_t rows_pending = n;

    /*** write pages **********************************************************/

    // Are these fixed values?
    auto subheader_off =
      100 +
      792 + totalvarnamesize + totalvarformatssize + totalvarlabelssize + // case 1
      24 +
      600 +
      68 +
      28 + k * 8 +
      28 + k * 16 +
      64 * k
    ;

    if (k > 1) subheader_off += 54; // case 8

    if (bit32 == 1) {
      subheader_off =
        100 +
        422 + totalvarnamesize + totalvarformatssize + totalvarlabelssize + // case 1 // kleiner
        12 +          // case 4
        304 +         // case 2
        64 +          // case 5
        24 + k * 8 +  // case 6            // 44 vielleicht kleiner
        24 + k * 12 + // case 7            // 56 vielleicht kleiner
        52 * k        // case 3
      ;

      if (k > 1) subheader_off += 50; // case 8 // something is wrong here
    }

    auto addextra = 0;
    if (((subheader_off) % 8) > 0) {
      // get the required offset from the end of the file so that the content
      // can be written completely
      subheader_off = ceil(8 * ceil((double)subheader_off / 8));
      addextra = 1;
    }

    if (debug)
      Rcout << "SUBHEADER_OFFSET: " << subheader_off << std::endl;


    // end of page 1
    auto pos_at_end_of_page1 = headersize + pagesize - subheader_off;

    // case1: 800 containing + k * 4 (varnames) + k * 4 (formats)
    // case4: 24
    // case2: 600
    // casex: 68
    // casex: 36
    // casex: 44
    // casex: 64 * k (attributes)

    if (debug)
      Rcout << "pos_at_end_of_page1 = " << pos_at_end_of_page1 << std::endl;

    // Rcout << subheader_off << std::endl;
    // Rcout << (headersize + pagecount*pagesize) << std::endl;


    checkUserInterrupt();

    pageseqnum32++;
    unk1 = 0;
    unk2 = 0;
    unk3 = 0; // 62950; // 62676; // 2* pagesize - sas.tellg();

    /* unk3 something like max number of something per page? */
    if (debug) Rcout << "pageseqnum32: " << pageseqnum32 << std::endl;
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


    // not really required, but maybe useful to track positions and length?
    std::vector<PO_Tab> potabs(SUBHEADER_COUNT);

    // all 3 values need to be defined depending on the input
    writebin(PAGE_TYPE, sas, swapit);
    // if n is larger than one page, Block count is off
    block_count_pos1 = sas.tellg();
    writebin(BLOCK_COUNT, sas, swapit);
    writebin(SUBHEADER_COUNT, sas, swapit);
    writebin(unk16, sas, swapit);

    // write entries for all subheaders
    // write them as dummies first, once their position is known fill sh_off.
    //  requires to track the position of each
    pos_SH_C = sas.tellg();
    for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
      if (u64 == 4) {
        if (debug) Rcout << "pos is " <<  sas.tellg() << std::endl;

        writebin(potabs[i].SH_OFF, sas, swapit);
        writebin(potabs[i].SH_LEN, sas, swapit);
        writebin(potabs[i].COMPRESSION, sas, swapit);
        writebin(potabs[i].SH_TYPE, sas, swapit);

        writebin(zero16, sas, swapit);
        writebin(zero16, sas, swapit);
        writebin(zero16, sas, swapit);
      } else {
        writebin((int32_t)potabs[i].SH_OFF, sas, swapit);
        writebin((int32_t)potabs[i].SH_LEN, sas, swapit);
        writebin(potabs[i].COMPRESSION, sas, swapit);
        writebin(potabs[i].SH_TYPE, sas, swapit);

        writebin(zero16, sas, swapit);
      }

    }

    if (debug) Rcout << " ---------- writing data ------------- " << std::endl;

    if (debug)
      Rcout << sas.tellg() << std::endl;

    auto data_pos_page1 = sas.tellg();
    // Rcout << pos_at_end_of_page1 << " < " << data_pos_page1 << " < " << rowlength << std::endl;
    int64_t rows_on_page1 = floor((pos_at_end_of_page1 - data_pos_page1) / (double)rowlength);
    // Rcout << rows_on_page1 << std::endl;

    if (n <= rows_on_page1) {
      rows_on_page1 = n;
    }

    pos = sas.tellg();
    // Rcout << "pos: " << pos << std::endl;
    // Rcout << "block_count_pos1: " << block_count_pos1 << std::endl;
    sas.seekg(block_count_pos1, sas.beg);
    BLOCK_COUNT = rows_on_page1 + SUBHEADER_COUNT;
    writebin(BLOCK_COUNT, sas, swapit);
    sas.seekg(pos, sas.beg);

    bool firstpage = 0;
    if (compress == 0) {
      for (uint64_t i = 0; i < rows_on_page1; ++i) {

        for (auto j = 0; j < k; ++j) {

          // auto ord = ordered[j];
          auto wid = colwidth[j];
          auto typ = vartypes[j];

          if ((wid < 8) && (typ == 1)) {

            double val_d = 0.0;
            val_d = as<NumericVector>(dat[j])[i];

            if ( val_d == NA_REAL  || R_IsNA(val_d) || R_IsNaN(val_d) || std::isinf(val_d) ) {
              uint32_t NA1 = 0;
              uint32_t NA2 = 4294966784;
              writebin(NA1, sas, 0);
              writebin(NA2, sas, 0);
            } else {
              writebin(val_d, sas, swapit);
            }


          }

          if ((wid == 8) && (typ == 1)) {

            double val_d = 0.0;
            val_d = as<NumericVector>(dat[j])[i];

            if ( val_d == NA_REAL  || R_IsNA(val_d) || R_IsNaN(val_d) || std::isinf(val_d) ) {
              uint32_t NA1 = 0;
              uint32_t NA2 = 4294966784;
              writebin(NA1, sas, 0);
              writebin(NA2, sas, 0);
            } else {
              writebin(val_d, sas, swapit);
            }
          }

          if ((wid > 0) && (typ == 2)) {

            int32_t const len = wid;

            CharacterVector cv_s = NA_STRING;
            cv_s = Rcpp::as<CharacterVector>(dat[j])[i];
            std::string val_s = Rcpp::as<std::string>(cv_s);

            if (cv_s[0] == NA_STRING) {
              val_s.clear();
            }

            // rcout << val_s << std::endl;
            writestr(val_s, len, sas);

          }

        }
      }
    }
    if (debug)
      Rcout << sas.tellg() << std::endl;

    rows_written = rows_on_page1;
    rows_pending = n - rows_written;


    // Rcout << "pos_at_end_of_page1: " << pos_at_end_of_page1 << std::endl;
    sas.seekg(pos_at_end_of_page1 , sas.beg);



    /************************************************************************/

    auto shc = SUBHEADER_COUNT;
    uint64_t pre_shlen = 0, post_shlen = 0;
    // while(shc < SUBHEADER_COUNT) {

    int16_t swlen = 0, proclen = 8, comprlen = 0, textoff = 28,
      addtextoff = 7, fmtkey = 0, fmtkey2 = 0,
      fmt32 = 0, fmt322 = 0, ifmt32 = 0, ifmt322 = 0;

    int16_t lenremain16 = 0;


    // SUBHEADER_COUNT - 1 is first pos;
    shc--;

    // final subheader indication where the data begins
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize; /* TODO: Check this */

    if (debug) Rcout << potabs[shc].SH_OFF << std::endl;


    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;
    // on the last subheader
    potabs[shc].COMPRESSION = 1;
    /************************************************************************/

    /**** case 3 ************************************************************/

    for (auto i = 0; i < k; ++i)
    {
      shc--;

      idxofflen fmts, lbls, unks;

      if (debug) Rcout << shc << std::endl;
      pre_shlen = sas.tellg();
      potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

      auto offsetpos = 32;
      if (bit32 == 1) offsetpos = 28; // no clue at all

      // everything is in reverse order - 1 for c index
      auto idx = k-i -1;

      // write offsets in opposite order. calc offsets -1. Presumably SAS
      //  writes the cases to a buffer internally.
      //  Stacks them and writes the entire stack as one
      for (auto z = 0; z < (k - i); ++z) {
        std::string nams = varnames[z];
        offsetpos += nams.size();

        std::string lblz = varlabels[z];
        offsetpos += lblz.size();
        lbls.OFF = offsetpos;

        std::string fmtz = varformats[z];
        offsetpos += fmtz.size();
        fmts.OFF = offsetpos;
      }

      // calc length of len3
      auto len3 = 0;
      for (auto ii = 0; ii < k; ++ii) {
        len3 += 64;
      }

      if (debug)
        Rcout << "-------- case 3 "<< sas.tellg() << std::endl;

      uint32_t case31 = 4294966270;
      uint32_t case32 = 4294967295;
      // int64_t case3 = 4294966270;

      writebin(case31, sas, 0);
      if (bit32 == 0)
        writebin(case32, sas, 0);

      writebin(unk16, sas, swapit);           // 1
      writebin(unk16, sas, swapit);           // 2
      writebin(unk16, sas, swapit);           // 3
      writebin(unk16, sas, swapit);           // 4
      if (vartypes[idx] == 1) {
        fmt32 = width[idx];
      } else {
        fmt32 = 0;
      }
      writebin(fmt32, sas, swapit);          // 5  : formatlen for decimals
      fmt322 = 1; // fix one decimal
      // character or integer
      if (vartypes[idx] == 2) {
        fmt322 = 0;
      }
      writebin(fmt322, sas, swapit);          // 6 : decimal for numerics
      writebin(ifmt32, sas, swapit);          // 7
      writebin(ifmt322, sas, swapit);         // 8
      if (vartypes[idx] == 2) {
        fmtkey = width[idx];
      } else {
        fmtkey = 0;
      }
      writebin(fmtkey, sas, swapit);           // 9  : formatlen for characters
      writebin(fmtkey2, sas, swapit);         // 10 : decimal for characters?
      writebin(unk16, sas, swapit);           // 11
      writebin(unk16, sas, swapit);           // 12
      writebin(unk16, sas, swapit);           // 13
      writebin(unk16, sas, swapit);           // 14 off + len
      writebin(unk16, sas, swapit);           // 15 1 w char

      if (u64 == 4) {
        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);
      }

      // fmt32s.push_back(  fmt32  + (double)fmt322/10);
      // ifmt32s.push_back( ifmt32 + (double)ifmt322/10);
      // fmtkeys.push_back( fmtkey + (double)fmtkey2/10);

      int16_t fmtslen = varformats[idx].size();

      // for characters there is no offset and no format


      if (debug)
        Rcout << varformats[idx] << "___" << offsetpos << "____" << fmtslen << std::endl;

      fmts.IDX = 0, fmts.LEN = fmtslen;
      if (fmtslen == 0) fmts.OFF = 0;

      writebin(fmts.IDX, sas, swapit);
      writebin(fmts.OFF, sas, swapit);
      writebin(fmts.LEN, sas, swapit);

      if (debug)
        Rcout << "fmts: "  << fmts.IDX << ", " << fmts.OFF <<
          ", " << fmts.LEN << std::endl;

      // fmt.push_back(fmts);

      int16_t lblslen = varlabels[idx].size();
      lbls.IDX = 0, lbls.LEN = lblslen;
      if (lblslen == 0) fmts.OFF = 0;

      writebin(lbls.IDX, sas, swapit);
      writebin(lbls.OFF, sas, swapit);
      writebin(lbls.LEN, sas, swapit);

      if (debug)
        Rcout << "lbls: " << lbls.IDX << ", " << lbls.OFF <<
          ", " << lbls.LEN << std::endl;

      // lbl.push_back(lbls);

      writebin(unks.IDX, sas, swapit);
      writebin(unks.OFF, sas, swapit);
      writebin(unks.LEN, sas, swapit);

      // unk.push_back(unks);

      if ((unks.IDX != 0) || (unks.OFF != 0) || (unks.LEN != 0)) {
        warning("case3: unk is not 0 as expected, but %d %d %d\n",
                unks.IDX, unks.OFF, unks.LEN);
        // Rcout << unks.IDX << ", " << unks.OFF <<
        //   ", " << unks.LEN << std::endl;
      }

      // break;
      post_shlen = sas.tellg();
      potabs[shc].SH_LEN = post_shlen - pre_shlen;
    }

    /************************************************************************/


    /**** case 8 ************************************************************/

    // only with more than one variable. maybe related to index
    if (k > 1) {
      shc--;
      if (debug) Rcout << shc << std::endl;
      pre_shlen = sas.tellg();
      potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

      if (debug)
        Rcout << "-------- case 8 "<< sas.tellg() << std::endl;

      uint32_t case81 = 4294967294;
      uint32_t case82 = 4294967295;

      writebin(case81, sas, 0);
      if (bit32 == 0)
        writebin(case82, sas, 0);

      int16_t cls = k;
      int64_t lenremain = 14 +
        cls * 2 + 8 +
        12; // 14 below, (k+2) * 2 and double?

      int32_t ptk32 = 2143813666;
      writebin(ptk32, sas, swapit); // unkown large number
      writebin(unk16, sas, swapit); // 0
      writebin(unk16, sas, swapit); // 0

      if (u64 == 4) {  // lenremain
        writebin(lenremain, sas, swapit);
      } else {
        writebin((int32_t)lenremain, sas, swapit);
      }

      if (debug)
        Rcout << "lenremain "<< lenremain << std::endl;

      writebin((int16_t)k, sas, swapit);  // number of varnames?
      writebin(cls, sas, swapit);    // counter for unk loop below
      int16_t ptk16 = 1;
      writebin(ptk16, sas, swapit);  // 1
      writebin((int16_t)k, sas, swapit);  // number of varnames?
      writebin(unk16, sas, swapit);  // 3233
      writebin(unk16, sas, swapit);  // 3233
      writebin(unk16, sas, swapit);  // 3233

      lenremain -= 14;

      // Rcout << lenremain << " " << cls << std::endl;

      int16_t res = 0;
      // for (auto cl = 0; cl < cls; ++cl) {
      ptk16 = -1;
      writebin(ptk16, sas, swapit);
      writebin((int16_t)k, sas, swapit);
      writebin(unk16, sas, swapit);
      writebin(unk16, sas, swapit);
      // }

      // to be on par with lenremain
      writebin(unkdub, sas, 0);

      // // padding? Not in lenremain
      writebin(unk16, sas, swapit); // 0
      // writebin(unk16, sas, swapit); // 0
      // writebin(unk16, sas, swapit); // 0
      // writebin(unk16, sas, swapit); // 0
      // writebin(unk16, sas, swapit); // 0


      // break;
      post_shlen = sas.tellg();
      potabs[shc].SH_LEN = post_shlen - pre_shlen;
      potabs[shc].SH_TYPE = 1;
    }

    /************************************************************************/


    /**** case 7 ************************************************************/

    /* Column Attributes */

    shc--;
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

    if (debug)
      Rcout << "-------- case 7 "<< sas.tellg()  << std::endl;

    uint32_t case71 = 4294967292;
    uint32_t case72 = 4294967295;

    writebin(case71, sas, 0);
    if (bit32 == 0)
      writebin(case72, sas, 0);

    int8_t divs = 16;
    if (u64 != 4) divs = 12;

    lenremain16 =
      k * divs + 8
    ;
    if (debug) Rprintf("lenremain16 %d \n", lenremain16);

    writebin(lenremain16, sas, swapit);
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0

    /* Column Attributes Pointers
     * should be filled from R or somewhere above
     */
    std::vector<CN_Att> capois(k);

    auto prevoffset = 0;
    for (auto i = 0; i < k; ++i) {

      capois[i].CN_OFF = prevoffset;

      if (u64 == 4) {
        writebin(capois[i].CN_OFF, sas, swapit);
      } else {
        writebin((int32_t)capois[i].CN_OFF,
                 sas, swapit);
      }

      capois[i].CN_WID = colwidth[i];
      capois[i].NM_FLAG = 1024; // ?
      capois[i].CN_TYP = vartypes[i];

      writebin(capois[i].CN_WID, sas, swapit);
      writebin(capois[i].NM_FLAG, sas, swapit);
      writebin(capois[i].CN_TYP, sas, swapit);
      writebin(capois[i].UNK8, sas, swapit);

      prevoffset += capois[i].CN_WID;
    }

    writebin(unk32, sas, swapit); // 0
    writebin(unk32, sas, swapit); // 0
    writebin(unk32, sas, swapit); // 0


    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;
    potabs[shc].SH_TYPE = 1;

    /************************************************************************/


    /**** case 6 ************************************************************/

    /* Column Name */

    shc--;
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

    if (debug)
      Rcout << "-------- case 6 "<< sas.tellg() << std::endl;

    uint32_t case61 = 4294967295;
    uint32_t case62 = 4294967295;

    writebin(case61, sas, 0);
    if (bit32 == 0)
      writebin(case62, sas, 0);

    lenremain16 = k * 8 + 8; // empty double at the end?
    if (debug) Rprintf("lenremain16 %d \n", lenremain16);

    writebin(lenremain16, sas, swapit);
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0

    lenremain16 -= 8;
    int8_t div = 8;

    auto cmax = lenremain16 / div;

    /* Column Name Pointers */
    std::vector<CN_Poi> cnpoi(k);



    inioff = 36;

    auto prevlen = 0 + inioff;

    for (auto cn = 0; cn < k; ++cn) {
      cnpoi[cn].CN_IDX = pg;
      cnpoi[cn].CN_OFF = prevlen;
      if (debug) Rcout << "as<std::string>(nvarnames[cn]).size();" << as<std::string>(nvarnames[cn]).size() << std::endl;
      if (debug) Rcout << "varnames[cn].size();" << varnames[cn].size() << std::endl;
      // varnames size < 4 are padded. This len field
      // tells SAS the actual size of the varname. max value is 32
      int8_t max_varname_len = 32L;
      cnpoi[cn].CN_LEN = std::min((int8_t)as<std::string>(nvarnames[cn]).size(), max_varname_len);
      prevlen += varnames[cn].size();
      prevlen += varlabels[cn].size();
      prevlen += varformats[cn].size();

      if(debug) Rcout << "prevlen: " << prevlen << std::endl;
    }


    for (auto cn = 0; cn < k; ++cn) {
      writebin(cnpoi[cn].CN_IDX, sas, swapit);
      writebin(cnpoi[cn].CN_OFF, sas, swapit);
      writebin(cnpoi[cn].CN_LEN, sas, swapit);
      writebin(cnpoi[cn].zeros,  sas, swapit);
    }
    writebin(unkdub, sas, 0);

    // padding? Not in lenremain
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0


    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;
    potabs[shc].SH_TYPE = 1;

    /************************************************************************/


    /**** case 5 ************************************************************/

    /* Column Text */

    shc--;
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

    if (debug)
      Rcout << "-------- case 5 "<< sas.tellg() << std::endl;

    uint32_t case51 = 4294967293;
    uint32_t case52 = 4294967295;

    writebin(case51, sas, 0);
    if (bit32 == 0)
      writebin(case52, sas, 0);


    int16_t len =
      6 +                 // 3 * unk16
      16 +                // empty string
      8 +                 // proc
      totalvarnamesize +  // varnamesize %% 4
      totalvarlabelssize +
      totalvarformatssize +
      2                    // int32 at end
      ;
    if (k > 1) len += 4;

    auto c5first = 0;
    auto c5typ = 0;

    writebin(len, sas, swapit);
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0
    writebin(unk16, sas, swapit); // 0

    if ((PAGE_TYPE != 1024) && (c5first == 0)) {
      writebin(unk16, sas, swapit); // 0 |     0 | 27977
      int16_t ptk16 = 5120; // related to visual representation 0 = fills full space
      writebin(ptk16, sas, swapit); // 0 | 15872 | 30064
    }

    // len starting here

    if ((c5typ == 0) && (pg == 0)) {

      // compression
      if (compress == 0) {
        std::string none = " ";
        writestr(none, 16, sas);
      }

      std::string proc = "DATASTEP";
      writestr(proc, 8, sas);

      for (auto i = 0; i < k; ++i) {
        writestr(varnames[i], varnames[i].size(), sas);
        writestr(varlabels[i], varlabels[i].size(), sas);
        writestr(varformats[i], varformats[i].size(), sas);
      }

      // padding
      // if(addextra == 1)
      writebin(unk32, sas, 0); // presumably wrong?

      // handling of labels is similar, but labels may exceed size of 8
      // still must be dividable by 8 (maybe 4) if 14, add 2

      writebin(unk32, sas, 0);
      writebin(unk32, sas, 0);

      // // additional software string
      // if (swlen > 0) {
      //   sw.resize(swlen, '\0');
      //   sw = readstring(sw, sas);
      // }

      ++c5typ;
    }


    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;
    potabs[shc].SH_TYPE = 1;

    /************************************************************************/


    /**** case 2 ************************************************************/

    shc--;
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

    if (debug)
      Rcout << "-------- case 2 "<< sas.tellg() << std::endl;

    uint32_t case21 = 4294966272;
    uint32_t case22 = 4294967295;

    writebin(case21, sas, 0);
    if (bit32 == 0)
      writebin(case22, sas, 0);

    inioff = 36;
    if (bit32 == 1) inioff = 32;

    int64_t off = inioff + totalvarnamesize + totalvarformatssize; // offset with format

    int64_t ptk64 = 3;
    if (k > 1) ptk64++;

    if (u64 == 4) {
      writebin(off, sas, swapit);
      writebin(ptk64, sas, swapit);
    } else {
      writebin((int32_t)off, sas, swapit);
      writebin(unk32, sas, swapit);
    }

    int16_t num_nonzero = 7;
    writebin(num_nonzero, sas, swapit);

    int8_t unklen = 94; // should be 94
    if (u64 != 4) unklen = 50;


    unklen = unklen/2 - 4;
    for (int jj = 0; jj < unklen; ++jj) {
      writebin(unk16, sas, swapit);
      // 4th from the end is 1804 meaning is unknown
    }

    // int16_t ptk16 = 1804;
    int16_t ptk16 = 0;
    writebin(ptk16, sas, swapit);
    writebin(unk16, sas, swapit);
    writebin(unk16, sas, swapit);
    writebin(unk16, sas, swapit);

    std::vector<SCV> scv(12);

    for (int8_t i = 0; i < 12; ++i) {


      // SIG -4; FIRST 1; F_POS 6; LAST 1; L_POS 6
      if (i == 0) {
        scv[i].SIG = -4;
        scv[i].FIRST = 1;
        scv[i].F_POS = 6;
        scv[i].LAST = 1;
        scv[i].L_POS = 6;
      }

      // SIG -3; FIRST 1; F_POS 4; LAST 1; L_POS 4
      if (i == 1) {
        scv[i].SIG = -3;
        scv[i].FIRST = 1;
        scv[i].F_POS = 4;
        scv[i].LAST = 1;
        scv[i].L_POS = 4;
      }

      // SIG -1; FIRST 1; F_POS 5; LAST 1; L_POS 5
      if (i == 2) {
        scv[i].SIG = -1;
        scv[i].FIRST = 1;
        scv[i].F_POS = 5;
        scv[i].LAST = 1;
        scv[i].L_POS = 5;
      }

      // SIG -2; FIRST 1; F_POS 7; LAST 1; L_POS 7
      if (i == 3) {
        scv[i].SIG = -2;
        scv[i].FIRST = 1;
        scv[i].F_POS = 7;
        scv[i].LAST = 1;
        scv[i].L_POS = 7;
      }

      // SIG -5; FIRST 0; F_POS 0; LAST 0; L_POS 0
      if (i == 4) {
        scv[i].SIG = -5;
        scv[i].FIRST = 0;
        scv[i].F_POS = 0;
        scv[i].LAST = 0;
        scv[i].L_POS = 0;
      }
      // SIG -6; FIRST 0; F_POS 0; LAST 0; L_POS 0
      if (i == 5) {
        scv[i].SIG = -6;
        scv[i].FIRST = 0;
        scv[i].F_POS = 0;
        scv[i].LAST = 0;
        scv[i].L_POS = 0;
      }
      // SIG -7; FIRST 0; F_POS 0; LAST 0; L_POS 0
      if (i == 6) {
        scv[i].SIG = -7;
        scv[i].FIRST = 0;
        scv[i].F_POS = 0;
        scv[i].LAST = 0;
        scv[i].L_POS = 0;
      }

      if (u64 == 4) {
        writebin(scv[i].SIG, sas, swapit);
        writebin(scv[i].FIRST, sas, swapit);
        writebin(scv[i].F_POS, sas, swapit);

        if ((i == 0) && (scv[i].SIG != -4))
          warning("first SIG is not -4");

        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);

        writebin(scv[i].LAST, sas, swapit);
        writebin(scv[i].L_POS, sas, swapit);

        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);
        writebin(unk16, sas, swapit);

      } else {
        writebin((int32_t)scv[i].SIG, sas, swapit);
        writebin((int32_t)scv[i].FIRST, sas, swapit);
        writebin(scv[i].F_POS, sas, swapit);

        writebin(unk16, sas, swapit);

        writebin((int32_t)scv[i].LAST, sas, swapit);
        writebin(scv[i].L_POS, sas, swapit);

        writebin(unk16, sas, swapit);
      }

      if (debug)
        Rprintf("SIG %d; FIRST %d; F_POS %d; LAST %d; L_POS %d\n",
                scv[i].SIG, scv[i].FIRST, scv[i].F_POS,
                scv[i].LAST, scv[i].L_POS);

    }


    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;

    /************************************************************************/


    /**** case 4 ************************************************************/

    /* Column Size */

    if (debug)
      Rcout << "-------- case 4 "<< sas.tellg() << std::endl;

    shc--;
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

    uint32_t case41 = 4143380214;
    uint32_t case42 = 0;

    writebin(case41, sas, 0);
    if (bit32 == 0)
      writebin(case42, sas, 0);

    uint64_t uunk64 = 0;
    int64_t colnum = k;

    if (u64 == 4) {
      writebin(colnum, sas, swapit);
      writebin(uunk64, sas, swapit);
    } else {
      writebin((int32_t)colnum, sas, swapit);
      writebin((int32_t)uunk64, sas, swapit);
    }
    if (debug)
      Rprintf("colnum %d; uunk64 %d\n",
              colnum, uunk64);

    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;

    /************************************************************************/

    /**** case 1 ************************************************************/

    shc--;
    if (debug) Rcout << shc << std::endl;
    pre_shlen = sas.tellg();
    potabs[shc].SH_OFF = pre_shlen - headersize - pg*pagesize;

    uint32_t case11 = 4160223223;
    uint32_t case12 = 0; //4294967295;

    writebin(case11, sas, 0);
    if (bit32 == 0)
      writebin(case12, sas, 0);

    /* Row Size */

    int16_t pgwpossh = 0, pgwpossh2 = 0, numzeros = 37,
      sh_num = 0, cn_maxlen = 0, l_maxlen = 0, todata = 0,
      rowsonpg = 0;
    int32_t sasvers = 8;

    if (debug)
      Rcout << "-------- case 1 "<< sas.tellg() << std::endl;

    int64_t rowcount = 0, delobs = 0;

    rowcount = n;

    int64_t pgsize = pagesize, pgc = 0,
      rcmix = rows_on_page1, // looks like number of rows on rcmix pages, but is just a larger number?
      pgwsh = 0, pgwsh2 = 0; //

    int64_t colf_p1 = k, colf_p2 = 0;
    int16_t dataoffset = 256;

    int64_t
    ptk64_01 = 1,
      ptk64_02 = 2,
      ptk64_03 = 1,
      ptk64_04 = 7,
      ptk64_05 = 1,
      ptk64_06 = 9,
      ptk64_07 = 1,
      ptk64_08 = BLOCK_COUNT,
      ptk64_09 = 1,
      ptk64_10 = 7;

    if (k > 1) {
      ptk64_04 = 7+k; // namestr
      ptk64_06 += k;
      ptk64_10++;
    }

    if (u64 == 4) {

      /* */
      // pkt64_2 is +2 with k * 2
      int64_t pkt64_1 = 240, pkt64_2 = 21, pkt64_3 = 0, pkt64_4 = 2240529;
      writebin(pkt64_1, sas, swapit);
      writebin(pkt64_2, sas, swapit);
      writebin(pkt64_3, sas, swapit);
      writebin(pkt64_4, sas, swapit);

      writebin(rowlength, sas, swapit); // rowlength
      rowcount_pos = sas.tellg();
      writebin(rowcount, sas, swapit);  // rowcount
      writebin(delobs, sas, swapit);
      writebin(unk64, sas, swapit);

      writebin(colf_p1, sas, swapit); // k related
      writebin(colf_p2, sas, swapit);
      int64_t ptk64 = 0; // changing this makes SAS complain about missing tables
      writebin(ptk64, sas, swapit); // p3 and p4?
      ptk64 = 34; // was 8 for k = 1
      writebin(ptk64, sas, swapit);
      writebin(pgsize, sas, swapit);
      writebin(unk64, sas, swapit);
      rcmix_pos = sas.tellg();
      writebin(rcmix, sas, swapit);

      uint64_t uunk64 = 0;
      /* next two indicate the end of the initial header ? */

      ptk64 = -1;
      writebin(ptk64, sas, swapit);
      writebin(ptk64, sas, swapit);

      for (int z = 0; z < numzeros; ++z) {
        writebin(unk64, sas, swapit);
      }

      writebin(pgidx, sas, swapit);

      // padding? 68 bytes: zeros
      for (int z = 0; z < 8; ++z) {
        writebin(unk64, sas, swapit);
      }

      writebin(unk32, sas, swapit); // padding

      ptk64_07 = pagecount;
      writebin(ptk64_01, sas, swapit); // val 1?
      writebin(ptk64_02, sas, swapit); // val 2?
      writebin(ptk64_03, sas, swapit); // val 1?
      writebin(ptk64_04, sas, swapit); // val 7?
      writebin(ptk64_05, sas, swapit); // val 7?
      writebin(ptk64_06, sas, swapit); // val 7?
      pagecount_pos2 = sas.tellg();
      writebin(ptk64_07, sas, swapit); // val 7?
      writebin(ptk64_08, sas, swapit); // val 7?
      writebin(ptk64_09, sas, swapit); // val 7?
      writebin(ptk64_10, sas, swapit); // val 7?

      // writebin(pgwsh, sas, swapit);
      // writebin(pgwpossh, sas, swapit);
      // writebin(unk16, sas, swapit); // padding
      // writebin(unk16, sas, swapit); // padding

      // writebin(unk16, sas, swapit); // padding
      // writebin(pgwsh2, sas, swapit);
      // writebin(pgwpossh2, sas, swapit);
      // writebin(unk16, sas, swapit); // padding

      // writebin(unk16, sas, swapit); // padding
      // writebin(unk16, sas, swapit); // padding
      //
      // writebin(pgc, sas, swapit);
      //
      // writebin(unk16, sas, swapit); // val ?
      // writebin(unk16, sas, swapit); // padding
      // writebin(unk16, sas, swapit); // padding
      // writebin(unk16, sas, swapit); // padding
      //
      // writebin(unk64, sas, swapit); // val 1?
      //
      // writebin(addtextoff, sas, swapit); // val 7 | 8?
      // writebin(unk16, sas, swapit); // padding
      // writebin(unk16, sas, swapit); // padding
      // writebin(unk16, sas, swapit); // padding

      for (int z = 0; z < 10; ++z) {
        writebin(unk64, sas, swapit); // 0
      }

      int16_t ptk16_1 = 0, ptk16_2 = 8, ptk16_3 = 4, ptk16_4 = 0;
      writebin(ptk16_1, sas, swapit); // val
      writebin(ptk16_2, sas, swapit); // val 0|8 ?
      writebin(ptk16_3, sas, swapit); // val 4
      writebin(ptk16_4, sas, swapit); // val 0
      writebin(todata, sas, swapit); // val 12,32|0?
      // if (todata == 12) hasproc = false;

      ptk16_1 = 0, ptk16_2 = 0, ptk16_3 = 20, ptk16_4 = 8;
      writebin(ptk16_1, sas, swapit);
      writebin(ptk16_2, sas, swapit); // val 0
      writebin(ptk16_3, sas, swapit); // val 20 | 28
      writebin(ptk16_4, sas, swapit); // 8

      writebin(unk16, sas, swapit); // 0
      writebin(unk16, sas, swapit); // 12
      writebin(unk16, sas, swapit); // 8
      writebin(unk16, sas, swapit); // 0

      ptk16_1 = 12, ptk16_2 = 8, ptk16_3 = 0, ptk16_4 = 0;
      writebin(ptk16_1, sas, swapit); // 12 | 20
      writebin(ptk16_2, sas, swapit); // 8
      writebin(ptk16_3, sas, swapit); // 0
      writebin(textoff, sas, swapit); // 28 | 36
      writebin(proclen, sas, swapit);

      for (int z = 0; z < 8; ++z) {
        writebin(unk32, sas, swapit); // 0
      }

      writebin(unk16, sas, swapit); // 4
      writebin(unk16, sas, swapit); // 1


      sh_num = 4, cn_maxlen = 1, l_maxlen = 1;
      writebin(sh_num, sas, swapit);
      writebin(cn_maxlen, sas, swapit);
      writebin(l_maxlen, sas, swapit);

      /* maybe SAS version information at o131018 ? */
      writebin(sasvers, sas, swapit); // 1
      writebin(unk32, sas, swapit); // 2
      writebin(unk32, sas, swapit); // 3

      writebin(rowsonpg, sas, swapit);


      writebin(unk16, sas, swapit); // 1

      int32_t pkt32 = 0;
      writebin(pkt32, sas, swapit); // 2
      writebin(unk16, sas, swapit); // 4
      writebin(unk16, sas, swapit); // 5
      writebin(unk16, sas, swapit); // 8
      // pkt32 = 1;
      // writebin(pkt32, sas, swapit); // 10 // width?
      // writebin(unk32, sas, swapit); // 10
      // ptk64 = -1;
      writebin(n, sas, swapit); // 6
      // writebin(unk32, sas, swapit); // 10
      writebin(unk16, sas, swapit); // 9
      writebin(unk16, sas, swapit); // 12
      writebin(unk16, sas, swapit); // 13
      writebin(unk16, sas, swapit); // 14
      writebin(unk16, sas, swapit); // 15
      writebin(dataoffset, sas, swapit); // 16
      writebin(unk16, sas, swapit); // 17
      writebin(unk16, sas, swapit);// 18
      writebin(unk16, sas, swapit); // 19
      writebin(unk16, sas, swapit); // 20

      writebin(unk32, sas, swapit); // extra
      /* */

    } else {

      writebin(unk32, sas, swapit);
      writebin(unk32, sas, swapit);
      writebin(unk32, sas, swapit);
      writebin(unk32, sas, swapit);

      writebin((int32_t)rowlength, sas, swapit);
      writebin((int32_t)rowcount, sas, swapit);
      writebin((int32_t)delobs, sas, swapit); // deleted obs?
      writebin(unk32, sas, swapit);
      writebin((int32_t)colf_p1, sas, swapit);
      writebin((int32_t)colf_p2, sas, swapit);
      writebin(unk32, sas, swapit);
      writebin(unk32, sas, swapit);
      writebin((int32_t)pgsize, sas, swapit);
      writebin(unk32, sas, swapit);
      writebin((int32_t)rcmix, sas, swapit);
      writebin(uunk32, sas, swapit);
      writebin(uunk32, sas, swapit);

      for (int z = 0; z < numzeros; ++z) {
        writebin((int32_t)unk64, sas, swapit);
      }

      writebin(pgidx, sas, swapit);

      for (int z = 0; z < 10; ++z) {
        writebin((int32_t)unk64, sas, swapit);;
      }

      ptk64_07 = pagecount;
      writebin((int32_t)ptk64_01, sas, swapit); // val 1?
      writebin((int32_t)ptk64_02, sas, swapit); // val 2?
      writebin((int32_t)ptk64_03, sas, swapit); // val 1?
      writebin((int32_t)ptk64_04, sas, swapit); // val 7?
      writebin((int32_t)ptk64_05, sas, swapit); // val 7?
      writebin((int32_t)ptk64_06, sas, swapit); // val 7?
      pagecount_pos2 = sas.tellg();
      writebin((int32_t)ptk64_07, sas, swapit); // val 7?
      writebin((int32_t)ptk64_08, sas, swapit); // val 7?
      writebin((int32_t)ptk64_09, sas, swapit); // val 7?
      writebin((int32_t)ptk64_10, sas, swapit); // val 7?

      // // padding?
      // writebin(unk32, sas, swapit);
      // writebin(unk32, sas, swapit);
      //
      // writebin(unk32, sas, swapit); // val 1?
      // writebin(unk16, sas, swapit); // val 2?
      //
      // writebin(unk16, sas, swapit); // padding
      // writebin((int32_t)pgwsh, sas, swapit);
      // writebin(pgwpossh, sas, swapit);
      //
      // writebin(unk16, sas, swapit); // padding
      //
      // writebin((int32_t)pgwsh2, sas, swapit);
      // writebin(pgwpossh2, sas, swapit);
      // writebin(unk16, sas, swapit); // padding
      // writebin((int32_t)pgc, sas, swapit);
      //
      // writebin(unk16, sas, swapit); // val ?
      // writebin(unk16, sas, swapit); // padding
      //
      // writebin((int32_t)unk64, sas, swapit); // val 1?
      //
      // writebin(addtextoff, sas, swapit); // val 7 | 8?
      // writebin(unk16, sas, swapit); // padding

      for (int z = 0; z < 10; ++z) {
        writebin((int32_t)unk64, sas, swapit); // 0
      }

      int16_t ptk16_1 = 0, ptk16_2 = 8, ptk16_3 = 4, ptk16_4 = 0;
      writebin(ptk16_1, sas, swapit); // val
      writebin(ptk16_2, sas, swapit); // val 0|8 ?
      writebin(ptk16_3, sas, swapit); // val 4
      writebin(ptk16_4, sas, swapit); // val 0
      writebin(todata, sas, swapit); // val 12,32|0? //

      ptk16_1 = 0, ptk16_2 = 0, ptk16_3 = 20, ptk16_4 = 8;
      writebin(swlen, sas, swapit);
      writebin(ptk16_2, sas, swapit); // val 0?
      writebin(ptk16_3, sas, swapit); // val 20?
      writebin(ptk16_4, sas, swapit); //

      writebin(unk16, sas, swapit); // 0
      writebin(unk16, sas, swapit); // 12
      writebin(unk16, sas, swapit); // 8 compr. code length?
      writebin(unk16, sas, swapit); // 0

      ptk16_1 = 12, ptk16_2 = 8, ptk16_3 = 0, ptk16_4 = 0;
      writebin(ptk16_1, sas, swapit); // 12
      writebin(ptk16_2, sas, swapit); // 8
      writebin(ptk16_3, sas, swapit); // 0
      writebin(textoff, sas, swapit); // 28
      writebin(proclen, sas, swapit);

      writebin(unk32, sas, swapit); // 1
      writebin(unk32, sas, swapit); // 2
      writebin(unk32, sas, swapit); // 3
      writebin(unk32, sas, swapit); // 4
      writebin(unk32, sas, swapit); // 5
      writebin(unk32, sas, swapit); // 6
      writebin(unk32, sas, swapit); // 7
      writebin(unk32, sas, swapit); // 8
      writebin(unk16, sas, swapit); // 4
      writebin(unk16, sas, swapit); // 1

      sh_num = 4, cn_maxlen = 1, l_maxlen = 1;
      writebin(sh_num, sas, swapit);
      writebin(cn_maxlen, sas, swapit);
      writebin(l_maxlen, sas, swapit);

      writebin(sasvers, sas, swapit); // 1
      writebin(unk32, sas, swapit); // 2
      writebin(unk32, sas, swapit); // 3

      writebin(rowsonpg, sas, swapit);

      writebin(unk16, sas, swapit); // 1

      int32_t pkt32 = 0;
      writebin(pkt32, sas, swapit); // 2
      writebin(unk16, sas, swapit); // 4
      writebin(unk16, sas, swapit); // 5
      writebin(unk16, sas, swapit); // 6
      writebin(unk32, sas, swapit); // 7
      writebin(unk16, sas, swapit); // 9
      writebin(unk32, sas, swapit); // 10
      writebin(unk16, sas, swapit); // 12
      writebin(unk16, sas, swapit); // 13
      writebin(unk16, sas, swapit); // 14
      writebin(unk16, sas, swapit); // 15
      writebin(dataoffset, sas, swapit); // 16
      writebin(unk16, sas, swapit); // 17
      writebin(unk16, sas, swapit);// 18
      writebin(unk16, sas, swapit); // 19
      writebin(unk16, sas, swapit); // 20
    }


    // if (debug)
    //   Rprintf("swlen = %d, todata %d, textoff %d\n",
    //           swlen, todata, textoff);
    //
    //
    // if (!((dataoffset == 256) || (dataoffset == 1280)) && debug)
    //   warning("debug: dataoffset is unexpectedly %d\n",
    //           dataoffset);

    post_shlen = sas.tellg();
    potabs[shc].SH_LEN = post_shlen - pre_shlen;

    /************************************************************************/


    auto exp_endpos = headersize + pagesize; // pagecount *
    auto got_endpos = sas.tellg();

    if (debug) Rcout << "END OF FILE REACHED " << got_endpos << std::endl;
    if (got_endpos > exp_endpos)
      warning("Unexpected position at end of file. %d and not %d",
              got_endpos, exp_endpos);


    /*write correct offsets *************************************************/

    // Rcout << "pos_SH_C: " << pos_SH_C << std::endl;
    sas.seekg(pos_SH_C, sas.beg);
    if (debug) Rcout << pos_SH_C << std::endl;
    if (debug) Rcout << "Writing subheader pos at " << sas.tellg() << std::endl;
    // Rcout << "pos is " << sas.tellg() << " and " << pos_SH_C << std::endl;
    for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
      if (u64 == 4) {
        writebin(potabs[i].SH_OFF, sas, 0);
        writebin(potabs[i].SH_LEN, sas, 0);
        writebin(potabs[i].COMPRESSION, sas, 0);
        writebin(potabs[i].SH_TYPE, sas, 0);

        writebin(zero16, sas, 0);
        writebin(zero16, sas, 0);
        writebin(zero16, sas, 0);

      } else {

        writebin((int32_t)potabs[i].SH_OFF, sas, swapit);
        writebin((int32_t)potabs[i].SH_LEN, sas, swapit);
        writebin(potabs[i].COMPRESSION, sas, swapit);
        writebin(potabs[i].SH_TYPE, sas, swapit);

        writebin(zero16, sas, swapit);

      }
    }


    // write additional pages
    while (rows_written < n) {

      // pos at begin of page
      pos = headersize + (double)pagecount * pagesize;

      ++pagecount;

      if (debug)
        Rcout << "@@@@ write additoinal page @@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

      // Rcout << "pos: " << pos << std::endl;
      sas.seekg(pos, sas.beg);

      // write empty page
      for (auto i = 0; i < (pagesize/sizeof(double)); ++i) {
        writebin(unkdub, sas, 0);
      }
      // Rcout << "pos: " << pos << std::endl;
      sas.seekg(pos, sas.beg);

      int64_t pos_beg = (int64_t)pos + 40;
      int64_t pos_end = headersize + ((double)pagecount) * pagesize;

      // Rcout << pos_at_end_of_page1 << " < " << data_pos_page1 << " < " << rowlength << std::endl;

      // potential number of rows on this page
      BLOCK_COUNT = floor((pos_end - pos_beg) / (double)rowlength);
      if (BLOCK_COUNT > rows_pending) {
        BLOCK_COUNT = rows_pending;
      } else {
        rows_pending -= BLOCK_COUNT;
      }

      auto rows_on_page = BLOCK_COUNT + rows_written;
      if (debug) Rcout << rows_on_page << std::endl;

      // if (n <= (rows_on_page)) {
      //   rows_on_page = n;
      // }

      ++pageseqnum32;
      unk32 = 0, unk1 = 0, unk2 = 0, unk3 = 0;

      // Page Offset Table
      if (debug) Rcout << "pageseqnum32: " << pageseqnum32 << std::endl;
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

      // not really required, but maybe useful to track positions and length?
      std::vector<PO_Tab> potabs(SUBHEADER_COUNT);

      PAGE_TYPE = 256;
      // BLOCK_COUNT = 1;
      SUBHEADER_COUNT = 0;

      // all 3 values need to be defined depending on the input
      writebin(PAGE_TYPE, sas, swapit);
      writebin(BLOCK_COUNT, sas, swapit);
      writebin(SUBHEADER_COUNT, sas, swapit);
      writebin(unk16, sas, swapit);


      if (debug) Rcout << " ---------- writing data ------------- " << std::endl;

      // rows_written;

      if (debug) {
        Rcout << sas.tellg() << std::endl;
        Rcout << "from i = " << rows_written << std::endl;
        Rcout << "to i < " << rows_on_page << std::endl;
      }

      bool firstpage = 0;
      if (compress == 0) {

        for (uint64_t i = rows_written; i < rows_on_page; ++i) {

          for (auto j = 0; j < k; ++j) {

            // auto ord = ordered[j];
            auto wid = colwidth[j];
            auto typ = vartypes[j];

            if ((wid < 8) && (typ == 1)) {

              double val_d = 0.0;
              val_d = as<NumericVector>(dat[j])[i];

              if ( val_d == NA_REAL  || R_IsNA(val_d) || R_IsNaN(val_d) || std::isinf(val_d) ) {
                uint32_t NA1 = 0;
                uint32_t NA2 = 4294966784;
                writebin(NA1, sas, 0);
                writebin(NA2, sas, 0);
              } else {
                writebin(val_d, sas, swapit);
              }


            }

            if ((wid == 8) && (typ == 1)) {

              double val_d = 0.0;
              val_d = as<NumericVector>(dat[j])[i];

              if ( val_d == NA_REAL  || R_IsNA(val_d) || R_IsNaN(val_d) || std::isinf(val_d) ) {
                uint32_t NA1 = 0;
                uint32_t NA2 = 4294966784;
                writebin(NA1, sas, 0);
                writebin(NA2, sas, 0);
              } else {
                writebin(val_d, sas, swapit);
              }
            }

            if ((wid > 0) && (typ == 2)) {

              int32_t const len = wid;

              CharacterVector cv_s = NA_STRING;
              cv_s = Rcpp::as<CharacterVector>(dat[j])[i];
              std::string val_s = Rcpp::as<std::string>(cv_s);

              if (cv_s[0] == NA_STRING) {
                val_s.clear();
              }

              // rcout << val_s << std::endl;
              writestr(val_s, len, sas);

            }

          }

        }
      }
      if (debug)
        Rcout << sas.tellg() << std::endl;

      rows_written = rows_on_page;

    }

    // update pagecount in file header
    // Rcout << "pagecount_pos1: " << pagecount_pos1 << std::endl;
    sas.seekg(pagecount_pos1, sas.beg);
    if (u64 == 4) {
      writebin(pagecount, sas, swapit);
    } else {
      writebin((int32_t)pagecount, sas, swapit);
    }

    sas.seekg(pagecount_pos2, sas.beg);
    // Rcout << "pagecount_pos2: "<< pagecount_pos2 << std::endl;
    if (u64 == 4) {
      writebin(pagecount, sas, swapit);
    } else {
      writebin((int32_t)pagecount, sas, swapit);
    }

  }
}
