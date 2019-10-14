/*
 * Copyright (C) 2019 Jan Marvin Garbuszus
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
#include <stdint.h>
#include <string>
#include <fstream>
#include <streambuf>
#include <regex>

#include "sas.h"

using namespace Rcpp;


//' Reads SAS data files
//'
//' @param filePath The full systempath to the sas7bdat file you want to import.
//' @param debug print debug information
//' @import Rcpp
//' @export
// [[Rcpp::export]]
Rcpp::List readsas(const char * filePath, const bool debug)
{
  std::ifstream sas(filePath, std::ios::in | std::ios::binary);
  if (sas) {


    auto k = 0;
    int compr = 0;
    auto ctr = 0, varoffset = 0, offset = 0;

    bool hasattributes = 0, hasproc = 1, swapit = 0, c5first = 0;
    int8_t  unk8  = 0;
    int8_t ALIGN_1_CHECKER_VALUE = 0;
    int8_t ALIGN_2_VALUE = 0;
    int8_t ENDIANNESS = 0;
    int16_t unk16 = 0;
    int32_t unk32 = 0;
    int64_t unk64 = 0;
    double unkdub = 0;

    int8_t u64 = 0;
    int16_t dataoffset = 0;
    int16_t swlen = 0, proclen = 0, comprlen = 0, textoff = 0, todata = 0,
      addtextoff = 0, fmtkey = 0, fmtkey2 = 0,
      fmt32 = 0, fmt322 = 0, ifmt32 = 0, ifmt322 = 0;
    int16_t PAGE_TYPE = 0, BLOCK_COUNT = 0, SUBHEADER_COUNT = 0;

    uint32_t pageseqnum32 = 0;
    uint64_t pageseqnum64 = 0;

    double created = 0, created2 = 0; // 8
    double modified = 0, modified2 = 0; // 16

    std::string compression = "";
    std::string proc = "";
    std::string sw = "";
    std::string sasfile  (8, '\0');
    std::string filetype (8, '\0');
    std::string sasrel   (8, '\0');
    std::string sasserv  (16, '\0');
    std::string osver    (16, '\0');
    std::string osmaker  (16, '\0');
    std::string osname   (16, '\0');
    std::string dataset  (64, '\0');

    std::vector<idxofflen> fmt;
    std::vector<idxofflen> lbl;
    std::vector<idxofflen> unk;

    std::vector<int64_t> data_pos;
    std::vector<uint64_t> varname_pos;

    Rcpp::NumericVector fmt32s;
    Rcpp::NumericVector ifmt32s;
    Rcpp::NumericVector fmtkeys;
    Rcpp::IntegerVector vartyps;
    Rcpp::IntegerVector colwidth;
    Rcpp::IntegerVector coloffset;
    Rcpp::IntegerVector page_type;

    Rcpp::CharacterVector labels;
    Rcpp::CharacterVector formats;
    Rcpp::CharacterVector varnames; // (colnum)

    // read 2* 4*4 = 32
    // 0 - 31: Magic Number
    int32_t mn1 = 0, mn2 = 0, mn3 = 0, mn4 = 0;
    int32_t mn5 = 0, mn6 = 0, mn7 = 0, mn8 = 0;

    mn1 = readbin(mn1, sas, 0);
    mn2 = readbin(mn2, sas, 0);
    mn3 = readbin(mn3, sas, 0);
    mn4 = readbin(mn4, sas, 0);
    mn5 = readbin(mn5, sas, 0);
    mn6 = readbin(mn6, sas, 0);
    mn7 = readbin(mn7, sas, 0);
    mn8 = readbin(mn8, sas, 0);

    if (debug)
      Rprintf("Magicnumber: %d, %d, %d, %d, %d, %d, %d, %d  \n",
              mn1, mn2, mn3, mn4, mn5, mn6, mn7, mn8);

    if (mn1 != 0)
      Rcpp::warning("mn1 != 0");
    // End Magicnumber

    /*
     * Most likely the following blocks of 4 are handled as fixed values. eg as
     *  int32 with a specific meaning.
     *
     */

    if (debug) Rcout << sas.tellg() << std::endl;

    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    ALIGN_1_CHECKER_VALUE = readbin(ALIGN_1_CHECKER_VALUE, sas, 0); // 51 or 34
    if (ALIGN_1_CHECKER_VALUE == 51) u64 = 4;// 51 is b'3'
    if (debug) Rprintf("ALIGN_1_CHECKER_VALUE: %d \n", ALIGN_1_CHECKER_VALUE);

    unk8 = readbin(unk8, sas, 0);
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, 0);
    if (debug) Rprintf("%d\n", unk8);
    int8_t U64_BYTE_CHECKER_VALUE = 0;
    U64_BYTE_CHECKER_VALUE = readbin(U64_BYTE_CHECKER_VALUE, sas, 0);
    if (U64_BYTE_CHECKER_VALUE == 51) ALIGN_2_VALUE = 4;
    if (debug) Rprintf("U64_BYTE_CHECKER_VALUE: %d \n", U64_BYTE_CHECKER_VALUE);
    // end block of 4

    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, 0);
    if (debug) Rprintf("%d\n", unk8);
    ENDIANNESS = readbin(ENDIANNESS, sas, 0);
    if (debug) Rprintf("ENDIANNESS: %d \n", ENDIANNESS);
    if (ENDIANNESS == 0) swapit = 1;
    // if (ENDIANNESS < 1) stop("Endiannes < 1 found");
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    // (1) 49 Unix (2) 50 Win
    uint8_t PLATFORM = 0;
    PLATFORM = readbin(PLATFORM, sas, swapit);
    if (debug) Rprintf("PLATFORM: %d \n", PLATFORM);
    // end block of 4


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // 1|4
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // ?
    if (debug) Rprintf("%d\n", unk8);


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 3
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 1
    if (debug) Rprintf("%d\n", unk8);


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);


    int8_t U64_BYTE_CHECKER_VALUE2 = 0;

    /* SAS written files repeat the first two blocks here */

    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("ALIGN_1_CHECKER_VALUE2 %d\n", unk8);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    U64_BYTE_CHECKER_VALUE2 = readbin(U64_BYTE_CHECKER_VALUE2, sas, swapit);
    if (debug) Rprintf("U64_BYTE_CHECKER_VALUE2 %d\n", U64_BYTE_CHECKER_VALUE2);

    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("%d\n", unk8);
    ENDIANNESS = readbin(ENDIANNESS, sas, swapit);
    if (debug) Rprintf("ENDIANNESS: %d \n", ENDIANNESS);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("unk2 %d\n", unk8);
    unk8 = readbin(unk8, sas, swapit);
    if (debug) Rprintf("Platform2 %d\n", unk8);

    // o64

    if (debug) Rcout << " ---- block ---- " << sas.tellg() << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // 4
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 51
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 1
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 35
    if (debug) Rprintf("%d\n", unk8);


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // interpreted as sas release ?
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // data representation 1 linux ?
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // file encoding?
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // sys encoding?
    if (debug) Rprintf("%d\n", unk8);


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 16 | 32
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 3
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 1
    if (debug) Rprintf("%d\n", unk8);


    // o80
    unk32 = readbin(unk32, sas, swapit); // 0
    // Rprintf("%d\n", unk32);

    unk32 = readbin(unk32, sas, swapit); // 0
    // Rprintf("%d\n", unk32);


    // o84 SAS FILE
    sasfile = readstring(sasfile, sas);
    if (debug) Rcout << sasfile << std::endl;

    // o92 dataset name
    dataset = readstring(dataset, sas);
    if (debug) Rcout << dataset << std::endl;

    // o156 filetype 'DATA    '
    filetype = readstring(filetype, sas);
    if (debug) Rcout << filetype << std::endl;

    if (ALIGN_2_VALUE == 4) {
      unk32 = readbin(unk32, sas, swapit);
      if (debug) Rcout << unk32 << std::endl;
    }

    created = readbin(created, sas, swapit);
    if (debug) Rcout << created << std::endl;

    modified = readbin(modified, sas, swapit);
    if (debug) Rcout << modified << std::endl;

    created2 = readbin(created2, sas, swapit);
    if (debug) Rcout << created2 << std::endl;

    modified2 = readbin(modified2, sas, swapit);
    if (debug) Rcout << modified2 << std::endl;

    int32_t headersize = 0;
    headersize = readbin(headersize, sas, swapit);
    if (debug) Rprintf("headersize: %d \n", headersize);
    if (headersize <= 0) stop("headersize <= 0");

    int32_t pagesize = 0;
    pagesize = readbin(pagesize, sas, swapit);
    if (debug) Rprintf("pagesize: %d \n", pagesize);
    if (pagesize <= 0) stop("pagesize <= 0");

    int64_t pagecount = 0;
    if (u64 == 4) {
      pagecount = readbin(pagecount, sas, swapit);
    } else {
      pagecount = readbin((int32_t)pagecount, sas, swapit);
    }
    if (debug) Rprintf("pagecount: %d \n", pagecount);

    std::vector<int64_t> rowsperpage(pagecount);

    unkdub = readbin(unkdub, sas, swapit); // 0
    if (debug) Rcout << unkdub << std::endl;

    sasrel = readstring(sasrel, sas);
    if (debug) Rcout << "SAS release: " << sasrel << std::endl;

    sasserv = readstring(sasserv, sas);
    if (debug) Rcout << "SAS server: " << sasserv << std::endl;

    // osversion
    osver = readstring(osver, sas);
    if (debug) Rcout << "OS ver: " <<  osver << std::endl;

    // osmaker
    osmaker = readstring(osmaker, sas);
    if (debug) Rcout << "OS maker: " << osmaker << std::endl; // eg WIN

    // osname
    osname = readstring(osname, sas); // x86_64
    if (debug) Rcout << "OS name: " << osname << std::endl;

    uint32_t uunk32 = 0;

    // unk
    uunk32 = readbin(uunk32, sas, swapit);
    if (debug) Rcout << uunk32 << std::endl;

    // three identical unks
    uunk32 = readbin(uunk32, sas, swapit);
    if (debug) Rcout << uunk32 << std::endl;
    uunk32 = readbin(uunk32, sas, swapit);
    if (debug) Rcout << uunk32 << std::endl;
    uunk32 = readbin(uunk32, sas, swapit);
    if (debug) Rcout << uunk32 << std::endl;

    unkdub = readbin(unkdub, sas, swapit); // 0
    if (debug) Rcout << unkdub << std::endl;
    unkdub = readbin(unkdub, sas, swapit); // 0
    if (debug) Rcout << unkdub << std::endl;

    // Rcout << sas.tellg() << std::endl;

    // page seq num at 320|328
    int32_t PAGE_SIZE = 0, PAGE_COUNT = 0;

    pageseqnum32 = readbin(pageseqnum32, sas, swapit);
    if (debug) Rcout << "pageseqnum: " << pageseqnum32 << std::endl;

    unk32 = readbin(unk32, sas, swapit); // 0 padding?
    if (debug) Rprintf("unk32: %d \n", unk32);

    // 3rd timestamp ? 0
    double thrdts = 0;
    thrdts = readbin(thrdts, sas, swapit);
    if (debug) Rcout << "3. TS " << thrdts << std::endl;

    // rest is filled with zeros. read it anyway.
    uint64_t num_zeros = headersize - sas.tellg();

    // Rcout << num_zeros << std::endl;

    for (int i = 0; i < num_zeros; ++i) {
      int8_t zero = 0;
      zero = readbin(zero, sas, swapit);
      if (zero!=0)
        warning("Error. Expected 0, is %d", zero);
    }

    // debug
    int64_t pagestart = sas.tellg();
    if (debug) Rprintf("position: %d\n", pagestart);

    // end of Header ---------------------------------------------------------//


    int8_t alignval = 8;
    if (u64 != 4) alignval = 4;

    int64_t rowlength = 0, rowcount = 0;
    int64_t colf_p1 = 0, colf_p2 = 0;
    int64_t colnum = 0;
    std::vector<std::string> stringvec(pagecount) ;

    auto totalrows = 0;
    std::vector<int32_t> totalrowsvec(pagecount);


    // begin reading pages ---------------------------------------------------//
    for (auto pg = 0; pg < pagecount; ++pg) {
      checkUserInterrupt();

      /* should already be at this position for pg == 1 */
      if (pagecount > 0) {
        auto pagenumx = headersize + pg * pagesize;
        sas.seekg(pagenumx, sas.beg);
        // Rcout << pagenumx << std::endl;
      }


      int64_t unk1 = 0, unk2 = 0, unk3 = 0;

      // Page Offset Table
      if (u64 == 4) {
        pageseqnum32 = readbin(pageseqnum32, sas, swapit);
        unk32 = readbin(unk32, sas, swapit);
        unk1 = readbin(unk1, sas, swapit);
        unk2 = readbin(unk2, sas, swapit);
        unk3 = readbin(unk3, sas, swapit);
      } else {
        pageseqnum32 = readbin(pageseqnum32, sas, swapit);
        unk1 = readbin((int32_t)unk1, sas, swapit);
        unk2 = readbin((int32_t)unk2, sas, swapit);
        unk3 = readbin((int32_t)unk3, sas, swapit);
      }

      // Rcout << "pageseqnum: " << pageseqnum32 << std::endl;
      // Rcout << unk1 << " " << unk2 << " " << unk3 << std::endl;

      PAGE_TYPE = readbin(PAGE_TYPE, sas, swapit);
      BLOCK_COUNT = readbin(BLOCK_COUNT, sas, swapit);
      SUBHEADER_COUNT = readbin(SUBHEADER_COUNT, sas, swapit);
      unk16 = readbin(unk16, sas, swapit);
      if (debug) Rprintf("unk16: %d \n", unk16);

      page_type.push_back(PAGE_TYPE);

      rowsperpage[pg] = BLOCK_COUNT - SUBHEADER_COUNT;
      totalrows += rowsperpage[pg];
      totalrowsvec[pg] = totalrows;

      if (debug)
        Rprintf("PAGE_TYPE: %d ; BLOCK_COUNT: %d ; SUBHEADER_COUNT: %d ---- \n",
                PAGE_TYPE, BLOCK_COUNT, SUBHEADER_COUNT);


      int16_t zero = 0;
      int64_t sh_tot_len = 0;
      uint64_t dataoff = 0;

      std::vector<PO_Tab> potabs(SUBHEADER_COUNT);


      if (( PAGE_TYPE == 1024 || PAGE_TYPE == 640 || PAGE_TYPE == 512 ||
          PAGE_TYPE == 256 || PAGE_TYPE == 0))
      {
        for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
          if (u64 == 4) {

            potabs[i].SH_OFF = readbin(potabs[i].SH_OFF, sas, swapit);
            potabs[i].SH_LEN = readbin(potabs[i].SH_LEN, sas, swapit);
            potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, swapit);
            potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, swapit);

            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;
            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;
            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;

            if (debug)
              Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPR.: %d ; SH_TYPE: %d \n",
                      potabs[i].SH_OFF, potabs[i].SH_LEN,
                      potabs[i].COMPRESSION, potabs[i].SH_TYPE);

          } else {

            potabs[i].SH_OFF = readbin((int32_t)potabs[i].SH_OFF, sas, swapit);
            potabs[i].SH_LEN = readbin((int32_t)potabs[i].SH_LEN, sas, swapit);
            potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, swapit);
            potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, swapit);

            zero = readbin(zero, sas, swapit);

            if (debug)
              Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPR.: %d ; SH_TYPE: %d \n",
                      potabs[i].SH_OFF, potabs[i].SH_LEN,
                      potabs[i].COMPRESSION, potabs[i].SH_TYPE);
          }

          sh_tot_len += potabs[i].SH_LEN;
          dataoff = potabs[i].SH_OFF - (rowsperpage[pg] * rowlength);
        }

        if (debug) {
          Rcout << "data offset " << dataoff << std::endl;
          Rcout << "data offset ------------------------------- " << std::endl;
        }


        uint64_t pos = sas.tellg();


        auto sh_end_pos = 0;

        if (PAGE_TYPE != 0) sh_end_pos = sas.tellg();
        data_pos.push_back( sh_end_pos );

        if (debug)
          Rprintf("sh_end_pos: %d\n", sh_end_pos);


        auto pg_vars = 0;


        // from now on, we will seek to every position inside the sas file
        for (auto sc = 0; sc < SUBHEADER_COUNT; ++sc)
        {

          auto pagepos = (headersize + pg * pagesize) + potabs[sc].SH_OFF;
          sas.seekg(pagepos, sas.beg);


          int64_t sas_offset = alignval;
          if (u64 == 4) {
            sas_offset = readbin(sas_offset, sas, swapit);
          } else {
            sas_offset = readbin((int32_t)sas_offset, sas, swapit);
          }

          std::string sas_hex = int_to_hex(sas_offset);

          auto sas_offset_table = 0;
          if (sas_hex.compare("f7f7f7f7") == 0 ||
              sas_hex.compare("fffffffff7f7f7f7") == 0 ||
              sas_hex.compare("f7f7f7f700000000") == 0 ||
              sas_hex.compare("f7f7f7f7fffffbfe") == 0 )
            sas_offset_table = 1;
          if (sas_hex.compare("fffffc00") == 0 ||
              sas_hex.compare("fffffffffffffc00") == 0 )
            sas_offset_table = 2;
          if (sas_hex.compare("fffffbfe") == 0 ||
              sas_hex.compare("fffffffffffffbfe") == 0)
            sas_offset_table = 3;
          if (sas_hex.compare("f6f6f6f6") == 0 ||
              sas_hex.compare("fffffffff6f6f6f6") == 0 ||
              sas_hex.compare("f6f6f6f600000000") == 0 ||
              sas_hex.compare("f6f6f6f6fffffbfe") == 0 )
            sas_offset_table = 4;
          if (sas_hex.compare("fffffffd") == 0 ||
              sas_hex.compare("fffffffffffffffd") == 0)
            sas_offset_table = 5;
          if (sas_hex.compare("ffffffff") == 0 ||
              sas_hex.compare("ffffffffffffffff") == 0)
            sas_offset_table = 6;
          if (sas_hex.compare("fffffffc") == 0 ||
              sas_hex.compare("fffffffffffffffc") == 0)
            sas_offset_table = 7;
          if (sas_hex.compare("fffffffe") == 0 ||
              sas_hex.compare("fffffffffffffffe") == 0)
            sas_offset_table = 8;
          if (sas_hex.compare("0") == 0)
            sas_offset_table = 9;

          switch(sas_offset_table)
          {


            // new offset --------------------------------------------------- //
          case 1:
            {
              /* Row Size */

              int16_t pgwpossh = 0, pgwpossh2 = 0, numzeros = 37,
                sh_num = 0, cn_maxlen = 0, l_maxlen = 0,
                rowsonpg = 0;
              int32_t pgidx = 0;
              int64_t pgsize = 0, pgc = 0, rcmix = 0, pgwsh = 0, pgwsh2 = 0;


              // Rcout << "-------- case 1 "<< sas.tellg() << std::endl;



              if (u64 == 4) {

                /* */


                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;

                rowlength = readbin(rowlength, sas, swapit);
                if (debug) Rcout << "rowlength " << rowlength << std::endl;
                rowcount = readbin(rowcount, sas, swapit);
                if (debug) Rcout << "rowcount " << rowcount << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;

                colf_p1 = readbin(colf_p1, sas, swapit);
                if (debug) Rcout << colf_p1 << std::endl;
                colf_p2 = readbin(colf_p2, sas, swapit);
                if (debug) Rcout << colf_p2 << std::endl;
                unk64 = readbin(unk64, sas, swapit); // p3 and p4?
                if (debug) Rcout << unk64 << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl;
                pgsize = readbin(pgsize, sas, swapit);
                unk64 = readbin(unk64, sas, swapit);
                rcmix =  readbin(rcmix, sas, swapit);

                uint64_t uunk64 = 0;
                /* next two indicate the end of the initial header ? */
                uunk64 = readbin(uunk64, sas, swapit);
                uunk64 = readbin(uunk64, sas, swapit);

                for (int z = 0; z < numzeros; ++z) {
                  unk64 = readbin(unk64, sas, swapit);
                  if (unk64 != 0)
                    warning("val is %d. expected a zero", unk64);
                }

                pgidx = readbin(pgidx, sas, swapit);

                // padding? 68 bytes: zeros
                for (int z = 0; z < 8; ++z) {
                  unk64 = readbin(unk64, sas, swapit);
                  if (unk64 != 0)
                    warning("val0 is %d. expected a zero", unk64);
                }
                unk32 = readbin(unk32, sas, swapit);

                unk64 = readbin(unk64, sas, swapit); // val 1?
                unk16 = readbin(unk16, sas, swapit); // val 2?

                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding

                pgwsh = readbin(pgwsh, sas, swapit);
                pgwpossh = readbin(pgwpossh, sas, swapit);

                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding

                pgwsh2 = readbin(pgwsh2, sas, swapit);
                pgwpossh2 = readbin(pgwpossh2, sas, swapit);

                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding

                pgc = readbin(pgc, sas, swapit);

                unk16 = readbin(unk16, sas, swapit); // val ?
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding

                unk64 = readbin(unk64, sas, swapit); // val 1?

                addtextoff = readbin(addtextoff, sas, swapit); // val 7 | 8?
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding
                unk16 = readbin(unk16, sas, swapit); // padding

                for (int z = 0; z < 10; ++z) {
                  unk64 = readbin(unk64, sas, swapit); // 0
                  if (unk64 != 0)
                    warning("val1 is %d. expected a zero", unk64);
                }

                if (debug) Rcout << "###############" << std::endl;

                unk16 = readbin(unk16, sas, swapit); // val
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 0|8 ?
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 4
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 0
                if (debug) Rcout << unk16 << std::endl;
                todata = readbin(unk16, sas, swapit); // val 12,32|0?
                if (debug) Rcout << todata << std::endl;

                if (todata == 12)
                  hasproc = false;

                if (debug) Rcout << "###############" << std::endl;

                swlen = readbin(swlen, sas, swapit);
                if (debug) Rcout << swlen << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 0
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 20 | 28
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 8
                if (debug) Rcout << unk16 << std::endl;
                // unk64 = readbin(unk64, sas, swapit); // ?
                // Rcout << unk64 << std::endl;

                if (debug) Rcout << "###############" << std::endl;

                unk16 = readbin(unk16, sas, swapit); // 0
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 12
                if (debug) Rcout << unk16 << std::endl;
                comprlen = readbin(unk16, sas, swapit); // 8
                if (debug) Rcout << comprlen << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 0
                if (debug) Rcout << unk16 << std::endl;

                if (debug) Rcout << "###############" << std::endl;

                unk16 = readbin(unk16, sas, swapit); // 12 | 20
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 8
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 0
                if (debug) Rcout << unk16 << std::endl;
                textoff = readbin(textoff, sas, swapit); // 28 | 36
                if (debug) Rcout << textoff << std::endl;
                proclen = readbin(proclen, sas, swapit);
                if (debug) Rcout << "proclen " << proclen << std::endl;

                if (debug) Rcout << "###############" << std::endl;

                for (int z = 0; z < 8; ++z) {
                  unk32 = readbin(unk32, sas, swapit); // 0
                  if (unk64 != 0)
                    warning("val2 is %d. expected a zero", unk64);
                }

                unk16 = readbin(unk16, sas, swapit); // 4
                unk16 = readbin(unk16, sas, swapit); // 1

                sh_num = readbin(sh_num, sas, swapit);
                cn_maxlen = readbin(cn_maxlen, sas, swapit);
                l_maxlen = readbin(l_maxlen, sas, swapit);

                /* maybe SAS version information at o131018 ? */
                for (int z = 0; z < 3; ++z) {
                  unk32 = readbin(unk32, sas, swapit); // 0
                  if (unk32 != 0)
                    warning("val3 is %d at %d. expected a zero",
                            unk32, sas.tellg());
                }

                rowsonpg = readbin(rowsonpg, sas, swapit);

                int datofs = 0;
                for (int z = 0; z < 20; ++z) {
                  unk16 = readbin(unk16, sas, swapit); // 0
                  if (unk16 != 0) {
                    warning("val4 is %d at %d. expected a zero",
                            unk16, sas.tellg());
                    datofs = unk16;
                  }
                }

                dataoffset = datofs;

                /* */

              } else {
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;

                rowlength = readbin((int32_t)rowlength, sas, swapit);
                if (debug) Rcout << "rowlength "<< rowlength << std::endl;
                rowcount = readbin((int32_t)rowcount, sas, swapit);
                if (debug) Rcout << "rowcount "<< rowcount << std::endl;
                unk32 = readbin(unk32, sas, swapit); // deleted variables?
                if (debug) Rcout << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                colf_p1 = readbin((int32_t)colf_p1, sas, swapit);
                if (debug) Rcout << "colfp1 " << colf_p1 << std::endl;
                colf_p2 = readbin((int32_t)colf_p2, sas, swapit);
                if (debug) Rcout << "colfp2 " << colf_p2 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                pgsize = readbin((int32_t)pgsize, sas, swapit);
                if (debug) Rcout << "pgsize " << pgsize << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                rcmix =  readbin((int32_t)rcmix, sas, swapit);
                if (debug) Rcout << "rcmix " << rcmix << std::endl;
                uunk32 = readbin(uunk32, sas, swapit);
                if (debug) Rcout << uunk32 << std::endl;
                uunk32 = readbin(uunk32, sas, swapit);
                if (debug) Rcout << uunk32 << std::endl;

                for (int z = 0; z < numzeros; ++z) {
                  unk64 = readbin((int32_t)unk64, sas, swapit);
                  if (unk64 != 0)
                    warning("val1 is %d. expected a zero", unk64);
                }

                pgidx = readbin(pgidx, sas, swapit);


                for (int z = 0; z < 8; ++z) {
                  unk64 = readbin((int32_t)unk64, sas, swapit);
                  if (debug) Rcout << unk64 << std::endl;
                }

                // padding?
                unk32 = readbin(unk32, sas, swapit);
                unk32 = readbin(unk32, sas, swapit);

                unk32 = readbin(unk32, sas, swapit); // val 1?
                if (debug) Rcout << unk32 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 2?
                if (debug) Rcout << unk16 << std::endl;

                unk16 = readbin(unk16, sas, swapit); // padding

                pgwsh = readbin((int32_t)pgwsh, sas, swapit);
                if (debug) Rcout << "pgwsh " << pgwsh << std::endl;
                pgwpossh = readbin(pgwpossh, sas, swapit);
                if (debug) Rcout << "pgwpossh " << pgwpossh << std::endl;

                unk16 = readbin(unk16, sas, swapit); // padding

                pgwsh2 = readbin((int32_t)pgwsh2, sas, swapit);
                if (debug) Rcout << "pgwsh2 " << pgwsh2 << std::endl;
                pgwpossh2 = readbin(pgwpossh2, sas, swapit);
                if (debug) Rcout << "pgwpossh2 " << pgwpossh2 << std::endl;

                unk16 = readbin(unk16, sas, swapit); // padding

                pgc = readbin((int32_t)pgc, sas, swapit);
                if (debug) Rcout << "pgc " << pgc << std::endl;

                unk16 = readbin(unk16, sas, swapit); // val ?
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // padding

                unk64 = readbin((int32_t)unk64, sas, swapit); // val 1?
                if (debug) Rcout << unk64 << std::endl;

                addtextoff = readbin(addtextoff, sas, swapit); // val 7 | 8?
                if (debug) Rcout << addtextoff << std::endl;
                unk16 = readbin(unk16, sas, swapit); // padding

                for (int z = 0; z < 10; ++z) {
                  unk64 = readbin((int32_t)unk64, sas, swapit); // 0
                  if (unk64 != 0)
                    warning("val2 is %d. expected a zero", unk64);
                }

                if (debug) Rcout << "###############" << std::endl;

                unk16 = readbin(unk16, sas, swapit); // val
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 0|8 ?
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 4
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 0
                if (debug) Rcout << unk16 << std::endl;
                todata = readbin(todata, sas, swapit); // val 12,32|0? //
                if (debug) Rcout << todata << std::endl;

                if (todata == 12)
                  hasproc = false;

                if (debug) Rcout << "###############" << std::endl;

                swlen = readbin(swlen, sas, swapit);
                if (debug) Rcout << "swlen " << swlen << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 0?
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 20?
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); //
                if (debug) Rcout << unk16 << std::endl;

                if (debug) Rcout << "###############" << std::endl;

                unk16 = readbin(unk16, sas, swapit); // 0
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 12
                if (debug) Rcout << unk16 << std::endl;
                comprlen = readbin(unk16, sas, swapit); // 8 compr. code length?
                if (debug) Rcout << comprlen << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 0
                if (debug) Rcout << unk16 << std::endl;

                if (debug) Rcout << "###############" << std::endl;

                unk16 = readbin(unk16, sas, swapit); // 12
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 8
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 0
                if (debug) Rcout << unk16 << std::endl;
                textoff = readbin(textoff, sas, swapit); // 28
                if (debug) Rcout << textoff << std::endl;
                proclen = readbin(proclen, sas, swapit);
                if (debug) Rcout << "proclen " << proclen << std::endl;

                if (debug) Rcout << "###############" << std::endl;

                for (int z = 0; z < 8; ++z) {
                  unk32 = readbin(unk32, sas, swapit); // 0
                  if (unk32 != 0)
                    warning("val3 is %d. expected a zero", unk32);
                }

                unk16 = readbin(unk16, sas, swapit); // 4
                unk16 = readbin(unk16, sas, swapit); // 1

                sh_num = readbin(sh_num, sas, swapit);
                cn_maxlen = readbin(cn_maxlen, sas, swapit);
                l_maxlen = readbin(l_maxlen, sas, swapit);

                for (int z = 0; z < 3; ++z) {
                  unk32 = readbin(unk32, sas, swapit); // 0
                }

                rowsonpg = readbin(rowsonpg, sas, swapit);

                int datofs = 0;
                for (int z = 0; z < 20; ++z) {
                  unk16 = readbin(unk16, sas, swapit); // 0
                  if (unk16 != 0) {
                    warning("val4 is %d at %d. expected a zero",
                            unk16, sas.tellg());
                    datofs = unk16;
                  }
                }

                dataoffset = datofs;
              }

              if (debug)
                Rprintf("swlen = %d, todata %d, textoff %d\n",
                        swlen, todata, textoff);


              if (!((dataoffset == 256) | (dataoffset == 1280)))
                warning("dataoffset is unexpectedly %d\n",
                        dataoffset);

              break;
            }



            // new offset --------------------------------------------------- //
          case 2:
            {

              // Rcout << "-------- case 2 "<< sas.tellg() << std::endl;

              int64_t off = 0;

              if (u64 == 4) {
                off = readbin(off, sas, swapit);
                // Rcout << off << std::endl;
                unk64 = readbin(unk64, sas, swapit);
                // Rcout << unk64 << std::endl;
              } else {
                off = readbin((int32_t)off, sas, swapit);
                // Rcout << off << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                // Rcout << unk32 << std::endl;
              }

              int16_t num_nonzero = 0;
              num_nonzero = readbin(num_nonzero, sas, swapit);
              // Rcout << "numzeros " <<  num_nonzero << std::endl;

              int8_t unklen = 94; // should be 94
              if (u64 != 4) unklen = 50;
              for (int jj = 0; jj < unklen/2; ++jj) {
                unk16 = readbin(unk16, sas, swapit);
                // 4th from the end is 1804 meaning is unknown
              }

              std::vector<SCV> scv(12);

              for (int8_t i = 0; i < 12; ++i) {

                if (u64 == 4) {
                  scv[i].SIG = readbin(scv[i].SIG, sas, swapit);
                  scv[i].FIRST = readbin(scv[i].FIRST, sas, swapit);
                  scv[i].F_POS = readbin(scv[i].F_POS, sas, swapit);

                  if ((i == 0) & (scv[i].SIG != -4))
                    warning("first SIG is not -4");

                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;

                  scv[i].LAST = readbin(scv[i].LAST, sas, swapit);
                  scv[i].L_POS = readbin(scv[i].L_POS, sas, swapit);

                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;

                } else {
                  scv[i].SIG = readbin((int32_t)scv[i].SIG, sas, swapit);
                  scv[i].FIRST = readbin((int32_t)scv[i].FIRST, sas, swapit);
                  scv[i].F_POS = readbin(scv[i].F_POS, sas, swapit);

                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;

                  scv[i].LAST = readbin((int32_t)scv[i].LAST, sas, swapit);
                  scv[i].L_POS = readbin(scv[i].L_POS, sas, swapit);

                  unk16 = readbin(unk16, sas, swapit);
                  // Rcout << unk16 << std::endl;
                }

                if (debug)
                  Rprintf("SIG %d; FIRST %d; F_POS %d; LAST %d; L_POS %d\n",
                          scv[i].SIG, scv[i].FIRST, scv[i].F_POS,
                          scv[i].LAST, scv[i].L_POS);

              }

              break;
            }


            // new offset --------------------------------------------------- //
          case 3:
            {

              // Rcout << "-------- case 3 "<< sas.tellg() << std::endl;

              hasattributes = 1;

              int64_t unk64 = 0;

              unk16 = readbin(unk16, sas, swapit);           // 1
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 2
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 3
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 4
              // Rcout << unk16 << std::endl;
              fmt32 = readbin(fmt32, sas, swapit);           // 5
              // Rcout << fmt32 << std::endl;
              fmt322 = readbin(fmt322, sas, swapit);         // 6
              // Rcout << fmt322 << std::endl;
              ifmt32 = readbin(ifmt32, sas, swapit);         // 7
              // Rcout << ifmt32 << std::endl;
              ifmt322 = readbin(ifmt322, sas, swapit);       // 8
              // Rcout << ifmt322 << std::endl;
              fmtkey = readbin(fmtkey, sas, swapit);         // 9
              // Rcout << fmtkey << std::endl;
              fmtkey2 = readbin(fmtkey2, sas, swapit);       // 10
              // Rcout << fmtkey2 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 11
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 12
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 13
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 14 off + len
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);           // 15 1 w char
              // Rcout << unk16 << std::endl;

              if (u64 == 4) {
                unk16 = readbin(unk16, sas, swapit);
                // Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit);
                // Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit);
                // Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit);
                // Rcout << unk16 << std::endl;
              }

              // Rcout << "-------- case 3 "<< sas.tellg()  << std::endl;

              fmt32s.push_back(  fmt32  + (double)fmt322/10);
              ifmt32s.push_back( ifmt32 + (double)ifmt322/10);
              fmtkeys.push_back( fmtkey + (double)fmtkey2/10);

              idxofflen fmts, lbls, unks;

              fmts.IDX = readbin(fmts.IDX, sas, swapit);
              fmts.OFF = readbin(fmts.OFF, sas, swapit);
              fmts.LEN = readbin(fmts.LEN, sas, swapit);

              if (debug)
                Rcout << fmts.IDX << ", " << fmts.OFF <<
                  ", " << fmts.LEN << std::endl;

              fmt.push_back(fmts);

              lbls.IDX = readbin(lbls.IDX, sas, swapit);
              lbls.OFF = readbin(lbls.OFF, sas, swapit);
              lbls.LEN = readbin(lbls.LEN, sas, swapit);

              if (debug)
                Rcout << lbls.IDX << ", " << lbls.OFF <<
                  ", " << lbls.LEN << std::endl;

              lbl.push_back(lbls);

              unks.IDX = readbin(unks.IDX, sas, swapit);
              unks.OFF = readbin(unks.OFF, sas, swapit);
              unks.LEN = readbin(unks.LEN, sas, swapit);


              unk.push_back(unks);

              if (unks.IDX != 0 | unks.OFF != 0 | unks.LEN != 0) {
                warning("case3: unk is not 0 as expected, but %d %d %d\n",
                        unks.IDX, unks.OFF, unks.LEN);
                // Rcout << unks.IDX << ", " << unks.OFF <<
                //   ", " << unks.LEN << std::endl;
              }

              break;
            }


            // new offset --------------------------------------------------- //
          case 4:
            {
              /* Column Size */

              // Rcout << "-------- case 4 "<< sas.tellg() << std::endl;

              uint64_t uunk64 = 0;

              if (u64 == 4) {
                colnum = readbin(colnum, sas, swapit);
                uunk64 = readbin(uunk64, sas, swapit);
              } else {
                colnum = readbin((int32_t)colnum, sas, swapit);
                uunk64 = readbin((int32_t)uunk64, sas, swapit);
              }

              if (debug)
                Rprintf("colnum %d; uunk64 %d\n",
                        colnum, uunk64);

              break;
            }

            // new offset --------------------------------------------------- //

          case 5:
            {
              /* Column Text */

              // Rcout << "-------- case 5 "<< sas.tellg() << std::endl;

              int16_t len = 0;

              varname_pos.push_back( sas.tellg() );

              len = readbin(len, sas, swapit);
              // Rcout << len << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0 vars on p1 ?
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0 vars on p2 ?
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              // not sure yet. pg == 0 is required
              if ((PAGE_TYPE != 1024) & (c5first == 0)) {
                unk16 = readbin(unk16, sas, swapit); // 0
                // Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 0
                // Rcout << unk16 << std::endl;
              }

              if (debug)
                Rprintf("SH_LEN %d; len %d; newlen: %d\n",
                        potabs[sc].SH_LEN, len);

              c5first = 1;

              break;
            }


            // new offset --------------------------------------------------- //
          case 6:
            {
              /* Column Name */

              // Rcout << "-------- case 6 "<< sas.tellg() << std::endl;


              int16_t lenremain = 0;
              int8_t tmp = 16; // always 16?
              if (!hasproc) tmp = 0;

              lenremain = readbin(lenremain, sas, swapit);
              if (debug) Rprintf("lenremain %d \n", lenremain);
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              if (varname_pos.size() > 0) {

                if (!(PAGE_TYPE == 1024)) {

                  uint64_t pos_beg = sas.tellg();
                  int8_t tmp = 16; // always 16?
                  if (!hasproc) tmp = 0;

                  if (debug)
                    Rprintf("%d, %d, %d, %d; %d\n",
                            comprlen, tmp, proclen, swlen, varname_pos[0]);

                  uint64_t txtpos = varname_pos[0] + 12;

                  sas.seekg(txtpos, sas.beg);

                  // compression
                  if (comprlen > 0) {
                    compression.resize(comprlen, '\0');
                    compression = readstring(compression, sas);

                    if (compression.compare("SASYZCRL") == 0)
                      compr = 1;

                    if (compression.compare("SASYZCR2") == 0)
                      compr = 2;
                    // Rcout << compression << std::endl;
                  }

                  // 16 whitespaces
                  std::string empty (tmp, '\0');
                  if (tmp > 0) {
                    empty = readstring(empty, sas);
                    if (!(empty.compare("                ") == 0))
                      warning("non empty 'empty' string found %s \n",
                              empty);
                  }

                  // proc that created the file
                  // if (proclen > 0) {
                  //   proc.resize(proclen, '\0');
                  //   proc = readstring(proc, sas);
                  // }

                  // additional software string
                  if (swlen > 0) {
                    sw.resize(swlen, '\0');
                    sw = readstring(sw, sas);
                  }


                  if (debug)
                    Rcout << "here we go!\n" <<
                      compression << "\n" <<
                        empty << "\n" <<
                          proc << "\n" <<
                            sw << std::endl;

                  // Rprintf("compr %d \n", compr  );

                  sas.seekg(pos_beg, sas.beg);
                }


                // stop("stop");


                /* Column Name Pointers */
                auto cmax = (lenremain + alignval)/8;

                // Rcout << cmax << std::endl;

                std::vector<CN_Poi> cnpois(cmax);

                for (auto i = 0; i < cmax; ++i) {

                  auto idx    = readbin(cnpois[i].CN_IDX, sas, swapit);
                  auto off    = readbin(cnpois[i].CN_OFF, sas, swapit);
                  auto len    = readbin(cnpois[i].CN_LEN, sas, swapit);
                  auto zeros  = readbin(cnpois[i].zeros,  sas, swapit);

                  pos = sas.tellg();


                  if (len > 0) {

                    if (debug)
                      Rprintf("CN_IDX %d; CN_OFF %d; CN_LEN %d; zeros %d \n",
                              idx, off, len, zeros);

                    int64_t vpos = (varname_pos[pg_vars] + off);
                    sas.seekg(vpos, sas.beg);

                    std::string varname(len, '\0');
                    varname = readstring(varname, sas);
                    varnames.push_back(varname);

                    // Rcout << vpos << " : " << varname << std::endl;

                    sas.seekg(pos, sas.beg);

                  }

                }

                ++pg_vars;
              }

              break;
            }


            // new offset --------------------------------------------------- //
          case 7:
            {
              /* Column Attributes */

              // Rcout << "-------- case 7 "<< sas.tellg()  << std::endl;

              int16_t lenremain = 0;
              lenremain = readbin(lenremain, sas, swapit);
              if (debug) Rprintf("lenremain %d \n", lenremain);

              // zeros as padding?
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              // auto cmax = (lenremain + alignval) / (alignval);
              // if (debug) Rcout << cmax << std::endl;

              auto cmax = colf_p1; if (pg == 2) cmax = colf_p2;

              /* Column Attributes Pointers */
              std::vector<CN_Att> capois(cmax);

              for (auto i = 0; i < cmax; ++i) {

                if (u64 == 4) {
                  capois[i].CN_OFF     = readbin(capois[i].CN_OFF, sas, swapit);
                } else {
                  capois[i].CN_OFF     = readbin((int32_t)capois[i].CN_OFF,
                                                 sas, swapit);
                }
                capois[i].CN_WID     = readbin(capois[i].CN_WID, sas, swapit);
                capois[i].NM_FLAG    = readbin(capois[i].NM_FLAG, sas, swapit);
                capois[i].CN_TYP     = readbin(capois[i].CN_TYP, sas, swapit);
                capois[i].UNK8       = readbin(capois[i].UNK8, sas, swapit);


                if (debug)
                  Rprintf("OFF %d; WID: %d; FLAG %d; TYP %d; UNK8 %d\n",
                          capois[i].CN_OFF, capois[i].CN_WID, capois[i].NM_FLAG,
                          capois[i].CN_TYP, capois[i].UNK8 );

                if ((capois[i].CN_TYP >= 1) & (capois[i].CN_TYP <= 2)) {
                  coloffset.push_back( capois[i].CN_OFF );
                  colwidth.push_back( capois[i].CN_WID );
                  vartyps.push_back( capois[i].CN_TYP );
                }

              }

              break;
            }

          case 8:
            {

              // Rcout << "-------- case 8 "<< sas.tellg() << std::endl;

              int16_t cls = 0;

              unk32 = readbin(unk32, sas, swapit); // ?
              // Rcout << unk32 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // padding? 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // padding? 0
              // Rcout << unk16 << std::endl;

              if (u64 == 4) {  // lenremain
                unk64 = readbin(unk64, sas, swapit);
              } else {
                unk64 = readbin((int32_t)unk64, sas, swapit); //
              }

              if (debug)
                Rcout << "lenremain "<< unk64 << std::endl; // 92

              unk16 = readbin(unk16, sas, swapit);  // 25
              // Rcout << unk16 << std::endl;
              cls = readbin(cls, sas, swapit);      // 37
              // Rcout << cls << std::endl;
              unk16 = readbin(unk16, sas, swapit);  // 1
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);  // 25
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);  // 3233
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);  // 3233
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit);  // 3233
              // Rcout << unk16 << std::endl;

              // 2 * CL
              for (auto cl = 0; cl < cls; ++cl) {
                // Rcout << "------" << std::endl;
                unk1 = readbin((int8_t)unk1, sas, swapit); //
                unk2 = readbin((int8_t)unk2, sas, swapit); //
                if (debug) Rcout << unk1 <<  " : " << unk2 << std::endl;
              }


              // 8
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              break;

            }

          case 9:
            {

              /* not sure about the purpose. appears in a few SAS 9 files */
              auto unklen = potabs[sc].SH_LEN - alignval;
              for (auto ul = 0; ul < unklen; ++ul) {
                unk8 = readbin(unk8, sas, 0);
                if (debug) Rprintf("%d\n", unk8);
              }

              break;
            }

            // not implemented ---------------------------------------------- //
          default:
            {

              // else it is padding?
              if ((potabs[sc].SH_LEN > alignval) &
                  (potabs[sc].COMPRESSION != 4))
            {
              Rcout << "---- unimplemented "<< sas.tellg() << std::endl;
              Rcout << "SAS HEX STRING: "  << sas_hex << std::endl;

              auto unklen = potabs[sc].SH_LEN - alignval;
              Rcout << "unklen is " << unklen << std::endl;


              sas.seekg(unklen, sas.cur);

            }

              break;
            }

          }
        }

      } else{
        Rcout << "found unimplemented PAGE_TYPE " << PAGE_TYPE << std::endl;
      }
    }


    if (hasattributes) {

      auto len = fmt.size();

      for (auto i = 0; i < len; ++i) {

        /* read formats and labels */
        std::string format = "";
        if (fmt[i].LEN > 0) {
          uint64_t fpos = (varname_pos[fmt[i].IDX] + fmt[i].OFF);
          sas.seekg(fpos, sas.beg);
          format.resize(fmt[i].LEN, '\0');
          format = readstring(format, sas);
        }

        std:: string label = "";
        if (lbl[i].LEN > 0) {
          uint64_t lpos = (varname_pos[lbl[i].IDX] + lbl[i].OFF);
          sas.seekg(lpos, sas.beg);
          label.resize(lbl[i].LEN, '\0');
          label = readstring(label, sas);
        }

        if (debug)
          Rcout << format << " : " << label << std::endl;

        formats.push_back( format );
        labels.push_back( label );
      }
    }

    if (compr != 0)
      warning("File contains unhandled compression. No data read. %d\n",
              compression);


    // ---------------------------------------------------------------------- //

    // 1. Create Rcpp::List
    Rcpp::List df(colnum);
    for (auto i=0; i<colnum; ++i)
    {
      int32_t const type = vartyps[i];

      // Rcout << type << std::endl;

      switch(type)
      {
      case 1:
        SET_VECTOR_ELT(df, i, NumericVector(no_init(rowcount)));
        break;

      default:
        SET_VECTOR_ELT(df, i, CharacterVector(no_init(rowcount)));
      break;
      }
    }


    // new offset ----------------------------------------------------------- //

    if (compr == 0) {

      auto page = 0;
      sas.seekg(data_pos[0], sas.beg);

      auto ii = 0;
      for (auto i = 0; i < rowcount; ++i) {

        if (pagecount>0) {
          while (totalrowsvec[page] == 0) {
            ++page;
            ii = 0;

            if (page == pagecount)
              break;
          }

          if (totalrowsvec[page] == i) {
            ++page;
            ii = 0;
          }
        }

        // Rcout << page << std::endl;
        // Rcout << i << " " << ii << std::endl;
        // Rcout << totalrowsvec[page] << std::endl;

        auto pp = data_pos[page];
        auto pos = pp + rowlength * ii; /* + proclen; */


        // auto pt = page_type[page];

        // Rcout << pt << std::endl;
        // Rcout << i << " " << dataoffset << " " << (int)alignval << std::endl;

        /* unknown */
        if (!(dataoffset == 256) & (page == 0)) {
          pos += alignval;
        }

        // if (page_type[page] == 256)
        //   pos -= 4;


        // if (!(dataoffset == 1280) & (page == 1))

        sas.seekg(pos, sas.beg);

        uint64_t tempoffs = sas.tellg();

        pos = 0;

        for (auto j = 0; j < colnum; ++j) {

          auto wid = colwidth[j];
          auto typ = vartyps[j];


          // Rcout  << varnames[j] << std::endl;

          uint64_t off = coloffset[j];
          uint64_t readpos = tempoffs + off;

          if (wid < 8 & typ == 1) {

            double val_d = 0.0;

            sas.seekg(readpos, sas.beg);

            val_d = readbinlen(val_d, sas, 0, wid);

            // Rcout << val_d << std::endl;

            if (std::isnan(val_d))
              REAL(VECTOR_ELT(df,j))[i] = NA_REAL;
            else
              REAL(VECTOR_ELT(df,j))[i] = val_d;

          }

          if (wid == 8 & typ == 1) {

            double val_d = 0.0;

            sas.seekg(readpos, sas.beg);

            val_d = readbin(val_d, sas, swapit);

            // Rcout << val_d << std::endl;

            if (std::isnan(val_d))
              REAL(VECTOR_ELT(df,j))[i] = NA_REAL;
            else
              REAL(VECTOR_ELT(df,j))[i] = val_d;

          }

          if (wid > 0 & typ == 2) {

            std::string val_s(wid, ' ');

            sas.seekg(readpos, sas.beg);

            val_s = readstring(val_s, sas);

            val_s = std::regex_replace(val_s,
                                       std::regex(" +$"), "$1");

            // Rcout << val_s << std::endl;

            as<CharacterVector>(df[j])[i] = val_s;

          }

        }

        ++ii;
      }
    }


    if (debug)
      Rprintf("%d %d \n", rowcount, colnum);


    Rcpp::IntegerVector rvec = seq(1, rowcount);
    Rcpp::IntegerVector cvec = seq(1, colnum);

    // 3. Create a data.frame
    df.attr("row.names") = rvec;
    if (compr == 0)
      df.attr("names") = varnames;
    else
      df.attr("names") = cvec;
    df.attr("class") = "data.frame";

    // close file
    sas.close();

    df.attr("varnames") = varnames;
    df.attr("labels") = labels;
    df.attr("formats") = formats;
    df.attr("created") = created;
    df.attr("created2") = created2;
    df.attr("modified") = modified;
    df.attr("modified2") = modified2;
    df.attr("thrdts") = thrdts;

    df.attr("sasfile") = sasfile;
    df.attr("dataset") = dataset;
    df.attr("filetype") = filetype;
    df.attr("compression") = compression;
    df.attr("proc") = proc;
    df.attr("sw") = sw;
    df.attr("sasrel") = sasrel;
    df.attr("sasserv") = sasserv;
    df.attr("osver") = osver;
    df.attr("osmaker") = osmaker;
    df.attr("osname") = osname;
    df.attr("fmtkeys") = fmtkeys;
    df.attr("fmt32") = fmt32s;
    df.attr("ifmt32") = ifmt32s;

    df.attr("colwidth") = colwidth;
    df.attr("coloffset") = coloffset;
    df.attr("vartyps") = vartyps;

    return(df);

  } else {
    return (-1);
  }
}
