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

    bool hasattributes = 0, hasproc = 1;
    int8_t  unk8  = 0;
    int16_t unk16 = 0;
    int32_t unk32 = 0;
    int64_t unk64 = 0;
    double unkdub = 0;

    int16_t swlen = 0, proclen = 0, comprlen = 0, textoff = 0, todata = 0,
      addtextoff = 0;
    int16_t PAGE_TYPE = 0, BLOCK_COUNT = 0, SUBHEADER_COUNT = 0;

    std::vector<idxofflen> fmt;
    std::vector<idxofflen> lbl;
    std::vector<idxofflen> unk;

    Rcpp::IntegerVector vartyps;
    Rcpp::IntegerVector colwidth;
    Rcpp::IntegerVector coloffset;
    Rcpp::IntegerVector data_pos;
    Rcpp::CharacterVector varnames; // (colnum)
    Rcpp::CharacterVector formats;
    Rcpp::CharacterVector labels;

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
      Rcpp::stop("mn1 != 0");
    // End Magicnumber

    // created block ?
    // o32 (char?)
    uint8_t ALIGN_1_CHECKER_VALUE = 0;
    ALIGN_1_CHECKER_VALUE = readbin(ALIGN_1_CHECKER_VALUE, sas, 0); // 51 or 34

    int8_t ALIGN_1_VALUE = 0;
    if (ALIGN_1_CHECKER_VALUE == 51) { // 51 is b'3'
      ALIGN_1_VALUE = 4;
    }

    if (debug) Rprintf("ALIGN_1_CHECKER_VALUE: %d \n", ALIGN_1_CHECKER_VALUE);

    readbin(unk8, sas, 0);
    unk8 = readbin(unk8, sas, 0);

    if (unk8 != 0)
      Rcpp::warning("Expected 0");

    // o35
    int8_t U64_BYTE_CHECKER_VALUE = 0;
    U64_BYTE_CHECKER_VALUE = readbin(U64_BYTE_CHECKER_VALUE, sas, 0);

    int8_t ALIGN_2_VALUE = 0;
    if (U64_BYTE_CHECKER_VALUE == 51) {
      ALIGN_2_VALUE = 4;
    }

    readbin(unk8, sas, 0);

    // o37 1 = little
    int8_t ENDIANNESS = 0;
    ENDIANNESS = readbin(ENDIANNESS, sas, 0);
    Rprintf("ENDIANNESS: %d \n", ENDIANNESS);

    if (ENDIANNESS != 1)
      stop("Big Endian found");

    if (!(PAGE_TYPE == 512 | PAGE_TYPE == 256 | PAGE_TYPE == 0))
      stop("Unhandled page Type %d detected", PAGE_TYPE);


    readbin(unk8, sas, 0);

    // o39 (char?) (1) 49 Unix (2) 50 Win
    uint8_t PLATFORM = 0;
    PLATFORM = readbin(PLATFORM, sas, 0);
    if (debug) Rprintf("PLATFORM: %d \n", PLATFORM);

    // o40
    unkdub = readbin(unk64, sas, 0);
    if (debug) Rcout << unk64 << std::endl;

    // o48
    unkdub = readbin(unk64, sas, 0);
    if (debug) Rcout << unk64 << std::endl;

    // modified block ?
    // repeat 32-39?
    // First block initial file/os type, this block last file/os type?
    // // o56 (char?)
    // uint8_t ALIGN_1_CHECKER_VALUE = 0;
    // ALIGN_1_CHECKER_VALUE = readbin(ALIGN_1_CHECKER_VALUE, sas, 0);
    //
    // int8_t ALIGN_1_VALUE = 0;
    // if (ALIGN_1_CHECKER_VALUE == 33) {
    //   ALIGN_1_VALUE = 4;
    //   Rcpp::stop("ALIGN_1_VALUE == 4!");
    // }
    // Rprintf("ALIGN_1_CHECKER_VALUE: %d \n", ALIGN_1_CHECKER_VALUE);

    readbin(unk8, sas, 0);
    readbin(unk8, sas, 0);
    readbin(unk8, sas, 0);

    // // o59
    // int8_t U64_BYTE_CHECKER_VALUE = 0;
    // U64_BYTE_CHECKER_VALUE = readbin(U64_BYTE_CHECKER_VALUE, sas, 0);
    //
    // int8_t ALIGN_2_VALUE = 0;
    // if (U64_BYTE_CHECKER_VALUE == 33) {
    //   ALIGN_2_VALUE = 4;
    //   Rcpp::stop("ALIGN_2_VALUE == 4!");
    // }
    readbin(unk8, sas, 0);
    readbin(unk8, sas, 0);

    // // o61 1 = little
    // int8_t ENDIANNESS = 0;
    // ENDIANNESS = readbin(ENDIANNESS, sas, 0);
    // Rprintf("ENDIANNESS: %d \n", ENDIANNESS);

    readbin(unk8, sas, 0);
    readbin(unk8, sas, 0);

    // // o62 (char?) (1) 49 Unix (2) 50 Win
    // uint8_t PLATFORM = 0;
    // PLATFORM = readbin(PLATFORM, sas, 0);
    // Rprintf("PLATFORM: %d \n", PLATFORM);
    readbin(unk8, sas, 0);

    // o64
    unk32 = readbin(unk32, sas, 0); // 1st byte: yes (?)
    if (debug) Rcout << unk32 << std::endl;
    unk32 = readbin(unk32, sas, 0); // 4th byte: length of string
    if (debug) Rcout << unk32 << std::endl;
    unk32 = readbin(unk32, sas, 0); // 1. and 2. int16 = 1?
    if (debug) Rcout << unk32 << std::endl;
    unk32 = readbin(unk32, sas, 0); // 0
    if (debug) Rcout << unk32 << std::endl;
    unk32 = readbin(unk32, sas, 0); // 0
    if (debug) Rcout << unk32 << std::endl;

    // o84 SAS FILE
    std::string sasfile (8, '\0');
    readstring(sasfile, sas);
    if (debug) Rcout << sasfile << std::endl;

    // o92 dataset name
    std::string dataset (64, '\0');
    readstring(dataset, sas);
    if (debug) Rcout << dataset << std::endl;

    // o156 filetype 'DATA    '
    std::string filetype (8, '\0');
    readstring(filetype, sas);
    if (debug) Rcout << filetype << std::endl;

    double created = 0;
    created = readbin(created, sas, 0);
    if (debug) Rcout << created << std::endl;

    double modified = 0;
    modified = readbin(modified, sas, 0);
    if (debug) Rcout << modified << std::endl;

    unkdub = readbin(unkdub, sas, 0);
    if (debug) Rcout << unkdub << std::endl;

    unkdub = readbin(unkdub, sas, 0);
    if (debug) Rcout << unkdub << std::endl;

    if (ALIGN_1_VALUE == 4) {
      unk32 = readbin(unk32, sas, 0); // padding for u64
    }

    int32_t headersize = 0;
    headersize = readbin(headersize, sas, 0);
    Rprintf("headersize: %d \n", headersize);

    int32_t pagesize = 0;
    pagesize = readbin(pagesize, sas, 0);
    if (debug) Rprintf("pagesize: %d \n", pagesize);

    int64_t pagecount = 0;
    if (ALIGN_2_VALUE == 4) {
      pagecount = readbin(pagecount, sas, 0);
    } else {
      pagecount = readbin((int32_t)pagecount, sas, 0);
    }
    Rprintf("pagecount: %d \n", pagecount);

    std::vector<int64_t> rowsperpage(pagecount);


    unkdub = readbin(unkdub, sas, 0);
    if (debug) Rcout << unkdub << std::endl;

    std::string sasrel (8, '\0');
    readstring(sasrel, sas);
    Rcout << "SAS release: " << sasrel << std::endl;

    std::string sasserv (16, '\0');
    readstring(sasserv, sas);
    if (debug) Rcout << "SAS server: " << sasserv << std::endl;

    // osversion
    std::string osver (16, '\0');
    readstring(osver, sas);
    if (debug) Rcout << osver << std::endl;

    // osmaker
    std::string osmaker (16, '\0');
    readstring(osmaker, sas);
    if (debug) Rcout << osmaker << std::endl;

    // osname
    std::string osname (16, '\0');
    readstring(osname, sas); // x86_64
    if (debug) Rcout << osname << std::endl;

    // 32 unknown byte
    unkdub = readbin(unkdub, sas, 0);
    if (debug) Rcout << unkdub << std::endl;
    unkdub = readbin(unkdub, sas, 0);
    if (debug) Rcout << unkdub << std::endl;
    unkdub = readbin(unkdub, sas, 0); // 0
    if (debug) Rcout << unkdub << std::endl;
    unkdub = readbin(unkdub, sas, 0); // 0
    if (debug) Rcout << unkdub << std::endl;


    int32_t PAGE_SIZE = 0, PAGE_COUNT = 0;
    PAGE_SIZE = readbin(PAGE_SIZE, sas, 0);
    if (debug) Rprintf("PAGE_SIZE: %d \n", PAGE_SIZE);

    PAGE_COUNT = readbin(PAGE_COUNT, sas, 0);

    // 3rd timestamp
    double thrdts = 0;
    thrdts = readbin(thrdts, sas, 0);
    if (debug) Rcout << thrdts << std::endl;

    // rest is filled with zeros. read it anyway.
    int32_t num_zeros = headersize - sas.tellg();

    Rcout << num_zeros << std::endl;

    for (int i = 0; i < num_zeros; ++i) {
      int8_t zero = 0;
      zero = readbin(zero, sas, 0);
      if (zero!=0)
        warning("Error. Expected 0, is %d", zero);
    }

    // debug
    int64_t pagestart = sas.tellg();
    if (debug) Rprintf("position: %d\n", pagestart);
    // end of Header ---------------------------------------------------------//



    int8_t alignval = 8;
    if (ALIGN_2_VALUE != 4) alignval = 4;

    int64_t rowlength = 0, rowcount = 0;
    int64_t colf_p1 = 0, colf_p2 = 0;
    int64_t colnum = 0;
    std::vector<std::string> stringvec(pagecount) ;

    auto totalrows = 0;
    std::vector<int32_t> totalrowsvec(pagecount);


    // begin reading pages ---------------------------------------------------//
    for (auto pg = 0; pg < pagecount; ++pg) {

      Rcout << "########### PAGE" << pg << " ###############"  << std::endl;
      Rcout << "########### " << sas.tellg() << " ###############"  << std::endl;


      if (pagecount > 0) {
        auto pagenumx = headersize + pg * pagesize;
        sas.seekg(pagenumx, sas.beg);
      }

      // Page Offset Table

      int64_t pageseqnum = 0;

      if (ALIGN_2_VALUE == 4) {
        pageseqnum = readbin(pageseqnum, sas, 0);
        unk64 = readbin(unk64, sas, 0);
        if (debug) Rcout << unk64 << std::endl;
        unk64 = readbin(unk64, sas, 0);
        if (debug) Rcout << unk64 << std::endl;
        unk64 = readbin(unk64, sas, 0);
        if (debug) Rcout << unk64 << std::endl;
      } else {
        pageseqnum = readbin((int32_t)pageseqnum, sas, 0);
        unk32 = readbin(unk32, sas, 0);
        if (debug) Rcout << unk32 << std::endl;
        unk32 = readbin(unk32, sas, 0);
        if (debug) Rcout << unk32 << std::endl;
        unk32 = readbin(unk32, sas, 0);
        if (debug) Rcout << unk32 << std::endl;
      }
      if (debug) Rprintf("pageseqnum: %d \n", pageseqnum);

      PAGE_TYPE = readbin(PAGE_TYPE, sas, 0);
      BLOCK_COUNT = readbin(BLOCK_COUNT, sas, 0);
      SUBHEADER_COUNT = readbin(SUBHEADER_COUNT, sas, 0);
      unk16 = readbin(unk16, sas, 0);

      rowsperpage[pg] = BLOCK_COUNT - SUBHEADER_COUNT;
      totalrows += rowsperpage[pg];
      totalrowsvec[pg] = totalrows;

      // if (debug)
      Rprintf("PAGE_TYPE: %d ; BLOCK_COUNT: %d ; SUBHEADER_COUNT: %d \n",
              PAGE_TYPE, BLOCK_COUNT, SUBHEADER_COUNT);


      if (debug) Rprintf("unk16: %d \n", unk16);

      auto SL = 24; if(ALIGN_2_VALUE != 4) SL = 12;

      std::vector<PO_Tab> potabs(SUBHEADER_COUNT);

      int8_t zero = 0;

      for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
        if (ALIGN_2_VALUE == 4) {

          potabs[i].SH_OFF = readbin(potabs[i].SH_OFF, sas, 0);
          potabs[i].SH_LEN = readbin(potabs[i].SH_LEN, sas, 0);
          potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, 0);
          potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, 0);

          zero = readbin(zero, sas, 0);
          zero = readbin(zero, sas, 0);
          zero = readbin(zero, sas, 0);
          zero = readbin(zero, sas, 0);
          zero = readbin(zero, sas, 0);
          zero = readbin(zero, sas, 0);

          // if (debug)
          Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPRESSION: %d ; SH_TYPE: %d \n",
                  potabs[i].SH_OFF, potabs[i].SH_LEN,
                  potabs[i].COMPRESSION, potabs[i].SH_TYPE);

        } else {

          potabs[i].SH_OFF = readbin((int32_t)potabs[i].SH_OFF, sas, 0);
          potabs[i].SH_LEN = readbin((int32_t)potabs[i].SH_LEN, sas, 0);
          potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, 0);
          potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, 0);

          zero = readbin(zero, sas, 0);
          zero = readbin(zero, sas, 0);

          // if (debug)
          Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPRESSION: %d ; SH_TYPE: %d \n",
                  potabs[i].SH_OFF, potabs[i].SH_LEN,
                  potabs[i].COMPRESSION, potabs[i].SH_TYPE);
        }
      }

      // debug

      if (pg == 0) {
        while (sas.tellg() % 8 != 0) {
          readbin(unk32, sas, 0); // padding?
        }
      }

      auto sh_end_pos = 0;

      if (PAGE_TYPE != 0)
        sh_end_pos = sas.tellg();

      Rprintf("sh_end_pos: %d\n", sh_end_pos);
      data_pos.push_back( sh_end_pos );

      for (auto sc = 0; sc < SUBHEADER_COUNT; ++sc)
      {
        // if ((PAGE_TYPE == 256) | (PAGE_TYPE == 512)) {
        //   auto pos = pagestart + potabs[sc].SH_OFF;
        //   sas.seekg(pos, sas.beg);
        // } else if (PAGE_TYPE == 0 || PAGE_TYPE == 1024) {

        //   auto pos = headersize + (pg+1) * potabs[sc].SH_OFF;
        //
        //   if (pg != 0)
        //     pos = headersize + (pg+1) * potabs[sc].SH_OFF;
        //
        //   sas.seekg(pos, sas.beg);
        //
        // } else if (PAGE_TYPE == 1024) {
        auto pos = (headersize + pg * pagesize) + potabs[sc].SH_OFF;
        sas.seekg(pos, sas.beg);
        // }
        // if (debug) Rprintf("%d \n", pos);


        int64_t sas_offset = alignval;
        if (ALIGN_1_VALUE == 4) {
          sas_offset = readbin(sas_offset, sas, 0);
        } else {
          sas_offset = readbin((int32_t)sas_offset, sas, 0);
        }

        std::string sas_hex = int32_to_hex(sas_offset);

        auto sas_offset_table = 0;
        if (sas_hex.compare("f7f7f7f7") == 0 ||
            sas_hex.compare("fffffffff7f7f7f7") == 0)
          sas_offset_table = 1;
        if (sas_hex.compare("fffffc00") == 0)
          sas_offset_table = 2;
        if (sas_hex.compare("fffffbfe") == 0)
          sas_offset_table = 3;
        if (sas_hex.compare("f6f6f6f6") == 0)
          sas_offset_table = 4;
        if (sas_hex.compare("fffffffd") == 0)
          sas_offset_table = 5;
        if (sas_hex.compare("ffffffff") == 0)
          sas_offset_table = 6;
        if (sas_hex.compare("fffffffc") == 0)
          sas_offset_table = 7;

        const int8_t SOT = sas_offset_table;


        switch(sas_offset_table)
        {


          // new offset ----------------------------------------------------- //
        case 1:
        { /* Row Size */


          Rcout << "-------- case 1 "<< sas.tellg() << " --------" << std::endl;

          int16_t pgwpossh = 0, pgwpossh2 = 0, numzeros = 37,
            sh_num = 0, cn_maxlen = 0, l_maxlen = 0,
            rowsonpg = 0;
          int32_t pgidx = 0;
          int64_t pgsize = 0, pgc = 0, rcmix = 0, pgwsh = 0, pgwsh2 = 0;


          if (ALIGN_2_VALUE == 4) {
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;

            rowlength = readbin(rowlength, sas, 0);
            Rcout << rowlength << std::endl;
            rowcount = readbin(rowcount, sas, 0);
            if (debug) Rcout << rowcount << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;

            colf_p1 = readbin(colf_p1, sas, 0);
            if (debug) Rcout << colf_p1 << std::endl;
            colf_p2 = readbin(colf_p2, sas, 0);
            if (debug) Rcout << colf_p2 << std::endl;
            unk64 = readbin(unk64, sas, 0); // p3 and p4?
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            pgsize = readbin(pgsize, sas, 0);
            unk64 = readbin(unk64, sas, 0);
            rcmix =  readbin(rcmix, sas, 0);

            unk64 = readbin(unk64, sas, 0); /* end of initial header ? */
            unk64 = readbin(unk64, sas, 0); /*                         */

            for (int z = 0; z < numzeros; ++z) {
              unk64 = readbin(unk64, sas, 0);
              if (unk64 != 0)
                warning("val is %d. expected a zero", unk64);
            }

            pgidx = readbin(pgidx, sas, 0);

            // padding? 68 bytes: zeros
            for (int z = 0; z < 8; ++z) {
              unk64 = readbin(unk64, sas, 0);
            }
            unk32 = readbin(unk32, sas, 0);

            unk64 = readbin(unk64, sas, 0); // val 1?
            unk16 = readbin(unk16, sas, 0); // val 2?

            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding

            pgwsh = readbin(pgwsh, sas, 0);
            pgwpossh = readbin(pgwpossh, sas, 0);

            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding

            pgwsh2 = readbin(pgwsh2, sas, 0);
            pgwpossh2 = readbin(pgwpossh2, sas, 0);

            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding

            pgc = readbin(pgc, sas, 0);

            unk16 = readbin(unk16, sas, 0); // val ?
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding

            unk64 = readbin(unk64, sas, 0); // val 1?

            addtextoff = readbin(addtextoff, sas, 0); // val 7 | 8?
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding
            unk16 = readbin(unk16, sas, 0); // padding

            for (int z = 0; z < 10; ++z) {
              unk64 = readbin(unk64, sas, 0); // 0
            }

            Rcout << "###############" << std::endl;

            unk16 = readbin(unk16, sas, 0); // val
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 0|8 ?
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 4
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 0
            Rcout << unk16 << std::endl;
            todata = readbin(unk16, sas, 0); // val 12,32|0?
            Rcout << todata << std::endl;

            if (todata == 12)
              hasproc = false;

            Rcout << "###############" << std::endl;

            swlen = readbin(swlen, sas, 0);
            Rcout << swlen << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 0
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 20 | 28
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 8
            Rcout << unk16 << std::endl;
            // unk64 = readbin(unk64, sas, 0); // ?
            // Rcout << unk64 << std::endl;

            Rcout << "###############" << std::endl;

            unk16 = readbin(unk16, sas, 0); // 0
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 12
            Rcout << unk16 << std::endl;
            comprlen = readbin(unk16, sas, 0); // 8
            Rcout << comprlen << std::endl;
            unk16 = readbin(unk16, sas, 0); // 0
            Rcout << unk16 << std::endl;

            Rcout << "###############" << std::endl;

            unk16 = readbin(unk16, sas, 0); // 12 | 20
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 8
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 0
            Rcout << unk16 << std::endl;
            textoff = readbin(textoff, sas, 0); // 28 | 36
            Rcout << textoff << std::endl;
            proclen = readbin(proclen, sas, 0);
            Rcout << proclen << std::endl;

            Rcout << "###############" << std::endl;

            for (int z = 0; z < 8; ++z) {
              unk32 = readbin(unk32, sas, 0); // 0
            }

            unk16 = readbin(unk16, sas, 0); // 4
            unk16 = readbin(unk16, sas, 0); // 1

            sh_num = readbin(sh_num, sas, 0);
            cn_maxlen = readbin(cn_maxlen, sas, 0);
            l_maxlen = readbin(l_maxlen, sas, 0);

            for (int z = 0; z < 3; ++z) {
              unk32 = readbin(unk32, sas, 0); // 0
            }

            rowsonpg = readbin(rowsonpg, sas, 0);

            for (int z = 0; z < 10; ++z) {
              unk32 = readbin(unk32, sas, 0); // 0
            }

          } else {
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;

            rowlength = readbin((int32_t)rowlength, sas, 0);
            Rcout << "rowlength "<< rowlength << std::endl;
            rowcount = readbin((int32_t)rowcount, sas, 0);
            Rcout << "rowcount "<< rowcount << std::endl;
            if (debug) Rcout << rowcount << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;

            colf_p1 = readbin((int32_t)colf_p1, sas, 0);
            Rcout << colf_p1 << std::endl;
            colf_p2 = readbin((int32_t)colf_p2, sas, 0);
            Rcout << colf_p2 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            pgsize = readbin((int32_t)pgsize, sas, 0);
            unk32 = readbin(unk32, sas, 0);
            rcmix =  readbin((int32_t)rcmix, sas, 0);
            unk32 = readbin(unk32, sas, 0);
            unk32 = readbin(unk32, sas, 0);

            for (int z = 0; z < numzeros; ++z) {
              unk64 = readbin((int32_t)unk64, sas, 0);
              if (unk64 != 0)
                warning("val is %d. expected a zero", unk64);
            }

            pgidx = readbin(pgidx, sas, 0);


            for (int z = 0; z < 8; ++z) {
              unk64 = readbin((int32_t)unk64, sas, 0);
            }

            // padding?
            unk32 = readbin(unk32, sas, 0);
            unk32 = readbin(unk32, sas, 0);

            unk32 = readbin(unk32, sas, 0); // val 1?
            unk16 = readbin(unk16, sas, 0); // val 2?

            unk16 = readbin(unk16, sas, 0); // padding

            pgwsh = readbin((int32_t)pgwsh, sas, 0);
            pgwpossh = readbin(pgwpossh, sas, 0);

            unk16 = readbin(unk16, sas, 0); // padding

            pgwsh2 = readbin((int32_t)pgwsh2, sas, 0);
            pgwpossh2 = readbin(pgwpossh2, sas, 0);

            unk16 = readbin(unk16, sas, 0); // padding

            pgc = readbin((int32_t)pgc, sas, 0);

            unk16 = readbin(unk16, sas, 0); // val ?
            unk16 = readbin(unk16, sas, 0); // padding

            unk64 = readbin((int32_t)unk64, sas, 0); // val 1?
            Rcout << unk64 << std::endl;

            addtextoff = readbin(addtextoff, sas, 0); // val 7 | 8?
            Rcout << addtextoff << std::endl;
            unk16 = readbin(unk16, sas, 0); // padding

            for (int z = 0; z < 10; ++z) {
              unk64 = readbin((int32_t)unk64, sas, 0); // 0
            }

            Rcout << "###############" << std::endl;

            unk16 = readbin(unk16, sas, 0); // val
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 0|8 ?
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 4
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 0
            Rcout << unk16 << std::endl;
            todata = readbin(todata, sas, 0); // val 12,32|0? //
            Rcout << todata << std::endl;

            if (todata == 12)
              hasproc = false;

            Rcout << "###############" << std::endl;

            swlen = readbin(swlen, sas, 0);
            Rcout << "swlen " << swlen << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 0?
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // val 20?
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); //
            Rcout << unk16 << std::endl;

            Rcout << "###############" << std::endl;

            unk16 = readbin(unk16, sas, 0); // 0
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 12
            Rcout << unk16 << std::endl;
            comprlen = readbin(unk16, sas, 0); // 8 compression code length?
            Rcout << comprlen << std::endl;
            unk16 = readbin(unk16, sas, 0); // 0
            Rcout << unk16 << std::endl;

            Rcout << "###############" << std::endl;

            unk16 = readbin(unk16, sas, 0); // 12
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 8
            Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 0
            Rcout << unk16 << std::endl;
            textoff = readbin(textoff, sas, 0); // 28
            Rcout << textoff << std::endl;
            proclen = readbin(proclen, sas, 0);
            Rcout << proclen << std::endl;

            Rcout << "###############" << std::endl;

            for (int z = 0; z < 8; ++z) {
              unk32 = readbin(unk32, sas, 0); // 0
            }

            unk16 = readbin(unk16, sas, 0); // 4
            unk16 = readbin(unk16, sas, 0); // 1

            sh_num = readbin(sh_num, sas, 0);
            cn_maxlen = readbin(cn_maxlen, sas, 0);
            l_maxlen = readbin(l_maxlen, sas, 0);

            for (int z = 0; z < 3; ++z) {
              unk32 = readbin(unk32, sas, 0); // 0
            }

            rowsonpg = readbin(rowsonpg, sas, 0);

            for (int z = 0; z < 10; ++z) {
              unk32 = readbin(unk32, sas, 0); // 0
            }

          }


          if (!hasproc)
            proclen = 0;

          // if (addtextoff == 8)
          //   textoff = 0;

          if (todata > 0) {
            offset = todata + swlen;
          } else {
            offset = textoff + proclen + swlen;
          }

          // if (debug)
          Rprintf("swlen = %d, proclen = %d, todata %d, textoff %d, offset %d \n",
                  swlen, proclen, todata, textoff, offset);

          break;
        }


          // new offset ----------------------------------------------------- //
        case 2:
        {

          Rcout << "-------- case 2 "<< sas.tellg() << " --------" << std::endl;

          int64_t off = 0;

          if (ALIGN_2_VALUE == 4) {
            off = readbin(off, sas, 0);
            // Rcout << off << std::endl;
            unk64 = readbin(unk64, sas, 0);
            // Rcout << unk64 << std::endl;
          } else {
            off = readbin((int32_t)off, sas, 0);
            // Rcout << off << std::endl;
            unk32 = readbin(unk32, sas, 0);
            // Rcout << unk32 << std::endl;
          }

          int16_t num_nonzero = 0;
          num_nonzero = readbin(num_nonzero, sas, 0);
          // Rcout << num_nonzero << std::endl;

          int8_t unklen = 94; if (ALIGN_2_VALUE != 4) unklen = 50;
          std::string unkstr(unklen, '\0');
          unkstr = readstring(unkstr, sas);
          // Rcout << unkstr << std::endl;

          std::vector<SCV> scv(12);

          for (int8_t i = 0; i < 12; ++i) {

            if (ALIGN_2_VALUE == 4) {
              scv[i].SIG = readbin(scv[i].SIG, sas, 0);
              scv[i].FIRST = readbin(scv[i].FIRST, sas, 0);
              scv[i].F_POS = readbin(scv[i].F_POS, sas, 0);

              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;

              scv[i].LAST = readbin(scv[i].LAST, sas, 0);
              scv[i].L_POS = readbin(scv[i].L_POS, sas, 0);

              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;

            } else {
              scv[i].SIG = readbin((int32_t)scv[i].SIG, sas, 0);
              scv[i].FIRST = readbin((int32_t)scv[i].FIRST, sas, 0);
              scv[i].F_POS = readbin(scv[i].F_POS, sas, 0);

              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;

              scv[i].LAST = readbin((int32_t)scv[i].LAST, sas, 0);
              scv[i].L_POS = readbin(scv[i].L_POS, sas, 0);

              unk16 = readbin(unk16, sas, 0);
              // Rcout << unk16 << std::endl;
            }

            if (debug)
              Rprintf("SIG %d; FIRST %d; F_POS %d; LAST %d; L_POS %d\n",
                      scv[i].SIG, scv[i].FIRST, scv[i].F_POS,
                      scv[i].LAST, scv[i].L_POS);

          }

          break;
        }


          // new offset ----------------------------------------------------- //
        case 3:
        {

          Rcout << "-------- case 3 "<< sas.tellg() << " --------" << std::endl;

          hasattributes = 1;

          int64_t unk64 = 0;

          int8_t unklen = 38;if (ALIGN_2_VALUE != 4) unklen = 30;
          std::string unkstr(unklen, '\0');
          unkstr = readstring(unkstr, sas);
          // Rcout << unkstr << std::endl;

          idxofflen fmts;
          idxofflen lbls;
          idxofflen unks;

          fmts.IDX = readbin(fmts.IDX, sas, 0);
          fmts.OFF = readbin(fmts.OFF, sas, 0);
          fmts.LEN = readbin(fmts.LEN, sas, 0);

          fmt.push_back(fmts);

          lbls.IDX = readbin(lbls.IDX, sas, 0);
          lbls.OFF = readbin(lbls.OFF, sas, 0);
          lbls.LEN = readbin(lbls.LEN, sas, 0);

          lbl.push_back(lbls);

          unks.IDX = readbin(unks.IDX, sas, 0);
          unks.OFF = readbin(unks.OFF, sas, 0);
          unks.LEN = readbin(unks.LEN, sas, 0);

          unk.push_back(unks);

          if (unks.IDX != 0 | unks.OFF != 0 | unks.LEN != 0)
            warning("unk not 0 as expected\n");

          break;
        }


          // new offset ----------------------------------------------------- //
        case 4:
        { /* Column Size */

          Rcout << "-------- case 4 "<< sas.tellg() << " --------" << std::endl;

          int64_t unk1 = 0;

          if (ALIGN_2_VALUE == 4) {
            colnum = readbin(colnum, sas, 0);
            unk1 = readbin(unk1, sas, 0);
          } else {
            colnum = readbin((int32_t)colnum, sas, 0);
            unk1 = readbin((int32_t)unk1, sas, 0);
          }

          if (debug)
            Rprintf("colnum %d; unk1 %d\n",
                    colnum, unk1);

          break;
        }

          // new offset ----------------------------------------------------- //

        case 5:
        { /* Column Text */

          Rcout << "-------- case 5 "<< sas.tellg() << " --------" << std::endl;
          Rcout << sas.tellg() << std::endl;

          int16_t len = 0;
          int8_t tmp = 16; // always 16?
          if (!hasproc) tmp = 0;

          len = readbin(len, sas, 0);
          Rcout << len << std::endl;
          unk16 = readbin(unk16, sas, 0); // 0 vars on p1 ?
          Rcout << unk16 << std::endl;
          unk16 = readbin(unk16, sas, 0); // 0 vars on p2 ?
          Rcout << unk16 << std::endl;
          unk16 = readbin(unk16, sas, 0); // 0
          // Rcout << unk16 << std::endl;

          if (PAGE_TYPE != 1024) {
            unk16 = readbin(unk16, sas, 0); // 0
            // Rcout << unk16 << std::endl;
            unk16 = readbin(unk16, sas, 0); // 0
            // Rcout << unk16 << std::endl;
          }

          std::string compression = "";
          std::string proc = "";
          std::string sw = "";

          Rprintf("%d %d %d %d \n", comprlen, tmp, proclen, swlen);
          varoffset = comprlen + tmp + proclen + swlen;

          if (PAGE_TYPE != 1024) {

            // compression
            if (comprlen > 0) {
              compression.resize(comprlen, '\0');
              compression = readstring(compression, sas);

              if (compression.compare("SASYZCRL") == 0)
                compr = 1;

              if (compression.compare("SASYZCR2") == 0)
                compr = 2;
            }

            // 16 whitespaces
            std::string empty (tmp, '\0');
            if (tmp > 0) {
              empty = readstring(empty, sas);
              if (!(empty.compare("                ") == 0))
                warning("non empty 'empty' string found");
            }

            // proc that created the file
            if (proclen > 0) {
              proc.resize(proclen, '\0');
              proc = readstring(proc, sas);
            }

            // additional software string
            if (swlen > 0) {
              sw.resize(swlen, '\0');
              sw = readstring(sw, sas);
            }


            Rcout << compression << "\n" <<
              empty << "\n" <<
                proc << "\n" << sw << std::endl;

          }

          int32_t newlen = len - varoffset;

          std::string CN_IDX_STR (newlen, '\0');
          CN_IDX_STR = readstring(CN_IDX_STR, sas);
          Rcout << CN_IDX_STR << std::endl;

          if (debug)
            Rprintf("SH_LEN %d; len %d; newlen: %d\n",
                    potabs[sc].SH_LEN, len, newlen);




          stringvec[pg] = CN_IDX_STR;

          break;
        }


          // new offset ----------------------------------------------------- //
        case 6:
        { /* Column Name */

          Rcout << "-------- case 6 "<< sas.tellg() << " --------" << std::endl;

          int16_t lenremain = 0;
          lenremain = readbin(lenremain, sas, 0);
          if (debug) Rprintf("lenremain %d \n", lenremain);

          unk16 = readbin(unk16, sas, 0);
          // Rcout << unk16 << std::endl;
          unk16 = readbin(unk16, sas, 0);
          // Rcout << unk16 << std::endl;
          unk16 = readbin(unk16, sas, 0);
          // Rcout << unk16 << std::endl;

          auto cmax = (lenremain + alignval)/8;

          // Rcout << offset << std::endl;

          /* Column Name Pointers */
          std::vector<CN_Poi> cnpois(cmax);

          for (auto i = 0; i < cmax; ++i) {
            auto idx    = readbin(cnpois[i].CN_IDX, sas, 0);
            auto off    = readbin(cnpois[i].CN_OFF, sas, 0);
            auto len    = readbin(cnpois[i].CN_LEN, sas, 0);
            auto zeros  = readbin(cnpois[i].zeros,  sas, 0);


            if (debug)
              Rprintf("CN_IDX %d; CN_OFF %d; CN_LEN %d; zeros %d \n",
                      idx, off, len, zeros);

            if (len > 0) {
              off -= offset;

              std::string varname = stringvec[idx].substr(off, len);
              // Rcout << varname << std::endl;
              varnames.push_back(varname);
            }

          }

          break;
        }


          // new offset ----------------------------------------------------- //
        case 7:
        { /* Column Attributes */

          Rcout << "-------- case 7 "<< sas.tellg() << " --------" << std::endl;

          int16_t lenremain = 0;
          lenremain = readbin(lenremain, sas, 0);
          if (debug) Rprintf("lenremain %d \n", lenremain);

          // zeros as padding?
          unk16 = readbin(unk16, sas, 0);
          // Rcout << unk16 << std::endl;
          unk16 = readbin(unk16, sas, 0);
          // Rcout << unk16 << std::endl;
          unk16 = readbin(unk16, sas, 0);
          // Rcout << unk16 << std::endl;

          auto cmax = (lenremain + alignval) / (alignval+8);
          if (debug) Rcout << cmax << std::endl;

          /* Column Attributes Pointers */
          std::vector<CN_Att> capois(cmax);

          for (auto i = 0; i < cmax; ++i) {

            if (ALIGN_1_VALUE == 4) {
              capois[i].CN_OFF     = readbin(capois[i].CN_OFF, sas, 0);
            } else {
              capois[i].CN_OFF     = readbin((int32_t)capois[i].CN_OFF, sas, 0);
            }
            capois[i].CN_WID     = readbin(capois[i].CN_WID, sas, 0);
            capois[i].NM_FLAG    = readbin(capois[i].NM_FLAG, sas, 0);
            capois[i].CN_TYP     = readbin(capois[i].CN_TYP, sas, 0);
            capois[i].UNK8       = readbin(capois[i].UNK8, sas, 0);


            // if (debug)
            Rprintf("OFF %d; WID: %d; FLAG %d; TYP %d; UNK8 %d\n",
                    capois[i].CN_OFF, capois[i].CN_WID, capois[i].NM_FLAG,
                    capois[i].CN_TYP, capois[i].UNK8 );

            if (capois[i].CN_TYP > 0) {
              coloffset.push_back( capois[i].CN_OFF );
              colwidth.push_back( capois[i].CN_WID );
              vartyps.push_back( capois[i].CN_TYP );
            }

          }

          break;
        }

          // not implemented ------------------------------------------------ //
        default:
        {

          Rcout << "---- unimplemented "<< sas.tellg() << " ----" << std::endl;

          // else it is padding?
          if (potabs[sc].SH_LEN > alignval)
          {
            sas_hex = int32_to_hex(sas_offset);
            // Rcout << sas_hex << std::endl;
            // Rcout << "at offset " << sas.tellg() << std::endl;

            auto unklen = potabs[sc].SH_LEN - alignval;

            std::string unkstr (unklen, '\0');

            unkstr = readstring(unkstr, sas);
            // std::cout << unkstr << std::endl;

          }

          break;
        }

        }
      }
    }

    if (hasattributes) {
      Rcout << "------------ attributes -------------" << std::endl;

      auto idx = 0;

      if (addtextoff == 8 & PAGE_TYPE == 1024) {
        idx = (stringvec.size() -1);
        offset = 8;
      }

      // read format and labels
      for (auto var = 0; var < colnum; ++var) {

        // std::cout << stringvec[pg] << std::endl;

        std::string format = "";
        int64_t fidx = 0, foff = fmt[var].OFF, flen = fmt[var].LEN;
        // Rprintf("fidx %d; foff %d; flen %d \n", fidx, foff, flen);
        if (flen > 0) {
          foff -= offset;
          format = stringvec[idx].substr(foff, flen);
        }

        std:: string label = "";
        int64_t lidx = 0, loff = lbl[var].OFF, llen = lbl[var].LEN;
        // Rprintf("lidx %d; loff %d; llen %d \n", lidx, loff, llen);
        if (llen > 0) {
          loff -= offset;
          label = stringvec[idx].substr(loff, llen);
        }

        // Rcout << format << " : " << label << std::endl;

        formats.push_back( format );
        labels.push_back( label );
      }

    }

    if (compr != 0)
      stop("File contains unhandled compression");


    // ---------------------------------------------------------------------- //

    // 1. Create Rcpp::List
    Rcpp::List df(colnum);
    for (auto i=0; i<colnum; ++i)
    {
      int32_t const type = vartyps[i];

      if (debug) Rcout << type << std::endl;

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


    // // new offset ---------------------------------------------------------//


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

      // Rcout << totalrowsvec[page] << std::endl;
      // Rprintf("ii: %d\n", ii);

      auto pp = data_pos[page];
      Rcout << pp << std::endl;

      sas.seekg((pp + rowlength * ii), sas.beg);
      Rcout << "---- data "<< sas.tellg() << " ----" << std::endl;
      // Rcout << "page: " << page << std::endl;
      int64_t tempoffs = sas.tellg();
      // Rcout << tempoffs << std::endl;

      auto pos = 0;

      for (auto j = 0; j < colnum; ++j) {

        auto wid = colwidth[j];
        auto typ = vartyps[j];

        int64_t off = coloffset[j];
        int64_t readpos = tempoffs + off;
        // Rcout << wid << std::endl;
        // Rcout << off << std::endl;
        // Rcout << readpos << std::endl;

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

          val_d = readbin(val_d, sas, 0);

          // Rcout << val_d << std::endl;

          if (std::isnan(val_d))
            REAL(VECTOR_ELT(df,j))[i] = NA_REAL;
          else
            REAL(VECTOR_ELT(df,j))[i] = val_d;

        }

        if (wid > 0 & typ == 2) {

          std::string val_s(wid, '\0');

          sas.seekg(readpos, sas.beg);

          val_s = readstring(val_s, sas);

          val_s = std::regex_replace(val_s,
                                     std::regex(" +$"), "$1");

          as<CharacterVector>(df[j])[i] = val_s;

        }

        // std::string val_str = row.substr(pos, wid);

        // Rcout << wid << std::endl;

        // as<CharacterVector>(df[j])[i] = val_str;
        //
        // Rcout << val_str << std::endl;
        //
        // pos += wid;
      }

      ++ii;
    }




    IntegerVector rvec = seq(1, rowcount);

    // 3. Create a data.frame
    df.attr("row.names") = rvec;
    df.attr("names") = varnames;
    df.attr("class") = "data.frame";

    // close file
    sas.close();

    df.attr("labels") = labels;
    df.attr("formats") = formats;

    return(df);

  } else {
    return (-1);
  }
}
