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


// RLE_COMPRESSION = b'SASYZCRL'
// RDC_COMPRESSION = b'SASYZCR2'


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

    int8_t  unk8  = 0;
    int16_t unk16 = 0;
    int32_t unk32 = 0;
    int64_t unk64 = 0;
    double unkdub = 0;

    Rcpp::IntegerVector vartyps;
    Rcpp::IntegerVector colwidth;
    Rcpp::IntegerVector data_pos;
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
    if (debug) Rprintf("ENDIANNESS: %d \n", ENDIANNESS);

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
    if (debug) Rcout << "SAS release: " << sasrel << std::endl;

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

    for (int i = 0; i < num_zeros; ++i) {
      int8_t zero = 0;
      zero = readbin(zero, sas, 0);
      if (zero!=0)
        Rcpp::stop("Error. Expected 0, is %d", zero);
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
    std::vector<std::string> stringvec ;


    // begin reading pages ---------------------------------------------------//
    for (auto pg = 0; pg < pagecount; ++pg) {

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

      int16_t PAGE_TYPE = 0, BLOCK_COUNT = 0, SUBHEADER_COUNT = 0;
      PAGE_TYPE = readbin(PAGE_TYPE, sas, 0);
      BLOCK_COUNT = readbin(BLOCK_COUNT, sas, 0);
      SUBHEADER_COUNT = readbin(SUBHEADER_COUNT, sas, 0);
      unk16 = readbin(unk16, sas, 0);

      rowsperpage[pg] = BLOCK_COUNT - SUBHEADER_COUNT;

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

      //   if (PAGE_TYPE == 256 || PAGE_TYPE == 512) {
      //     sh_end_pos = sas.tellg();
      //     // debug
      //     Rprintf("position: %d\n", sh_end_pos);
      //
      //     data_pos.push_back( sh_end_pos );
      //
      //   }

      // debug
      auto sh_end_pos = sas.tellg();
      Rprintf("sh_end_pos: %d\n", sh_end_pos);
      data_pos.push_back( sh_end_pos );

      for (auto sc = 0; sc < SUBHEADER_COUNT; ++sc)
      {

        auto pos = pagestart + potabs[sc].SH_OFF;
        sas.seekg(pos, sas.beg);
        if (debug) Rprintf("%d \n", pos);


        int64_t sas_offset = alignval;
        if (ALIGN_1_VALUE == 4) {
          sas_offset = readbin(sas_offset, sas, 0);
        } else {
          sas_offset = readbin((int32_t)sas_offset, sas, 0);
        }

        // Rcout << std::hex << sas_offset << std::endl;

        std::string sas_hex = int32_to_hex(sas_offset);

        // new offset ----------------------------------------------------------//

        if (sas_hex.compare("f7f7f7f7") == 0 ||
            sas_hex.compare("fffffffff7f7f7f7") == 0) { // F7F7F7F7

          if (ALIGN_2_VALUE == 4) {
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
          } else {
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
          }


          if (ALIGN_2_VALUE == 4) {
            rowlength = readbin(rowlength, sas, 0);
            Rcout << rowlength << std::endl;
            rowcount = readbin(rowcount, sas, 0);
            if (debug) Rcout << rowcount << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
          } else {
            rowlength = readbin((int32_t)rowlength, sas, 0);
            Rcout << rowlength << std::endl;
            rowcount = readbin((int32_t)rowcount, sas, 0);
            if (debug) Rcout << rowcount << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
          }


          if (ALIGN_2_VALUE == 4) {
            colf_p1 = readbin(colf_p1, sas, 0);
            if (debug) Rcout << colf_p1 << std::endl;
            colf_p2 = readbin(colf_p2, sas, 0);
            if (debug) Rcout << colf_p2 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
            unk64 = readbin(unk64, sas, 0);
            if (debug) Rcout << unk64 << std::endl;
          } else {
            colf_p1 = readbin((int32_t)colf_p1, sas, 0);
            if (debug) Rcout << colf_p1 << std::endl;
            colf_p2 = readbin((int32_t)colf_p2, sas, 0);
            if (debug) Rcout << colf_p2 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
            unk32 = readbin(unk32, sas, 0);
            if (debug) Rcout << unk32 << std::endl;
          }

        } else

          // new offset ----------------------------------------------------------//


          if (sas_hex.compare("fffffc00") == 0) { // f6f6f6f6

            int64_t off = 0;

            if (ALIGN_2_VALUE == 4) {
              off = readbin(off, sas, 0);
              Rcout << off << std::endl;
              unk64 = readbin(unk64, sas, 0);
              Rcout << unk64 << std::endl;
            } else {
              off = readbin((int32_t)off, sas, 0);
              Rcout << off << std::endl;
              unk32 = readbin(unk32, sas, 0);
              Rcout << unk32 << std::endl;
            }

            int16_t num_nonzero = 0;
            num_nonzero = readbin(num_nonzero, sas, 0);
            Rcout << num_nonzero << std::endl;

            int8_t unklen = 94; if (ALIGN_2_VALUE != 4) unklen = 50;
            std::string unkstr(unklen, '\0');
            unkstr = readstring(unkstr, sas);
            Rcout << unkstr << std::endl;

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

              Rprintf("SIG %d; FIRST %d; F_POS %d; LAST %d; L_POS %d\n",
                      scv[i].SIG, scv[i].FIRST, scv[i].F_POS,
                      scv[i].LAST, scv[i].L_POS);

            }


          } else

            // new offset ----------------------------------------------------------//


            if (sas_hex.compare("fffffbfe") == 0) { // f6f6f6f6

              int64_t sig = 0;

              if (ALIGN_2_VALUE == 4) {
                sig = readbin(sig, sas, 0);
              } else {
                sig = readbin((int32_t)sig, sas, 0);
              }

              int8_t unklen = 38; if (ALIGN_2_VALUE != 4) unklen = 30;
              std::string unkstr(unklen, '\0');
              unkstr = readstring(unkstr, sas);

              int16_t colfsig = 0, colfoff = 0, colflen = 0, collidx = 0, colllen = 0;

              colfsig = readbin(colfsig, sas, 0);
              colfoff = readbin(colfoff, sas, 0);
              colflen = readbin(colflen, sas, 0);
              collidx = readbin(collidx, sas, 0);
              colllen = readbin(colllen, sas, 0);
              unk16 = readbin(unk16, sas, 0);
              unk16 = readbin(unk16, sas, 0);
              unk16 = readbin(unk16, sas, 0);

            } else


              // new offset ----------------------------------------------------------//


              if (sas_hex.compare("f6f6f6f6") == 0) { // f6f6f6f6


                if (ALIGN_2_VALUE == 4) {
                  colnum = readbin(colnum, sas, 0);
                  if (debug) Rcout << colnum << std::endl;
                  unk64 = readbin(unk64, sas, 0);
                  if (debug) Rcout << unk64 << std::endl;
                  unk64 = readbin(unk64, sas, 0);
                  if (debug) Rcout << unk64 << std::endl;
                } else {
                  colnum = readbin((int32_t)colnum, sas, 0);
                  if (debug) Rcout << colnum << std::endl;
                  unk32 = readbin(unk32, sas, 0);
                  if (debug) Rcout << unk32 << std::endl;
                  unk32 = readbin(unk32, sas, 0);
                  if (debug) Rcout << unk32 << std::endl;
                }

              } else

                // new offset ----------------------------------------------------------//

                if (sas_hex.compare("fffffffd") == 0) { // f6f6f6f6

                  int16_t len = 0;
                  len = readbin(len, sas, 0);
                  if (debug) Rprintf("%d\n", len);

                  unk16 = readbin(unk16, sas, 0);
                  if (debug) Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, 0);
                  if (debug) Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, 0);
                  if (debug) Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, 0);
                  if (debug) Rcout << unk16 << std::endl;
                  unk16 = readbin(unk16, sas, 0);
                  if (debug) Rcout << unk16 << std::endl;

                  std::string CN_IDX_STR (len, '\0');
                  CN_IDX_STR = readstring(CN_IDX_STR, sas);

                  // std::cout << CN_IDX_STR << std::endl;

                  stringvec.push_back(CN_IDX_STR);

                } else


                  // new offset ----------------------------------------------------------//


                  if (sas_hex.compare("ffffffff") == 0) { // f6f6f6f6

                    int16_t lenremain = 0;
                    lenremain = readbin(lenremain, sas, 0);
                    if (debug) Rprintf("lenremain %d \n", lenremain);

                    unk16 = readbin(unk16, sas, 0);
                    if (debug) Rcout << unk16 << std::endl;
                    unk16 = readbin(unk16, sas, 0);
                    if (debug) Rcout << unk16 << std::endl;
                    unk16 = readbin(unk16, sas, 0);
                    if (debug) Rcout << unk16 << std::endl;

                    auto cmax = (lenremain + alignval)/8;

                    std::vector<CN_Poi> cnpois(cmax);

                    for (auto i = 0; i < cmax; ++i) {
                      auto idx    = readbin(cnpois[i].CN_IDX, sas, 0);
                      auto off    = readbin(cnpois[i].CN_OFF, sas, 0);
                      auto len    = readbin(cnpois[i].CN_LEN, sas, 0);
                      auto zeros  = readbin(cnpois[i].zeros,  sas, 0);


                      if (!(len <= 0)) {
                        off -= 12; // reduce off

                        if (debug)
                          Rprintf("CN_IDX %d; CN_OFF %d; CN_LEN %d; zeros %d \n",
                                  idx, off,
                                  len, zeros);

                        std::string varname = stringvec[idx].substr(off, len);

                        if (debug) Rcout << varname << std::endl;

                        varnames.push_back(varname);
                      }

                    }

                  } else


                    // new offset ----------------------------------------------------------//

                    if (sas_hex.compare("fffffffc") == 0) { // f6f6f6f6

                      int16_t lenremain = 0;
                      lenremain = readbin(lenremain, sas, 0);
                      if (debug) Rprintf("lenremain %d \n", lenremain);

                      unk16 = readbin(unk16, sas, 0);
                      // Rcout << unk16 << std::endl;
                      unk16 = readbin(unk16, sas, 0);
                      // Rcout << unk16 << std::endl;
                      unk16 = readbin(unk16, sas, 0);
                      // Rcout << unk16 << std::endl;

                      auto cmax = (lenremain + alignval) / (alignval+8);
                      if (debug) Rcout << cmax << std::endl;

                      std::vector<CN_Att> cnpois(cmax);

                      for (auto i = 0; i < cmax; ++i) {

                        if (ALIGN_1_VALUE == 4) {
                          cnpois[i].CN_OFF     = readbin(cnpois[i].CN_OFF, sas, 0);
                        } else {
                          cnpois[i].CN_OFF     = readbin((int32_t)cnpois[i].CN_OFF, sas, 0);
                        }
                        cnpois[i].CN_WID     = readbin(cnpois[i].CN_WID, sas, 0);
                        cnpois[i].NM_FLAG    = readbin(cnpois[i].NM_FLAG, sas, 0);
                        cnpois[i].CN_TYP     = readbin(cnpois[i].CN_TYP, sas, 0);
                        cnpois[i].UNK8       = readbin(cnpois[i].UNK8, sas, 0);

                        auto wid = cnpois[i].CN_WID;
                        auto typ = cnpois[i].CN_TYP;

                        // Rprintf("WID: %d\n", wid  );
                        // Rprintf("typ: %d\n", typ );

                        if (typ > 0) {
                          colwidth.push_back( wid );
                          vartyps.push_back( typ );
                        }

                      }
                    }

                    else
                    {
                      sas_hex = int32_to_hex(sas_offset);
                      Rcout << sas_hex << std::endl;

                      std::string unkstr (potabs[sc].SH_LEN, '\0');

                      unkstr = readstring(unkstr, sas);

                    }

      }
    }

    // ---------------------------------------------------------------------- //

    // if (pagecount > 1) {
    //
    //   auto page2 = headersize + pagesize;
    //   sas.seekg(page2, sas.beg);
    //
    //   if (debug) Rcout << "page 2: " << page2 << std::endl;
    //
    //   // Page Offset Table
    //
    //   // int64_t pageseqnum = 0;
    //
    //   if (ALIGN_2_VALUE == 4) {
    //     pageseqnum = readbin(pageseqnum, sas, 0);
    //     unk64 = readbin(unk64, sas, 0);
    //     if (debug) Rcout << unk64 << std::endl;
    //     unk64 = readbin(unk64, sas, 0);
    //     if (debug) Rcout << unk64 << std::endl;
    //     unk64 = readbin(unk64, sas, 0);
    //     if (debug) Rcout << unk64 << std::endl;
    //   } else {
    //     pageseqnum = readbin((int32_t)pageseqnum, sas, 0);
    //     unk32 = readbin(unk32, sas, 0);
    //     if (debug) Rcout << unk32 << std::endl;
    //     unk32 = readbin(unk32, sas, 0);
    //     if (debug) Rcout << unk32 << std::endl;
    //     unk32 = readbin(unk32, sas, 0);
    //     if (debug) Rcout << unk32 << std::endl;
    //   }
    //   if (debug) Rprintf("pageseqnum: %d \n", pageseqnum);
    //
    //   PAGE_TYPE = 0, BLOCK_COUNT = 0, SUBHEADER_COUNT = 0;
    //   PAGE_TYPE = readbin(PAGE_TYPE, sas, 0);
    //   BLOCK_COUNT = readbin(BLOCK_COUNT, sas, 0);
    //   SUBHEADER_COUNT = readbin(SUBHEADER_COUNT, sas, 0);
    //   unk16 = readbin(unk16, sas, 0);
    //
    //   rowsperpage[1] = BLOCK_COUNT - SUBHEADER_COUNT;
    //
    //   // if (debug)
    //   Rprintf("PAGE_TYPE: %d ; BLOCK_COUNT: %d ; SUBHEADER_COUNT: %d \n",
    //           PAGE_TYPE, BLOCK_COUNT, SUBHEADER_COUNT);
    //
    //   if (debug) Rprintf("unk16: %d \n", unk16);
    //
    //
    //   for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
    //     if (ALIGN_2_VALUE == 4) {
    //
    //       potabs[i].SH_OFF = readbin(potabs[i].SH_OFF, sas, 0);
    //       potabs[i].SH_LEN = readbin(potabs[i].SH_LEN, sas, 0);
    //       potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, 0);
    //       potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, 0);
    //
    //       zero = readbin(zero, sas, 0);
    //       zero = readbin(zero, sas, 0);
    //       zero = readbin(zero, sas, 0);
    //       zero = readbin(zero, sas, 0);
    //       zero = readbin(zero, sas, 0);
    //       zero = readbin(zero, sas, 0);
    //
    //       // if (debug)
    //       Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPRESSION: %d ; SH_TYPE: %d \n",
    //               potabs[i].SH_OFF, potabs[i].SH_LEN,
    //               potabs[i].COMPRESSION, potabs[i].SH_TYPE);
    //
    //     } else {
    //
    //       potabs[i].SH_OFF = readbin((int32_t)potabs[i].SH_OFF, sas, 0);
    //       potabs[i].SH_LEN = readbin((int32_t)potabs[i].SH_LEN, sas, 0);
    //       potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, 0);
    //       potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, 0);
    //
    //       zero = readbin(zero, sas, 0);
    //       zero = readbin(zero, sas, 0);
    //
    //       // if (debug)
    //       Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPRESSION: %d ; SH_TYPE: %d \n",
    //               potabs[i].SH_OFF, potabs[i].SH_LEN,
    //               potabs[i].COMPRESSION, potabs[i].SH_TYPE);
    //     }
    //   }
    //
    //
    //
    //   if (PAGE_TYPE == 256 || PAGE_TYPE == 512) {
    //     sh_end_pos = sas.tellg();
    //     // debug
    //     Rprintf("position: %d\n", sh_end_pos);
    //
    //     data_pos.push_back( sh_end_pos );
    //
    //   }
    //
    // }


    // "FDFFFFFF"

    // auto dl = (position + 7) % 8 * 8;

    // Rcout << dl << std::endl;

    //
    //       int32_t pos = sas.tellg();
    //       Rprintf("pos: %d\n", pos);
    //
    //       // B = 16
    //       // DL = (B+8+SC*SL+PL+7) % 8 * 8
    //       int32_t dl = 0;
    //       dl = (16 + 8+sc*12 + pagesize + 7) % 8 * 8;
    //
    //       Rprintf("dl: %d \n", dl);
    //
    //       // sas.seekg(dl, ios::cur);
    //
    //       Rprintf("pos: %d \n", sas.tellg());
    // }
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


    // // new offset ------------------------------------------------------------//
    // pos = pagestart + potabs[].SH_OFF;


    auto page = 0;
    sas.seekg(data_pos[0], sas.beg);

    std::string row(rowlength, '\0');
    for (auto i = 0; i < rowcount; ++i) {

      // row = readstring(row, sas);

      // auto rowsonpage = rowsperpage[page];
      // auto rowsread = rowcount - i;

      // if ((pagecount > 0) & (page < pagecount)) {
      //   auto rowsonpage = rowsperpage[page];
      //   sas.seekg(data_pos[page], sas.beg);
      //
      //   ++page;
      //
      //   Rcout << page << std::endl;
      // }

      auto pos = 0;
      for (auto j = 0; j < colnum; ++j) {

        auto wid = colwidth[j];
        auto typ = vartyps[j];

        if (wid == 8 & typ == 1) {

          double val_d = 0.0;

          val_d = readbin(val_d, sas, 0);

          // Rcout << val_d << std::endl;

          REAL(VECTOR_ELT(df,j))[i] = val_d;

        }

        if (wid > 0 & typ == 2) {

          std::string val_s(wid, '\0');

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

      // Rcout << row << std::endl;



      if (i == (rowsperpage[page] -1) ) {
        ++page;
        sas.seekg(data_pos[page], sas.beg);
      }

    }




    IntegerVector rvec = seq(1, rowcount);

    // 3. Create a data.frame
    df.attr("row.names") = rvec;
    df.attr("names") = varnames;
    df.attr("class") = "data.frame";


    // close file
    sas.close();

    return(df);

  } else {
    return (-1);
  }
}
