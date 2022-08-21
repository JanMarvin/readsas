/*
 * Copyright (C) 2019, 2022 Jan Marvin Garbuszus
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
#include <bitset>

#include "sas.h"
#include "uncompress.h"

using namespace Rcpp;


//' Reads SAS data files
//'
//' @param filePath The full systempath to the sas7bdat file you want to import.
//' @param debug print debug information
//' @param selectrows integer vector of selected rows
//' @param selectcols character vector of selected rows
//' @import Rcpp
//' @export
// [[Rcpp::export]]
Rcpp::List readsas(const char * filePath,
                   const bool debug,
                   Nullable<IntegerVector> selectrows_,
                   Nullable<CharacterVector> selectcols_)
{
  std::ifstream sas(filePath, std::ios::in | std::ios::binary | std::ios::ate);
  auto sas_size = sas.tellg();
  if (sas) {

    sas.seekg(0, std::ios_base::beg);

    const std::string tempstr = ".readsas_unc_tmp_file";
    std::fstream out (tempstr, std::ios::out | std::ios::binary);


    int compr = 0;

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
    uint8_t encoding = 0;
    int16_t dataoffset = 0;
    int16_t swlen = 0, proclen = 0, comprlen = 0, textoff = 0, todata = 0,
      addtextoff = 0, fmtkey = 0, fmtkey2 = 0,
      fmt32 = 0, fmt322 = 0, ifmt32 = 0, ifmt322 = 0;
    int16_t PAGE_TYPE = 0, BLOCK_COUNT = 0, SUBHEADER_COUNT = 0;

    uint32_t pageseqnum32 = 0;

    double created = 0, created2 = 0;   // 8
    double modified = 0, modified2 = 0; // 16

    std::string compression = "";
    std::string compstr = "";
    std::string proc = "";
    std::string enc = "";
    std::string sw = "";
    std::string sasfile  (8, '\0');
    std::string filetype (8, '\0');
    std::string sasrel   (8, '\0');
    std::string sasserv  (16, '\0');
    std::string osver    (16, '\0');
    std::string osmaker  (16, '\0');
    std::string osname   (16, '\0');
    std::string dataset  (64, '\0');

    std::vector<CN_Poi> cnpois;
    std::vector<idxofflen> fmt;
    std::vector<idxofflen> lbl;
    std::vector<idxofflen> unk;
    std::vector<int16_t> c8vec;

    Rcpp::NumericVector fmt32s;
    Rcpp::NumericVector ifmt32s;
    Rcpp::NumericVector fmtkeys;
    Rcpp::IntegerVector vartyps;
    Rcpp::IntegerVector colwidth;
    Rcpp::IntegerVector page_type;
    Rcpp::IntegerVector coloffset;
    Rcpp::IntegerVector cnidx;
    Rcpp::IntegerVector cnoff;
    Rcpp::IntegerVector cnlen;
    Rcpp::IntegerVector cnzer;

    Rcpp::CharacterVector labels;
    Rcpp::CharacterVector formats;
    Rcpp::CharacterVector varnames; // (k)

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
    encoding = readbin(encoding, sas, swapit); // file encoding?
    if (debug) Rprintf("%d\n", encoding);
    enc = SASEncoding(encoding);
    unk8 = readbin(unk8, sas, swapit); // sys encoding?
    if (debug) Rprintf("%d\n", unk8);


    if (debug) Rcout << " ---- block ---- " << sas.tellg()  << std::endl;
    /* begin block of 4 ----------------------------------------------------- */
    unk8 = readbin(unk8, sas, swapit); // 0
    if (debug) Rprintf("%d\n", unk8);
    unk8 = readbin(unk8, sas, swapit); // 16 | 32: PAGE_BIT_OFFSET?
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
    dataset = std::regex_replace(dataset, std::regex(" +$"), "$1");
    if (debug) Rcout << dataset << std::endl;

    // o156 filetype 'DATA    '
    filetype = readstring(filetype, sas);
    filetype = std::regex_replace(filetype, std::regex(" +$"), "$1");
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

    uint32_t headersize = 0;
    headersize = readbin(headersize, sas, swapit);
    if (debug) Rprintf("headersize: %d \n", headersize);
    if (headersize <= 0) stop("headersize <= 0");

    uint32_t pagesize = 0;
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

    /*
     * theoretically every page contains data and/or varnames. practically
     * this must not be true. 1024 does not contain data, only varnames. 512
     * might contain both or only varnames.
     */

    std::vector<uint64_t> data_pos(pagecount, 0);
    std::vector<uint64_t> varname_pos;
    std::vector<uint64_t> label_pos;
    std::vector<int64_t>  rowsperpage(pagecount, 0);

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

    for (uint64_t i = 0; i < num_zeros; ++i) {
      int8_t zero = 0;
      zero = readbin(zero, sas, swapit);
      if (zero!=0 && debug)
        warning("debug: expected 0, is %d", zero);
    }

    // debug
    int64_t pagestart = sas.tellg();
    if (debug) Rprintf("position: %d\n", pagestart);

    // end of Header ---------------------------------------------------------//


    uint8_t alignval = 8;
    if (u64 != 4) alignval = 4;

    uint64_t rowlength = 0, n = 0, delobs = 0;
    int64_t colf_p1 = 0, colf_p2 = 0;
    int64_t k = 0;
    std::vector<std::string> stringvec(pagecount) ;

    auto totalrows = 0;
    std::vector<uint32_t> totalrowsvec(pagecount);
    std::vector<uint32_t> pageseqnum(pagecount);

    std::vector<std::string> pagedelmarker(pagecount);

    int8_t PAGE_BIT_OFFSET = 0;
    int8_t SUBHEADER_POINTER_LENGTH = 0;
    int8_t SUBHEADER_POINTERS_OFFSET = 8;

    if (u64 == 4) {
      PAGE_BIT_OFFSET = 32;
      SUBHEADER_POINTER_LENGTH = 24;
    } else {
      PAGE_BIT_OFFSET = 16;
      SUBHEADER_POINTER_LENGTH = 12;
    }

    uint64_t pre_pagenumx = 0;
    uint64_t pagenumx = 0;

    // begin reading pages ---------------------------------------------------//
    for (auto pg = 0; pg < pagecount; ++pg) {
      checkUserInterrupt();

      // Rcout << "--- new page ------------------------------------" << std::endl;

      /* should already be at this position for pg == 1 */
      if (pagecount > 0) {
        pre_pagenumx = pagenumx;
        pagenumx = headersize + (double)pg * pagesize;

        sas.seekg(pagenumx, sas.beg);

        if (pagenumx <= pre_pagenumx)
          stop("pagenumx did not increase");
      }

      int64_t unk1 = 0, unk2 = 0;
      int32_t c5typ = 0; /* check if variable name, format or label */
      int64_t PAGE_DELETED_POINTER_LENGTH = 0;

      // Page Offset Table
      if (u64 == 4) {
        pageseqnum32 = readbin(pageseqnum32, sas, swapit);
        unk32 = readbin(unk32, sas, swapit);
        unk1 = readbin(unk1, sas, swapit);
        unk2 = readbin(unk2, sas, swapit);
        PAGE_DELETED_POINTER_LENGTH = readbin(PAGE_DELETED_POINTER_LENGTH, sas, swapit);
      } else {
        pageseqnum32 = readbin(pageseqnum32, sas, swapit);
        unk1 = readbin((int32_t)unk1, sas, swapit);
        unk2 = readbin((int32_t)unk2, sas, swapit);
        PAGE_DELETED_POINTER_LENGTH = readbin((int32_t)PAGE_DELETED_POINTER_LENGTH, sas, swapit);
      }

      pageseqnum[pg] = pageseqnum32;

      PAGE_TYPE       = readbin(PAGE_TYPE, sas, swapit);
      BLOCK_COUNT     = readbin(BLOCK_COUNT, sas, swapit);
      SUBHEADER_COUNT = readbin(SUBHEADER_COUNT, sas, swapit);
      unk16           = readbin(unk16, sas, swapit);

      page_type.push_back(PAGE_TYPE);

      rowsperpage[pg] = BLOCK_COUNT - SUBHEADER_COUNT;
      totalrows += rowsperpage[pg];
      totalrowsvec[pg] = totalrows;

      if (debug)
        Rprintf("PAGE_TYPE: %d ; BC: %d ; SC: %d ; unk16: %d ---- \n",
                PAGE_TYPE, BLOCK_COUNT, SUBHEADER_COUNT, unk16);

      int16_t zero = 0;
      // int64_t sh_tot_len = 0;
      uint64_t dataoff = 0;

      std::vector<PO_Tab> potabs(SUBHEADER_COUNT);


      if ((
          PAGE_TYPE == 16384 ||                   // PAGE_META_TYPE_2
            PAGE_TYPE == 1024 ||                    // PAGE_AMD_TYPE
            PAGE_TYPE == 640 || PAGE_TYPE == 512 || // PAGE_MIX_TYPE_2   PAGE_MIX_TYPE_1
            PAGE_TYPE == 384 || PAGE_TYPE == 256 || // PAGE_DATA_TYPE_2  PAGE_DATA_TYPE
            PAGE_TYPE == 128 ||                     // PAGE_CMETA_TYPE
            PAGE_TYPE == 0))                        // PAGE_META_TYPE_1
      {

        for (auto i = 0; i < SUBHEADER_COUNT; ++i) {
          if (u64 == 4) {

            potabs[i].SH_OFF = readbin(potabs[i].SH_OFF, sas, swapit);           // 8
            potabs[i].SH_LEN = readbin(potabs[i].SH_LEN, sas, swapit);           // 16
            potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, swapit); // 20
            potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, swapit);         // 24

            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;
            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;
            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;

          } else {

            potabs[i].SH_OFF = readbin((uint32_t)potabs[i].SH_OFF, sas, swapit); // 4
            potabs[i].SH_LEN = readbin((uint32_t)potabs[i].SH_LEN, sas, swapit); // 8
            potabs[i].COMPRESSION = readbin(potabs[i].COMPRESSION, sas, swapit); // 10
            potabs[i].SH_TYPE = readbin(potabs[i].SH_TYPE, sas, swapit);         // 12

            zero = readbin(zero, sas, swapit);
            // Rcout << zero << std::endl;
          }

          if (debug)
            Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPR.: %d ; SH_TYPE: %d \n",
                    potabs[i].SH_OFF, potabs[i].SH_LEN,
                    potabs[i].COMPRESSION, potabs[i].SH_TYPE);

          dataoff = potabs[i].SH_OFF - ((double)rowsperpage[pg] * rowlength);

        }

        if (debug) {
          Rcout << "data offset " << dataoff << std::endl;
          Rcout << "data offset ------------------------------- " << std::endl;
        }

        uint64_t sh_end_pos = 0;

        if (PAGE_TYPE != 0) sh_end_pos = sas.tellg();

        data_pos[pg] = sh_end_pos;

        if (debug)
          // Rprintf("sh_end_pos: %d\n", sh_end_pos);
          Rcout << "sh_end_pos: " << sh_end_pos << std::endl;


        // from now on, we will seek to every position inside the sas file
        for (auto sc = 0; sc < SUBHEADER_COUNT; ++sc)
        {

          if (debug)
            Rcout << "Subheader Count: " << sc << std::endl;

          uint64_t pagepos = (headersize + (double)pg * pagesize) + potabs[sc].SH_OFF;

          // 2 files, where this is a problem
          if(potabs[sc].SH_OFF == 0 || potabs[sc].SH_LEN == 0)
            break;

          sas.seekg(pagepos, sas.beg);

          // not sure yet, whats the right thing to do here
          bool page0not = 0;
          if ((pg == 0) & (sc != 3))
            page0not = (PAGE_TYPE == 0) &  (potabs[sc].SH_LEN == rowlength);

          int64_t sas_offset = alignval;
          if (! ((potabs[sc].COMPRESSION == 4) |
              (PAGE_TYPE == -28672) | page0not ) ) {
            if (u64 == 4) {
              sas_offset = readbin(sas_offset, sas, swapit);
            } else {
              sas_offset = readbin((int32_t)sas_offset, sas, swapit);
            }
          }

          std::string sas_hex = int_to_hex(sas_offset);

          if (debug)
            Rcout << "SAS Hex: " << sas_hex << std::endl;

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
          if (potabs[sc].COMPRESSION == 4)
            sas_offset_table = 9;
          if (page0not)
            sas_offset_table = 10;


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

              if (debug)
                Rcout << "-------- case 1 "<< sas.tellg() << std::endl;



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
                n = readbin(n, sas, swapit);
                if (debug) Rcout << "n " << n << std::endl;
                delobs = readbin(delobs, sas, swapit);
                if (debug) Rcout << "delobs " << delobs << std::endl;
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

                /* next two indicate the end of the initial header ? */
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl; // -1
                unk64 = readbin(unk64, sas, swapit);
                if (debug) Rcout << unk64 << std::endl; // -1

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
                  if (unk64 != 0 && debug)
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
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 1
                if (debug) Rcout << unk16 << std::endl;

                sh_num = readbin(sh_num, sas, swapit);
                if (debug) Rcout << sh_num << std::endl;
                cn_maxlen = readbin(cn_maxlen, sas, swapit);
                if (debug) Rcout << cn_maxlen << std::endl;
                l_maxlen = readbin(l_maxlen, sas, swapit);
                if (debug) Rcout << l_maxlen << std::endl;

                /* maybe SAS version information at o131018 ? */
                unk32 = readbin(unk32, sas, swapit); // 1
                // Rcout << "1 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 2
                // if (unk32 != 0) stop("unk32 1. expected 0 is %d", unk32);
                unk32 = readbin(unk32, sas, swapit); // 3
                // if (unk32 != 0) stop("unk32 1. expected 0 is %d", unk32);

                rowsonpg = readbin(rowsonpg, sas, swapit);


                unk16 = readbin(unk16, sas, swapit); // 1
                if (unk16 != 0) stop("unk16 01. expected 0 is %d", unk16);
                unk32 = readbin(unk32, sas, swapit); // 2
                // Rcout << "2 " << unk32 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 4
                if (unk16 != 0) stop("unk16 04. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 5
                if (unk16 != 0) stop("unk16 05. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 6
                if (unk16 != 0) stop("unk16 06. expected 0 is %d", unk16);
                unk32 = readbin(unk32, sas, swapit); // 7
                // Rcout << "7 " << unk32 << std::endl; // nrows
                unk16 = readbin(unk16, sas, swapit); // 9
                if (unk16 != 0) stop("unk16 09. expected 0 is %d", unk16);
                unk32 = readbin(unk32, sas, swapit); // 10
                // Rcout << "10 "<< unk32 << std::endl; // delobs
                unk16 = readbin(unk16, sas, swapit); // 12
                if (unk16 != 0) stop("unk16 12. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 13
                if (unk16 != 0) stop("unk16 13. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 14
                if (unk16 != 0) stop("unk16 14. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 15
                if (unk16 != 0) stop("unk16 15. expected 0 is %d", unk16);
                dataoffset = readbin(dataoffset, sas, swapit); // 16
                // Rcout << dataoffset << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 17
                if (unk16 != 0) stop("unk16 17. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit);// 18
                if (unk16 != 0) stop("unk16 18. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 19
                if (unk16 != 0) stop("unk16 19. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 20
                if (unk16 != 0) stop("unk16 20. expected 0 is %d", unk16);


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
                if (debug) Rcout << "rowlength " << rowlength << std::endl;
                n = readbin((int32_t)n, sas, swapit);
                if (debug) Rcout << "rowcount " << n << std::endl;
                delobs = readbin((int32_t)delobs, sas, swapit); // deleted obs?
                if (debug) Rcout << "delobs " << delobs << std::endl;
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
                  if (unk64 != 0 && debug)
                    warning("val1 is %d. expected a zero", unk64);
                }

                pgidx = readbin(pgidx, sas, swapit);


                for (int z = 0; z < 8; ++z) {
                  unk64 = readbin((int32_t)unk64, sas, swapit);
                  if (debug) Rcout << unk64 << std::endl;
                }

                // padding?
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit);
                if (debug) Rcout << unk32 << std::endl;

                unk32 = readbin(unk32, sas, swapit); // val 1?
                if (debug) Rcout << unk32 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // val 2?
                if (debug) Rcout << unk16 << std::endl;

                unk16 = readbin(unk16, sas, swapit); // padding
                if (debug) Rcout << unk16 << std::endl;

                pgwsh = readbin((int32_t)pgwsh, sas, swapit);
                if (debug) Rcout << "pgwsh " << pgwsh << std::endl;
                pgwpossh = readbin(pgwpossh, sas, swapit);
                if (debug) Rcout << "pgwpossh " << pgwpossh << std::endl;

                unk16 = readbin(unk16, sas, swapit); // padding
                if (debug) Rcout << unk16 << std::endl;

                pgwsh2 = readbin((int32_t)pgwsh2, sas, swapit);
                if (debug) Rcout << "pgwsh2 " << pgwsh2 << std::endl;
                pgwpossh2 = readbin(pgwpossh2, sas, swapit);
                if (debug) Rcout << "pgwpossh2 " << pgwpossh2 << std::endl;

                unk16 = readbin(unk16, sas, swapit); // padding
                if (debug) Rcout << unk16 << std::endl;

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


                unk32 = readbin(unk32, sas, swapit); // 1
                // Rcout << "1 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 2
                // Rcout << "2 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 3
                // Rcout << "3 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 4
                // Rcout << "4 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 5
                // Rcout << "5 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 6
                // Rcout << "6 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 7
                // Rcout << "7 " << unk32 << std::endl;
                unk32 = readbin(unk32, sas, swapit); // 8
                // Rcout << "8 " << unk32 << std::endl;

                unk16 = readbin(unk16, sas, swapit); // 4
                if (debug) Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 1
                if (debug) Rcout << unk16 << std::endl;

                sh_num = readbin(sh_num, sas, swapit);
                if (debug) Rcout << sh_num << std::endl;
                cn_maxlen = readbin(cn_maxlen, sas, swapit);
                if (debug) Rcout << cn_maxlen << std::endl;
                l_maxlen = readbin(l_maxlen, sas, swapit);
                if (debug) Rcout << l_maxlen << std::endl;

                for (int z = 0; z < 3; ++z) {
                  unk32 = readbin(unk32, sas, swapit); // 0
                }

                rowsonpg = readbin(rowsonpg, sas, swapit);


                unk16 = readbin(unk16, sas, swapit); // 1
                if (unk16 != 0) stop("unk16 01. expected 0 is %d", unk16);
                unk32 = readbin(unk32, sas, swapit); // 2
                // Rcout << "2 " << unk32 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 4
                if (unk16 != 0) stop("unk16 04. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 5
                if (unk16 != 0) stop("unk16 05. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 6
                if (unk16 != 0) stop("unk16 06. expected 0 is %d", unk16);
                unk32 = readbin(unk32, sas, swapit); // 7
                // Rcout << "7 " << unk32 << std::endl; // nrows
                unk16 = readbin(unk16, sas, swapit); // 9
                if (unk16 != 0) stop("unk16 09. expected 0 is %d", unk16);
                unk64 = readbin((int32_t)unk64, sas, swapit); // 10
                if (debug) Rcout << "delobs "<< unk64 << std::endl; // delobs?
                // if (unk16 != 0) stop("unk16 11. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 12
                if (unk16 != 0) stop("unk16 12. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 13
                if (unk16 != 0) stop("unk16 13. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 14
                if (unk16 != 0) stop("unk16 14. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 15
                if (unk16 != 0) stop("unk16 15. expected 0 is %d", unk16);
                dataoffset = readbin(dataoffset, sas, swapit); // 16
                // Rcout << dataoffset << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 17
                if (unk16 != 0) stop("unk16 17. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit);// 18
                if (unk16 != 0) stop("unk16 18. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 19
                if (unk16 != 0) stop("unk16 19. expected 0 is %d", unk16);
                unk16 = readbin(unk16, sas, swapit); // 20
                if (unk16 != 0) stop("unk16 20. expected 0 is %d", unk16);
              }

              if (debug)
                Rprintf("swlen = %d, todata %d, textoff %d\n",
                        swlen, todata, textoff);


              if (!((dataoffset == 1) ||(dataoffset == 256) || (dataoffset == 1280)))
                warning("debug: dataoffset is unexpectedly %d\n",
                        dataoffset);

              break;
            }



            // new offset --------------------------------------------------- //
          case 2:
            {

              if (debug)
                Rcout << "-------- case 2 "<< sas.tellg() << std::endl;

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

              if (debug)
                Rcout << "-------- case 3 "<< sas.tellg() << std::endl;

              hasattributes = 1;

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

              if ((unks.IDX != 0) | (unks.OFF != 0) | (unks.LEN != 0)) {
                warning("case3: unk is not 0 as expected, but %d %d %d\n",
                        unks.IDX, unks.OFF, unks.LEN);
              }

              break;
            }


            // new offset --------------------------------------------------- //
          case 4:
            {
              /* Column Size */

              if (debug)
                Rcout << "-------- case 4 "<< sas.tellg() << std::endl;

              uint64_t uunk64 = 0;

              if (u64 == 4) {
                k = readbin(k, sas, swapit);
                uunk64 = readbin(uunk64, sas, swapit);
              } else {
                k = readbin((int32_t)k, sas, swapit);
                uunk64 = readbin((int32_t)uunk64, sas, swapit);
              }

              if (debug)
                Rprintf("k %d; uunk64 %d\n",
                        k, uunk64);

              break;
            }

            // new offset --------------------------------------------------- //

          case 5:
            {
              /* Column Text */

              if (debug)
                Rcout << "-------- case 5 "<< sas.tellg() << std::endl;

              int16_t len = 0;

              varname_pos.push_back( sas.tellg() );

              len = readbin(len, sas, swapit);
              // Rcout << len << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              if ((PAGE_TYPE != 1024) & (c5first == 0)) {
                unk16 = readbin(unk16, sas, swapit); // 0 |     0 | 27977
                // Rcout << unk16 << std::endl;
                unk16 = readbin(unk16, sas, swapit); // 0 | 15872 | 30064
                // Rcout << unk16 << std::endl;
              }

              if ((c5typ == 0) & (pg == 0)) {

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
                if (proclen > 0) {
                  proc.resize(proclen, '\0');
                  proc = readstring(proc, sas);
                }

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

                sas.seekg(pos_beg, sas.beg);
              }

              if (debug)
                Rprintf("SH_LEN %d; len %d;\n",
                        potabs[sc].SH_LEN, len);

              ++c5typ;

              break;
            }


            // new offset --------------------------------------------------- //
          case 6:
            {
              /* Column Name */

              if (debug)
                Rcout << "-------- case 6 "<< sas.tellg() << std::endl;


              int16_t lenremain = 0;

              lenremain = readbin(lenremain, sas, swapit);
              if (debug) Rprintf("lenremain %d \n", lenremain);

              int8_t div = 8;
              lenremain -= 8;

              auto cmax = lenremain / div;


              unk16 = readbin(unk16, sas, swapit); // 0
              if (unk16 != 0) stop("unk16 1. expected 0 is %d", unk16);
              unk16 = readbin(unk16, sas, swapit); // 0
              if (unk16 != 0) stop("unk16 2. expected 0 is %d", unk16);
              unk16 = readbin(unk16, sas, swapit); // 0
              if (unk16 != 0) stop("unk16 3. expected 0 is %d", unk16);

              /* Column Name Pointers */
              CN_Poi cnpoi;

              for (auto cn = 0; cn < cmax; ++cn) {

                cnpoi.CN_IDX    = readbin(cnpoi.CN_IDX, sas, swapit);
                cnpoi.CN_OFF    = readbin(cnpoi.CN_OFF, sas, swapit);
                cnpoi.CN_LEN    = readbin(cnpoi.CN_LEN, sas, swapit);
                cnpoi.zeros     = readbin(cnpoi.zeros,  sas, swapit);

                cnpois.push_back( cnpoi );

                cnidx.push_back( cnpoi.CN_IDX );
                cnoff.push_back( cnpoi.CN_OFF );
                cnlen.push_back( cnpoi.CN_LEN );
                cnzer.push_back( cnpoi.zeros );

                if (debug) {
                  Rprintf("CN_IDX %d; CN_OFF %d; CN_LEN %d; zeros %d; len %d;\n",
                          cnpois[cn].CN_IDX,
                          cnpois[cn].CN_OFF,
                          cnpois[cn].CN_LEN,
                          cnpois[cn].zeros,
                          cn);
                }


              }

              break;
            }


            // new offset --------------------------------------------------- //
          case 7:
            {
              /* Column Attributes */

              if (debug)
                Rcout << "-------- case 7 "<< sas.tellg()  << std::endl;

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

              int8_t divs = 16;
              if (u64 != 4) divs = 12;

              lenremain -= 8;

              auto cmax = lenremain / divs;

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


                if ((capois[i].CN_TYP >= 1) & (capois[i].CN_TYP <= 2) &
                    (capois[i].CN_WID >= 0) & // just > ?
                    ((uint32_t)capois[i].CN_WID <= pagesize)) {
                  if (debug)
                    Rprintf("OFF %d; WID: %d; FLAG %d; TYP %d; UNK8 %d\n",
                            capois[i].CN_OFF, capois[i].CN_WID,
                            capois[i].NM_FLAG,
                            capois[i].CN_TYP, capois[i].UNK8 );

                  coloffset.push_back( capois[i].CN_OFF );
                  colwidth.push_back( capois[i].CN_WID );
                  vartyps.push_back( capois[i].CN_TYP );
                }

              }

              break;
            }

          case 8:
            {

              if (debug)
                Rcout << "-------- case 8 "<< sas.tellg() << std::endl;

              int16_t cls = 0;
              int64_t lenremain = 0;

              unk32 = readbin(unk32, sas, swapit); // unkown large number
              // Rcout << unk32 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin(unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              if (u64 == 4) {  // lenremain
                lenremain = readbin(lenremain, sas, swapit);
              } else {
                lenremain = readbin((int32_t)lenremain, sas, swapit);
              }

              if (debug)
                Rcout << "lenremain "<< lenremain << std::endl; // 92

              unk16 = readbin(unk16, sas, swapit);
              // Rcout << unk16 << std::endl; // number of varnames?
              cls = readbin(cls, sas, swapit);
              // Rcout << cls << std::endl;   // counter for unk loop below
              unk16 = readbin(unk16, sas, swapit);
              // Rcout << unk16 << std::endl; // 1
              unk16 = readbin(unk16, sas, swapit);
              // Rcout << unk16 << std::endl; // number of varnames?
              unk16 = readbin(unk16, sas, swapit);  // 3233
              // Rcout << unk16 << std::endl; // 0
              unk16 = readbin(unk16, sas, swapit);  // 3233
              // Rcout << unk16 << std::endl; // 0
              unk16 = readbin(unk16, sas, swapit);  // 3233
              // Rcout << unk16 << std::endl; // 0

              lenremain -= 14;

              // Rcout << lenremain << " " << cls << std::endl;

              for (auto cl = 0; cl < cls; ++cl) {
                int16_t res = 0;
                res = readbin(res, sas, swapit);
                c8vec.push_back(res);
              }


              // Rcout << "---------------------------" << std::endl;

              // 8
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;
              unk16 = readbin((int8_t)unk16, sas, swapit); // 0
              // Rcout << unk16 << std::endl;

              // Rcout << "---------------------------" << std::endl;

              break;

            }

          case 9:
            {

              if (debug)
                Rcout << "-------- case 9 "<< sas.tellg()  << std::endl;


              auto clen = potabs[sc].SH_LEN;
              std::string cstr(clen, '\0');
              std::string ustr(rowlength, '\0');
              cstr = readstring(cstr, sas);

              if (compr == 1)
                ustr = SASYZCRL(clen, rowlength, cstr, debug);

              if (compr == 2)
                ustr = SASYZCR2(clen, rowlength, cstr, debug);

              writestr(ustr, ustr.size(), out);

              break;
            }

          case 10:
            {
              if (debug)
                Rcout << "-------- case 10 "<< sas.tellg()  << std::endl;

              if ((potabs[sc].SH_LEN > alignval) &
                  (potabs[sc].SH_LEN < pagesize))
              {
                // uncompressed row containing data
                auto unklen = potabs[sc].SH_LEN;
                std::string unkstr(unklen, '\0');
                unkstr = readstring(unkstr, sas);
                writestr(unkstr, unkstr.size(), out);
              }

              break;
            }

            // not implemented ---------------------------------------------- //
          default:
            {

              if (debug) {
              Rcout << "---- unimplemented "<< sas.tellg() << std::endl;
              Rcout << "SAS HEX STRING: "  << sas_hex << std::endl;
              Rcout << "rowlength is " << rowlength << std::endl;

              Rprintf("SH_OFF: %d ; SH_LEN: %d ; COMPR.: %d ; SH_TYPE: %d \n",
                      potabs[sc].SH_OFF, potabs[sc].SH_LEN,
                      potabs[sc].COMPRESSION, potabs[sc].SH_TYPE);

            }

              // some subheaders are pointers to positions inside the file.
              // their use for SAS is unknown and they are not required for R.
              // if SC_LEN == alignval it is just padding?
              if ((potabs[sc].SH_LEN > alignval) &
                  (potabs[sc].SH_LEN < pagesize) &
                  (potabs[sc].COMPRESSION != 4) )
              {
                auto unklen = potabs[sc].SH_LEN;
                std::string unkstr(unklen, '\0');
              }

              break;
            }

          }
        }



        /*
         * The code to detect removed rows is based on the parso library
         * Licensed under the Apache License, Version 2.0. Copyright 2015 EPAM
         *
         */

        // check for deleted rows
        if (PAGE_TYPE == 384 || PAGE_TYPE == 640 || PAGE_TYPE == 1024) {

          uint64_t start_pos = sas.tellg();

          // TODO this calculation should be replaced with alignval assignment
          auto alignCorrection = (
            (PAGE_BIT_OFFSET + 8) +
              SUBHEADER_POINTERS_OFFSET +
              (int64_t)((double)SUBHEADER_COUNT * SUBHEADER_POINTER_LENGTH)
          ) % 8;

          if (debug)
            Rcout << "alignCorrection: " << alignCorrection << std::endl;

          uint64_t deletedMapOffset = (PAGE_BIT_OFFSET + 8) +
            PAGE_DELETED_POINTER_LENGTH +
            alignCorrection + // alignval?
            ((double)SUBHEADER_COUNT * SUBHEADER_POINTER_LENGTH) +
            ((double)rowsperpage[pg] * rowlength);

          uint64_t dmo_pos = (headersize + (double)pg * pagesize) + deletedMapOffset;

          if (debug) {
            Rcout << "SUBHEADER_COUNT " << SUBHEADER_COUNT << std::endl;
            Rcout << "rowlength " << rowlength << std::endl;
            Rcout << "dMO " << deletedMapOffset << std::endl;
            Rcout << "read from: " << dmo_pos << std::endl;
          }

          sas.seekg(dmo_pos, sas.beg);


          // every byte contains information for 8 rows
          int32_t dm_len = (int32_t)std::ceil((double)rowsperpage[pg] / 8);
          if (debug) Rcout << "dm_len: " << dm_len << std::endl;

          std::string delmarker = "";
          for (auto dm = 0; dm < dm_len; ++dm) {
            unk8 = readbin(unk8, sas, swapit);
            if (debug) Rprintf("unk8: %d\n", unk8);

            delmarker += std::bitset<8>(unk8).to_string(); //to binary
            if (debug) Rcout << delmarker << std::endl;
          }

          pagedelmarker[pg] = delmarker;

          sas.seekg(start_pos, sas.beg);

        }

      }
    }

    out.close();


    if (debug)
      Rcout << "varnames ----------------------------" << std::endl;

    for (size_t i = 0; i < cnpois.size(); ++i) {

      if (debug)
        Rprintf("CN_IDX %d; CN_OFF %d; CN_LEN %d; zeros %d \n",
                cnpois[i].CN_IDX, cnpois[i].CN_OFF,
                cnpois[i].CN_LEN, cnpois[i].zeros);

      uint64_t vpos = (varname_pos[cnpois[i].CN_IDX] + cnpois[i].CN_OFF);
      sas.seekg(vpos, sas.beg);

      std::string varname(cnpois[i].CN_LEN, '\0');
      varname = readstring(varname, sas);

      if ((int16_t)varname.size() == cnpois[i].CN_LEN)
        varnames.push_back(varname);

      if (debug)
        Rcout << i << " : " << vpos << " : " << varname << std::endl;

    }




    if (hasattributes) {

      for (auto i = 0; i < k; ++i) {

        /* read formats and labels */
        std::string format = "";
        if ((size_t)i < fmt.size() && fmt[i].LEN > 0) {
          uint64_t fpos = (varname_pos[fmt[i].IDX] + fmt[i].OFF);
          sas.seekg(fpos, sas.beg);
          format.resize(fmt[i].LEN, '\0');
          format = readstring(format, sas);
        }

        std:: string label = "";
        if ((size_t)i < lbl.size() && lbl[i].LEN > 0) {
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

    if ((compr < 0) | (compr > 2))
      warning("File contains unhandled compression. No data imported. %d\n",
              compression);

    if (debug)
      Rcout << "---- data part ----" << std::endl;


    // ---------------------------------------------------------------------- //

    uint64_t nmin = 0, nmax = 0;
    uint64_t nn   = 0;

    // if  selectrows is c(0,0) use full data
    if (selectrows_.isNull()) {
      nmin = 1;
      nmax = n;
    } else {
      IntegerVector selectrows(selectrows_);
      if (selectrows.size() != 2) stop("selected rows must be vector of two.");
      nmin = selectrows(0);
      nmax = selectrows(1);
    }

    // make sure that n is not greater than nmax or nmin
    if (n < nmax)
      nmax = n;
    if (n < nmin)
      nmin = n;

    if (debug)
      Rcout << "reading nmin/nmax: " << nmin << " / " << nmax << std::endl;

    // sequences of column and row
    IntegerVector cvec = seq(0, (k-1));
    IntegerVector rvec = seq(nmin, nmax);
    nn = rvec.size();


    /* this is dangerous, the length is correct, but the position in the file
     * differs. right now it is not always exactly this. To make this working
     * modifications to the import process are required.
     */
    // rowlength
    IntegerVector rlen = rowlength;

    // check if vars are selected
    bool selectvars = selectcols_.isNotNull();

    // select vars: either select every var or only matched cases. This will
    // return index positions of the selected variables. If non are selected the
    // index position is cvec
    IntegerVector select = cvec, nselect;
    if (selectvars) {
      CharacterVector selectcols(selectcols_);
      select = choose(selectcols, varnames);
    }

    IntegerVector sels = select;

    // separate the selected from the not selected cases
    LogicalVector ll = is_na(select);
    nselect = cvec[ll == 1];
    select = cvec[ll == 0];

    // get new cvec_
    cvec = cvec_(ll);

    uint32_t kk = select.size();

    // shrink variables to selected size
    CharacterVector varnames_kk = clone(varnames)[select];
    IntegerVector vartyps_kk = clone(vartyps)[select];


    // 2. fill it with data

    // 1. Create Rcpp::List
    Rcpp::List df(kk);
    for (uint32_t i = 0; i < kk; ++i)
    {
      int32_t const type = vartyps_kk[i];

      switch(type)
      {
      case 1:
        SET_VECTOR_ELT(df, i, NumericVector(no_init(nn)));
        break;

      default:
        SET_VECTOR_ELT(df, i, CharacterVector(no_init(nn)));
      break;
      }
    }

    // IntegerVector coloffset_kk = clone(coloffset);
    // coloffset_kk[ll == 1] = -1;
    //
    // Rf_PrintValue(coloffset_kk);
    Rcpp::IntegerVector ordered = order_(coloffset);

    if (debug)
    {
      Rcout << "ll "        << ll           << std::endl;
      Rcout << "vartyps "   << vartyps      << std::endl;
      Rcout << "vartyps "   << vartyps_kk   << std::endl;
      Rcout << "coloffset " << coloffset    << std::endl;
      // Rcout << "coloffset " << coloffset_kk << std::endl;
      Rcout << "ordered "   << ordered      << std::endl;
    }

    // new offset ----------------------------------------------------------- //

    std::vector<bool> deleted(n);
    std::vector<bool> valid(n);

    bool firstpage = 0;
    bool keepr = 0;

    // sas provides two modes, compressed and uncompressed data. compressed
    // data has to be uncompressed and consists of always single rows. un-
    // compressed data is unordered aka might have gaps. So we have to be more
    // careful and handle certain offsets.
    // uncompressed data might contain deleted rows. most likely these are not
    // available in compressed data. BUT THIS IS A GUESS.
    if (compr == 0) {

      Rcout << "no compression" << std::endl;

      auto page = 0;
      sas.seekg(data_pos[0], sas.beg);

      auto i = -1;  // counter output data frame
      uint64_t ii = 0;  // row on the selected page
      for (uint64_t iii = 0; iii < n; ++iii) {

        /* nmin is not a c vector starting at 0. i is initialized at -1 so will
         * be 0 once its bigger than nmin. This allows to import only the
         * selected rows. Once nmax is reached, import will stop.
         */
        if (iii >= (nmin-1) ) {
          keepr = 1;
          ++i;
        }

        // Rcout << "---------------------" << std::endl;
        // Rcout << iii << " " << nmin << std::endl;
        if (debug && i == 0)
          Rcout << "row i / ii / iii / keepr: " << i << " " << ii << " " << iii <<" " << keepr << std::endl;

        if (iii >= nmax) break;

        if (pagecount > 0) {

          while (totalrowsvec[page] == 0) {
            ++page;
            ii = 0;

            if (page == pagecount)
              break;
          }

          if (debug && i == 0)
            Rcout << "page: " << page << std::endl;

          // beg handle delmarker
          std::string page_pagedelmarker = pagedelmarker[page];

          // Rcout << page_pagedelmarker << " : " << ii << std::endl;


          std::string delmarked = "";
          if (!page_pagedelmarker.empty() && (uint64_t)page_pagedelmarker.size() > 0 && ii < page_pagedelmarker.size()) {
            delmarked = page_pagedelmarker[ii];
            // Rcout << page_pagedelmarker[ii] << std::endl;
          }
          if (delmarked == "1")
            deleted[iii] = true;
          else
            deleted[iii] = false;

          valid[iii] = true;
          // end handle delmarker


          if (totalrowsvec[page] == ii) {
            ++page;
            ii = 0;
            firstpage = 1;
          }
        }

        uint64_t pp = data_pos[page];
        uint64_t pos = pp + (double)rowlength * ii;

        /* unknown */
        if (!(dataoffset == 1 || dataoffset == 256) & (firstpage == 0)) {
          pos += alignval;
        }


        // check if position is equal to expected position
        if (pos != (uint64_t)sas.tellg()) {
          // Rcout << "shift" << std::endl;
          sas.seekg(pos, sas.beg);
        }

        pos = 0;

        for (auto j = 0; j < k; ++j) {

          auto ord = ordered[j];
          auto wid = colwidth[ord];
          auto typ = vartyps[ord];
          auto col = cvec[ord];

          if (debug && i == 0)
            Rcout << "ord/wid/typ: " << ord << " : " << wid << " : " << typ << std::endl;

          bool keepc = col >= 0;

          // Rcout << "---------------------" << std::endl;
          if (debug && i == 0)
            Rcout << "  --col j / ord / keepc: " << j << " " << ord << " "<< keepc << std::endl;

          if (debug && i == 0)
            Rcout << ord << std::endl;

          if ((wid < 8) & (typ == 1)) {

            double val_d = 0.0;

            val_d = readbinlen(val_d, sas, 0, wid);

            if (debug && i == 0)
              Rcout << val_d << std::endl;

            if (debug && i == 0)
              Rcout << "writing: " << i << "/" << col << std::endl;

            if (keepr) {
              if (keepc) {
                if (std::isnan(val_d))
                  REAL(VECTOR_ELT(df,col))[i] = NA_REAL;
                else
                  REAL(VECTOR_ELT(df,col))[i] = val_d;
              }
            }

          }

          if ((wid == 8) & (typ == 1)) {

            double val_d = 0.0;

            val_d = readbin(val_d, sas, swapit);

            if (debug && i == 0)
              Rcout << val_d << std::endl;

            if (debug && i == 0)
              Rcout << "writing: " << i << "/" << col << std::endl;

            if (keepr) {
              if (keepc) {
                if (std::isnan(val_d))
                  REAL(VECTOR_ELT(df,col))[i] = NA_REAL;
                else
                  REAL(VECTOR_ELT(df,col))[i] = val_d;
              }
            }

          }

          if ((wid > 0) & (typ == 2)) {

            std::string val_s(wid, ' ');

            val_s = readstring(val_s, sas);

            val_s = std::regex_replace(val_s,
                                       std::regex(" +$"), "$1");

            if (debug && i == 0)
              Rcout << val_s << std::endl;

            if (debug && i == 0)
              Rcout << "writing: " << i << "/" << col << std::endl;

            if (keepr)
              if (keepc) {
                as<CharacterVector>(df[col])[i] = val_s;
              }

          }

        }

        // check if eof is reached sas.eof() did not work
        if ((sas_size - sas.tellg()) == 0) {
          // Rcout << "eof reached" << std::endl;
          break;
        }

        ++ii;
      }
    }

    // close file. compressed data is imported from a different file
    sas.close();


    if ((compr == 1) || (compr == 2)) {

      Rcout << "compression" << std::endl;

      std::ifstream sas(tempstr, std::ios::in | std::ios::binary | std::ios::ate);

      if (sas) {

        sas_size = sas.tellg();
        sas.seekg(0, std::ios_base::beg);

        auto i = -1;
        auto ii = 0;
        for (uint64_t iii = 0; iii < n; ++iii) {

          if (debug && i == 0)
            Rcout << "row: " << i << " --------------------" <<std::endl;

          valid[iii] = true;

          /* nmin is not a c vector starting at 0. i is initialized at -1 so will
           * be 0 once its bigger than nmin. This allows to import only the
           * selected rows. Once nmax is reached, import will stop.
           */
          if (iii >= (nmin-1) ) {
            keepr = 1;
            ++i;
          }

          for (auto j = 0; j < k; ++j) {

            auto ord = ordered[j];
            auto wid = colwidth[ord];
            auto typ = vartyps[ord];
            auto col = cvec[ord];

            bool keepc = col >= 0;

            if (debug && i == 0)
              Rcout << ord << " : " << wid << " : " << typ << std::endl;


            if ((wid > 0) & (wid < 8) & (typ == 1)) {

              double val_d = 0.0;

              val_d = readbinlen(val_d, sas, 0, wid);

              if (debug && i == 0)
                Rcout << val_d << std::endl;

              if (keepc) {
                if (std::isnan(val_d))
                  REAL(VECTOR_ELT(df,col))[i] = NA_REAL;
                else
                  REAL(VECTOR_ELT(df,col))[i] = val_d;
              }

            }

            if ((wid == 8) & (typ == 1)) {

              double val_d = 0.0;

              val_d = readbin(val_d, sas, swapit);

              if (debug && i == 0)
                Rcout << val_d << std::endl;

              if (keepc) {
                if (std::isnan(val_d))
                  REAL(VECTOR_ELT(df,col))[i] = NA_REAL;
                else
                  REAL(VECTOR_ELT(df,col))[i] = val_d;
              }

            }

            if ((wid > 0) & (typ == 2)) {

              std::string val_s(wid, ' ');

              val_s = readstring(val_s, sas);

              val_s = std::regex_replace(val_s,
                                         std::regex(" +$"), "$1");

              if (debug && i == 0)
                Rcout << val_s << std::endl;

              if (keepc)
                as<CharacterVector>(df[col])[i] = val_s;

            }
          }

          // check if eof is reached sas.eof() did not work
          if ((sas_size - sas.tellg()) == 0) {
            // Rcout << "eof reached" << std::endl;
            break;
          }

          ++ii;
        }

      } else {
        warning("could not read from compressed pages");
      }

      sas.close();

    }

    // Rf_PrintValue(df);

    if (!debug) {
      if ( remove(tempstr.c_str()) != 0 )
        warning("Could not remove temporary file containing",
                "uncompressed data");
    }

    if (debug) {
      Rprintf("%d %d \n", nn, kk);
    }

    if (nn > 0)
      rvec = seq(1, nn);

    if (varnames.size() > kk)
      cvec = seq(1, kk);

    // 3. Create a data.frame
    df.attr("row.names") = rvec;

    if (varnames_kk.size() == kk)
      df.attr("names") = varnames_kk;
    else
      df.attr("names") = cvec;
    df.attr("class") = "data.frame";

    if (varnames.size() > kk)
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
    df.attr("encoding") = enc;
    df.attr("fmtkeys") = fmtkeys;
    df.attr("fmt32") = fmt32s;
    df.attr("ifmt32") = ifmt32s;

    df.attr("rowcount") = nn;
    df.attr("rowlength") = rowlength;
    df.attr("deleted_rows") = delobs;
    df.attr("colwidth") = colwidth;
    // df.attr("coloffset") = coloffset;
    df.attr("vartyps") = vartyps;
    df.attr("c8vec") = c8vec;

    df.attr("headersize") = headersize;
    df.attr("pagesize") = pagesize;

    // if (kk >= 0 && (uint64_t)kk == nn) {
    //   // deleted needs to be reduced if we ever skip the first row
    //   deleted.resize(nn);
    // }

    df.attr("cvec") = cvec;
    df.attr("rvec") = rvec;
    df.attr("deleted") = deleted;
    df.attr("valid") = valid;

    if (debug) {
      df.attr("cnidx") = cnidx;
      df.attr("cnoff") = cnoff;
      df.attr("cnlen") = cnlen;
      df.attr("cnzer") = cnzer;
    }


    return(df);

  } else {
    return (-1);
  }
}
