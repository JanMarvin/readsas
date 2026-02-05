#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H

/*
 * Both functions are c++ conversions of java functions from the parso library
 * Licensed under the Apache License, Version 2.0. Copyright 2015 EPAM
 *
 */

#include <math.h>
#include "sas.h"
#include "uncompress.h"

std::string SASYZCRL(uint64_t rowlen, uint64_t reslen, const std::string& rowstr, bool debug) {
    std::string res;
    res.reserve(reslen);

    const uint8_t* row = reinterpret_cast<const uint8_t*>(rowstr.data());
    uint32_t rowoff = 0;

    while (rowoff < rowlen) {
        uint8_t control = row[rowoff];
        uint8_t cbyte = control & 0xF0;
        uint8_t ebyte = control & 0x0F;

        rowoff++;

        switch (cbyte) {
            case 0x00: case 0x10: case 0x20: case 0x30: { // Large Literal Copy
                if (rowoff < rowlen) {
                    int32_t len = (row[rowoff] & 0xFF) + 64 + (control << 8);
                    rowoff++;
                    if (rowoff + len <= rowlen) {
                        res.append(rowstr, rowoff, len);
                        rowoff += len;
                    }
                }
                break;
            }
            case 0x40: { // RLE: Repeat row[rowoff+1]
                if (rowoff < rowlen) {
                    int32_t count = (ebyte << 8) + (row[rowoff] & 0xFF) + 18;
                    rowoff++;
                    if (rowoff < rowlen) {
                        res.append(count, (char)row[rowoff]);
                        rowoff++;
                    }
                }
                break;
            }
            case 0x60: { // Repeat Space (0x20)
                int32_t count = (ebyte << 8) + (row[rowoff] & 0xFF) + 17;
                res.append(count, ' ');
                rowoff++;
                break;
            }
            case 0x70: { // Repeat Zero (0x00)
                int32_t count = (ebyte << 8) + (row[rowoff] & 0xFF) + 17;
                res.append(count, '\0');
                rowoff++;
                break;
            }
            case 0x80: case 0x90: case 0xA0: case 0xB0: { // Small Literal Copy
                int32_t len = (control - 0x7F);
                if (rowoff + len <= rowlen) {
                    res.append(rowstr, rowoff, len);
                    rowoff += len;
                }
                break;
            }
            case 0xC0: { // Short RLE
                int32_t count = ebyte + 3;
                if (rowoff < rowlen) {
                    res.append(count, (char)row[rowoff]);
                    rowoff++;
                }
                break;
            }
            case 0xD0: { // Short Repeat '@' (0x40)
                res.append(ebyte + 2, '@');
                break;
            }
            case 0xE0: { // Short Repeat Space (0x20)
                res.append(ebyte + 2, ' ');
                break;
            }
            case 0xF0: { // Short Repeat Zero (0x00)
                res.append(ebyte + 2, '\0');
                break;
            }
            default:
                break;
        }
    }

    if (res.size() != reslen)
        Rcpp::warning("SASYZCRL: string size %zu != %zu\n", res.size(), reslen);

    return res;
}

std::string SASYZCR2(uint64_t rowlen, uint64_t reslen, const std::string& rowstr, bool debug) {
    std::string res;
    res.reserve(reslen);

    const uint8_t* row = reinterpret_cast<const uint8_t*>(rowstr.data());
    uint32_t rowoff = 0;
    uint32_t cbit = 0, cmsk = 0;

    while (rowoff < rowlen && res.size() < reslen) {
        cmsk >>= 1;
        if (cmsk == 0) {
            if (rowoff + 1 >= rowlen) break;
            // Read 16-bit control mask
            cbit = (static_cast<uint16_t>(row[rowoff]) << 8) | row[rowoff + 1];
            rowoff += 2;
            cmsk = 0x8000;
        }

        // Control bit 0: Raw Literal Byte
        if ((cbit & cmsk) == 0) {
            if (rowoff < rowlen) {
                res += static_cast<char>(row[rowoff++]);
            }
            continue;
        }

        if (rowoff >= rowlen) break;
        uint8_t ctrl = row[rowoff++];
        uint8_t cmd = (ctrl >> 4) & 0x0F;
        uint8_t len_nibble = ctrl & 0x0F;

        switch (cmd) {
            case 0: { // Short RLE
                uint16_t count = len_nibble + 3;
                if (rowoff < rowlen) {
                    char b = static_cast<char>(row[rowoff++]);
                    uint32_t can_add = std::min<uint32_t>(count, reslen - res.size());
                    res.append(can_add, b);
                }
                break;
            }
            case 1: { // Long RLE
                if (rowoff < rowlen) {
                    uint16_t count = (static_cast<uint16_t>(row[rowoff++]) << 4) + len_nibble + 19;
                    if (rowoff < rowlen) {
                        char b = static_cast<char>(row[rowoff++]);
                        uint32_t can_add = std::min<uint32_t>(count, reslen - res.size());
                        res.append(can_add, b);
                    }
                }
                break;
            }
            case 2: { // Long LZ77 (Pattern)
                if (rowoff + 1 < rowlen) {
                    uint32_t ofs = len_nibble + 3 + (static_cast<uint32_t>(row[rowoff++]) << 4);
                    uint16_t count = static_cast<uint16_t>(row[rowoff++]) + 16;

                    if (ofs <= res.size()) {
                        uint32_t pos = res.size() - ofs;
                        uint32_t max_copy = std::min<uint32_t>(count, reslen - res.size());
                        for (uint32_t i = 0; i < max_copy; ++i)
                            res += res[pos + i];
                    }
                }
                break;
            }
            default: { // Short LZ77 (cmd >= 3)
                if (rowoff < rowlen) {
                    uint32_t ofs = len_nibble + 3 + (static_cast<uint32_t>(row[rowoff++]) << 4);
                    uint16_t count = cmd;

                    if (ofs <= res.size()) {
                        uint32_t pos = res.size() - ofs;
                        uint32_t max_copy = std::min<uint32_t>(count, reslen - res.size());
                        for (uint32_t i = 0; i < max_copy; ++i)
                            res += res[pos + i];
                    }
                }
                break;
            }
        }
    }

    if (res.size() != reslen) {
        Rcpp::warning("SASYZCR2 mismatch: %zu != %zu", res.size(), reslen);
        if (res.size() < reslen) res.append(reslen - res.size(), '\0');
        else if (res.size() > reslen) res.resize(reslen);
    }

    return res;
}

#endif
