#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H

std::string SASYZCRL(uint64_t rowlen, uint64_t reslen, std::string rowstr,
                     bool debug);

std::string SASYZCR2(uint64_t rowlen, uint64_t reslen, std::string rowstr,
                     bool debug);

#endif
