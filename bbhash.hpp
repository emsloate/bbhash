//
//  rank_support.hpp
//  
//
//  Created by Elliott Sloate on 10/22/20.
//

#ifndef rank_support_hpp
#define rank_support_hpp

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "MurmurHash3.h"
#include <sdsl/bit_vectors.hpp>

class bbhash
{
 
public:
    std::vector<std::string> * keys;
    int num_input;
    double gamma;

    sdsl::bit_vector A0;
    sdsl::rank_support_v<1> A0_rank;
    int rank_A0;

    sdsl::bit_vector A1;
    sdsl::rank_support_v<1> A1_rank;
    int rank_A1;
    
    sdsl::bit_vector A2;
    sdsl::rank_support_v<1> A2_rank;
    int rank_A2;

    int seed_0 = 777;
    int seed_1 = 888;
    int seed_2 = 999;
    std::unordered_map<std::string,uint64_t> ht;

    sdsl::bit_vector* bitvecs[3] = {&A0,&A1,&A2};
    sdsl::rank_support_v<1>* rank_sups[3] = {&A0_rank,&A1_rank,&A2_rank};


    bbhash();
    bbhash(std::vector<std::string> *keys_in, int num_keys, double gamma);
    uint64_t query(std::string key);
    std::tuple<sdsl::bit_vector,int> arrayfill(int size, int seed);
    std::unordered_map<std::string,uint64_t> htfill(int key_rem,int num_keys);
    uint64_t size_bv();
    uint64_t size_ht();
    void save(std::string& fname);
    void load(std::string& fname);
};
#endif /* rank_support_hpp */
