#include <iostream>
#include "MurmurHash3.h"
#include <sdsl/bit_vectors.hpp>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <time.h>
#include <unistd.h>
#include <chrono>
#include "bbhash.hpp"

using namespace std;
using namespace sdsl;


bbhash::bbhash(){}
bbhash::bbhash(vector<string>* keys_in, int num_keys, double gamma){
    /*
        inputs: 
        keys_in - a pointer to vector of strings, which is the list of keys for our table
        num_keys - the number of keys in the keyset
        gamma - tuning parameter that effects size of bitarrays
        
        output: a BBhash object, with 3 bitarrays & a hash table to store our keys
    */


    gamma = gamma;
    keys = keys_in;
    num_input = num_keys;
    
    // size of bitvectors is effected by gamma parameter
    int size = int (gamma * num_input);
    //fill our first bitvector, remove strings that will map to A0
    tuple<bit_vector,int> fillzero = arrayfill(size, seed_0);
    A0 = get<0>(fillzero);
    A0_rank = rank_support_v<1> (&A0);
    //the "rank" of vector A0
    rank_A0 = A0_rank(size);

    //repeat this process two more times
    
    int keys1 = get<1>(fillzero);
    size = int (gamma * keys1);
    tuple<bit_vector,int> fillone = arrayfill(size, seed_1);
    A1 = get<0>(fillone);
    A1_rank = rank_support_v<1> (&A1);
    //the "rank" of vector A1
    rank_A1 = A1_rank(size);

    int keys2 = get<1>(fillone);
    size = int (gamma * keys2);
    tuple<bit_vector,int> filltwo = arrayfill(size, seed_2);
    A2 = get<0>(filltwo);
    A2_rank = rank_support_v<1> (&A2);
    //the "rank" of vector A2
    rank_A2 = A2_rank(size);

    //all keys left are put in a hash table
    int keys_rem = get<1>(filltwo);
    ht = htfill(keys_rem,num_input);
}
unordered_map<string,uint64_t> bbhash::htfill(int key_rem,int num_keys){
    /*
        inputs:
        
        keys_rem - how many keys are left
        num_keys - the number of keys we started with when creating the class
        
        output:
        ht - a dictionary that maps each key to a value in the range (num_keys − keys_rem + 1, num_keys) inclusive
    */

    unordered_map<string,uint64_t> ht;
    //starting point
    int index = num_keys - key_rem + 1;

    //just loop through each key left and map it to a number in the range specified
    for (string key: (*keys)){ 
        pair<string,uint64_t> p(key,uint64_t(index));
        ht.insert(p);
        index ++;
    }
    return ht;

}
tuple<bit_vector,int> bbhash::arrayfill(int size, int seed){
    /*
        inputs:
        size - desired size of the bitarrays
        seed - seed for xxhash 
        
        output:
        a bitvector array that has been "filled", the rank of the bitvector, the remaining number of keys left
     */
    //arr, the bitarray that will be returned, values should be 1 iff a key has been mapped to that location with no collisions
    bit_vector arr = bit_vector(size);
    //cout << size << endl;
    //cout << arr.size() << endl;
    //c - the bitarray which tells us if there is a collision at the same point in arr
    bit_vector c = bit_vector(size);
    //number of keys that have a collision
    int keys_left = 0;
    
    int filled = 0;
    int alc = 0;
    int newc= 0;
    uint64_t * out = new uint64_t[2];
    for(vector<string>::iterator it = keys->begin() ; it != keys->end(); it++){
        string key = *it;
        //hashing the key
        int len = sizeof(key.c_str());
        MurmurHash3_x64_128(key.c_str(),len,seed,out);
        //position in arr
        int pos = int( out[0] % size);
        //if c[pos] == 1, there is a collision at pos already
        if (c[pos] == 1){
            arr[pos] = 0;
            alc++;
            //keys_left++;
        }
        //otherwise, if arr[pos] = 0 & c[pos] = 0, there has not been a key mapped here yet
        else if(arr[pos] == 0 && c[pos] == 0){
            //cout << "Mapping to " << pos << endl;
            arr[pos] = 1;
            filled++;
        }
        //there is already a key mapped here, now we have a collision
        else if(arr[pos] == 1 && c[pos] == 0){
            arr[pos] = 0;
            c[pos] = 1;
            newc++;
            //keys_left += 2;
        }
    }
    //cout << filled << "\t" << newc << "\t" << alc << endl;
    int kl = 2 * newc + alc;
    //loop through keys again, removing keys that have no collisions
    uint64_t * out2 =  new uint64_t[2];
    for (vector<string>::iterator it2 = keys->begin() ; it2 != keys->end();){
        //hashing the key
        string key = *it2;
        int len = sizeof(key.c_str());
        //const void * key_p = key.c_str();
        
        MurmurHash3_x64_128(key.c_str(),len,seed,out2);
        //position in arr
        int pos = int( out2[0] % size);
        //if it has no collisions, value in arr should be 1
        if(arr[pos] == 1){
            //vector<int>::iterator position = find(keys.begin(), keys.end(), 8);
            keys->erase(it2);
            //rem++;
        }
        //otherwise the key should remain, increment pointer to vector 
        else{
            keys_left++;
            it2++;
        }
    }
    //cout << rem << endl;
    //number of keys left
    int num_left = keys->size();
    //cout << num_left << endl;
    return tuple<bit_vector,int> (arr,num_left);
}

uint64_t bbhash::query(string key){
    /*
    input: string key, the key which we are searching for
    output: uint64_t, the value our key maps to. Will be in (1,N) for keys in the set. Will be n + 1234567 for keys not in the set.
    */

    //cummulative rank of the tables we've passed
    int cum_rank = 0;
    
    int i = 0;
    vector<int> seeds = {seed_0,seed_1,seed_2};
    vector<int> sizes = {int(A0.size()), int(A1.size()), int(A2.size())};
    vector<int> ranks = {rank_A0,rank_A1,rank_A2};

    uint64_t out[2]= {0};
    int len = sizeof(key.c_str());
    while(i <= 2){
        //hash the key with the ith seed
        int seed = seeds[i];
        MurmurHash3_x86_32(key.c_str(),len,seed,out);
        //get position in current table
        int pos = out[0] % sizes[i];

        //bitvecs is array of pointers to our bitvectors, if Ai[position]  == 1, we found a "match"
        if( (*bitvecs[i])[pos] == 1 ){
            //return the sum of ranks of vectors already visited + rank in current vector
            // rank function is not inclusive of current positition
            return uint64_t(cum_rank + (*rank_sups[i])(pos + 1));
        }
        else{
            cum_rank += ranks[i];
            i++;
        }

    }

    //if not in tables, return value in our hash table
    std::unordered_map<string,uint64_t>::const_iterator search = ht.find(key);
    if(search == ht.end()){
        return num_input + 1234567;
    }
    else{
        return search->second; 
    }
}

uint64_t bbhash::size_bv(){
    return A0.size() + A1.size() + A2.size();
}

uint64_t bbhash::size_ht(){
    return ht.size();
}

void testfile(string file, string outfile_size, string outfile_speed, int num_keys){
    /*
    inputs: 
    string file - the file name that contains the keys
    string outfile_size - file to write size results to
    string outfile_speed - file to write speed results to
    int num_keys - the number of keys in the file

    */
    
    
    //read file lines into vector
    vector<string> keys;
    //create a list of randomly selected indices (for vector testkeys)
    vector<int> indices;

    string line;
    ifstream in(file); 
    while (getline(in, line)){
        if (line.size() > 0){
            keys.push_back(line);
            indices.push_back(rand() % num_keys);
        }   
    }

    in.close();
    //values of gamma to test
    double gammas[5] = {1,2,4,6,8};
    vector<string> testkeys;
    uint64_t times = 0;
    ofstream speed_results;
    ofstream size_results;
    //test each value of gamma by runing num_keys * 10 queries on a bbhash object made from testkeys
    
    for (double gamma: gammas){
        testkeys = keys;
        bbhash bb = bbhash(&testkeys,num_keys,gamma);
        size_results.open(outfile_size,ios_base::app);
        //Size of bitvectors, hashtable for gamma and num_keys
        size_results << gamma << "," << num_keys << "," << bb.size_bv() << "," << bb.size_ht() << endl;
        size_results.close();
        // querying randomly chosen keys in our keyset
        for(int i = 0;i<10;i++){
            for(int index: indices){
                uint64_t start = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                bb.query(keys[index]);
                uint64_t end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                uint64_t time = end-start;
                times = times + time;
            }
        }
        //convert time to milliseconds
        double avg_time = double(times) / 1000000;
        double num_queries = 10.0 * num_keys;
        
        speed_results.open(outfile_speed,ios_base::app);
        //Avg lookup time for gamma, num_keys
        speed_results << gamma << ", " << num_keys << "," << avg_time / num_queries << endl;
        speed_results.close();
        times = 0;
    }
    

}
void testalien(string file, string alienfile,string outfile, int num_keys){
    //read file lines into vector
    vector<string> keys;
    //create a list of randomly selected indices (for vector testkeys)
    vector<int> indices;

    string line;
    ifstream in(file); 
    while (getline(in, line)){
        if (line.size() > 0){
            keys.push_back(line);
            indices.push_back(rand() % num_keys);
        }   
    }
    in.close();

    //read alien keys into a vector
    vector<string> aliens;
    ifstream in2(alienfile); 
    while (getline(in2, line)){
        if (line.size() > 0){
            aliens.push_back(line);
        }   
    }
    //values of gamma to test
    double gammas[5] = {1,2,4,6,8};
    vector<string> testkeys;
    uint64_t times = 0;
    ofstream results;

    int num_alien = 0;
    for(double gamma: gammas){
        testkeys = keys;
        bbhash bb = bbhash(&testkeys,num_keys,gamma);
        for(string alien_key:aliens){
            //since we return n + 1234567 for alien keys, just have to check if we return a value > n
            if (bb.query(alien_key) > num_keys){
                num_alien ++;
            }
        }
        results.open(outfile,ios_base::app);
        //write % of alien keys detected
        results << gamma << "," << num_keys << "," << num_alien / 1000.0 << endl;
        num_alien = 0;
        results.close();

    }
}
int main(int argc, char *argv[]){
    //get size of index
    string command = argv[1];
    int num_keys = stoi(argv[2]);
    //get file name
    string file = argv[3];
    //first command executes speed and size tests
    //second executes alien key test
    if(command == "speed_size"){
        //files for results
        string out_size = argv[4];
        string out_speed = argv[5];
        testfile(file,out_size,out_speed,num_keys);
    }
    else if(command == "alien"){
        string aliens = argv[4];
        string outfile = argv[5];
        testalien(file,aliens,outfile,num_keys);
    }
    

    return 0;
    
}