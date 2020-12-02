# bbhash

Implemenation of BBHash algorithm - CMSC858D.
To run, you should have [sdsl-lite](https://github.com/simongog/sdsl-lite) installed on your machine. Credit for MurmurHash3 algorithm implentation goes to [Austin Appleby](https://github.com/aappleby/smhasher/tree/master/src). If you wish to change things and recompile, recompile with ```g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib bbhash.cpp MurmurHash3.cpp -o bbhash -lsdsl -ldivsufsort -ldivsufsort64```. Otherwise, you can run the speed and size tests or alien key tests.

Speed and size tests:

``` ./bbhash speed_size <num_keys> <text file containing keys> <output file for size results> <output file for speed results> ```

Alien key test:

``` ./bbhash alien <num_keys> <text file containing keys> <text file containing alien keys (assumed to have 1000 keys)> <output file for alien key detection results> ```

