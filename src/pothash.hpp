#ifndef POTHASH_HPP
#define POTHASH_HPP

#include <map>
#include <string>

using namespace std;

class PotHash
{
private:
    int paramD;
    int subseqSizeK;
    map<char, int> alphabet;

    int *tableThetaDegrees;

    int *tableAx;
    int *tableAy;
    int *tableB;

public:
    PotHash();
    PotHash(const int, const int, const map<char, int>&);
    ~PotHash();

    void generateTables();

    void solveDP(const string&);
};

#endif