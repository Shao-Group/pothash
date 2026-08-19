#ifndef POTHASH_HPP
#define POTHASH_HPP

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace std;

inline constexpr float NEG_INF = numeric_limits<float>::lowest();
#define PI acos(-1.0)

struct DpTableCell
{
    float dpTableCellScore;
    uint64_t dpTableCellSeed;
};

struct Seed
{
    int seedParamD;
    float seedScore;
    string seedSubsequence;
};

class PotHash
{
    int paramD;
    int subseqSizeK;
    map<char, int> alphabet;

    // ThetaA[k, sigma]: angle in degrees
    float *tableA;
    // ThetaB[k, sigma, d]: angle in degrees
    float *tableB;

    map<int, char> reverseAlphabet;

    // P(r, theta, thetaPrime) = r * cos(theta - thetaPrime); angles in degrees
    float projectScalar(const float, const float, const float) const;

    int indexA(const int, const int) const;
    int indexB(const int, const int, const int) const;
    // DP cell T[n,k,d,lastBase]
    int indexDP(const int, const int, const int, const int) const;

    string generateSeedSubsequence(const uint64_t&);
    
public:
    PotHash();
    PotHash(const int, const int, const map<char, int>&);
    ~PotHash();

    void generateTables(const bool, const int);
    void loadTables(const int);

    vector<Seed> solveDP(const string&, const uint8_t&);
};

#endif
