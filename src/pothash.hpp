#ifndef POTHASH_HPP
#define POTHASH_HPP

#include <cmath>
#include <cstdint>
#include <map>
#include <string>

using namespace std;

#define NEG_INF numeric_limits<float>::lowest()
#define PI acos(-1.0)

struct Vector2D
{
    float vectorComponentX = 0;
    float vectorComponentY = 0;

    float magnitude() const
    {
        return sqrt(this->vectorComponentX * this->vectorComponentX + this->vectorComponentY * this->vectorComponentY);
    }

    Vector2D normalize() const
    {
        return Vector2D{this->vectorComponentX / this->magnitude(), this->vectorComponentY / this->magnitude()};
    }

    Vector2D operator+(const Vector2D& other) const
    {
        return Vector2D{this->vectorComponentX + other.vectorComponentX, this->vectorComponentY + other.vectorComponentY};
    }
};

struct DpTableCell
{
    Vector2D dpTableCellScoreVector;
    uint64_t dpTableCellSeed;
};

struct Seed
{
    int seedParamD;
    Vector2D seedScoreVector;
    string seedSubsequence;
};

class PotHash
{
    int paramD;
    int subseqSizeK;
    map<char, int> alphabet;

    int *tableThetaDegrees;

    Vector2D *tableA;
    int *tableB;

    Vector2D projectVectorOnRay(const Vector2D&, const int);
    
public:
    PotHash();
    PotHash(const int, const int, const map<char, int>&);
    ~PotHash();

    void generateTables();

    Seed solveDP(const string&);
};

#endif