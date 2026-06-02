#ifndef POTHASH_HPP
#define POTHASH_HPP

#include <cmath>
#include <map>
#include <string>

using namespace std;

#define NEG_INF numeric_limits<float>::lowest()
#define PI acos(-1)

struct Vector2D
{
public:
    float vectorComponentX;
    float vectorComponentY;

    Vector2D() : Vector2D(0, 0)
    {
        // Calling parameterized constructor from default constructor with default arguments
    }

    Vector2D(const float vectorComponentX, const float vectorComponentY)
    {
        this->vectorComponentX = vectorComponentX;
        this->vectorComponentY = vectorComponentY;
    }

    Vector2D operator+(const Vector2D& other) const
    {
        return Vector2D(this->vectorComponentX + other.vectorComponentX, this->vectorComponentY + other.vectorComponentY);
    }

    float magnitude() const
    {
        return sqrt(this->vectorComponentX * this->vectorComponentX + this->vectorComponentY * this->vectorComponentY);
    }

    Vector2D normalize() const
    {
        return Vector2D(this->vectorComponentX / this->magnitude(), this->vectorComponentY / this->magnitude());
    }
};

class PotHash
{
private:
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

    void solveDP(const string&);
};

#endif