#include "pothash.hpp"

#include <algorithm>
#include <chrono>
#include <random>

PotHash::PotHash() : PotHash(64, 16, {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}})
{
    // Calling parameterized constructor from default constructor with default arguments
}

PotHash::PotHash(const int paramD, const int subseqSizeK, const map<char, int>& alphabet)
{
    this->paramD = (paramD % 2 == 0) ? paramD : paramD + 1;
    this->subseqSizeK = subseqSizeK;
    this->alphabet = alphabet;

    this->tableThetaDegrees = new int[paramD];

    this->tableAx = new int[subseqSizeK * alphabet.size()];
    this->tableAy = new int[subseqSizeK * alphabet.size()];
    this->tableB = new int[subseqSizeK * alphabet.size()];
}

PotHash::~PotHash()
{
    delete[] this->tableThetaDegrees;

    delete[] this->tableAx;
    delete[] this->tableAy;
    delete[] this->tableB;
}

void PotHash::generateTables()
{
    random_device seedSource;

    mt19937 pseudoRandomNumberEngine(seedSource());

    uniform_int_distribution<int> uniformIntegerDistributionTableTheta(0, 360);

    vector<int> randomThetaDegrees;

    for (int i = 0; i < this->paramD / 2; i++)
    {
        randomThetaDegrees.push_back(uniformIntegerDistributionTableTheta(pseudoRandomNumberEngine));
    }

    for (int i = this->paramD / 2; i < this->paramD; i++)
    {
        randomThetaDegrees.push_back(360 - randomThetaDegrees[i - this->paramD / 2]);
    }

    unsigned int currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();

    shuffle(randomThetaDegrees.begin(), randomThetaDegrees.end(), default_random_engine(currentTimeBasedSeed));

    for (int d = 0; d < this->paramD; d++)
    {
        // tableThetaDegrees contains random ray angles in degrees to avoid storing floating-point values on file; should be converted to radians upon usage
        this->tableThetaDegrees[d] = randomThetaDegrees[d];
    }

    uniform_int_distribution<int> uniformIntegerDistributionTablesA(- (1 << 10), 1 << 10);

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            // tableAx and tableAy contain X- and Y-coefficients of unnormalized 2D vectors to avoid storing floating-point values on file; should be normalized upon usage
            this->tableAx[k * this->alphabet.size() + sigma] = uniformIntegerDistributionTablesA(pseudoRandomNumberEngine);
            this->tableAy[k * this->alphabet.size() + sigma] = uniformIntegerDistributionTablesA(pseudoRandomNumberEngine);
        }
    }

    vector<int> paramDValues;

    for (int i = 0; i < this->paramD; i++)
    {
        paramDValues.push_back(i);
    }

    if (this->alphabet.size() > this->paramD)
    {
        for (int i = this->paramD; i < this->alphabet.size(); i++)
        {
            paramDValues.push_back(i % this->paramD);
        }
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();

        shuffle(paramDValues.begin(), paramDValues.end(), default_random_engine(currentTimeBasedSeed));

        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            this->tableB[k * this->alphabet.size() + sigma] = paramDValues[sigma];
        }
    }

    return;
}

void PotHash::solveDP(const string& sequence)
{
    return;
}