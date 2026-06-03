#include "pothash.hpp"

#include <algorithm>
#include <chrono>
#include <random>

Vector2D PotHash::projectVectorOnRay(const Vector2D& vector, const int rayAngleDegree)
{
    float projectionOfVectorOnRay = vector.vectorComponentX * cos(rayAngleDegree * PI / 180) + vector.vectorComponentY * sin(rayAngleDegree * PI / 180);

    return Vector2D{(float) (projectionOfVectorOnRay * cos(rayAngleDegree * PI / 180)), (float) (projectionOfVectorOnRay * sin(rayAngleDegree * PI / 180))};
}

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

    this->tableA = new Vector2D[subseqSizeK * alphabet.size()];
    this->tableB = new int[subseqSizeK * alphabet.size()];
}

PotHash::~PotHash()
{
    delete[] this->tableThetaDegrees;

    delete[] this->tableA;
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

    uniform_int_distribution<int> uniformIntegerDistributionTableA(- (1 << 10), 1 << 10);

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            // tableA contains integer X- and Y-coefficients of unnormalized 2D vectors to avoid storing floating-point values on file; should be normalized upon usage
            this->tableA[k * this->alphabet.size() + sigma] = Vector2D{(float) uniformIntegerDistributionTableA(pseudoRandomNumberEngine), (float) uniformIntegerDistributionTableA(pseudoRandomNumberEngine)};
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

Seed PotHash::solveDP(const string& sequence)
{
    int seqSizeN = sequence.length();

    DpTableCell *dpTableT = new DpTableCell[(seqSizeN + 1) * (this->subseqSizeK + 1) * this->paramD];

    for (int n = 0; n <= seqSizeN; n++)
    {
        // Initialization: no matter what the prefix sequence is, null vector for empty (k = 0) subsequence when d = 0 at the beginning
        // Note that Score(z1) = Project(Score(0) + A[1, z1], ("0" + B[1, z1]) mod D); "0" in the formulation is the reason why we have null vector for k = 0 and d = 0 at the beginning
        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + 0].dpTableCellScoreVector = Vector2D{0, 0};
        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + 0].dpTableCellSeed = 0;
    }

    for (int n = 0; n <= seqSizeN; n++)
    {
        for (int d = 1; d < this->paramD; d++)
        {
            // Initialization: no matter what the prefix sequence is, invalid vector for empty (k = 0) subsequence when d > 0 at the beginning
            dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + d].dpTableCellScoreVector = Vector2D{NEG_INF, NEG_INF};
            dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + d].dpTableCellSeed = UINT64_MAX;
        }
    }

    for (int k = 1; k <= this->subseqSizeK; k++)
    {
        for (int d = 0; d < this->paramD; d++)
        {
            // Initialization: when prefix sequence is empty, no valid subsequence can be formed unless the subsequence itself is empty
            dpTableT[0 * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = Vector2D{NEG_INF, NEG_INF};
            dpTableT[0 * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = UINT64_MAX;
        }
    }

    for (int n = 1; n <= seqSizeN; n++)
    {
        for (int k = 1; k <= this->subseqSizeK; k++)
        {
            for (int d = 0; d < this->paramD; d++)
            {
                // Case-1: do not add n-th base of sequence as k-th base of subsequence
                Vector2D case1ScoreVector = dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector;

                // Case-2: add n-th base of sequence as k-th base of subsequence
                int dPrevious = (d - this->tableB[(k - 1) * this->alphabet.size() + this->alphabet[sequence[n - 1]]] + this->paramD) % this->paramD;

                Vector2D case2ScoreVector = Vector2D{NEG_INF, NEG_INF};

                if (dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrevious].dpTableCellScoreVector.vectorComponentX > NEG_INF && dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrevious].dpTableCellScoreVector.vectorComponentY > NEG_INF)
                {
                    case2ScoreVector = this->projectVectorOnRay(dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrevious].dpTableCellScoreVector + this->tableA[(k - 1) * this->alphabet.size() + this->alphabet[sequence[n - 1]]].normalize(), this->tableThetaDegrees[d]);
                }

                // Update the DP table cell at location [n, k, d] with a score vector and a prefix of seed, based on the magnitude of score vectors from above two cases
                if ((case1ScoreVector.vectorComponentX > NEG_INF && case1ScoreVector.vectorComponentY > NEG_INF) && (case2ScoreVector.vectorComponentX > NEG_INF && case2ScoreVector.vectorComponentY > NEG_INF))
                {
                    if (case1ScoreVector.magnitude() > case2ScoreVector.magnitude())
                    {
                        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = case1ScoreVector;
                        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed;
                    }
                    else
                    {
                        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = case2ScoreVector;
                        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = (dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrevious].dpTableCellSeed << 2) | (uint64_t) this->alphabet[sequence[n - 1]];
                    }
                }
                else if ((case1ScoreVector.vectorComponentX > NEG_INF && case1ScoreVector.vectorComponentY > NEG_INF) && (case2ScoreVector.vectorComponentX == NEG_INF && case2ScoreVector.vectorComponentY == NEG_INF))
                {
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = case1ScoreVector;
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed;
                }
                else if ((case1ScoreVector.vectorComponentX == NEG_INF && case1ScoreVector.vectorComponentY == NEG_INF) && (case2ScoreVector.vectorComponentX > NEG_INF && case2ScoreVector.vectorComponentY > NEG_INF))
                {
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = case2ScoreVector;
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = (dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrevious].dpTableCellSeed << 2) | (uint64_t) this->alphabet[sequence[n - 1]];
                }
                else
                {
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = Vector2D{NEG_INF, NEG_INF};
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = UINT64_MAX;
                }
            }
        }
    }

    Seed seed;

    seed.seedParamD = -1;

    for (int d = 0; d < this->paramD; d++)
    {
        if (dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentX > NEG_INF && dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentY > NEG_INF)
        {
            seed.seedParamD = d;

            break;
        }
    }

    seed.seedScoreVector = (seed.seedParamD == -1) ? Vector2D{NEG_INF, NEG_INF} : dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + seed.seedParamD].dpTableCellScoreVector;

    delete[] dpTableT;

    return seed;
}