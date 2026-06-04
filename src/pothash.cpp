#include "pothash.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>

Vector2D PotHash::projectVectorOnRay(const Vector2D &vector, const int rayAngleDegree)
{
    float projectionOfVectorOnRay = vector.vectorComponentX * cos(rayAngleDegree * PI / 180) + vector.vectorComponentY * sin(rayAngleDegree * PI / 180);

    return Vector2D{(float) (projectionOfVectorOnRay * cos(rayAngleDegree * PI / 180)), (float) (projectionOfVectorOnRay * sin(rayAngleDegree * PI / 180))};
}

string PotHash::generateSeedSubsequence(const uint64_t &seedSubsequence)
{
    string seed = "";
    
    uint64_t seedSubsequenceCopy = seedSubsequence;

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        seed = this->reverseAlphabet[seedSubsequenceCopy & 3] + seed;

        seedSubsequenceCopy = seedSubsequenceCopy >> 2;
    }

    return seed;
}

PotHash::PotHash() : PotHash(16, 16, {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}})
{
    // Calling parameterized constructor from default constructor with default arguments
}

PotHash::PotHash(const int paramD, const int subseqSizeK, const map<char, int> &alphabet)
{
    this->paramD = (paramD % 2 == 0) ? paramD : paramD + 1;
    this->subseqSizeK = subseqSizeK;
    this->alphabet = alphabet;

    this->tableThetaDegrees = new int[paramD];

    this->tableA = new Vector2D[subseqSizeK * alphabet.size()];
    this->tableB = new int[subseqSizeK * alphabet.size()];

    for (auto const &pair: alphabet)
    {
        this->reverseAlphabet[pair.second] = pair.first;
    }
}

PotHash::~PotHash()
{
    delete[] this->tableThetaDegrees;

    delete[] this->tableA;
    delete[] this->tableB;
}

void PotHash::generateTables(const bool shouldSaveTables, const int tablesFileVersion)
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

    if (shouldSaveTables)
    {
        string tablesFileFolderPath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(this->alphabet.size());

        filesystem::path tablesFileFolder(tablesFileFolderPath);

        if (!filesystem::exists(tablesFileFolder))
        {
            if (!filesystem::create_directories(tablesFileFolder))
            {
                exit(EXIT_FAILURE);
            }
        }

        string tablesFilePath = tablesFileFolderPath + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(this->alphabet.size()) + string("_Version") + to_string(tablesFileVersion);

        ofstream tablesFile(tablesFilePath);

        if (!tablesFile.is_open())
        {
            exit(EXIT_FAILURE);
        }

        for (int d = 0; d < this->paramD; d++)
        {
            tablesFile << this->tableThetaDegrees[d] << " ";
        }

        tablesFile << "\n" << endl;

        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
            {
                tablesFile << (int) this->tableA[k * this->alphabet.size() + sigma].vectorComponentX << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;

        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
            {
                tablesFile << (int) this->tableA[k * this->alphabet.size() + sigma].vectorComponentY << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;

        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
            {
                tablesFile << this->tableB[k * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile.close();
    }

    return;
}

void PotHash::loadTables(const int tablesFileVersion)
{
    string tablesFilePath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(this->alphabet.size()) + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(this->alphabet.size()) + string("_Version") + to_string(tablesFileVersion);

    ifstream tablesFile(tablesFilePath);

    if (!tablesFile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    for (int d = 0; d < this->paramD; d++)
    {
        tablesFile >> this->tableThetaDegrees[d];
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            tablesFile >> this->tableA[k * this->alphabet.size() + sigma].vectorComponentX;
        }
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            tablesFile >> this->tableA[k * this->alphabet.size() + sigma].vectorComponentY;
        }
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            tablesFile >> this->tableB[k * this->alphabet.size() + sigma];
        }
    }

    tablesFile.close();

    return;
}

vector<Seed> PotHash::solveDP(const string &sequence, const uint8_t &seedGenerationMode)
{
    int seqSizeN = sequence.length();

    DpTableCell *dpTableT = new DpTableCell[(seqSizeN + 1) * (this->subseqSizeK + 1) * this->paramD];

    for (int n = 0; n <= seqSizeN; n++)
    {
        // Initialization: no matter what the prefix sequence is, null vector for empty (k = 0) subsequence when d = 0 at the beginning
        // Note that Score(z1) = Project(Score(0) + A[1, z1], ("0" + B[1, z1]) mod D); "0" in the formulation is the reason why we have null vector for k = 0 and d = 0 at the beginning
        dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + 0] = DpTableCell{Vector2D{0, 0}, 0};
    }

    for (int n = 0; n <= seqSizeN; n++)
    {
        for (int d = 1; d < this->paramD; d++)
        {
            // Initialization: no matter what the prefix sequence is, invalid vector for empty (k = 0) subsequence when d > 0 at the beginning
            dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + d] = DpTableCell{Vector2D{NEG_INF, NEG_INF}, UINT64_MAX};
        }
    }

    for (int k = 1; k <= this->subseqSizeK; k++)
    {
        for (int d = 0; d < this->paramD; d++)
        {
            // Initialization: when prefix sequence is empty, no valid subsequence can be formed unless the subsequence itself is empty
            dpTableT[0 * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d] = DpTableCell{Vector2D{NEG_INF, NEG_INF}, UINT64_MAX};
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
                bool isCase1Better = false, isCase2Better = false;

                if ((case1ScoreVector.vectorComponentX > NEG_INF && case1ScoreVector.vectorComponentY > NEG_INF) && (case2ScoreVector.vectorComponentX > NEG_INF && case2ScoreVector.vectorComponentY > NEG_INF))
                {
                    if (case1ScoreVector.magnitude() > case2ScoreVector.magnitude())
                    {
                        isCase1Better = true;
                    }
                    else
                    {
                        isCase2Better = true;
                    }
                }
                else if (case1ScoreVector.vectorComponentX > NEG_INF && case1ScoreVector.vectorComponentY > NEG_INF)
                {
                    isCase1Better = true;
                }
                else if (case2ScoreVector.vectorComponentX > NEG_INF && case2ScoreVector.vectorComponentY > NEG_INF)
                {
                    isCase2Better = true;
                }

                if (isCase1Better)
                {
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = case1ScoreVector;
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed;
                }
                else if (isCase2Better)
                {
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellScoreVector = case2ScoreVector;
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d].dpTableCellSeed = (dpTableT[(n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrevious].dpTableCellSeed << 2) | (uint64_t) this->alphabet[sequence[n - 1]];
                }
                else
                {
                    dpTableT[n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d] = DpTableCell{Vector2D{NEG_INF, NEG_INF}, UINT64_MAX};
                }
            }
        }
    }

    vector<Seed> seeds;

    if (seedGenerationMode == 0)
    {
        // Extract a single valid seed from the dpTableT[N, K] column having the lowest d value associated with it
        Seed seed = Seed{-1, Vector2D{NEG_INF, NEG_INF}, "X"};

        for (int d = 0; d < this->paramD; d++)
        {
            if (dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentX > NEG_INF && dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentY > NEG_INF)
            {
                seed.seedParamD = d;

                break;
            }
        }

        if (seed.seedParamD != -1)
        {
            seed.seedScoreVector = dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + seed.seedParamD].dpTableCellScoreVector;

            seed.seedSubsequence = this->generateSeedSubsequence(dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + seed.seedParamD].dpTableCellSeed);
        }

        seeds.push_back(seed);
    }
    else if (seedGenerationMode == 1)
    {
        // Extract a single valid seed from the dpTableT[N, K] column having the highest magnitude of score vector associated with it
        Seed seed = Seed{-1, Vector2D{NEG_INF, NEG_INF}, "X"};

        for (int d = 0; d < this->paramD; d++)
        {
            if (dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentX > NEG_INF && dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentY > NEG_INF)
            {
                if (seed.seedScoreVector.vectorComponentX > NEG_INF && seed.seedScoreVector.vectorComponentY > NEG_INF)
                {
                    if (dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.magnitude() > seed.seedScoreVector.magnitude())
                    {
                        seed.seedParamD = d;
                        seed.seedScoreVector = dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector;
                    }
                }
                else
                {
                    seed.seedParamD = d;
                    seed.seedScoreVector = dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector;
                }
            }
        }

        if (seed.seedParamD != -1)
        {
            seed.seedSubsequence = this->generateSeedSubsequence(dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + seed.seedParamD].dpTableCellSeed);
        }

        seeds.push_back(seed);
    }
    else if (seedGenerationMode == 2)
    {
        // Extract all D seeds from the dpTableT[N, K] column regardless of them being valid or invalid
        for (int d = 0; d < this->paramD; d++)
        {
            Seed seed = Seed{-1, Vector2D{NEG_INF, NEG_INF}, "X"};

            if (dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentX > NEG_INF && dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector.vectorComponentY > NEG_INF)
            {
                seed.seedParamD = d;
                seed.seedScoreVector = dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScoreVector;
                seed.seedSubsequence = this->generateSeedSubsequence(dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellSeed);
            }

            seeds.push_back(seed);
        }
    }
    else
    {
        exit(EXIT_FAILURE);
    }

    delete[] dpTableT;

    return seeds;
}