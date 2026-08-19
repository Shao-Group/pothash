#include "pothash.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>

float PotHash::projectScalar(const float r, const float thetaDegrees, const float thetaPrimeDegrees) const
{
    float thetaRadian = thetaDegrees * PI / 180.0;
    float thetaPrimeRadian = thetaPrimeDegrees * PI / 180.0;

    return r * cos(thetaRadian - thetaPrimeRadian);
}

int PotHash::indexA(const int k, const int sigma) const
{
    return k * (int) this->alphabet.size() + sigma;
}

int PotHash::indexB(const int k, const int sigma, const int d) const
{
    return k * (int) this->alphabet.size() * this->paramD + sigma * this->paramD + d;
}

int PotHash::indexDP(const int n, const int k, const int d, const int lastBase) const
{
    const int alphabetSize = (int) this->alphabet.size();

    return n * (this->subseqSizeK + 1) * this->paramD * alphabetSize + k * this->paramD * alphabetSize + d * alphabetSize + lastBase;
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
    // August 16 formulation: D need not be even (no global complementary Theta set)
    this->paramD = paramD;
    this->subseqSizeK = subseqSizeK;
    this->alphabet = alphabet;

    const int alphabetSize = (int) alphabet.size();

    this->tableA = new float[this->subseqSizeK * alphabetSize];
    this->tableB = new float[this->subseqSizeK * alphabetSize * this->paramD];

    for (auto const &pair: alphabet)
    {
        this->reverseAlphabet[pair.second] = pair.first;
    }
}

PotHash::~PotHash()
{
    delete[] this->tableA;
    delete[] this->tableB;
}

void PotHash::generateTables(const bool shouldSaveTables, const int tablesFileVersion)
{
    random_device seedSource;

    mt19937 pseudoRandomNumberEngine(seedSource());

    const int alphabetSize = (int) this->alphabet.size();

    // ThetaA and ThetaB store angles in degrees in [0, 360), matching PDF [0, 2π)
    uniform_real_distribution<float> uniformAngleDegrees(0.0f, 360.0f);

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            this->tableA[this->indexA(k, sigma)] = uniformAngleDegrees(pseudoRandomNumberEngine);

            for (int d = 0; d < this->paramD; d++)
            {
                this->tableB[this->indexB(k, sigma, d)] = uniformAngleDegrees(pseudoRandomNumberEngine);
            }
        }
    }

    if (shouldSaveTables)
    {
        string tablesFileFolderPath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize);

        filesystem::path tablesFileFolder(tablesFileFolderPath);

        if (!filesystem::exists(tablesFileFolder))
        {
            if (!filesystem::create_directories(tablesFileFolder))
            {
                exit(EXIT_FAILURE);
            }
        }

        string tablesFilePath = tablesFileFolderPath + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize) + string("_Version") + to_string(tablesFileVersion);

        ofstream tablesFile(tablesFilePath);

        if (!tablesFile.is_open())
        {
            exit(EXIT_FAILURE);
        }

        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                tablesFile << this->tableA[this->indexA(k, sigma)] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;

        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                for (int d = 0; d < this->paramD; d++)
                {
                    tablesFile << this->tableB[this->indexB(k, sigma, d)] << " ";
                }

                tablesFile << endl;
            }
        }

        tablesFile.close();
    }

    return;
}

void PotHash::loadTables(const int tablesFileVersion)
{
    const int alphabetSize = (int) this->alphabet.size();

    string tablesFilePath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize) + filesystem::path::preferred_separator + string("tables_D") + to_string(this->paramD) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize) + string("_Version") + to_string(tablesFileVersion);

    ifstream tablesFile(tablesFilePath);

    if (!tablesFile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            tablesFile >> this->tableA[this->indexA(k, sigma)];
        }
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            for (int d = 0; d < this->paramD; d++)
            {
                tablesFile >> this->tableB[this->indexB(k, sigma, d)];
            }
        }
    }

    tablesFile.close();

    return;
}

vector<Seed> PotHash::solveDP(const string &sequence, const uint8_t &seedGenerationMode)
{
    const int alphabetSize = (int) this->alphabet.size();
    const int seqSizeN = (int) sequence.length();
    const int dpPlane = (this->subseqSizeK + 1) * this->paramD * alphabetSize;

    // T[n,k,d,σ]: best score for length-k subsequences of s[1..n] with index d and last base σ
    DpTableCell *dpTableT = new DpTableCell[(seqSizeN + 1) * dpPlane];

    // Empty (k=0): score 0 for every d; lastBase unused — fill all σ slots with 0 for simple skip/take wiring
    for (int n = 0; n <= seqSizeN; n++)
    {
        for (int d = 0; d < this->paramD; d++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                dpTableT[this->indexDP(n, 0, d, sigma)] = DpTableCell{0.0f, 0};
            }
        }
    }

    // n=0 cannot form positive-length subsequences
    for (int k = 1; k <= this->subseqSizeK; k++)
    {
        for (int d = 0; d < this->paramD; d++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                dpTableT[this->indexDP(0, k, d, sigma)] = DpTableCell{NEG_INF, UINT64_MAX};
            }
        }
    }

    for (int n = 1; n <= seqSizeN; n++)
    {
        for (int k = 1; k <= this->subseqSizeK; k++)
        {
            for (int d = 0; d < this->paramD; d++)
            {
                int currentBase = this->alphabet[sequence[n - 1]];

                float angleA = this->tableA[this->indexA(k - 1, currentBase)];
                float angleBCurrent = this->tableB[this->indexB(k - 1, currentBase, d)];
                float addProjection = this->projectScalar(1.0f, angleA, angleBCurrent);

                // Precompute best previous projection over (d', σ_prev) when taking sn (k >= 2)
                float maxPrevProjection = NEG_INF;
                uint64_t bestPrevSeed = UINT64_MAX;
                bool hasValidPrev = false;

                if (k >= 2)
                {
                    for (int dPrev = 0; dPrev < this->paramD; dPrev++)
                    {
                        for (int sigmaPrev = 0; sigmaPrev < alphabetSize; sigmaPrev++)
                        {
                            int takePreviousIndex = this->indexDP(n - 1, k - 1, dPrev, sigmaPrev);
                            float prevScore = dpTableT[takePreviousIndex].dpTableCellScore;

                            if (prevScore > NEG_INF)
                            {
                                float angleBPrev = this->tableB[this->indexB(k - 2, sigmaPrev, dPrev)];
                                float prevProjection = this->projectScalar(prevScore, angleBPrev, angleBCurrent);

                                if (!hasValidPrev || prevProjection > maxPrevProjection)
                                {
                                    maxPrevProjection = prevProjection;
                                    bestPrevSeed = dpTableT[takePreviousIndex].dpTableCellSeed;
                                    hasValidPrev = true;
                                }
                            }
                        }
                    }
                }

                float takeScore = NEG_INF;
                uint64_t takeSeed = UINT64_MAX;

                if (k == 1)
                {
                    takeScore = (addProjection > 0.0f) ? addProjection : 0.0f;
                    takeSeed = (uint64_t) currentBase;
                }
                else if (hasValidPrev)
                {
                    float clampedPrev = (maxPrevProjection > 0.0f) ? maxPrevProjection : 0.0f;
                    float combined = addProjection + clampedPrev;

                    takeScore = (combined > 0.0f) ? combined : 0.0f;
                    takeSeed = (bestPrevSeed << 2) | (uint64_t) currentBase;
                }

                for (int sigma = 0; sigma < alphabetSize; sigma++)
                {
                    int currentIndex = this->indexDP(n, k, d, sigma);
                    int skipIndex = this->indexDP(n - 1, k, d, sigma);

                    // Case-1: sn not picked — keep same last base σ
                    float case1Score = dpTableT[skipIndex].dpTableCellScore;
                    uint64_t case1Seed = dpTableT[skipIndex].dpTableCellSeed;

                    // Case-2: sn picked — only valid when new last base equals sn
                    float case2Score = NEG_INF;
                    uint64_t case2Seed = UINT64_MAX;

                    if (sigma == currentBase)
                    {
                        case2Score = takeScore;
                        case2Seed = takeSeed;
                    }

                    bool isCase1Better = false, isCase2Better = false;

                    if (case1Score > NEG_INF && case2Score > NEG_INF)
                    {
                        if (case1Score > case2Score)
                        {
                            isCase1Better = true;
                        }
                        else
                        {
                            isCase2Better = true;
                        }
                    }
                    else if (case1Score > NEG_INF)
                    {
                        isCase1Better = true;
                    }
                    else if (case2Score > NEG_INF)
                    {
                        isCase2Better = true;
                    }

                    if (isCase1Better)
                    {
                        dpTableT[currentIndex] = DpTableCell{case1Score, case1Seed};
                    }
                    else if (isCase2Better)
                    {
                        dpTableT[currentIndex] = DpTableCell{case2Score, case2Seed};
                    }
                    else
                    {
                        dpTableT[currentIndex] = DpTableCell{NEG_INF, UINT64_MAX};
                    }
                }
            }
        }
    }

    auto bestOverLastBase = [&](const int d, float &outScore, uint64_t &outSeed) -> bool
    {
        outScore = NEG_INF;
        outSeed = UINT64_MAX;
        bool found = false;

        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            int terminalIndex = this->indexDP(seqSizeN, this->subseqSizeK, d, sigma);
            float score = dpTableT[terminalIndex].dpTableCellScore;

            if (score > NEG_INF && (!found || score > outScore))
            {
                outScore = score;
                outSeed = dpTableT[terminalIndex].dpTableCellSeed;
                found = true;
            }
        }

        return found;
    };

    vector<Seed> seeds;

    if (seedGenerationMode == 0)
    {
        Seed seed = Seed{-1, NEG_INF, "X"};

        for (int d = 0; d < this->paramD; d++)
        {
            float score;
            uint64_t packedSeed;

            if (bestOverLastBase(d, score, packedSeed))
            {
                seed.seedParamD = d;
                seed.seedScore = score;
                seed.seedSubsequence = this->generateSeedSubsequence(packedSeed);

                break;
            }
        }

        seeds.push_back(seed);
    }
    else if (seedGenerationMode == 1)
    {
        Seed seed = Seed{-1, NEG_INF, "X"};

        for (int d = 0; d < this->paramD; d++)
        {
            float score;
            uint64_t packedSeed;

            if (bestOverLastBase(d, score, packedSeed))
            {
                if (seed.seedParamD == -1 || score > seed.seedScore)
                {
                    seed.seedParamD = d;
                    seed.seedScore = score;
                    seed.seedSubsequence = this->generateSeedSubsequence(packedSeed);
                }
            }
        }

        seeds.push_back(seed);
    }
    else if (seedGenerationMode == 2)
    {
        for (int d = 0; d < this->paramD; d++)
        {
            Seed seed = Seed{-1, NEG_INF, "X"};
            float score;
            uint64_t packedSeed;

            if (bestOverLastBase(d, score, packedSeed))
            {
                seed.seedParamD = d;
                seed.seedScore = score;
                seed.seedSubsequence = this->generateSeedSubsequence(packedSeed);
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
