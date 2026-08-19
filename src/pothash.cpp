#include "pothash.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>

float PotHash::projectScalar(const float r, const float thetaDegrees, const int rayIndex) const
{
    float thetaRadian = thetaDegrees * PI / 180.0;
    float rayAngleRadian = this->tableThetaDegrees[rayIndex] * PI / 180.0;

    return r * cos(thetaRadian - rayAngleRadian);
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

    this->tableThetaDegrees = new int[this->paramD];
    this->tableA = new float[this->subseqSizeK * alphabet.size()];

    for (auto const &pair: alphabet)
    {
        this->reverseAlphabet[pair.second] = pair.first;
    }
}

PotHash::~PotHash()
{
    delete[] this->tableThetaDegrees;
    delete[] this->tableA;
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
        // tableThetaDegrees contains random ray angles in degrees; converted to radians upon usage
        this->tableThetaDegrees[d] = randomThetaDegrees[d];
    }

    // tableA stores angles in degrees in [0, 360), matching PDF angles in [0, 2π)
    uniform_real_distribution<float> uniformRealDistributionTableA(0.0f, 360.0f);

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < this->alphabet.size(); sigma++)
        {
            this->tableA[k * this->alphabet.size() + sigma] = uniformRealDistributionTableA(pseudoRandomNumberEngine);
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
                tablesFile << this->tableA[k * this->alphabet.size() + sigma] << " ";
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
            tablesFile >> this->tableA[k * this->alphabet.size() + sigma];
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
        for (int d = 0; d < this->paramD; d++)
        {
            // f(∅, d) = 0 for every d (August 4 formulation; no B/ψ constraint on the empty subsequence)
            dpTableT[n * (this->subseqSizeK + 1) * this->paramD + 0 * this->paramD + d] = DpTableCell{0.0f, 0};
        }
    }

    for (int k = 1; k <= this->subseqSizeK; k++)
    {
        for (int d = 0; d < this->paramD; d++)
        {
            // Initialization: when prefix sequence is empty, no valid subsequence of length k > 0 can be formed
            dpTableT[0 * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d] = DpTableCell{NEG_INF, UINT64_MAX};
        }
    }

    for (int n = 1; n <= seqSizeN; n++)
    {
        for (int k = 1; k <= this->subseqSizeK; k++)
        {
            for (int d = 0; d < this->paramD; d++)
            {
                int currentIndex = n * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d;
                int skipIndex = (n - 1) * (this->subseqSizeK + 1) * this->paramD + k * this->paramD + d;
                int currentBase = this->alphabet[sequence[n - 1]];

                // Case-1: do not add n-th base of sequence as k-th base of subsequence
                float case1Score = dpTableT[skipIndex].dpTableCellScore;
                uint64_t case1Seed = dpTableT[skipIndex].dpTableCellSeed;

                // Case-2: add n-th base; maximize over all previous ray indices d'
                float maxPrevProjection = NEG_INF;
                uint64_t bestPrevSeed = UINT64_MAX;
                bool hasValidPrev = false;

                for (int dPrev = 0; dPrev < this->paramD; dPrev++)
                {
                    int takePreviousIndex = (n - 1) * (this->subseqSizeK + 1) * this->paramD + (k - 1) * this->paramD + dPrev;
                    float prevScore = dpTableT[takePreviousIndex].dpTableCellScore;

                    if (prevScore > NEG_INF)
                    {
                        float prevProjection = this->projectScalar(prevScore, (float) this->tableThetaDegrees[dPrev], d);

                        if (!hasValidPrev || prevProjection > maxPrevProjection)
                        {
                            maxPrevProjection = prevProjection;
                            bestPrevSeed = dpTableT[takePreviousIndex].dpTableCellSeed;
                            hasValidPrev = true;
                        }
                    }
                }

                float case2Score = NEG_INF;
                uint64_t case2Seed = UINT64_MAX;

                if (hasValidPrev)
                {
                    float angleA = this->tableA[(k - 1) * this->alphabet.size() + currentBase];
                    float addProjection = this->projectScalar(1.0f, angleA, d);
                    float clampedPrev = (maxPrevProjection > 0.0f) ? maxPrevProjection : 0.0f;
                    float combined = addProjection + clampedPrev;

                    case2Score = (combined > 0.0f) ? combined : 0.0f;
                    case2Seed = (bestPrevSeed << 2) | (uint64_t) currentBase;
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

    vector<Seed> seeds;

    if (seedGenerationMode == 0)
    {
        // Extract a single valid seed from the dpTableT[N, K] column having the lowest d value associated with it
        Seed seed = Seed{-1, NEG_INF, "X"};

        for (int d = 0; d < this->paramD; d++)
        {
            if (dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScore > NEG_INF)
            {
                seed.seedParamD = d;

                break;
            }
        }

        if (seed.seedParamD != -1)
        {
            int terminalIndex = seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + seed.seedParamD;

            seed.seedScore = dpTableT[terminalIndex].dpTableCellScore;
            seed.seedSubsequence = this->generateSeedSubsequence(dpTableT[terminalIndex].dpTableCellSeed);
        }

        seeds.push_back(seed);
    }
    else if (seedGenerationMode == 1)
    {
        // Extract a single valid seed with maximum scalar score (PDF: final score = max_d f(z, d))
        Seed seed = Seed{-1, NEG_INF, "X"};

        for (int d = 0; d < this->paramD; d++)
        {
            float terminalScore = dpTableT[seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d].dpTableCellScore;

            if (terminalScore > NEG_INF)
            {
                if (seed.seedParamD == -1 || terminalScore > seed.seedScore)
                {
                    seed.seedParamD = d;
                    seed.seedScore = terminalScore;
                }
            }
        }

        if (seed.seedParamD != -1)
        {
            int terminalIndex = seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + seed.seedParamD;

            seed.seedSubsequence = this->generateSeedSubsequence(dpTableT[terminalIndex].dpTableCellSeed);
        }

        seeds.push_back(seed);
    }
    else if (seedGenerationMode == 2)
    {
        // Extract all D seeds from the dpTableT[N, K] column regardless of them being valid or invalid
        for (int d = 0; d < this->paramD; d++)
        {
            Seed seed = Seed{-1, NEG_INF, "X"};
            int terminalIndex = seqSizeN * (this->subseqSizeK + 1) * this->paramD + this->subseqSizeK * this->paramD + d;

            if (dpTableT[terminalIndex].dpTableCellScore > NEG_INF)
            {
                seed.seedParamD = d;
                seed.seedScore = dpTableT[terminalIndex].dpTableCellScore;
                seed.seedSubsequence = this->generateSeedSubsequence(dpTableT[terminalIndex].dpTableCellSeed);
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
