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

int PotHash::indexA(const int k, const int sigma, const int db) const
{
    return k * (int) this->alphabet.size() * this->paramDb + sigma * this->paramDb + db;
}

int PotHash::indexB(const int k, const int sigma, const int da, const int db) const
{
    return k * (int) this->alphabet.size() * this->paramDa * this->paramDb + sigma * this->paramDa * this->paramDb + da * this->paramDb + db;
}

int PotHash::indexC(const int k, const int sigma) const
{
    return k * (int) this->alphabet.size() + sigma;
}

int PotHash::indexDP(const int n, const int k, const int da, const int db, const int lastBase) const
{
    const int alphabetSize = (int) this->alphabet.size();

    return n * (this->subseqSizeK + 1) * this->paramDa * this->paramDb * alphabetSize
         + k * this->paramDa * this->paramDb * alphabetSize
         + da * this->paramDb * alphabetSize
         + db * alphabetSize
         + lastBase;
}

string PotHash::generateSeedSubsequence(const uint64_t &seedSubsequence) const
{
    string seed = "";
    uint64_t seedSubsequenceCopy = seedSubsequence;

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        seed = this->reverseAlphabet.at(seedSubsequenceCopy & 3) + seed;
        seedSubsequenceCopy = seedSubsequenceCopy >> 2;
    }

    return seed;
}

PotHash::PotHash() : PotHash(17, 17, 16, {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}})
{
}

PotHash::PotHash(const int paramDa, const int paramDb, const int subseqSizeK, const map<char, int> &alphabet)
{
    this->paramDa = paramDa;
    this->paramDb = paramDb;
    this->subseqSizeK = subseqSizeK;
    this->alphabet = alphabet;

    const int alphabetSize = (int) alphabet.size();

    this->tableA = new float[this->subseqSizeK * alphabetSize * this->paramDb];
    this->tableB = new float[this->subseqSizeK * alphabetSize * this->paramDa * this->paramDb];
    this->tableC = new int[this->subseqSizeK * alphabetSize];

    for (auto const &pair: alphabet)
    {
        this->reverseAlphabet[pair.second] = pair.first;
    }
}

PotHash::~PotHash()
{
    delete[] this->tableA;
    delete[] this->tableB;
    delete[] this->tableC;
}

void PotHash::generateTables(const bool shouldSaveTables, const int tablesFileVersion)
{
    random_device seedSource;
    mt19937 pseudoRandomNumberEngine(seedSource());

    const int alphabetSize = (int) this->alphabet.size();
    uniform_real_distribution<float> uniformAngleDegrees(0.0f, 360.0f);
    uniform_int_distribution<int> uniformC(0, this->paramDb - 1);

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            this->tableC[this->indexC(k, sigma)] = uniformC(pseudoRandomNumberEngine);

            for (int db = 0; db < this->paramDb; db++)
            {
                this->tableA[this->indexA(k, sigma, db)] = uniformAngleDegrees(pseudoRandomNumberEngine);

                for (int da = 0; da < this->paramDa; da++)
                {
                    this->tableB[this->indexB(k, sigma, da, db)] = uniformAngleDegrees(pseudoRandomNumberEngine);
                }
            }
        }
    }

    if (shouldSaveTables)
    {
        string tablesFileFolderPath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_Da") + to_string(this->paramDa) + string("_Db") + to_string(this->paramDb) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize);

        filesystem::path tablesFileFolder(tablesFileFolderPath);

        if (!filesystem::exists(tablesFileFolder))
        {
            if (!filesystem::create_directories(tablesFileFolder))
            {
                exit(EXIT_FAILURE);
            }
        }

        string tablesFilePath = tablesFileFolderPath + filesystem::path::preferred_separator + string("tables_Da") + to_string(this->paramDa) + string("_Db") + to_string(this->paramDb) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize) + string("_Version") + to_string(tablesFileVersion);

        ofstream tablesFile(tablesFilePath);

        if (!tablesFile.is_open())
        {
            exit(EXIT_FAILURE);
        }

        // C: K lines of Sigma ints
        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                tablesFile << this->tableC[this->indexC(k, sigma)] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;

        // A: for each k, sigma: Db angles
        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                for (int db = 0; db < this->paramDb; db++)
                {
                    tablesFile << this->tableA[this->indexA(k, sigma, db)] << " ";
                }

                tablesFile << endl;
            }
        }

        tablesFile << endl;

        // B: for each k, sigma, da: Db angles
        for (int k = 0; k < this->subseqSizeK; k++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                for (int da = 0; da < this->paramDa; da++)
                {
                    for (int db = 0; db < this->paramDb; db++)
                    {
                        tablesFile << this->tableB[this->indexB(k, sigma, da, db)] << " ";
                    }

                    tablesFile << endl;
                }
            }
        }

        tablesFile.close();
    }

    return;
}

void PotHash::loadTables(const int tablesFileVersion)
{
    const int alphabetSize = (int) this->alphabet.size();

    string tablesFilePath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_Da") + to_string(this->paramDa) + string("_Db") + to_string(this->paramDb) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize) + filesystem::path::preferred_separator + string("tables_Da") + to_string(this->paramDa) + string("_Db") + to_string(this->paramDb) + string("_K") + to_string(this->subseqSizeK) + string("_Sigma") + to_string(alphabetSize) + string("_Version") + to_string(tablesFileVersion);

    ifstream tablesFile(tablesFilePath);

    if (!tablesFile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            tablesFile >> this->tableC[this->indexC(k, sigma)];
        }
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            for (int db = 0; db < this->paramDb; db++)
            {
                tablesFile >> this->tableA[this->indexA(k, sigma, db)];
            }
        }
    }

    for (int k = 0; k < this->subseqSizeK; k++)
    {
        for (int sigma = 0; sigma < alphabetSize; sigma++)
        {
            for (int da = 0; da < this->paramDa; da++)
            {
                for (int db = 0; db < this->paramDb; db++)
                {
                    tablesFile >> this->tableB[this->indexB(k, sigma, da, db)];
                }
            }
        }
    }

    tablesFile.close();

    return;
}

Seed PotHash::solveDP(const string &sequence)
{
    const int alphabetSize = (int) this->alphabet.size();
    const int seqSizeN = (int) sequence.length();
    const int dpPlane = (this->subseqSizeK + 1) * this->paramDa * this->paramDb * alphabetSize;

    // T[n,k,da,db,σ]: best ω_da among length-k subsequences of s[1..n] with psi=db and last base σ
    DpTableCell *dpTableT = new DpTableCell[(seqSizeN + 1) * dpPlane];

    // Empty: psi=0, omega=0 for every da; lastBase unused — fill all σ slots
    for (int n = 0; n <= seqSizeN; n++)
    {
        for (int da = 0; da < this->paramDa; da++)
        {
            for (int db = 0; db < this->paramDb; db++)
            {
                for (int sigma = 0; sigma < alphabetSize; sigma++)
                {
                    int cellIdx = this->indexDP(n, 0, da, db, sigma);

                    if (db == 0)
                    {
                        dpTableT[cellIdx] = DpTableCell{0.0f, 0};
                    }
                    else
                    {
                        dpTableT[cellIdx] = DpTableCell{NEG_INF, UINT64_MAX};
                    }
                }
            }
        }
    }

    // n=0 cannot form positive-length subsequences
    for (int k = 1; k <= this->subseqSizeK; k++)
    {
        for (int da = 0; da < this->paramDa; da++)
        {
            for (int db = 0; db < this->paramDb; db++)
            {
                for (int sigma = 0; sigma < alphabetSize; sigma++)
                {
                    dpTableT[this->indexDP(0, k, da, db, sigma)] = DpTableCell{NEG_INF, UINT64_MAX};
                }
            }
        }
    }

    for (int n = 1; n <= seqSizeN; n++)
    {
        for (int k = 1; k <= this->subseqSizeK; k++)
        {
            for (int da = 0; da < this->paramDa; da++)
            {
                for (int db = 0; db < this->paramDb; db++)
                {
                    int currentBase = this->alphabet[sequence[n - 1]];
                    int cVal = this->tableC[this->indexC(k - 1, currentBase)];
                    int dbPrev = (db - cVal + this->paramDb) % this->paramDb;

                    float angleA = this->tableA[this->indexA(k - 1, currentBase, db)];
                    float angleBCurrent = this->tableB[this->indexB(k - 1, currentBase, da, db)];
                    float addProjection = this->projectScalar(1.0f, angleA, angleBCurrent);

                    float maxPrevProjection = NEG_INF;
                    uint64_t bestPrevSeed = UINT64_MAX;
                    bool hasValidPrev = false;

                    if (k == 1)
                    {
                        // Empty previous requires dbPrev == 0; projection of 0 is 0
                        if (dbPrev == 0)
                        {
                            maxPrevProjection = 0.0f;
                            bestPrevSeed = 0;
                            hasValidPrev = true;
                        }
                    }
                    else
                    {
                        for (int daPrev = 0; daPrev < this->paramDa; daPrev++)
                        {
                            for (int sigmaPrev = 0; sigmaPrev < alphabetSize; sigmaPrev++)
                            {
                                int takePreviousIndex = this->indexDP(n - 1, k - 1, daPrev, dbPrev, sigmaPrev);
                                float prevScore = dpTableT[takePreviousIndex].dpTableCellScore;

                                if (prevScore > NEG_INF)
                                {
                                    float angleBPrev = this->tableB[this->indexB(k - 2, sigmaPrev, daPrev, dbPrev)];
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

                    if (hasValidPrev)
                    {
                        float clampedPrev = (maxPrevProjection > 0.0f) ? maxPrevProjection : 0.0f;
                        float combined = addProjection + clampedPrev;

                        takeScore = (combined > 0.0f) ? combined : 0.0f;
                        takeSeed = (bestPrevSeed << 2) | (uint64_t) currentBase;
                    }

                    for (int sigma = 0; sigma < alphabetSize; sigma++)
                    {
                        int currentIndex = this->indexDP(n, k, da, db, sigma);
                        int skipIndex = this->indexDP(n - 1, k, da, db, sigma);

                        float case1Score = dpTableT[skipIndex].dpTableCellScore;
                        uint64_t case1Seed = dpTableT[skipIndex].dpTableCellSeed;

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
    }

    // Final: omega_db = max_{da,σ} T[N,K,da,db,σ]; pick smallest valid psi=db
    Seed seed = Seed{-1, NEG_INF, "X"};

    for (int db = 0; db < this->paramDb; db++)
    {
        float bestOmegaForDb = NEG_INF;
        int bestDaForDb = -1;
        uint64_t bestSeedForDb = UINT64_MAX;

        for (int da = 0; da < this->paramDa; da++)
        {
            for (int sigma = 0; sigma < alphabetSize; sigma++)
            {
                int terminalIndex = this->indexDP(seqSizeN, this->subseqSizeK, da, db, sigma);
                float score = dpTableT[terminalIndex].dpTableCellScore;

                if (score > NEG_INF && (bestDaForDb == -1 || score > bestOmegaForDb))
                {
                    bestOmegaForDb = score;
                    bestDaForDb = da;
                    bestSeedForDb = dpTableT[terminalIndex].dpTableCellSeed;
                }
            }
        }

        if (bestDaForDb != -1)
        {
            seed.seedPsi = db;
            seed.seedOmega = bestOmegaForDb;
            seed.seedSubsequence = this->generateSeedSubsequence(bestSeedForDb);
            break;
        }
    }

    delete[] dpTableT;

    return seed;
}
