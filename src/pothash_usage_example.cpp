#include "pothash.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    int paramD = 16, subseqSizeK = 16;

    map<char, int> defaultAlphabet = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    cout << "Instantiating a PotHash object with paramD = " << paramD << " and subseqSizeK = " << subseqSizeK << "..." << endl;

    PotHash potHash(paramD, subseqSizeK, defaultAlphabet);

    bool shouldSaveTables = false;

    int tablesFileVersion = 0;

    cout << "\nGenerating tables with PotHash object..." << endl;

    potHash.generateTables(shouldSaveTables, tablesFileVersion);

    string sequence = "ACGTACGTACGTACGTACGT";

    cout << "\nGenerating seeds for sequence \"" << sequence << "\" with PotHash object..." << endl;

    vector<Seed> seeds = potHash.solveDP(sequence, 0); // seedGenerationMode = 0 means that valid seed with minimum d will be picked
    
    cout << "\nSeeds picked in seedGenerationMode = 0: {seedParamD, seedScoreVector, magnitude(seedScoreVector), seedSubsequence}\n" << endl;

    for (int seedIdx = 0; seedIdx < seeds.size(); seedIdx++)
    {
        if (seeds[seedIdx].seedParamD != -1)
        {
            cout << "-> " << seeds[seedIdx].seedParamD << ", <" << seeds[seedIdx].seedScoreVector.vectorComponentX << ", " << seeds[seedIdx].seedScoreVector.vectorComponentY << ">, " << seeds[seedIdx].seedScoreVector.magnitude() << ", " << seeds[seedIdx].seedSubsequence << endl;
        }
        else
        {
            cout << "-> \"No valid seed found\"" << endl;
        }
    }

    seeds = potHash.solveDP(sequence, 1); // seedGenerationMode = 1 means that valid seed with maximum magnitude of score vector will be picked

    cout << "\nSeeds picked in seedGenerationMode = 1: {seedParamD, seedScoreVector, magnitude(seedScoreVector), seedSubsequence}\n" << endl;

    for (int seedIdx = 0; seedIdx < seeds.size(); seedIdx++)
    {
        if (seeds[seedIdx].seedParamD != -1)
        {
            cout << "-> " << seeds[seedIdx].seedParamD << ", <" << seeds[seedIdx].seedScoreVector.vectorComponentX << ", " << seeds[seedIdx].seedScoreVector.vectorComponentY << ">, " << seeds[seedIdx].seedScoreVector.magnitude() << ", " << seeds[seedIdx].seedSubsequence << endl;
        }
        else
        {
            cout << "-> \"No valid seed found\"" << endl;
        }
    }

    seeds = potHash.solveDP(sequence, 2); // seedGenerationMode = 2 means that all seeds will be picked regardless of them being valid or invalid

    cout << "\nSeeds picked in seedGenerationMode = 2: {seedParamD, seedScoreVector, magnitude(seedScoreVector), seedSubsequence}\n" << endl;

    for (int seedIdx = 0; seedIdx < seeds.size(); seedIdx++)
    {
        if (seeds[seedIdx].seedParamD != -1)
        {
            cout << "-> " << seeds[seedIdx].seedParamD << ", <" << seeds[seedIdx].seedScoreVector.vectorComponentX << ", " << seeds[seedIdx].seedScoreVector.vectorComponentY << ">, " << seeds[seedIdx].seedScoreVector.magnitude() << ", " << seeds[seedIdx].seedSubsequence << endl;
        }
        else
        {
            cout << "-> \"No valid seed found\"" << endl;
        }
    }
    
    cout << "\nDone!" << endl;

    return EXIT_SUCCESS;
}