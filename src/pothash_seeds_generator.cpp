#include "pothash.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    if (argc != 6)
    {
        cerr << "Invalid number of command-line arguments provided!\nCorrect usage: ./pothash_seeds_generator [datasetFileName] [paramD] [subseqSizeK] [numRepeats] [seedGenerationMode]" << endl;

        return EXIT_FAILURE;
    }

    string datasetFileName = argv[1];

    string datasetFilePath = string("..") + filesystem::path::preferred_separator + string("data") + filesystem::path::preferred_separator + datasetFileName;

    ifstream datasetFile(datasetFilePath);

    if (!datasetFile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    string seedsFileFolderPath = string("..") + filesystem::path::preferred_separator + string("seeds");

    filesystem::path seedsFileFolder(seedsFileFolderPath);

    if (!filesystem::exists(seedsFileFolder))
    {
        if (!filesystem::create_directories(seedsFileFolder))
        {
            exit(EXIT_FAILURE);
        }
    }

    int paramD = stoi(argv[2]), subseqSizeK = stoi(argv[3]), numRepeats = stoi(argv[4]), seedGenerationMode = stoi(argv[5]);

    map<char, int> defaultAlphabet = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    string seedsFilePath = seedsFileFolderPath + filesystem::path::preferred_separator + string("seeds_") + datasetFileName + string("_D") + to_string(paramD) + string("_K") + to_string(subseqSizeK) + string("_Sigma") + to_string(defaultAlphabet.size()) + string("_numRepeats") + to_string(numRepeats) + string("_seedGenerationMode") + to_string(seedGenerationMode);

    ofstream seedsFile(seedsFilePath);

    if (!seedsFile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    cout << "Instantiating a PotHash object with paramD = " << paramD << " and subseqSizeK = " << subseqSizeK << "..." << endl;

    PotHash potHash(paramD, subseqSizeK, defaultAlphabet);

    string sequence;

    cout << "Generating seeds for dataset \"" << datasetFileName << "\" with PotHash object having paramD = " << paramD << ", subseqSizeK = " << subseqSizeK << " and numRepeats = " << numRepeats << " in seedGenerationMode = " << seedGenerationMode << "..." << endl;

    while (datasetFile >> sequence)
    {
        seedsFile << sequence << endl;

        for (int repeat = 1; repeat <= numRepeats; repeat++)
        {
            potHash.loadTables(repeat);

            vector<Seed> seeds = potHash.solveDP(sequence, seedGenerationMode);

            for (int seedIdx = 0; seedIdx < seeds.size(); seedIdx++)
            {
                if (seeds[seedIdx].seedParamD != -1)
                {
                    seedsFile << repeat << "," << seeds[seedIdx].seedParamD << "," << seeds[seedIdx].seedScoreVector.vectorComponentX << "," << seeds[seedIdx].seedScoreVector.vectorComponentY << "," << seeds[seedIdx].seedScoreVector.magnitude() << "," << seeds[seedIdx].seedSubsequence << endl;
                }
                else
                {
                    seedsFile << repeat << "," << seeds[seedIdx].seedParamD << "," << -1 << "," << -1 << "," << -1 << "," << seeds[seedIdx].seedSubsequence << endl;
                }
            }
        }
    }

    datasetFile.close();

    seedsFile.close();

    cout << "Done!" << endl;

    return EXIT_SUCCESS;
}