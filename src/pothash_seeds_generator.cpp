#include "pothash.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    if (argc != 6)
    {
        cerr << "Invalid number of command-line arguments provided!\nCorrect usage: ./pothash_seeds_generator [datasetFileName] [paramDa] [paramDb] [subseqSizeK] [numRepeats]" << endl;

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

    int paramDa = stoi(argv[2]), paramDb = stoi(argv[3]), subseqSizeK = stoi(argv[4]), numRepeats = stoi(argv[5]);

    map<char, int> defaultAlphabet = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    string seedsFilePath = seedsFileFolderPath + filesystem::path::preferred_separator + string("seeds_") + datasetFileName + string("_Da") + to_string(paramDa) + string("_Db") + to_string(paramDb) + string("_K") + to_string(subseqSizeK) + string("_Sigma") + to_string(defaultAlphabet.size()) + string("_numRepeats") + to_string(numRepeats);

    ofstream seedsFile(seedsFilePath);

    if (!seedsFile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    cout << "Instantiating a PotHash object with paramDa = " << paramDa << ", paramDb = " << paramDb << " and subseqSizeK = " << subseqSizeK << "..." << endl;

    PotHash potHash(paramDa, paramDb, subseqSizeK, defaultAlphabet);

    string sequence;

    cout << "Generating seeds for dataset \"" << datasetFileName << "\" with PotHash object having numRepeats = " << numRepeats << " (SubseqHash-style min valid psi)..." << endl;

    while (datasetFile >> sequence)
    {
        seedsFile << sequence << endl;

        for (int repeat = 1; repeat <= numRepeats; repeat++)
        {
            potHash.loadTables(repeat);

            Seed seed = potHash.solveDP(sequence);

            if (seed.seedPsi != -1)
            {
                seedsFile << repeat << "," << seed.seedPsi << "," << seed.seedOmega << "," << seed.seedSubsequence << endl;
            }
            else
            {
                seedsFile << repeat << "," << seed.seedPsi << "," << -1 << "," << seed.seedSubsequence << endl;
            }
        }
    }

    datasetFile.close();

    seedsFile.close();

    cout << "Done!" << endl;

    return EXIT_SUCCESS;
}
