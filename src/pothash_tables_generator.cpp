#include "pothash.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cerr << "Invalid number of command-line arguments provided!\nCorrect usage: ./pothash_tables_generator [paramD] [subseqSizeK] [numTablesFileVersions]" << endl;

        return EXIT_FAILURE;
    }

    int paramD = stoi(argv[1]), subseqSizeK = stoi(argv[2]);

    map<char, int> defaultAlphabet = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    cout << "Instantiating a PotHash object with paramD = " << paramD << " and subseqSizeK = " << subseqSizeK << "..." << endl;

    PotHash potHash(paramD, subseqSizeK, defaultAlphabet);
    
    int numTablesFileVersions = stoi(argv[3]);

    bool shouldSaveTables = true;

    cout << "Generating and saving " << numTablesFileVersions << " sets of tables with PotHash object..." << endl;

    for (int tablesFileVersion = 1; tablesFileVersion <= numTablesFileVersions; tablesFileVersion++)
    {
        potHash.generateTables(shouldSaveTables, tablesFileVersion);
    }

    cout << "Done!" << endl;

    return EXIT_SUCCESS;
}