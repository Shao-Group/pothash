#include "pothash.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    int paramDa = 17, paramDb = 17, subseqSizeK = 16;

    map<char, int> defaultAlphabet = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    cout << "Instantiating a PotHash object with paramDa = " << paramDa << ", paramDb = " << paramDb << " and subseqSizeK = " << subseqSizeK << "..." << endl;

    PotHash potHash(paramDa, paramDb, subseqSizeK, defaultAlphabet);

    bool shouldSaveTables = false;

    int tablesFileVersion = 0;

    cout << "\nGenerating tables with PotHash object..." << endl;

    potHash.generateTables(shouldSaveTables, tablesFileVersion);

    string sequence = "ACGTACGTACGTACGTACGT";

    cout << "\nGenerating seed for sequence \"" << sequence << "\" with PotHash object..." << endl;

    Seed seed = potHash.solveDP(sequence);

    cout << "\nSeed picked (min valid psi): {psi, omega, subsequence}\n" << endl;

    if (seed.seedPsi != -1)
    {
        cout << "-> " << seed.seedPsi << ", " << seed.seedOmega << ", " << seed.seedSubsequence << endl;
    }
    else
    {
        cout << "-> \"No valid seed found\"" << endl;
    }
    
    cout << "\nDone!" << endl;

    return EXIT_SUCCESS;
}
