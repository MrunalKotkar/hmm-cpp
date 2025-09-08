#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

const int SHIFT = 18;        // Caesar shift key
const string ALPHABET = "abcdefghijklmnopqrstuvwxyz";
char c;
int count = 0;
const int MAX_CHARS = 50000;

char shiftChar(char c, int shift) {
    int pos = c - 'a';
    int newPos = (pos + shift) % 26;
    return 'a' + newPos;
}

int main() {

    std::ifstream in("BrownCorpus");       // input file
    if (!in.is_open()) {
        std::cerr << "Error opening Brown Corpus File\n";
        return 1;
    }

    ofstream plainFile("plaintext.txt");
    ofstream cipherFile("ciphertext.txt");

    if (!plainFile || !cipherFile) {
        cerr << "Error: cannot open output files." << endl;
        return 1;
    }

    //Read from Brown Corpus file and shift by 18 to generate cipher text file
    while (in.get(c) && count < MAX_CHARS) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            plainFile.put(c);
            char enc = shiftChar(c, SHIFT);
            cipherFile.put(enc);
            count++;
        }
    }

    plainFile.close();
    cipherFile.close();

    cout << "Generated plaintext.txt (50k chars) and ciphertext.txt" << endl;
    return 0;
}

