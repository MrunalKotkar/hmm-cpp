#include <bits/stdc++.h>
using namespace std;

const int M = 26;      // alphabet size
char c;
const int MAX_CHARS = 1000000;
static inline int idx(char c) { return c - 'a'; }

int main() {
    int count = 0;
    vector<char> letters(MAX_CHARS);
    
    //Read 1000000 characters from Brown Corpus
    std::ifstream in("BrownCorpus");       // input file
    if (!in.is_open()) {
        std::cerr << "Error opening Brown Corpus File\n";
        return 1;
    }

    while (in.get(c) && count < MAX_CHARS) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            letters[count] = c;
            count++;
        }
    }

    // Count digraphs
    static long long counts[M][M] = {{0}};
    for (int t = 0; t < MAX_CHARS - 1; t++) {
        int i = idx(letters[t]);
        int j = idx(letters[t + 1]);
        counts[i][j]++;
    }

    // Add-5 smoothing
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            counts[i][j] += 5;
        }
    }

    // Final row-stochastic matrix
    double A[M][M];
    for (int i = 0; i < M; i++) {
        long double rowSum = 0.0L;
        for (int j = 0; j < M; j++) rowSum += counts[i][j];
        for (int j = 0; j < M; j++) {
            A[i][j] = static_cast<double>(counts[i][j] / rowSum);
        }
    }

    // Print entire A matrix
    cout << fixed << setprecision(6);
    cout << "{\n";
    for (int i = 0; i < M; i++) {
        cout << "  {";
        for (int j = 0; j < M; j++) {
            cout << A[i][j];
            if (j < M - 1) cout << ", ";
        }
        cout << "}";
        if (i < M - 1) cout << ",";
        cout << "\n";
    }
    cout << "};\n";

    return 0;
}
