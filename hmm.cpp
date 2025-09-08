#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <iomanip>

using namespace std;

const int N = 27;   // number of hidden states
const int M = 27;  // 26 letters + space
const string ALPHABET = "abcdefghijklmnopqrstuvwxyz ";
const int MAX_T = 50000; //maximum character to read from the file

// Step struct: one time step of the sequence
struct Step {
    int obs;
    vector<double> alpha, beta, gamma;
    vector<vector<double>> diGamma;
    double c;

    Step() {
        alpha.assign(N, 0.0);
        beta.assign(N, 0.0);
        gamma.assign(N, 0.0);
        diGamma.assign(N, vector<double>(N, 0.0));
        c = 0.0;
    }
};

// Printing matrices
void printPi(const vector<double>& pi) {
    double sum = 0;
    for (auto val : pi) { cout << val << " "; sum += val; }
    cout << ", sum=" << sum << endl;
}

void printA(const vector<vector<double>>& A) {
    for (auto& row : A) {
        double sum = 0;
        for (auto val : row) { cout << val << " "; sum += val; }
        cout << ", sum=" << sum << endl;
    }
}

void printBT(const vector<vector<double>>& B) {
    for (int j = 0; j < M; j++) {
        cout << ALPHABET[j] << " ";
        for (int i = 0; i < N; i++) cout << B[i][j] << " ";
        cout << endl;
    }
}

// Initialize pi, A, B randomly
void initMatrices(vector<double>& pi, vector<vector<double>>& A, vector<vector<double>>& B) {
    srand((unsigned)time(0));

    double prob, temp;
    // pi
    prob = 1.0 / N;
    temp = 0;
    for (int i = 0; i < N; i++) {
        pi[i] = prob + ((rand() % 2 == 0 ? 1 : -1) * (double)(rand() % 8) / 8.0 * prob / 10.0);
        temp += pi[i];
    }
    for (int i = 0; i < N; i++) pi[i] /= temp;

    // A
    for (int i = 0; i < N; i++) {
        temp = 0;
        for (int j = 0; j < N; j++) {
            A[i][j] = prob + ((rand() % 2 == 0 ? 1 : -1) * (double)(rand() % 8) / 8.0 * prob / 10.0);
            temp += A[i][j];
        }
        for (int j = 0; j < N; j++) A[i][j] /= temp;
    }

    // B
    prob = 1.0 / M;
    for (int i = 0; i < N; i++) {
        temp = 0;
        for (int j = 0; j < M; j++) {
            B[i][j] = prob + ((rand() % 2 == 0 ? 1 : -1) * (double)(rand() % 8) / 8.0 * prob / 10.0);
            temp += B[i][j];
        }
        for (int j = 0; j < M; j++) B[i][j] /= temp;
    }
}

// Read observations from text file
int getObservations(const string& fname, vector<Step>& steps) {
    ifstream in(fname);
    if (!in.is_open()) {
        cerr << "Error opening " << fname << endl;
        exit(1);
    }

    string line;
    int num = 0;

    while (getline(in, line) && num < MAX_T) {
        // Handle Windows line endings (\r\n)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        bool isLastSpace = true; // treat "start of line" as if last was a space

        for (char c : line) {
            c = tolower(c);
            size_t pos = ALPHABET.find(c);

            if (pos != string::npos) {
                // if this is a space, only add it if previous wasn't space
                if (c == ' ') {
                    if (!isLastSpace) {
                        Step s;
                        s.obs = (int)pos;
                        steps.push_back(s);
                        num++;
                        isLastSpace = true;
                    }
                } else {
                    Step s;
                    s.obs = (int)pos;
                    steps.push_back(s);
                    num++;
                    isLastSpace = false;
                }

                if (num >= MAX_T) break;
            }
        }

        if (num >= MAX_T) break;

        // Add a single space for the newline, but only if last wasn't space
        if (!isLastSpace) {
            Step s;
            s.obs = (int)ALPHABET.find(' ');
            steps.push_back(s);
            num++;
        }
    }

    return num; // number of observations
}


// Forward pass
void alphaPass(vector<Step>& steps, const vector<double>& pi, const vector<vector<double>>& A, const vector<vector<double>>& B, int T) {
    double sum = 0.0;
    //compute α0(i)
    for (int i = 0; i < N; i++) {
        steps[0].alpha[i] = pi[i] * B[i][steps[0].obs];
        sum += steps[0].alpha[i];
    }
    //scale α0(i)
    steps[0].c = 1.0 / sum;
    for (int i = 0; i < N; i++) {
        steps[0].alpha[i] *= steps[0].c;
    }
    //compute αt(i)
    for (int t = 1; t < T; t++) {
        sum = 0.0;
        for (int i = 0; i < N; i++) {
            steps[t].alpha[i] = 0.0;
            for (int j = 0; j < N; j++) {
                steps[t].alpha[i] += steps[t-1].alpha[j] * A[j][i];
            }
            steps[t].alpha[i] *= B[i][steps[t].obs];
            sum += steps[t].alpha[i];
        }
        //scale αt(i)
        steps[t].c = 1.0 / sum;
        for (int i = 0; i < N; i++) {
            steps[t].alpha[i] *= steps[t].c;
        }
    }
}

// Backward pass
void betaPass(vector<Step>& steps, const vector<vector<double>>& A, const vector<vector<double>>& B, int T) {
    //beta(T-1)[i] = 1 is scaled by c(T-1)
    for (int i = 0; i < N; i++) {
        steps[T-1].beta[i] = steps[T-1].c;
    }

    for (int t = T-2; t >= 0; t--) {
        for (int i = 0; i < N; i++) {
            steps[t].beta[i] = 0.0;
            for (int j = 0; j < N; j++) {
                steps[t].beta[i] += A[i][j] * B[j][steps[t+1].obs] * steps[t+1].beta[j];
            }
            //scale beta(t)[i] with same factor as alpha(t)[i]
            steps[t].beta[i] *= steps[t].c;
        }
    }
}

// Gamma & DiGamma
void computeGammas(vector<Step>& steps, const vector<vector<double>>& A, const vector<vector<double>>& B, int T) {
    for (int t = 0; t < T-1; t++) {
        double denom = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                denom += steps[t].alpha[i] * A[i][j] * B[j][steps[t+1].obs] * steps[t+1].beta[j];
            }
        }
        for (int i = 0; i < N; i++) {
            steps[t].gamma[i] = 0.0;
            for (int j = 0; j < N; j++) {
                steps[t].diGamma[i][j] = (steps[t].alpha[i] * A[i][j] * B[j][steps[t+1].obs] * steps[t+1].beta[j]) / denom;
                steps[t].gamma[i] += steps[t].diGamma[i][j];
            }
        }
    }
    //special case for gamma(T-1)[i]
    double denom = 0.0;
    for (int i = 0; i < N; i++) {
        denom += steps[T-1].alpha[i];
    }
    for (int i = 0; i < N; i++) {
        steps[T-1].gamma[i] = steps[T-1].alpha[i] / denom;
    }
}

// Re-estimate pi
void reEstimatePi(vector<Step>& steps, vector<double>& piBar) {
    for (int i = 0; i < N; i++) {
        piBar[i] = steps[0].gamma[i];
    }
}

// Re-estimate A
void reEstimateA(vector<Step>& steps, vector<vector<double>>& Abar, int T) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double numer = 0, denom = 0;
            for (int t = 0; t < T-1; t++) {
                numer += steps[t].diGamma[i][j];
                denom += steps[t].gamma[i];
            }
            Abar[i][j] = numer / denom;
        }
    }
}

// Re-estimate B
void reEstimateB(vector<Step>& steps, vector<vector<double>>& Bbar, int T) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            double numer = 0, denom = 0;
            for (int t = 0; t < T-1; t++) {
                if (steps[t].obs == j)
                    numer += steps[t].gamma[i];
                denom += steps[t].gamma[i];
            }
            Bbar[i][j] = numer / denom;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " filename" << endl;
        return 1;
    }

    string fname = argv[1];
    int maxIters = 100;

    vector<Step> steps;
    int T = getObservations(fname, steps);
    cout << "T = " << T << endl;

    vector<double> pi(N), piBar(N);
    vector<vector<double>> A(N, vector<double>(N));
    vector<vector<double>> Abar(N, vector<double>(N));
    vector<vector<double>> B(N, vector<double>(M));
    vector<vector<double>> Bbar(N, vector<double>(M));

    initMatrices(pi, A, B);

    cout << "Initial pi:" << endl; printPi(pi);
    cout << "Initial A:" << endl; printA(A);
    cout << "Initial B^T:" << endl; printBT(B);

    double logProb = -numeric_limits<double>::infinity(), newLogProb = 0.0;

    for (int iter = 0; iter < maxIters; iter++) {
        alphaPass(steps, pi, A, B, T);
        betaPass(steps, A, B, T);
        computeGammas(steps, A, B, T);
        reEstimatePi(steps, piBar);
        reEstimateA(steps, Abar, T);
        reEstimateB(steps, Bbar, T);

        pi = piBar;
        A = Abar;
        B = Bbar;

        newLogProb = 0;
        for (int t = 0; t < T; t++) newLogProb += log(steps[t].c);
        newLogProb = -newLogProb;

        cout << fixed << setprecision(5);
        cout << "Iteration " << iter+1 << " logProb=" << newLogProb << endl;
        if (newLogProb <= logProb) break;
        logProb = newLogProb;
    }

    
    cout << "Final pi:" << endl; printPi(pi);
    cout << "Final A:" << endl; printA(A);
    cout << "Final B^T:" << endl; printBT(B);
    cout << fixed << setprecision(5);
    cout << "Final logProb=" << newLogProb << endl;

    return 0;
}
