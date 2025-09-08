# Hidden Markov Model in C++ (from scratch)

An implementation of Hidden Markov Models (HMMs) in C++ from scratch, as part of a machine learning course.
Includes experiments on English text and substitution cipher breaking.

# Features

- Discrete Hidden Markov Model implemented in modern C++
- Forward & Backward algorithms
- Viterbi decoding for most likely state sequence
- Baum–Welch re-estimation for training
- Config-driven experiments (N = 2, 3, 4, 27 states and cipher text cracking)

# 📂 Repository Structure
.hmm-cpp/

├── hmm.cpp                   # HMM code in C++ from scratch

├── ciphertest_problem/       # Application of HMM on ciphertext

├── datasets/                 # All required dataset files to run the code

├── docs/                     # Explanation documents

└── examples/                 # Sample outputs


# Quick Start
1. Clone & build
git clone https://github.com/your-username/hmm-cpp.git
cd hmm-cpp
g++ hmm.cpp -o hmm
hmm.exe browncorpus
