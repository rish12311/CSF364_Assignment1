Maximum Clique Algorithm Implementations
This repository contains implementations of three algorithms for solving the Maximum Clique Problem: Tomita, Eppstein, and Chiba & Nishizeki. These implementations are written in C++ and optimized for performance.
Setup Requirements

C++ compiler with C++17 support
Unix-like environment (Linux, macOS, WSL, etc.)

Preprocessing
Before running the algorithms, you need to preprocess your input data:

Compile the preprocessing script:
bashCopyg++ -std=c++17 -O3 -o preprocess preprocess.cpp

Run the preprocessing script:
bashCopy./preprocess <raw_input_file>.txt <processed_output_file>.txt

Use the generated <processed_output_file>.txt as input to the algorithm executables.

Compilation
Each algorithm is implemented in a separate C++ file. Compile them using:
bashCopy# Compile Tomita Algorithm
g++ -std=c++17 -O3 -o tomita Tomita.cpp

# Compile Eppstein Algorithm
g++ -std=c++17 -O3 -o eppstein Eppstein.cpp

# Compile Chiba and Nishizeki Algorithm
g++ -std=c++17 -O3 -o chiba Chiba.cpp
Execution
Run each algorithm with your preprocessed input file:
bashCopy# Run Tomita Algorithm
./tomita <input_file>.txt

# Run Eppstein Algorithm
./eppstein <input_file>.txt

# Run Chiba and Nishizeki Algorithm
./chiba <input_file>.txt
Replace <input_file>.txt with the path to your preprocessed data file.
Algorithm Descriptions

Tomita Algorithm: An exact algorithm for finding the maximum clique in a graph using a branch and bound approach with pivot selection.
Eppstein Algorithm: An algorithm for finding maximal cliques using a recursive approach with vertex ordering.
Chiba and Nishizeki Algorithm: An algorithm designed to efficiently find all maximal cliques in sparse graphs.

Input Format
The input file should be in the following format:
Each line should have two integers u and v representing an edge between vertices u and v.

Output
Each algorithm will output:

The histogram of each clique size
The vertices in the maximum clique
The execution time
