#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <deque>
#include <algorithm>
#include <chrono>
#include <string>
#include <sstream>
#include <iomanip>  // Added for setw formatting

// Custom hash function for std::set to use in unordered_set
struct SetHash {
    template <typename T>
    std::size_t operator()(const std::set<T>& s) const {
        std::size_t seed = s.size();
        for (const auto& i : s) {
            seed ^= std::hash<T>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Custom equality function for std::set to use in unordered_set
struct SetEqual {
    template <typename T>
    bool operator()(const std::set<T>& lhs, const std::set<T>& rhs) const {
        return lhs == rhs;
    }
};

/**
 * Read a graph from a text file where each line contains an edge "u v".
 * Convert it to an undirected graph with zero-based vertex indices.
 * 
 * @param filename Path to the text file containing edge data
 * @param graph Output parameter: A map representing the undirected graph with zero-based indices
 * @param idxToVertex Output parameter: A map from zero-based indices to original vertex IDs
 * @return True if successful, false otherwise
 */
bool readGraphFromFile(
    const std::string& filename,
    std::unordered_map<int, std::vector<int>>& graph,
    std::unordered_map<int, std::string>& idxToVertex
) {
    // Step 1: Read the edges and identify all vertices
    std::vector<std::pair<std::string, std::string>> edges;
    std::set<std::string> vertices;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        
        // Parse the edge
        std::istringstream iss(line);
        std::string u, v;
        if (iss >> u >> v) {
            edges.push_back({u, v});
            vertices.insert(u);
            vertices.insert(v);
        }
    }
    
    // Step 2: Map the original vertex IDs to zero-based indices
    std::unordered_map<std::string, int> vertexToIdx;
    int index = 0;
    for (const auto& v : vertices) {
        vertexToIdx[v] = index;
        idxToVertex[index] = v;
        index++;
    }
    
    // Step 3: Create the undirected graph with zero-based indices
    int n = vertices.size();
    for (int i = 0; i < n; i++) {
        graph[i] = std::vector<int>();
    }
    
    for (const auto& edge : edges) {
        int uIdx = vertexToIdx[edge.first];
        int vIdx = vertexToIdx[edge.second];
        
        // Add edges for both directions (undirected graph)
        if (std::find(graph[uIdx].begin(), graph[uIdx].end(), vIdx) == graph[uIdx].end()) {
            graph[uIdx].push_back(vIdx);
        }
        if (std::find(graph[vIdx].begin(), graph[vIdx].end(), uIdx) == graph[vIdx].end()) {
            graph[vIdx].push_back(uIdx);
        }
    }
    
    return true;
}

/**
 * Iterative implementation of the CLIQUE algorithm for finding all unique maximal cliques.
 * 
 * @param graph A map representing the graph where keys are vertices (0 to n-1)
 *              and values are lists of adjacent vertices.
 * @return A vector of unique maximal cliques of maximum size (each represented as a vector of vertices).
 */
std::vector<std::vector<int>> findUniqueMaxCliques(const std::unordered_map<int, std::vector<int>>& graph) {
    // Sort vertices by degree (non-decreasing order)
    std::vector<int> vertices;
    for (const auto& entry : graph) {
        vertices.push_back(entry.first);
    }
    
    std::sort(vertices.begin(), vertices.end(), [&graph](int a, int b) {
        return graph.at(a).size() < graph.at(b).size();
    });
    
    // Map to 1-based indexing for the algorithm (as it expects 1-based)
    std::unordered_map<int, int> idx0to1;
    std::unordered_map<int, int> idx1to0;
    for (int v : vertices) {
        idx0to1[v] = v + 1;
        idx1to0[v + 1] = v;
    }
    
    // Convert graph to use 1-based indices for internal algorithm
    std::unordered_map<int, std::vector<int>> indexedGraph;
    for (int v : vertices) {
        std::vector<int> neighbors;
        for (int u : graph.at(v)) {
            neighbors.push_back(idx0to1[u]);
        }
        indexedGraph[idx0to1[v]] = neighbors;
    }
    
    int n = vertices.size();
    
    // Use a set to store unique cliques
    std::unordered_set<std::set<int>, SetHash, SetEqual> uniqueCliquesSet;
    
    auto neighbors = [&indexedGraph](int v) -> const std::vector<int>& {
        return indexedGraph.at(v);
    };
    
    // Stack for iterative implementation
    // Each item contains (i, C, state)
    // state: 0 = first branch, 1 = second branch
    std::deque<std::tuple<int, std::set<int>, int>> stack;
    stack.push_back(std::make_tuple(1, std::set<int>(), 0));
    
    while (!stack.empty()) {
        auto [i, C, state] = stack.back();
        stack.pop_back();
        
        // If we've already processed the first branch and now coming back
        // for the second branch, we need to decide if vertex i should be included
        if (state == 1) {
            // Get neighbors of i
            std::set<int> Ni(neighbors(i).begin(), neighbors(i).end());
            
            // C ∩ N(i)
            std::set<int> C_intersect_Ni;
            for (int v : C) {
                if (Ni.find(v) != Ni.end()) {
                    C_intersect_Ni.insert(v);
                }
            }
            
            // Initialize T and S
            std::unordered_map<int, int> T, S;
            for (int v = 1; v <= n; v++) {
                T[v] = 0;
                S[v] = 0;
            }
            
            // Compute T[y] = |N(y) ∩ C ∩ N(i)| for all vertices y
            for (int x : C_intersect_Ni) {
                for (int y : neighbors(x)) {
                    if (C.find(y) == C.end() && y != i) {
                        T[y]++;
                    }
                }
            }
            
            // Compute S[y] = |N(y) ∩ (C - N(i))| for all vertices y
            std::set<int> C_minus_Ni;
            for (int v : C) {
                if (Ni.find(v) == Ni.end()) {
                    C_minus_Ni.insert(v);
                }
            }
            
            for (int x : C_minus_Ni) {
                for (int y : neighbors(x)) {
                    if (C.find(y) == C.end() && y != i) {
                        S[y]++;
                    }
                }
            }
            
            // Apply the maximality and lexicographic tests
            bool canIncludeI = true;
            
            // Maximality test
            for (int y : Ni) {
                if (C.find(y) != C.end() && y < i && T[y] == C_intersect_Ni.size()) {
                    canIncludeI = false;
                    break;
                }
            }
            
            if (canIncludeI) {
                // Lexicographic test
                std::vector<int> C_minus_Ni_list(C_minus_Ni.begin(), C_minus_Ni.end());
                std::sort(C_minus_Ni_list.begin(), C_minus_Ni_list.end());
                int p = C_minus_Ni_list.size();
                
                // Process vertices y with S[y] > 0
                for (int y = 1; y <= n; y++) {
                    if (y < i && C.find(y) == C.end() && T[y] == C_intersect_Ni.size() && S[y] > 0) {
                        // Find k such that y is adjacent to exactly k vertices in C - N(i)
                        int adjacencyCount = 0;
                        for (int k = 0; k < p; k++) {
                            int vertex = C_minus_Ni_list[k];
                            const auto& vertexNeighbors = neighbors(vertex);
                            if (std::find(vertexNeighbors.begin(), vertexNeighbors.end(), y) != vertexNeighbors.end()) {
                                adjacencyCount++;
                            }
                        }
                        
                        // Lexicographic test condition
                        if (adjacencyCount > 0) {
                            // Find first vertex in C_minus_Ni that y is not adjacent to
                            for (int k = 0; k < p; k++) {
                                int vertex = C_minus_Ni_list[k];
                                const auto& vertexNeighbors = neighbors(vertex);
                                if (std::find(vertexNeighbors.begin(), vertexNeighbors.end(), y) == vertexNeighbors.end()) {
                                    if (y > vertex) {
                                        canIncludeI = false;
                                        break;
                                    }
                                }
                            }
                            if (!canIncludeI) break;
                        }
                    }
                }
                
                // Process vertices y with S[y] = 0
                if (canIncludeI && !C_intersect_Ni.empty()) {
                    for (int y = 1; y <= n; y++) {
                        if (y < i && C.find(y) == C.end() && y != i && 
                            T[y] == C_intersect_Ni.size() && S[y] == 0) {
                            // For this case, we just need to check if y is larger than
                            // the largest vertex in C - N(i)
                            if (p > 0 && y > C_minus_Ni_list[p-1]) {
                                canIncludeI = false;
                                break;
                            }
                        }
                    }
                }
            }
            
            // If tests pass, include vertex i in the clique
            if (canIncludeI) {
                std::set<int> newC;
                for (int v : C) {
                    if (Ni.find(v) != Ni.end()) {
                        newC.insert(v);
                    }
                }
                newC.insert(i);
                stack.push_back(std::make_tuple(i + 1, newC, 0));
            }
            
            continue;
        }
        
        // First branch: continue without adding vertex i
        if (i <= n) {
            // Push the second branch onto the stack to be handled later
            stack.push_back(std::make_tuple(i, C, 1));
            
            // Then continue with the first branch (without vertex i)
            stack.push_back(std::make_tuple(i + 1, C, 0));
        } else {
            // We've considered all vertices, so C is a maximal clique
            std::set<int> cliqueVertices;
            for (int v : C) {
                cliqueVertices.insert(idx1to0[v]);
            }
            uniqueCliquesSet.insert(cliqueVertices);
        }
    }
    
    // Convert set of sets to vector of vectors
    std::vector<std::vector<int>> allCliques;
    for (const auto& clique : uniqueCliquesSet) {
        std::vector<int> cliqueVec(clique.begin(), clique.end());
        std::sort(cliqueVec.begin(), cliqueVec.end());
        allCliques.push_back(cliqueVec);
    }
    
    return allCliques;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <graph_file>" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    
    try {
        auto startTime = std::chrono::high_resolution_clock::now();
        
        std::unordered_map<int, std::vector<int>> graph;
        std::unordered_map<int, std::string> idxToVertex;
        
        if (!readGraphFromFile(filename, graph, idxToVertex)) {
            std::cerr << "Error: File '" << filename << "' not found" << std::endl;
            return 1;
        }
        
        int n = graph.size();
        std::cout << "Read graph with " << n << " vertices (mapped to 0-" << n-1 << ")" << std::endl;
        
        auto readTime = std::chrono::high_resolution_clock::now();
        auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(readTime - startTime).count() / 1000.0;
        std::cout << "Time to read graph: " << readDuration << " seconds" << std::endl;
        
        std::cout << "Finding maximal cliques..." << std::endl;
        auto cliqueStartTime = std::chrono::high_resolution_clock::now();
        
        std::vector<std::vector<int>> maxCliques = findUniqueMaxCliques(graph);
        
        auto cliqueEndTime = std::chrono::high_resolution_clock::now();
        auto cliqueDuration = std::chrono::duration_cast<std::chrono::milliseconds>(cliqueEndTime - cliqueStartTime).count() / 1000.0;
        
        if (!maxCliques.empty()) {
            // Count cliques by size
            std::unordered_map<int, int> sizeCount;
            int maxCliqueSize = 0;
            int totalCliques = maxCliques.size();
            
            for (const auto& clique : maxCliques) {
                int size = clique.size();
                sizeCount[size]++;
                maxCliqueSize = std::max(maxCliqueSize, size);
            }
            
            // Print the histogram header
            std::cout << "Total Maximal Cliques: " << totalCliques << std::endl;
            std::cout << "Largest Clique Size: " << maxCliqueSize << std::endl;
            std::cout << std::endl;
            std::cout << "Clique Size Histogram:" << std::endl;
            std::cout << "--------------------" << std::endl;
            std::cout << "Size | Count | Histogram" << std::endl;
            std::cout << "--------------------" << std::endl;
            
            // Calculate scale factor - how many cliques per asterisk
            // This can be adjusted based on your preference (5, 50, 500, 5000)
            int scaleBase;
            if (totalCliques > 50000) scaleBase = 5000;
            else if (totalCliques > 5000) scaleBase = 500;
            else if (totalCliques > 500) scaleBase = 50;
            else scaleBase = 5;
            
            // Find the maximum count to determine scale
            int maxCount = 0;
            for (int size = 2; size <= maxCliqueSize; size++) {
                if (sizeCount.find(size) != sizeCount.end()) {
                    maxCount = std::max(maxCount, sizeCount[size]);
                }
            }
            
            // Calculate appropriate scale factor
            int scale = std::max(1, static_cast<int>(std::ceil(static_cast<double>(maxCount) / 70.0 / scaleBase)) * scaleBase);
            
            // Print histogram
            for (int size = 2; size <= maxCliqueSize; size++) {
                int count = sizeCount.find(size) != sizeCount.end() ? sizeCount[size] : 0;
                
                std::cout << std::setw(3) << size << " | " 
                          << std::setw(6) << count << " | ";
                
                // Print asterisks
                int stars = std::ceil(static_cast<double>(count) / scale);
                for (int i = 0; i < stars; i++) {
                    std::cout << "*";
                }
                std::cout << std::endl;
            }
            
            std::cout << "--------------------" << std::endl;
            std::cout << "Scale: Each * represents approximately " << scale << " clique(s)" << std::endl;
            std::cout << std::endl;
            std::cout << "Execution Time: " << cliqueDuration << " seconds" << std::endl;
            
        } else {
            std::cout << "No cliques found in the graph." << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}