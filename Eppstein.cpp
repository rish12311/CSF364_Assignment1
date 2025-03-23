#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <queue>
#include <iomanip>
#include <chrono>

// Global vector to store all maximal cliques
std::vector<std::vector<int>> cliques{};

// Function to compute the degeneracy ordering of a graph
std::vector<int> compute_degeneracy_ordering(
    const std::unordered_map<int, std::set<int>>& graph) {
    
    // Create a copy of the graph that we can modify
    std::unordered_map<int, std::set<int>> temp_graph = graph;
    std::vector<int> ordering;
    ordering.reserve(graph.size());
    
    // Create a priority queue to select vertices with minimum degree
    using DegreeVertex = std::pair<int, int>;
    std::priority_queue<DegreeVertex, std::vector<DegreeVertex>, std::greater<DegreeVertex>> pq;
    
    // Initialize with all vertices and their degrees
    for (const auto& [vertex, neighbors] : temp_graph) {
        pq.push({static_cast<int>(neighbors.size()), vertex});
    }
    
    // Process vertices in order of minimum degree
    while (!pq.empty()) {
        auto [degree, vertex] = pq.top();
        pq.pop();
        
        // Check if this vertex has already been processed or its degree has changed
        if (temp_graph.find(vertex) == temp_graph.end() || temp_graph[vertex].size() != degree) {
            continue;
        }
        
        // Add vertex to the ordering
        ordering.push_back(vertex);
        
        // Update degrees of neighbors
        std::vector<int> neighbors(temp_graph[vertex].begin(), temp_graph[vertex].end());
        temp_graph.erase(vertex);
        
        for (const auto& neighbor : neighbors) {
            if (temp_graph.find(neighbor) != temp_graph.end()) {
                temp_graph[neighbor].erase(vertex);
                pq.push({static_cast<int>(temp_graph[neighbor].size()), neighbor});
            }
        }
    }
    
    return ordering;
}

// Fixed version of Bron-Kerbosch with pivot
void bron_kerbosch_pivot(
    std::set<int> R,    // Current clique
    std::set<int> P,    // Candidate vertices
    std::set<int> X,    // Excluded vertices
    const std::unordered_map<int, std::set<int>>& graph) {
    
    // If both P and X are empty, R is a maximal clique
    if (P.empty() && X.empty()) {
        if (R.size() >= 2) {  // Only include cliques of size >= 2
            std::vector<int> clique(R.begin(), R.end());
            std::sort(clique.begin(), clique.end());  // Sort for consistent output
            cliques.push_back(clique);
        }
        return;
    }
    
    // Choose a pivot vertex u from P ∪ X to maximize |P ∩ N(u)|
    std::set<int> P_union_X;
    std::set_union(P.begin(), P.end(), X.begin(), X.end(), 
                   std::inserter(P_union_X, P_union_X.begin()));
    
    int pivot = -1;
    size_t max_intersection = 0;
    
    for (int u : P_union_X) {
        const auto& neighbors = graph.at(u);
        size_t intersection_size = 0;
        
        for (int v : P) {
            if (neighbors.find(v) != neighbors.end()) {
                intersection_size++;
            }
        }
        
        if (intersection_size > max_intersection) {
            max_intersection = intersection_size;
            pivot = u;
        }
    }
    
    // If no suitable pivot found, use first vertex in P
    if (pivot == -1 && !P.empty()) {
        pivot = *P.begin();
    }
    
    // Vertices in P that are not neighbors of pivot
    std::set<int> P_minus_N_pivot;
    
    if (pivot != -1) {
        const auto& pivot_neighbors = graph.at(pivot);
        
        for (int v : P) {
            if (pivot_neighbors.find(v) == pivot_neighbors.end()) {
                P_minus_N_pivot.insert(v);
            }
        }
    } else {
        // If no pivot, use all vertices in P
        P_minus_N_pivot = P;
    }
    
    // For each vertex v not connected to pivot
    for (int v : P_minus_N_pivot) {
        const auto& v_neighbors = graph.at(v);
        
        // R ∪ {v}
        std::set<int> R_new = R;
        R_new.insert(v);
        
        // P ∩ N(v)
        std::set<int> P_new;
        for (int p : P) {
            if (p != v && v_neighbors.find(p) != v_neighbors.end()) {
                P_new.insert(p);
            }
        }
        
        // X ∩ N(v)
        std::set<int> X_new;
        for (int x : X) {
            if (v_neighbors.find(x) != v_neighbors.end()) {
                X_new.insert(x);
            }
        }
        
        // Recursive call
        bron_kerbosch_pivot(R_new, P_new, X_new, graph);
        
        // Move v from P to X
        P.erase(v);
        X.insert(v);
    }
}

// Corrected degeneracy-based Bron-Kerbosch algorithm
void bron_kerbosch_degeneracy(const std::unordered_map<int, std::set<int>>& graph) {
    std::vector<int> degeneracy_order = compute_degeneracy_ordering(graph);
    
    // Process vertices in degeneracy order
    for (size_t i = 0; i < degeneracy_order.size(); ++i) {
        int v = degeneracy_order[i];
        
        // Create set P of neighbors that come later in the ordering
        std::set<int> P;
        for (int neighbor : graph.at(v)) {
            auto it = std::find(degeneracy_order.begin() + i + 1, degeneracy_order.end(), neighbor);
            if (it != degeneracy_order.end()) {
                P.insert(neighbor);
            }
        }
        
        // Create set X of neighbors that come earlier in the ordering
        std::set<int> X;
        for (int neighbor : graph.at(v)) {
            auto it = std::find(degeneracy_order.begin(), degeneracy_order.begin() + i, neighbor);
            if (it != degeneracy_order.end()) {
                X.insert(neighbor);
            }
        }
        
        // Call Bron-Kerbosch with singleton clique {v}
        std::set<int> R = {v};
        bron_kerbosch_pivot(R, P, X, graph);
    }
}

// Function to load and create a graph from a file
std::unordered_map<int, std::set<int>> load_graph(const std::string& filename) {
    std::unordered_map<int, std::set<int>> graph;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return graph;
    }
    
    std::string line;
    
    // Skip comment lines that start with #
    while (std::getline(file, line)) {
        if (line.empty() || line[0] != '#') {
            break;
        }
    }
    
    // Process the data lines
    do {
        if (!line.empty() && line[0] != '#') {
            std::istringstream iss(line);
            int from, to;
            
            if (iss >> from >> to) {
                graph[from].insert(to);
                graph[to].insert(from); // Make the graph undirected
            }
        }
    } while (std::getline(file, line));
    
    file.close();
    
    // Ensure all vertices have entry in graph
    for (const auto& [vertex, neighbors] : graph) {
        for (int neighbor : neighbors) {
            if (graph.find(neighbor) == graph.end()) {
                graph[neighbor] = std::set<int>();
            }
        }
    }
    
    return graph;
}

// Function to create a simple test graph
std::unordered_map<int, std::set<int>> create_test_graph() {
    std::unordered_map<int, std::set<int>> graph;
    
    // Define the edges for a simple test case
    std::vector<std::pair<int, int>> edges = {
        {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}, // Clique of size 4
        {4, 5}, {5, 6}, {5, 7}, {6, 7},                  // Triangle 5-6-7 connected to clique
        {7, 8}, {8, 9}, {9, 10}, {10, 7}                 // 4-cycle 7-8-9-10
    };
    
    // Add edges to create an undirected graph
    for (const auto& [u, v] : edges) {
        graph[u].insert(v);
        graph[v].insert(u);
    }
    
    return graph;
}

// Print clique information and histogram
void print_clique_stats() {
    if (cliques.empty()) {
        std::cout << "No maximal cliques found." << std::endl;
        return;
    }
    
    // Find the maximum clique size
    size_t max_size = 0;
    for (const auto& clique : cliques) {
        max_size = std::max(max_size, clique.size());
    }
    
    // Count cliques by size
    std::vector<int> size_counts(max_size + 1, 0);
    for (const auto& clique : cliques) {
        size_counts[clique.size()]++;
    }
    
    // Print summary
    std::cout << "Total Maximal Cliques: " << cliques.size() << std::endl;
    std::cout << "Largest Clique Size: " << max_size << std::endl;
    std::cout << std::endl;
    
    // Print histogram
    std::cout << "Clique Size Histogram:" << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << "Size | Count  | Histogram" << std::endl;
    std::cout << "----------------------" << std::endl;
    
    // Calculate scale for histogram
    int max_count = *std::max_element(size_counts.begin(), size_counts.end());
    int scale_factor = (max_count > 50) ? (max_count / 50) : 1;
    
    for (size_t i = 2; i <= max_size; ++i) {  // Start from size 2
        int count = size_counts[i];
        int stars = count / scale_factor;
        if (count > 0 && stars == 0) stars = 1;  // Ensure at least one star for non-zero counts
        
        std::cout << std::setw(4) << i << " | " << std::setw(6) << count << " | ";
        for (int j = 0; j < stars; ++j) {
            std::cout << "*";
        }
        std::cout << std::endl;
    }
    
    std::cout << "----------------------" << std::endl;
    std::cout << "Scale: Each * represents approximately " << scale_factor << " clique(s)" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string filename;
    std::unordered_map<int, std::set<int>> graph;
    
    // Check if filename was provided as command line argument
    if (argc > 1) {
        filename = argv[1];
        graph = load_graph(filename);
        
        if (graph.empty()) {
            std::cout << "Error: Could not load graph from file " << filename << std::endl;
            std::cout << "Using test graph instead." << std::endl;
            graph = create_test_graph();
        } else {
            std::cout << "Using graph from " << filename << " with " << graph.size() << " vertices." << std::endl;
        }
    } else {
        // No filename provided, try default or use test graph
        filename = "email-Enron.txt";
        graph = load_graph(filename);
        
        if (graph.empty()) {
            std::cout << "Default file not found. Using test graph instead." << std::endl;
            graph = create_test_graph();
        } else {
            std::cout << "Using graph from " << filename << " with " << graph.size() << " vertices." << std::endl;
        }
    }
    
    // Clear any previous results
    cliques.clear();
    
    // Measure execution time
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Run the algorithm
    bron_kerbosch_degeneracy(graph);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    
    // Print results
    print_clique_stats();
    
    std::cout << "Execution Time: " << std::fixed << std::setprecision(3) 
              << elapsed_seconds.count() << " seconds" << std::endl;
    
    // Find the maximum clique(s)
    size_t max_clique_size = 0;
    for (const auto& clique : cliques) {
        max_clique_size = std::max(max_clique_size, clique.size());
    }
    
    // Print all maximal cliques of the maximum size
    std::cout << "\nMaximal Cliques (size " << max_clique_size << "):" << std::endl;
    for (const auto& clique : cliques) {
        if (clique.size() == max_clique_size) {
            std::cout << "[ ";
            for (int v : clique) {
                std::cout << v << " ";
            }
            std::cout << "]" << std::endl;
        }
    }
    
    // Optionally print all maximal cliques if there aren't too many
    const size_t MAX_CLIQUES_TO_PRINT = 20;
    if (cliques.size() <= MAX_CLIQUES_TO_PRINT && cliques.size() > 0) {
        std::cout << "\nAll Maximal Cliques:" << std::endl;
        for (const auto& clique : cliques) {
            std::cout << "[ ";
            for (int v : clique) {
                std::cout << v << " ";
            }
            std::cout << "]" << std::endl;
        }
    }
    
    return 0;
}