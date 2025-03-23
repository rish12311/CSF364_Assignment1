#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <memory>
#include <queue>
#include <string>
#include <limits>

// Memory-efficient graph representation
class Graph {
private:
    int n; // Number of vertices
    std::vector<std::unordered_set<int>> adjacencyList; // Using unordered_set for faster lookups

public:
    Graph(int numVertices) : n(numVertices) {
        adjacencyList.resize(n);
        if (n > 1000000) {
            int estimatedAvgDegree = 30;
            for (auto& neighbors : adjacencyList) {
                neighbors.reserve(estimatedAvgDegree);
            }
        }
    }

    void addEdge(int u, int v) {
        if (u >= 0 && v >= 0 && u < n && v < n && u != v) {
            adjacencyList[u].insert(v);
            adjacencyList[v].insert(u);
        }
    }

    int getNumVertices() const { return n; }

    bool hasEdge(int u, int v) const {
        if (u >= 0 && v >= 0 && u < n && v < n) {
            return adjacencyList[u].find(v) != adjacencyList[u].end();
        }
        return false;
    }

    const std::unordered_set<int>& getNeighbors(int v) const { return adjacencyList[v]; }

    size_t getDegree(int v) const {
        if (v >= 0 && v < n) {
            return adjacencyList[v].size();
        }
        return 0;
    }

    static Graph loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::string line;
        int maxVertexId = -1;
        size_t numEdges = 0;
        std::unordered_set<int> uniqueVertices;
        
        std::cout << "Scanning file to determine graph size..." << std::endl;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            int u, v;
            if (iss >> u >> v) {
                uniqueVertices.insert(u);
                uniqueVertices.insert(v);
                maxVertexId = std::max(maxVertexId, std::max(u, v));
                numEdges++;
            }
        }
        
        int numVertices = maxVertexId + 1;  // Fixed: Use maxVertexId + 1 instead of uniqueVertices.size()
        std::cout << "Detected " << uniqueVertices.size() << " vertices with max id " << maxVertexId 
                  << " and approximately " << numEdges << " edges" << std::endl;
        
        Graph graph(numVertices);
        
        file.clear();
        file.seekg(0);
        
        std::cout << "Loading edges into graph..." << std::endl;
        
        size_t edgeCount = 0;
        size_t reportInterval = numEdges / 10;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            int u, v;
            if (iss >> u >> v) {
                graph.addEdge(u, v);
                edgeCount++;
                
                if (reportInterval > 0 && edgeCount % reportInterval == 0) {
                    std::cout << "Loaded " << edgeCount << " edges (" 
                              << (edgeCount * 100 / numEdges) << "%)" << std::endl;
                }
            }
        }
        
        return graph;
    }

    Graph pruneGraph(int minDegree) const {
        std::vector<int> vertexMap(n, -1);
        int newVertexCount = 0;
        
        for (int i = 0; i < n; i++) {
            if (getDegree(i) >= minDegree) {
                vertexMap[i] = newVertexCount++;
            }
        }
        
        Graph prunedGraph(newVertexCount);
        
        for (int oldU = 0; oldU < n; oldU++) {
            int newU = vertexMap[oldU];
            if (newU != -1) {
                for (int oldV : adjacencyList[oldU]) {
                    int newV = vertexMap[oldV];
                    if (newV != -1 && newU < newV) {
                        prunedGraph.addEdge(newU, newV);
                    }
                }
            }
        }
        
        return prunedGraph;
    }
};

class MaximalCliqueEnumerator {
private:
    const Graph& graph;
    std::vector<int> Q;
    std::vector<std::vector<int>> maximalCliques;
    std::vector<int> largestClique;
    size_t maximalCliqueCount;
    size_t maxCliqueSize;
    bool storeCliques;
    bool verbose;
    std::unordered_map<int, int> sizeDistribution;
    
    // Fixed: Improved pivot selection function
    int findPivotVertex(const std::vector<int>& SUBG, const std::vector<int>& CAND) {
        int maxConnections = -1;
        int pivotVertex = -1;

        // Try to find a vertex in SUBG that has the most connections to CAND
        for (int v : SUBG) {
            const auto& neighbors = graph.getNeighbors(v);
            int connections = 0;
            
            for (int u : CAND) {
                if (neighbors.find(u) != neighbors.end()) {
                    connections++;
                }
            }
            
            if (connections > maxConnections) {
                maxConnections = connections;
                pivotVertex = v;
            }
        }
        
        return pivotVertex;
    }

    void EXPAND(std::vector<int>& SUBG, std::vector<int>& CAND) {
        // Fixed: Check if both CAND and SUBG are empty, then Q is a maximal clique
        if (CAND.empty() && SUBG.empty()) {
            maximalCliqueCount++;
            sizeDistribution[Q.size()]++;
            
            if (Q.size() > maxCliqueSize) {
                maxCliqueSize = Q.size();
                largestClique = Q;
            }
            
            if (storeCliques) {
                maximalCliques.push_back(Q);
            }
            
            if (verbose && maximalCliqueCount % 10000 == 0) {
                std::cout << "Found " << maximalCliqueCount << " maximal cliques so far. "
                          << "Current max size: " << maxCliqueSize << std::endl;
            }
            return;
        }
        
        // If CAND is empty but SUBG is not, then Q is not maximal
        if (CAND.empty()) {
            return;
        }
        
        int u = findPivotVertex(SUBG, CAND);
        
        std::vector<int> EXTu;
        if (u != -1) {
            const auto& neighborsOfU = graph.getNeighbors(u);
            
            // Find vertices in CAND that are not connected to the pivot
            for (int v : CAND) {
                if (neighborsOfU.find(v) == neighborsOfU.end()) {
                    EXTu.push_back(v);
                }
            }
        } else {
            // If no pivot is found, use all of CAND
            EXTu = CAND;
        }
        
        // If no candidates remain after pivot optimization, check if current clique is maximal
        if (EXTu.empty()) {
            // If there are no vertices in SUBG that can extend Q, then Q is maximal
            bool isMaximal = true;
            for (int v : SUBG) {
                const auto& neighbors = graph.getNeighbors(v);
                bool canExtend = true;
                
                for (int q : Q) {
                    if (neighbors.find(q) == neighbors.end()) {
                        canExtend = false;
                        break;
                    }
                }
                
                if (canExtend) {
                    isMaximal = false;
                    break;
                }
            }
            
            if (isMaximal && !Q.empty()) {
                maximalCliqueCount++;
                sizeDistribution[Q.size()]++;
                
                if (Q.size() > maxCliqueSize) {
                    maxCliqueSize = Q.size();
                    largestClique = Q;
                }
                
                if (storeCliques) {
                    maximalCliques.push_back(Q);
                }
            }
            return;
        }
        
        while (!EXTu.empty()) {
            int q = EXTu.back();
            EXTu.pop_back();
            
            // Add vertex q to the current clique
            Q.push_back(q);
            
            // Create new candidate and subgraph sets for recursive call
            std::vector<int> SUBGq;
            std::vector<int> CANDq;
            
            const auto& neighborsOfQ = graph.getNeighbors(q);
            
            // Vertices in SUBG connected to q go into SUBGq
            for (int v : SUBG) {
                if (neighborsOfQ.find(v) != neighborsOfQ.end()) {
                    SUBGq.push_back(v);
                }
            }
            
            // Vertices in CAND connected to q go into CANDq
            for (int v : CAND) {
                if (v != q && neighborsOfQ.find(v) != neighborsOfQ.end()) {
                    CANDq.push_back(v);
                }
            }
            
            // Recursive call
            EXPAND(SUBGq, CANDq);
            
            // Remove q from candidates
            CAND.erase(std::remove(CAND.begin(), CAND.end(), q), CAND.end());
            
            // Add q to SUBG for future iterations
            SUBG.push_back(q);
            
            // Remove q from the current clique
            Q.pop_back();
        }
    }

public:
    MaximalCliqueEnumerator(const Graph& g, bool storeAllCliques = false, bool verboseOutput = false)
        : graph(g), maximalCliqueCount(0), maxCliqueSize(0), storeCliques(storeAllCliques), 
          verbose(verboseOutput) {}

    void enumerateMaximalCliques(bool pruneGraph = true, int degreeThreshold = 5) {
        auto startTime = std::chrono::high_resolution_clock::now();
        
        const Graph* workingGraph = &graph;
        std::unique_ptr<Graph> prunedGraphPtr;
        
        if (pruneGraph && graph.getNumVertices() > 10000) {
            std::cout << "Pruning vertices with degree < " << degreeThreshold << "..." << std::endl;
            prunedGraphPtr = std::make_unique<Graph>(graph.pruneGraph(degreeThreshold));
            workingGraph = prunedGraphPtr.get();
            std::cout << "Pruned graph has " << workingGraph->getNumVertices() << " vertices" << std::endl;
        }
        
        std::vector<int> SUBG;
        std::vector<int> CAND;
        for (int i = 0; i < workingGraph->getNumVertices(); i++) {
            // Only include vertices that actually exist in the graph
            if (workingGraph->getDegree(i) > 0) {
                CAND.push_back(i);
            }
        }
        
        std::cout << "Starting maximal clique enumeration..." << std::endl;
        
        EXPAND(SUBG, CAND);
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        double seconds = duration.count() / 1000.0;

        std::cout << "\nTotal Maximal Cliques: " << maximalCliqueCount << std::endl;
        std::cout << "Largest Clique Size: " << maxCliqueSize << std::endl;
        
        std::cout << "\nClique Size Histogram:" << std::endl;
        std::cout << "--------------------" << std::endl;
        std::cout << "Size | Count | Histogram" << std::endl;
        std::cout << "--------------------" << std::endl;
        
        int maxCount = 0;
        for (const auto& entry : sizeDistribution) {
            maxCount = std::max(maxCount, entry.second);
        }
        
        int scaleFactor = 1;  // Fixed: Start with scale factor of 1
        if (maxCount > 100) {
            scaleFactor = 5;
        }
        if (maxCount > 500) {
            scaleFactor = 50;
        }
        if (maxCount > 5000) {
            scaleFactor = 500;
        }

        std::vector<std::pair<int, int>> sortedDist(sizeDistribution.begin(), sizeDistribution.end());
        std::sort(sortedDist.begin(), sortedDist.end());
        
        for (const auto& entry : sortedDist) {
            int size = entry.first;
            int count = entry.second;
            int stars = (count + scaleFactor - 1) / scaleFactor;  // Round up
            
            std::cout << std::setw(4) << size << " | " 
                      << std::setw(6) << count << " | ";
            
            for (int i = 0; i < stars; i++) {
                std::cout << "*";
            }
            std::cout << std::endl;
        }
        
        std::cout << "--------------------" << std::endl;
        std::cout << "Scale: Each * represents " << scaleFactor << " clique(s)" << std::endl;

        std::cout << "\nExecution Time: " << std::fixed << std::setprecision(3) << seconds << " seconds" << std::endl;
        
        if (!largestClique.empty()) {
            std::cout << "\nLargest Clique Vertices: ";
            for (int vertex : largestClique) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::vector<int> getLargestClique() const { return largestClique; }
    
    size_t getLargestCliqueSize() const { return maxCliqueSize; }
    
    // Added: Function to get the clique count for a specific size
    int getCliqueCountForSize(int size) const {
        auto it = sizeDistribution.find(size);
        if (it != sizeDistribution.end()) {
            return it->second;
        }
        return 0;
    }
    
    // Added: Function to get all clique counts
    const std::unordered_map<int, int>& getCliqueDistribution() const {
        return sizeDistribution;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <graph_file> [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --store-cliques    Store all cliques in memory" << std::endl;
        std::cerr << "  --verbose          Show detailed progress information" << std::endl;
        std::cerr << "  --degree-threshold=N  Prune vertices with degree less than N (default: 5)" << std::endl;
        std::cerr << "  --no-prune         Don't prune low-degree vertices" << std::endl;
        return 1;
    }
    
    bool storeCliques = true;
    bool verbose = false;
    bool pruneGraph = true;
    int degreeThreshold = 5;
    
    for (int i = 2; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--store-cliques") {
            storeCliques = true;
        } else if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--no-prune") {
            pruneGraph = false;
        } else if (arg.find("--degree-threshold=") == 0) {
            degreeThreshold = std::stoi(arg.substr(18));
        }
    }
    
    try {
        std::cout << argv[0] << std::endl;
        
        // Check if input is from a file or standard input
        if (std::string(argv[1]) == "-") {
            // Create a graph from standard input
            std::vector<std::pair<int, int>> edges;
            int maxVertex = -1;
            std::string line;
            
            std::cout << "Reading graph from standard input..." << std::endl;
            
            while (std::getline(std::cin, line)) {
                if (line.empty() || line[0] == '#') continue;
                
                std::istringstream iss(line);
                int u, v;
                if (iss >> u >> v) {
                    edges.push_back({u, v});
                    maxVertex = std::max(maxVertex, std::max(u, v));
                }
            }
            
            Graph graph(maxVertex + 1);
            for (const auto& edge : edges) {
                graph.addEdge(edge.first, edge.second);
            }
            
            std::cout << "Graph created with " << maxVertex + 1 << " vertices and " 
                      << edges.size() << " edges" << std::endl;
            
            MaximalCliqueEnumerator enumerator(graph, storeCliques, verbose);
            enumerator.enumerateMaximalCliques(pruneGraph, degreeThreshold);
        } else {
            // Load graph from file
            std::cout << "Loading graph from " << argv[1] << "..." << std::endl;
            
            auto loadStart = std::chrono::high_resolution_clock::now();
            Graph graph = Graph::loadFromFile(argv[1]);
            auto loadEnd = std::chrono::high_resolution_clock::now();
            
            std::cout << "Graph loaded with " << graph.getNumVertices() << " vertices" << std::endl;
            std::cout << "Loading time: " 
                      << std::chrono::duration_cast<std::chrono::seconds>(loadEnd - loadStart).count() 
                      << " seconds" << std::endl;
            
            MaximalCliqueEnumerator enumerator(graph, storeCliques, verbose);
            enumerator.enumerateMaximalCliques(pruneGraph, degreeThreshold);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}