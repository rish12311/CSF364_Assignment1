# CSF364 Assignment 1 by Rishabh Goyal, Yash Gupta, Soham Mangle and Vani Jain

This repository contains implementations of three algorithms for solving the **Maximum Clique Problem**: 
- **Tomita Algorithm**
- **Eppstein Algorithm**
- **Chiba & Nishizeki Algorithm**

These implementations are written in **C++** and optimized for performance.

---
## 📌 Setup Requirements

To run these algorithms, ensure you have:
- A **C++ compiler** with C++17 support
- A **Unix-like environment** (Linux, macOS, WSL, etc.)

---
## ⚙️ Preprocessing

Before executing the algorithms, you need to preprocess your input data.

### **1️⃣ Compile the Preprocessing Script**
```bash
g++ -std=c++17 -O3 -o preprocess preprocess.cpp
```

### **2️⃣ Run the Preprocessing Script**
```bash
./preprocess <raw_input_file>.txt <processed_output_file>.txt
```

Use the generated `<processed_output_file>.txt` as input for the algorithms.

---
## 🚀 Compilation

Each algorithm is implemented in a separate C++ file. Compile them as follows:

### **Tomita Algorithm**
```bash
g++ -std=c++17 -O3 -o tomita Tomita.cpp
```

### **Eppstein Algorithm**
```bash
g++ -std=c++17 -O3 -o eppstein Eppstein.cpp
```

### **Chiba & Nishizeki Algorithm**
```bash
g++ -std=c++17 -O3 -o chiba Chiba.cpp
```

---
## ▶️ Execution

Run each algorithm using the preprocessed input file:

### **Tomita Algorithm**
```bash
./tomita <input_file>.txt
```

### **Eppstein Algorithm**
```bash
./eppstein <input_file>.txt
```

### **Chiba & Nishizeki Algorithm**
```bash
./chiba <input_file>.txt
```

Replace `<input_file>.txt` with the path to your processed data file.

---
## 📖 Algorithm Descriptions

### **1️⃣ Tomita Algorithm**
An exact algorithm for finding the **maximum clique** in a graph using a **branch and bound** approach with **pivot selection**.

### **2️⃣ Eppstein Algorithm**
A recursive algorithm for finding **maximal cliques** using an optimized **vertex ordering** strategy.

### **3️⃣ Chiba & Nishizeki Algorithm**
A specialized algorithm designed to efficiently find **all maximal cliques** in **sparse graphs**.

---
## 📂 Input Format

The input file should contain edges of the graph in the following format:
```
u v  # An edge between vertex u and vertex v
```
Each line should contain two integers **u** and **v**, representing an **undirected edge** between vertices **u** and **v**.

---
## 📊 Output

Each algorithm outputs the following:
- **Histogram** of clique sizes
- **Vertices** in the maximum clique
- **Execution time**
