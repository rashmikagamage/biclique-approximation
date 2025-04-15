

# P-Q Biclique Counting

  

This is the implementation of the (p-q)-biclique approximation algorithm.

  

---

  

## 🛠 Build Instructions

  

1. Ensure you have `cmake` and `make` installed.

2. From the project root directory, run:

  

```bash

cmake  .

make

```

  

---

  

## 🚀 How to Run

  
```
bash
  /bin/run -f {datafile} -zstar3 -p 4 -q 5 -e 0.01 -d 0.01
```


### Parameters 
-  `-f <file>`: Path to the input graph file. **(Required)**

-  `-p <int>` and `-q <int>`: Size of the left and right parts of the biclique.  **(Required)**

-  `-e <float>`: Error parameter for approximation (default: **0.05**).

-  `-d <float>`: Optional parameter for failure probability (default: **0.05**, can be set inside the code).

```
bash

./biclique_counter  -f  graph.txt  -zstar3  -p  4  -q  5

./biclique_counter  -f  graph.txt  -zstar4  -p  4  -q  4  -e 0.01 -d -0.05

./biclique_counter  -f  graph.txt  -zstar5  -p  4  -q  4  -e  0.01
```
  

---

  

## 📥 Input Format

  

The input graph should be in the following format:

  

```

n m e

u1 v1

u2 v2

...

```

  

Where:

-  `n`: Number of vertices in the **left** partition (U)

-  `m`: Number of vertices in the **right** partition (V)

-  `e`: Total number of edges between U and V

- Each line after the first describes one undirected edge from U to V

  

### 📌 Example Input Graph


```
5 5 6
0 0
0 1
1 2
1 4
3 4
4 5
```

  

This describes a ** bipartite graph** with sets U = {0, 1, 2, 3, 4} and V = {0, 1, 2, 3, 4}, with 6 edges between them.

 ---

