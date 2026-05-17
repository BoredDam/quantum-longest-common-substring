# A Quantum Longest Common Substring Implementation
This repository is related to an upcoming paper for the [QUASAR26](https://sites.google.com/view/quasar26) conference, 
which is about an actual implementation of the Quantum Longest Common Substring algorithm, from the paper

> *Domenico Cantone, Simone Faro, Arianna Pavone, Caterina Viola,
  Quantum algorithms for longest common and palindromic substrings in the circuit model*.
>
> Abstract:
> 
> The Longest Common Substring (LCS) and Longest Palindromic Substring (LPS) problems are fundamental challenges in string processing,
  traditionally solved in linear time using classical computation through suffix trees. Recent breakthroughs by Le Gall and Seddighin
  [1] introduced sublinear quantum query algorithms, while Akmal and Jin [2] further improved the LCS complexity to $\tilde{\mathcal{O}}(n2/3)$.
  While these results are remarkable in the quantum query model, their practical implementation on real quantum hardware remains elusive.
  In this paper, we bridge this gap by presenting the first $\tilde{\mathcal{O}}(n)$ quantum algorithms for both LCS and LPS in the circuit model of computation.
  Our circuits are explicitly constructed and analyzed in terms of size and depth, achieving polylogarithmic overheads while preserving the $\tilde{\mathcal{O}}(n)$ depth bound.
  This provides, for the first time, concrete circuit-level blueprints and resource estimates for quantum solutions to LCS and LPS.
  Keywords: Quantum computing; Text processing; Sequence analysis

The code in this project provides an accurate implementation of the classical-quantum hybrid algorithm of the paper, 
while also applying some optimizations and some strategies (like a BBHT-inspired approach to the Grover search and 
2 auxillary qubits to optimize the search on large solution spaces).

## Folder structure

```
quantum-longest-common-substring/
├── .gitignore
├── LICENSE
├── README.md                      you are here! :)
├── data/
│   ├── 16char-2bit.csv            benchmarks on 16 chars strings - 2 bit encoding
│   ├── 8char-2bit.csv             benchmarks on 8 chars strings - 2 bit encoding
│   └── analysis.ipynb             
├── playground/
│   └── dna.ipynb                  quantum-LCS in a bioinformatics <3 inspired scenario
└── src/
    ├── QLCS.py                    the actual code!
    ├── QuantumLCS.ipynb
    └── for-testing.ipynb
```

