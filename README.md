
# Cuckoo Search Optimization with Lévy Flight & Adam Optimizer for the Salomon Problem

[![University](https://img.shields.io/badge/University-University%20of%20Ioannina-blue.svg)](https://www.uoi.gr)
[![Department](https://img.shields.io/badge/Department-Computer%20Science%20%26%20Telecommunications-lightgrey.svg)](https://www.dit.uoi.gr/)
[![Framework](https://img.shields.io/badge/Framework-XOPTIMUS-orange.svg)](https://github.com/itsoulos/GlobalOptimus)
[![Language](https://img.shields.io/badge/Language-C++-blueviolet.svg)]()

A robust C++ implementation of the **Cuckoo Search (CS)** metaheuristic algorithm enhanced with **Lévy Flights** and hybridized with the **Adam Optimizer**, applied to benchmark and minimize the multi-modal **Salomon Problem**. Developed as a semester project for the optimization course at the Department of Computer Science and Telecommunications, University of Ioannina.

---

## 📚 Table of Contents
1. [Overview](#-overview)
2. [Theoretical Background](#-theoretical-background)
   - [Cuckoo Search Algorithm](#cuckoo-search-algorithm)
   - [Lévy Flight & Random Walk](#lévy-flight--random-walk)
   - [Salomon Benchmark Problem](#salomon-benchmark-problem)
   - [Adam Optimizer Hybridization](#adam-optimizer-hybridization)
3. [Project Structure & Code Architecture](#-project-structure--code-architecture)
4. [Experimental Results & Analysis](#-experimental-results--analysis)
5. [Getting Started & Installation](#-getting-started--installation)
6. [Author & Credits](#-author--credits)

---

## 🔍 Overview

Optimization problems frequently suffer from complex error surfaces characterized by numerous local minima, flat regions, and sharp discontinuities. This project implements a sophisticated hybrid framework within the **XOPTIMUS** environment:
* **Global Exploration:** Uses **Cuckoo Search** combined with **Lévy Flights** to escape local traps and thoroughly scan multi-dimensional spaces.
* **Local Exploitation:** Integrates a modified top-$k$ differential local search and a final **Adam Optimizer** pass to precisely locate the floor of local or global basins.

---

## 📐 Theoretical Background

### Cuckoo Search Algorithm
Developed by Xin-She Yang and Suash Deb in 2009, Cuckoo Search mimics the obligate brood parasitism of certain cuckoo species that lay their eggs in the nests of other host birds. 
* **Rule 1:** Each cuckoo lays one egg at a time, deposited in a randomly chosen nest.
* **Rule 2:** The best nests with high-quality eggs (solutions) carry over to the next generations.
* **Rule 3:** The number of available host nests is fixed, and a host bird can discover an alien egg with probability $p_a \in [0, 1]$, subsequently abandoning the nest or building a new one.

### Lévy Flight & Random Walk
To answer *“where does the cuckoo fly?”*, the algorithm utilizes a random walk whose step lengths follow a Lévy probability distribution with heavy tails:
$$x_i^{t+1} = x_i^t + \alpha \otimes L(\lambda)$$
Where $\alpha$ is the step size scaling factor, and the step size $s$ is computed using Mantegna's algorithm:
$$s = \frac{U}{|V|^{\frac{1}{\lambda}}}$$
with $U \sim N(0, \sigma_u^2)$ and $V \sim N(0, \sigma_v^2)$.

### Salomon Benchmark Problem
The Salomon function is a multi-modal mathematical benchmark test defined over a hypercube domain (typically bounded within $[-100, 100]$ for dimensions $d = 4$):
$$f(x) = 1 - \cos\left(2\pi\sqrt{\sum_{i=1}^d x_i^2}\right) + 0.1\left(\sqrt{\sum_{i=1}^d x_i^2}\right)$$
* **Global Minimum:** $f(x) = 0$ at $x_i = 0$.
* **Challenge:** The function features a global minimum at the center surrounded by concentric ripples/rings acting as deceptive local traps.

### Adam Optimizer Hybridization
Once the Cuckoo Search phase identifies a promising region (ring), an **Adam Optimizer** instance is instantiated. Leveraging momentum and RMSprop adaptively:
* It smooths out noisy gradients caused by the trigonometric components of the Salomon function.
* It scales step sizes independently per variable, allowing it to rapidly traverse flat regions and brake over steep gradients to sink directly to the floor of the local ring.

---

## 🛠️ Project Structure & Code Architecture

The implementation relies on four core C++ files built on top of the XOPTIMUS framework classes (`Optimizer` and `Problem`):

| File | Description |
| :--- | :--- |
| **`Usermethod.h`** | Header file declaring class attributes (`var_num_nests`, `var_Pa`, `var_top_k`, nesting vectors, etc.) and core method interfaces (`init`, `step`, `levy_search`, `nest_abandonment`, `local_search`, `done`). |
| **`Usermethod.cpp`** | Implements the main Cuckoo Search lifecycle, Lévy flight generation, top-$k$ differential mutation local search, abandonment logic, and the post-processing Adam optimization hook. |
| **`Userproblem.h`** | Header file extending the framework's `Problem` class, declaring the objective function and gradient computations. |
| **`Userproblem.cpp`** | Defines the Salomon benchmark function calculation and analytical gradient derivatives (`gradient`) utilizing the chain rule for the Adam optimizer. |

---

## 📊 Experimental Results & Observations

* **Ring Traps:** Due to the concentric wave structure of the Salomon function, running smaller populations (e.g., 20 nests) often traps the global search in intermediate rings, yielding recurring fitness values such as $\approx 1.099$, $\approx 2.799$, or $\approx 3.199$ (corresponding to the exact depth of those rings).
* **Hybrid Success:** The Adam optimizer perfectly drives the solutions down to the exact floor of whatever ring the Cuckoo Search detects.
* **Scalability:** Increasing the maximum iterations (e.g., from 30 to 100) or population size allows the algorithm to bypass intermediate rings and approach the true global minimum at zero.

---

## 🚀 Getting Started & Execution

1. **Prerequisites:** Ensure you have the **XOPTIMUS** C++ optimization framework installed and configured in your development environment.
2. **Setup:** Place `Usermethod.h`, `Usermethod.cpp`, `Userproblem.h`, and `Userproblem.cpp` into your working directory.
3. **Execution:** Load the problem (`userproblem`) and method (`UserMethod`) through the XOPTIMUS interface, configure parameters (`popsize`, `max_iters`, `Pa`, etc.) as desired, and run.

---

## 📝 Author & Academic Use

* **Student:** Maria Alexandra Ndontou (ΜΑΡΙΑ ΑΛΕΞΑΝΔΡΑ ΝΤΟΝΤΟΥ)
* **Student ID:** 2799
* **Institution:** University of Ioannina — Department of Computer Science and Telecommunications
* **Course:** Optimization (Βελτιστοποίηση) — Semester Project

*This document is provided for academic reference and educational purposes.*
