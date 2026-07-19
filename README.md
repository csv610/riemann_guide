# Computational Analysis of Riemann's Discoveries

This repository contains a self-contained, publication-ready mathematical and computational package that explores and verifies several major discoveries of **Bernhard Riemann (1826–1866)**.

The project combines historical mathematical theory with modern software engineering practices, providing a formal academic paper, a student-friendly laboratory guide, and an independent Python verification suite with zero external dependencies.

---

## 📂 Project Directory Structure

*   [riemann.tex](file:///Users/csv610/Projects/Riemann/riemann.tex) / [riemann.pdf](file:///Users/csv610/Projects/Riemann/riemann.pdf) — The formal academic paper, structured in standard AMS-LaTeX (`amsart`) format.
*   [riemann_student_guide.tex](file:///Users/csv610/Projects/Riemann/riemann_student_guide.tex) / [riemann_student_guide.pdf](file:///Users/csv610/Projects/Riemann/riemann_student_guide.pdf) — An accessible, educational companion guide for high school and undergraduate students, featuring 8 computational lab exercises and selected answers.
*   [riemann_codes.py](file:///Users/csv610/Projects/Riemann/riemann_codes.py) — The core Python module containing all simulations, analytic continuations, Riemann integration methods, prime counting, and spherical geodesic calculations.
*   [test_riemann.py](file:///Users/csv610/Projects/Riemann/test_riemann.py) — A 6-test unit suite verifying the mathematical accuracy of all Python implementations.

---

## 🧬 Topics Verified & Explored

1.  **The Riemann Zeta Function \& Riemann Hypothesis:** Evaluating $\zeta(s)$ on the complex critical line $s = 1/2 + it$ via Dirichlet eta analytic continuation, and locating the first non-trivial zero crossing at $t \approx 14.1347$.
2.  **Riemann Integration (Riemann Sums):** Approximating definite integrals using Left, Right, and Midpoint Riemann sums, and verifying the quadratic $O(1/n^2)$ error convergence rate of the midpoint method.
3.  **Prime Counting \& the Logarithmic Integral $\operatorname{Li}(x)$:** Demonstrating the Prime Number Theorem by comparing exact prime counts $\pi(x)$ with Riemann's logarithmic integral approximation.
4.  **Riemannian Geometry \& Spherical Geodesics:** Computing shortest paths (great-circle distances) on curved 2D spheres using a spherical Riemannian metric tensor.

---

## 🚀 Getting Started

### Prerequisites
You only need a standard installation of **Python 3.x** and a **LaTeX distribution** (like TeX Live, MacTeX, or MiKTeX) to compile the documents. No third-party Python packages are required.

### 1. Run the Computational Demos
To run the full suite of simulations and print out the numerical verifications:
```bash
python3 riemann_codes.py
```

### 2. Run the Unit Test Suite
To verify the math functions against standard reference values:
```bash
python3 -m unittest test_riemann.py
```

### 3. Compile the Documents
To compile the formal paper and the student guide to PDFs:
```bash
pdflatex riemann.tex
pdflatex riemann_student_guide.tex
```
*(Note: It is recommended to compile twice to resolve cross-references correctly.)*

---

## ✍️ Author
*   **Chaman Singh Verma** — Independent Researcher
*   Email: [csv610@gmail.com](mailto:csv610@gmail.com)
