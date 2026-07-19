# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-07-19

### Added
- Formal academic paper (`riemann.tex` / `riemann.pdf`) in AMS-LaTeX (`amsart`) format
- Student companion guide (`riemann_student_guide.tex` / `riemann_student_guide.pdf`) with 8 lab exercises and selected answers
- Pure-Python computational module (`riemann_codes.py`) with zero external dependencies:
  - Riemann Zeta function via Dirichlet eta analytic continuation
  - Left/Right/Midpoint Riemann integration
  - Prime counting π(x) and logarithmic integral Li(x) via midpoint Riemann sum
  - Spherical geodesic distances using Riemannian metric tensor
- Unit test suite (`test_riemann.py`) — 6 tests covering all four mathematical domains
- README with build/test instructions
- MIT License
- GitHub Actions CI workflow (compile + test)

### Fixed (Student Guide)
- Replaced markdown `**bold**` with LaTeX `\textbf{}` throughout
- Removed six `---` horizontal-rule artifacts that rendered as em-dashes
- Corrected Topic 1 Exercise 1: "magnitude |ζ| zero crossing" → "real part sign crossing" (matches implementation)
- Rewrote Topic 3 Exercise 2 (Chebyshev bound): now asks student to explain why finite-x ratio can exceed asymptotic bound
- Cleaned up draft sentence "? Wait..." in Topic 3 Answer 2
- Adjusted audience claim from "high school students" / "high-school algebra" to "advanced high-school / early undergraduate" with complex-numbers + trigonometry prerequisites
- Fixed "Cryptography & prime distribution" → "Prime distribution & primality testing" in roadmap table
- Added REPL import example to Getting Started section
- Added accuracy caveat for 5000-term eta series near critical line (±0.02 in t, |ζ| ∼ 0.003 at zeros)
- Title formatting: proper `\subtitle` macro instead of `\\` inside `\title{}`

### Fixed (Formal Paper)
- Title formatting: added `\subtitle` macro for cleaner AMS-LaTeX compliance