"""
riemann_codes.py

A student-friendly Python exploration of Bernhard Riemann's mathematical work.

Topics included:
1. The Riemann Zeta Function & Critical Line Zeros
2. Riemann Integration (Left, Right, and Midpoint Riemann Sums)
3. Prime Counting Function & Logarithmic Integral Li(x)
4. Riemannian Geometry (Geodesics on a Spherical Metric)
5. General Riemannian Geometry (Curvature & Geodesics on Manifolds)
6. Riemann Explicit Formula (Primes from Zeta Zeros)

Pure Python, no third-party libraries required.
"""

import math
from zeta_zeros import ZETA_ZEROS_T

# ============================================================
# Special Functions
# ============================================================

def exponential_integral_Ei(z):
    """
    Compute the exponential integral Ei(z) for complex z.
    Uses series for small |z|, continued fraction for large |z|.
    Branch cut along negative real axis.
    """
    # For small |z|, use series: Ei(z) = γ + ln(z) + Σ_{k=1}∞ z^k/(k·k!)
    if abs(z) < 8:
        euler_gamma = 0.57721566490153286060651209008240243104215933593992
        log_z = complex(math.log(abs(z)), math.atan2(z.imag, z.real))
        result = complex(euler_gamma, 0) + log_z
        term = z
        for k in range(1, 100):
            result += term / k
            term *= z / (k + 1)
            if abs(term / k) < 1e-15:
                break
        return result
    
    # For large |z|, use continued fraction:
    # Ei(z) = e^z/z * (1 + 1/(z + 2/(z + 3/(z + ...))))
    # Using Lentz's algorithm
    tiny = 1e-30
    f = z
    c = f
    d = 0.0
    h = 1.0
    
    for n in range(1, 100):
        a = n
        b = z
        d = b + a * d
        if abs(d) < tiny:
            d = tiny
        c = b + a / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = c * d
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    
    exp_z = complex(math.cos(z.imag), math.sin(z.imag)) * math.exp(z.real)
    return exp_z * h / z


def complex_logarithmic_integral(x, rho, use_ei=True):
    """
    Compute Li(x^ρ) = Ei(ρ·ln(x)) for x > 1, complex ρ with Re(ρ) > 0.
    This is the proper way to compute the complex logarithmic integral.
    """
    if x <= 1:
        return complex(0.0, 0.0)
    
    log_x = math.log(x)
    z = complex(rho.real * log_x, rho.imag * log_x)
    
    if use_ei:
        # Li(x^ρ) = Ei(ρ·ln(x))
        return exponential_integral_Ei(z)
    else:
        # Fallback to direct integration (for verification)
        return complex_logarithmic_integral_direct(x, rho)


def complex_logarithmic_integral_direct(x, rho, steps=2000):
    """Direct complex line integration - kept for verification."""
    if x <= 1:
        return complex(0.0, 0.0)
    
    # Upper limit: x^ρ = exp(ρ * ln(x))
    upper = complex(math.cos(rho.imag * math.log(x)) * math.exp(rho.real * math.log(x)),
                    math.sin(rho.imag * math.log(x)) * math.exp(rho.real * math.log(x)))
    lower = complex(2.0, 0.0)
    
    # Straight-line integration in complex plane
    dt = (upper - lower) / steps
    result = complex(0.0, 0.0)
    for k in range(steps):
        t = lower + (k + 0.5) * dt  # midpoint
        if abs(t) > 1e-10:
            # 1/ln(t) in complex plane
            log_t = complex(math.log(abs(t)), math.atan2(t.imag, t.real))
            if abs(log_t) > 1e-10:
                result += dt / log_t
    return result


def find_zeros_on_critical_line(num_zeros, t_start=0.0, t_step=0.05, max_t=100.0):
    """
    Find the first num_zeros non-trivial zeros of ζ(s) on the critical line s = 1/2 + it.
    Returns list of t values where ζ(1/2 + it) = 0.
    """
    zeros = []
    t = t_start
    if t < 10:
        t = 10.0
    
    prev_t = t
    prev_val = riemann_zeta(complex(0.5, prev_t)).real
    
    while len(zeros) < num_zeros and t <= max_t:
        t += t_step
        val = riemann_zeta(complex(0.5, t)).real
        if prev_val * val < 0:
            zero_t = prev_t - prev_val * (t - prev_t) / (val - prev_val)
            zeros.append(zero_t)
        prev_t = t
        prev_val = val
    
    return zeros


def explicit_formula_pi(x, zeros, use_ei=True):
    """
    Compute Riemann's explicit formula for π(x):
    π(x) ≈ Li(x) - Σ_ρ Li(x^ρ) - log(2) + ∫_x^∞ dt/(t(t²-1)log t)
    
    For simplicity, we omit the integral term and -log(2).
    zeros is a list of t values where ρ = 1/2 + it.
    """
    if x <= 2:
        return 0
    
    # Main term: Li(x)
    main_term = logarithmic_integral_li(x)
    
    # Sum over zeros using Ei-based complex integration
    zero_sum = complex(0.0, 0.0)
    for t in zeros:
        rho = complex(0.5, t)
        li_rho = complex_logarithmic_integral(x, rho, use_ei=use_ei)
        li_conj = complex_logarithmic_integral(x, complex(0.5, -t), use_ei=use_ei)
        zero_sum += li_rho + li_conj
    
    result = main_term - zero_sum.real
    return result


def demo_explicit_formula():
    print("=" * 70)
    print("6. Riemann Explicit Formula (Primes from Zeta Zeros)")
    print("=" * 70)
    print("Riemann's explicit formula expresses π(x) in terms of zeta zeros:")
    print("  π(x) ≈ Li(x) - Σ_ρ Li(x^ρ) - log(2) + ∫_x^∞ dt/(t(t²-1)log t)")
    print("where the sum is over non-trivial zeros ρ = 1/2 + iγ.")
    print()
    print("Using precomputed high-precision zeros (Odlyzko tables, ~15 digits):")
    print()
    
    # Use first 6 high-precision zeros
    zeros = ZETA_ZEROS_T[:6]
    print(f"First 6 high-precision zeros: {[f'{z:.6f}' for z in zeros]}")
    print()
    
    print(f"{'x':>5} | {'π(x) exact':>10} | {'Li(x)':>10} | {'Explicit (1 zero)':>16} | {'Explicit (3 zeros)':>16} | {'Explicit (6 zeros)':>16}")
    print("-" * 85)
    for x in [10, 20, 50, 100, 200]:
        pi_exact = prime_counting_pi(x)
        li_x = logarithmic_integral_li(x)
        
        # Try with 1, 3, 6 zeros
        exp_1 = explicit_formula_pi(x, zeros[:1], use_ei=True)
        exp_3 = explicit_formula_pi(x, zeros[:3], use_ei=True)
        exp_6 = explicit_formula_pi(x, zeros[:6], use_ei=True)
        
        print(f"{x:5d} | {pi_exact:10d} | {li_x:10.4f} | {exp_1:16.4f} | {exp_3:16.4f} | {exp_6:16.4f}")
    
    print()
    print("The sum over zeros corrects Li(x) toward the exact π(x).")
    print("With high-precision zeros and Ei-based Li(x^ρ), the formula converges to π(x).")


# ============================================================
# 1. The Riemann Zeta Function & Critical Line Zeros
# ============================================================

def riemann_zeta(s, terms=5000):
    """
    Compute the Riemann Zeta function zeta(s) for a complex number s.
    Uses the Dirichlet Eta function (alternating zeta series) analytic continuation:
    eta(s) = sum_{n=1}^inf (-1)^{n-1} / n^s
    zeta(s) = eta(s) / (1 - 2^{1-s})
    This converges for Re(s) > 0, which covers the critical strip.
    s is represented as Python's native complex type.
    """
    if s.real <= 0:
        raise ValueError("This implementation only converges for Re(s) > 0.")
    
    # Calculate Dirichlet Eta sum
    eta_val = complex(0.0, 0.0)
    for n in range(1, terms + 1):
        # term = (-1)^{n-1} * n^{-s}
        # n^{-s} = exp(-s * ln(n))
        exponent = -s * math.log(n)
        term = complex(math.cos(exponent.imag), math.sin(exponent.imag)) * math.exp(exponent.real)
        if n % 2 == 0:
            eta_val -= term
        else:
            eta_val += term
           
    # Calculate divisor: 1 - 2^{1-s}
    two_pow = math.exp((1.0 - s.real) * math.log(2.0))
    two_pow_complex = two_pow * complex(math.cos(-s.imag * math.log(2.0)), math.sin(-s.imag * math.log(2.0)))
    divisor = 1.0 - two_pow_complex
    
    if abs(divisor) < 1e-12:
        return complex(float('nan'), float('nan'))  # Singularity at s = 1
        
    return eta_val / divisor


def find_zeta_zero_on_critical_line(t_start, t_end, step=0.01):
    """
    Scan a range of t values on the critical line s = 0.5 + it
    and find sign changes of the real part of zeta(s) to locate zeros.
    """
    zeros = []
    prev_t = t_start
    prev_val = riemann_zeta(complex(0.5, prev_t)).real
    
    t = t_start + step
    while t <= t_end:
        val = riemann_zeta(complex(0.5, t)).real
        # Check for zero crossing (sign change)
        if prev_val * val < 0:
            # Linear interpolation for zero crossing
            zero_t = prev_t - prev_val * (t - prev_t) / (val - prev_val)
            zeros.append(zero_t)
        prev_t = t
        prev_val = val
        t += step
    return zeros


def demo_zeta_function():
    print("=" * 70)
    print("1. Riemann Zeta Function & Critical Line Zeros")
    print("=" * 70)
    print("Evaluating Riemann Zeta on the critical line s = 0.5 + it:")
    # Evaluate at known non-zero points
    print(f"zeta(0.5 + 10i) = {riemann_zeta(complex(0.5, 10.0))}")
    
    # Scan for the first non-trivial zero (known to be near t = 14.1347)
    print("\nScanning for non-trivial zeros on s = 0.5 + it for t in [13.0, 16.0]...")
    zeros = find_zeta_zero_on_critical_line(13.0, 16.0)
    for z in zeros:
        print(f"-> Non-trivial zero located at: s = 0.5 + {z:.4f}i")
        print(f"   Value at zero: |zeta(s)| = {abs(riemann_zeta(complex(0.5, z))):.6f}")
    print()


# ============================================================
# 2. Riemann Integration (Riemann Sums)
# ============================================================

def riemann_sum(f, a, b, n, method="midpoint"):
    """
    Calculate the Riemann Sum approximation of the integral of f from a to b with n subdivisions.
    Methods supported: 'left', 'right', 'midpoint'.
    """
    dx = (b - a) / n
    total_sum = 0.0
    for i in range(n):
        if method == "left":
            x = a + i * dx
        elif method == "right":
            x = a + (i + 1) * dx
        elif method == "midpoint":
            x = a + (i + 0.5) * dx
        else:
            raise ValueError("Method must be 'left', 'right', or 'midpoint'")
        total_sum += f(x)
    return total_sum * dx


def demo_riemann_integration():
    print("=" * 70)
    print("2. Riemann Integration (Riemann Sums)")
    print("=" * 70)
    # Test function: f(x) = x^2, integral from 0 to 1 is 1/3 ~ 0.333333
    f = lambda x: x * x
    print("Integrating f(x) = x^2 from 0 to 1 (Exact Integral = 1/3):")
    for n in [10, 100, 1000]:
        left = riemann_sum(f, 0.0, 1.0, n, "left")
        right = riemann_sum(f, 0.0, 1.0, n, "right")
        mid = riemann_sum(f, 0.0, 1.0, n, "midpoint")
        print(f"Subdivisions n = {n:4d}:")
        print(f"  * Left Riemann Sum:     {left:.6f}  (Error: {abs(left - 1/3):.6f})")
        print(f"  * Right Riemann Sum:    {right:.6f}  (Error: {abs(right - 1/3):.6f})")
        print(f"  * Midpoint Riemann Sum: {mid:.6f}  (Error: {abs(mid - 1/3):.6f})")
    print()


# ============================================================
# 3. Prime Counting Function & Logarithmic Integral Li(x)
# ============================================================

def is_prime(n):
    """Trial division helper to check if n is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def prime_counting_pi(x):
    """Calculate the number of primes less than or equal to x."""
    count = 0
    for i in range(2, int(x) + 1):
        if is_prime(i):
            count += 1
    return count


def logarithmic_integral_li(x, n=5000):
    """
    Compute the Logarithmic Integral Li(x) = integral_2^x (1 / ln(t)) dt
    using Midpoint Riemann Sum integration.
    """
    if x < 2:
        raise ValueError("Li(x) is only computed for x >= 2 in this context.")
    f = lambda t: 1.0 / math.log(t)
    return riemann_sum(f, 2.0, x, n, "midpoint")


def demo_prime_distribution():
    print("=" * 70)
    print("3. Prime Counting Function & Logarithmic Integral")
    print("=" * 70)
    print("Comparing pi(x) (exact prime count) with Riemann's logarithmic integral Li(x):")
    print(f"{'x':>6} | {'pi(x)':>8} | {'Li(x)':>10} | {'Absolute Error':>14}")
    print("-" * 47)
    for x in [10, 50, 100, 500, 1000]:
        pi_val = prime_counting_pi(x)
        li_val = logarithmic_integral_li(x)
        print(f"{x:6d} | {pi_val:8d} | {li_val:10.4f} | {abs(pi_val - li_val):14.4f}")
    print()


# ============================================================
# 4. Riemannian Geometry (Geodesics on a Spherical Metric)
# ============================================================

def spherical_distance(pt1, pt2, R=6371.0):
    """
    Calculate the geodesic distance between two points on a sphere of radius R.
    Coordinates are represented as (theta, phi) in radians:
    theta: colatitude (0 to pi, measured from the North Pole)
    phi: longitude (0 to 2*pi)
    Uses the Riemannian metric ds^2 = R^2 d_theta^2 + R^2 sin^2(theta) d_phi^2.
    """
    theta1, phi1 = pt1
    theta2, phi2 = pt2
    
    # Boundary checks
    if not (0.0 <= theta1 <= math.pi) or not (0.0 <= theta2 <= math.pi):
        raise ValueError("Colatitude theta must be in range [0, pi] radians.")
       
    # Formula for geodesic distance on sphere:
    # d = R * arccos(cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2))
    cos_arg = math.cos(theta1)*math.cos(theta2) + math.sin(theta1)*math.sin(theta2)*math.cos(phi1 - phi2)
    # Clamp cos_arg to [-1.0, 1.0] to prevent floating point domain errors
    cos_arg = max(-1.0, min(1.0, cos_arg))
    return R * math.acos(cos_arg)


def demo_riemannian_geometry():
    print("=" * 70)
    print("4. Riemannian Geometry (Geodesics on a Sphere)")
    print("=" * 70)
    # Earth radius ~ 6371.0 km
    # Coordinate 1: London, UK (Latitude 51.5074 N, Longitude 0.1278 W)
    # Colatitude theta = 90 - latitude = 38.4926 degrees -> convert to radians
    # Longitude phi = -0.1278 degrees -> convert to radians
    pt_london = (math.radians(90.0 - 51.5074), math.radians(-0.1278))
    
    # Coordinate 2: New York, USA (Latitude 40.7128 N, Longitude 74.0060 W)
    # Colatitude theta = 90 - latitude = 49.2872 degrees -> convert to radians
    # Longitude phi = -74.0060 degrees -> convert to radians
    pt_ny = (math.radians(90.0 - 40.7128), math.radians(-74.0060))
    
    dist = spherical_distance(pt_london, pt_ny)
    print("Using Earth's spherical Riemannian metric tensor, calculate geodesic distance:")
    print("From: London, UK")
    print("To:   New York, USA")
    print(f"-> Geodesic Great-Circle Distance = {dist:.2f} km")
    print()


# ============================================================
# 5. General Riemannian Geometry (n-D Manifolds, Curvature, Geodesics)
# ============================================================

def christoffel_symbols(metric_func, coords, h=1e-5):
    """
    Compute Christoffel symbols of the second kind Gamma^i_{jk} for a Riemannian metric.
    
    Parameters:
    - metric_func: function taking coordinates (list/tuple) -> metric tensor g_ij (n x n matrix)
    - coords: list/tuple of coordinate values (length n)
    - h: finite difference step
    
    Returns:
    - Gamma: 3D list Gamma[i][j][k] = Gamma^i_{jk} (n x n x n)
    
    Formula: Gamma^i_{jk} = 0.5 * g^{il} * (dg_{lj}/dx^k + dg_{lk}/dx^j - dg_{jk}/dx^l)
    """
    n = len(coords)
    # Get metric and its inverse at the point
    g = metric_func(coords)
    # Compute inverse metric g^{ij}
    g_inv = invert_matrix(g)
    
    # Compute partial derivatives of metric components: dg_{ab}/dx^c
    dg = [[[0.0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                coords_plus = list(coords)
                coords_minus = list(coords)
                coords_plus[c] += h
                coords_minus[c] -= h
                g_plus = metric_func(coords_plus)
                g_minus = metric_func(coords_minus)
                dg[a][b][c] = (g_plus[a][b] - g_minus[a][b]) / (2 * h)
    
    # Compute Christoffel symbols
    Gamma = [[[0.0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                val = 0.0
                for l in range(n):
                    val += g_inv[i][l] * (dg[l][j][k] + dg[l][k][j] - dg[j][k][l])
                Gamma[i][j][k] = 0.5 * val
    return Gamma


def invert_matrix(m):
    """Invert a square matrix using Gaussian elimination (pure Python)."""
    n = len(m)
    # Augment with identity
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(m)]
    # Forward elimination
    for col in range(n):
        # Find pivot
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1e-12:
            raise ValueError("Matrix is singular")
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        # Normalize pivot row
        piv_val = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= piv_val
        # Eliminate other rows
        for r in range(n):
            if r != col:
                factor = aug[r][col]
                for j in range(2 * n):
                    aug[r][j] -= factor * aug[col][j]
    # Extract inverse
    return [row[n:] for row in aug]


def riemann_curvature_tensor(metric_func, coords, h=1e-5):
    """
    Compute Riemann curvature tensor R^i_{jkl} for a Riemannian metric.
    
    Formula: R^i_{jkl} = d_k Gamma^i_{jl} - d_l Gamma^i_{jk} 
                         + Gamma^i_{km} Gamma^m_{jl} - Gamma^i_{lm} Gamma^m_{jk}
    """
    n = len(coords)
    Gamma = christoffel_symbols(metric_func, coords, h)
    
    # Compute derivatives of Christoffel symbols
    dGamma = [[[[0.0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for l in range(n):
                for k in range(n):
                    coords_plus = list(coords)
                    coords_minus = list(coords)
                    coords_plus[k] += h
                    coords_minus[k] -= h
                    Gamma_plus = christoffel_symbols(metric_func, coords_plus, h)
                    Gamma_minus = christoffel_symbols(metric_func, coords_minus, h)
                    dGamma[i][j][l][k] = (Gamma_plus[i][j][l] - Gamma_minus[i][j][l]) / (2 * h)
    
    # Compute Riemann tensor
    R = [[[[0.0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    val = dGamma[i][j][l][k] - dGamma[i][j][k][l]
                    for m in range(n):
                        val += Gamma[i][k][m] * Gamma[m][j][l] - Gamma[i][l][m] * Gamma[m][j][k]
                    R[i][j][k][l] = val
    return R


def sectional_curvature(metric_func, coords, u, v, h=1e-5):
    """
    Compute sectional curvature K(u,v) for two tangent vectors u, v at coords.
    
    K(u,v) = R(u,v,u,v) / (g(u,u)*g(v,v) - g(u,v)^2)
    where R(u,v,u,v) = R_{ijkl} u^i v^j u^k v^l = g_{im} R^m_{jkl} u^i v^j u^k v^l
    """
    n = len(coords)
    R_mixed = riemann_curvature_tensor(metric_func, coords, h)  # R^m_{jkl}
    g = metric_func(coords)
    
    # Compute R(u,v,u,v) = g_{im} R^m_{jkl} u^i v^j u^k v^l
    Ru_v_u_v = 0.0
    for i in range(n):
        for m in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        Ru_v_u_v += g[i][m] * R_mixed[m][j][k][l] * u[i] * v[j] * u[k] * v[l]
    
    # Compute denominator: |u|^2 |v|^2 - <u,v>^2
    def norm_sq(vec):
        s = 0.0
        for i in range(n):
            for j in range(n):
                s += g[i][j] * vec[i] * vec[j]
        return s
    
    def inner(vec1, vec2):
        s = 0.0
        for i in range(n):
            for j in range(n):
                s += g[i][j] * vec1[i] * vec2[j]
        return s
    
    denom = norm_sq(u) * norm_sq(v) - inner(u, v)**2
    if abs(denom) < 1e-12:
        return float('nan')
    return Ru_v_u_v / denom


def geodesic_ode(metric_func, state, t, h=1e-5):
    """
    RHS of geodesic equation: d^2 x^i/dt^2 + Gamma^i_{jk} (dx^j/dt)(dx^k/dt) = 0
    
    state = [x^1, ..., x^n, v^1, ..., v^n] where v^i = dx^i/dt
    Returns derivative of state.
    """
    n = len(state) // 2
    coords = state[:n]
    vel = state[n:]
    Gamma = christoffel_symbols(metric_func, coords, h)
    
    # dv^i/dt = -Gamma^i_{jk} v^j v^k
    accel = [0.0] * n
    for i in range(n):
        for j in range(n):
            for k in range(n):
                accel[i] -= Gamma[i][j][k] * vel[j] * vel[k]
    
    return vel + accel


def rk4_step(f, state, t, dt, *args):
    """Single RK4 step."""
    k1 = f(state, t, *args)
    k2 = f([state[i] + 0.5 * dt * k1[i] for i in range(len(state))], t + 0.5 * dt, *args)
    k3 = f([state[i] + 0.5 * dt * k2[i] for i in range(len(state))], t + 0.5 * dt, *args)
    k4 = f([state[i] + dt * k3[i] for i in range(len(state))], t + dt, *args)
    return [state[i] + dt/6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(len(state))]


def solve_geodesic(metric_func, initial_coords, initial_vel, t_span, dt=0.01, h=1e-5):
    """
    Solve geodesic ODE using RK4.
    
    Parameters:
    - metric_func: function(coords) -> metric tensor
    - initial_coords: list of starting coordinates
    - initial_vel: list of initial velocities (dx^i/dt)
    - t_span: (t0, t1) integration interval
    - dt: time step
    - h: finite difference step for Christoffel symbols
    
    Returns:
    - trajectory: list of (t, coords, vel) tuples
    """
    t0, t1 = t_span
    state = initial_coords + initial_vel
    t = t0
    trajectory = [(t, state[:len(initial_coords)], state[len(initial_coords):])]
    
    def ode(state, t):
        return geodesic_ode(metric_func, state, t, h)
    
    while t < t1:
        if t + dt > t1:
            dt = t1 - t
        state = rk4_step(ode, state, t, dt)
        t += dt
        trajectory.append((t, state[:len(initial_coords)], state[len(initial_coords):]))
    return trajectory


def sphere_metric(coords, R=1.0):
    """
    Metric tensor for a 2-sphere of radius R.
    coords = [theta, phi] (colatitude, longitude)
    g = [[R^2, 0], [0, R^2 sin^2(theta)]]
    """
    theta, phi = coords
    return [[R*R, 0.0], [0.0, R*R * math.sin(theta)**2]]


def hyperbolic_metric(coords, R=1.0):
    """
    Metric tensor for 2D hyperbolic plane (Poincare half-plane model).
    coords = [x, y] with y > 0
    g = [[R^2/y^2, 0], [0, R^2/y^2]]
    """
    x, y = coords
    if y <= 0:
        raise ValueError("y must be > 0 for hyperbolic metric")
    factor = (R/y)**2
    return [[factor, 0.0], [0.0, factor]]


def torus_metric(coords, R=2.0, r=1.0):
    """
    Metric tensor for a 2-torus (embedding in R^3).
    coords = [u, v] where u,v in [0, 2*pi)
    Major radius R, minor radius r
    """
    u, v = coords
    # Embedding: x = (R + r*cos(v))*cos(u), y = (R + r*cos(v))*sin(u), z = r*sin(v)
    # Induced metric: g_uu = (R + r*cos(v))^2, g_vv = r^2, g_uv = 0
    guu = (R + r * math.cos(v))**2
    gvv = r * r
    return [[guu, 0.0], [0.0, gvv]]


def demo_general_riemannian_geometry():
    print("=" * 70)
    print("5. General Riemannian Geometry (Curvature & Geodesics on Manifolds)")
    print("=" * 70)
    
    # ---- Sphere ----
    print("\n--- 2-Sphere (R=1) ---")
    coords = [math.pi/2, 0.0]  # Equator, prime meridian
    g = sphere_metric(coords)
    print(f"Metric at equator (theta=pi/2, phi=0): {g}")
    
    Gamma = christoffel_symbols(sphere_metric, coords)
    print(f"Non-zero Christoffel symbols:")
    for i in range(2):
        for j in range(2):
            for k in range(2):
                if abs(Gamma[i][j][k]) > 1e-10:
                    print(f"  Gamma^{i}_{{{j}{k}}} = {Gamma[i][j][k]:.6f}")
    
    R = riemann_curvature_tensor(sphere_metric, coords)
    print(f"Riemann tensor component R^0_101 (should be 1 for unit sphere): {R[0][1][0][1]:.6f}")
    
    # Sectional curvature for sphere = 1/R^2 = 1
    u = [1.0, 0.0]  # theta direction
    v = [0.0, 1.0]  # phi direction
    K = sectional_curvature(sphere_metric, coords, u, v)
    print(f"Sectional curvature K(d/dtheta, d/dphi) = {K:.6f} (expected 1.0)")
    
    # ---- Geodesic on sphere: great circle from equator going "north" ----
    print("\n--- Geodesic on Sphere ---")
    init_coords = [math.pi/2, 0.0]  # Start at equator
    init_vel = [-0.1, 0.0]  # Move toward north pole (decreasing theta)
    traj = solve_geodesic(sphere_metric, init_coords, init_vel, (0, 10), dt=0.05)
    print(f"Start: theta={init_coords[0]:.4f}, phi={init_coords[1]:.4f}")
    print(f"Initial velocity: dtheta/dt={init_vel[0]:.4f}, dphi/dt={init_vel[1]:.4f}")
    print(f"After t=10: theta={traj[-1][1][0]:.6f}, phi={traj[-1][1][1]:.6f}")
    print(f"  (Should approach north pole: theta=0, phi=0)")
    
    # ---- Hyperbolic plane ----
    print("\n--- Hyperbolic Plane (Poincare half-plane) ---")
    coords_h = [0.0, 1.0]  # x=0, y=1
    g_h = hyperbolic_metric(coords_h)
    print(f"Metric at (x=0, y=1): {g_h}")
    K_h = sectional_curvature(hyperbolic_metric, coords_h, [1.0, 0.0], [0.0, 1.0])
    print(f"Sectional curvature = {K_h:.6f} (expected -1.0)")
    
    # ---- Torus ----
    print("\n--- 2-Torus (R=2, r=1) ---")
    coords_t = [0.0, 0.0]
    g_t = torus_metric(coords_t)
    print(f"Metric at (u=0, v=0): {g_t}")
    Gamma_t = christoffel_symbols(torus_metric, coords_t)
    print(f"Non-zero Christoffel symbols at (0,0):")
    for i in range(2):
        for j in range(2):
            for k in range(2):
                if abs(Gamma_t[i][j][k]) > 1e-10:
                    print(f"  Gamma^{i}_{{{j}{k}}} = {Gamma_t[i][j][k]:.6f}")


# ============================================================
# 6. Riemann Explicit Formula (Primes from Zeta Zeros)
# ============================================================
# [Already defined above - see special functions section]


# ============================================================
# Main Runner
# ============================================================

def run_all():
    demo_zeta_function()
    demo_riemann_integration()
    demo_prime_distribution()
    demo_riemannian_geometry()
    demo_general_riemannian_geometry()
    demo_explicit_formula()


if __name__ == "__main__":
    run_all()