"""
Riemann Visualization Module
Generates publication-quality plots for Riemann's discoveries.

Requires: matplotlib, numpy
Install: pip install matplotlib numpy
"""

import math
try:
    import numpy as np
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    np = None
    plt = None

from riemann_codes import (
    riemann_zeta, find_zeta_zero_on_critical_line, 
    riemann_sum, prime_counting_pi, logarithmic_integral_li,
    spherical_distance, christoffel_symbols, riemann_curvature_tensor,
    sectional_curvature, solve_geodesic, sphere_metric, 
    hyperbolic_metric, torus_metric,
    find_zeros_on_critical_line, explicit_formula_pi
)
from zeta_zeros import ZETA_ZEROS_T


def plot_zeta_critical_line(save_path=None):
    """Plot |ζ(1/2 + it)| on the critical line."""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available. Install with: pip install matplotlib numpy")
        return
    
    # Scan a range
    t_vals = np.linspace(0, 50, 1000)
    zeta_vals = []
    for t in t_vals:
        z = riemann_zeta(complex(0.5, t), terms=2000)
        zeta_vals.append(abs(z))
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Magnitude
    ax1.plot(t_vals, zeta_vals, 'b-', linewidth=0.8, label=r'$|\zeta(1/2 + it)|$')
    ax1.set_ylabel(r'$|\zeta(1/2 + it)|$')
    ax1.set_title('Riemann Zeta Function on the Critical Line')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_yscale('log')
    
    # Real and imaginary parts
    real_vals = []
    imag_vals = []
    for t in t_vals:
        z = riemann_zeta(complex(0.5, t), terms=2000)
        real_vals.append(z.real)
        imag_vals.append(z.imag)
    
    ax2.plot(t_vals, real_vals, 'r-', linewidth=0.8, label=r'$\Re(\zeta)$')
    ax2.plot(t_vals, imag_vals, 'g-', linewidth=0.8, label=r'$\Im(\zeta)$')
    ax2.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax2.set_xlabel('t (imaginary height)')
    ax2.set_ylabel(r'$\Re(\zeta), \Im(\zeta)$')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Mark known zeros
    for t_zero in ZETA_ZEROS_T[:10]:
        ax1.axvline(x=t_zero, color='orange', linestyle='--', alpha=0.5, linewidth=0.8)
        ax2.axvline(x=t_zero, color='orange', linestyle='--', alpha=0.5, linewidth=0.8)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()


def plot_prime_staircase(save_path=None):
    """Plot π(x) vs Li(x) and the explicit formula corrections."""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available.")
        return
    
    x_vals = list(range(2, 201))
    pi_vals = [prime_counting_pi(x) for x in x_vals]
    li_vals = [logarithmic_integral_li(x) for x in x_vals]
    
    # Explicit formula with 1, 3, 6 zeros
    zeros_6 = ZETA_ZEROS_T[:6]
    explicit_1 = [explicit_formula_pi(x, zeros_6[:1]) for x in x_vals]
    explicit_3 = [explicit_formula_pi(x, zeros_6[:3]) for x in x_vals]
    explicit_6 = [explicit_formula_pi(x, zeros_6[:6]) for x in x_vals]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    # Prime staircase vs Li(x)
    ax1.step(x_vals, pi_vals, 'k-', where='post', linewidth=1.5, label=r'$\pi(x)$ (exact)')
    ax1.plot(x_vals, li_vals, 'b--', linewidth=1.5, label=r'$\operatorname{Li}(x)$')
    ax1.set_ylabel('Count')
    ax1.set_title('Prime Counting Function π(x) vs Approximations')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0, 50)
    
    # Explicit formula corrections
    ax2.step(x_vals, pi_vals, 'k-', where='post', linewidth=1.5, label=r'$\pi(x)$ (exact)')
    ax2.plot(x_vals, explicit_1, 'r-', linewidth=1.2, label='Explicit (1 zero)')
    ax2.plot(x_vals, explicit_3, 'g-', linewidth=1.2, label='Explicit (3 zeros)')
    ax2.plot(x_vals, explicit_6, 'b-', linewidth=1.2, label='Explicit (6 zeros)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('Count')
    ax2.set_title('Riemann Explicit Formula: Convergence to π(x)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_ylim(-10, 50)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()


def plot_geodesics_on_manifolds(save_path=None):
    """Visualize geodesics on sphere, hyperbolic plane, and torus."""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available.")
        return
    
    fig = plt.figure(figsize=(15, 5))
    
    # 1. Sphere geodesic (meridian)
    ax1 = fig.add_subplot(131, projection='3d')
    init_coords = [math.pi/2, 0.0]
    init_vel = [-0.5, 0.0]
    traj = solve_geodesic(sphere_metric, init_coords, init_vel, (0, 6), dt=0.02)
    
    # Convert to 3D
    xs, ys, zs = [], [], []
    for _, coords, _ in traj:
        theta, phi = coords[0], coords[1]
        xs.append(math.sin(theta) * math.cos(phi))
        ys.append(math.sin(theta) * math.sin(phi))
        zs.append(math.cos(theta))
    
    # Sphere wireframe
    u = np.linspace(0, 2*np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    u, v = np.meshgrid(u, v)
    x = np.sin(v) * np.cos(u)
    y = np.sin(v) * np.sin(u)
    z = np.cos(v)
    ax1.plot_surface(x, y, z, alpha=0.1, rstride=2, cstride=2)
    ax1.plot(xs, ys, zs, 'r-', linewidth=2, label='Geodesic (meridian)')
    ax1.set_title('Sphere: Meridian Geodesic')
    ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_zlabel('z')
    
    # 2. Hyperbolic plane (Poincaré half-plane) - vertical geodesic
    ax2 = fig.add_subplot(132)
    init_coords_h = [0.0, 1.0]
    init_vel_h = [0.0, -0.2]
    traj_h = solve_geodesic(hyperbolic_metric, init_coords_h, init_vel_h, (0, 5), dt=0.02)
    xs_h = [p[1][0] for p in traj_h]
    ys_h = [p[1][1] for p in traj_h]
    ax2.plot(xs_h, ys_h, 'b-', linewidth=2, label='Geodesic (vertical)')
    
    # Show boundary
    ax2.axhline(y=0, color='k', linewidth=2)
    ax2.fill_between([-2, 2], 0, 5, alpha=0.1, color='blue')
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(0, 1.5)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y > 0')
    ax2.set_title('Hyperbolic Plane: Vertical Geodesic')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    # 3. Torus geodesic
    ax3 = fig.add_subplot(133, projection='3d')
    init_coords_t = [0.0, 0.0]
    init_vel_t = [0.3, 0.5]
    traj_t = solve_geodesic(torus_metric, init_coords_t, init_vel_t, (0, 20), dt=0.02)
    
    # Torus embedding: x = (R + r*cos(v))*cos(u), etc.
    R, r = 2.0, 1.0
    u_vals = np.linspace(0, 2*np.pi, 30)
    v_vals = np.linspace(0, 2*np.pi, 20)
    u_vals, v_vals = np.meshgrid(u_vals, v_vals)
    x_t = (R + r*np.cos(v_vals)) * np.cos(u_vals)
    y_t = (R + r*np.cos(v_vals)) * np.sin(u_vals)
    z_t = r * np.sin(v_vals)
    ax3.plot_surface(x_t, y_t, z_t, alpha=0.1, rstride=2, cstride=2)
    
    # Geodesic path
    xs_t, ys_t, zs_t = [], [], []
    for _, coords, _ in traj_t:
        u, v = coords[0], coords[1]
        xs_t.append((R + r*math.cos(v)) * math.cos(u))
        ys_t.append((R + r*math.cos(v)) * math.sin(u))
        zs_t.append(r * math.sin(v))
    ax3.plot(xs_t, ys_t, zs_t, 'g-', linewidth=2, label='Geodesic')
    ax3.set_title('Torus: Geodesic')
    ax3.set_xlabel('x'); ax3.set_ylabel('y'); ax3.set_zlabel('z')
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()


def plot_sectional_curvature(save_path=None):
    """Plot sectional curvature on different manifolds."""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available.")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # 1. Sphere: constant K=1
    ax = axes[0]
    thetas = np.linspace(0.1, math.pi-0.1, 100)
    K_sphere = [sectional_curvature(sphere_metric, [t, 0], [1,0], [0,1]) for t in thetas]
    ax.plot(thetas, K_sphere, 'b-', linewidth=2)
    ax.axhline(y=1, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel(r'$\theta$ (colatitude)')
    ax.set_ylabel('Sectional Curvature K')
    ax.set_title('Sphere: K = 1 (constant)')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.5, 1.5)
    
    # 2. Hyperbolic: constant K=-1
    ax = axes[1]
    ys = np.linspace(0.2, 3, 100)
    K_hyper = [sectional_curvature(hyperbolic_metric, [0, y], [1,0], [0,1]) for y in ys]
    ax.plot(ys, K_hyper, 'r-', linewidth=2)
    ax.axhline(y=-1, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('y')
    ax.set_ylabel('Sectional Curvature K')
    ax.set_title('Hyperbolic Plane: K = -1 (constant)')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1.5, -0.5)
    
    # 3. Torus: variable curvature
    ax = axes[2]
    vs = np.linspace(0, 2*math.pi, 200)
    K_torus = [sectional_curvature(torus_metric, [0, v], [1,0], [0,1]) for v in vs]
    ax.plot(vs, K_torus, 'g-', linewidth=2)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('v (minor angle)')
    ax.set_ylabel('Sectional Curvature K')
    ax.set_title('Torus: Variable Curvature')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()


def plot_riemann_sums(save_path=None):
    """Plot Riemann sum convergence for f(x) = x^2 on [0,1]."""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available.")
        return
    
    f = lambda x: x * x
    exact = 1.0 / 3.0
    ns = np.arange(1, 101)
    left = [riemann_sum(f, 0.0, 1.0, n, "left") for n in ns]
    right = [riemann_sum(f, 0.0, 1.0, n, "right") for n in ns]
    mid = [riemann_sum(f, 0.0, 1.0, n, "midpoint") for n in ns]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Convergence plot
    ax1.plot(ns, left, 'r-', label='Left', linewidth=1.5)
    ax1.plot(ns, right, 'b-', label='Right', linewidth=1.5)
    ax1.plot(ns, mid, 'g-', label='Midpoint', linewidth=1.5)
    ax1.axhline(y=exact, color='k', linestyle='--', alpha=0.5, label='Exact (1/3)')
    ax1.set_ylabel('Integral Approximation')
    ax1.set_title('Riemann Sums for f(x) = x² on [0,1]')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0.2, 0.45)
    
    # Error plot (log-log)
    err_left = [abs(l - exact) for l in left]
    err_right = [abs(r - exact) for r in right]
    err_mid = [abs(m - exact) for m in mid]
    
    ax2.loglog(ns, err_left, 'r-', label='Left Error', linewidth=1.5)
    ax2.loglog(ns, err_right, 'b-', label='Right Error', linewidth=1.5)
    ax2.loglog(ns, err_mid, 'g-', label='Midpoint Error', linewidth=1.5)
    ax2.loglog(ns, 1.0/ns, 'k--', alpha=0.5, label='O(1/n)')
    ax2.loglog(ns, 1.0/ns**2, 'k:', alpha=0.5, label='O(1/n²)')
    ax2.set_xlabel('n (subdivisions)')
    ax2.set_ylabel('Absolute Error')
    ax2.set_title('Convergence Rate (log-log)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()


def plot_spherical_geodesics(save_path=None):
    """Plot spherical geodesics - London to NY great circle and spherical triangle."""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib not available.")
        return
    
    fig = plt.figure(figsize=(12, 5))
    
    # 1. London to NY great circle
    ax1 = fig.add_subplot(121, projection='3d')
    # London and NY coordinates
    pt_london = (math.radians(90.0 - 51.5074), math.radians(-0.1278))
    pt_ny = (math.radians(90.0 - 40.7128), math.radians(-74.0060))
    
    # Geodesic path using solve_geodesic
    # For great circle, we can just interpolate
    def great_circle_path(p1, p2, steps=100):
        """Return points along great circle."""
        theta1, phi1 = p1
        theta2, phi2 = p2
        
        # Convert to 3D vectors
        v1 = np.array([math.sin(theta1)*math.cos(phi1), math.sin(theta1)*math.sin(phi1), math.cos(theta1)])
        v2 = np.array([math.sin(theta2)*math.cos(phi2), math.sin(theta2)*math.sin(phi2), math.cos(theta2)])
        
        # Spherical linear interpolation (slerp)
        dot = np.dot(v1, v2)
        dot = max(-1.0, min(1.0, dot))
        omega = math.acos(dot)
        sin_omega = math.sin(omega)
        
        points = []
        for i in range(steps + 1):
            t = i / steps
            if sin_omega > 1e-10:
                v = math.sin((1-t)*omega)/sin_omega * v1 + math.sin(t*omega)/sin_omega * v2
            else:
                v = v1
            points.append(v)
        return points
    
    gc_points = great_circle_path(pt_london, pt_ny, 50)
    
    # Sphere wireframe
    u = np.linspace(0, 2*np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    u, v = np.meshgrid(u, v)
    x = np.sin(v) * np.cos(u)
    y = np.sin(v) * np.sin(u)
    z = np.cos(v)
    ax1.plot_surface(x, y, z, alpha=0.1, rstride=2, cstride=2)
    
    # Great circle path
    xs = [p[0] for p in gc_points]
    ys = [p[1] for p in gc_points]
    zs = [p[2] for p in gc_points]
    ax1.plot(xs, ys, zs, 'r-', linewidth=2, label='London → NY')
    
    # Mark London and NY
    london_3d = np.array([math.sin(pt_london[0])*math.cos(pt_london[1]), 
                          math.sin(pt_london[0])*math.sin(pt_london[1]), math.cos(pt_london[0])])
    ny_3d = np.array([math.sin(pt_ny[0])*math.cos(pt_ny[1]), 
                      math.sin(pt_ny[0])*math.sin(pt_ny[1]), math.cos(pt_ny[0])])
    ax1.scatter([london_3d[0]], [london_3d[1]], [london_3d[2]], color='red', s=50, label='London')
    ax1.scatter([ny_3d[0]], [ny_3d[1]], [ny_3d[2]], color='blue', s=50, label='New York')
    
    ax1.set_title('Great Circle: London → New York (5570 km)')
    ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_zlabel('z')
    ax1.legend()
    
    # 2. Spherical triangle with 3 right angles
    ax2 = fig.add_subplot(122, projection='3d')
    
    # Triangle vertices: North Pole, Equator at 0°, Equator at 90°
    v1 = np.array([0, 0, 1])  # North pole
    v2 = np.array([1, 0, 0])  # Equator at prime meridian
    v3 = np.array([0, 1, 0])  # Equator at 90°E
    
    # Sphere wireframe
    ax2.plot_surface(x, y, z, alpha=0.1, rstride=2, cstride=2)
    
    # Triangle edges (great circles)
    def slerp_points(a, b, n=50):
        dot = np.dot(a, b)
        dot = max(-1.0, min(1.0, dot))
        omega = math.acos(dot)
        sin_omega = math.sin(omega)
        pts = []
        for i in range(n + 1):
            t = i / n
            if sin_omega > 1e-10:
                v = math.sin((1-t)*omega)/sin_omega * a + math.sin(t*omega)/sin_omega * b
            else:
                v = a
            pts.append(v)
        return pts
    
    edge1 = slerp_points(v1, v2)
    edge2 = slerp_points(v2, v3)
    edge3 = slerp_points(v3, v1)
    
    for edge, color in zip([edge1, edge2, edge3], ['red', 'blue', 'green']):
        xs = [p[0] for p in edge]
        ys = [p[1] for p in edge]
        zs = [p[2] for p in edge]
        ax2.plot(xs, ys, zs, color=color, linewidth=2)
    
    # Mark vertices
    ax2.scatter([v1[0], v2[0], v3[0]], [v1[1], v2[1], v3[1]], [v1[2], v2[2], v3[2]], 
                color='black', s=60)
    
    ax2.set_title('Spherical Triangle: 3 Right Angles (Sum = 270°)')
    ax2.set_xlabel('x'); ax2.set_ylabel('y'); ax2.set_zlabel('z')
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()


def generate_all_plots(output_dir="."):
    """Generate all visualization plots."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("Generating zeta critical line plot...")
    plot_zeta_critical_line(os.path.join(output_dir, "zeta_critical_line.png"))
    
    print("Generating prime staircase plot...")
    plot_prime_staircase(os.path.join(output_dir, "prime_staircase.png"))
    
    print("Generating geodesics on manifolds...")
    plot_geodesics_on_manifolds(os.path.join(output_dir, "geodesics_manifolds.png"))
    
    print("Generating sectional curvature plot...")
    plot_sectional_curvature(os.path.join(output_dir, "sectional_curvature.png"))
    
    print("Generating Riemann sums convergence plot...")
    plot_riemann_sums(os.path.join(output_dir, "riemann_sums.png"))
    
    print("Generating spherical geodesics plot...")
    plot_spherical_geodesics(os.path.join(output_dir, "spherical_geodesics.png"))
    
    print(f"All plots saved to {output_dir}/")


if __name__ == "__main__":
    if MATPLOTLIB_AVAILABLE:
        generate_all_plots("plots")
    else:
        print("Matplotlib not installed. Run: pip install matplotlib numpy")