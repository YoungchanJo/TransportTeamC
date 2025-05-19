#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numbers>

// --------------------------------------------------
// Fixed-Source 1D Slab Transport with Consistent DSA
// --------------------------------------------------
// Implements discrete ordinates transport (Diamond-Difference)
// with fixed source Q(x) and consistent DSA acceleration.
// Based on Termproject.pdf (Group C requirement), DOM - DD.pdf, and DSA.pdf.
// --------------------------------------------------

// Physical and numerical parameters (Termproject.pdf specifies L, N, M)
constexpr double L = 10.0;    // slab length [cm]
constexpr int    N = 200;     // number of spatial cells (mesh resolution)
constexpr int    M = 16;      // number of discrete angles (Gauss-Legendre points)
constexpr double tol = 1e-8;  // convergence tolerance for source iteration + DSA
constexpr int    maxIter = 5000;    // maximum number of SI+DSA iterations
constexpr double dx = L / N;       // cell width Δx = L/N
constexpr int    N1 = N / 2;       // grid index separating Region1 and Region2
constexpr double PI = M_PI;        // π constant (used in quadrature)

// Region-specific cross-sections and sources (Termproject.pdf Table 1)
// Region 1: 0 ≤ x < L/2
constexpr double sigma_t1 = 1.0;   // total cross section [1/cm]
constexpr double sigma_s1 = 0.9;   // scattering cross section [1/cm]
constexpr double Q1 = 1.0;   // fixed source term
// Region 2: L/2 ≤ x ≤ L
constexpr double sigma_t2 = 0.5;
constexpr double sigma_s2 = 0.3;
constexpr double Q2 = 0.0;

// Step 1: Gauss–Legendre Quadrature Initialization
// (See DOM - DD.pdf Section 2.1: Gauss-Legendre Quadrature)
// Computes M-point roots of P_M(z) and weights w_i = 2/[(1-z_i^2)(P'_M(z_i))^2]
void initializeQuadrature(double mu[M], double w[M]) {
    int half = (M + 1) / 2;  // use symmetry: only compute half the roots
    for (int i = 0; i < half; ++i) {
        // Initial guess: Eq.(2.3) in DOM - DD.pdf
        double z = std::cos(PI * (i + 0.75) / (M + 0.5));
        double zPrev, P0, P1, dP;
        // Newton-Raphson: refine root of P_M(z) (Eq.2.6 in DOM - DD.pdf)
        do {
            P0 = 1.0;       // P_0(z)
            P1 = z;         // P_1(z)
            // Recurrence for P_j(z) up to j=M (Eq.2.5)
            for (int j = 2; j <= M; ++j) {
                double P2 = ((2 * j - 1) * z * P1 - (j - 1) * P0) / j;
                P0 = P1;
                P1 = P2;
            }
            // Derivative P'_M(z) (Eq.2.7)
            dP = M * (z * P1 - P0) / (z * z - 1);
            zPrev = z;
            z = zPrev - P1 / dP;
        } while (std::abs(z - zPrev) > 1e-14);

        // Assign symmetric directions μ and weights
        mu[i] = -z;
        mu[M - 1 - i] = z;
        double weight = 2.0 / ((1 - z * z) * dP * dP);
        w[i] = weight;
        w[M - 1 - i] = weight;
    }
}

// Step 2: Diamond-Difference transport sweep
// (See DOM - DD.pdf Section 3.2: Diamond-Difference Discretization)
// Solves discrete ordinates transport: μ∂ψ/∂x + Σ_t ψ = Σ_s φ/2 + Q/2
// using diamond-difference: ψ_i = (ψ_{i-1/2}+ψ_{i+1/2})/2
void transportSweep(
    const double phi[N],              // current scalar flux φ
    const double sigma_s[N],          // scattering cross-sections
    const double sigma_t[N],          // total cross-sections
    const double mu[M],               // discrete directions
    const double w[M],                // quadrature weights
    const double Q[N],                // fixed source term
    double phi_half[N]                // output: half-step scalar flux
) {
    // Initialize φ^(half) to zero
    std::fill(phi_half, phi_half + N, 0.0);

    // Loop over all angles n = 1..M
    for (int n = 0; n < M; ++n) {
        double psi_edge = 0.0;        // boundary angular flux (vacuum BC, Eq.3.4)
        double mu_n = mu[n];
        double absMu = std::abs(mu_n);
        double w_n = w[n];

        if (mu_n > 0) {
            // sweep left-to-right (DOM - DD.pdf Fig.3.1)
            for (int i = 0; i < N; ++i) {
                // scattering source: Σ_s φ Δx / 2
                double scat_term = 0.5 * sigma_s[i] * dx * phi[i];
                // fixed source: Q Δx / 2
                double src_term = 0.5 * Q[i] * dx;
                // numerator and denominator from diamond-difference (Eq.3.7)
                double num = 2 * mu_n * psi_edge + scat_term + src_term;
                double den = 2 * mu_n + sigma_t[i] * dx;
                double psi = num / den;  // ψ_{i+1/2}
                phi_half[i] += w_n * psi; // accumulate weighted flux
                psi_edge = 2 * psi - psi_edge; // update for next cell (Eq.3.8)
            }
        }
        else {
            // sweep right-to-left for μ < 0
            for (int i = N - 1; i >= 0; --i) {
                double scat_term = 0.5 * sigma_s[i] * dx * phi[i];
                double src_term = 0.5 * Q[i] * dx;
                double num = 2 * absMu * psi_edge + scat_term + src_term;
                double den = 2 * absMu + sigma_t[i] * dx;
                double psi = num / den;
                phi_half[i] += w_n * psi;
                psi_edge = 2 * psi - psi_edge;
            }
        }
    }
}

// Step 3: Consistent DSA (Eq.3.28–3.31 in DSA.pdf)
// Builds low-order diffusion equation to accelerate convergence
void buildConsistentDSA(
    const double sigma_t[N],
    const double sigma_s[N],
    const double D[N],               // diffusion coefficient D=1/(3Σ_t)
    const double phi[N],              // old scalar flux
    const double phi_half[N],         // high-order flux
    double f_edge[N + 1],             // output: edge flux corrections
    double delta_phi[N]               // output: scalar flux correction
) {
    static double r[N];          // residual (Eq.3.28)
    static double A[N + 1], B[N + 1], C[N + 1], RHS[N + 1]; // tridiagonal matrix
    static double cp[N + 1], dp[N + 1]; // Thomas solver arrays

    // 1) Compute scattering-only residual r_j = 0.5 Σ_s (φ_half - φ) Δx
    for (int j = 0; j < N; ++j) {
        r[j] = 0.5 * sigma_s[j] * (phi_half[j] - phi[j]) * dx;
    }

    // 2) Apply Dirichlet BCs for f_edge at j=0 and j=N (Eq.3.29–3.30)
    A[0] = 0; B[0] = 1;     C[0] = 0;     RHS[0] = phi_half[0];
    A[N] = 0; B[N] = 1;     C[N] = 0;     RHS[N] = phi_half[N - 1];

    // 3) Fill interior tridiagonal coefficients (Eq.3.28)
    for (int j = 1; j < N; ++j) {
        A[j] = -D[j - 1] / dx;                         // lower diag
        B[j] = (D[j - 1] + D[j]) / dx + sigma_t[j] * dx;  // main diag
        C[j] = -D[j] / dx;                         // upper diag
        RHS[j] = r[j];                                 // RHS
    }

    // 4) Solve tridiagonal system via Thomas algorithm (DSA.pdf Sec.3.3)
    cp[0] = C[0] / B[0];
    dp[0] = RHS[0] / B[0];
    for (int i = 1; i <= N; ++i) {
        double m = B[i] - A[i] * cp[i - 1];
        cp[i] = C[i] / m;
        dp[i] = (RHS[i] - A[i] * dp[i - 1]) / m;
    }
    f_edge[N] = dp[N];
    for (int i = N - 1; i >= 0; --i) {
        f_edge[i] = dp[i] - cp[i] * f_edge[i + 1];
    }

    // 5) Compute scalar flux correction δφ_j = (f_{j-1/2}+f_{j+1/2})/2 (Eq.3.31)
    for (int j = 0; j < N; ++j) {
        delta_phi[j] = 0.5 * (f_edge[j] + f_edge[j + 1]);
    }
}

// Step 4: Update φ and check convergence criterion
// |φ^(k+1) - φ^(k)|_∞ < tol
// (Implements Eq.3.31 update and ∞-norm check)
double updateFlux(
    double phi[N],
    const double phi_half[N],
    const double delta_phi[N]
) {
    double maxErr = 0.0;
    for (int i = 0; i < N; ++i) {
        double newPhi = phi_half[i] + delta_phi[i];
        maxErr = std::max(maxErr, std::abs(newPhi - phi[i]));
        phi[i] = newPhi;
    }
    return maxErr;
}

int main() {
    // Allocate quadrature arrays and initialize (Step 1)
    double mu[M], w[M];
    initializeQuadrature(mu, w);

    // Allocate solution and coefficient arrays
    double phi[N], phi_half[N], sigma_t[N], sigma_s[N], D[N], Q[N];
    double f_edge[N + 1], delta_phi[N];

    // Initialize region properties and diffusion coefficients
    for (int i = 0; i < N; ++i) {
        phi[i] = 1.0;  // initial guess for scalar flux
        if (i < N1) {
            sigma_t[i] = sigma_t1; sigma_s[i] = sigma_s1; Q[i] = Q1;  // Region 1
        }
        else {
            sigma_t[i] = sigma_t2; sigma_s[i] = sigma_s2; Q[i] = Q2;  // Region 2
        }
        D[i] = 1.0 / (3.0 * sigma_t[i]);  // diffusion coefficient
    }

    // Main Source Iteration + DSA Acceleration Loop (Termproject.pdf workflow)
    for (int iter = 1; iter <= maxIter; ++iter) {
        transportSweep(phi, sigma_s, sigma_t, mu, w, Q, phi_half);      // high-order
        buildConsistentDSA(sigma_t, sigma_s, D, phi, phi_half, f_edge, delta_phi); // low-order
        double err = updateFlux(phi, phi_half, delta_phi);             // update & error
        if (err < tol) break;  // convergence reached
    }

    // Output converged flux at cell centers (x = (i+0.5)Δx)
    for (int i = 0; i < N; ++i) {
        double x = (i + 0.5) * dx;
        std::cout << x << " " << phi[i] << "\n";
    }
    return 0;
}
