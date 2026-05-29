# PHEMICS

**PHEMICS** — **Parallel Hot-plasma EMIC Solver for Multi-ion Anisotropic Plasmas**  
Pronounced: **FEE-mix**

PHEMICS is a MATLAB solver for the linear kinetic dispersion relation of **parallel-propagating EMIC waves** in a multi-ion plasma. It calculates the complex normalized wave frequency

\[
x = \frac{\omega}{\Omega_p}
  = \frac{\omega_r}{\Omega_p}
  + i\frac{\gamma}{\Omega_p},
\]

where \(\omega_r\) is the real wave frequency and \(\gamma\) is the linear growth or damping rate. A positive imaginary part indicates wave growth, while a negative imaginary part indicates damping.

The code supports both cold ion components and finite-temperature anisotropic ion components. Finite-temperature ions are treated kinetically through the plasma dispersion function.

---

## Features

- Solves the parallel-propagating left-hand EMIC dispersion relation.
- Supports multi-ion plasmas, including H\(^+\), He\(^+\), O\(^+\), or user-defined ion species.
- Allows both cold-fluid ion components and finite-temperature kinetic ion components.
- Includes temperature anisotropy through

  \[
  A = \frac{T_\perp}{T_\parallel} - 1.
  \]

- Tracks a dispersion branch over a user-defined normalized wave-number grid.
- Uses a rational multi-pole approximation to evaluate the plasma dispersion function.

---

## File structure

The main files are:

```text
PHEMICS.m
solveBranch.m
realImagEquation.m
dispersionRelation.m
plasmaDispersionFunction.m
```

Their calling sequence is:

```text
PHEMICS
  └── solveBranch
        └── realImagEquation
              └── dispersionRelation
                    └── plasmaDispersionFunction
```

### File descriptions

| File | Description |
|---|---|
| `PHEMICS.m` | Main user interface. Converts input ion parameters into the solver parameter structure and calls `solveBranch`. |
| `solveBranch.m` | Tracks one complex dispersion branch over a wave-number grid using `fsolve`. |
| `realImagEquation.m` | Converts the complex dispersion equation into two real equations. |
| `dispersionRelation.m` | Evaluates the parallel L-mode kinetic dispersion relation. |
| `plasmaDispersionFunction.m` | Evaluates the plasma dispersion function using a rational multi-pole approximation. |

---

## Requirements

- MATLAB
- Optimization Toolbox, because the solver uses `fsolve`

All `.m` files should be placed in the same directory, or the directory should be added to the MATLAB path.

---

## Quick start

```matlab
clear; clc;

% Background plasma parameters
B0 = 31200/5.5^3 * 1e-9;   % background magnetic field, T
ne = 50e6;                 % electron density, m^-3

% Normalized wave number
% K = k c / (omega_pe * sqrt(eps)), where eps = me/mp
K = logspace(-2, 1, 2000);

% Cold ion populations
% dd_cold_ion = [ns/ne;
%                qs/qe;
%                ms/mp]
dd_cold_ion = [
    0.8   0.1
    1     1
    1     4
];

% Warm/hot ion populations
% dd_warm_ion = [ns/ne;
%                qs/qe;
%                ms/mp;
%                Tperp/Tpara - 1;
%                Tpara_keV]
dd_warm_ion = [
    0.05   0.05
    1      1
    1      4
    1.5    0
    10     0.01
];

% Initial guesses for different branches
xstart_H  = 0.7 + 1e-4i;   % H-band initial guess
xstart_He = 0.1 + 1e-4i;   % He-band initial guess

% Solve branches
[K, wH]  = PHEMICS(K, xstart_H,  B0, ne, dd_cold_ion, dd_warm_ion);
[K, wHe] = PHEMICS(K, xstart_He, B0, ne, dd_cold_ion, dd_warm_ion);
```

---

## Plot example

```matlab
figure('Color','w')

subplot(2,1,1)
plot(K, real(wH),  'r--', 'LineWidth', 1.5); hold on
plot(K, real(wHe), 'r-',  'LineWidth', 1.5);
set(gca, 'XScale', 'log')
xlabel('K = k c / (\omega_{pe}\epsilon^{1/2})')
ylabel('\omega_r / \Omega_p')
legend('H band', 'He band')
grid on

subplot(2,1,2)
plot(K, imag(wH),  'r--', 'LineWidth', 1.5); hold on
plot(K, imag(wHe), 'r-',  'LineWidth', 1.5);
yline(0, 'k-')
set(gca, 'XScale', 'log')
xlabel('K = k c / (\omega_{pe}\epsilon^{1/2})')
ylabel('\gamma / \Omega_p')
legend('H band', 'He band')
grid on
```

---

## Input format

### `K`

Normalized wave number:

\[
K = \frac{k c}{\omega_{pe}\sqrt{\epsilon}},
\qquad
\epsilon = \frac{m_e}{m_p}.
\]

### `fstart`

Initial guess of the complex normalized frequency:

\[
f_\mathrm{start} = \frac{\omega}{\Omega_p}
= \frac{\omega_r}{\Omega_p} + i\frac{\gamma}{\Omega_p}.
\]

Although the input name is `fstart` in the current code, it is dimensionless and represents the normalized complex angular frequency, not frequency in Hz.

### `B0`

Background magnetic field in tesla.

### `ne`

Electron number density in m\(^{-3}\).  
For reference:

\[
1~\mathrm{cm}^{-3} = 10^6~\mathrm{m}^{-3}.
\]

### `dd_cold_ion`

Cold ion parameter matrix:

```matlab
dd_cold_ion = [ns/ne;
               qs/qe;
               ms/mp];
```

Each column represents one cold ion species. For example, cold H\(^+\) and cold He\(^+\):

```matlab
dd_cold_ion = [
    0.8   0.1
    1     1
    1     4
];
```

This input can be empty:

```matlab
dd_cold_ion = [];
```

### `dd_warm_ion`

Warm/hot ion parameter matrix:

```matlab
dd_warm_ion = [ns/ne;
               qs/qe;
               ms/mp;
               Tperp/Tpara - 1;
               Tpara_keV];
```

Each column represents one finite-temperature ion species. This input cannot be empty.

The first warm ion population is used as the reference species for the kinetic normalization. In typical EMIC applications, this should be the hot H\(^+\) population.

---

## Output

```matlab
[K, w] = PHEMICS(K, fstart, B0, ne, dd_cold_ion, dd_warm_ion);
```

- `K`: normalized wave-number array.
- `w`: complex normalized frequency array.

The real and imaginary parts are:

```matlab
omega_r_norm = real(w);   % omega_r / Omega_p
gamma_norm   = imag(w);   % gamma / Omega_p
```

---

## Notes on branch tracking

The dispersion equation may have multiple roots. The branch obtained by the solver depends on the initial guess `fstart`.

Typical initial guesses are:

```matlab
xstart_H  = 0.5  + 1e-4i;   % H band
xstart_H  = 0.7  + 1e-4i;   % H band, alternative
xstart_He = 0.1  + 1e-4i;   % He band
```

If many `NaN` values appear in the output, try:

- using a better initial guess,
- refining the `K` grid,
- scanning `K` in the opposite direction,
- avoiding points too close to ion cyclotron resonances,
- adjusting the `fsolve` options in `solveBranch.m`.

---

## Physical assumptions

PHEMICS currently assumes:

- parallel propagation,
- left-hand polarized EMIC waves,
- cold electrons,
- cold-fluid response for cold ion species,
- kinetic response for finite-temperature ion species,
- bi-Maxwellian temperature anisotropy for warm/hot ions.

Charge neutrality should approximately satisfy:

\[
\sum_j Z_j \frac{n_j}{n_e} \approx 1,
\]

where the sum includes all cold and warm ion species.

---

## Numerical cautions

- Very small `K` can be difficult because the plasma dispersion function argument can become large.
- Roots close to ion cyclotron frequencies can be numerically delicate.
- Branch tracking is sensitive to the initial guess.
- A species with zero parallel temperature should be treated as a cold ion, not as a warm/hot ion.
- The rational approximation used for the plasma dispersion function should be tested against a high-accuracy reference if very high precision is required.

---

## Recommended repository layout

```text
PHEMICS/
├── README.md
├── PHEMICS.m
├── solveBranch.m
├── realImagEquation.m
├── dispersionRelation.m
├── plasmaDispersionFunction.m
└── examples/
    └── example_basic.m
```

---

## Author

**Xue Z.X.**

Code development assistance: ChatGPT, OpenAI.  
Final implementation, verification, and scientific interpretation by Xue Z.X.

---

## License

Please add a license file before public release. For an open-source GitHub repository, common choices include MIT, BSD-3-Clause, GPL-3.0, or Apache-2.0.

---

## Citation

If you use PHEMICS in scientific work, please cite the relevant paper, thesis, or repository release associated with this code.

A suggested software citation format is:

```text
Xue, Z. X. (2026). PHEMICS: Parallel Hot-plasma EMIC Solver for Multi-ion Anisotropic Plasmas. GitHub repository.
```
