# GHOST — Graph Hamiltonian Open-System Transport

GHOST propagates neutrino flavour states (two or three flavours, vacuum or
matter, with or without decoherence) by computing `|ψ(t)⟩ = exp(−iHt)|ψ(0)⟩`
with a Chebyshev expansion and a sparse matrix–vector product. The propagator
does not know it is working with neutrinos. It accepts any sparse Hermitian `H`,
and flavour oscillation is one small instance of that.

What it computes:

- two- and three-flavour oscillation in vacuum and in constant-density matter,
  including the CP phase;
- propagation through a varying matter profile: a power law `ρ ∝ r⁻ⁿ` (the
  supernova envelope) or a tabulated profile `ρ(r)`, `Y_e(r)`, with optional
  density noise;
- energy scans for neutrinos and antineutrinos, e.g. a DUNE-like spectrum of 200
  energies in one call;
- decoherence (pure dephasing in the flavour basis) through quantum
  trajectories;
- the magnetic sector: spin precession of a Dirac neutrino with a magnetic
  moment along a tabulated trajectory `ρ, Y_e, B(s)`, with an optional
  neutrino background and collision rate (Majorana neutrinos and antineutrinos
  are not supported yet);
- coherent flavour probabilities, mass-basis populations and the incoherent
  (phase-averaged) flavour probabilities. For astrophysical sources the
  incoherent ones are the observable.

GHOST also plugs into [SNEWPY](https://github.com/SNEWS2/snewpy) as a
`FlavorTransformation`, so the full SNEWPY → SNOwGLoBES chain can use a
propagated, time- and energy-dependent conversion instead of the adiabatic or
non-adiabatic limits.

> Code comments, the user-facing CLI messages and the Manual
> (`docs/apostila/apostila.pdf`) are in Portuguese.

## Repository layout

| path | contents |
|---|---|
| `core/` | the numerical core, C99 that also compiles as C++17 (prefix `dae_`) |
| `native/` | the `ghost` command-line program, the acceptance test suite, benchmarks |
| `wasm/` | the WebAssembly build of the same core |
| `src/`, `index.html`, `public/` | the browser app (React + Vite), running the core in a Web Worker |
| `templates/` | C++, Python and Wolfram Language templates used by the app's export |
| `tools/neutrinos/snewpy_ghost.py` | the SNEWPY integration |
| `tools/neutrinos/trajetoria_ns.py` | ray tables for the magnetic sector |
| `tools/` | core hash, amalgamation, smoke tests and CSV/dump comparison |
| `specs/` | example input files |
| `docs/apostila/` | the Manual (LaTeX source and PDF) |
| `CONVENTIONS.md` | the rules the core follows (in Portuguese) |

## Requirements

- a C99 compiler (GCC or Clang) and `make`;
- Node.js 20 or later, for the web app and the build tools;
- [Emscripten](https://emscripten.org/) (`emsdk`), only for the WebAssembly
  build;
- Python 3 with NumPy and SNEWPY, only for the SNEWPY integration;
- XeLaTeX, only to rebuild the Manual.

## Build and test

```sh
npm install
make -C native test          # builds the core and the CLI, runs the acceptance tests
make -C native cxx-check     # the core must also compile as C++17
make -C native asan          # AddressSanitizer + UndefinedBehaviorSanitizer

source ~/emsdk/emsdk_env.sh  # wherever emsdk is installed
npm run wasm                 # core → WebAssembly
npm test                     # web app unit tests (Vitest)
npm run build                # static site in dist/
npm run dev                  # local development server
```

`npm run ci` runs the core tests, the C++17 check, the type check, the unit
tests, the build and a smoke test in one go. The GitHub workflow also runs the
sanitizers and the full suite under WebAssembly.

## Running a calculation

A calculation is described by a `spec.json` file. The spec is parsed only by the
C core, so the CLI, the web app and the SNEWPY wrapper read it the same way.

```sh
./native/build/bin/ghost run specs/neutrino-msw.json --saida out.csv
./native/build/bin/ghost traj <spec.json> --saida rho.csv     # quantum trajectories
./native/build/bin/ghost --version
```

`ghost run` also accepts `--estado` (write the final state) and
`--realizacoes N` (number of realisations).

A two-flavour MSW example (`specs/neutrino-msw.json`):

```json
{
  "format_version": 3,
  "seed": 1,
  "neutrino": {
    "n_flavors": 2,
    "theta12_deg": 10.0,
    "dm21_sq_eV2": 7.42e-5,
    "energy_MeV": 20.0,
    "flavor": "e",
    "rho_g_cm3": 46.1, "Ye": 0.5
  },
  "time": { "t1": 2500.0, "nt": 500 },
  "observables": { "target": 1, "population": false },
  "realizations": 1
}
```

Main keys:

- `neutrino`: mixing angles and `delta_cp_deg`, mass splittings,
  `energy_MeV`, initial `flavor` (`e`, `mu`, `tau`), `antineutrino`, and a
  constant density `rho_g_cm3`, `Ye`.
- `matter`: a varying profile, with `profile` set to `constant`, `power_law`
  (`rho_ref_g_cm3`, `r_ref_km`, `exponent`, `r_start_km`) or `table` (`r_km`,
  `rho_table`, `ye_table`), plus `tolerance` and the optional `noise_*` keys.
- `scan`: an energy scan (`energy_min_MeV`, `energy_max_MeV`, `points`, `log`,
  `both_signs`).
- `magnetic` and `trajectory`: the magnetic sector (`nature`, `mu_nu_bohr`,
  `initial_helicity`, `integrator`, `coherence`, and the tabulated `s_km`,
  `rho_table`, `ye_table`, `b1_gauss`, `b2_gauss`). The tables along a ray
  through a static spherical star, in flat, Schwarzschild or tabulated
  metrics, can be generated with `tools/neutrinos/trajetoria_ns.py`.
- `trajectories`: quantum-trajectory settings (`n_traj`, `gamma_deph`,
  `rho_stride`, `output_mode`).

### Output

The output is a CSV file. Its header lines (`#!`) carry the program version, the
**core hash** and the full canonical spec, so every result records the code and
the input that produced it. For a neutrino run the columns are:

```
t, norm, ipr, coh_l1, p_target, p_e, p_mu, [p_tau,] m_1, m_2, [m_3,] inc_e, inc_mu, [inc_tau]
```

`p_*` are coherent flavour probabilities, `m_*` mass-basis populations, and
`inc_*` the incoherent flavour probabilities.

The core hash is a hash of the core sources (`npm run hash`), and the binary
reports it:

```
$ ghost --version
ghost 0.1.0 (nucleo <hash>)
```

### Example specs

| file | what it runs |
|---|---|
| `specs/neutrino-vacuo.json` | two flavours in vacuum, written directly as a two-node graph |
| `specs/neutrino-msw.json` | two flavours in constant-density matter (MSW resonance) |
| `specs/neutrino-dune.json` | three flavours with CP phase, DUNE baseline and density |
| `specs/neutrino-espectro-dune.json` | DUNE-like energy scan, neutrinos and antineutrinos |
| `specs/neutrino-supernova.json` | three flavours through a supernova envelope `ρ ∝ r⁻³` |
| `specs/magnetico/` | the magnetic sector: one crossing of `Y_e = 1/3` with gradient `g` = 0.01, 0.1 and 1 per km, to compare with Landau–Zener |
| `specs/microtubulo.json`, `specs/sonda/`, `specs/oraculo/` | graph quantum walks, and reference outputs used by the tests |

## SNEWPY integration

`tools/neutrinos/snewpy_ghost.py` defines `GhostMSW`, a SNEWPY
`FlavorTransformation` that runs the `ghost` binary through the density profile
of each time step and returns incoherent probabilities:

```python
from snewpy_ghost import GhostMSW, perfil_lei_potencia

xform = GhostMSW(perfil_lei_potencia(), energias_MeV=[5, 10, 20, 40],
                 tempos_s=[0.0])
xform.prob_ee(0.0, 20.0)      # P(νe → νe), propagated
```

The profile is a callable `t_s → (r_km, rho_g_cm3)`. A static profile can be
passed as the tuple itself. `perfil_choque` builds a power-law envelope with a
moving shock front. Build the CLI first (`make -C native`): the wrapper looks for
`native/build/bin/ghost` relative to its own location, or takes `binario=`.

## The web app

`npm run dev` (or the static build in `dist/`) opens a browser front end that
runs the same core, compiled to WebAssembly, in a Web Worker. It builds graphs
(path, cycle, grid, hypercube, complete, stochastic block model, microtubule
lattice), propagates the walk, plots the observables and exports the calculation
as a self-contained C++, Python or Wolfram Language program, carrying the spec
and the core hash. It can re-import a CSV produced elsewhere. The interface is
available in Portuguese, English, French and Italian.

## How it works

Two techniques let a graph propagator handle neutrino flavour:

- **Self-loops are the diagonal.** Mass terms and the matter potential enter as
  self-loops `[i, i, w]` in the edge list, and mixing enters as off-diagonal
  edges.
- **The real `2n × 2n` embedding.** The core uses only real arithmetic (no
  `_Complex`, so that it stays valid C++17). A complex Hermitian `H = X + iY`
  (X symmetric, Y antisymmetric) is propagated through the real symmetric
  matrix

  ```
  H̃ = [ X  −Y ]
      [ Y   X ]
  ```

  whose invariant subspace `{(ψ, −iψ)}` reproduces the dynamics of `H`. This
  holds for unitary evolution only. A collapse in the trajectory method leaves
  that subspace.

The spectral bound used by the Chebyshev expansion is guaranteed, not
estimated. The core has no mutable global state and is bit-for-bit
reproducible, and it is never built with `-ffast-math`. `CONVENTIONS.md` lists
these rules and the reasons for them.

## Licence

MIT. See [`LICENSE`](LICENSE). Algorithm credits and the licences of bundled
fonts and packages are in [`NOTICE`](NOTICE).
