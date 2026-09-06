# SALT — Steady-state After Last Transients

`SALT` solves the periodic steady state of a transmission network directly,
using physics-based EMT device models for a single-harmonic, balanced
system. It solves for the complex bus voltages, the internal device states,
and the **system frequency**, which deviates from 60 Hz nominal whenever
there is a power mismatch.

`EMT` integrates the same device models in the time domain. Running EMT
from the SALT solution and watching it stay put is the verification that
the steady state is correct.

## Quick start

```matlab
run('<...>/SALT/scripts/addpathSALT.m')   % or addpath(genpath('<...>/SALT'))
run_SALT                                  % solve, verify, plot, report
```

Set `caseName` at the top of `run_SALT.m` to any name in the case table
just below it.

MATPOWER must be on the path for `psse2mpc`, `runpf`, `loadcase` and
`define_constants`; `addpathSALT` adds it if it is not already there.
`addpathSALT` also checks that no SALT function is shadowed by another
folder and errors if one is — several file names are shared with the legacy
solver directory, and MATLAB gives the current folder priority over the
path.

## Layout

```
run_SALT.m           entry point: case selection and solver settings
solver/              SALT.m  steady-state solver (parser + stamping)
                     EMT.m   time-domain solver
models/indices/      variable layouts, exposed as globals
models/ss/           steady-state stamp kernels + triplet helpers
models/td/           time-domain stamp kernels
models/init/         operating-point initialization, SALT -> EMT transfer
models/network/      Park transforms and three-phase helpers
scripts/             simulate / post-process / path setup
data/                network case files and device workbooks, with workbook
                     format designed from consulting ParaEmt [1],[2]
utils/               numeric and plotting helpers
docs/                device_models.tex - the analytical model equations
results/             generated plots and <case>_results.txt
```

## Cases

| Case | Network | Workbook | Load scale | Gen outage | Devices |
|---|---|---|---|---|---|
| `4bus`            | `4bus.raw`                    | `4bus_SALT.xlsx`                     | 0.75 | -    | 1 GENROU, 1 GFM IBR, 2 composite loads |
| `39bus`           | `case39_modified.m`           | `39bus_SALT.xlsx`                    | 0.97 | -    | 3 GENROU, 7 GFM IBR, 21 composite loads |
| `39bus_highPenetrationIBR` | `case39_modified.m`  | `39bus_highPenetrationIBR_SALT.xlsx` | 0.98 | bus 34 | 1 GENROU, 9 GFM IBR, 21 composite loads |
| `case118`         | `case118_modified.m`          | `case118_SALT.xlsx`                  | 1    | bus 6 | 27 GENROU, 27 GFM IBR, 99 composite loads |
| `case300`         | `case300_modified.m`          | `case300_SALT.xlsx`                  | 1    | bus 20 | 35 GENROU, 34 GFM IBR, 199 composite loads |
| `case1354pegase`  | `case1354pegase_modified.m`   | `case1354pegase_SALT.xlsx`           | 1    | bus 221 | 130 GENROU, 130 GFM IBR, 673 composite loads |
| `case3375wp`      | `case3375wp_modified.m`       | `case3375wp_SALT.xlsx`               | 1    | bus 10080 | 196 GENROU, 196 GFM IBR, 2424 composite loads |
| `case2kbus`       | `case_ACTIVSg2000_modified.m` | `case2kbus_SALT.xlsx`                | 1    | bus 1009 | 295 GENROU, 97 GFM IBR, 1125 composite loads |
| `case13659pegase` | `case13659pegase_modified.m`  | `case13659pegase_SALT.xlsx`          | 1    | bus 18 | 2046 GENROU, 2046 GFM IBR, 5499 composite loads |

The load scale DIVIDES the load power, so below 1 is a heavier load and a
lower system frequency. Every large case loses its THIRD generating unit,
which is what the timing driver `script_largeSystemTiming.m` does; the
device counts above are before that outage.

The two 39-bus rows are the same network at different inverter penetration,
matching the two legacy drivers: `run_39bus_SS_vsPF.m` keeps machines at
buses 30, 31 and 39, and `run_39bus_SS_concatenate.m` keeps only bus 30.

The six large systems are the cases in the timing driver
`script_largeSystemTiming.m`, brought over so the two codebases describe the
same systems. Their networks carry that driver's adaptations — transformer
ratios forced to 1 and phase shifts to 0, consistent with how ParaEMT models
a transformer [1] — and their workbooks carry its device parameters: every
inverter on a 200 MVA / 16.5 kV base, the WECC composite-load split, and
Kundur's small industrial motor.

Which generators become inverters follows the same rule as the driver: the
wind and solar machines where the case carries a `genfuel` field
(ACTIVSg2000), and the first half of the generator list otherwise.

An EMT verification run costs about an hour per two periods at 2000 buses,
so the case table sets `verify` false for the six large cases; they solve
and report the operating point without the time-domain check.

## Input files

Two files describe a system.

**Network file** — a PSS/E `.raw` or any MATPOWER-loadable `.m`/`.mat`
case. Supplies the buses, transmission lines, transformers, loads and the
generator dispatch.

**Device workbook** — one sheet per dynamic model. Column 1 holds the
parameter names, each further column is one device. A `comments` row
records where the parameter set came from and is ignored by the parser. The
sheet format is designed after consulting ParaEmt [1],[2]

| Sheet | `type` | Model |
|---|---|---|
| `gen`             | `GENROU`      | round-rotor synchronous machine |
| `exc`             | `EXDC1`       | DC exciter |
| `gov`             | `TGOV1`       | governor / turbine |
| `gfm ibr generic` | `GFM_GENERIC` | droop-controlled grid-forming inverter |
| `zip load`        | `ZIP`         | static ZIP load |
| `im load kundur`  | `IM_KUNDUR`   | induction motor |

The constructor checks that every sheet's `type` row names a model SALT
implements, and that the machines and inverters between them account for
exactly the generators in the power-flow case.

## Constructing a solver

```matlab
salt = SALT(networkFile, deviceFile, ...
    maxDeltaT   = 1e-5, ...   % sample grid for the derived period
    loadScale   = 0.5,  ...   % <1 = heavier load = lower frequency
    lineOutage  = false, ...  % or indices into the transmission-line list
    genOutage   = false, ...  % or bus numbers
    loadOutage  = false, ...  % or bus numbers
    rerunPfInit = true);
```

All outage options default to `false`. The constructor parses both files,
applies the outages, sets every device parameter through the normal setters
(`setGenrouParams`, `setExcDc1Params`, …), computes the operating point and
builds the stamp maps.

### The `rerunPfInit` flag

SALT does not take generator setpoints or controller references as inputs.
It back-calculates them from a **base case**: given each bus voltage and
each device's power injection, it inverts the device models to recover the
field voltage, mechanical torque, exciter reference and governor reference
that would hold the machine at that operating point.

`rerunPfInit` chooses what that base case is.

**`rerunPfInit = true`** — the power flow is solved first (MATPOWER by
default) and its solution becomes the base case. Because the setpoints are
derived from a consistent operating point, the system has no power
mismatch, and SALT reproduces the power-flow solution exactly at nominal
60 Hz. This is the natural starting point: it confirms the model chain is
self-consistent, and it gives contingency studies a clean pre-disturbance
reference.

**`rerunPfInit = false`** — the bus voltages and injections recorded in the
network file are used as they stand. Nothing guarantees they satisfy the
power-flow equations, so the residual mismatch is real, and the frequency
settles wherever generation and load actually balance — away from 60 Hz.

## Loads and bus shunts

A bus counts as a load if it draws real **or** reactive power, and each load
bus is then routed by what it physically is:

- **Real power present** → composite load: an induction motor plus a ZIP.
  How its states are initialised, and where the motor power factor does and
  does not enter, is documented at the top of
  `models/init/compLoadValInitCalculations_assumeW0.m`.
- **Reactive power only** → constant impedance. The composite load sizes its
  motor on real power, so a zero-P bus would give the motor a zero MVA base;
  such a bus is a shunt reactor or capacitor and is modelled as one.
- **Fixed shunts** (`Gs`, `Bs` in the bus table) → the same
  constant-impedance path.

Unlike a lossless element in series between two buses, a lossless shunt
reactor has the finite admittance `1/(jwL)` and needs no MNA branch current,
so both capacitors and reactors stamp directly.

Shunts are not a detail at scale: neither the 4-bus nor the 39-bus case has
any, but ACTIVSg2000 carries 18.3 GVAr across 149 buses — 96% of its load
reactive demand, so ignoring them moves bus voltages by up to 0.23 pu and
the system frequency by 0.04 Hz — and case13659pegase has 8799
constant-impedance elements against 5499 composite loads.

The row layout is phase-major: all bus and MNA real/imaginary rows for phase
A, then phase B, then phase C, then the device blocks, then the single
system-frequency row. `getRowIdxAbc_fromBusNums` returns indices in that
order, so a gather over several buses runs bus-fastest within each phase.

## How stamping works

Both solvers use the same contract:

- `buildStampMaps` runs **once** and fills `stampRow` / `stampCol` — the
  row and column of every Jacobian entry — via `saltStampTriplets`.
- Each Newton iteration refills `stampVal` only. There is one vectorized
  call per device *model*, not per device: states are gathered as
  `[nState × 1 × nDevice]` pages and the Jacobian blocks come out of
  `pagemtimes` for the whole fleet at once.
- `saltAssembleY` turns the three arrays into the sparse matrix. Duplicate
  `(row,col)` pairs accumulate, which is what the network stamps rely on.

Newton's method is stamped in the form
`J(x[k])·x[k+1] = J(x[k])·x[k] − F(x[k])`, so a linear equation contributes
its coefficients to `J` and its constant to the right-hand side directly.

Model variable layouts live in globals, so no index vector is recomputed
inside the Newton loop. There are two sets, populated separately by the
class that owns them: `saltInitIndices` fills the steady-state layout
(`idxGenrou_Idq`, …) from SALT's constructor, and `emtInitIndices` fills the
time-domain layout (`idxLoadCompTD_Iabc`, …) from EMT's. A steady-state
solve that never builds an EMT twin leaves the TD set untouched.

The analytical equations behind each stamp kernel are in
`docs/device_models.tex`.

## Post-processing

`scriptPostprocess_SALT` builds the EMT twin with `salt.buildEmt`, runs it
period by period until consecutive periods line up in voltage, current and
speed, then plots each requested bus with the SALT waveform dashed on top
of the EMT waveform.

The command window gets the system frequency; the bus table is one line per
bus — 13659 of them on the largest case — so it is written to
`results/<case>_results.txt`:

```
================ SALT steady-state solution ================
Case             : 4bus
System frequency : 59.239469 Hz  (w = 372.212560 rad/s)
Deviation from 60 Hz : -0.760531 Hz

Bus      |V| [PU]         angle(V) [deg]
---      --------         --------------
1        1.002894          +0.0000
2        0.904885         -12.1482
3        0.877921         -20.1393
4        0.924537         -17.4942
============================================================
```

Plots are written to `results/` as `<networkFileName>_bus<N>.png`. Set the
plot list for a case in the case table in `run_SALT.m`.

# Artificial intelligence (AI) disclaimer

The use of AI (Anthropic, Opus 5. 2026, August 28 - September 1) was applied in
helping with documentation generation, code beautification, and code
conversion from MATLAB to Python.

# References
[1] M. Xiong, B. Wang, D. Vaidhynathan, J. Maack, M. Reynolds, A. Hoke, 
    K. Sun, J. Tan, "ParaEMT: an open source, parallelizable, and 
    HPC-compatible EMT simulator for large-scale IBR-rich power grids," 
    IEEE Trans. Power Del., vol. 39, no. 2, pp. 911-921, Apr. 2024.
[2] M. Xiong, B. Wang, D. Vaidhynathan, J. Maack, M. Reynolds, A. Hoke, K.
    Sun, D. Ramasubramanian, V. Verma, J. Tan, "An open-source parallel EMT 
    simulation framework," Electric Power Syst. Res., vol. 235, 2024, Art. 
    no. 110734.