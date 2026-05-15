# Building ngspice as a Shared Library (Windows)

## Prerequisites

- Visual Studio 2022 (Community or higher) with C++ workload
- CMake 3.20+ (optional, VS solution provided)

## Steps

### Option A: Visual Studio Solution (recommended)

1. Open `ref/ngspice/visualc/sharedspice.sln` in VS2022
2. Set configuration: Release | x64
3. In project properties → C\C++ → Preprocessor:
   - Ensure `SIMULATOR` and `HAS_WINGUI` are defined
4. Build → Build Solution
5. Output: `ref/ngspice/visualc/sharedspice/x64/Release/ngspice.dll`

### Option B: CMake

```bash
cd ref/ngspice
mkdir build && cd build
cmake .. -G "Visual Studio 17 2022" -A x64 \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON
cmake --build . --config Release
```

### Verifying Instrumentation

After building, verify both exports are present:

```bash
dumpbin /exports ngspice.dll | findstr ni_
```

Should show:
- `ni_instrument_register`- per-NR-iteration callback (extended: time, dt, mode, pre-solve RHS)
- `ni_topology_register`- one-time topology callback (node names, device state base offsets)

### Callback Signatures

The iteration callback signature (extended from original):

```c
typedef void (*NI_InstrumentCallback)(
    int iteration,       /* 0-based NR iteration */
    int matrixSize,      /* CKTmaxEqNum + 1 */
    double *rhs,         /* post-solve voltages (CKTrhs) */
    double *rhsOld,      /* previous iteration voltages (CKTrhsOld) */
    double *preSolveRhs, /* RHS after CKTload, before SMPsolve */
    double *state0,      /* CKTstate0 (device state vector) */
    int numStates,       /* CKTnumStates */
    int noncon,          /* CKTnoncon */
    int converged,       /* 1 if converged */
    double simTime,      /* CKTtime */
    double dt,           /* CKTdelta */
    int cktMode          /* CKTmode flags */
);
```

The topology callback fires once before the first NR iteration.
String arrays are passed as single pipe-delimited buffers (`char*`)
because koffi cannot marshal `char**` across FFI callbacks:

```c
typedef void (*NI_TopologyCallback)(
    const char *nodeNamesJoined,   /* "0|net1|net2|..." */
    int *nodeNumbers, int nodeCount,
    const char *devNamesJoined,    /* "q1|q2|r1|..." */
    const char *devTypesJoined,    /* "BJT|BJT|Resistor|..." */
    int *devStateBases, int devCount,
    int matrixSize, int numStates
);
```
