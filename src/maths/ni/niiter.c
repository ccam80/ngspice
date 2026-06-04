/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2001 AlansFixes
**********/

/*
 * NIiter(ckt,maxIter)
 *
 *  This subroutine performs the actual numerical iteration.
 *  It uses the sparse matrix stored in the circuit struct
 *  along with the matrix loading program, the load data, the
 *  convergence test function, and the convergence parameters
 */

#include "ngspice/ngspice.h"
#include "ngspice/trandefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/fteext.h"
#include "ngspice/devdefs.h"
#include "ngspice/const.h"
#include "../../spicelib/devices/bjt/bjtdefs.h"
#include "../../spicelib/devices/dio/diodefs.h"
#include "../../spicelib/devices/mos1/mos1defs.h"
#include "../../spicelib/devices/jfet/jfetdefs.h"
#include "../../spicelib/devices/cap/capdefs.h"
#include "../../spicelib/devices/ind/inddefs.h"
#include "../../spicelib/devices/res/resdefs.h"

/*
 * Minimal mirror of the sparse matrix element and frame structs.
 * Only the fields needed for CSC matrix extraction are declared here.
 * Layout matches MatrixElement and MatrixFrame in maths/sparse/spdefs.h
 * (with INITIALIZE=NO, which is the ngspice build configuration).
 */
typedef struct NiMatrixElement {
    double                   Real;
    double                   Imag;
    int                      Row;
    int                      Col;
    struct NiMatrixElement  *NextInRow;
    struct NiMatrixElement  *NextInCol;
} NiMatrixElement;

typedef struct NiMatrixFrame {
    double                AbsThreshold;
    int                   AllocatedSize;
    int                   AllocatedExtSize;
    int                   Complex;
    int                   CurrentSize;
    NiMatrixElement     **Diag;
    int                  *DoCmplxDirect;
    int                  *DoRealDirect;
    int                   Elements;
    int                   Error;
    int                   ExtSize;
    int                  *ExtToIntColMap;
    int                  *ExtToIntRowMap;
    int                   Factored;
    int                   Fillins;
    NiMatrixElement     **FirstInCol;
    NiMatrixElement     **FirstInRow;
} NiMatrixFrame;

/* ---- NR iteration instrumentation (extended) ---- */

/*
 * Struct-based per-iteration callback data.
 * Passed by pointer to ni_instrument_cb_v2 after every NR iteration.
 * Using a struct eliminates parameter ordering bugs in FFI.
 */
typedef struct {
    int iteration;
    int matrixSize;       /* CKTmaxEqNum + 1- matrix dimension reported to harness */
    int rhsBufSize;       /* SMPmatSize(CKTmatrix) + 1- actual rhs/rhsOld/preSolveRhs slot count;
                             may be < matrixSize when devices stamp into ground via TrashCan */
    double *rhs;
    double *rhsOld;
    double *preSolveRhs;
    double *state0;
    double *state1;
    double *state2;
    double *state3;
    int numStates;
    int noncon;
    int converged;
    double simTime;
    double dt;
    int cktMode;
    double ag0;
    double ag1;
    int integrateMethod;
    int order;
    int *matrixColPtr;
    int *matrixRowIdx;
    double *matrixVals;
    int matrixNnz;
    int *devConvFailed;
    int devConvCount;
    int numLimitEvents;
    int *limitDevIdx;
    int *limitJunctionId;
    double *limitVBefore;
    double *limitVAfter;
    int *limitWasLimited;
    /* Timestep-alignment fields (added for phase-aware capture) */
    double simTimeStart; /* CKTtime BEFORE the current NR solve began (set by dctran/niiter) */
    double phaseGmin;    /* live ckt->CKTdiagGmin (diagonal conductance in spFactor) */
    double phaseSrcFact; /* current source factor during src stepping, else 1 */
    int    phaseFlags;   /* bit0=inGminDynamic, bit1=inSrcSweep, bit2=inGminSpice3 */
} NiIterationData;

typedef void (*ni_instrument_cb_v2)(NiIterationData *data);

/*
 * Outer-loop callback: fired from dctran.c once per outer timestep-loop
 * iteration (after NIiter returns, after LTE check), ONLY during transient.
 * Exactly one of accepted/lteRejected/nrFailedRetry/finalFailure will be 1.
 */
typedef struct {
    double simTimeStart; /* CKTtime BEFORE the just-finished NR solve */
    double dt;           /* CKTdelta used */
    int    lteRejected;  /* 1 if LTE test failed and dt was cut */
    int    nrFailed;     /* 1 if NR did not converge */
    int    accepted;     /* 1 if the step was accepted and simTime advanced */
    int    finalFailure; /* 1 if dt <= delmin and giving up */
    double newDt;        /* next CKTdelta (after any trunc-error cut) */
} NiOuterData;

typedef void (*ni_outer_cb_t)(NiOuterData *data);

/*
 * Topology callback: fires once after the first CKTload, before first solve.
 * Provides node names and device state base offsets for state0 unpacking.
 *
 * String arrays are passed as single null-delimited buffers (char* with '\0'
 * between entries) because koffi cannot marshal char** across FFI callbacks.
 * The JS side splits on '\0' to recover individual strings.
 *
 * nodeNamesJoined: null-delimited node name strings ("0\0net1\0net2\0...")
 * nodeNumbers:     node number for each name
 * nodeCount:       number of entries
 * devNamesJoined:  null-delimited device instance names ("q1\0q2\0r1\0...")
 * devTypesJoined:  null-delimited device type names ("BJT\0BJT\0Resistor\0...")
 * devStateBases:   device state base offsets in CKTstate0
 * devCount:        number of devices
 * matrixSize:      CKTmaxEqNum + 1
 * numStates:       CKTnumStates
 */
typedef void (*NI_TopologyCallback)(
    const char *nodeNamesJoined,
    int *nodeNumbers,
    int nodeCount,
    const char *devNamesJoined,
    const char *devTypesJoined,
    int *devStateBases,
    int devCount,
    int matrixSize,
    int numStates,
    int *devNodeIndicesFlat,
    int *devNodeCounts
);

static ni_instrument_cb_v2 ni_instrument_cb = NULL;
static ni_outer_cb_t ni_outer_cb = NULL;
static NI_TopologyCallback ni_topology_cb = NULL;
static int ni_topology_sent = 0;
static double *ni_preSolveRhs = NULL;  /* scratch buffer for pre-solve RHS copy */
static int ni_preSolveRhsSize = 0;

/* Pre-LU matrix snapshot buffers (CSC format, populated after CKTload, before factorization) */
static int    *ni_mxColPtr     = NULL;
static int     ni_mxColPtrSize = 0;    /* allocated length of ni_mxColPtr (in ints) */
static int    *ni_mxRowIdx     = NULL;
static double *ni_mxVals       = NULL;
static int     ni_mxNnzBufSize = 0;    /* allocated length of ni_mxRowIdx / ni_mxVals */
static int     ni_mxNnz        = 0;    /* NNZ populated in the last pre-LU snapshot */

/* Phase flags set by cktop.c around each gmin/src ladder branch */
int ni_phase_flags = 0;     /* bit0=inGminDynamic, bit1=inSrcSweep, bit2=inGminSpice3 */
double ni_phase_gmin = 0.0;
double ni_phase_src_fact = 1.0;

/* ---- Item 9: Voltage limiting event collector ---- */

#define NI_MAX_LIMIT_EVENTS 256
#define NI_MAX_DEVICES 512

static int ni_limit_count = 0;
static int ni_limit_devIdx[NI_MAX_LIMIT_EVENTS];
static int ni_limit_junctionId[NI_MAX_LIMIT_EVENTS];
static double ni_limit_vBefore[NI_MAX_LIMIT_EVENTS];
static double ni_limit_vAfter[NI_MAX_LIMIT_EVENTS];
static int ni_limit_wasLimited[NI_MAX_LIMIT_EVENTS];

static GENinstance *ni_dev_map[NI_MAX_DEVICES];
static int ni_dev_map_count = 0;

__declspec(dllexport) void
ni_limit_reset(void)
{
    ni_limit_count = 0;
}

__declspec(dllexport) void
ni_limit_record(int devIdx, int junctionId, double vBefore, double vAfter)
{
    if (ni_limit_count >= NI_MAX_LIMIT_EVENTS) return;
    ni_limit_devIdx[ni_limit_count] = devIdx;
    ni_limit_junctionId[ni_limit_count] = junctionId;
    ni_limit_vBefore[ni_limit_count] = vBefore;
    ni_limit_vAfter[ni_limit_count] = vAfter;
    ni_limit_wasLimited[ni_limit_count] = (vBefore != vAfter) ? 1 : 0;
    ni_limit_count++;
}

__declspec(dllexport) int
ni_get_dev_index(GENinstance *inst)
{
    int i;
    for (i = 0; i < ni_dev_map_count; i++) {
        if (ni_dev_map[i] == inst) return i;
    }
    return -1;
}

/* Called from shared-lib consumer to register callbacks. */
__declspec(dllexport) void ni_instrument_register(ni_instrument_cb_v2 cb) {
    ni_instrument_cb = cb;
    ni_topology_sent = 0;  /* reset so topology fires on next run */
}

__declspec(dllexport) void ni_outer_register(ni_outer_cb_t cb) {
    ni_outer_cb = cb;
}

__declspec(dllexport) void ni_topology_register(NI_TopologyCallback cb) {
    ni_topology_cb = cb;
    ni_topology_sent = 0;
}

/* Called by cktop.c to set/clear phase flags before each gmin/src ladder call */
__declspec(dllexport) void ni_set_phase_flags(int flags, double gmin, double srcFact) {
    ni_phase_flags    = flags;
    ni_phase_gmin     = gmin;
    ni_phase_src_fact = srcFact;
}

/* Called by dctran.c to fire the outer-loop result callback */
__declspec(dllexport) void ni_fire_outer_cb(double simTimeStart, double dt,
                                             int lteRejected, int nrFailed,
                                             int accepted, int finalFailure,
                                             double newDt) {
    if (!ni_outer_cb) return;
    NiOuterData d;
    d.simTimeStart = simTimeStart;
    d.dt           = dt;
    d.lteRejected  = lteRejected;
    d.nrFailed     = nrFailed;
    d.accepted     = accepted;
    d.finalFailure = finalFailure;
    d.newDt        = newDt;
    ni_outer_cb(&d);
}

/* ============================================================================
 * AC sweep instrumentation (extended)
 *
 * Sibling of the DC/transient NR instrumentation (NiIterationData / NIiter
 * pre-LU CSC snapshot + pre-solve RHS + post-solve callback), specialised
 * for the AC linear path. AC has no NR loop: NIacIter performs a single
 * CKTacLoad -> SMPcLUfac (or SMPcReorder) -> SMPcSolve per frequency, with
 * the complex matrix split across .Real/.Imag in each spMatrix element and
 * the complex RHS split across CKTrhs/CKTirhs. After SMPcSolve, NIacIter
 * SWAPs CKTrhs<->CKTrhsOld and CKTirhs<->CKTirhsOld so the solution lives
 * in CKTrhsOld/CKTirhsOld for the rest of acan.c (consumed by CKTacDump).
 *
 * The hook fires once per frequency point. Capture sequencing inside
 * NIacIter (acan.c:297-298 invokes NIacIter once per freq):
 *   1. After CKTacLoad, before SMPcReorder/SMPcLUfac:
 *        ni_ac_capture_matrix(ckt)
 *      Snapshot the loaded complex Jacobian as external-coords CSC, both
 *      .Real and .Imag, before in-place LU overwrites them.
 *   2. After factor success, before SMPcSolve:
 *        ni_ac_capture_loaded_rhs(ckt)
 *      memcpy CKTrhs/CKTirhs- the per-frequency loaded RHS that drives
 *      the solve.
 *   3. After both SWAPs:
 *        ni_ac_capture_solution_and_fire(ckt)
 *      Solution lives in CKTrhsOld/CKTirhsOld; fill NiAcData, fire cb.
 *
 * Each helper short-circuits on `if (!ni_ac_cb) return;`- zero cost when
 * the bridge has not called ni_ac_register.
 *
 * cite: ref/ngspice/src/spicelib/analysis/acan.c:297-298- NIacIter call site
 * cite: ref/ngspice/src/maths/ni/niaciter.c:108-156- single-shot complex solve
 * ========================================================================== */

typedef struct {
    int    matrixSize;     /* CKTmaxEqNum + 1 */
    int    rhsBufSize;     /* SMPmatSize(CKTmatrix) + 1 */
    int    nnz;            /* CSC non-zero count for this frequency */
    int   *colPtr;         /* length matrixSize+1, CSC col offsets (external) */
    int   *rowIdx;         /* length nnz, external row index per entry */
    double *valsRe;        /* length nnz, loaded complex matrix Real */
    double *valsIm;        /* length nnz, loaded complex matrix Imag */
    double *rhsRe;         /* length rhsBufSize, loaded complex RHS Real */
    double *rhsIm;         /* length rhsBufSize, loaded complex RHS Imag */
    double *solRe;         /* length rhsBufSize, solution Real (= CKTrhsOld) */
    double *solIm;         /* length rhsBufSize, solution Imag (= CKTirhsOld) */
    double  omega;         /* CKTomega for this frequency point (rad/s) */
    double  freq;          /* omega / (2*pi), Hz */
} NiAcData;

typedef void (*ni_ac_cb_t)(NiAcData *data);

static ni_ac_cb_t ni_ac_cb = NULL;

/* Staging buffers- amortised across the frequency sweep. */
static int    *ni_ac_colPtr        = NULL;
static int     ni_ac_colPtrSize    = 0;     /* allocated length in ints */
static int    *ni_ac_rowIdx        = NULL;
static double *ni_ac_valsRe        = NULL;
static double *ni_ac_valsIm        = NULL;
static int     ni_ac_nnzBufSize    = 0;     /* allocated length of rowIdx/vals* */
static double *ni_ac_rhsRe         = NULL;
static double *ni_ac_rhsIm         = NULL;
static int     ni_ac_rhsBufSize    = 0;
static int     ni_ac_capturedNnz   = 0;     /* nnz emitted by ni_ac_capture_matrix this freq */
static int     ni_ac_capturedMsz   = 0;     /* matrixSize emitted by ni_ac_capture_matrix this freq */

__declspec(dllexport) void ni_ac_register(ni_ac_cb_t cb) {
    ni_ac_cb = cb;
}

/* Step 1: capture loaded complex matrix as external-coords CSC.
 * Faithful sibling of the NIiter pre-LU snapshot (this file, ~line 706-843):
 * same TrashCan-safe column cap (min(CKTmaxEqNum+1, AllocatedSize+1)), same
 * Ext<->Int inverse-map construction, same prefix-sum CSC build, additionally
 * emitting ep->Imag in parallel to ep->Real. */
__declspec(dllexport) void ni_ac_capture_matrix(CKTcircuit *ckt) {
    /* C89/MSVC: locals declared at top of block. */
    NiMatrixFrame   *mx;
    NiMatrixElement *ep;
    int   msz, colMax, nnz, intToExtSize;
    int  *IntToExtCol = NULL;
    int  *IntToExtRow = NULL;
    int  *tmpExtCol   = NULL;
    int  *tmpExtRow   = NULL;
    double *tmpValsRe = NULL;
    double *tmpValsIm = NULL;
    int  *cursor      = NULL;
    int   col, idx, ext, in, c, k, ec, pos, ii;
    int   extCol, extRow;

    if (!ni_ac_cb) return;
    mx = (NiMatrixFrame *)ckt->CKTmatrix;
    if (!mx) return;

    msz = ckt->CKTmaxEqNum + 1;
    colMax = msz;
    if (mx->AllocatedSize + 1 < colMax)
        colMax = mx->AllocatedSize + 1;

    /* Pass 1: count */
    nnz = 0;
    for (col = 1; col < colMax; col++) {
        for (ep = mx->FirstInCol[col]; ep != NULL; ep = ep->NextInCol) {
            nnz++;
        }
    }

    /* Realloc CSC output buffers if needed */
    if (ni_ac_colPtrSize < msz + 1) {
        FREE(ni_ac_colPtr);
        ni_ac_colPtr     = TMALLOC(int, msz + 1);
        ni_ac_colPtrSize = msz + 1;
    }
    if (ni_ac_nnzBufSize < nnz) {
        FREE(ni_ac_rowIdx);
        FREE(ni_ac_valsRe);
        FREE(ni_ac_valsIm);
        ni_ac_rowIdx     = TMALLOC(int,    nnz > 0 ? nnz : 1);
        ni_ac_valsRe     = TMALLOC(double, nnz > 0 ? nnz : 1);
        ni_ac_valsIm     = TMALLOC(double, nnz > 0 ? nnz : 1);
        ni_ac_nnzBufSize = nnz > 0 ? nnz : 1;
    }

    /* Build inverse maps; if absent (no preorder), identity. */
    intToExtSize = mx->AllocatedSize + 1;
    if (mx->ExtToIntColMap && mx->ExtToIntRowMap) {
        IntToExtCol = TMALLOC(int, intToExtSize + 1);
        IntToExtRow = TMALLOC(int, intToExtSize + 1);
        for (ii = 0; ii <= intToExtSize; ii++) {
            IntToExtCol[ii] = 0;
            IntToExtRow[ii] = 0;
        }
        for (ext = 1; ext <= mx->ExtSize; ext++) {
            in = mx->ExtToIntColMap[ext];
            if (in >= 1 && in <= intToExtSize) IntToExtCol[in] = ext;
        }
        for (ext = 1; ext <= mx->ExtSize; ext++) {
            in = mx->ExtToIntRowMap[ext];
            if (in >= 1 && in <= intToExtSize) IntToExtRow[in] = ext;
        }
    }

    /* Collect (extCol, extRow, valRe, valIm) triples. */
    tmpExtCol = TMALLOC(int,    nnz > 0 ? nnz : 1);
    tmpExtRow = TMALLOC(int,    nnz > 0 ? nnz : 1);
    tmpValsRe = TMALLOC(double, nnz > 0 ? nnz : 1);
    tmpValsIm = TMALLOC(double, nnz > 0 ? nnz : 1);
    idx = 0;
    for (col = 1; col < colMax; col++) {
        extCol = IntToExtCol ? IntToExtCol[col] : col;
        for (ep = mx->FirstInCol[col]; ep != NULL; ep = ep->NextInCol) {
            extRow = IntToExtRow ? IntToExtRow[ep->Row] : ep->Row;
            tmpExtCol[idx] = extCol;
            tmpExtRow[idx] = extRow;
            tmpValsRe[idx] = ep->Real;
            tmpValsIm[idx] = ep->Imag;
            idx++;
        }
    }

    /* Build CSC: zero colPtr, count-per-col, prefix-sum, scatter. */
    for (c = 0; c < msz + 1; c++) ni_ac_colPtr[c] = 0;
    for (k = 0; k < nnz; k++) {
        ec = tmpExtCol[k];
        if (ec >= 1 && ec < msz) ni_ac_colPtr[ec + 1]++;
    }
    for (c = 1; c <= msz; c++) ni_ac_colPtr[c] += ni_ac_colPtr[c - 1];

    cursor = TMALLOC(int, msz + 1);
    for (c = 0; c <= msz; c++) cursor[c] = ni_ac_colPtr[c];
    for (k = 0; k < nnz; k++) {
        ec = tmpExtCol[k];
        if (ec >= 1 && ec < msz) {
            pos = cursor[ec]++;
            ni_ac_rowIdx[pos] = tmpExtRow[k];
            ni_ac_valsRe[pos] = tmpValsRe[k];
            ni_ac_valsIm[pos] = tmpValsIm[k];
        }
    }

    FREE(cursor);
    FREE(tmpExtCol);
    FREE(tmpExtRow);
    FREE(tmpValsRe);
    FREE(tmpValsIm);
    FREE(IntToExtCol);
    FREE(IntToExtRow);

    ni_ac_capturedNnz = nnz;
    ni_ac_capturedMsz = msz;
}

/* Step 2: snapshot loaded complex RHS before SMPcSolve overwrites. */
__declspec(dllexport) void ni_ac_capture_loaded_rhs(CKTcircuit *ckt) {
    int sz;
    if (!ni_ac_cb) return;
    sz = SMPmatSize(ckt->CKTmatrix) + 1;
    if (ni_ac_rhsBufSize < sz) {
        FREE(ni_ac_rhsRe);
        FREE(ni_ac_rhsIm);
        ni_ac_rhsRe      = TMALLOC(double, sz);
        ni_ac_rhsIm      = TMALLOC(double, sz);
        ni_ac_rhsBufSize = sz;
    }
    memcpy(ni_ac_rhsRe, ckt->CKTrhs,  (size_t)sz * sizeof(double));
    memcpy(ni_ac_rhsIm, ckt->CKTirhs, (size_t)sz * sizeof(double));
}

/* Step 3: solution lives in CKTrhsOld / CKTirhsOld after both SWAPs. Fire cb. */
__declspec(dllexport) void ni_ac_capture_solution_and_fire(CKTcircuit *ckt) {
    NiAcData d;
    if (!ni_ac_cb) return;
    d.matrixSize = ckt->CKTmaxEqNum + 1;
    d.rhsBufSize = SMPmatSize(ckt->CKTmatrix) + 1;
    d.nnz        = ni_ac_capturedNnz;
    d.colPtr     = ni_ac_colPtr;
    d.rowIdx     = ni_ac_rowIdx;
    d.valsRe     = ni_ac_valsRe;
    d.valsIm     = ni_ac_valsIm;
    d.rhsRe      = ni_ac_rhsRe;
    d.rhsIm      = ni_ac_rhsIm;
    d.solRe      = ckt->CKTrhsOld;
    d.solIm      = ckt->CKTirhsOld;
    d.omega      = ckt->CKTomega;
    d.freq       = ckt->CKTomega / (2.0 * M_PI);
    ni_ac_cb(&d);
}

/*
 * Helper: join an array of strings into a single pipe-delimited buffer.
 * Returns a TMALLOC'd buffer; caller must FREE it.
 * Example: ["hello", "world"] -> "hello|world"
 * Pipe is used because koffi reads 'str' up to first \0.
 */
static char *
ni_join_strings(const char **strings, int count)
{
    int totalLen = 0, i;
    char *buf, *p;
    for (i = 0; i < count; i++) {
        totalLen += (int)strlen(strings[i] ? strings[i] : "") + 1; /* +1 for pipe or null term */
    }
    buf = TMALLOC(char, totalLen + 1);
    p = buf;
    for (i = 0; i < count; i++) {
        const char *s = strings[i] ? strings[i] : "";
        int len = (int)strlen(s);
        if (i > 0) { *p = '|'; p++; }
        memcpy(p, s, len);
        p += len;
    }
    *p = '\0';
    return buf;
}

/*
 * Pack node indices for a single device instance into devNodeIndicesFlat
 * at the given offset. Returns the number of nodes written.
 * Uses the device type name from DEVpublic.name to select the correct cast.
 */
static int
ni_pack_dev_nodes(GENinstance *inst, const char *typeName, int *flat, int offset)
{
    if (!inst || !typeName) return 0;

    if (strcmp(typeName, "BJT") == 0) {
        BJTinstance *d = (BJTinstance *)inst;
        flat[offset + 0] = d->BJTcolNode;
        flat[offset + 1] = d->BJTbaseNode;
        flat[offset + 2] = d->BJTemitNode;
        flat[offset + 3] = d->BJTsubstNode;
        return 4;
    }
    if (strcmp(typeName, "Diode") == 0) {
        DIOinstance *d = (DIOinstance *)inst;
        flat[offset + 0] = d->DIOposNode;
        flat[offset + 1] = d->DIOnegNode;
        return 2;
    }
    if (strcmp(typeName, "Mos1") == 0) {
        MOS1instance *d = (MOS1instance *)inst;
        flat[offset + 0] = d->MOS1dNode;
        flat[offset + 1] = d->MOS1gNode;
        flat[offset + 2] = d->MOS1sNode;
        flat[offset + 3] = d->MOS1bNode;
        return 4;
    }
    if (strcmp(typeName, "JFET") == 0) {
        JFETinstance *d = (JFETinstance *)inst;
        flat[offset + 0] = d->JFETdrainNode;
        flat[offset + 1] = d->JFETgateNode;
        flat[offset + 2] = d->JFETsourceNode;
        return 3;
    }
    if (strcmp(typeName, "Capacitor") == 0) {
        CAPinstance *d = (CAPinstance *)inst;
        flat[offset + 0] = d->CAPposNode;
        flat[offset + 1] = d->CAPnegNode;
        return 2;
    }
    if (strcmp(typeName, "Inductor") == 0) {
        INDinstance *d = (INDinstance *)inst;
        flat[offset + 0] = d->INDposNode;
        flat[offset + 1] = d->INDnegNode;
        return 2;
    }
    if (strcmp(typeName, "Resistor") == 0) {
        RESinstance *d = (RESinstance *)inst;
        flat[offset + 0] = d->RESposNode;
        flat[offset + 1] = d->RESnegNode;
        return 2;
    }
    return 0;
}

/* ---- Item 8: Per-element convergence detail ---- */

static void
ni_convTestAll(CKTcircuit *ckt, int *failedIndices, int maxDevices, int *failedCount)
{
    int i, devIdx;
    GENmodel *genmod;
    GENinstance *geninst;

    *failedCount = 0;
    devIdx = 0;

    for (i = 0; i < DEVmaxnum; i++) {
        if (!DEVices[i]) continue;
        for (genmod = ckt->CKThead[i]; genmod; genmod = genmod->GENnextModel) {
            for (geninst = genmod->GENinstances; geninst; geninst = geninst->GENnextInstance) {

                const char *tn = DEVices[i]->DEVpublic.name;
                int failed = 0;

                if (tn && strcmp(tn, "BJT") == 0) {
                    BJTmodel *bmod = (BJTmodel *)genmod;
                    BJTinstance *here = (BJTinstance *)geninst;
                    double vt = CONSTKoverQ * here->BJTtemp;
                    double vbe, vbc, delvbe, delvbc, cchat, cbhat, cc, cb, tol;
                    vbe = bmod->BJTtype * (
                        *(ckt->CKTrhsOld + here->BJTbasePrimeNode) -
                        *(ckt->CKTrhsOld + here->BJTemitPrimeNode));
                    vbc = bmod->BJTtype * (
                        *(ckt->CKTrhsOld + here->BJTbasePrimeNode) -
                        *(ckt->CKTrhsOld + here->BJTcolPrimeNode));
                    delvbe = vbe - *(ckt->CKTstate0 + here->BJTvbe);
                    delvbc = vbc - *(ckt->CKTstate0 + here->BJTvbc);
                    cchat = *(ckt->CKTstate0 + here->BJTcc) +
                        (*(ckt->CKTstate0 + here->BJTgm) +
                         *(ckt->CKTstate0 + here->BJTgo)) * delvbe -
                        (*(ckt->CKTstate0 + here->BJTgo) +
                         *(ckt->CKTstate0 + here->BJTgmu)) * delvbc;
                    cbhat = *(ckt->CKTstate0 + here->BJTcb) +
                        *(ckt->CKTstate0 + here->BJTgpi) * delvbe +
                        *(ckt->CKTstate0 + here->BJTgmu) * delvbc;
                    cc = *(ckt->CKTstate0 + here->BJTcc);
                    cb = *(ckt->CKTstate0 + here->BJTcb);
                    tol = ckt->CKTreltol * MAX(fabs(cchat), fabs(cc)) + ckt->CKTabstol;
                    if (fabs(cchat - cc) > tol) {
                        failed = 1;
                    } else {
                        tol = ckt->CKTreltol * MAX(fabs(cbhat), fabs(cb)) + ckt->CKTabstol;
                        if (fabs(cbhat - cb) > tol) {
                            failed = 1;
                        }
                    }
                } else if (tn && strcmp(tn, "Diode") == 0) {
                    DIOinstance *here = (DIOinstance *)geninst;
                    double vd, delvd, cdhat, cd, tol;
                    vd = *(ckt->CKTrhsOld + here->DIOposPrimeNode) -
                         *(ckt->CKTrhsOld + here->DIOnegNode);
                    delvd = vd - *(ckt->CKTstate0 + here->DIOvoltage);
                    cdhat = *(ckt->CKTstate0 + here->DIOcurrent) +
                            *(ckt->CKTstate0 + here->DIOconduct) * delvd;
                    cd = *(ckt->CKTstate0 + here->DIOcurrent);
                    tol = ckt->CKTreltol * MAX(fabs(cdhat), fabs(cd)) + ckt->CKTabstol;
                    if (fabs(cdhat - cd) > tol) {
                        failed = 1;
                    }
                } else if (tn && strcmp(tn, "Mos1") == 0) {
                    MOS1model *mmod = (MOS1model *)genmod;
                    MOS1instance *here = (MOS1instance *)geninst;
                    double vbs, vgs, vds, vbd, vgd, vgdo;
                    double delvbs, delvbd, delvgs, delvds, delvgd;
                    double cbhat, cdhat, tol;
                    vbs = mmod->MOS1type * (
                        *(ckt->CKTrhs + here->MOS1bNode) -
                        *(ckt->CKTrhs + here->MOS1sNodePrime));
                    vgs = mmod->MOS1type * (
                        *(ckt->CKTrhs + here->MOS1gNode) -
                        *(ckt->CKTrhs + here->MOS1sNodePrime));
                    vds = mmod->MOS1type * (
                        *(ckt->CKTrhs + here->MOS1dNodePrime) -
                        *(ckt->CKTrhs + here->MOS1sNodePrime));
                    vbd = vbs - vds;
                    vgd = vgs - vds;
                    vgdo = *(ckt->CKTstate0 + here->MOS1vgs) -
                           *(ckt->CKTstate0 + here->MOS1vds);
                    delvbs = vbs - *(ckt->CKTstate0 + here->MOS1vbs);
                    delvbd = vbd - *(ckt->CKTstate0 + here->MOS1vbd);
                    delvgs = vgs - *(ckt->CKTstate0 + here->MOS1vgs);
                    delvds = vds - *(ckt->CKTstate0 + here->MOS1vds);
                    delvgd = vgd - vgdo;
                    if (here->MOS1mode >= 0) {
                        cdhat = here->MOS1cd -
                            here->MOS1gbd * delvbd +
                            here->MOS1gmbs * delvbs +
                            here->MOS1gm * delvgs +
                            here->MOS1gds * delvds;
                    } else {
                        cdhat = here->MOS1cd -
                            (here->MOS1gbd - here->MOS1gmbs) * delvbd -
                            here->MOS1gm * delvgd +
                            here->MOS1gds * delvds;
                    }
                    cbhat = here->MOS1cbs + here->MOS1cbd +
                        here->MOS1gbd * delvbd +
                        here->MOS1gbs * delvbs;
                    tol = ckt->CKTreltol * MAX(fabs(cdhat), fabs(here->MOS1cd)) +
                          ckt->CKTabstol;
                    if (fabs(cdhat - here->MOS1cd) >= tol) {
                        failed = 1;
                    } else {
                        tol = ckt->CKTreltol *
                              MAX(fabs(cbhat), fabs(here->MOS1cbs + here->MOS1cbd)) +
                              ckt->CKTabstol;
                        if (fabs(cbhat - (here->MOS1cbs + here->MOS1cbd)) > tol) {
                            failed = 1;
                        }
                    }
                }

                if (failed && *failedCount < maxDevices) {
                    failedIndices[*failedCount] = devIdx;
                    (*failedCount)++;
                }
                devIdx++;
            }
        }
    }
}

/*
 * Send topology data once. Called from the NR loop after the first
 * successful CKTload.
 */
static void
ni_send_topology(CKTcircuit *ckt)
{
    CKTnode *node;
    GENmodel *model;
    GENinstance *inst;
    int i, nodeCount, devCount;
    const char **nodeNamePtrs;
    int *nodeNumbers;
    const char **devNamePtrs;
    const char **devTypePtrs;
    int *devStateBases;
    int *devNodeCounts;
    int *devNodeIndicesFlat;
    int totalNodeIndices;
    char *nodeNamesJoined, *devNamesJoined, *devTypesJoined;

    if (!ni_topology_cb || ni_topology_sent) return;
    ni_topology_sent = 1;

    /* Count nodes */
    nodeCount = 0;
    for (node = ckt->CKTnodes; node; node = node->next) nodeCount++;

    nodeNamePtrs = TMALLOC(const char*, nodeCount);
    nodeNumbers = TMALLOC(int, nodeCount);
    i = 0;
    for (node = ckt->CKTnodes; node; node = node->next) {
        nodeNamePtrs[i] = node->name ? node->name : "";
        nodeNumbers[i] = node->number;
        i++;
    }

    /* Count devices */
    devCount = 0;
    for (i = 0; i < DEVmaxnum; i++) {
        if (!DEVices[i]) continue;
        for (model = ckt->CKThead[i]; model; model = model->GENnextModel) {
            for (inst = model->GENinstances; inst; inst = inst->GENnextInstance) {
                devCount++;
            }
        }
    }

    devNamePtrs = TMALLOC(const char*, devCount);
    devTypePtrs = TMALLOC(const char*, devCount);
    devStateBases = TMALLOC(int, devCount);
    devNodeCounts = TMALLOC(int, devCount);

    /* Build instance-to-index map for limiting event recording */
    ni_dev_map_count = 0;

    {
        int di = 0;
        for (i = 0; i < DEVmaxnum; i++) {
            if (!DEVices[i]) continue;
            for (model = ckt->CKThead[i]; model; model = model->GENnextModel) {
                for (inst = model->GENinstances; inst; inst = inst->GENnextInstance) {
                    const char *tn = DEVices[i]->DEVpublic.name ? DEVices[i]->DEVpublic.name : "";
                    devNamePtrs[di]  = inst->GENname ? inst->GENname : "";
                    devTypePtrs[di]  = tn;
                    devStateBases[di] = inst->GENstate;
                    if (strcmp(tn, "BJT") == 0)        devNodeCounts[di] = 4;
                    else if (strcmp(tn, "Diode") == 0)     devNodeCounts[di] = 2;
                    else if (strcmp(tn, "Mos1") == 0)      devNodeCounts[di] = 4;
                    else if (strcmp(tn, "JFET") == 0)      devNodeCounts[di] = 3;
                    else if (strcmp(tn, "Capacitor") == 0) devNodeCounts[di] = 2;
                    else if (strcmp(tn, "Inductor") == 0)  devNodeCounts[di] = 2;
                    else if (strcmp(tn, "Resistor") == 0)  devNodeCounts[di] = 2;
                    else                                    devNodeCounts[di] = 0;
                    if (di < NI_MAX_DEVICES) {
                        ni_dev_map[di] = inst;
                    }
                    di++;
                }
            }
        }
        ni_dev_map_count = di < NI_MAX_DEVICES ? di : NI_MAX_DEVICES;
    }

    /* Compute total flat array size and allocate */
    totalNodeIndices = 0;
    for (i = 0; i < devCount; i++) totalNodeIndices += devNodeCounts[i];
    devNodeIndicesFlat = TMALLOC(int, totalNodeIndices > 0 ? totalNodeIndices : 1);

    /* Pack node indices into flat array */
    {
        int di = 0, flatOffset = 0;
        for (i = 0; i < DEVmaxnum; i++) {
            if (!DEVices[i]) continue;
            for (model = ckt->CKThead[i]; model; model = model->GENnextModel) {
                for (inst = model->GENinstances; inst; inst = inst->GENnextInstance) {
                    ni_pack_dev_nodes(inst, devTypePtrs[di], devNodeIndicesFlat, flatOffset);
                    flatOffset += devNodeCounts[di];
                    di++;
                }
            }
        }
    }

    /* Build pipe-delimited string buffers */
    nodeNamesJoined = ni_join_strings(nodeNamePtrs, nodeCount);
    devNamesJoined  = ni_join_strings(devNamePtrs,  devCount);
    devTypesJoined  = ni_join_strings(devTypePtrs,  devCount);

    ni_topology_cb(nodeNamesJoined, nodeNumbers, nodeCount,
                   devNamesJoined, devTypesJoined,
                   devStateBases, devCount,
                   ckt->CKTmaxEqNum + 1, ckt->CKTnumStates,
                   devNodeIndicesFlat, devNodeCounts);

    FREE(nodeNamesJoined);
    FREE(devNamesJoined);
    FREE(devTypesJoined);
    FREE(nodeNamePtrs);
    FREE(nodeNumbers);
    FREE(devNamePtrs);
    FREE(devTypePtrs);
    FREE(devStateBases);
    FREE(devNodeCounts);
    FREE(devNodeIndicesFlat);
}

/* Limit the number of 'singular matrix' warnings */
static int msgcount = 0;

/* NIiter() - return value is non-zero for convergence failure */

int
NIiter(CKTcircuit *ckt, int maxIter)
{
    double startTime, *OldCKTstate0 = NULL;
    int error, i, j;

    int iterno = 0;
    int ipass = 0;

    /* some convergence issues that get resolved by increasing max iter */
    if (maxIter < 100)
        maxIter = 100;

    if ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)) {
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        error = CKTload(ckt);
        if (error)
            return(error);
        return(OK);
    }

#ifdef WANT_SENSE2
    if (ckt->CKTsenInfo) {
        error = NIsenReinit(ckt);
        if (error)
            return(error);
    }
#endif

    if (ckt->CKTniState & NIUNINITIALIZED) {
        error = NIreinit(ckt);
        if (error) {
#ifdef STEPDEBUG
            printf("re-init returned error \n");
#endif
            return(error);
        }
    }

    /* OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1); */

    for (;;) {

        ckt->CKTnoncon = 0;

        /* Item 9: Reset limiting event collector before each NR iteration */
        ni_limit_reset();

#ifdef NEWPRED
        if (!(ckt->CKTmode & MODEINITPRED))
#endif
        {

            error = CKTload(ckt);
            /* printf("loaded, noncon is %d\n", ckt->CKTnoncon); */
            /* fflush(stdout); */
            iterno++;
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                printf("load returned error \n");
#endif
                FREE(OldCKTstate0);
                return (error);
            }

            /* printf("after loading, before solving\n"); */
            /* CKTdump(ckt); */

            /* Pre-LU matrix snapshot: capture MNA M immediately after CKTload,
             * before SMPpreOrder/SMPreorder/SMPluFac overwrites .Real in place.
             * Gated on ni_instrument_cb- zero cost when instrumentation is off.
             *
             * Column loop bound: we iterate FirstInCol[col] for col in
             * [1, min(CKTmaxEqNum + 1, AllocatedSize + 1)). FirstInCol is
             * allocated with exactly AllocatedSize + 1 entries
             * (spalloc.c:220, SP_CALLOC(..., SizePlusOne)), but CKTmaxEqNum
             * is tracked independently by the circuit layer and can exceed
             * AllocatedSize when a device stamps into Row 0 or Col 0 via
             * spGetElement's TrashCan branch (spbuild.c:272) without
             * triggering EnlargeMatrix. This happens whenever a BJT has a
             * terminal on ground node 0 (the bjtsetup.c TSTALLOC for
             * BJTcolColPtr, BJTcolColPrimePtr, etc. routes to TrashCan
             * instead of expanding the sparse frame). Iterating up to msz
             * in that case reads one-or-more entries past the end of
             * FirstInCol and dereferences garbage- confirmed segfault on
             * pnp-cc-harness.dts, root-caused 2026-04-12. Columns in
             * [AllocatedSize + 1, msz) are guaranteed to contain no stamped
             * elements (their equations exist only in the TrashCan slot),
             * so reporting them as empty in the CSC output is exact, not
             * an approximation. */
            if (ni_instrument_cb) {
                int ni_pre_msz = ckt->CKTmaxEqNum + 1;
                int ni_pre_nnz = 0;
                int ni_pre_col, ni_pre_idx;
                NiMatrixElement *ni_pre_ep;
                NiMatrixFrame *ni_pre_mx = (NiMatrixFrame *)ckt->CKTmatrix;
                int ni_pre_colMax;

                if (ni_pre_mx) {
                    /* Locals declared at top of block for C89/MSVC compatibility */
                    int     ni_intToExtSize;
                    int    *ni_IntToExtCol;
                    int    *ni_IntToExtRow;
                    int    *ni_tmpExtCol;
                    int    *ni_tmpExtRow;
                    double *ni_tmpVals;
                    int    *ni_cursor;
                    int     ni_extCol, ni_extRow;
                    int     ni_ii, ni_ext, ni_int, ni_c, ni_k, ni_ec, ni_pos;

                    ni_pre_colMax = ni_pre_msz;
                    if (ni_pre_mx->AllocatedSize + 1 < ni_pre_colMax)
                        ni_pre_colMax = ni_pre_mx->AllocatedSize + 1;

                    /* Pass 1: count non-zeros in the pre-factor sparsity pattern */
                    for (ni_pre_col = 1; ni_pre_col < ni_pre_colMax; ni_pre_col++) {
                        for (ni_pre_ep = ni_pre_mx->FirstInCol[ni_pre_col];
                             ni_pre_ep != NULL;
                             ni_pre_ep = ni_pre_ep->NextInCol) {
                            ni_pre_nnz++;
                        }
                    }

                    /* Realloc colPtr buffer if needed */
                    if (ni_mxColPtrSize < ni_pre_msz + 1) {
                        FREE(ni_mxColPtr);
                        ni_mxColPtr     = TMALLOC(int, ni_pre_msz + 1);
                        ni_mxColPtrSize = ni_pre_msz + 1;
                    }
                    /* Realloc rowIdx / vals buffers if needed */
                    if (ni_mxNnzBufSize < ni_pre_nnz) {
                        FREE(ni_mxRowIdx);
                        FREE(ni_mxVals);
                        ni_mxRowIdx     = TMALLOC(int,    ni_pre_nnz);
                        ni_mxVals       = TMALLOC(double, ni_pre_nnz);
                        ni_mxNnzBufSize = ni_pre_nnz;
                    }

                    /* Build inverse maps: IntToExtCol[intIdx] = extIdx,
                     * IntToExtRow[intIdx] = extIdx.
                     * ExtToIntColMap[ext] = int, so we invert by scanning.
                     * If the maps are NULL (before any SMPpreOrder), internal
                     * indices equal external indices- use identity mapping. */
                    ni_intToExtSize = ni_pre_mx->AllocatedSize + 1;
                    ni_IntToExtCol  = NULL;
                    ni_IntToExtRow  = NULL;
                    if (ni_pre_mx->ExtToIntColMap && ni_pre_mx->ExtToIntRowMap) {
                        ni_IntToExtCol = TMALLOC(int, ni_intToExtSize + 1);
                        ni_IntToExtRow = TMALLOC(int, ni_intToExtSize + 1);
                        /* initialise to 0 (unmapped entries stay 0 = ground) */
                        for (ni_ii = 0; ni_ii <= ni_intToExtSize; ni_ii++) {
                            ni_IntToExtCol[ni_ii] = 0;
                            ni_IntToExtRow[ni_ii] = 0;
                        }
                        /* ExtToIntColMap is indexed [1..ExtSize] */
                        for (ni_ext = 1; ni_ext <= ni_pre_mx->ExtSize; ni_ext++) {
                            ni_int = ni_pre_mx->ExtToIntColMap[ni_ext];
                            if (ni_int >= 1 && ni_int <= ni_intToExtSize)
                                ni_IntToExtCol[ni_int] = ni_ext;
                        }
                        for (ni_ext = 1; ni_ext <= ni_pre_mx->ExtSize; ni_ext++) {
                            ni_int = ni_pre_mx->ExtToIntRowMap[ni_ext];
                            if (ni_int >= 1 && ni_int <= ni_intToExtSize)
                                ni_IntToExtRow[ni_int] = ni_ext;
                        }
                    }

                    /* Scratch buffer: collect (extCol, extRow, val) triples,
                     * then scatter into CSC order by external column. */
                    ni_tmpExtCol = TMALLOC(int,    ni_pre_nnz > 0 ? ni_pre_nnz : 1);
                    ni_tmpExtRow = TMALLOC(int,    ni_pre_nnz > 0 ? ni_pre_nnz : 1);
                    ni_tmpVals   = TMALLOC(double, ni_pre_nnz > 0 ? ni_pre_nnz : 1);

                    /* Collect pass: walk FirstInCol by internal column, remap. */
                    ni_pre_idx = 0;
                    for (ni_pre_col = 1; ni_pre_col < ni_pre_colMax; ni_pre_col++) {
                        ni_extCol = ni_IntToExtCol
                            ? ni_IntToExtCol[ni_pre_col]
                            : ni_pre_col;
                        for (ni_pre_ep = ni_pre_mx->FirstInCol[ni_pre_col];
                             ni_pre_ep != NULL;
                             ni_pre_ep = ni_pre_ep->NextInCol) {
                            ni_extRow = ni_IntToExtRow
                                ? ni_IntToExtRow[ni_pre_ep->Row]
                                : ni_pre_ep->Row;
                            ni_tmpExtCol[ni_pre_idx] = ni_extCol;
                            ni_tmpExtRow[ni_pre_idx] = ni_extRow;
                            ni_tmpVals[ni_pre_idx]   = ni_pre_ep->Real;
                            ni_pre_idx++;
                        }
                    }

                    /* Build CSC: count entries per external column, then fill.
                     * Pass 2a: zero colPtr, then count per-column nnz. */
                    for (ni_c = 0; ni_c < ni_pre_msz + 1; ni_c++)
                        ni_mxColPtr[ni_c] = 0;
                    for (ni_k = 0; ni_k < ni_pre_nnz; ni_k++) {
                        ni_ec = ni_tmpExtCol[ni_k];
                        if (ni_ec >= 1 && ni_ec < ni_pre_msz)
                            ni_mxColPtr[ni_ec + 1]++;
                    }
                    /* Prefix-sum to get colPtr[0..msz]. */
                    for (ni_c = 1; ni_c <= ni_pre_msz; ni_c++)
                        ni_mxColPtr[ni_c] += ni_mxColPtr[ni_c - 1];

                    /* Pass 2b: scatter entries into rowIdx/vals using a
                     * running cursor array. */
                    ni_cursor = TMALLOC(int, ni_pre_msz + 1);
                    for (ni_c = 0; ni_c <= ni_pre_msz; ni_c++)
                        ni_cursor[ni_c] = ni_mxColPtr[ni_c];
                    for (ni_k = 0; ni_k < ni_pre_nnz; ni_k++) {
                        ni_ec = ni_tmpExtCol[ni_k];
                        if (ni_ec >= 1 && ni_ec < ni_pre_msz) {
                            ni_pos = ni_cursor[ni_ec]++;
                            ni_mxRowIdx[ni_pos] = ni_tmpExtRow[ni_k];
                            ni_mxVals[ni_pos]   = ni_tmpVals[ni_k];
                        }
                    }
                    FREE(ni_cursor);

                    FREE(ni_tmpExtCol);
                    FREE(ni_tmpExtRow);
                    FREE(ni_tmpVals);
                    FREE(ni_IntToExtCol);
                    FREE(ni_IntToExtRow);

                    ni_mxNnz = ni_pre_nnz;
                }
            }

            if (!(ckt->CKTniState & NIDIDPREORDER)) {
                error = SMPpreOrder(ckt->CKTmatrix);
                if (error) {
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("pre-order returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error); /* badly formed matrix */
                }
                ckt->CKTniState |= NIDIDPREORDER;
            }

            if ((ckt->CKTmode & MODEINITJCT) ||
                ((ckt->CKTmode & MODEINITTRAN) && (iterno == 1)))
            {
                ckt->CKTniState |= NISHOULDREORDER;
            }

            if (ckt->CKTniState & NISHOULDREORDER) {
                startTime = SPfrontEnd->IFseconds();
                error = SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                   ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
                ckt->CKTstat->STATreorderTime +=
                    SPfrontEnd->IFseconds() - startTime;
                if (error) {
                    /* new feature - we can now find out something about what is
                     * wrong - so we ask for the troublesome entry
                     * Limit the number of messages to 6, if not 'set ngdebug'.
                     */
                    if (ft_ngdebug || msgcount < 6) {
                        SMPgetError(ckt->CKTmatrix, &i, &j);
                        if(eq(NODENAME(ckt, i), NODENAME(ckt, j)))
                            SPfrontEnd->IFerrorf(ERR_WARNING, "singular matrix:  check node %s\n", NODENAME(ckt, i));
                        else
                            SPfrontEnd->IFerrorf(ERR_WARNING, "singular matrix:  check nodes %s and %s\n", NODENAME(ckt, i), NODENAME(ckt, j));
                        msgcount += 1;
                    }
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("reorder returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error); /* can't handle these errors - pass up! */
                }
                ckt->CKTniState &= ~NISHOULDREORDER;
            } else {
                startTime = SPfrontEnd->IFseconds();
                error = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                 ckt->CKTdiagGmin);
                ckt->CKTstat->STATdecompTime +=
                    SPfrontEnd->IFseconds() - startTime;
                if (error) {
                    if (error == E_SINGULAR) {
                        ckt->CKTniState |= NISHOULDREORDER;
                        DEBUGMSG(" forced reordering....\n");
                        continue;
                    }
                    /* CKTload(ckt); */
                    /* SMPprint(ckt->CKTmatrix, stdout); */
                    /* seems to be singular - pass the bad news up */
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("lufac returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error);
                }
            }

            /* moved it to here as if xspice is included then CKTload changes
               CKTnumStates the first time it is run */
            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));

            /* Send topology once after first successful load */
            ni_send_topology(ckt);

            /* Capture pre-solve RHS (loaded stamp contributions) before SMPsolve overwrites.
             * Buffer size MUST be SMPmatSize+1, not CKTmaxEqNum+1: ckt->CKTrhs is allocated
             * to SMPmatSize+1 doubles in nireinit.c. CKTmaxEqNum can exceed SMPmatSize when
             * devices stamp into ground row/col via TrashCan (spbuild.c)- using
             * CKTmaxEqNum+1 reads OOB from CKTrhs. */
            if (ni_instrument_cb) {
                int sz = SMPmatSize(ckt->CKTmatrix) + 1;
                if (ni_preSolveRhsSize < sz) {
                    FREE(ni_preSolveRhs);
                    ni_preSolveRhs = TMALLOC(double, sz);
                    ni_preSolveRhsSize = sz;
                }
                memcpy(ni_preSolveRhs, ckt->CKTrhs, sz * sizeof(double));
            }

            startTime = SPfrontEnd->IFseconds();
            SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
            ckt->CKTstat->STATsolveTime +=
                SPfrontEnd->IFseconds() - startTime;
#ifdef STEPDEBUG
            /*XXXX*/
            if (ckt->CKTrhs[0] != 0.0)
                printf("NIiter: CKTrhs[0] = %g\n", ckt->CKTrhs[0]);
            if (ckt->CKTrhsSpare[0] != 0.0)
                printf("NIiter: CKTrhsSpare[0] = %g\n", ckt->CKTrhsSpare[0]);
            if (ckt->CKTrhsOld[0] != 0.0)
                printf("NIiter: CKTrhsOld[0] = %g\n", ckt->CKTrhsOld[0]);
            /*XXXX*/
#endif
            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

            if (iterno > maxIter) {
                ckt->CKTstat->STATnumIter += iterno;
                /* we don't use this info during transient analysis */
                if (ckt->CKTcurrentAnalysis != DOING_TRAN) {
                    FREE(errMsg);
                    errMsg = copy("Too many iterations without convergence");
#ifdef STEPDEBUG
                    fprintf(stderr, "too many iterations without convergence: %d iter's (max iter == %d)\n",
                    iterno, maxIter);
#endif
                }
                FREE(OldCKTstate0);
                return(E_ITERLIM);
            }

            if ((ckt->CKTnoncon == 0) && (iterno != 1))
                ckt->CKTnoncon = NIconvTest(ckt);
            else
                ckt->CKTnoncon = 1;

#ifdef STEPDEBUG
            printf("noncon is %d\n", ckt->CKTnoncon);
#endif
            /* Instrumentation callback- fires after solve + convergence check */
            if (ni_instrument_cb) {
                int ni_msz = ckt->CKTmaxEqNum + 1;
                int ni_rhs_buf_size = SMPmatSize(ckt->CKTmatrix) + 1;
                {
                    NiIterationData ni_data;
                    int ni_convFailed[NI_MAX_DEVICES];
                    int ni_convFailedCount = 0;

                    ni_data.iteration       = iterno - 1;
                    ni_data.matrixSize      = ni_msz;
                    ni_data.rhsBufSize      = ni_rhs_buf_size;
                    ni_data.rhs             = ckt->CKTrhs;
                    ni_data.rhsOld          = ckt->CKTrhsOld;
                    ni_data.preSolveRhs     = ni_preSolveRhs;
                    ni_data.state0          = ckt->CKTstate0;
                    ni_data.state1          = ckt->CKTstate1;
                    ni_data.state2          = ckt->CKTstate2;
                    ni_data.state3          = ckt->CKTstate3;
                    ni_data.numStates       = ckt->CKTnumStates;
                    ni_data.noncon          = ckt->CKTnoncon;
                    ni_data.converged       = (ckt->CKTnoncon == 0 && iterno != 1) ? 1 : 0;
                    ni_data.simTime         = ckt->CKTtime;
                    ni_data.dt              = ckt->CKTdelta;
                    ni_data.cktMode         = ckt->CKTmode;
                    ni_data.ag0             = ckt->CKTag[0];
                    ni_data.ag1             = ckt->CKTag[1];
                    ni_data.integrateMethod = ckt->CKTintegrateMethod;
                    ni_data.order           = ckt->CKTorder;
                    ni_data.matrixColPtr    = ni_mxColPtr;
                    ni_data.matrixRowIdx    = ni_mxRowIdx;
                    ni_data.matrixVals      = ni_mxVals;
                    ni_data.matrixNnz       = ni_mxNnz;
                    /* Timestep-alignment fields */
                    ni_data.simTimeStart    = ckt->CKTsimTimeStart;
                    /* The live diagonal conductance added in spFactor (lines
                       above), not the phase-flag value: during gmin stepping the
                       two agree, but OPtran's pseudo-transient runs with the gmin
                       gillespie_src left on the diagonal (cktop.c:647) while the
                       phase-flag gmin is 0, so capture the actual value. */
                    ni_data.phaseGmin       = ckt->CKTdiagGmin;
                    ni_data.phaseSrcFact    = ni_phase_src_fact;
                    ni_data.phaseFlags      = ni_phase_flags;

                    /* Item 8: Per-element convergence detail */
                    ni_convTestAll(ckt, ni_convFailed, NI_MAX_DEVICES, &ni_convFailedCount);
                    ni_data.devConvFailed = ni_convFailedCount > 0 ? ni_convFailed : NULL;
                    ni_data.devConvCount  = ni_convFailedCount;

                    /* Item 9: Limiting events from this iteration */
                    ni_data.numLimitEvents  = ni_limit_count;
                    ni_data.limitDevIdx     = ni_limit_count > 0 ? ni_limit_devIdx : NULL;
                    ni_data.limitJunctionId = ni_limit_count > 0 ? ni_limit_junctionId : NULL;
                    ni_data.limitVBefore    = ni_limit_count > 0 ? ni_limit_vBefore : NULL;
                    ni_data.limitVAfter     = ni_limit_count > 0 ? ni_limit_vAfter : NULL;
                    ni_data.limitWasLimited = ni_limit_count > 0 ? ni_limit_wasLimited : NULL;

                    ni_instrument_cb(&ni_data);
                }

            }
        }

        if ((ckt->CKTnodeDamping != 0) && (ckt->CKTnoncon != 0) &&
            ((ckt->CKTmode & MODETRANOP) || (ckt->CKTmode & MODEDCOP)) &&
            (iterno > 1))
        {
            CKTnode *node;
            double diff, maxdiff = 0;
            for (node = ckt->CKTnodes->next; node; node = node->next)
                if (node->type == SP_VOLTAGE) {
                    diff = fabs(ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number]);
                    if (maxdiff < diff)
                        maxdiff = diff;
                }

            if (maxdiff > 10) {
                double damp_factor = 10 / maxdiff;
                if (damp_factor < 0.1)
                    damp_factor = 0.1;
                for (node = ckt->CKTnodes->next; node; node = node->next) {
                    diff = ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number];
                    ckt->CKTrhs[node->number] =
                        ckt->CKTrhsOld[node->number] + (damp_factor * diff);
                }
                for (i = 0; i < ckt->CKTnumStates; i++) {
                    diff = ckt->CKTstate0[i] - OldCKTstate0[i];
                    ckt->CKTstate0[i] = OldCKTstate0[i] + (damp_factor * diff);
                }
            }
        }

        if (ckt->CKTmode & MODEINITFLOAT) {
            if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset) {
                if (ipass)
                    ckt->CKTnoncon = ipass;
                ipass = 0;
            }
            if (ckt->CKTnoncon == 0) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
                return(OK);
            }
        } else if (ckt->CKTmode & MODEINITJCT) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
            ckt->CKTniState |= NISHOULDREORDER;
        } else if (ckt->CKTmode & MODEINITFIX) {
            if (ckt->CKTnoncon == 0)
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
            ipass = 1;
        } else if (ckt->CKTmode & MODEINITSMSIG) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITTRAN) {
            if (iterno <= 1)
                ckt->CKTniState |= NISHOULDREORDER;
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITPRED) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else {
            ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
            printf("bad initf state \n");
#endif
            FREE(OldCKTstate0);
            return(E_INTERN);
            /* impossible - no such INITF flag! */
        }

        /* build up the lvnim1 array from the lvn array */
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        /* printf("after loading, after solving\n"); */
        /* CKTdump(ckt); */
    }
    /*NOTREACHED*/
}

void NIresetwarnmsg(void) {
    msgcount = 0;
};
