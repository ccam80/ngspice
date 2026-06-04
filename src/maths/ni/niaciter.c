/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified 1999 Emmanuel Rouat
**********/

    /*
     * NIacIter(ckt)
     *
     *  This subroutine performs the actual numerical iteration.
     *  It uses the sparse matrix stored in the NIstruct by NIinit,
     *  along with the matrix loading program, the load data, the
     *  convergence test function, and the convergence parameters
     * - return value is non-zero for convergence failure 
     */

#include "ngspice/ngspice.h"
#include "ngspice/trandefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/sperror.h"

/* AC sweep instrumentation- defined in niiter.c. Each helper is a no-op
 * (single null-pointer test) when ni_ac_register has not been called, so
 * stock ngspice users pay zero cost. Capture sequencing matches the
 * comments at the top of the AC instrumentation block in niiter.c. */
extern void ni_ac_capture_matrix(CKTcircuit *ckt);
extern void ni_ac_capture_loaded_rhs(CKTcircuit *ckt);
extern void ni_ac_capture_solution_and_fire(CKTcircuit *ckt);

#ifdef RFSPICE
     // We don't need to reload the AC matrix for every port analysis
     // So we split the NIacIter in two functions
int NIspPreload(CKTcircuit* ckt)
{

    int error;
    int ignore;
    double startTime;

retry:
    ckt->CKTnoncon = 0;

    error = CKTacLoad(ckt);

    if (error) return(error);

    if (ckt->CKTniState & NIACSHOULDREORDER) {
        startTime = SPfrontEnd->IFseconds();
        error = SMPcReorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
            ckt->CKTpivotRelTol, &ignore);
        ckt->CKTstat->STATreorderTime +=
            SPfrontEnd->IFseconds() - startTime;
        ckt->CKTniState &= ~NIACSHOULDREORDER;
        if (error != 0) {
            /* either singular equations or no memory, in either case,
             * let caller handle problem
             */
            return(error);
        }
    }
    else {
        startTime = SPfrontEnd->IFseconds();
        error = SMPcLUfac(ckt->CKTmatrix, ckt->CKTpivotAbsTol);
        ckt->CKTstat->STATdecompTime +=
            SPfrontEnd->IFseconds() - startTime;
        if (error != 0) {
            if (error == E_SINGULAR) {
                /* the problem is that the matrix can't be solved with the
                 * current LU factorization.  Maybe if we reload and
                 * try to reorder again it will help...
                 */
                ckt->CKTniState |= NIACSHOULDREORDER;
                goto retry;
            }
            return(error); /* can't handle E_BADMATRIX, so let caller */
        }
    }

    return (OK);
};

int NIspSolve(CKTcircuit* ckt)
{
    double startTime;
    startTime = SPfrontEnd->IFseconds();
    SMPcSolve(ckt->CKTmatrix, ckt->CKTrhs,
        ckt->CKTirhs, ckt->CKTrhsSpare,
        ckt->CKTirhsSpare);
    ckt->CKTstat->STATsolveTime += SPfrontEnd->IFseconds() - startTime;

    ckt->CKTrhs[0] = 0;
    ckt->CKTrhsSpare[0] = 0;
    ckt->CKTrhsOld[0] = 0;
    ckt->CKTirhs[0] = 0;
    ckt->CKTirhsSpare[0] = 0;
    ckt->CKTirhsOld[0] = 0;

    SWAP(double*, ckt->CKTirhs, ckt->CKTirhsOld);

    SWAP(double*, ckt->CKTrhs, ckt->CKTrhsOld);
    return (OK);
};

#endif

int
NIacIter(CKTcircuit *ckt)
{
    int error;
    int ignore;
    double startTime;

retry:
    ckt->CKTnoncon=0;

    error = CKTacLoad(ckt);
    if(error) return(error);

    /* Step 1: snapshot the loaded complex Jacobian as CSC (re+im) before
     * SMPcLUfac overwrites .Real/.Imag with L/U. No-op without registration. */
    ni_ac_capture_matrix(ckt);

    if(ckt->CKTniState & NIACSHOULDREORDER) {
	startTime = SPfrontEnd->IFseconds();
        error = SMPcReorder(ckt->CKTmatrix,ckt->CKTpivotAbsTol,
                ckt->CKTpivotRelTol,&ignore);
	ckt->CKTstat->STATreorderTime +=
		SPfrontEnd->IFseconds()- startTime;
        ckt->CKTniState &= ~NIACSHOULDREORDER;
        if(error != 0) {
            /* either singular equations or no memory, in either case,
             * let caller handle problem
             */
            return(error);
        }
    } else {
	startTime = SPfrontEnd->IFseconds();
        error = SMPcLUfac(ckt->CKTmatrix,ckt->CKTpivotAbsTol);
	ckt->CKTstat->STATdecompTime += 
		SPfrontEnd->IFseconds()-startTime;
        if(error != 0) {
            if(error == E_SINGULAR) {
                /* the problem is that the matrix can't be solved with the
                 * current LU factorization.  Maybe if we reload and
                 * try to reorder again it will help...
                 */
                ckt->CKTniState |= NIACSHOULDREORDER;
                goto retry;
            }
            return(error); /* can't handle E_BADMATRIX, so let caller */
        }
    }

    /* Step 2: snapshot loaded complex RHS before SMPcSolve overwrites
     * CKTrhs/CKTirhs with the solution. No-op without registration. */
    ni_ac_capture_loaded_rhs(ckt);

    startTime = SPfrontEnd->IFseconds();
    SMPcSolve(ckt->CKTmatrix,ckt->CKTrhs,
            ckt->CKTirhs, ckt->CKTrhsSpare,
            ckt->CKTirhsSpare);
    ckt->CKTstat->STATsolveTime += SPfrontEnd->IFseconds() - startTime;

    ckt->CKTrhs[0] = 0;
    ckt->CKTrhsSpare[0] = 0;
    ckt->CKTrhsOld[0] = 0;
    ckt->CKTirhs[0] = 0;
    ckt->CKTirhsSpare[0] = 0;
    ckt->CKTirhsOld[0] = 0;

    SWAP(double *, ckt->CKTirhs, ckt->CKTirhsOld);

    SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);

    /* Step 3: solution lives in CKTrhsOld / CKTirhsOld after the SWAPs.
     * Fill NiAcData and fire the registered callback. No-op without registration. */
    ni_ac_capture_solution_and_fire(ckt);

    return(OK);
}
