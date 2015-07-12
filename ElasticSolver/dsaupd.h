/*
  In this header file, I have defined a simplified function
  call to the ARPACK solver routine for a simple, symmetric
  eigenvalue problem Av = lv.  The looping procedure and
  final extraction of eigenvalues and vectors is handled
  automatically.  Most of the parameters to the FORTRAN
  functions are hidden from the user, since most of them are
  determined from user input anyway.

  The remaining parameters to the function calls are as follows:

    dsaupd(int n, int nev, double *Evals)
    dsaupd(int n, int nev, doubel *Evals, double **Evecs)

    n: the order of the square matrix A
    nev: the number of eigenvalues to be found, starting at the
         bottom.  Note that the highest eigenvalues, or some
	 other choices, can be found.  For now, the choice of
	 the lowest nev eigenvalues is hard-coded.
    Evals: a one-dimensional array of length nev to hold the
           eigenvalues.
    Evecs: a two-dimensional array of size nev by n to hold the
           eigenvectors.  If this argument is not provided, the
	   eigenvectors are not calculated.  Note that the
	   vectors are stored as columns of Evecs, so that the
	   elements of vector i are the values Evecs[j][i].

  The function is overdefined, so that you can call it in either
  fashion and the appropriate version will be run.

  To use these function calls, there must be a function
  defined in the calling program of the form

    av(int n, double *in, double *out)

  where the function av finds out = A.in, and n is the order of
  the matrix A.  This function must be defined before the
  statement that includes this header file, as it needs to know
  about this function.  It is used in the looping procedure.

  Scot Shaw
  30 August 1999
*/

#define LOCAL_OGF_FORTRAN(x) x##_

// _________________________________________________________________
//
//              ARPACK Fortran routine prototypes
// _________________________________________________________________

typedef int ARint;
typedef int ARlogical;

extern "C" {

// double precision symmetric routines.

    void LOCAL_OGF_FORTRAN(dsaupd)(
        ARint *ido, char *bmat, ARint *n, char *which,
        ARint *nev, double *tol, double *resid,
        ARint *ncv, double *V, ARint *ldv,
        ARint *iparam, ARint *ipntr, double *workd,
        double *workl, ARint *lworkl, ARint *info
    ) ;

    void LOCAL_OGF_FORTRAN(dseupd)(
        ARlogical *rvec, char *HowMny, ARlogical *select,
        double *d, double *Z, ARint *ldz,
        double *sigma, char *bmat, ARint *n,
        char *which, ARint *nev, double *tol,
        double *resid, ARint *ncv, double *V,
        ARint *ldv, ARint *iparam, ARint *ipntr,
        double *workd, double *workl,
        ARint *lworkl, ARint *info
    ) ;

// double precision nonsymmetric routines.

    void LOCAL_OGF_FORTRAN(dnaupd)(
        ARint *ido, char *bmat, ARint *n, char *which,
        ARint *nev, double *tol, double *resid,
        ARint *ncv, double *V, ARint *ldv,
        ARint *iparam, ARint *ipntr, double *workd,
        double *workl, ARint *lworkl, ARint *info
    ) ;

    void LOCAL_OGF_FORTRAN(dneupd)(
        ARlogical *rvec, char *HowMny, ARlogical *select,
        double *dr, double *di, double *Z,
        ARint *ldz, double *sigmar,
        double *sigmai, double *workev,
        char *bmat, ARint *n, char *which,
        ARint *nev, double *tol, double *resid,
        ARint *ncv, double *V, ARint *ldv,
        ARint *iparam, ARint *ipntr,
        double *workd, double *workl,
        ARint *lworkl, ARint *info
    ) ;
}

#include <iostream>

void av(int n, double *in, double *out);

void dsaupd(int n, int nev, double *Evals);

void dsaupd(int n, int nev, double *Evals, double **Evecs);