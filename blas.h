#ifndef ERPA_BLAS_H
#define ERPA_BLAS_H

#include <psi4/libqt/qt.h>

typedef long int integer;
typedef double doublereal;

extern "C" {

extern int F_DGEEV(char *, char *, int *, double *, int *, double *, double *,
  double *, int *, double *, int *, double *, int *, int *);

extern int F_DGESVD(char *, char *, int *, int *, double *, int *,
  double *, double *, int *, double *, int *, double *, int *, int *);

}

namespace psi{ namespace fnocc {

/**
 * fortran-ordered dgemm
 */
void F_DGEMM(char transa,char transb, integer m, integer n, integer k,
            doublereal alpha,doublereal*A,integer lda,doublereal*B,integer ldb,
            doublereal beta,doublereal*C,integer ldc);

}}


namespace psi{ namespace erpa{

/**
  * diagonalize general matrix and keep eigenvectors
  */
void NonSymmetricEigenvalueEigenvector(long int dim, double * M, double * eigval, double * el, double * er);
/**
  * diagonalize general matrix, don't compute eigenvectors
  */
void GeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig, bool * energy_is_real);
int SymmetricGeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig);

/**
 * name mangling dggev
 */
extern "C" {
    void dggev(char&JOBVL, char&JOBVR, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double * ALPHAR, double * ALPHAI, double * BETA, double * VL, long int &LDVL, double * VR, long int &LDVR,
        double*WORK,long int&LWORK,long int&INFO);
};
inline void DGGEV(char&JOBVL, char&JOBVR, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double * ALPHAR, double * ALPHAI, double * BETA, double * VL, long int &LDVL, double * VR, long int &LDVR,
        double*WORK,long int&LWORK,long int&INFO){
    dggev(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
};



/**
 * name mangling dcopy
 */
extern "C" {
    void dcopy(integer&n,doublereal*dx,integer&incx,doublereal*dy,
         integer&incy);
};
inline void DCOPY(integer&n,doublereal*dx,integer&incx,doublereal*dy,
            integer&incy){
    dcopy(n,dx,incx,dy,incy);
};
/**
 * name mangling dnrm2
 */
extern"C"{
    double dnrm2(integer&N,doublereal*X,integer&INCX);
};
inline double DNRM2(integer&N,doublereal*X,integer&INCX){
    return dnrm2(N,X,INCX);
};
/**
 * name mangling dgesv
 */
extern"C" {
    void dgesv(integer &N,integer &NRHS,doublereal*A,integer &LDA,integer*IPIV,doublereal*B,integer &LDB,integer &INFO);
};
inline void DGESV(integer &N,integer &NRHS,doublereal*A,integer &LDA,integer*IPIV,doublereal*B,integer &LDB,integer &INFO){
    dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO);
};
/**
 * name mangling ddot
 */
extern "C" {
    double ddot(integer&n,doublereal*dx,integer&incx,doublereal*dy,integer&incy);
};
inline double DDOT(integer&n,doublereal*dx,integer&incx,doublereal*dy,integer&incy){
    return ddot(n,dx,incx,dy,incy);
}

/**
 * diagonalize a real symmetric matrix
 */
void Diagonalize(integer N,doublereal*A,doublereal*W);
/**
 * name mangling dsyev
 */
extern "C" {
    void dsyev(char&JOBZ,char&UPLO,integer&N,doublereal*A,integer&LDA,doublereal*W,doublereal*WORK,integer&LWORK,integer&INFO);
};
inline void DSYEV(char&JOBZ,char&UPLO,integer&N,doublereal*A,integer&LDA,doublereal*W,doublereal*WORK,integer&LWORK,integer&INFO){
    dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO);
};
/**
 * diagonalize a real symmetric packed matrix
 */
void Diagonalize2(integer N,doublereal*AP,doublereal*W,doublereal*Z);
/**
 * name mangling dspev
 */
extern "C" {
    void dspev(char&JOBZ,char&UPLO,integer&N,doublereal*AP,doublereal*W,doublereal*Z,integer&LDZ,doublereal*WORK,integer&INFO);
};
inline void DSPEV(char&JOBZ,char&UPLO,integer&N,doublereal*AP,doublereal*W,doublereal*Z,integer&LDZ,doublereal*WORK,integer&INFO){
    dspev(JOBZ,UPLO,N,AP,W,Z,LDZ,WORK,INFO);
};

/**
 *  General SVD
 */
void SVD(integer M,integer N,doublereal*A,doublereal*U,doublereal*VT,doublereal*S);
/**
 * name mangling dgesvd
 */
extern "C" {
    void dgesvd(char&JOBU,char&JOBVT,integer&M,integer&N,doublereal*A,integer&LDA,doublereal*S,doublereal*U,integer&LDU,doublereal*VT,integer&LDVT,doublereal*WORK,integer&LWORK,integer&INFO);
};
inline void DGESVD(char&JOBU,char&JOBVT,integer&M,integer&N,doublereal*A,integer&LDA,doublereal*S,doublereal*U,integer&LDU,doublereal*VT,integer&LDVT,doublereal*WORK,integer&LWORK,integer&INFO){
    dgesvd(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO);
};

}} // end of namespace

#endif
