/*
 *@BEGIN LICENSE
 *
 * ERPA, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 *
 *@END LICENSE
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<psi4/libmints/writer.h>
#include<psi4/libmints/writer_file_prefix.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libplugin/plugin.h>
#include<psi4/psi4-dec.h>
#include<psi4/liboptions/liboptions.h>

#include<psi4/libpsio/psio.hpp>
#include<psi4/libmints/wavefunction.h>
#include<psi4/psifiles.h>
#include<psi4/libmints/mintshelper.h>

#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>
#include<psi4/libiwl/iwl.h>

#include<time.h>

#include <psi4/libmints/molecule.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/basisset.h>

#include<psi4/libqt/qt.h>

#include "erpa_solver.h"
#include "fortran.h"

#include <psi4/libqt/lapack_intfc_mangle.h>
#include "blas.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace psi;
using namespace fnocc;
using namespace std;

namespace psi { namespace erpa {

void ERPASolver::IPExtendedRPA() {

    double * Aa  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Ab  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Ba  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Bb  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * cc  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * eig = (double*)malloc(nmo_*nmo_*sizeof(double));

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)Aa,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Ab,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Ba,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Bb,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)cc,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)eig,'\0',nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals

    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    long int n = nmo_*nmo_;

    // build two-electron integrals

    // note from AED (9/11/17): since the ERPA code below was designed for use 
    // with 3-index integrals, I had to play some tricks to make it work with
    // 4-index integrals.  If you want to use 4-index integrals (which is a bad
    // idea), the code will perform an SVD of the ERI tensor to yield factorized 
    // integrals that are appropriate for the 3-index code.  The factorization 
    // isn't symmetric, so there are now two buffers for 3-index integrals: Qmo_ 
    // and Qmo2.  In the case of DF, Qmo2 = Qmo_.  I say using the 4-index code
    // is a bad idea because it requires more storage.

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    double * Qmo2;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
        Qmo2 = Qmo_;
    }else {
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        // now, compute SVD of tei_full_sym_ and store result in Qmo_, Qmo2
        double * U  = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * VT = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * S  = (double*)malloc(nn1o2*sizeof(double));
        memset((void*)U,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)VT,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)S,'\0',nn1o2*sizeof(double));

        Qmo_ = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        memset((void*)Qmo_,'\0',nn1o2*nn1o2*sizeof(double));

        int lwork      = 5 * nn1o2;
        double * work  = (double*)malloc(lwork*sizeof(double));

        //info = C_DGESVD('A','A',nn1o2,nn1o2,tei,nn1o2,S,VT,nn1o2,U,nn1o2,work,lwork);
        char jobu  = 'A';
        char jobvt = 'A';
        int info;
        int dim = (int)nn1o2;
        info = F_DGESVD(&jobvt,&jobu,&dim,&dim,tei,&dim,S,U,&dim,VT,&dim,work,&lwork,&info);
        if ( info != 0 ) {
            throw PsiException("something went wrong in the SVD of the 4-index ERI tensor.",__FILE__,__LINE__);
        }

        // repack two-electron integrals (which were destroyed by SVD)
        memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));
        offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        double err = 0.0;
        for (int Q = 0; Q < nn1o2; Q++) {
            for (int kl = 0; kl < nn1o2; kl++) {
                Qmo_[kl*nn1o2 + Q] = U[Q*nn1o2 + kl] * S[Q];
            }
        }
        free(U);
        free(S);
        Qmo2 = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        C_DCOPY(nn1o2*nn1o2,VT,1,Qmo2,1);
        free(VT);
        nQ_ = nn1o2;
        F_DGEMM('t','n',nn1o2,nn1o2,nQ_,1.0,Qmo_,nQ_,Qmo2,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }

    // read TDPM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b,tei);

    // transform everything to NO basis
    //NOTransformation(D1a,D1b,D2aa,D2bb,D2ab,ha,hb,QmoA,QmoB);

    // check energy
    double e2 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2 +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }
    printf("%20.12lf %20.12lf %20.12lf %20.12lf\n",e1,e2,e1+e2,e1+e2+enuc_);
    Process::environment.globals["CURRENT ENERGY"] = e1 + e2 + enuc_;

    // build B(i,j) = <|[ j*, i ]|>
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {
            long int ij = i*nmo_+j;

            Ba[ij] = 2.0 * D1a[j*nmo_+i] - (double)(i==j);
            Bb[ij] = 2.0 * D1b[j*nmo_+i] - (double)(i==j);

        }
    }

    // build A1(i,j) = <|[ j*, [H1, i ] ]|>
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {

            long int ij = i*nmo_+j;

            double duma = h1[i*nmo_+j];
            double dumb = h1[i*nmo_+j];

            for (int q = 0; q < nmo_; q++) {
                duma    -= 2.0 * D1a[j*nmo_+q] * h1[i*nmo_+q];
                dumb    -= 2.0 * D1b[j*nmo_+q] * h1[i*nmo_+q];
            }

            Aa[ij] = duma;
            Ab[ij] = dumb;

        }
    }

    // build A2(i,j) = <| [ j*, [H2, i ] ] |>
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {
            long int ij = i*nmo_+j;

            double duma = 0.0;
            double dumb = 0.0;
            for (long int q = 0; q < nmo_; q++) {
                for (long int s = 0; s < nmo_; s++) {
                    double e1 = tei[INDEX(i,j) * nn1o2 + INDEX(q,s)];//TEI(i,j,q,s,SymmetryPair(symmetry[q],symmetry[s]));
                    double e2 = tei[INDEX(i,s) * nn1o2 + INDEX(q,j)];//TEI(i,s,q,j,SymmetryPair(symmetry[q],symmetry[j]));
                    duma += D1a[q*nmo_+s] * ( e1 - e2 ) + D1b[q*nmo_+s] * e1;
                    dumb += D1b[q*nmo_+s] * ( e1 - e2 ) + D1a[q*nmo_+s] * e1;
                }
            }

            /*for (long int q = 0; q < nmo_; q++) {

                long int jq = j*nmo_+q;
                long int qj = q*nmo_+j;

                for (long int s = 0; s < nmo_; s++) {
                    for (long int r = 0; r < nmo_; r++) {

                        long int rs = r*nmo_+s;
                        long int sr = s*nmo_+r;

                        double eint = TEI(q,s,i,r,SymmetryPair(symmetry[q],symmetry[s]));

                        duma -= D2aa[jq*nmo_*nmo_+rs] * eint;
                        duma -= D2ab[jq*nmo_*nmo_+rs] * eint;

                        dumb -= D2bb[jq*nmo_*nmo_+rs] * eint;
                        dumb -= D2ab[qj*nmo_*nmo_+sr] * eint;

                    }
                }
            }*/
            Aa[ij] += duma;
            Ab[ij] += dumb;
        }
    }

    //for (long int p = 0; p < nmo_; p++) {
    //    for (long int q = p; q < nmo_; q++) {
    //        double dum = A[p*nmo_+q]+A[q*nmo_+p];
    //        A[p*nmo_+q] = A[q*nmo_+p] = 0.5 * dum;

    //        dum = A[p*nmo_+q+nmo_*nmo_]+A[q*nmo_+p+nmo_*nmo_];
    //        A[p*nmo_+q+nmo_*nmo_] = A[q*nmo_+p+nmo_*nmo_] = 0.5 * dum;
    //    }
    //}

    double * newB = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * newA = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)newA,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)newB,'\0',nmo_*nmo_*sizeof(double));

    if ( options_.get_bool("OPTIMIZE_ERPA") ) {
        throw PsiException("erpa optimization not implemented for IP-ERPA",__FILE__,__LINE__);
        if ( !options_["TARGET_STATE"].has_changed() ) {
            throw PsiException("for erpa optimizations, please specify target state symmetries",__FILE__,__LINE__);
        }
    }


    // do this without symmetry first...
    //SymmetricGeneralizedEigenvalueProblem(2*nmo_*nmo_,A,B,cc,eig);
//    GeneralizedEigenvalueProblem(nmo_,Aa,Ba,cc,eig);
    //GeneralizedEigenvalueProblem(nmo_*nmo_,Ab,Bb,cc,eig);

    bool * energy_is_real = (bool*)malloc(nmo_*sizeof(bool));
    for (int h = 0; h < nirrep_; h++) {
        if ( nmopi_[h] == 0 ) continue;
        for (int i = 0; i < nmopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            for (int j = 0; j < nmopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                newA[i*nmopi_[h]+j] = Aa[ii*nmo_+jj];
                newB[i*nmopi_[h]+j] = Ba[ii*nmo_+jj];
            }
        }
        GeneralizedEigenvalueProblem(nmopi_[h],newA,newB,cc,eig,energy_is_real);
    }
    free(energy_is_real);

    /*long int count = 0;
    outfile->Printf("    state");
    outfile->Printf("          energy (Eh)");
    outfile->Printf("       ex energy (Eh)");
    outfile->Printf("       ex energy (eV)\n");
    //outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",0,Process::environment.globals["CURRENT ENERGY"],0.0,0.0);
    for (long int i = nmo_*nmo_-1; i >= 0; i--) {

        eig[i] = 1.0 / eig[i];

        if ( eig[i] > 0.0 && 27.21138 * eig[i] < 50000.0) {

            outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",count+1,eig[i]+Process::environment.globals["CURRENT ENERGY"],eig[i],eig[i] * 27.21138);
            count++;
        }
    }*/

    free(Aa);
    free(Ba);
    free(Ab);
    free(Bb);
    free(newA);
    free(newB);
    free(cc);
    free(eig);

    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);
}


}} // End namespaces

