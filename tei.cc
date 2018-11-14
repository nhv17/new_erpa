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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include<psi4/libpsio/psio.hpp>
#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libmints/wavefunction.h>
#include<psi4/libmints/mintshelper.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>

// #include<../bin/fnocc/blas.h>
#include<time.h>

#include"erpa_solver.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace std;
using namespace psi;
//using namespace fnocc;

namespace psi{ namespace erpa{

void ERPASolver::GetIntegrals() {

    // one-electron integrals:  
    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    SharedMatrix K1 (new Matrix(mints->so_potential()));
    K1->add(mints->so_kinetic());
    K1->transform(Ca_);

    // size of the tei buffer
    if ( is_df_ ) {

        // size of the 3-index integral buffer
        tei_full_dim_ = (long int) nQ_ * (long int) ( nmo_ - nfrzv_ ) * ( (long int) ( nmo_ - nfrzv_ ) + 1L ) / 2L ;

        // just point to 3-index integral buffer
        tei_full_sym_      = Qmo_;

    }else {

        // size of the 4-index integral buffer
        tei_full_dim_ = 0;
        for (int h = 0; h < nirrep_; h++) {
            tei_full_dim_ += (long int)gems_full[h] * ( (long int)gems_full[h] + 1L ) / 2L;
        }

        tei_full_sym_ = (double*)malloc(tei_full_dim_*sizeof(double));
        memset((void*)tei_full_sym_,'\0',tei_full_dim_*sizeof(double));

    }

    // allocate memory for oei tensor, blocked by symmetry, excluding frozen virtuals
    oei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim_ += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));

    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    if ( !is_df_ ) {
        // read tei's from disk
        GetTEIFromDisk();
    }

    // evaluate frozen core energy
    FrozenCoreEnergy();

}

void ERPASolver::GetTEIFromDisk() {

    double * temptei = (double*)malloc((long int)nmo_*(long int)nmo_*(long int)nmo_*(long int)nmo_*sizeof(double));
    memset((void*)temptei,'\0',(long int)nmo_*(long int)nmo_*(long int)nmo_*(long int)nmo_*sizeof(double));

    // read two-electron integrals from disk
    ReadIntegrals(temptei,(long int)nmo_);
  
    // load tei_full_sym_
    long int offset = 0;
    long int n2 = (long int)nmo_*(long int)nmo_;
    long int n3 = n2 * (long int)nmo_;
    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_full[h]; ij++) {
            long int i = bas_really_full_sym[h][ij][0];
            long int j = bas_really_full_sym[h][ij][1];
            for (long int kl = ij; kl < gems_full[h]; kl++) {
                long int k = bas_really_full_sym[h][kl][0];
                long int l = bas_really_full_sym[h][kl][1];
                tei_full_sym_[offset + INDEX(ij,kl)] = temptei[i*n3+j*n2+k*(long int)nmo_+l];
            }
        }
        offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
    }

    free(temptei);
}

void ERPASolver::FrozenCoreEnergy() {

    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    long int offset = 0;
    long int offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {

            int ifull = i + offset;

            efzc_ += 2.0 * oei_full_sym_[offset3 + INDEX(i,i)];

            long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < rstcpi_[h2] + frzcpi_[h2]; j++) {

                    int jfull = j + offset2;

                    int hij = SymmetryPair(h,h2);

                    double dum1 = TEI(ifull,ifull,jfull,jfull, 0);
                    double dum2 = TEI(ifull,jfull,ifull,jfull, hij);
                    efzc_ += 2.0 * dum1 - dum2;

                }
                offset2 += nmopi_[h2] - frzvpi_[h2];
            }
        }
        offset += nmopi_[h] - frzvpi_[h];
        offset3 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // adjust one-electron integrals for core repulsion contribution
    offset = 0;
    offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {

            int ifull = i + offset;

            for (int j = i; j < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; j++) {

                int jfull = j + offset;

                double dum = 0.0;

                long int offset2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int k = 0; k < rstcpi_[h2] + frzcpi_[h2]; k++) {

                        int kfull = k + offset2;

                        int hik = SymmetryPair(h,h2);

                        double dum1 = TEI(ifull,jfull,kfull,kfull,0);
                        double dum2 = TEI(ifull,kfull,jfull,kfull,hik);
                        dum += 2.0 * dum1 - dum2;

                    }
                    offset2 += nmopi_[h2] - frzvpi_[h2];
                }
                // TODO: verify
                oei_full_sym_[offset3+INDEX(i,j)] += dum;
            }
        }
        offset += nmopi_[h] - frzvpi_[h];
        offset3 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

}

double ERPASolver::TEI(int i, int j, int k, int l, int h) {
    double dum = 0.0;

    if ( is_df_ ) {

        dum = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,j),1,Qmo_+nQ_*INDEX(k,l),1);

    }else {

        int myoff = 0;
        for (int myh = 0; myh < h; myh++) {
            myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
        }

        int ij    = ibas_full_sym[h][i][j];
        int kl    = ibas_full_sym[h][k][l];

        dum = tei_full_sym_[myoff + INDEX(ij,kl)];
    }

    return dum;
}



}}
