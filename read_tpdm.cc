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

#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>

#include "erpa_solver.h"

using namespace psi;

namespace psi{namespace erpa{

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};


// quick and dirty way of reading in CI RDMs ... from ASCII files!!
void ERPASolver::ReadTPDM_CI(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b){

    FILE * fp_aa = fopen("junk.aa","r");
    FILE * fp_bb = fopen("junk.bb","r");
    FILE * fp_ab = fopen("junk.ab","r");

    char * dum = (char*)malloc(9999*sizeof(char));

    // aa
    fscanf(fp_aa,"%s",dum); 
    fscanf(fp_aa,"%s",dum); 
    fscanf(fp_aa,"%s",dum); 
    fscanf(fp_aa,"%s",dum); 
    fscanf(fp_aa,"%s",dum); 

    int naa = 0;
    fscanf(fp_aa,"%i",&naa);

    for (int n = 0; n < naa; n++) {
        int ik,jl;
        double val;
        fscanf(fp_aa,"%i %i %lf\n",&ik,&jl,&val);

        // ik = i * nmo + k;
        int k = ik % nmo_;
        int i = (ik - k) / nmo_;

        // jl = j * nmo + l;
        int l = jl % nmo_;
        int j = (jl - l) / nmo_;

        D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = val;
    }

    // bb
    fscanf(fp_bb,"%s",dum); 
    fscanf(fp_bb,"%s",dum); 
    fscanf(fp_bb,"%s",dum); 
    fscanf(fp_bb,"%s",dum); 
    fscanf(fp_bb,"%s",dum); 

    int nbb = 0;
    fscanf(fp_bb,"%i",&nbb);

    for (int n = 0; n < nbb; n++) {
        int ik,jl;
        double val;
        fscanf(fp_bb,"%i %i %lf\n",&ik,&jl,&val);

        // ik = i * nmo + k;
        int k = ik % nmo_;
        int i = (ik - k) / nmo_;

        // jl = j * nmo + l;
        int l = jl % nmo_;
        int j = (jl - l) / nmo_;

        D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = val;
    }

    // ab
    fscanf(fp_ab,"%s",dum); 
    fscanf(fp_ab,"%s",dum); 
    fscanf(fp_ab,"%s",dum); 
    fscanf(fp_ab,"%s",dum); 
    fscanf(fp_ab,"%s",dum); 

    int nab = 0;
    fscanf(fp_ab,"%i",&nab);

    for (int n = 0; n < nab; n++) {
        int ik,jl;
        double val;
        fscanf(fp_ab,"%i %i %lf\n",&ik,&jl,&val);

        // ik = i * nmo + k;
        int k = ik % nmo_;
        int i = (ik - k) / nmo_;

        // jl = j * nmo + l;
        int l = jl % nmo_;
        int j = (jl - l) / nmo_;

        D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = val;
    }

    free(dum);

}

void ERPASolver::CountNonExternalOrbitals() {

    // number of non external orbitals
    int nne = 0;

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D2AA) ) throw PsiException("No D2aa on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D2BB) ) throw PsiException("No D2bb on disk",__FILE__,__LINE__);


    // determine active orbitals from CASSCF:
    int ref_amo = 0;
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_V2RDM_D2AB,"NUMBER ACTIVE ORBITALS",(char*)&ref_amo,sizeof(int));

    active_reference_orbitals_ = (bool*)malloc(nmo_*sizeof(bool));
    memset((void*)active_reference_orbitals_,'\0',nmo_*sizeof(bool));

    psio_address addr = PSIO_ZERO;
    for (int i = 0; i < ref_amo; i++) {
        int t = 0;
        psio->read(PSIF_V2RDM_D2AB,"ACTIVE ORBITALS",(char*)&t,sizeof(int),addr,&addr);
        active_reference_orbitals_[t] = true;
    }
    psio->close(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    // determine inactive orbitals from CASSCF:
    int ref_inact = 0;
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_V2RDM_D2AB,"NUMBER INACTIVE ORBITALS",(char*)&ref_inact,sizeof(int));

    inactive_reference_orbitals_ = (bool*)malloc(nmo_*sizeof(bool));
    memset((void*)inactive_reference_orbitals_,'\0',nmo_*sizeof(bool));

    addr = PSIO_ZERO;
    for (int i = 0; i < ref_inact; i++) {
        int t = 0;
        psio->read(PSIF_V2RDM_D2AB,"INACTIVE ORBITALS",(char*)&t,sizeof(int),addr,&addr);
        inactive_reference_orbitals_[t] = true;
    }
    psio->close(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    // make new list of non-external orbitals
    full_to_no_ext_ = (int*)malloc(nmo_*sizeof(int));
    for (int i = 0; i < nmo_; i++) {
        full_to_no_ext_[i] = -999;
    }
    nmo_no_ext_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            bool i_external = ( !active_reference_orbitals_[ii] && !inactive_reference_orbitals_[ii] );
            if ( i_external ) continue;
            full_to_no_ext_[ii] = nmo_no_ext_;
            nmo_no_ext_++;
        }
    }
}

void ERPASolver::ReadTPDMLowMemory(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b) {

    int nne = nmo_no_ext_;

    memset((void*)D2aa,'\0',nne*nne*nne*nne*sizeof(double));
    memset((void*)D2bb,'\0',nne*nne*nne*nne*sizeof(double));
    memset((void*)D2ab,'\0',nne*nne*nne*nne*sizeof(double));

    if ( ref_rdm_ == "CI" ) {
        throw PsiException("there is no low-memory option for ref_rdm = CI",__FILE__,__LINE__);
    }else if ( ref_rdm_ == "CCSD" ) {
        throw PsiException("there is no low-memory option for ref_rdm = CCSD",__FILE__,__LINE__);
    }else {
        
        std::shared_ptr<PSIO> psio (new PSIO());

        if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
        if ( !psio->exists(PSIF_V2RDM_D2AA) ) throw PsiException("No D2aa on disk",__FILE__,__LINE__);
        if ( !psio->exists(PSIF_V2RDM_D2BB) ) throw PsiException("No D2bb on disk",__FILE__,__LINE__);

        psio_address addr_aa = PSIO_ZERO;
        psio_address addr_bb = PSIO_ZERO;
        psio_address addr_ab = PSIO_ZERO;

        // ab
        psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

        long int nab;
        psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

        for (int n = 0; n < nab; n++) {
            tpdm d2;
            psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
            int i = full_to_no_ext_[d2.i];
            int j = full_to_no_ext_[d2.j];
            int k = full_to_no_ext_[d2.k];
            int l = full_to_no_ext_[d2.l];
            long int id = i*nne*nne*nne+j*nne*nne+k*nne+l;
            D2ab[id] = d2.val;
        }
        psio->close(PSIF_V2RDM_D2AB,1);

        // aa
        psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

        long int naa;
        psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));

        for (int n = 0; n < naa; n++) {
            tpdm d2;
            psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
            int i = full_to_no_ext_[d2.i];
            int j = full_to_no_ext_[d2.j];
            int k = full_to_no_ext_[d2.k];
            int l = full_to_no_ext_[d2.l];
            long int id = i*nne*nne*nne+j*nne*nne+k*nne+l;
            D2aa[id] = d2.val;
        }
        psio->close(PSIF_V2RDM_D2AA,1);

        // bb
        psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

        long int nbb;
        psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));

        for (int n = 0; n < nbb; n++) {
            tpdm d2;
            psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
            int i = full_to_no_ext_[d2.i];
            int j = full_to_no_ext_[d2.j];
            int k = full_to_no_ext_[d2.k];
            int l = full_to_no_ext_[d2.l];
            long int id = i*nne*nne*nne+j*nne*nne+k*nne+l;
            D2bb[id] = d2.val;
        }
        psio->close(PSIF_V2RDM_D2BB,1);

    }

    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));

    double tra = 0.0;
    double trb = 0.0;

    for (int i = 0; i < nmo_; i++) {

        bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
        if ( i_external ) continue; 
        int ii = full_to_no_ext_[i];

        for (int j = 0; j < nmo_; j++) {

            bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
            if ( j_external ) continue; 
            int jj = full_to_no_ext_[j];

            double duma = 0.0;
            double dumb = 0.0;

            for (int k = 0; k < nmo_; k++) {

                bool k_external = ( !active_reference_orbitals_[k] && !inactive_reference_orbitals_[k] );
                if ( k_external ) continue; 
                int kk = full_to_no_ext_[k];

                duma += D2ab[ii*nne*nne*nne+kk*nne*nne+jj*nne+kk];
                duma += D2aa[ii*nne*nne*nne+kk*nne*nne+jj*nne+kk];

                dumb += D2ab[kk*nne*nne*nne+ii*nne*nne+kk*nne+jj];
                dumb += D2bb[ii*nne*nne*nne+kk*nne*nne+jj*nne+kk];
            }
            D1a[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * duma;
            D1b[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * dumb;

        }

        tra += D1a[i*nmo_+i];
        trb += D1b[i*nmo_+i];

    }
}

void ERPASolver::ReadTPDM(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b) {

    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    if ( ref_rdm_ == "CI" ) {

        ReadTPDM_CI(D2aa,D2bb,D2ab,D1a,D1b);

        // this will only work for closed shells
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                for (int k = 0; k < nmo_; k++) {
                    for (int l = 0; l < nmo_; l++) {
                        int ik = i*nmo_+k;
                        int jl = j*nmo_+l;
                        //double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                        double aa_1 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double aa_2 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double aa_3 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double aa_4 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];

                        double ab_1 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double ab_2 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double ab_3 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double ab_4 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                        double val = 0.5 * (ab_1-ab_2-ab_3+ab_4);
                        double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-val;
                        
                        // print here check against ci
                        //printf("%5i %5i %20.12lf\n",ik,jl,D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
                        //printf("%5i %5i %5i %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",i,j,k,l,ab_1,ab_2,ab_3,ab_4);
                        //if ( fabs(dum) > 1e-6 ) {
                        //    //printf("%5i %5i %20.12lf %20.12lf\n",ik,jl,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l],D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
                        //    printf("%5i %5i %5i %5i %20.12lf %20.12lf\n",i,j,k,l,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l],val);
                        //}
fflush(stdout);
                    }
                }
            }
        }

    }else if ( ref_rdm_ == "CCSD" ) {


        struct iwlbuf Bufab;
        iwl_buf_init(&Bufab,PSIF_MO_TPDM,1e-14,1,1);
        //iwl_buf_init(&Bufab,PSIF_MO_AB_TPDM,0.0,1,1);
        ReadTPDM_IWL(&Bufab,D2ab,nmo_);
        iwl_buf_close(&Bufab,1);

/*
        struct iwlbuf Bufaa;
        iwl_buf_init(&Bufaa,PSIF_MO_AA_TPDM,0.0,1,1);
        ReadTPDM_IWL(&Bufaa,D2aa,nmo_);
        iwl_buf_close(&Bufaa,1);
        C_DSCAL(nmo_*nmo_*nmo_*nmo_,4.0,D2aa,1);

        struct iwlbuf Bufbb;
        iwl_buf_init(&Bufbb,PSIF_MO_BB_TPDM,0.0,1,1);
        ReadTPDM_IWL(&Bufbb,D2bb,nmo_);
        iwl_buf_close(&Bufbb,1);
        C_DSCAL(nmo_*nmo_*nmo_*nmo_,4.0,D2bb,1);
*/

        // symmetrize.  is this necessary?
        for (int i = 0; i < nmo_*nmo_; i++) {
            for (int j = i+1; j < nmo_*nmo_; j++) {
                double dum = D2ab[i*nmo_*nmo_+j] + D2ab[j*nmo_*nmo_+i];
                D2ab[i*nmo_*nmo_+j] = 0.5 * dum;
                D2ab[j*nmo_*nmo_+i] = 0.5 * dum;

                dum = D2aa[i*nmo_*nmo_+j] + D2aa[j*nmo_*nmo_+i];
                D2aa[i*nmo_*nmo_+j] = 0.5 * dum;
                D2aa[j*nmo_*nmo_+i] = 0.5 * dum;

                dum = D2bb[i*nmo_*nmo_+j] + D2bb[j*nmo_*nmo_+i];
                D2bb[i*nmo_*nmo_+j] = 0.5 * dum;
                D2bb[j*nmo_*nmo_+i] = 0.5 * dum;
            }
        }

        // unpack spatial density into its spin blocks
        // this will only work for closed shells
/*
        for (int i = 0; i < doccpi_[0]; i++) {
            for (int b = doccpi_[0]; b < nmo_; b++) {
                for (int a = doccpi_[0]; a < nmo_; a++) {
                    for (int j = 0; j < doccpi_[0]; j++) {
printf("%5i %5i %5i %5i %20.12lf %20.12lf\n",i,b,j,a,D2ab[i*nmo_*nmo_*nmo_+a*nmo_*nmo_+b*nmo_+j],D2ab[i*nmo_*nmo_*nmo_+a*nmo_*nmo_+j*nmo_+b]);
                        double dum = D2ab[i*nmo_*nmo_*nmo_+a*nmo_*nmo_+j*nmo_+b];
                        D2ab[i*nmo_*nmo_*nmo_+a*nmo_*nmo_+b*nmo_+j] = dum;
                        D2ab[b*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+a] = dum;
                        //D2ab[i*nmo_*nmo_*nmo_+a*nmo_*nmo_+j*nmo_+b] = dum;
                        //D2ab[b*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+a] = D2ab[a*nmo_*nmo_*nmo_+i*nmo_*nmo_+b*nmo_+j];
                    }
                }
            }
        }
*/
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                for (int k = 0; k < nmo_; k++) {
                    for (int l = 0; l < nmo_; l++) {
                        int ik = i*nmo_+k;
                        int jl = j*nmo_+l;
                        //double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                        double aa_1 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double aa_2 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double aa_3 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double aa_4 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                        double aa_val = (aa_1-aa_2-aa_3+aa_4);

                        double bb_1 = D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double bb_2 = D2bb[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double bb_3 = D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double bb_4 = D2bb[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                        double bb_val = (bb_1-bb_2-bb_3+bb_4);

                        double ab_1 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double ab_2 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double ab_3 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double ab_4 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                        double ab_val = 0.5 * (ab_1-ab_2-ab_3+ab_4);
                        double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-ab_val;
                        
                        //D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = ab_val;//0.25 * aa_val;
                        //D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = ab_val;//0.25 * bb_val;

                    }
                }
            }
        }

/*
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                for (int k = 0; k < nmo_; k++) {
                    for (int l = 0; l < nmo_; l++) {
                        double dum1 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double dum2 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l] = 1.0 / 3.0 * (dum1 + 2.0 * dum2);
                    }
                }
            }
        }
        C_DCOPY(nmo_*nmo_*nmo_*nmo_,D2aa,1,D2ab,1);
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                for (int k = 0; k < nmo_; k++) {
                    for (int l = 0; l < nmo_; l++) {
                        double dum1 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double dum2 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double dum3 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double dum4 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                        D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = 0.5 * ( dum1 - dum2 - dum3 + dum4 );
                        D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = 0.5 * ( dum1 - dum2 - dum3 + dum4 );
                    }
                }
            }
        }
*/

        //for (int i = 0; i < nmo_; i++) {
        //    for (int j = 0; j < nmo_; j++) {
        //        for (int k = 0; k < nmo_; k++) {
        //            for (int l = 0; l < nmo_; l++) {
        //                int ik = i * nmo_ + k;
        //                int jl = j * nmo_ + l;
        //                printf("%5i %5i %20.12lf\n",ik,jl,D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
        //                //printf("%5i %5i %20.12lf\n",ik,jl,D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
        //                //printf("%5i %5i %20.12lf\n",ik,jl,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
        //            }
        //        }
        //    }
        //}

        //memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
        //memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
        

    }else {
        
        std::shared_ptr<PSIO> psio (new PSIO());

        if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
        if ( !psio->exists(PSIF_V2RDM_D2AA) ) throw PsiException("No D2aa on disk",__FILE__,__LINE__);
        if ( !psio->exists(PSIF_V2RDM_D2BB) ) throw PsiException("No D2bb on disk",__FILE__,__LINE__);

        //Ca_->print();

        psio_address addr_aa = PSIO_ZERO;
        psio_address addr_bb = PSIO_ZERO;
        psio_address addr_ab = PSIO_ZERO;

        // ab
        psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

        long int nab;
        psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

        for (int n = 0; n < nab; n++) {
            tpdm d2;
            psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
            int i = d2.i;
            int j = d2.j;
            int k = d2.k;
            int l = d2.l;
            long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
            D2ab[id] = d2.val;
        }
        psio->close(PSIF_V2RDM_D2AB,1);

        // aa
        psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

        long int naa;
        psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));

        for (int n = 0; n < naa; n++) {
            tpdm d2;
            psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
            int i = d2.i;
            int j = d2.j;
            int k = d2.k;
            int l = d2.l;
            long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
            D2aa[id] = d2.val;
        }
        psio->close(PSIF_V2RDM_D2AA,1);

        // bb
        psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

        long int nbb;
        psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));

        for (int n = 0; n < nbb; n++) {
            tpdm d2;
            psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
            int i = d2.i;
            int j = d2.j;
            int k = d2.k;
            int l = d2.l;
            long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
            D2bb[id] = d2.val;
        }
        psio->close(PSIF_V2RDM_D2BB,1);

        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                for (int k = 0; k < nmo_; k++) {
                    for (int l = 0; l < nmo_; l++) {
                        int ik = i*nmo_+k;
                        int jl = j*nmo_+l;
                        //double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                        double aa_1 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double aa_2 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double aa_3 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double aa_4 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];

                        double ab_1 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double ab_2 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                        double ab_3 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                        double ab_4 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                        double val = 0.5 * (ab_1-ab_2-ab_3+ab_4);
                        double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-val;

                        // print here check against ci
                        //printf("%5i %5i %20.12lf\n",ik,jl,D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
                        //printf("%5i %5i %5i %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",i,j,k,l,ab_1,ab_2,ab_3,ab_4);
                        //if ( fabs(dum) > 1e-6 ) {
                        //    //printf("%5i %5i %20.12lf %20.12lf\n",ik,jl,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l],D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
                        //    printf("%5i %5i %5i %5i %20.12lf %20.12lf\n",i,j,k,l,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l],val);
                        //}
//fflush(stdout);
                    }
                }
            }
        }

        // determine active orbitals from CASSCF:
        int ref_amo = 0;
        psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_V2RDM_D2AB,"NUMBER ACTIVE ORBITALS",(char*)&ref_amo,sizeof(int));

        active_reference_orbitals_ = (bool*)malloc(nmo_*sizeof(bool));
        memset((void*)active_reference_orbitals_,'\0',nmo_*sizeof(bool));

        psio_address addr = PSIO_ZERO;
        for (int i = 0; i < ref_amo; i++) {
            int t = 0;
            psio->read(PSIF_V2RDM_D2AB,"ACTIVE ORBITALS",(char*)&t,sizeof(int),addr,&addr);
            active_reference_orbitals_[t] = true;
        }
        psio->close(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

        // determine inactive orbitals from CASSCF:
        int ref_inact = 0;
        psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_V2RDM_D2AB,"NUMBER INACTIVE ORBITALS",(char*)&ref_inact,sizeof(int));

        inactive_reference_orbitals_ = (bool*)malloc(nmo_*sizeof(bool));
        memset((void*)inactive_reference_orbitals_,'\0',nmo_*sizeof(bool));

        addr = PSIO_ZERO;
        for (int i = 0; i < ref_inact; i++) {
            int t = 0;
            psio->read(PSIF_V2RDM_D2AB,"INACTIVE ORBITALS",(char*)&t,sizeof(int),addr,&addr);
            inactive_reference_orbitals_[t] = true;
        }
        psio->close(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    }

    // check traces:
    double traa = 0.0;
    double trbb = 0.0;
    double trab = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            traa += D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trbb += D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trab += D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }
    //printf("  tr(d2aa) = %20.12lf\n",traa);
    //printf("  tr(d2bb) = %20.12lf\n",trbb);
    //printf("  tr(d2ab) = %20.12lf\n",trab);

    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));

    double tra = 0.0;
    double trb = 0.0;

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {

            double duma = 0.0;
            double dumb = 0.0;
            for (int k = 0; k < nmo_; k++) {
                duma += D2ab[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
                duma += D2aa[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];

                dumb += D2ab[k*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+j];
                dumb += D2bb[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
            }
            D1a[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * duma;
            D1b[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * dumb;

        }

        tra += D1a[i*nmo_+i];
        trb += D1b[i*nmo_+i];

    }

}

void ERPASolver::ReadTPDM_IWL(iwlbuf *Buf,double*d2,int nmo) {

    memset((void*)d2,'\0',nmo*nmo*nmo*nmo*sizeof(double));

    unsigned long int lastbuf;
    Label *lblptr;
    Value *valptr;
    
    unsigned long int idx, p, q, r, s;
    
    lblptr = Buf->labels;
    valptr = Buf->values;
    lastbuf = Buf->lastbuf;
    
    outfile->Printf("\n");
    outfile->Printf("        Read TPDM ......");

    /**
      * first buffer (read in when Buf was initialized)
      */
    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
        p = (unsigned long int) lblptr[idx++];
        r = (unsigned long int) lblptr[idx++];
        q = (unsigned long int) lblptr[idx++];
        s = (unsigned long int) lblptr[idx++];

        double val = (double)valptr[Buf->idx];

        int pr = p*nmo + r;
        int qs = q*nmo + s;

        //printf("hey %5i %5i %5i %5i %5i %5i %20.12lf\n",p,q,r,s,pr,qs,val);
//if ( (p == 0 || p == 1) &&
//     (q == 0 || q == 1) &&
//     (r == 0 || r == 1) &&
//     (s == 0 || s == 1) ) {
//      printf("hey (%5i,%5i) %5i %5i %5i %5i %20.12lf\n",pr,qs,p,q,r,s,val);
//}

        //int pr = p*nmo + r;
        //int qs = q*nmo + s;

        //d2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += val;
        if ( pr != qs ) {
            d2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += val;
            d2[q*nmo*nmo*nmo+p*nmo*nmo+s*nmo+r] += val;
        }else {
            d2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += val;
        }
    }

    /**
      * now do the same for the rest of the buffers
      */
    while(!lastbuf){
        iwl_buf_fetch(Buf);
        lastbuf = Buf->lastbuf;
        for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
            p = (unsigned long int) lblptr[idx++];
            r = (unsigned long int) lblptr[idx++];
            q = (unsigned long int) lblptr[idx++];
            s = (unsigned long int) lblptr[idx++];
    
            double val = (double)valptr[Buf->idx];
            //printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,val);

            int pr = p*nmo + r;
            int qs = q*nmo + s;

//if ( (p == 0 || p == 1) &&
//     (q == 0 || q == 1) &&
//     (r == 0 || r == 1) &&
//     (s == 0 || s == 1) ) {
//        printf("hey (%i,%i) %5i %5i %5i %5i %20.12lf\n",pr,qs,p,q,r,s,val);
//}

            //d2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += val;
            if ( pr != qs ) {
                d2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += val;// * 0.5;
                d2[q*nmo*nmo*nmo+p*nmo*nmo+s*nmo+r] += val;// * 0.5;
            }else {
                d2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += val;
            }

        }
    
    }
    outfile->Printf("done.\n\n");
}

}} //end namespaces


