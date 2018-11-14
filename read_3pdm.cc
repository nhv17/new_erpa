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
#include <psi4/libqt/qt.h>
#include "erpa_solver.h"

using namespace psi;

namespace psi{namespace erpa{

struct dm3 {
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double val;
};


void ERPASolver::Read3PDM(double * D3aaa, double * D3aab, double * D3bba, double * D3bbb){

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D3AAA) ) throw PsiException("No D3aaa on disk.",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D3AAB) ) throw PsiException("No D3aab on disk.",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D3BBA) ) throw PsiException("No D3bba on disk.",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D3BBB) ) throw PsiException("No D3bbb on disk.",__FILE__,__LINE__);

    //Ca_->print();

    memset((void*)D3aaa,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D3aab,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D3bba,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D3bbb,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));

    psio_address addr_aaa = PSIO_ZERO;
    psio_address addr_bbb = PSIO_ZERO;
    psio_address addr_aab = PSIO_ZERO;
    psio_address addr_bba = PSIO_ZERO;

    // aab
    psio->open(PSIF_V2RDM_D3AAB,PSIO_OPEN_OLD);

    long int naab;
    psio->read_entry(PSIF_V2RDM_D3AAB,"length",(char*)&naab,sizeof(long int));

    for (int count = 0; count < naab; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3AAB,"D3aab",(char*)&d3,sizeof(dm3),addr_aab,&addr_aab);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3aab[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3AAB,1);

    // bba
    psio->open(PSIF_V2RDM_D3BBA,PSIO_OPEN_OLD);

    long int nbba;
    psio->read_entry(PSIF_V2RDM_D3BBA,"length",(char*)&nbba,sizeof(long int));

    for (int count = 0; count < nbba; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3BBA,"D3bba",(char*)&d3,sizeof(dm3),addr_bba,&addr_bba);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3bba[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3BBA,1);

    // aaa
    psio->open(PSIF_V2RDM_D3AAA,PSIO_OPEN_OLD);

    long int naaa;
    psio->read_entry(PSIF_V2RDM_D3AAA,"length",(char*)&naaa,sizeof(long int));

    for (int count = 0; count < naaa; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3AAA,"D3aaa",(char*)&d3,sizeof(dm3),addr_aaa,&addr_aaa);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3aaa[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3AAA,1);

    // bbb
    psio->open(PSIF_V2RDM_D3BBB,PSIO_OPEN_OLD);

    long int nbbb;
    psio->read_entry(PSIF_V2RDM_D3BBB,"length",(char*)&nbbb,sizeof(long int));

    for (int count = 0; count < nbbb; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3BBB,"D3bbb",(char*)&d3,sizeof(dm3),addr_bbb,&addr_bbb);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3bbb[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3BBB,1);

    // check traces:
    /*double traaa = 0.0;
    double trbbb = 0.0;
    double traab = 0.0;
    double trbba = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                traaa += D3aaa[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
                traab += D3aab[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
                trbba += D3bba[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
                trbbb += D3bbb[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
            }
        }
    }
    printf("  tr(d3aaa) = %20.12lf\n",traaa);
    printf("  tr(d3aab) = %20.12lf\n",traab);
    printf("  tr(d3bba) = %20.12lf\n",trbba);
    printf("  tr(d3bbb) = %20.12lf\n",trbbb);*/

}


}} //end namespaces


