/*
 *@BEGIN LICENSE
 *
 * ERPA
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

#include "erpa_solver.h"

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libpsio/psio.hpp>
#include<psi4/libciomr/libciomr.h>

#include <psi4/libiwl/iwl.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include "fortran.h"

//#include "david.h"

INIT_PLUGIN

//using namespace boost;
using namespace psi;

namespace psi{ namespace erpa {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "ERPA"|| options.read_globals()) {

        /*- Do include active/active excitations in the ERPA excited-state wave functions? -*/
        options.add_bool("ACTIVE_ACTIVE_EXCITATIONS", true);

        /*- Do spin adapt the ERPA equations? -*/
        options.add_bool("SPIN_ADAPT_ERPA", false);

        /*- Do compute correlation energy using adiabatic connection? -*/
        options.add_bool("ADIABATIC_CONNECTION", false);

        /*- Do compute correlation energy using ERPA? -*/
        options.add_bool("ERPA_CORRELATION", false);

        /*- Array containing target state symmetries for ERPA optimization. -*/
        options.add("TARGET_STATE", new ArrayType());

        /*- Array containing point-group symmetry of orbitals in energy order.
            Sometimes when degerenate orbitals are present, just sorting the orbitals
            gives a different ordering than DETCI used. -*/
        options.add("SYMMETRY_ENERGY_ORDER", new ArrayType());

        /*- Do optimize excited state energy? -*/
        options.add_bool("OPTIMIZE_ERPA",false);

        /*- Reference RDM type. -*/
        options.add_str("REF_RDM","V2RDM", "V2RDM CI CCSD");
        
        /*- Do reconstruct D3 from D1/D2 for Linear-response? -*/
        options.add_bool("RECONSTRUCT_D3",false);

        /*- Threshold below which elements of the A matrix are pruned -*/
        options.add_double("PRUNING_THRESHOLD", 1e-6);

        /*- Number of expansion coefficients to print -*/
        options.add_int("NUMBER_COEFFICIENTS", 10 );

        /*- Magnitude above which to print excitation coeffients -*/
        options.add_double("COEFFICIENT_THRESHOLD", 0.01);
    }

    return true;
}

extern "C" PSI_API
SharedWavefunction erpa(SharedWavefunction ref_wfn, Options& options)
{

    tstart();

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    ERPA                                             *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Extended random phase approximation              *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Eugene DePrince                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf("\n");

    std::shared_ptr<ERPASolver > erpa (new ERPASolver(ref_wfn,options));

    //outfile->Printf("\n");
    //outfile->Printf("    ==> ERPA expansion: aa/bb <==\n");
    //outfile->Printf("\n");

    // are we doing AC or ERPA excitaiton energies?
    if ( options.get_bool("ADIABATIC_CONNECTION") ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Alpha ERPA integration <==\n");
        outfile->Printf("\n");
        outfile->Printf("    ");
        outfile->Printf("               alpha");
        outfile->Printf("            w(alpha)\n");
       // double e_ac0 = erpa->AlphaERPA(0.0,true,false);
        //e_ac0       += erpa->AlphaERPA(0.0,false,false);
        double e_ac0_s = erpa->DOCIAlphaERPA(0.0,true);
        //e_ac0       += erpa->DOCIAlphaERPA(0.0,false);
        double e_ac0_t = erpa->DOCIAlphaERPA(0.0,false);
     //   double e_ac01_s = erpa->DOCIAlphaERPA(0.1,true);//using linear relationship between w and alpha (10/22/18)
     //   double e_ac01_t = erpa->DOCIAlphaERPA(0.1,false);//10/22/18
 
        double sum = 0.0;
        for (int i = 1; i < 10; i++) {
           // sum += 2.0 * erpa->AlphaERPA(0.1 * i,true,false);
            //sum += 2.0 * erpa->AlphaERPA(0.1 * i,false,false);
              sum += 2.0 * erpa->DOCIAlphaERPA(0.1 * i,true); 
              sum += 2.0 * erpa->DOCIAlphaERPA(0.1 * i,false);
    //        sum += 2.0 * (e_ac01_s - e_ac0_s) * i + e_ac0_s; //this part used to get data based on 2 initial points
    //        sum += 2.0 * (e_ac01_t - e_ac0_t) * i + e_ac0_t; 
           
        }
    //        sum += (e_ac01_s - e_ac0_s) * 10 + e_ac0_s; 
    //        sum += (e_ac01_t - e_ac0_t) * 10 + e_ac0_t; 
       // double e_ac1 = erpa->AlphaERPA(1.0,true,false);
       // e_ac1       += erpa->AlphaERPA(1.0,false);
       double e_ac1_s = erpa->DOCIAlphaERPA(1.0,true);
        //e_ac1       += erpa->DOCIAlphaERPA(1.0,false);
       double e_ac1_t = erpa->DOCIAlphaERPA(1.0,false);
        outfile->Printf("\n");

       double ac_corr  = (e_ac0_s + e_ac1_s + e_ac0_t + e_ac1_t + sum) * 0.1 / 2.0;
       //   double ac_corr  = (e_ac0_s + 2.0 * e_ac01_s + e_ac0_t + 2.0 * e_ac01_t + sum) * 0.1 / 2.0;
              
      // double ac1_corr = (e_ac0 + e_ac1) / 2.0;
        //double ac_corr  = (e_ac1 + sum) * 0.1 / 2.0;
        //double ac1_corr = (e_ac1) / 2.0;

        outfile->Printf("    Reference energy             %20.12lf\n",Process::environment.globals["v2RDM TOTAL ENERGY"]);
        outfile->Printf("\n");
       // outfile->Printf("    * AC1 correlation energy     %20.12lf\n",ac1_corr);
       // outfile->Printf("    * AC1 total energy           %20.12lf\n",Process::environment.globals["v2RDM TOTAL ENERGY"] + ac1_corr);
        outfile->Printf("\n");
        outfile->Printf("    * AC correlation energy      %20.12lf\n",ac_corr);
        outfile->Printf("    * AC total energy            %20.12lf\n",Process::environment.globals["v2RDM TOTAL ENERGY"] + ac_corr);
        outfile->Printf("\n");

    }else if ( options.get_bool("ERPA_CORRELATION") ) {

        outfile->Printf("\n");
        outfile->Printf("    ==> ERPA <==\n");
        outfile->Printf("\n");
        double e_ac1 = erpa->AlphaERPA(1.0,true,false);
        //double e_ac1 = erpa->DOCIAlphaERPA(1.0,true);
        //e_ac1       += erpa->DOCIAlphaERPA(1.0,false);

        outfile->Printf("    Reference energy             %20.12lf\n",Process::environment.globals["v2RDM TOTAL ENERGY"]);
        outfile->Printf("\n");
        outfile->Printf("    * ERPA correlation energy    %20.12lf\n",e_ac1 - Process::environment.globals["v2RDM TOTAL ENERGY"]);
        outfile->Printf("    * ERPA total energy          %20.12lf\n",e_ac1);
        outfile->Printf("\n");

    }else {

        if ( options.get_bool("SPIN_ADAPT_ERPA") ) {
            //erpa->SpinAdaptedExtendedRPA();
            erpa->AlphaERPA(1.0,true,true);
        }else {
            erpa->NewExtendedRPA("AA");
        }
    }


    tstop();
    return ref_wfn;
}

void ERPASolver::FCIDensity(double*D1a,double*D1b,double*D2aa,double*D2bb,double*D2ab){

    double * newaa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * newbb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * newab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * newa  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * newb  = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)newaa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)newbb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)newab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)newa,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)newb,'\0',nmo_*nmo_*sizeof(double));
    ReadDetciDensity(newaa,newbb,newab,newa,newb,nmo_);

    // stupid qt order
    int * qt_to_energy_order = (int*)malloc(nmo_*sizeof(int));
    int * energy_to_qt_order = (int*)malloc(nmo_*sizeof(int));
    memset((void*)qt_to_energy_order,'\0',nmo_*sizeof(int));
    memset((void*)energy_to_qt_order,'\0',nmo_*sizeof(int));

    int ndocc = 0;
    int nsocc = 0;
    for (int h = 0; h < nirrep_; h++){
        ndocc    += doccpi_[h];
        nsocc    += soccpi_[h];
    }

    int count_qt = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < ndocc; i++) {
            if ( symmetry_energy_order[i] == h+1 ) {
                energy_to_qt_order[i] = count_qt;
                qt_to_energy_order[count_qt++] = i;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int i = ndocc; i < ndocc+nsocc; i++) {
            if ( symmetry_energy_order[i] == h+1 ) {
                energy_to_qt_order[i] = count_qt;
                qt_to_energy_order[count_qt++] = i;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int i = ndocc+nsocc; i < nmo_; i++) {
            if ( symmetry_energy_order[i] == h+1 ) {
                energy_to_qt_order[i] = count_qt;
                qt_to_energy_order[count_qt++] = i;
            }
        }
    }
    int count_pitzer = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < doccpi_[h]; i++) {
            count_pitzer++;
        }
        for (int i = doccpi_[h]; i < nmopi_[h]; i++) {
            count_pitzer++;
        }
    }

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {
                    if ( i > j && k < l ) newaa[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = -newaa[(j*nmo_+i)*nmo_*nmo_ + k*nmo_+l];
                    if ( i < j && k > l ) newaa[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = -newaa[(i*nmo_+j)*nmo_*nmo_ + l*nmo_+k];
                    if ( i > j && k > l ) newaa[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] =  newaa[(j*nmo_+i)*nmo_*nmo_ + l*nmo_+k];

                    if ( i > j && k < l ) newbb[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = -newbb[(j*nmo_+i)*nmo_*nmo_ + k*nmo_+l];
                    if ( i < j && k > l ) newbb[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = -newbb[(i*nmo_+j)*nmo_*nmo_ + l*nmo_+k];
                    if ( i > j && k > l ) newbb[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] =  newbb[(j*nmo_+i)*nmo_*nmo_ + l*nmo_+k];

                    long int ie = pitzer_to_energy_order[i];
                    long int je = pitzer_to_energy_order[j];
                    long int ke = pitzer_to_energy_order[k];
                    long int le = pitzer_to_energy_order[l];

                    long int ii = energy_to_qt_order[ie];
                    long int jj = energy_to_qt_order[je];
                    long int kk = energy_to_qt_order[ke];
                    long int ll = energy_to_qt_order[le];

                    D2aa[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = newaa[(ie*nmo_+je)*nmo_*nmo_ + ke*nmo_+le];
                    D2bb[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = newbb[(ie*nmo_+je)*nmo_*nmo_ + ke*nmo_+le];
                    D2ab[(i*nmo_+j)*nmo_*nmo_ + k*nmo_+l] = newab[(ie*nmo_+je)*nmo_*nmo_ + ke*nmo_+le];
                }
            }
        }
    }
    free(qt_to_energy_order);
    free(energy_to_qt_order);
    free(newaa);
    free(newbb);
    free(newab);
    free(newa);
    free(newb);

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            double duma = 0.0;
            double dumb = 0.0;
            for (int k = 0; k < nmo_; k++) {
                duma += D2ab[(i*nmo_+k)*nmo_*nmo_+(j*nmo_+k)];
                dumb += D2ab[(k*nmo_+i)*nmo_*nmo_+(k*nmo_+j)];
            }
            D1a[i*nmo_+j] = duma / nbeta_;
            D1b[i*nmo_+j] = dumb / nalpha_;
        }
    }
}

}} // End namespaces

