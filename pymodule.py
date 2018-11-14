#
#@BEGIN LICENSE
#
# erpa by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

import psi4
import re
import os
import math
import warnings
from psi4.driver.procrouting import proc_util
import psi4.driver.p4util as p4util 

#import psi4
#import re
#import os
#import inputparser
#import math
#import warnings
#import driver
#from molutil import *
#import p4util
#from p4util.exceptions import *
#from procedures import *
#from psi4.driver import util
#from psi4.driver import p4util

def run_erpa(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    erpa can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('erpa')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

#    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')
    erpa_corr_global = psi4.core.get_option('ERPA', 'ERPA_CORRELATION')
    erpa_corr_local  = psi4.core.get_global_option('ERPA_CORRELATION')
    if not erpa_corr_global and not erpa_corr_local:
        nat_orbs_global = psi4.core.get_option('V2RDM_CASSCF', 'NAT_ORBS')
        nat_orbs_local  = psi4.core.get_global_option('NAT_ORBS')
        if not nat_orbs_global and not nat_orbs_local:
            raise psi4.ValidationError("""ERPA requires that the V2RDM_CASSCF option "nat_orbs" be set to true.""")


    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    returnvalue = psi4.core.plugin('erpa.so', ref_wfn)

 #   optstash.restore()

    return returnvalue


# Integration with driver routines
#psi4.driver['energy']['erpa'] = run_erpa
psi4.driver.procedures['energy']['erpa'] = run_erpa

def exampleFN():
    # Your Python code goes here
    pass
