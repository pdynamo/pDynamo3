"""Defines a QC/MM electrostatic model compatible with the PySCF program."""

from pCore             import NotInstalledError
from pMolecule.NBModel import QCMMElectrostaticModel
from pScientific       import Units

try:
    import numpy
except:
    raise NotInstalledError ( "numpy not installed." )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelPySCF ( QCMMElectrostaticModel ):

    _classLabel  = "PySCF QC/MM Electrostatic Model"

    def QCMMGradients ( self, target ):
        # . This works for mean-field (R/U)(HF/KS) methods. For post-HF methods, the orbital response needs to be implemented.
        # . Gradients for ROHF with qmmm charges do not work in PySCF (1.7.0).
        # . gradients3B/M are handled by NBModel.PairwiseInteraction* in pyrex/csource.
        # . The gradients are in atomic units!
        if target.scratch.doGradients:
               coords     = target.qcState.mf.mm_mol.atom_coords()
               charges    = target.qcState.mf.mm_mol.atom_charges()
               qm_coords  = target.qcState.mole.atom_coords()
               qm_charges = target.qcState.mole.atom_charges()
               dmi        = target.qcState.mf.make_rdm1()
               if len(dmi.shape) > 2: 
                   dm = dmi[0] + dmi[1] # . Unrestricted/RestrictedOpen-shell DM, so we sum alpha and beta contribs.
               else: 
                   dm = dmi
               # . PySCF does not calculate gradients at the MM centers.
               # . From pyscf/examples/qmmm/30-force_on_mm_particles.py, Author: Qiming Sun <osirpt.sun@gmail.com>
               # . Gradient from the interaction between QM atoms and MM particles
               # \sum_K d/dR (1/|r_K-R|) = \sum_K (r_K-R)/|r_K-R|^3
               dr = qm_coords[:,None,:] - coords
               r  = numpy.linalg.norm(dr, axis=2)
               g  = numpy.einsum('r,R,rRx,rR->Rx', qm_charges, charges, dr, r**-3)

               # . Gradient from the interaction between electron density and MM particles
               # d/dR <i| (1/|r-R|) |j> = <i| d/dR (1/|r-R|) |j> = <i| -d/dr (1/|r-R|) |j>
               #   = <d/dr i| (1/|r-R|) |j> + <i| (1/|r-R|) |d/dr j>
               for i, q in enumerate(charges):
                   with target.qcState.mole.with_rinv_origin(coords[i]):
                       v = target.qcState.mole.intor('int1e_iprinv')
                   f =(numpy.einsum('ij,xji->x', dm, v) +
                       numpy.einsum('ij,xij->x', dm, v.conj())) * -q
                   g[i] += f

               n = len (coords)
               if n > 0:
                  factor      = Units.Length_Angstroms_To_Bohrs * Units.Energy_Hartrees_To_Kilojoules_Per_Mole
                  gradients3B = target.scratch.Get ( "bpGradients3", None )
                  gradients3M = target.scratch.gradients3
                  mmAtoms     = target.mmState.pureMMAtoms
                  nM          = len ( mmAtoms )
                  for i in range ( n ):
                      if i < nM:
                          s = mmAtoms[i]
                          for j in range(3): gradients3M[s,j] = g[i][j] * factor
                      else:
                          s = i-nM
                          for j in range(3): gradients3B[s,j] = g[i][j] * factor

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
