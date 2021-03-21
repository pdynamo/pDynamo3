"""The PySCF QC model."""

import glob, os

try:
    import pyscf
    _PYSCFFound = True
except:
    _PYSCFFound = False

from pCore             import NotInstalledError
from pMolecule.QCModel import QCModel           , \
                              QCModelError
from pScientific       import Units

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultFunctional = "blyp"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelPySCF ( QCModel ):
    """The PySCF QC model class."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "PySCF QC Model"
    _summarizable = dict ( QCModel._summarizable )
    _attributable.update ( { "deleteJobFiles" : True    ,
                             "functional"     : "blyp"  ,
                             "method"         : "RHF"   , 
                             "mf"             : None    , 
                             "mf_kwargs"      : dict    ,
                             "mole"           : None    ,
                             "mole_kwargs"    : dict    , 
                             "orbitalBasis"   : "3-21G" ,
                             "pyscf"          : None    ,
                             "pySCFscratch"   : None    } )
    _summarizable.update ( { "functional"     : "Functional"    ,
                             "method"         : "Method"        ,
                             "orbitalBasis"   : "Orbital Basis" } )

    def __del__ ( self ):
        """Deallocation."""
        self.DeleteJobFiles ( )

    def _CheckOptions ( self ):
        """Check options."""
        super ( QCModelPySCF, self )._CheckOptions ( )
        # . pySCF.
        if _PYSCFFound: self.pyscf = pyscf
        else: raise NotInstalledError ( "pySCF not installed." )
        # . pySCFscratch.
        scratch = os.getenv ( 'PYSCF_TMPDIR' )
        if scratch is None:
            scratch                    = os.getenv ( "PDYNAMO3_SCRATCH" )
            os.environ['PYSCF_TMPDIR'] = str ( scratch )
        self.pySCFscratch = scratch
        # . Functional.
        if ( 'KS' in self.method.upper ( ) ) and ( self.functional is None ):
            self.functional = _DefaultFunctional

    def CreateMole ( self, target, doQCMM ):
        """Create PySCF mole and mean-field objects"""
        state        = getattr ( target, self.__class__._stateName )
        n            = len ( state.atomicNumbers )
        coordinates3 = target.scratch.qcCoordinates3AU
        mole         = self.pyscf.gto.Mole()
        mole.atom    = [[state.atomicNumbers[i], (coordinates3[i][0], coordinates3[i][1], coordinates3[i][2])] for i in range(n)]
        mole.basis   = self.orbitalBasis
        mole.charge  = target.electronicState.charge
        mole.spin    = target.electronicState.multiplicity - 1 # 2S
        mole.unit    = 'Bohr'
        mole.verbose = 0
        mole.__dict__.update(self.mole_kwargs)
        state.mole = mole.build()
        state.mf   = state.mole.apply ( self.method, **self.mf_kwargs )
        if doQCMM:
            charges       = []
            chargesB      = getattr ( target.qcmmState, "bpCharges", None )
            chargesM      = target.mmState.charges
            coords        = []
            coordinates3B = target.scratch.Get ( "bpCoordinates3", None                ) 
            coordinates3M = target.scratch.Get ( "coordinates3NB", target.coordinates3 )
            mmAtoms       = target.mmState.pureMMAtoms
            qScale        = 1.0 / target.qcmmElectrostatic.dielectric
            nM            = len ( mmAtoms )
            if chargesB is None: nB = 0
            else:                nB = len ( chargesB )
            for i in mmAtoms:
               charges.append(qScale * chargesM[i])
               coords.append ([ coordinates3M[i,0], coordinates3M[i,1], coordinates3M[i,2] ])
            for i in range ( nB ):
               charges.append(qScale * chargesB[i])
               coords.append ([ coordinates3B[i,0], coordinates3B[i,1], coordinates3B[i,2] ])
            # . Coordinates for MM centers stored in Angstroms.
            state.mf = self.pyscf.qmmm.add_mm_charges(state.mf, coords, charges, unit='Angstrom')

    def DeleteJobFiles ( self ):
        """Delete job files."""
        if self.deleteJobFiles:
            try:
                jobFiles = glob.glob ( os.path.join ( self.pySCFscratch, "tmp????????" ) )
                for jobFile in jobFiles: os.remove ( jobFile )
            except:
                pass

    def Energy ( self, target ):
        """Calculate the quantum chemical energy and gradient."""
        state           = getattr ( target, self.__class__._stateName )
        xc              = None
        if 'KS' in self.method.upper(): xc = self.functional
        doGradients     = target.scratch.doGradients
        doQCMM          = ( len ( target.atoms ) > len ( state.qcAtoms ) )
        try: 
            self.CreateMole ( target, doQCMM )
        except:
            raise QCModelError ( "Error creating PySCF mole and mean-field objects." )

        # . Calculate energy.
        try:     
            state.mf.run(xc=xc)
        except:
           raise QCModelError ( "Error calculating PySCF energy." )

        if not state.mf.converged: raise QCModelError ( "SCF energy calculation did not converge." )
        target.scratch.energyTerms["PySCF QC"] = ( state.mf.energy_tot() * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        if doGradients:
            try: 
                g  = state.mf.nuc_grad_method()
                ga = g.kernel()
            except:
                raise QCModelError ( "Error calculating PySCF gradient." )  
            # . QC gradients
            for i in range ( len ( state.atomicNumbers ) ):
                for j in range ( 3 ):
                   target.scratch.qcGradients3AU[i,j] = ga[i,j]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

