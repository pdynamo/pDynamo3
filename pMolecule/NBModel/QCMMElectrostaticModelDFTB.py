"""Defines the QC/MM electrostatic model appropriate for the DFTB+ program."""

from  pCore                  import logFile , \
                                    LogFileActive
from  pScientific            import Units
from .QCMMElectrostaticModel import QCMMElectrostaticModel

# . Lattice derivatives?

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelDFTB ( QCMMElectrostaticModel ):

    _classLabel = "DFTB QC/MM Electrostatic Model"

    def QCMMGradients ( self, target ):
        """Process the MM gradients. These should already have been read and stored by the QC model."""
        if target.scratch.doGradients:
            gradientsMM = target.scratch.dftbOutputData["MM Gradients"]
            gradients3B = target.scratch.Get ( "bpGradients3", None )
            gradients3M = target.scratch.gradients3
            mmAtoms     = target.mmState.pureMMAtoms
            nM          = len ( mmAtoms )
            for i in range ( nM ):
                s = mmAtoms[i]
                gradients3M[s,0] += gradientsMM[i,0]
                gradients3M[s,1] += gradientsMM[i,1]
                gradients3M[s,2] += gradientsMM[i,2]
            if gradients3B is not None:
                for i in range ( gradients3B.rows ):
                    gradients3B[i,0] += gradientsMM[i+nM,0]
                    gradients3B[i,1] += gradientsMM[i+nM,1]
                    gradients3B[i,2] += gradientsMM[i+nM,2]

    def QCMMPotentials ( self, target ):
        """Write MM data to an external point-charge file."""
        # . X Y Z Q with coordinates in Angstroms.
        outPath = target.qcState.paths.get ( "PointCharge", None )
        if outPath is not None:
            state         = getattr ( target, self.__class__._stateName )
            chargesB      = state.bpCharges
            chargesM      = target.mmState.charges
            coordinates3B = target.scratch.Get ( "bpCoordinates3", None                )
            coordinates3M = target.scratch.Get ( "coordinates3NB", target.coordinates3 )
            mmAtoms       = target.mmState.pureMMAtoms
            qScale        = 1.0 / self.dielectric
            nM            = len ( mmAtoms  )
            if chargesB is None: nB = 0
            else:                nB = len ( chargesB )
            if ( nB + nM ) > 0:
                pcFile = open ( outPath, "w" )        
                for i in mmAtoms:
                    pcFile.write ( " {:25.15f}{:25.15f}{:25.15f}{:20.10f}\n".format ( coordinates3M[i,0]   ,
                                                                                      coordinates3M[i,1]   ,
                                                                                      coordinates3M[i,2]   ,
                                                                                      qScale * chargesM[i] ) )
                for i in range ( nB )                           :
                    pcFile.write ( " {:25.15f}{:25.15f}{:25.15f}{:20.10f}\n".format ( coordinates3B[i,0]   ,
                                                                                      coordinates3B[i,1]   ,
                                                                                      coordinates3B[i,2]   ,
                                                                                      qScale * chargesB[i] ) )
                pcFile.close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
