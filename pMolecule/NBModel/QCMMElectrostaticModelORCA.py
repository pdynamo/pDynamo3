"""Defines the QC/MM electrostatic model appropriate for ORCA."""

from  pCore                  import logFile                , \
                                    LogFileActive
from  pScientific            import Units
from .QCMMElectrostaticModel import QCMMElectrostaticModel

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelORCA ( QCMMElectrostaticModel ):

    # . Defaults.
    _classLabel  = "ORCA QC/MM Electrostatic Model"

    def QCMMGradients ( self, target ):
        """Read MM data from a point-charge gradient file in atomic units and convert to pDynamo units."""
        # . The gradients are in atomic units!
        if target.scratch.doGradients:
            inPath = target.qcState.paths.get ( "PCGrad", None )
            if inPath is not None:
                pcgFile = open ( inPath, "r" ) 
                n       = int ( next ( pcgFile ) )
                if n > 0:
                    factor      = Units.Length_Angstroms_To_Bohrs * Units.Energy_Hartrees_To_Kilojoules_Per_Mole
                    gradients3B = target.scratch.Get ( "bpGradients3", None )   
                    gradients3M = target.scratch.gradients3
                    mmAtoms     = target.mmState.pureMMAtoms
                    nM          = len ( mmAtoms )
                    for i in range ( n ):
                        tokens = next ( pcgFile ).split ( )
                        gX     = float ( tokens[0] ) * factor
                        gY     = float ( tokens[1] ) * factor
                        gZ     = float ( tokens[2] ) * factor
                        if i < nM:
                            s = mmAtoms[i]
                            gradients3M[s,0] += gX
                            gradients3M[s,1] += gY
                            gradients3M[s,2] += gZ
                        else:
                            s = i-nM
                            gradients3B[s,0] += gX
                            gradients3B[s,1] += gY
                            gradients3B[s,2] += gZ
                pcgFile.close ( )

    def QCMMPotentials ( self, target ):
        """Write MM data to an external point-charge file."""
        # . The coordinates are in Angstroms!
        outPath = target.qcState.paths.get ( "PC", None )
        if outPath is not None:
            state         = getattr ( target, self.__class__._stateName )
            chargesB      = state.bpCharges
            chargesM      = target.mmState.charges
            coordinates3B = target.scratch.Get ( "bpCoordinates3", None                )
            coordinates3M = target.scratch.Get ( "coordinates3NB", target.coordinates3 )
            mmAtoms       = target.mmState.pureMMAtoms
            qScale        = 1.0 / self.dielectric
            nM            = len ( mmAtoms )
            if chargesB is None: nB = 0
            else:                nB = len ( chargesB )
            if ( nB + nM ) > 0:
                pcFile = open ( outPath, "w" )        
                pcFile.write ( "{:10d}\n".format ( nB + nM ) )
                for i in mmAtoms:
                    pcFile.write ( " {:10.5f}{:20.10f}{:20.10f}{:20.10f}\n".format ( qScale * chargesM[i] ,
                                                                                     coordinates3M[i,0]   ,
                                                                                     coordinates3M[i,1]   ,
                                                                                     coordinates3M[i,2]   ) )
                for i in range ( nB ):
                    pcFile.write ( " {:10.5f}{:20.10f}{:20.10f}{:20.10f}\n".format ( qScale * chargesB[i] ,
                                                                                     coordinates3B[i,0]   ,
                                                                                     coordinates3B[i,1]   ,
                                                                                     coordinates3B[i,2]   ) )
                pcFile.close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
