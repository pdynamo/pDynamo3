"""Example 26 - charge perturbation."""

_Moves           = 100
_VolumeFrequency = 500

from Definitions import *

# . Header.
logFile.Header ( )

# . Set some parameters.
_Lambdas     =    21
_Name        = "example26.ptGeo"
_Solute      =     0
_Temperature = 300.0
_DLambda     =   1.0 / float ( _Lambdas - 1 )
_RT          = Constants.Molar_Gas * _Temperature / 1000.0

# . Read in the system.
solution = ImportSystem ( os.path.join ( pklPath, "water216_cubicBox_mc.pkl" ) )
solution.Summary ( )

# . Define a random number generator.
randomNumberGenerator = RandomNumberGenerator ( )

# . Initialize the dictionary that will hold the free energy values.
dG = {}

# . Perform simulations at different coupling constants.
for i in range ( _Lambdas - 1, -1, -1 ):

    # . Reset the random number generator.
    randomNumberGenerator.SetSeed ( 622199 + i )

    # . Get the value of the coupling parameter.
    Lambda = float ( i ) * _DLambda

    # . Scale the solute's charge parameters.
    MonteCarlo_ScaleIsolateInteractionParameters ( solution, _Solute, chargeScale = Lambda )

    # . Equilibration.
    MonteCarlo_SystemGeometry ( solution                                      ,
                                blocks                = 10                    ,
                                moves                 = _Moves                ,
                                randomNumberGenerator = randomNumberGenerator ,
                                temperature           = _Temperature          ,
                                volumeFrequency       = _VolumeFrequency      )

    # . Data-collection.
    mcData = ExportTrajectory ( os.path.join ( scratchPath, _Name ), solution )
    MonteCarlo_SystemGeometry ( solution                                      ,
                                blocks                = 20                    ,
                                moves                 = _Moves                ,
                                randomNumberGenerator = randomNumberGenerator ,
                                temperature           = _Temperature          ,
                                trajectories          = [ ( mcData, 100 ) ]   ,
                                volumeFrequency       = _VolumeFrequency      )

    # . Define a trajectory object for reading.
    mcData = ImportTrajectory ( os.path.join ( scratchPath, _Name ), solution )
    mcData.ReadHeader ( )

    # . Initialize the accumulators.
    gB = gF = 0.0

    # . Loop over the frames in the trajectory.
    while mcData.RestoreOwnerData ( ):

        # . Get the interaction energy at i.
        MonteCarlo_ScaleIsolateInteractionParameters ( solution, _Solute, chargeScale = Lambda )
        eI = MonteCarlo_IsolateInteractionEnergy ( solution, _Solute )

        # . Calculate the energy at i-1.
        if i > 0:
            MonteCarlo_ScaleIsolateInteractionParameters ( solution, _Solute, chargeScale = Lambda - _DLambda )
            eJ = MonteCarlo_IsolateInteractionEnergy ( solution, _Solute )
            gB += math.exp ( - ( eJ - eI ) / _RT )

        # . Calculate the energy at i+1.
        if i < ( _Lambdas - 1 ):
            MonteCarlo_ScaleIsolateInteractionParameters ( solution, _Solute, chargeScale = Lambda + _DLambda )
            eJ = MonteCarlo_IsolateInteractionEnergy ( solution, _Solute )
            gF += math.exp ( - ( eJ - eI ) / _RT )

    # . Scale and save the values.
    if len ( mcData ) == 0:
        gB = 1.0
        gF = 1.0
    else:
        gB /= float ( len ( mcData ) )
        gF /= float ( len ( mcData ) )
    if i > 0:                dG[(i,i-1)] = - _RT * math.log ( gB )
    if i < ( _Lambdas - 1 ): dG[(i,i+1)] = - _RT * math.log ( gF )

    # . Finish up.
    mcData.ReadFooter ( )
    mcData.Close ( )

# . Output the results.
table = logFile.GetTable ( columns = [ 12, 12, 16, 16, 16 ] )
table.Start   ( )
table.Title   ( "Calculated Free Energies" )
table.Heading ( "Lambda I"     )
table.Heading ( "Lambda J"     )
table.Heading ( "dG (I->J)"    )
table.Heading ( "dG (I<-J)"    )
table.Heading ( "dG (average)" )
dGijTot = dGjiTot = 0.0
for j in range ( _Lambdas - 2, -1, -1 ):
    i = j + 1
    dGij = dG[(i,j)]
    dGji = dG[(j,i)]
    dGijTot += dGij
    dGjiTot += dGji
    table.Entry ( "{:12.4f}".format ( float ( i ) * _DLambda ) )
    table.Entry ( "{:12.4f}".format ( float ( j ) * _DLambda ) )
    table.Entry ( "{:16.4e}".format ( dGij                   ) )
    table.Entry ( "{:16.4e}".format ( dGji                   ) )
    table.Entry ( "{:16.4e}".format ( 0.5 * ( dGij - dGji )  ) )
table.Entry ( "Total:", align = Align.Left, columnSpan = 2 )
table.Entry ( "{:16.3f}".format ( dGijTot ) )
table.Entry ( "{:16.3f}".format ( dGjiTot ) )
table.Entry ( "{:16.3f}".format ( 0.5 * ( dGijTot - dGjiTot ) ) )
table.Stop ( )

# . Save the final conformation.
ExportSystem ( os.path.join ( scratchPath, "Example26Q.pkl" ), solution )

# . Footer.
logFile.Footer ( )
