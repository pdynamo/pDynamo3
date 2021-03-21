"""Definition of simple QC/MM systems for testing."""

import math, os.path

from pBabel                       import ImportSystem
from pCore                        import logFile                  , \
                                         LogFileActive            , \
                                         Selection                , \
                                         TestScript_InputDataPath
from pMolecule                    import EnergyModelPriority
from pMolecule.MMModel            import MMModelOPLS
from pMolecule.NBModel            import NBModelFull
from pMolecule.NBModel.NBDefaults import _NonUpdatablePairLists
from pMolecule.QCModel            import ElectronicState
from pScientific                  import Units
from pScientific.Arrays           import Array                    , \
                                         StorageType

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Paths.
_dataPath = os.path.join ( TestScript_InputDataPath ( "pMolecule" ), "mol" )

# . Molecule definitions.
_keywordLabels = ( "label", "dataPath", "fileName", "fileOutName", "qcCharge", "multiplicity", "qcSelection" )
_moleculeData  = ( ( "Alanine Dipeptide"    , _dataPath, "bAla_c7eq"        , "bAla_c7eq"         ,  0, 1, Selection.FromIterable ( range ( 10, 14 ) ) ) ,
                   ( "Chloride-Chloride"    , _dataPath, "chlorideChloride" , "chlorideChloride"  , -1, 1, Selection.FromIterable ( [ 0 ] ) ) ,
                   ( "Chloride-Water 1"     , _dataPath, "chlorideWater"    , "chlorideWater1"    , -1, 1, Selection.FromIterable ( [ 3 ] ) ) ,
                   ( "Chloride-Water 2"     , _dataPath, "chlorideWater"    , "chlorideWater2"    ,  0, 1, Selection.FromIterable ( [ 0, 1, 2 ] ) ) ,
                   ( "Cyclohexane 6"        , _dataPath, "cyclohexane"      , "cyclohexane6"      ,  0, 1, Selection.FromIterable ( range (  6     ) ) ) ,
                   ( "Cyclohexane 9"        , _dataPath, "cyclohexane"      , "cyclohexane9"      ,  0, 1, Selection.FromIterable ( range (  9     ) ) ) ,
                   ( "Disulfide Bridge"     , _dataPath, "disulfideBridge"  , "disulfideBridge"   ,  0, 1, Selection.FromIterable ( [ 5, 6, 16, 17, 27, 28, 38, 39 ] ) ) ,
                   ( "Tyrosine Dipeptide 1" , _dataPath, "tyrosineDipeptide", "tyrosineDipeptide1",  0, 1, Selection.FromIterable ( list ( range (  6, 14 ) ) + list ( range ( 22, 29 ) ) ) ) ,
                   ( "Tyrosine Dipeptide 2" , _dataPath, "tyrosineDipeptide", "tyrosineDipeptide2",  1, 2, Selection.FromIterable ( list ( range (  6, 14 ) ) + list ( range ( 22, 29 ) ) ) ) ,
                   ( "Water Dimer 1"        , _dataPath, "waterDimer_cs"    , "waterDimer_cs1"    ,  0, 1, Selection.FromIterable ( range (  3     ) ) ) ,
                   ( "Water Dimer 2"        , _dataPath, "waterDimer_cs"    , "waterDimer_cs2"    ,  0, 1, Selection.FromIterable ( range (  3,  6 ) ) ) )

#===================================================================================================================================
# . Class for a QC/MM test system.
#===================================================================================================================================
class QCMMTestSystem:
    """QC/MM test system."""

    def __init__ ( self, **options ):
        """Constructor."""
        for ( attribute, value ) in options.items ( ): setattr ( self, attribute, value )

    # . qcmmModels is a dictionary typically of the form:
    #   { "qcmmElectrostatic" : QCMMElectrostaticModel ( ) ,
    #     "qcmmLennard-Jones  : QCMMLennardJonesModel  ( ) }
    def GetSystem ( self                             ,
                    doQCMM                 = True    ,
                    electronicStateOptions = None    ,
                    log                    = logFile ,
                    nbModel                = None    ,
                    qcmmModels             = None    ,
                    qcModel                = None    ):
        """Get the system with the energy model defined."""
        # . Basic setup.
        molecule       = ImportSystem ( os.path.join ( self.dataPath, self.fileName + ".mol" ) )
        molecule.label = self.label
        molecule.DefineMMModel ( self.mmModel )
        # . Set up the QC model.
        if qcModel is not None:
            esOptions = { "charge" : self.qcCharge, "isSpinRestricted" : ( self.multiplicity == 1 ), "multiplicity" : self.multiplicity }
            if electronicStateOptions is not None: esOptions.update ( electronicStateOptions )
            molecule.electronicState = ElectronicState.WithOptions ( **esOptions )
            if doQCMM: molecule.DefineQCModel ( qcModel, qcSelection = self.qcSelection )
            else:      molecule.DefineQCModel ( qcModel )
        # . Set up the NB model.
        hasQCMMModels = doQCMM and ( qcmmModels is not None ) and ( len ( qcmmModels ) > 0 )
        if ( qcModel is None ) or doQCMM:
            if nbModel is None: nbModel = self.nbModel
            molecule.DefineNBModel ( nbModel, assignQCMMModels = ( not hasQCMMModels ) )
        # . QC/MM models.
        if hasQCMMModels:
            for key in sorted ( qcmmModels.keys ( ) ):
                molecule.AddEnergyModel ( key, qcmmModels[key], priority = EnergyModelPriority.QCMMModel )
        # . Summary.
        if LogFileActive ( log ):
            molecule.Summary ( log = log )
            log.Paragraph ( "Formula = " + molecule.atoms.FormulaString ( ) + "." )
        # . Finish up.
        return molecule

    @property
    def mmModel ( self ): return MMModelOPLS.WithParameterSet ( "protein" )
    @property
    def nbModel ( self ): return NBModelFull.WithDefaults ( )

#===================================================================================================================================
# . Method for testing the QC/MM multipole interaction energy.
#===================================================================================================================================
def TestQCMMMultipoleEnergy ( target, energyTolerance = 1.0e-4, log = logFile ):
    """Test the QC/MM electrostatic energy.

    Target must be a valid QC(MNDO)/MM system with a current energy calculation.
    """
    # . Constants.
    sqrt3 = math.sqrt ( 3.0 )
    # . Basis function indices in pDynamo order.
    S     = 0
    Pz    = 1
    Px    = 2
    Py    = 3
    Dz2   = 4
    Dxz   = 5
    Dyz   = 6
    Dx2y2 = 7
    Dxy   = 8
    # . Nodes.
    qcmmModel = target.qcmmElectrostatic
    qcmmState = target.qcmmState
    scratch   = target.scratch
    pNode     = scratch.Get ( _NonUpdatablePairLists )
    # . Arrays.
    density        = scratch.onePDMP.density
    indices        = target.qcState.orbitalBases.centerFunctionPointers
    nuclearCharges = target.qcState.nuclearCharges
    parameters     = target.qcState.mndoParameters
    qcCoordinates3 = scratch.qcCoordinates3QCMM
    # . Charge separations.
    atomicNumbers     = parameters.atomicNumbers
    uniqueEntries     = parameters.uniqueEntries
    chargeSeparations = { key : value.chargeSeparations for ( key, value ) in uniqueEntries.items ( ) }
    # . Get BP and MM arrays.
    MMSets    = [ ( target.mmState.charges, scratch.Get ( "coordinates3NB", target.coordinates3 ), pNode.Get ( "qcmmElectrostatic", None ) ) ]
    bpCharges = qcmmState.bpCharges
    if bpCharges is not None:
        MMSets.append ( ( bpCharges, scratch.bpCoordinates3, pNode.Get ( "qcbpElectrostatic", None ) ) )
    # . Allocated arrays.
    B    = Array.WithExtent ( 9           , storageType = StorageType.Symmetric )
    fock = Array.WithExtent ( density.rows, storageType = StorageType.Symmetric )
    # . Loop over pairs.
    eElectronic = 0.0
    eNuclear    = 0.0
    fock.Set ( 0.0 )
    for ( mmCharges, mmCoordinates3, pairList ) in MMSets:
        for ( q, p ) in pairList:
            orbitals = indices[q+1] - indices[q]
            i0       = indices[q]
            # . Parameters.
            qM = mmCharges[p]
            dA = dB = qA = qB = qC = 0.0
            if qcmmModel.multipoleOrder > 0:
                ( _, ddp1, ddp2, ddp3, ddp4, ddp5 ) = chargeSeparations[atomicNumbers[q]]
                dA = ddp1
                dB = ddp4
                if qcmmModel.multipoleOrder > 1:
                    qA = 0.5 * ddp2**2
                    qB = 0.5 * ddp3**2
                    qC = 0.5 * ddp5**2
            # . Distance factors.
            dX  = ( mmCoordinates3[p,0] - qcCoordinates3[q,0] ) * Units.Length_Angstroms_To_Bohrs
            dY  = ( mmCoordinates3[p,1] - qcCoordinates3[q,1] ) * Units.Length_Angstroms_To_Bohrs
            dZ  = ( mmCoordinates3[p,2] - qcCoordinates3[q,2] ) * Units.Length_Angstroms_To_Bohrs
            r   = math.sqrt ( dX*dX + dY*dY + dZ*dZ )
            s1  = 1.0 / r
            s2  = s1**2
            s3  = s1**3
            dX /= r
            dY /= r
            dZ /= r
            # . Matrix elements.
            B.Set ( 0.0 )
            B[S,S] = s1
            if orbitals > 1:
                # . Operators.
                Tx  = dX * s2
                Ty  = dY * s2
                Tz  = dZ * s2
                Txx = ( 3.0 * dX * dX - 1.0 ) * s3
                Tyy = ( 3.0 * dY * dY - 1.0 ) * s3
                Tzz = ( 3.0 * dZ * dZ - 1.0 ) * s3
                Txy =   3.0 * dX * dY         * s3
                Txz =   3.0 * dX * dZ         * s3
                Tyz =   3.0 * dY * dZ         * s3
                # . Matrix elements.
                B[Px, S] =      dA * Tx
                B[Py, S] =      dA * Ty
                B[Pz, S] =      dA * Tz
                B[Px,Px] = s1 + qA * Txx
                B[Py,Py] = s1 + qA * Tyy
                B[Pz,Pz] = s1 + qA * Tzz
                B[Px,Py] =      qA * Txy
                B[Pz,Px] =      qA * Txz
                B[Pz,Py] =      qA * Tyz
                if orbitals > 4:
                    B[Dx2y2,    S] =      qB * ( Txx - Tyy ) * 0.5
                    B[Dxz  ,    S] =      qB * Txz
                    B[Dz2  ,    S] =      qB * Tzz * 0.5 * sqrt3
                    B[Dyz  ,    S] =      qB * Tyz
                    B[Dxy  ,    S] =      qB * Txy
                    B[Dx2y2,   Px] =      dB * Tx
                    B[Dx2y2,   Py] =    - dB * Ty
                    B[Dxz  ,   Px] =      dB * Tz
                    B[Dxz  ,   Pz] =      dB * Tx
                    B[Dz2  ,   Px] =    - dB * Tx       / sqrt3
                    B[Dz2  ,   Py] =    - dB * Ty       / sqrt3
                    B[Dz2  ,   Pz] =      dB * Tz * 2.0 / sqrt3
                    B[Dyz  ,   Py] =      dB * Tz
                    B[Dyz  ,   Pz] =      dB * Ty
                    B[Dxy  ,   Px] =      dB * Ty
                    B[Dxy  ,   Py] =      dB * Tx
                    B[Dx2y2,Dx2y2] = s1 - qC * Tzz
                    B[Dxz  ,  Dxz] = s1 - qC * Tyy
                    B[Dz2  ,  Dz2] = s1 + qC * Tzz
                    B[Dyz  ,  Dyz] = s1 - qC * Txx
                    B[Dxy  ,  Dxy] = s1 - qC * Tzz
                    B[Dxz  ,Dx2y2] =      qC * Txz
                    B[Dz2  ,Dx2y2] =    - qC * ( Txx - Tyy ) / sqrt3
                    B[Dyz  ,Dx2y2] =    - qC * Tyz
                    B[Dz2  ,  Dxz] =      qC * Txz / sqrt3
                    B[Dyz  ,  Dxz] =      qC * Txy
                    B[Dxy  ,  Dxz] =      qC * Tyz
                    B[Dyz  ,  Dz2] =      qC * Tyz / sqrt3
                    B[Dxy  ,  Dz2] =    - qC * Txy * 2.0 / sqrt3
                    B[Dxy  ,  Dyz] =      qC * Txz
            # . Energies and Fock.
            B.Scale ( -qM )
            eNuclear += ( qM * nuclearCharges[q] * s1 )
            for u in range ( orbitals ):
                for v in range ( u ):
                    eElectronic     += ( 2.0 * B[u,v] * density[u+i0,v+i0] )
                    fock[u+i0,v+i0] += B[u,v] 
                eElectronic += ( B[u,u] * density[u+i0,u+i0] )
                fock[u+i0,u+i0] += B[u,u]
    # . Check the results.
    eTest1 = ( eElectronic                     + eNuclear ) * Units.Energy_Hartrees_To_Kilojoules_Per_Mole
    eTest2 = ( fock.TraceOfProduct ( density ) + eNuclear ) * Units.Energy_Hartrees_To_Kilojoules_Per_Mole
    eQCMM  = scratch.energyTerms["QC/MM Electrostatic"]
    isOK   = ( ( math.fabs ( eQCMM - eTest1 ) <= energyTolerance ) and ( math.fabs ( eQCMM - eTest2 ) <= energyTolerance ) )
    if ( not isOK ) and LogFileActive ( log ):
        table = log.GetTable ( columns = [ 20, 20, 20 ] )
        table.Start   ( )
        table.Title   ( "QC/MM Test Energy Deviations" )
        table.Heading ( "Test Type"    )
        table.Heading ( "QC/MM Energy" )
        table.Heading ( "Type Energy"  )
        table.Heading ( "Deviation"    )
        table.Entry ( "Energy" )
        table.Entry ( "{:.4f}".format ( eQCMM  ) )
        table.Entry ( "{:.4f}".format ( eTest1 ) )
        table.Entry ( "{:.4f}".format ( math.fabs ( eQCMM - eTest1 ) ) )
        table.Entry ( "Fock" )
        table.Entry ( "{:.4f}".format ( eQCMM  ) )
        table.Entry ( "{:.4f}".format ( eTest2 ) )
        table.Entry ( "{:.4f}".format ( math.fabs ( eQCMM - eTest2 ) ) )
        table.Stop ( )
    return isOK

#===================================================================================================================================
# . Set up the test systems (as a dictionary).
#===================================================================================================================================
qcmmTestSystems = {}
for values in _moleculeData:
    options = { key : value for ( key, value ) in zip ( _keywordLabels, values ) }
    qcmmTestSystems[options["label"]] = QCMMTestSystem ( **options )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
