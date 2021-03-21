"""PySCF QC/MM potential using the PySCF Model classes."""

import os.path

from pBabel            import ImportSystem
from pCore             import logFile                     , \
                              NotInstalledError           , \
                              Selection                   , \
                              TestScript_InputDataPath    , \
                              TestScriptExit_NotInstalled
from pMolecule         import SystemGeometryObjectiveFunction
from pMolecule.MMModel import MMModelOPLS
from pMolecule.QCModel import ElectronicState

try:
    from addOns.pySCF import NBModelPySCF , \
                             QCModelPySCF
    import pyscf
except:
    TestScriptExit_NotInstalled ( )

# . Start.
logFile.Header ( )

#--------------------------------------------
# . PySCF input.
#--------------------------------------------
method = 'RHF' ; xc=None ; doGradients=True
basis  = '3-21G'
charge = 0
spin   = 0 

#--------------------------------------------
# . Loop over systems.
for ( name, qcSelec, testGradients ) in ( ( "waterDimer_cs", Selection.FromIterable ( [ 0, 1, 2 ]      ), True ) ,
                                          ( "bAla_c7eq"    , Selection.FromIterable ( range ( 10, 14 ) ), True ) ):

    # . Define pDynamo molecule.
    MOL                 = ImportSystem ( os.path.join ( TestScript_InputDataPath ( "pMolecule" ), "mol", "{:s}.mol".format ( name ) ) )
    MOL.electronicState = ElectronicState.WithOptions ( charge = charge, multiplicity = spin + 1 )

    # . Define the MM, NB and QC models.
    mmModel = MMModelOPLS.WithParameterSet ( "protein" ) 
    nbModel = NBModelPySCF.WithDefaults ( )
    qcModel = QCModelPySCF.WithOptions ( deleteJobFiles = True        ,
                                         functional     = xc          ,
                                         method         = method      ,
                                         mf_kwargs      = { 'diis'    : pyscf.scf.ADIIS ( ) }, # . Needs to be redefined ...
                                         mole_kwargs    = { 'verbose' : 0 } ,
                                         orbitalBasis   = basis       )

    # . Define the energy model.
    MOL.DefineMMModel ( mmModel )
    MOL.DefineQCModel ( qcModel , qcSelection = qcSelec )
    MOL.DefineNBModel ( nbModel ) 
    MOL.Summary ( )

    # . Calculate an energy and gradients.
    MOL.Energy ( doGradients = doGradients )      
    if testGradients:
        of = SystemGeometryObjectiveFunction.FromSystem ( MOL )
        of.TestGradients ( )

# . Stop.
logFile.Footer ( )
