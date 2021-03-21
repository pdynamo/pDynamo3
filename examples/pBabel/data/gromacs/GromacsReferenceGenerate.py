"""Write reference yaml files for Gromacs data."""

import os

from pCore import logFile       , \
                  LogFileActive , \
                  TestDataSet   , \
                  TestReal      , \
                  YAMLMappingFile_FromObject

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
#_SystemLabels = [ "ava", "1atp_peptide", "crambin" ] 
_SystemLabels = [ "ava", "1atp_peptide" ] 

# . Force field names.
#_ForceFields  = [ "CHARMM", "AMBER", "OPLS" ]
_ForceFields  = [ "CHARMM", "AMBER" ]

# . Options.
_AbsoluteErrorTolerance         = 0.5
_GradientAbsoluteErrorTolerance = 1.0e-03

# . Reference values.
pDynamoObserved = { \
      'ava-AMBER'          : {'Fourier Dihedral': 44.670834719049054, 'Potential Energy': -101.28190740795, 'MM/MM Lennard-Jones': -14.064513428536756, 'MM/MM Electrostatic': -1050.0422912611261, 'Harmonic Bond': 20.775629951402532, 'Fourier Improper': 0.019987121794708691, 'MM/MM 1-4 Electrostatic': 840.04593339686073, 'MM/MM 1-4 Lennard-Jones': 43.894106585814619, 'Gradient Error': 1.2700677785915104e-05, 'Harmonic Angle': 13.418405506791224}, \
      'ava-CHARMM'         : {'Fourier Dihedral': 13.72169044570512, 'Harmonic Improper': 0.037322600424997687, 'MM/MM Lennard-Jones': -12.581979424742224, 'Potential Energy': 107.90881710249582, 'Urey-Bradley': 0.92423796871695385, 'MM/MM Electrostatic': -1017.1123629775817, 'Harmonic Bond': 3.2071943938419158, 'MM/MM 1-4 Electrostatic': 1113.9924841186792, 'MM/MM 1-4 Lennard-Jones': 11.627072056416131, 'CMAP Dihedral': -16.421353658193567, 'Gradient Error': 1.3989508545364515e-05, 'Harmonic Angle': 10.514511579229001}, \
      '1atp_peptide-AMBER' : {'Fourier Dihedral': 659.62440280233182, 'Potential Energy': 2815.8698774217401, 'MM/MM Lennard-Jones': 848.57104030116125, 'MM/MM Electrostatic': -3440.7572268824629, 'Harmonic Bond': 2246.4445046382393, 'Fourier Improper': 1.4496707406907781, 'MM/MM 1-4 Electrostatic': 1469.475305857085, 'MM/MM 1-4 Lennard-Jones': 539.66062169405723, 'Harmonic Angle': 491.40155827063762}, \
      '1atp_peptide-CHARMM': {'Fourier Dihedral': 311.01789688872691, 'Harmonic Improper': 2.0763471355852237, 'MM/MM Lennard-Jones': 831.67641502232198, 'Potential Energy': 2758.554213523988, 'Urey-Bradley': 153.13246782399202, 'MM/MM Electrostatic': -2060.6027695430394, 'Harmonic Bond': 2566.5962344099862, 'MM/MM 1-4 Electrostatic': 119.79882361641397, 'MM/MM 1-4 Lennard-Jones': 478.20597361491809, 'CMAP Dihedral': -73.974179211310783, 'Harmonic Angle': 430.62700376639339}, \
                  }

# . U-B term in Gromacs includes harmonic angle term, but not in pDynamo. Sum should be equal.
gromacsObserved = { \
     'ava-AMBER'          : {'Fourier Dihedral': 44.671059, 'Fourier Improper': 0.019947, 'MM/MM Lennard-Jones': -14.064210, 'Potential Energy': -101.319580, 'MM/MM Electrostatic': -1050.042969, 'Harmonic Bond': 20.774445, 'MM/MM 1-4 Electrostatic': 840.012268, 'MM/MM 1-4 Lennard-Jones': 43.893955, 'Gradient Error': 1.4035588240091101e-05, 'Harmonic Angle': 13.415921}, \
     'ava-CHARMM'         : {'Fourier Dihedral': 13.721725, 'Harmonic Improper': 0.037332, 'MM/MM Lennard-Jones': -12.581890, 'Potential Energy': 107.907959, 'MM/MM Electrostatic': -1017.112183, 'Harmonic Bond': 3.207108, 'MM/MM 1-4 Electrostatic': 1113.991455, 'MM/MM 1-4 Lennard-Jones': 11.627161, 'CMAP Dihedral': -16.421263, 'Gradient Error': 1.4035588240091101e-05, 'Harmonic Angle + U-B': 11.438545}, \
     '1atp_peptide-AMBER' : {'Fourier Dihedral': 659.624268, 'Fourier Improper': 1.449597, 'MM/MM Lennard-Jones': 848.560913, 'Potential Energy': 2815.853516, 'MM/MM Electrostatic': -3440.765137, 'Harmonic Bond': 2246.516357, 'MM/MM 1-4 Electrostatic': 1469.416626, 'MM/MM 1-4 Lennard-Jones': 539.657837, 'Harmonic Angle': 491.392853}, \
     '1atp_peptide-CHARMM': {'Fourier Dihedral': 311.017822, 'Harmonic Improper': 2.076368, 'MM/MM Lennard-Jones': 831.665039, 'Potential Energy': 2758.623047, 'MM/MM Electrostatic': -2060.593994, 'Harmonic Bond': 2566.679688, 'MM/MM 1-4 Electrostatic': 119.798035, 'MM/MM 1-4 Lennard-Jones': 478.203125, 'CMAP Dihedral': -73.974022, 'Harmonic Angle + U-B': 583.750916}, \
     'crambin-AMBER'      : {'Fourier Dihedral': 1589.967041, 'Fourier Improper':8.644552, 'MM/MM Lennard-Jones': -1170.950806, 'Potential Energy': -3259.375732, 'MM/MM Electrostatic': -13067.493164, 'Harmonic Bond': 502.356781, 'MM/MM 1-4 Electrostatic': 7718.257812, 'MM/MM 1-4 Lennard-Jones': 771.424194, 'Harmonic Angle': 388.418030}, \
     'crambin-CHARMM'     : {'Fourier Dihedral': 0.0, 'Harmonic Improper': 0.0, 'MM/MM Lennard-Jones': 0.0, 'Potential Energy': 0.0, 'MM/MM Electrostatic': 0.0, 'Harmonic Bond': 0.0, 'MM/MM 1-4 Electrostatic': 0.0, 'MM/MM 1-4 Lennard-Jones': 0.0, 'Harmonic Angle': 0.0}, \
                  }

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Initialization.
gromacsReference = True
#
if gromacsReference: 
    referenceDataPath = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "examples", "pBabel", "data", "reference", "GromacsTopCrdRead_gromacsValues.yaml" )
    localObserved = gromacsObserved
else               :
    referenceDataPath = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "examples", "pBabel", "data", "reference", "GromacsTopCrdRead_pDynamoValues.yaml" )
    localObserved = pDynamoObserved
referenceData = TestDataSet.WithOptions ( label = "Gromacs Test" )

# . Loop over systems and force-fields.
for label in _SystemLabels:
    for ff in _ForceFields:

        # . Generate reference data.
        dataLabel = label + "-" + ff
        localData = TestDataSet.WithOptions ( label = dataLabel, parent = referenceData )
        for ( key, value ) in localObserved[dataLabel].items ( ):
            if key == "Gradient Error": localData.AddDatum ( TestReal.WithOptions ( label = key, value = value, parent = localData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance ) )
            else:                       localData.AddDatum ( TestReal.WithOptions ( label = key, value = value, parent = localData, absoluteErrorTolerance = _AbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
        referenceData.AddDatum ( localData )

# . Write reference data file.
YAMLMappingFile_FromObject ( referenceDataPath, "Gromacs Test", referenceData )
