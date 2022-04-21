"""Mathematical and physical unit conversions."""

#===================================================================================================================================
# . Declarations from C.
#===================================================================================================================================
cdef extern from "Units.h":

    cdef float Units_Angle_Degrees_To_Radians
    cdef float Units_Angle_Radians_To_Degrees

    cdef float Units_Dipole_Atomic_Units_To_Debyes

    cdef float Units_Energy_Calories_To_Joules
    cdef float Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole
    cdef float Units_Energy_Electron_Volts_To_Kilojoules_Per_Mole
    cdef float Units_Energy_Hartrees_To_Electron_Volts
    cdef float Units_Energy_Hartrees_To_Kilojoules_Per_Mole
    cdef float Units_Energy_Joules_To_Calories
    cdef float Units_Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
    cdef float Units_Energy_Kilojoules_Per_Mole_To_Kilocalories_Per_Mole

    cdef float Units_Length_Angstroms_To_Bohrs
    cdef float Units_Length_Angstroms_To_Meters
    cdef float Units_Length_Bohrs_To_Angstroms

    cdef float Units_Mass_AMU_To_Kg

    cdef float Units_Pressure_Atmospheres_To_Kilojoules_Per_Mole
    cdef float Units_Pressure_Atmospheres_To_Pascals

    cdef float Units_Quadrupole_Atomic_Units_To_Buckinghams

#===================================================================================================================================
# . Declarations into Python.
#===================================================================================================================================
Angle_Degrees_To_Radians                            = Units_Angle_Degrees_To_Radians
Angle_Radians_To_Degrees                            = Units_Angle_Radians_To_Degrees
Dipole_Atomic_Units_To_Debyes                       = Units_Dipole_Atomic_Units_To_Debyes
Energy_Calories_To_Joules                           = Units_Energy_Calories_To_Joules
Energy_E2Angstroms_To_Kilojoules_Per_Mole           = Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole
Energy_Electron_Volts_To_Kilojoules_Per_Mole        = Units_Energy_Electron_Volts_To_Kilojoules_Per_Mole
Energy_Hartrees_To_Electron_Volts                   = Units_Energy_Hartrees_To_Electron_Volts
Energy_Hartrees_To_Kilojoules_Per_Mole              = Units_Energy_Hartrees_To_Kilojoules_Per_Mole
Energy_Joules_To_Calories                           = Units_Energy_Joules_To_Calories
Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole = Units_Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
Energy_Kilojoules_Per_Mole_To_Kilocalories_Per_Mole = Units_Energy_Kilojoules_Per_Mole_To_Kilocalories_Per_Mole
Length_Angstroms_To_Bohrs                           = Units_Length_Angstroms_To_Bohrs
Length_Angstroms_To_Meters                          = Units_Length_Angstroms_To_Meters
Length_Bohrs_To_Angstroms                           = Units_Length_Bohrs_To_Angstroms
Mass_AMU_To_Kg                                      = Units_Mass_AMU_To_Kg
Pressure_Atmospheres_To_Kilojoules_Per_Mole         = Units_Pressure_Atmospheres_To_Kilojoules_Per_Mole
Pressure_Atmospheres_To_Pascals                     = Units_Pressure_Atmospheres_To_Pascals
Quadrupole_Atomic_Units_To_Buckinghams              = Units_Quadrupole_Atomic_Units_To_Buckinghams
