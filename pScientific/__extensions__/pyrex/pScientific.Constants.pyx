"""Fundamental mathematical and physical constants."""

#===================================================================================================================================
# . Declarations from C.
#===================================================================================================================================
cdef extern from "Constants.h":

    cdef float Constant_Ln10
    cdef float Constant_Pi

    cdef float Constant_Atomic_Mass
    cdef float Constant_Avogadro_Number
    cdef float Constant_Bohr_Radius
    cdef float Constant_Boltzmann
    cdef float Constant_Electron_Charge
    cdef float Constant_Molar_Gas
    cdef float Constant_Molar_Ideal_Gas_Volume
    cdef float Constant_Planck
    cdef float Constant_Speed_Of_Light
    cdef float Constant_Vacuum_Permittivity

#===================================================================================================================================
# . Declarations into Python.
#===================================================================================================================================
Ln10                   = Constant_Ln10
Pi                     = Constant_Pi

Atomic_Mass            = Constant_Atomic_Mass
Avogadro_Number        = Constant_Avogadro_Number
Bohr_Radius            = Constant_Bohr_Radius
Boltzmann              = Constant_Boltzmann
Electron_Charge        = Constant_Electron_Charge
Molar_Gas              = Constant_Molar_Gas
Molar_Ideal_Gas_Volume = Constant_Molar_Ideal_Gas_Volume
Planck                 = Constant_Planck
Speed_Of_Light         = Constant_Speed_Of_Light
Vacuum_Permittivity    = Constant_Vacuum_Permittivity
