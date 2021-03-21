# ifndef _UNITS
# define _UNITS

# include "Constants.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unit conversions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Angle. */
# define Units_Angle_Degrees_To_Radians                               ( Constant_Pi / 180.0e+00 ) 
# define Units_Angle_Radians_To_Degrees                               ( 180.0e+00 / Constant_Pi ) 

/* . Dipole. */
# define Units_Dipole_Atomic_Units_To_Debyes                          2.54176568e+00

/* . Energy. */
# define Units_Energy_Calories_To_Joules                              4.184e+00
# define Units_Energy_Joules_To_Calories                              ( 1.0e+00 / 4.184e+00 )
# define Units_Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole    Units_Energy_Calories_To_Joules
# define Units_Energy_Kilojoules_Per_Mole_To_Kilocalories_Per_Mole    Units_Energy_Joules_To_Calories

/* . Conversion from e^2/Angstroms to kJ mol^-1. */
# define Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole              ( ( 1.0e+7  * Constant_Avogadro_Number * Constant_Electron_Charge * Constant_Electron_Charge ) / \
                                                                        ( 4.0e+00 * Constant_Pi * Constant_Vacuum_Permittivity ) )

/* . Should be equivalent to Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole * Units_Length_Angstroms_To_Bohrs. */
# define Units_Energy_Hartrees_To_Kilojoules_Per_Mole                 2625.5e+00

# define Units_Energy_Hartrees_To_Electron_Volts                      27.2113845
# define Units_Energy_Electron_Volts_To_Kilojoules_Per_Mole           ( Units_Energy_Hartrees_To_Kilojoules_Per_Mole / Units_Energy_Hartrees_To_Electron_Volts )

/* . Length. */
# define Units_Length_Angstroms_To_Meters                             1.0e-10

# define Units_Length_Angstroms_To_Bohrs                              ( Units_Length_Angstroms_To_Meters / Constant_Bohr_Radius )
# define Units_Length_Bohrs_To_Angstroms                              ( Constant_Bohr_Radius / Units_Length_Angstroms_To_Meters )

/* . Mass. */
# define Units_Mass_AMU_To_Kg                                         Constant_Atomic_Mass

/* . Pressure. */
# define Units_Pressure_Atmospheres_To_Pascals                        1.013250e+5
# define Units_Pressure_Atmospheres_To_Kilojoules_Per_Mole            ( Units_Pressure_Atmospheres_To_Pascals * Constant_Avogadro_Number * 1.0e-33 )

# endif
