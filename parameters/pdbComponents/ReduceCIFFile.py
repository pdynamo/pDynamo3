"""Reduce a full PDB component cif file."""

import os, os.path

# . List of components to include in the reduced file.
_ToInclude = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "UNK", "VAL", \
               "HOH", \
               "ATP", "MG", "TPO", \
               "CL", "K", "NA" ]

# . Paths.
_InPath  = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), "fullComponents.cif" )
_OutPath = "components.cif"

# . Process the file.
inFile    = open ( _InPath , "r" )
outFile   = open ( _OutPath, "w" )
writeLine = False
for line in inFile:
    if line.startswith ( "data_" ):
        key = line[5:].strip ( )
        writeLine = ( key in _ToInclude )
    if writeLine: outFile.write ( line )

# . Finish up.
inFile.close  ( )
outFile.close ( )
