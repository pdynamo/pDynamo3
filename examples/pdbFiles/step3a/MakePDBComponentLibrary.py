"""Make the PDB component library."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions import pdbDataPath
from pBabel      import MakeDefaultPDBComponentLibrary , \
                        PDBComponentVariant
from pCore       import logFile

# . Header.
logFile.Header ( )

# . Make the default library.
library = MakeDefaultPDBComponentLibrary ( fullLibrary = True, libraryPaths = [ pdbDataPath ] )

# . Add a default variant to ATP.
atp = library.GetComponent ( "ATP" )
atp.variants = [ "Fully_Deprotonated" ]

# . Create the ATP variant.
atpVariant = PDBComponentVariant.WithOptions ( componentLabel = "ATP"                              ,
                                               label          = "Fully_Deprotonated"               ,
                                               atomsToDelete  = [ "HOG2", "HOG3", "HOB2", "HOA2" ] ,
                                               formalCharges  = { "O2G" : -1, "O3G" : -1, "O2B" : -1, "O2A" : -1 } )

# . Get Mg.
mg = library.GetComponent ( "MG" )

# . Create the TPO component (THR with a side chain PO3).
tpo = library.GetComponent ( "TPO" )
tpo.leftAtom         = "N"
tpo.leftLink         = "Peptide"
tpo.leftTermination  = "NTerminal"
tpo.isInChain        = True
tpo.isHeteroatom     = True
tpo.rightAtom        = "C"
tpo.rightLink        = "Peptide"
tpo.rightTermination = "CTerminal"
tpo.variants = [ "Fully_Deprotonated" ]

# . Create the TPO variant.
tpoVariant = PDBComponentVariant.WithOptions ( componentLabel = "TPO"                ,
                                               label          = "Fully_Deprotonated" ,
                                               atomsToDelete  = [ "HOP2", "HOP3" ]   ,
                                               formalCharges  = { "O2P" : -1, "O3P" : -1 } )

# . Save all the items.
library.AddItems ( ( atp, mg, tpo, atpVariant, tpoVariant ) )

# . Footer.
logFile.Footer ( )
