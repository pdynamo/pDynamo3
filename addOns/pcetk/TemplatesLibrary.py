import os, glob, collections

from  pCore         import logFile, LogFileActive, YAMLUnpickle
from .CEModelError  import CEModelError
from .Constants     import YAMLPATHIN
from .ESTFileReader import ESTFileReader

TSite     = collections.namedtuple ("TSite"     , "label  instances  center  atomLabels")
TInstance = collections.namedtuple ("TInstance" , "label  Gmodel  charges  nprotons")


class TemplatesLibrary:
    """A class to represent a library of sites."""

    def __init__ (self, customFiles=None, log=logFile):
        """Constructor."""
        # . Include standard files
        files = glob.glob (os.path.join (YAMLPATHIN, "sites", "*.yaml"))
        # . Include custom files
        if customFiles:
            for filename in customFiles:
                fileTruncated = os.path.basename (filename)
                if os.path.exists (filename):
                    if LogFileActive (log):
                        log.Paragraph ( "Including custom file: {:s}".format ( fileTruncated ) )
                    files.append (filename)
                else:
                    raise CEModelError ( "Custom file {:s} not found.".format ( fileTruncated ) )

                # . Check for a known extension
                stem, extension = os.path.splitext (filename)
                if not extension in ( ".yaml", ".YAML", ".est", ".EST" ):
                    raise CEModelError ( "Unknown extension of file: {:s}".format ( fileTruncated ) )
        # . Parse all files
        self.log   = log
        self.files = files
        self._Parse ()


    def _Parse (self):
        """Parse parameter files."""
        library = []
        for filename in self.files:
            stem, extension = os.path.splitext (filename)

            # . Load a YAML file
            if   extension in (".yaml", ".YAML"):
                yaml    = YAMLUnpickle (filename)
                collect = []
                for template in yaml["instances"]:
                    instance = TInstance (
                        label     =  template["label"]    ,
                        Gmodel    =  template["Gmodel"]   ,
                        charges   =  template["charges"]  ,
                        nprotons  =  template["protons"]  ,
                            )
                    collect.append (instance)
                templateSite = TSite (
                    atomLabels  =  yaml["atoms"]  ,
                    label       =  yaml["site"]   ,
                    center      =  None           ,
                    instances   =  collect        ,
                        )

            # . Load an EST file (extended-MEAD format)
            elif extension in (".est", ".EST"):
                reader    = ESTFileReader.FromPath (filename)
                reader.Parse ()
                collect   = []
                for template in reader.siteInstances:
                    instance = TInstance (
                        label     =  template["label"]    ,
                        Gmodel    =  template["Gmodel"]   ,
                        charges   =  template["charges"]  ,
                        nprotons  =  template["protons"]  ,
                            )
                    collect.append (instance)
                templateSite = TSite (
                    atomLabels  =  reader.siteAtoms   ,
                    label       =  reader.siteLabel   ,
                    center      =  reader.siteCenter  ,
                    instances   =  collect            ,
                        )
            # . Add or replace entry
            found = False
            for index, site in enumerate (library):
                if site.label == templateSite.label:
                    found = True
                    break
            if found:
                # . Entry exists in the library, overwrite it
                library[index] = templateSite
            else:
                # . Add a new entry
                library.append (templateSite)
        # . Finish up
        self.library = library


    def __getitem__ (self, key):
        """Find and return a site from the library."""
        for site in self.library:
            if site.label == key:
                return site
        raise CEModelError ( "Site {:s} not found in the library.".format ( key ) )

    def __contains__ (self, key):
        """Check if a site is in the library."""
        for site in self.library:
            if site.label == key:
                return True
        return False

    def __len__ (self):
        if hasattr (self, "library"):
            return len (self.library)
        return 0

    @property
    def nsites (self):
        if hasattr (self, "library"):
            return len (self.library)
        return 0


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
