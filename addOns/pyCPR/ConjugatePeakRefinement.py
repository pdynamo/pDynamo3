"""Python-based Conjugate peak refinement (PyCPR) module."""

import math, os, string

from  sys                                    import exit
from  time                                   import strftime
from  pBabel                                 import ExportTrajectory
from  pCore                                  import Clone                                     , \
                                                    logFile                                   , \
                                                    LogFileActive                             , \
                                                    PickleFileExtension
from  pMolecule                              import MultiLayerSystemGeometryObjectiveFunction , \
                                                    SystemGeometryObjectiveFunction
from  pScientific.Arrays                     import Array
from  pScientific.ObjectiveFunctionIterators import LBFGSMinimizer
from .CPRSaddlePointRefinement               import CPRSaddlePointRefinement

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Default restart file name.
_DefaultRestartFile = "cpr.restart"

#===================================================================================================================================
# . CPR objective function.
#===================================================================================================================================
class CPRObjectiveFunction ( MultiLayerSystemGeometryObjectiveFunction ):
    """The CPR objective function."""

    _attributable = dict ( MultiLayerSystemGeometryObjectiveFunction._attributable )
    _attributable.update ( { "imageTrajectory" : None ,
                             "numberOfImages"  : 0    } )

    def ActivateImage (self, iImage):
        """Activate an image."""
        self.imageTrajectory.RestoreOwnerData (index=iImage)
        # . A fudge.
        if self.freeAtoms is None: self.iCoordinates3 = self.system.coordinates3.iterator
        else:                      self.iCoordinates3 = self.system.coordinates3.RowIterator ( selection = self.freeAtoms )

    def InitializeImages (self, imageTrajectory):
        """Initialize image data, given a trajectory."""
        self.imageTrajectory = imageTrajectory
        self.numberOfImages  = len (imageTrajectory)

#===================================================================================================================================
# . CPR Refinement.
#===================================================================================================================================
class CPRRefinement:
    """CPR refinement."""

    # . Option attributes.
    optionAttributes = { "finiteDifferenceStep"          : 1.0e-3,
                         "finiteDifference"              : True  ,
                         "lineSearchInitialStep"         : 0.1   ,
                         "interpolationIncreasement"     :   3   ,
                         "rmsdMaximumTolerance"          :  1.5  ,
                         "inputRestartFile"              : None  ,
                         "outputRestartFile"             : None  ,
                         "stepsPerSegment"               :   3   ,
                         "upperStepLimit"                : None  ,
                         "lowestStepLimit"               : None  ,
                         "maximalMinimizationSteps"      : None  ,
                         "maxCPRrun"                     : 1000  ,
                         "numberOfSuccessiveLowGradients": None  ,
                         "initialOrthogonalRuns"         :   0   ,
                         "rmsGradientTolerance"          :  0.1  ,
                         "rmsdLoopPreventionX0"          : 1.0e-6,  # . bigger values can cause problems because different points are detected to be the same!!
                         "rmsdLoopPreventionS0"          : 1.0e-6,  # . bigger values can cause problems because different points are detected to be the same!!
                         "tauTolAdd"                     : 0.15  ,
                         "tauTolRefine"                  : 0.05  ,
                         "increaseTau"                   : False ,
                         "breakIfTauReached"             : True  ,
                         "finalUnrefinableRefinement"    : False ,
                         "scratchDirectory"              :  None ,
                         "betaType"                      : "PRP" }  # FR: Fletcher-Reeves; PRP: Polak-Ribiere and Polyak; HS: Hestenes-Stiefel  

    # . Option attribute names.
    #optionAttributeNames = { "Maximal Minimization Steps"     : "maximalMinimizationSteps"    ,
    #                         "RMS Gradient Tolerance" : "rmsGradientTolerance" }

    # . Objective function attributes.
    objectiveFunctionAttributes = { "CPRrun"                     :  0    ,
                                    "ds0"                        : None  ,
                                    "error"                      : None  ,
                                    "exitMessage"                : None  ,
                                    "energy"                     : None  ,
                                    "segmentsHighest"            : None  ,
                                    "segmentsLowest"             : None  ,
                                    "nFraction"                  :   6   ,
                                    "iImage"                     :  -1   ,
                                    "log"                        : None  ,
                                    "numberOfImages"             : None  ,
                                    "numberOfVariables"          : 0     ,
                                    "numberOfSaddle"             : 0     ,
                                    "numberOfUnrefinables"       : 0     ,
                                    "numberOfStationary"         : 0     ,
                                    "numberOfMaximaInSegment"    : 0     ,
                                    "increasedTau"               : None  ,
                                    "changeInterpolations"       : None  ,
                                    "objectiveFunction"          : None  ,
                                    "optimizer"                  : None  ,
                                    "optimizerState"             : None  ,
                                    "orthogonal"                 : False ,
                                    "isConverged"                : False ,
                                    "isSaddle"                   : None  ,
                                    "isStationary"               : None  ,
                                    "isUnrefinable"              : None  ,
                                    "finalUnrefRef"              : False ,
                                    "fullyRefinedPath"           : False ,
                                    "tangent"                    : None  ,
                                    "variables"                  : None  ,
                                    "listOfAddedStructures"      : None  }

    # . Public methods.
    def __init__ (self, **options):
        """Constructor."""
        # . Set defaults for all options.
        keyvalues = self.__class__.optionAttributes.items ( )
        for (key, value) in keyvalues: setattr (self, key, value)
        # . Set passed in parameters.
        unknowns = set ()
        for (key, value) in options.items ( ):
            if key in self.__class__.optionAttributes: setattr (self, key, value)
            else: unknowns.add (key)
        if len (unknowns) > 0: raise ValueError ("Invalid options: " + ", ".join (sorted (unknowns)) + ".")
        self._CheckOptions ( )

    def _CheckOptions ( self ):
        """Check options."""
        # . Set up the output restart file.
        # . A fudge to ensure that the restart file is put in scratch and not the current directory if no name is given.
        if self.outputRestartFile is None:
            self.outputRestartFile = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), _DefaultRestartFile )

    def Continue (self):
        """Check to see if iterations should continue."""
        if self.numberOfImages == 2:  # . Special initial case with no intermediate image. Has to be treated separately.
             maxunrefinable = 1
        else:  # . All path points except terminal points should become saddle points.
            maxunrefinable = self.numberOfImages - 2  
        notFinished = (self.error is None)                                                                \
                      and (self.CPRrun < self.maxCPRrun)                                                  \
                      and (not self.fullyRefinedPath)                                                     \
                      and ((self.numberOfSaddle + self.numberOfUnrefinables + self.numberOfStationary) < maxunrefinable or self.numberOfMaximaInSegment > 0)       
        if not (notFinished):
            if not (self.error is None):
                self.exitMessage = self.error
            elif self.fullyRefinedPath:
                self.exitMessage = "No further maxima within the segments or along the path. CPR finished successfully."
            elif not (self.CPRrun < self.maxCPRrun):
                self.exitMessage = "Maximum number of CPR runs ({:d}) reached".format(self.maxCPRrun)
            elif not ((self.numberOfSaddle + self.numberOfUnrefinables + self.numberOfStationary) < maxunrefinable or self.numberOfMaximaInSegment > 0):
                self.exitMessage = "All path points are not further refinable and no higher points are in the segments. CPR finished successfully."

        return notFinished

    def CountInSegment (self):
        """Count how many inSegment images are there."""
        self.numberOfMaximaInSegment = 0
        for i in range (0, self.numberOfImages - 1):
            if self.segmentsHighest[i][1] == "inSegment": self.numberOfMaximaInSegment += 1

    def CPRIterate ( self, objectiveFunction, system, imageTrajectory, log = logFile ):
        """Apply the algorithm to a function."""
        trajectory = imageTrajectory
        self.Header ( log = log )
        try:
            self.Initialize (objectiveFunction)
            self.DetermineStepSize (system)
            if self.inputRestartFile is not None:  # . Start from restart file.
               self.ReadInRestartFile (system, trajectory)
            else:  # . With no restart file all segments have to be scanned.
               for i in range (0, self.numberOfImages - 1):
                  self.SegmentInterpolation (system, trajectory, i)

            self.PathSummary ("initial")
            self.CountInSegment ()  # . Count inSegment images for self.Continue().

            while (self.Continue ()):
              if self.CPRrun < self.initialOrthogonalRuns:
                   self.optimizer.orthogonal = True
              else:
                   self.optimizer.orthogonal = self.orthogonal

              # . Find segment containing highest point and decide what to do.
              err_report = self.StructureToOptimizer (system, imageTrajectory)
              if err_report is not None:
                 logFile.Paragraph ("\nCPR cycle terminated!\n {:s}".format (err_report))
                 break
              # . Initialize again to update if the path is not yet fully refined.
              if not self.fullyRefinedPath:
                  self.InitializeAgain (imageTrajectory, system)
                  self.PathSummary ("other")
                  self.CountInSegment ()  # . Count inSegment images for self.Continue().
                  self.WriteOutRestartFile ()
                  self.CPRrun += 1
              
            if (self.finalUnrefinableRefinement and self.error is None):  # . Final unrefinable refinement
              self.finalUnrefRef = True
    
              # . Find segment containing highest point and decide what to do.
              err_report = self.StructureToOptimizer (system, imageTrajectory)
              # . Initialize again to update.
              self.InitializeAgain (imageTrajectory, system)
              self.PathSummary ("other")
              self.WriteOutRestartFile () 
              if err_report is not None:
                 logFile.Paragraph ("\nCPR cycle terminated!\n {:s}".format (err_report))
              else: logFile.Paragraph ("CPR finished successfully after trying to refine remaining points marked as unrefinable.")
            if (self.exitMessage is not None):
                logFile.Paragraph ("\nCPR terminated\n {:s}".format (self.exitMessage))
                self.Footer ( log = log )
        except Exception as error:
            logFile.Paragraph ( "\n*** Error : {:s}. ***".format ( error.args[0] ) )
        self.Finalize ( )
        return self.state

    def DetermineStepSize (self, system):
        """Estimate step size of the segment scanning steps."""
        # . Determine step size of initial scan.
        reactant = Array.WithExtent ( self.numberOfVariables ) ; reactant.Set ( 0.0 ) 
        product  = Array.WithExtent ( self.numberOfVariables ) ; product.Set  ( 0.0 )
        self.variables[ 0                       ].CopyTo ( reactant )  # . Start structure.
        self.variables[ self.numberOfImages - 1 ].CopyTo ( product  )  # . End structure.
        product.Add ( reactant, scale = -1.0 ) 
        displacement = product.Norm2 () / float (3 * self.numberOfVariables)
        # . Check for valid path.
        Rsmall = 1.0e-10
        Rbig = 1.0e+20
        if ( displacement < Rsmall and self.numberOfImages < 3 ):
           logFile.Paragraph ("\nIdentical reactant and product, but no intermediate path-point between them.")
           logFile.Paragraph ("Should have at least a different initial and final state or one intermediate point.")
           raise ValueError ("Need at least two different structures!")
        largestStepSize = displacement / 6.0
        largestStepSize *= math.sqrt (3 * self.numberOfVariables) * math.sqrt (3 * self.numberOfVariables)
        if (self.nFraction > 0 and displacement > Rsmall):
           largestStepSize = ((displacement * math.sqrt (3 * self.numberOfVariables) * math.sqrt (3 * self.numberOfVariables)) / float (self.nFraction + 1)) 
        if (self.upperStepLimit is not None):
           largestStepSize = self.upperStepLimit
        else:
           largestStepSize /= math.sqrt(3 * self.numberOfVariables)     
        largestStepSize *= math.sqrt(3 * self.numberOfVariables)  
        if (largestStepSize < Rsmall):          
           logFile.Paragraph ("\nSize of the step for the scan is too small, largestStepSize = {:f}".format (largestStepSize))
           if (displacement < Rsmall):
              logFile.Paragraph ("\nIdentical reactant and product structures.")
              logFile.Paragraph ("Give largest allowed scanning step explicitly.")
              logFile.Paragraph ("Adjust 'upperStepLimit' option to a value greater than {:f}".format (Rsmall))
           if (self.stepsPerSegment >= 2):
              largestStepSize = Rbig
              logFile.Paragraph ("\nStep size is ignored and set to {:f}".format (largestStepSize))
              logFile.Paragraph ("for further step size calculations.")
           else:
              logFile.Paragraph ("Step size cannot be near zero! Failed!")
              raise ValueError ("Adjust the step size!") 
        self.upperStepLimit = largestStepSize
        if (self.lowestStepLimit is not None):
           self.lowestStepLimit = self.lowestStepLimit
        else:
           self.lowestStepLimit = 0.0 
        self.lowestStepLimit *= (3 * self.numberOfVariables)

    def Finalize (self):
        """Clear up after working on a function."""
        # . The final state should correspond to the best path found.
        # . Finalize the optimizer.
        self.optimizerState.Finalize ()
        self.state["Exit Message"              ] = self.exitMessage    
        if self.error is not None: self.state["Error"] = self.error
        # . Delete all objective function attributes.
        keys = self.__class__.objectiveFunctionAttributes.keys ()
        for key in keys: delattr (self, key) 

    def Footer ( self, log = logFile ):
        """Footer."""
        if LogFileActive ( log ):
            userInformation = """                      
     
      Thank you for using PyCPR
     
      Florian J. Gisdon, Martin Culka, G. Matthias Ullmann (2016): PyCPR - A Python-based
      Implementation of the Conjugate Peak Refinement (CPR) Algorithm for Finding
      Transition State Structures. J. Mol. Model., 22: 242, 2016
      DOI:  10.1007/s00894-016-3116-8
                          
                              """
            log.Text ( userInformation )

    def Header ( self, log = logFile ):
        """Header."""
        if LogFileActive ( log ):
            userInformation = """
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
            
          Python-based Conjugate Peak Refinement (PyCPR)
                              
                      Program Version 1.0
                                       
      If you find PyCPR useful and use it 
      for publication, please cite: 
      
      Florian J. Gisdon, Martin Culka, G. Matthias Ullmann (2016): PyCPR - A Python-based
      Implementation of the Conjugate Peak Refinement (CPR) Algorithm for Finding
      Transition State Structures. J. Mol. Model., 22: 242, 2016
      DOI:  10.1007/s00894-016-3116-8
   
              
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
                              """
            log.Text ( userInformation )

    def HighestPoint (self):
        """Finding highest energy structure over all segments."""
        highestEnergy = None
        segment = None
        for i in range (0, self.numberOfImages - 1):
           if not ((self.segmentsHighest[i][1] == "saddle") or                        \
                   (self.segmentsHighest[i][1] == "stationary") or                    \
                   (self.segmentsHighest[i][1] == "unrefStart") or                    \
                   (self.segmentsHighest[i][1] == "unrefStop") or                     \
                   (self.segmentsHighest[i][1] == "initial") or                     \
                   (self.segmentsHighest[i][1] == "final") or                     \
                   ((i == 0) and (self.segmentsHighest[i][1] == "flankStart")) or     \
                   ((i == self.numberOfImages - 2) and (self.segmentsHighest[i][1] == "flankStop"))):
               if (highestEnergy is None):
                   highestEnergy = self.segmentsHighest[i][2]
                   segment = i
               elif (highestEnergy is not None):
                  if (self.segmentsHighest[i][2] > highestEnergy):
                     highestEnergy = self.segmentsHighest[i][2]
                     segment = i
        return (segment)
    
    def Initialize (self, objectiveFunction):
        """Set up before iterating."""
        # . Initialize the state.
        self.state = {}
        # . Set all objective function attributes.
        keyvalues = self.__class__.objectiveFunctionAttributes.items ( )
        for (key, value) in keyvalues: setattr (self, key, value)
        # . Check the objective function.
        if not self.IsGoodObjectiveFunction (objectiveFunction): raise TypeError ("Invalid objective function.")
        self.objectiveFunction = objectiveFunction
        # . Number of variables and images.
        self.numberOfImages = self.objectiveFunction.numberOfImages
        self.numberOfVariables = self.objectiveFunction.numberOfVariables
        # . Allocate space.
        self.energy = [None for iImage in range (self.numberOfImages)]
        self.rmsGradient = [None for iImage in range (self.numberOfImages)]
        self.isSaddle = [False for iImage in range (self.numberOfImages)]
        self.isStationary = [False for iImage in range (self.numberOfImages)]
        self.isUnrefinable = [False for iImage in range (self.numberOfImages)]
        self.ds0 = [ Array.WithExtent (self.numberOfVariables) for iImage in range (self.numberOfImages - 1)]
        self.segmentsHighest = [ (Array.WithExtent (self.numberOfVariables), None, None) for iImage in range (self.numberOfImages - 1) ]
        self.segmentsLowest = [ (Array.WithExtent (self.numberOfVariables), None, None) for iImage in range (self.numberOfImages - 1) ]
        self.variables = [ Array.WithExtent (self.numberOfVariables) for iImage in range (self.numberOfImages) ]
        self.listOfAddedStructures = []
        for iImage in range (self.numberOfImages - 1):
            self.segmentsHighest[iImage][0].Set (0.0)
            self.segmentsLowest[iImage][0].Set (0.0)
            self.ds0[iImage].Set (0.0)
        for iImage in range (self.numberOfImages):
            self.variables[iImage].Set (0.0)
        # . Read trajectory into self.variables.
#        from pScientific.Arrays import ArrayPrint, ArrayPrint2D
#        ArrayPrint2D ( self.objectiveFunction.system.coordinates3, itemFormat = "{:.3f}", title = "Coordinates3 Before" )
        for iImage in range (0, self.numberOfImages):
            self.objectiveFunction.ActivateImage (iImage)
            self.objectiveFunction.VariablesGet (self.variables[iImage])
#            ArrayPrint2D ( self.objectiveFunction.system.coordinates3, itemFormat = "{:.3f}", title = "Coordinates3 {:d}".format ( iImage ) )
#            ArrayPrint ( self.variables[iImage], itemFormat = "{:.3f}", title = "Variables {:d}".format ( iImage ) )
        # . Set up the optimizer.
        if (self.maximalMinimizationSteps == None): # . Set the default value for maximalMinimizationSteps
            self.maximalMinimizationSteps = (self.numberOfVariables - 1)
        if (self.numberOfSuccessiveLowGradients == None):   # . Set the default value for numberOfSuccessiveLowGradients
            self.numberOfSuccessiveLowGradients = int (math.sqrt (self.numberOfVariables / 3))
        if (self.maximalMinimizationSteps <= self.numberOfSuccessiveLowGradients):
            logFile.Paragraph ("The number of maximalMinimizationSteps is smaller then the numberOfSuccessiveLowGradients.")
            logFile.Paragraph ("The algorithm quits, because it is impossible to find saddle points. The number of maximalMinimizationSteps has to be greater than the numberOfSuccessiveLowGradients.")
            self.error = "ERROR!\nNumber of maximalMinimizationSteps and numberOfSuccessiveLowGradients don't allow to find saddle points.\nCPR terminates."
        elif (self.maximalMinimizationSteps < 2*self.numberOfSuccessiveLowGradients):
            logFile.Paragraph ("\nWARNING!")
            logFile.Paragraph ("The number of maximalMinimizationSteps is smaller then two times the numberOfSuccessiveLowGradients.")
            logFile.Paragraph ("This could complicate the procedure to find saddle points. It is advisable to increase the number of maximalMinimizationSteps.\n")
        self.optimizer = CPRSaddlePointRefinement.WithOptions ( maximalMinimizationSteps       = self.maximalMinimizationSteps       ,
                                                                numberOfSuccessiveLowGradients = self.numberOfSuccessiveLowGradients ,
                                                                rmsGradientTolerance           = self.rmsGradientTolerance           ,
                                                                rmsdMaximumTolerance           = self.rmsdMaximumTolerance           ,
                                                                finiteDifferenceStep           = self.finiteDifferenceStep           ,
                                                                finiteDifference               = self.finiteDifference               ,
                                                                initialStep                    = self.lineSearchInitialStep          ,
                                                                orthogonal                     = self.orthogonal                     ,
                                                                betaType                       = self.betaType                       ,
                                                                breakIfTauReached              = self.breakIfTauReached              )
        self.optimizerState = self.optimizer.StateFromObjectiveFunction (self.objectiveFunction)
        #self.optimizer.Initialize (self.optimizerState)  # . Unnecessary to initialize here, because it is initialized within each run of the optimizer
        # . Check options
        if (not self.breakIfTauReached and self.increaseTau):
            logFile.Paragraph ("\nWARNING!")
            logFile.Paragraph ("The option to break if the threshold for the value tau is reached was set to False and the option to increase tau was set to True.") 
            logFile.Paragraph ("Since 'increaseTau = True' has no influence in that case it was set to False.") 
            self.increaseTau = False
        # . Check if a restart file exists and rename if necessary
        if(os.path.isfile("{:s}".format(self.outputRestartFile))):
            logFile.Paragraph ("Restart file already exists, moving to {:s}.old-renamed{:s}".format (self.outputRestartFile, strftime("%y%m%d%H%M%S")))
            os.rename("{:s}".format (self.outputRestartFile), "{:s}.old-renamed{:s}".format (self.outputRestartFile,strftime("%y%m%d%H%M%S")))
        if self.error != None:
            logFile.Paragraph ( self.error )
            exit()

    def InitializeAgain (self, imageTrajectory, system):
        """Update before iterating."""
        trajectory = ExportTrajectory ( imageTrajectory.path, system, append = True )
        self.objectiveFunction.InitializeImages (trajectory)
        self.numberOfImages = len (self.variables)

    def IsGoodObjectiveFunction (self, objectiveFunction):
        """Return a flag indicating if the objective function is valid."""
        return isinstance (objectiveFunction, CPRObjectiveFunction) and (objectiveFunction.numberOfImages >= 2)

    def Optimization (self, system, imageTrajectory, segment):
        """Perform a thorough line maximization and minimize in conjugate directions."""
        # . Check if scratch directory is there and delete content for the next calculation.
        if self.scratchDirectory:
            logFile.Paragraph ("\nDeleted files:")
            for item in os.listdir(self.scratchDirectory):
                  scratchFolderContent = os.path.join (self.scratchDirectory, "{:s}".format (item))
                  os.remove (scratchFolderContent)
                  logFile.Paragraph ( "{:s}".format (scratchFolderContent) )
        # . Perform a thorough line maximization.
        # . Find the normalized direction of the maximization.
        tangent = Array.WithExtent (self.numberOfVariables); tangent.Set (0.0)
        tmp = Array.WithExtent (self.numberOfVariables); tmp.Set (0.0)
        # . Decide if the refinement was good for existing path-point refinement, for new point "elif".
        if (self.segmentsHighest[segment][1] != "inSegment"):
           if (self.segmentsHighest[segment][1] == "flankStart") or (self.segmentsHighest[segment][1] == "unrefStart"): 
               i = segment
               if segment == 0:
                   return ("\nImpossible to refine the initial path point!")
           elif (self.segmentsHighest[segment][1] == "flankStop") or (self.segmentsHighest[segment][1] == "unrefStop"): 
               i = segment + 1
               if segment == (self.numberOfImages - 2):
                   return ("\nImpossible to refine the final path point!")
           self.objectiveFunction.VariablesPut (self.variables[i])  # . Load in the image to optimizer.
           self.optimizerState = self.optimizer.StateFromObjectiveFunction (self.objectiveFunction)
           # . Determine the tangent.
           self.variables[i - 1].CopyTo (tmp)
           tmp.Add ( self.variables[i], scale = -1.0 )
           tmp.Normalize () 
           self.variables[i + 1].CopyTo (tangent)
           tangent.Add ( self.variables[i], scale = -1.0 )
           tangent.Normalize ()
           tangent.Add ( tmp, scale = -1.0 )
           # . Search for a nearby maximum.
           tangent.CopyTo (self.optimizerState.s0)  
           self.optimizer.onlyMax = True
           self.optimizer.tauTol = self.tauTolRefine
           report = self.optimizer.Iterate(self.optimizerState)
           if report["Converged"]:
               logFile.Paragraph ("\nThere is a nearby maximum, refine the path point {:d} along the tangent!".format (i))
               self.optimizer.onlyMax = False
               report = self.optimizer.Iterate(self.optimizerState)
               self.optimizerState.x.CopyTo (self.variables[i])
               self.rmsGradient[i] = self.optimizerState.rmsGradient
               # . Update the refined image.  
               if self.optimizerState.saddle: 
                   self.numberOfSaddle += 1
                   self.isSaddle[i] = self.optimizerState.saddle
               elif self.optimizerState.stationary:
                   self.numberOfStationary += 1
                   self.isStationary[i] = self.optimizerState.stationary
               temporary = Clone (system.coordinates3)
               self.objectiveFunction.VariablesPut (self.variables[i])
               imageTrajectory.WriteFrame (system.coordinates3, frame=i)
               self.energy[i] = self.optimizerState.f
               self.segmentsHighest[segment][2] = self.energy[i]
# . MJF.
               system.coordinates3 = temporary  # Restore temporary system coordinates
               if system.freeAtoms is None: self.objectiveFunction.iCoordinates3 = system.coordinates3.iterator
               else:                        self.objectiveFunction.iCoordinates3 = system.coordinates3.RowIterator ( selection = self.freeAtoms )
#
               #system.coordinates3 = temporary
               self.InitializeAgain (imageTrajectory, system)
               for k in  (i - 1, i):  # . Update segment information.
                    self.SegmentInterpolation (system, imageTrajectory, k)
           else:
             if self.finalUnrefRef == True:
               logFile.Paragraph ("\nThere is no nearby maximum along s0 for unrefinable point {:d}!".format (i))
               return
             else:
               logFile.Paragraph ("\nThere is no nearby maximum along s0, deleting the path point {:d}!".format (i))
               temporary = Clone (system.coordinates3)
               del self.energy[i]
               del self.rmsGradient[i]
               del self.isSaddle[i]
               del self.isStationary[i]
               del self.isUnrefinable[i]
               del self.variables[i]
               # . Delete written trajectory frames before writing them new.
               delFrame = os.path.join (imageTrajectory.path, "frame{:d}".format (self.numberOfImages - 1) + PickleFileExtension) 
               os.remove (delFrame)
               # logFile.Paragraph ("\nImage {:s} deleted from path.".format (delFrame))
               # . Write trajectory new.
               for j in range (0, self.numberOfImages - 1):
                    self.objectiveFunction.VariablesPut (self.variables[j])
                    imageTrajectory.WriteFrame (system.coordinates3, frame=j)
               del self.segmentsHighest[i]
               del self.segmentsLowest[i]
               del self.ds0[i]
# . MJF.
               system.coordinates3 = temporary  # Restore temporary system coordinates
               if system.freeAtoms is None: self.objectiveFunction.iCoordinates3 = system.coordinates3.iterator
               else:                        self.objectiveFunction.iCoordinates3 = system.coordinates3.RowIterator ( selection = self.freeAtoms )
#
#               system.coordinates3 = temporary
               self.InitializeAgain (imageTrajectory, system)
               self.SegmentInterpolation (system, imageTrajectory, i - 1)
               return
        # . Decide if the refinement was good for refinement of a new point.
        elif (self.segmentsHighest[segment][1] == "inSegment"):
           # . Copy the image to be optimized and the initial search direction to the optimizer state.
           self.objectiveFunction.VariablesPut (self.segmentsHighest[segment][0])
           self.optimizerState = self.optimizer.StateFromObjectiveFunction (self.objectiveFunction)
           self.ds0[segment].CopyTo (self.optimizerState.s0)

           ############################################################
           # . Test if structure was already added
           structureFound = None
           for structure in range(len(self.listOfAddedStructures)):
               # Calculate RMSD of structure
               rmsdAdded = -1
               #logFile.Paragraph structure
               self.listOfAddedStructures[structure][0].CopyTo (tmp)
               tmp.Add ( self.segmentsHighest[segment][0], scale = -1.0 )
               rmsdAdded = tmp.Norm2()
               #logFile.Paragraph "RMSD of added structure in comparison to all structures in the list: {:f}".format(rmsdAdded)
               if (0.0 <= rmsdAdded < self.rmsdLoopPreventionX0):
                   # Calculate RMSD of ds0
                   rmsdAdded = -1
                   self.listOfAddedStructures[structure][1].CopyTo (tmp)
                   tmp.Add ( self.ds0[segment], scale = -1.0 )
                   rmsdAdded = tmp.Norm2()
                   #logFile.Paragraph "RMSD of added ds0 in comparison to ds0 of hit: {:f}".format(rmsdAdded)
                   if (0.0 <= rmsdAdded < self.rmsdLoopPreventionS0):
                       structureFound = structure
                       break
           
           if structureFound != None: 
           #    change stored information and determine procedure
                procedure = self.listOfAddedStructures[structure][2]
                self.listOfAddedStructures[structure][2] += 1
           elif structureFound == None:
                # . Copy x0 and s0 to the list of structures and vectors for added structures
                tmpAddStruc = [ Array.WithExtent (self.numberOfVariables), Array.WithExtent (self.numberOfVariables), 0 ] # . x0,s0,state of optimization method
                for position in range (2):
                    tmpAddStruc[position].Set (0.0)
                self.segmentsHighest[segment][0].CopyTo (tmpAddStruc[0])
                self.ds0[segment].CopyTo (tmpAddStruc[1])
                tmpAddStruc[2] += 1
                # . Add temp to list
                self.listOfAddedStructures.append(tmpAddStruc)
                procedure = 0
           #logFile.Paragraph "Procedure {:d}".format(procedure)
           ############################################################
           
           unrefinable = False
           # . Do the maximization and subsequent conjugate minimizations.
           self.optimizer.onlyMax = False
           if (procedure == 0):
               self.optimizer.tauTol = self.tauTolAdd
           elif (procedure > 0):
               if(self.increaseTau):
                   if (self.interpolationIncreasement):
                       if (procedure == 1):
                           logFile.Paragraph ("\nFutile loop detected, interpolation number was increased by {:d}".format(self.interpolationIncreasement))
                           self.SegmentInterpolation (system, imageTrajectory, segment, modification = self.interpolationIncreasement)  # . Increased interpolation number
                           if ((self.segmentsHighest[segment][1] == "saddle") or                        \
                                   (self.segmentsHighest[segment][1] == "stationary") or                    \
                                   (self.segmentsHighest[segment][1] == "unrefStart") or                    \
                                   (self.segmentsHighest[segment][1] == "unrefStop") or                     \
                                   (self.segmentsHighest[segment][1] == "initial") or                     \
                                   (self.segmentsHighest[segment][1] == "final") or                     \
                                   ((segment == 0) and (self.segmentsHighest[segment][1] == "flankStart")) or     \
                                   ((segment == self.numberOfImages - 2) and (self.segmentsHighest[segment][1] == "flankStop"))):
                                    logFile.Paragraph ("Nothing to refine within the segment.")
                                    return
                           elif ((self.segmentsHighest[segment][1] == "flankStart") or                        \
                                   (self.segmentsHighest[segment][1] == "flankStop")):
                                    logFile.Paragraph ( "NOTE: Recurison occurred! Optimization called from within itself!" )
                                    self.Optimization (system, imageTrajectory, segment)
                                    logFile.Paragraph ( "WARNING! The function proceeds after the recursive call!" )
                                    return                     # . Possible solution for getting flanking points after interpolation increasement: Recursion
                           else:
                                    self.objectiveFunction.VariablesPut (self.segmentsHighest[segment][0])
                                    self.optimizerState = self.optimizer.StateFromObjectiveFunction (self.objectiveFunction)   
                                    self.ds0[segment].CopyTo (self.optimizerState.s0)
                                    self.optimizer.tauTol = self.tauTolAdd
                                    self.changeInterpolations = segment + 1 
                       elif (procedure == 2):
                           self.optimizer.tauTol = 2 * self.tauTolAdd
                           self.increasedTau = segment + 1
                           logFile.Paragraph ("\nFutile loop detected, tau was increased")
                       elif (procedure == 3):
                           self.optimizer.orthogonal = True  # . Orthogonal directions.
                           self.optimizer.tauTol = 2 * self.tauTolAdd
                           logFile.Paragraph ("\nFutile loop detected, switching to orthogonal directions for this run")
                       elif (procedure > 3):
                           unrefinable = True
                   elif (not self.interpolationIncreasement):  
                       if (procedure == 1):
                           self.optimizer.tauTol = 2 * self.tauTolAdd
                           self.increasedTau = segment + 1
                           logFile.Paragraph ("\nFutile loop detected, tau was increased")
                       elif (procedure == 2):
                           self.optimizer.orthogonal = True  # . Orthogonal directions.
                           self.optimizer.tauTol = 2 * self.tauTolAdd
                           logFile.Paragraph ("\nFutile loop detected, switching to orthogonal directions for this run")
                       elif (procedure > 2):
                           unrefinable = True
               elif (not self.increaseTau): 
                   if (self.interpolationIncreasement):
                       if (procedure == 1):
                           logFile.Paragraph ("\nFutile loop detected, interpolation number was increased by {:d}".format(self.interpolationIncreasement))
                           self.SegmentInterpolation (system, imageTrajectory, segment, modification = self.interpolationIncreasement)  # . Increased interpolation number
                           if ((self.segmentsHighest[segment][1] == "saddle") or                        \
                                   (self.segmentsHighest[segment][1] == "stationary") or                    \
                                   (self.segmentsHighest[segment][1] == "unrefStart") or                    \
                                   (self.segmentsHighest[segment][1] == "unrefStop") or                     \
                                   (self.segmentsHighest[segment][1] == "initial") or                     \
                                   (self.segmentsHighest[segment][1] == "final") or                     \
                                   ((segment == 0) and (self.segmentsHighest[segment][1] == "flankStart")) or     \
                                   ((segment == self.numberOfImages - 2) and (self.segmentsHighest[segment][1] == "flankStop"))):
                                    logFile.Paragraph ("Nothing to refine within the segment.")
                                    return
                           else:
                                    self.objectiveFunction.VariablesPut (self.segmentsHighest[segment][0])
                                    self.optimizerState = self.optimizer.StateFromObjectiveFunction (self.objectiveFunction)   
                                    self.ds0[segment].CopyTo (self.optimizerState.s0)
                                    self.optimizer.tauTol = self.tauTolAdd
                                    self.changeInterpolations = segment + 1 
                       elif (procedure == 2):
                           self.optimizer.orthogonal = True  # . Orthogonal directions.
                           logFile.Paragraph ("\nFutile loop detected, switching to orthogonal directions for this run")
                       elif (procedure > 2):
                           unrefinable = True               
                   elif (not self.interpolationIncreasement):
                       if (procedure == 1):
                           self.optimizer.orthogonal = True  # . Orthogonal directions.
                           logFile.Paragraph ("\nFutile loop detected, switching to orthogonal directions for this run")
                       elif (procedure > 1):
                           unrefinable = True 
           
           ############################################################  
           # . Start the optimizer (maximization and minimization)
           report = self.optimizer.Iterate(self.optimizerState)
           self.optimizer.orthogonal = self.orthogonal  # . Restore original setting for orthogonal directions.
           # . Copy back information from the optimizer.
           self.optimizerState.x.CopyTo (self.segmentsHighest[segment][0])
           self.segmentsHighest[segment][2] = self.optimizerState.f
           if self.optimizerState.saddle: self.numberOfSaddle += 1
           elif self.optimizerState.stationary: self.numberOfStationary += 1
           # . Save the current system coordinates.           
           temporary = Clone (system.coordinates3)
           # . Save the images after the new one which will be inserted
           for j in range (self.numberOfImages - 1, segment, -1):  # . Reverse order.
              self.objectiveFunction.VariablesPut (self.variables[j])
              imageTrajectory.WriteFrame (system.coordinates3, frame=j + 1) 
           # . New image
           self.objectiveFunction.VariablesPut (self.segmentsHighest[segment][0])
           imageTrajectory.WriteFrame (system.coordinates3, frame=segment + 1)
           # . Insert generated image and additional information
           self.variables.insert (segment + 1, self.segmentsHighest[segment][0])
           self.energy.insert (segment + 1, self.segmentsHighest[segment][2])
           self.rmsGradient.insert (segment + 1, self.optimizerState.rmsGradient)
           self.isSaddle.insert (segment + 1, self.optimizerState.saddle)
           self.isStationary.insert (segment + 1, self.optimizerState.stationary)
           if unrefinable:
               self.isUnrefinable.insert (segment + 1, True)
               self.numberOfUnrefinables += 1
           else:
               self.isUnrefinable.insert (segment + 1, False)
           tmpHighSeg = [ Array.WithExtent (self.numberOfVariables), None, None ]; tmpHighSeg[0].Set (0.0)
           tmpLowSeg = [ Array.WithExtent (self.numberOfVariables), None, None ]; tmpLowSeg[0].Set (0.0)
           self.segmentsHighest.insert (segment + 1, tmpHighSeg)
           self.segmentsLowest.insert (segment + 1, tmpLowSeg)
           tmpds0 = Array.WithExtent (self.numberOfVariables) ; tmpds0.Set (0.0)
           self.ds0.insert (segment + 1, tmpds0)
           self.InitializeAgain (imageTrajectory, system)
           for k in  (segment, segment + 1):
               self.SegmentInterpolation (system, imageTrajectory, k)
           # . Restore temporary system coordinates.
#           system.coordinates3 = temporary 
# . MJF.
           system.coordinates3 = temporary  # Restore temporary system coordinates
           if system.freeAtoms is None: self.objectiveFunction.iCoordinates3 = system.coordinates3.iterator
           else:                        self.objectiveFunction.iCoordinates3 = system.coordinates3.RowIterator ( selection = self.freeAtoms )
#
           return

    def PathSummary ( self, state, log = logFile ):
        """Output initial path summary."""
        if LogFileActive (log):
            # . Output.
            table = log.GetTable ( columns = [ 8, 22, 30, 12, 20 ] )
            table.Start ( )
            if   state == "initial"                              : title = "Path Summary before CPRrun {:d}".format ( self.CPRrun )
            elif state == "other" and ( not self.finalUnrefRef ) : title = "Path Summary after CPRrun {:d}".format  ( self.CPRrun )
            elif state == "other" and       self.finalUnrefRef   : title = "Path Summary after final unrefinable refinement"
            table.Title   ( title )
            table.Heading ( "Image"        )
            table.Heading ( "Energy"       )
            table.Heading ( "Rel. Energy"  )
            table.Heading ( "Type"         )
            table.Heading ( "RMS Gradient" )
            for i in range ( self.numberOfImages ):
                if   i == 0                       : imageType = "initial"
                elif i == self.numberOfImages - 1 : imageType = "final"
                elif self.isSaddle      [i]       : imageType = "saddle"
                elif self.isStationary  [i]       : imageType = "stationary"
                elif self.isUnrefinable [i]       : imageType = "unrefinable"
                else                              : imageType = "normal"
                table.Entry ( "{:d}".format     ( i                               ) )  
                table.Entry ( "{:.6f}".format   ( self.energy[i]                  ) )
                table.Entry ( "{:.2f}".format   ( self.energy[i] - self.energy[0] ) )
                table.Entry ( "{:10s}".format   ( imageType                       ) )
                if self.rmsGradient[i] is None: table.Entry ( "n/a" )
                else:                           table.Entry ( "{:12.6f}".format ( self.rmsGradient[i] ) )
            table.Stop ( )
            table = log.GetTable ( columns = [ 5, 20, 22, 12, 22, 12, 18 ] )
            table.Start ( )
            if   state == "initial"                              : title = "Segment overview before CPRrun {:d}".format ( self.CPRrun )
            elif state == "other" and ( not self.finalUnrefRef ) : title = "Segment overview after CPRrun {:d}".format  ( self.CPRrun )
            elif state == "other" and       self.finalUnrefRef   : title = "Path Summary after final unrefinable refinement"
            table.Title   ( title            )
            table.Heading ( "Segment"        )
            table.Heading ( "E."             )
            table.Heading ( "Rel. E."        )
            table.Heading ( "Highest"        )
            table.Heading ( "Rel. E."        )
            table.Heading ( "Lowest"         )
            table.Heading ( "Segment Length" )
            for i in range ( self.numberOfImages - 1 ):
                table.Entry ("{:d}".format   ( i                                           ) )
                table.Entry ("{:.5f}".format ( self.segmentsHighest[i][2]                  ) )
                table.Entry ("{:.2f}".format ( self.segmentsHighest[i][2] - self.energy[0] ) )
                table.Entry ("{:10s}".format ( self.segmentsHighest[i][1]                  ) )
                table.Entry ("{:.2f}".format ( self.segmentsLowest[i][2] - self.energy[0]  ) )
                table.Entry ("{:10s}".format ( self.segmentsLowest[i][1]                   ) )
                table.Entry ("{:.5f}".format ( self.ds0[i].Norm2 ( )                       ) )
            table.Stop ( )

    def ReadInRestartFile (self, system , imageTrajectory):
         """Read in the restart file."""
         images = True
         with  open(self.inputRestartFile, 'r') as restartFile:
             for line in restartFile:
                 if line.strip():  # . Ignore empty lines.
                     columns = string.split (line)
                     # . Reading in the parts and skip each header.
                     if (columns[0] == "Image"): i = 0; continue
                     if (columns[0] == "Segment"): i = 0; images = False; continue
                     if images:
                         # . 2nd column: image energy; 3rd column: image type; 4th column: image gradient (if available).
                         self.energy[i] = float (columns[1])
                         if columns[2] == "saddle": 
                             self.isSaddle[i] = True
                             self.numberOfSaddle += 1
                         elif columns[2] == "stationary": 
                             self.isStationary[i] = True
                             self.numberOfStationary += 1
                         elif columns[2] == "unrefinable": 
                             self.isUnrefinable[i] = True
                             self.numberOfUnrefinables += 1
                         if len (columns) == 4 and (isinstance(columns[3], int) or isinstance(columns[3], float)):
                             self.rmsGradient[i] = float (columns[3])
                     else:
                         # . 2nd column: segment energy; 3rd column: segment type.
                         tmpHighSeg = [ Array.WithExtent (self.numberOfVariables), None, None ]; tmpHighSeg[0].Set (0.0)
                         tmpHighSeg[2] = float (columns[1])
                         tmpHighSeg[1] = columns[2]
                         self.segmentsHighest[i] = tmpHighSeg
                         tmpLowSeg = [ Array.WithExtent (self.numberOfVariables), None, None ]; tmpLowSeg[0].Set (0.0)
                         tmpLowSeg[2] = float (columns[3])
                         tmpLowSeg[1] = columns[4]
                         self.segmentsLowest[i] = tmpLowSeg
                         if self.segmentsHighest[i][1] == "inSegment":  # . For type 'inSegment', scanning of that segment required.
                             self.SegmentInterpolation (system , imageTrajectory, i) 
                     i += 1

    def SegmentInterpolation (self, system, imageTrajectory, i, modification = False):
        """Scanning for local E max. along trajectory."""
        # . Check if scratch directory is there and delete content for the next calculation.
        if self.scratchDirectory:
            logFile.Paragraph ("\nDeleted files:")
            for item in os.listdir(self.scratchDirectory):
                  scratchFolderContent = os.path.join (self.scratchDirectory, "{:s}".format (item))
                  os.remove (scratchFolderContent)
                  logFile.Paragraph ( "{:s}".format (scratchFolderContent) )
        # . Check for fixed atoms. Superimposing all coordinates3 with reference.
        if system.freeAtoms is None:
            reference = imageTrajectory.ReadFrame (frame=0)
            system.coordinates3.Superimpose (reference)
        # . Calculate vector between two neighboring path-points
        logFile.Paragraph ("\nScanning for local energy maximum along trajectory")
        logFile.Paragraph ("Scanning segment {:d}".format (i))
        if (self.numberOfImages >= 2):
           # . Calculate the vector between two path-points iImage, iImage+1 ('tangent').
           # . This is the tangent for the max. energy structure along that path.
           start = Array.WithExtent (self.numberOfVariables);  start.Set (0.0) 
           self.variables[i    ].CopyTo ( start       )  # . Start structure.
           self.variables[i + 1].CopyTo ( self.ds0[i] )  # . End structure (just temporary).
           self.ds0[i].Add ( start, scale = -1.0 )  # . ds0[i] = Vector between End and Start of segment [i].
           # . Determine step size for each segment
           segmentLength = self.ds0[i].Norm2 ()
           try:
               1.0 / segmentLength
           except ZeroDivisionError:
               logFile.Paragraph ( """WARNING: Structures with number {0:d} and {1:d} are the same, 
         no action performed for this segment. 
         Structure {0:d} is set as highest energy structure 
         of this segment, structure {1:d} as lowest.""".format(i, i+1) )
               self.variables[i    ].CopyTo (self.segmentsHighest[i][0])  # . Start structure.
               self.objectiveFunction.VariablesPut (self.segmentsHighest[i][0])
               self.energy[i] = self.energy[i + 1] = system.Energy ( log = None )
               self.variables[i + 1].CopyTo (self.segmentsLowest[i][0])  # . Start structure.
               tmpHighSeg = [ Array.WithExtent (self.numberOfVariables), None, None ];  tmpHighSeg[0].Set (0.0)
               # . Highest energy image is flanking image start.
               if self.isSaddle[i]:
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}, and it is a saddle point".format (i))
                 tmpHighSeg[1] = "saddle"
               elif self.isStationary[i]:
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}, and it is a stationary point".format (i))
                 tmpHighSeg[1] = "stationary"
               elif self.isUnrefinable[i]:
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}, but it is an unrefinable point".format (i))
                 tmpHighSeg[1] = "unrefStart"
               else:   
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}".format (i))
                 logFile.Paragraph ( tmpHighSeg[1] )
                 tmpHighSeg[1] = "flankStart"
               tmpHighSeg[2] = self.energy[i]
               self.segmentsHighest[i] = tmpHighSeg
               tmpLowSeg = [ Array.WithExtent (self.numberOfVariables), None, None ];  tmpLowSeg[0].Set (0.0)
               logFile.Paragraph ("E (kJ/mol) = {:f}".format (self.energy[i]))
               # . Lowest energy image is flanking image stop.
               if self.isSaddle[i + 1]:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}, and it is a saddle point".format (i))
                  tmpLowSeg[1] = "saddle"
               elif self.isStationary[i + 1]:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}, and it is a stationary point".format (i))
                  tmpLowSeg[1] = "stationary"
               elif self.isUnrefinable[i + 1]:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}, and it is an unrefinable point".format (i))
                  tmpLowSeg[1] = "unrefStop"
               else:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}".format (i))
                  tmpLowSeg[1] = "flankStop"
               tmpLowSeg[2] = self.energy[i + 1]
               logFile.Paragraph ("E (kJ/mol) = {:f}".format (self.energy[i + 1]))
               self.segmentsLowest[i] = tmpLowSeg
               return
           
           if (self.stepsPerSegment >= 2):
              stepSize = min (float (self.upperStepLimit), float (segmentLength) / float (self.stepsPerSegment))
           else:
              stepSize = self.upperStepLimit
           if (stepSize < self.lowestStepLimit):
              stepSize = self.lowestStepLimit
           interpolations = max (int (float (segmentLength) / float (stepSize)), 1)
           if (modification):
               interpolations += self.interpolationIncreasement
           stepSize = float (segmentLength) / float (interpolations)
           logFile.Paragraph ("Interpolation number,            interpolations    =   {:d}".format (interpolations))
           logFile.Paragraph ("Size of the scanning step,       stepSize          =   {:f}".format (stepSize))
           # . Checking function for energy maximum along that path
           activeImage = Array.WithExtent (self.numberOfVariables);  activeImage.Set (0.0) 
           self.variables[i    ].CopyTo (activeImage)  # . Start structure.
           EmaxImage = Array.WithExtent (self.numberOfVariables);  EmaxImage.Set (0.0) 
           EminImage = Array.WithExtent (self.numberOfVariables);  EminImage.Set (0.0)
           highestEnergy = None
           lowestEnergy = None
           tempEnergy = None
           # . Save the current system coordinates.           
           temporary = Clone (system.coordinates3)
           # Checking generated images.
           for scan in range (1, interpolations):
              start.CopyTo (activeImage)  # . Start structure.
              if (highestEnergy is None and lowestEnergy is None):
                 activeImage.Add ( self.ds0[i], scale = ( float ( scan ) / float (interpolations) ) )
                 self.objectiveFunction.VariablesPut (activeImage)
                 highestEnergy = lowestEnergy = system.Energy ( log = None ) 
                 activeImage.CopyTo (EmaxImage)
                 activeImage.CopyTo (EminImage)
                 EmaxScan = "PeakinSegment"
                 EminScan = "ValleyinSegment"
              elif (highestEnergy is not None and lowestEnergy is None):
                  activeImage.Add ( self.ds0[i], scale = ( float ( scan ) / float (interpolations) ) )
                  self.objectiveFunction.VariablesPut (activeImage)
                  tempEnergy = lowestEnergy = system.Energy ( log = None )
                  activeImage.CopyTo (EminImage)
                  EminScan = "ValleyinSegment"
                  if (tempEnergy > highestEnergy):
                    highestEnergy = tempEnergy
                    activeImage.CopyTo (EmaxImage)
                    EmaxScan = "PeakinSegment"
              elif (highestEnergy is None and lowestEnergy is not  None):
                  activeImage.Add ( self.ds0[i], scale = ( float ( scan ) / float (interpolations) ) )
                  self.objectiveFunction.VariablesPut (activeImage)
                  tempEnergy = highestEnergy = system.Energy ( log = None )
                  activeImage.CopyTo (EmaxImage)
                  EmaxScan = "PeakinSegment"
                  if (tempEnergy < lowestEnergy):
                    lowestEnergy = tempEnergy
                    activeImage.CopyTo (EmaxImage)
                    EminScan = "ValleyinSegment"
              elif (highestEnergy is not None and lowestEnergy is not  None):
                 activeImage.Add ( self.ds0[i], scale = ( float ( scan ) / float (interpolations) ) )
                 self.objectiveFunction.VariablesPut (activeImage)
                 tempEnergy = system.Energy ( log = None )
                 if (tempEnergy > highestEnergy):
                    highestEnergy = tempEnergy
                    activeImage.CopyTo (EmaxImage)
                    EmaxScan = "PeakinSegment"
                 if (tempEnergy < lowestEnergy):
                    lowestEnergy = tempEnergy
                    activeImage.CopyTo (EminImage)
                    EminScan = "ValleyinSegment"
              self.objectiveFunction.VariablesPut (EmaxImage)
           # . Checking flanking images.
           if (self.energy[i] is None):
              self.objectiveFunction.VariablesPut (self.variables[i])
              self.energy[i] = system.Energy ( log = None )
           if (self.energy[i + 1] is None):
              self.objectiveFunction.VariablesPut (self.variables[i + 1])
              self.energy[i + 1] = system.Energy ( log = None )
           if (self.energy[i] >= self.energy[i + 1]):
              if (self.energy[i] >= highestEnergy):
                 self.variables[i].CopyTo (EmaxImage)
                 highestEnergy = self.energy[i]
                 EmaxScan = i
              if (self.energy[i + 1] <= lowestEnergy):
                 self.variables[i + 1].CopyTo (EminImage)
                 lowestEnergy = self.energy[i + 1]
                 EminScan = (i + 1)   
           elif (self.energy[i] < self.energy[i + 1]): 
              if (self.energy[i + 1] >= highestEnergy):
                 self.variables[i + 1].CopyTo (EmaxImage)
                 highestEnergy = self.energy[i + 1]
                 EmaxScan = (i + 1)
              if (self.energy[i] <= lowestEnergy):
                 self.variables[i].CopyTo (EminImage)
                 lowestEnergy = self.energy[i]
                 EminScan = i
           # . Store highest segment images.
           tmpHighSeg = [ Array.WithExtent (self.numberOfVariables), None, None ];  tmpHighSeg[0].Set (0.0)
           if (EmaxScan == i):  # . If highest energy image is flanking image start.
              if self.isSaddle[i]:
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}, and it is a saddle point".format (i))
                 tmpHighSeg[1] = "saddle"
              elif self.isStationary[i]:
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}, and it is a stationary point".format (i))
                 tmpHighSeg[1] = "stationary"
              elif self.isUnrefinable[i]:
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}, but it is an unrefinable point".format (i))
                 tmpHighSeg[1] = "unrefStart"
              else:   
                 logFile.Paragraph ("Highest image is initial flanking image of segment {:d}".format (i))
                 tmpHighSeg[1] = "flankStart"
              tmpHighSeg[2] = highestEnergy
              logFile.Paragraph ("E (kJ/mol) = {:f}".format (highestEnergy))
              self.segmentsHighest[i] = tmpHighSeg
           if (EmaxScan == i + 1):  # . If highest energy image is flanking image stop.
              if self.isSaddle[i + 1]:
                  logFile.Paragraph ("Highest image is final flanking image of segment {:d}, and it is a saddle point".format (i))
                  tmpHighSeg[1] = "saddle"
              elif self.isStationary[i + 1]:
                  logFile.Paragraph ("Highest image is final flanking image of segment {:d}, and it is a stationary point".format (i))
                  tmpHighSeg[1] = "stationary"
              elif self.isUnrefinable[i + 1]:
                  logFile.Paragraph ("Highest image is final flanking image of segment {:d}, but it is an unrefinable point".format (i))
                  tmpHighSeg[1] = "unrefStop"
              else:
                  logFile.Paragraph ("Highest image is final flanking image of segment {:d}".format (i))
                  tmpHighSeg[1] = "flankStop"
              tmpHighSeg[2] = highestEnergy
              logFile.Paragraph ("E (kJ/mol) = {:f}".format (highestEnergy))
              self.segmentsHighest[i] = tmpHighSeg
           elif (EmaxScan == "PeakinSegment"):  # . If highest energy image is within the segment.
              logFile.Paragraph ("Highest image is within segment {:d}".format (i))
              tmpHighSeg[1] = "inSegment"
              tmpHighSeg[2] = highestEnergy
              logFile.Paragraph ("E (kJ/mol) = {:f}".format (highestEnergy))
              tmpHighSeg[0] = EmaxImage
              self.segmentsHighest[i] = tmpHighSeg
              
           tmpLowSeg = [ Array.WithExtent (self.numberOfVariables), None, None ];  tmpLowSeg[0].Set (0.0)   
           if (EminScan == i):  # . If lowest energy image is flanking image start.
              if self.isSaddle[i]:
                 logFile.Paragraph ("Lowest image is initial flanking image of segment {:d}, and it is a saddle point".format (i))
                 tmpLowSeg[1] = "saddle"
              elif self.isStationary[i]:
                 logFile.Paragraph ("Lowest image is initial flanking image of segment {:d}, and it is a stationary point".format (i))
                 tmpLowSeg[1] = "stationary"
              elif self.isUnrefinable[i]:
                 logFile.Paragraph ("Lowest image is initial flanking image of segment {:d}, and it is an unrefinable point".format (i))
                 tmpLowSeg[1] = "unrefStart"
              else:   
                 logFile.Paragraph ("Lowest image is initial flanking image of segment {:d}".format (i))
                 tmpLowSeg[1] = "flankStart"
              tmpLowSeg[2] = lowestEnergy
              logFile.Paragraph ("E (kJ/mol) = {:f}".format (lowestEnergy))
              self.segmentsLowest[i] = tmpLowSeg
           if (EminScan == i + 1):  # . If lowest energy image is flanking image stop.
              if self.isSaddle[i + 1]:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}, and it is a saddle point".format (i))
                  tmpLowSeg[1] = "saddle"
              elif self.isStationary[i + 1]:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}, and it is a stationary point".format (i))
                  tmpLowSeg[1] = "stationary"
              elif self.isUnrefinable[i + 1]:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}, and it is an unrefinable point".format (i))
                  tmpLowSeg[1] = "unrefStop"
              else:
                  logFile.Paragraph ("Lowest image is final flanking image of segment {:d}".format (i))
                  tmpLowSeg[1] = "flankStop"
              tmpLowSeg[2] = lowestEnergy
              logFile.Paragraph ("E (kJ/mol) = {:f}".format (lowestEnergy))
              self.segmentsLowest[i] = tmpLowSeg
           elif (EminScan == "ValleyinSegment"):  # . If lowest energy image is within the segment.
              logFile.Paragraph ("Lowest image is within segment {:d}".format (i))
              tmpLowSeg[1] = "inSegment"
              tmpLowSeg[2] = lowestEnergy
              logFile.Paragraph ("E (kJ/mol) = {:f}".format (lowestEnergy))
              tmpLowSeg[0] = EminImage
              self.segmentsLowest[i] = tmpLowSeg
# . MJF.
           system.coordinates3 = temporary  # Restore temporary system coordinates
           if system.freeAtoms is None: self.objectiveFunction.iCoordinates3 = system.coordinates3.iterator
           else:                        self.objectiveFunction.iCoordinates3 = system.coordinates3.RowIterator ( selection = self.freeAtoms )
#

    def StructureToOptimizer (self, system, imageTrajectory, log=logFile):
        """Getting highest structure and passing it to the saddle optimizer."""
        segment = None
        if self.finalUnrefRef:
            i = 0
            unrefinable = 0
            while i < (self.numberOfImages - 1):
                if self.segmentsHighest[i][1] == "unrefStop" or self.segmentsHighest[i][1] == "unrefStart": 
                    table = log.GetTable (columns=[40])  
                    table.Start ()
                    table.Entry ("{:38s}".format ("=== Final unrefinable refinement ==="))
                    table.Stop ()
                    report = self.Optimization (system, imageTrajectory, i)
                    unrefinable += 1
                i += 1
            if unrefinable == 0:
                report = "CPR finished successfully. No points are marked as unrefinable."
        else: 
          segment = self.HighestPoint ()
          if segment == None: 
              self.fullyRefinedPath = True
              report = None
          else:
              table = log.GetTable (columns=[36])  
              table.Start ()
              table.Entry ("{:32s}".format ("==== Starting CPRrun {:d} ====".format (self.CPRrun)))
              table.Stop ()
              if LogFileActive (log):
                  table = log.GetTable (columns=[36])
                  table.Start ()
                  table.Entry   ("{:36s}".format ("Structure to refine:"))
                  if (self.segmentsHighest[segment][1] == "inSegment"):
                      table.Entry ("New path-point within segment {:d} ".format (segment))
                  elif (self.segmentsHighest[segment][1] == "flankStart"):
                      table.Entry ("Existing path point {:d} ".format (segment))
                  elif (self.segmentsHighest[segment][1] == "flankStop"):
                      table.Entry ("Existing path point {:d} ".format (segment + 1))
                  table.Stop ()
              report = self.Optimization (system, imageTrajectory, segment)
        return (report)

    def WriteOutRestartFile (self):
         """Write out the restart file."""
         with open(self.outputRestartFile, 'w') as restartFile:
             restartFile.write("{:>6s}     {:20s} {:8s} {:15s} \n".format ("Image", "Energy[kJ/mol]", "type", "rmsGradient[kJ/mol]"))
             for i in range (0, self.numberOfImages):
                 if self.rmsGradient[i] == None: RMSGradient = ""
                 else: RMSGradient = self.rmsGradient[i]
                 if (i == 0):
                     restartFile.write("{:4d} {:19.6f}   {:9s} {:.6}\n".format (i, self.energy[i], "initial", RMSGradient))
                 elif (i == self.numberOfImages - 1):
                     restartFile.write("{:4d} {:19.6f}   {:9s} {:.6}\n".format (i, self.energy[i], "final", RMSGradient))
                 else:
                     if self.isSaddle[i]: restartFile.write("{:4d} {:19.6f}   {:9s} {:11.6}\n".format (i, self.energy[i], "saddle", RMSGradient))
                     elif self.isStationary[i]:  restartFile.write("{:4d} {:19.6f}   {:9s} {:11.6}\n".format (i, self.energy[i], "stationary", RMSGradient))
                     elif self.isUnrefinable[i]:  restartFile.write("{:4d} {:19.6f}   {:9s} {:11.6}\n".format (i, self.energy[i], "unrefinable", RMSGradient))
                     else:                restartFile.write("{:4d} {:19.6f}   {:9s} {:11.6}\n".format (i, self.energy[i], "normal", RMSGradient))
             restartFile.write("{:7s}    {:19s} {:10s}  {:19s} {:8s}\n".format ("Segment", "Highest[kJ/mol]", "type", "Lowest[kJ/mol]", "type"))
             for i in range (0, self.numberOfImages - 1):
                 restartFile.write("{:4d} {:19.6f}    {:10s} {:19.6f}   {:10s}\n".format (i, self.segmentsHighest[i][2], self.segmentsHighest[i][1], self.segmentsLowest[i][2], self.segmentsLowest[i][1]))

#===================================================================================================================================
# . CPR refinement.
#===================================================================================================================================
# . CPR refinement options.
_CPRRefinementOptions = { "finiteDifferenceStep"          : 1.0e-3,
                          "finiteDifference"              : True  ,
                          "lineSearchInitialStep"         : 0.1   ,
                          "interpolationIncreasement"     :   3   , # . Set preferably odd number
                          "rmsdMaximumTolerance"          :  1.5  ,
                          "inputRestartFile"              : None  ,
                          "outputRestartFile"             : None  ,
                          "stepsPerSegment"               :   3   ,
                          "upperStepLimit"                : None  ,
                          "lowestStepLimit"               : None  ,
                          "maximalMinimizationSteps"      : None  , # . a value around 800 is maybe enough
                          "maxCPRrun"                     : 1000  ,
                          "numberOfSuccessiveLowGradients": None  ,
                          "initialOrthogonalRuns"         :   0   ,
                          "rmsGradientTolerance"          :  0.1  , # . kJ/mol
                          "rmsdLoopPreventionX0"          : 1.0e-6, # . bigger values can cause problems because different points are detected to be the same!!
                          "rmsdLoopPreventionS0"          : 1.0e-6, # . bigger values can cause problems because different points are detected to be the same!!
                          "tauTolAdd"                     :  0.15 ,
                          "tauTolRefine"                  :  0.05 , # . Sinclair and Fletcher tau tolerance (no dimensionality correction as stated in Fischer&Karplus (1992)) 
                          "increaseTau"                   : False ,
                          "breakIfTauReached"             : True  ,
                          "finalUnrefinableRefinement"    : False ,
                          "scratchDirectory"              : None  ,
                          "betaType"                      : "PRP" } # . "FR": Fletcher-Reeves; "PRP": Polak-Ribiere and Polyak; "HS": Hestenes-Stiefel.          

# . Function.
def ConjugatePeakRefinementOptimizePath (system, imageTrajectory, **options):
    """CPR path optimization."""
    # . Get some variables.
    log = options.pop ( "log", logFile )
    # . Create an objective function.
    if isinstance (system, CPRObjectiveFunction): objectiveFunction = system
    else: objectiveFunction = CPRObjectiveFunction.FromSystem (system)
    objectiveFunction.InitializeImages (imageTrajectory)
    if (objectiveFunction.numberOfImages < 2): raise ValueError ("Invalid number of path-points: " + repr ( objectiveFunction.numberOfImages ) + ", need 2 or more to run CPR!")
    # . Optimization.
    toPass = dict ( _CPRRefinementOptions )
    toPass.update ( options ) 
    optimizer = CPRRefinement (**options)
    state = optimizer.CPRIterate (objectiveFunction, system, imageTrajectory, log=log)
    return state

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
