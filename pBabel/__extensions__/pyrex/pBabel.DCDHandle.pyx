"""DCD handling functions."""

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Default message.
_DefaultMessage = "DCD trajectory file error."

# . Defined messages.
_Messages = { CDCDStatus_AtomNumberMismatch      : "Atom number mismatch."           ,
              CDCDStatus_BadFormat               : "Bad format."                     ,
              CDCDStatus_BadRead                 : "Read error."                     ,
              CDCDStatus_BadSeek                 : "File positioning error."         ,
              CDCDStatus_BadWrite                : "Write error."                    ,
              CDCDStatus_FileAccessFailure       : "Unable to access file."          ,
              CDCDStatus_InvalidDataObject       : "Invalid or missing data object." ,
              CDCDStatus_InvalidFrameIndex       : "Invalid frame  index."           ,
              CDCDStatus_OpenFailed              : "File opening error."             ,
              CDCDStatus_OutOfMemory             : "Out of memory."                  }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DCDError ( Exception ):
    pass

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
cdef DCDStatus_Check ( CDCDStatus status ):
    """Check the status flag and raise an error if it is not normal."""
    if status != CDCDStatus_Normal: raise DCDError ( _Messages.get ( status, _DefaultMessage ) )
