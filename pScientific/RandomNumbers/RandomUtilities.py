"""Random utilities."""

import random, string

#===================================================================================================================================
# . Methods.
#===================================================================================================================================
def RandomString ( characters = ( string.ascii_lowercase + string.digits ), size = 12, startingCharacters = None ):
    """Generate a random string."""
    if startingCharacters is None: startingCharacters = string.ascii_lowercase
    return "".join ( [ random.choice ( startingCharacters ) ] + [ random.choice ( characters ) for x in range ( size - 1 ) ] )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
