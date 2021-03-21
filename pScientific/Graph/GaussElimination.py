"""Gauss elimination - boolean vectors only."""

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def GaussElimination ( vectors, yPosition, xPosition ):
    """Gauss elimination."""
    numberVectors = len ( vectors    )
    vectorLength  = len ( vectors[0] )
    if ( yPosition < numberVectors ) and ( xPosition < vectorLength ):
        current = yPosition
        try:
            # . Search for the first vector with xPosition set.
            while ( current < numberVectors ) and ( not vectors[current][xPosition] ):
                current += 1
            # . None found.
            if ( current >= numberVectors ):
                xPosition += 1
                return GaussElimination ( vectors, yPosition, xPosition )
            # . Change positions of these two lines.
            if current != yPosition:
                vector = vectors[current]
                vectors[current]   = vectors[yPosition]
                vectors[yPosition] = vector
        except Exception as e:
            print ( "Gauss elimination error: {:s}.".format ( e[0] ) )
            del vectors[yPosition:]
            return len ( vectors )

        # . XOR all vectors with xPosition set.
        for current in range ( numberVectors ):
            vector = vectors[current]
            if vector[xPosition] and ( current != yPosition ):
                for i in range ( vectorLength ):
                    vector[i] ^= vectors[yPosition][i]
        yPosition += 1
        xPosition += 1
        return GaussElimination ( vectors, yPosition, xPosition )
    # . Finish.
    else:
        # . Remove all empty vectors.
        del vectors[yPosition:]
        return len ( vectors )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
