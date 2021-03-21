"""Classes for manipulating restraints.

This should be reimplemented in C for more efficiency.
"""

import math

from pCore                 import AttributableObject , \
                                  Clone              , \
                                  Selection
from pScientific           import Units
from pScientific.Geometry3 import Coordinates3       , \
                                  Vector3

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RestraintEnergyModel ( AttributableObject ):
    """Restraint energy model."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "highEquilibriumValue" : None ,
                             "highForceConstant"    : None ,
                             "highPower"            : None ,
                             "isPeriodic"           : None ,
                             "lowEquilibriumValue"  : None ,
                             "lowForceConstant"     : None ,
                             "lowPower"             : None ,
                             "period"               : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintEnergyModel, self )._CheckOptions ( )
        isOK = isinstance ( self.lowPower , int ) and isinstance ( self.lowEquilibriumValue , float ) and isinstance ( self.lowForceConstant , float ) and \
               isinstance ( self.highPower, int ) and isinstance ( self.highEquilibriumValue, float ) and isinstance ( self.highForceConstant, float ) and \
               ( ( self.lowPower > 0 ) or ( self.lowForceConstant == 0.0 ) ) and ( ( self.highPower > 0 ) or ( self.highForceConstant == 0.0 ) ) and \
               ( self.lowEquilibriumValue <= self.highEquilibriumValue ) and ( self.lowForceConstant >= 0.0 ) and ( self.highForceConstant >= 0.0 )
        if not isOK: raise TypeError ( "Invalid option." )
        self.isPeriodic = ( self.period is not None )

    def Energy ( self, value ):
        """Calculate an energy and a derivative given a value."""
        energy = 0.0
        dedv   = 0.0
        v      = value
        if self.period is None:
            higher = self.highEquilibriumValue
            lower  = self.lowEquilibriumValue
        else:
            # . Make sure the boundaries are OK.
            higher = self.highEquilibriumValue - math.floor ( self.highEquilibriumValue / self.period ) * self.period
            lower  = self.lowEquilibriumValue  - math.floor ( self.lowEquilibriumValue  / self.period ) * self.period
            if lower > higher: lower -= self.period
            # . Move the value just to the left of lower.
            v = v - math.floor ( v / self.period ) * self.period
            while v > lower: v -= self.period
            # . Three possible cases.
            v0 = v + self.period
            if ( ( v0 >= lower ) and ( v0 <= higher ) ) or ( ( v0 > higher ) and ( ( v0 - higher ) < ( lower - v ) ) ): v = v0
        if v < lower:
            disp   = lower - v
            dedv   = self.lowForceConstant * disp**( self.lowPower-1 )
            energy = dedv * disp
            dedv  *= - float ( self.lowPower )
        elif v > higher:
            disp   = v - higher
            dedv   = self.highForceConstant * disp**( self.highPower-1 )
            energy = dedv * disp
            dedv  *= float ( self.highPower )
        return ( energy, dedv )

    @classmethod
    def FlatBottomedHarmonic ( selfClass, lowEquilibriumValue, highEquilibriumValue, forceConstant, period = None ):
        """Define a flat-bottomed harmonic energy model."""
        return selfClass.WithOptions ( highEquilibriumValue = highEquilibriumValue ,
                                       highForceConstant    = forceConstant        ,
                                       highPower            = 2                    ,
                                       lowEquilibriumValue  = lowEquilibriumValue  ,
                                       lowForceConstant     = forceConstant        ,
                                       lowPower             = 2                    ,
                                       period               = period               )

    @classmethod
    def Harmonic ( selfClass, equilibriumValue, forceConstant, period = None ):
        """Define a harmonic energy model."""
        return selfClass.WithOptions ( highEquilibriumValue = equilibriumValue ,
                                       highForceConstant    = forceConstant    ,
                                       highPower            = 2                ,
                                       lowEquilibriumValue  = equilibriumValue ,
                                       lowForceConstant     = forceConstant    ,
                                       lowPower             = 2                ,
                                       period               = period           )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Restraint ( AttributableObject ):
    """Base class for restraint types."""

    _attributable = dict ( AttributableObject._attributable )
    _classLabel   = "Undefined"
    _period       = None # . Period of the class's restraint variable. None implies no periodicity.
    _attributable.update ( { "energyModel" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( Restraint, self )._CheckOptions ( )
        # . Check the type and period of the energy model.
        if isinstance ( self.energyModel, RestraintEnergyModel ):
            period     = self.__class__._period
            isPeriodic = ( period is not None )
            isOK       = ( ( period is None ) and ( not self.energyModel.isPeriodic ) ) or \
                         ( isPeriodic and self.energyModel.isPeriodic and ( period == self.energyModel.period ) )
            if not isOK: raise ValueError ( "The restraint and energy model have incompatible periodicites." )
        else: raise TypeError ( "Invalid restraint energy model." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        return ( 0.0, 0.0 )

    def Merge ( self, indexIncrement ):
        """Merging."""
        return None

    def Prune ( self, selection ):
        """Pruning."""
        return None

#===================================================================================================================================
# . Subclasses.
#===================================================================================================================================
class RestraintAngleDotProduct ( Restraint ):
    """An angle dot product restraint.

    The dot product of the angle is restrained so the reference value
    must be between -1 and 1.
    """

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Angle Dot Product"
    _attributable.update ( { "point1" : None ,
                             "point2" : None ,
                             "point3" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintAngleDotProduct, self )._CheckOptions ( )
        if not ( isinstance ( self.point1, int ) and \
                 isinstance ( self.point2, int ) and \
                 isinstance ( self.point3, int ) ):
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        v12 = coordinates3.Displacement ( self.point1, self.point2 ) ; r12 = v12.Norm2 ( )
        v32 = coordinates3.Displacement ( self.point3, self.point2 ) ; r32 = v32.Norm2 ( )
        v12.Scale ( 1.0 / r12 )
        v32.Scale ( 1.0 / r32 )
        value = v12.Dot ( v32 )
        if   value >  1.0: value =  1.0
        elif value < -1.0: value = -1.0
        ( f, df ) = self.energyModel.Energy ( value )
        if gradients3 is not None:
            v32.Add ( v12, scale = -value ) # . This is s12 * d(costheta)/d(p1) in d32.
            v12.Scale ( 1.0 - value**2 )
            v12.Add ( v32, scale = -value ) # . This is s32 * d(costheta)/d(p3) in d12.
            v32.Scale ( df / r12 )
            v12.Scale ( df / r32 )
            gradients3.AddScaledRow ( self.point1,  1.0, v32 )
            gradients3.AddScaledRow ( self.point2, -1.0, v12 )
            gradients3.AddScaledRow ( self.point2, -1.0, v32 )
            gradients3.AddScaledRow ( self.point3,  1.0, v12 )
        return ( f, value )

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        new.point3 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ) and \
           ( self.point3 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), \
                                      selection.Position ( self.point3 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class RestraintDihedral ( Restraint ):
    """A dihedral angle restraint."""

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Dihedral"
    _period       = 360.0
    _attributable.update ( { "point1" : None ,
                             "point2" : None ,
                             "point3" : None ,
                             "point4" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintDihedral, self )._CheckOptions ( )
        if not ( isinstance ( self.point1, int ) and \
                 isinstance ( self.point2, int ) and \
                 isinstance ( self.point3, int ) and \
                 isinstance ( self.point4, int ) ):
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        # . Displacements.
        v12 = coordinates3.Displacement ( self.point1, self.point2 )
        v32 = coordinates3.Displacement ( self.point3, self.point2 )
        v34 = coordinates3.Displacement ( self.point3, self.point4 )
        # . Get m and n.
        m = Clone ( v12 ) ; m.Cross ( v32 )
        n = Clone ( v32 ) ; n.Cross ( v34 )
        # . Get the sizes of m and n.
        msize = m.Norm2 ( )
        nsize = n.Norm2 ( )
        # . Normalize m and n.
        m.Scale ( 1.0 / msize )
        n.Scale ( 1.0 / nsize )
        # . Get the dot-product.
        dotFac = m.Dot ( n )
        if    dotFac > 1.0: dotFac =  1.0
        elif dotFac < -1.0: dotFac = -1.0
        # . Get the sign of the angle.
        sgnfac = 1.0
        if v12.Dot ( n ) < 0.0: sgnfac = -1.0
        # . Determine the dihedral.
        value = sgnfac * math.acos ( dotFac ) * Units.Angle_Radians_To_Degrees
        # . Get the energy.
        ( f, df ) = self.energyModel.Energy ( value )
        # . Derivatives.
        if gradients3 is not None:
            df *= Units.Angle_Radians_To_Degrees
            # . Calculate r32.
            r32 = v32.Norm2 ( )
            # . Calculate dedi and dedl in m and n respectively.
            m.Scale (   df * r32 / msize )
            n.Scale ( - df * r32 / nsize )
            # . Calculate some additional factors.
            fact12 = v12.Dot ( v32 ) / ( r32 * r32 )
            fact34 = v34.Dot ( v32 ) / ( r32 * r32 )
            # . Gradients for i and l.
            gradients3.AddScaledRow ( self.point1, 1.0, m )
            gradients3.AddScaledRow ( self.point4, 1.0, n )
            # . Calculate dedj and dedk in v12 and v32 respectively.
            m.CopyTo ( v12 ) ; v12.Scale ( fact12 - 1.0 ) ; v12.Add ( n, scale = -fact34 )
            n.CopyTo ( v32 ) ; v32.Scale ( fact34 - 1.0 ) ; v32.Add ( m, scale = -fact12 )
            # . calculate the gradients.
            gradients3.AddScaledRow ( self.point2, 1.0, v12 )
            gradients3.AddScaledRow ( self.point3, 1.0, v32 )
        return ( f, value )

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        new.point3 += indexIncrement
        new.point4 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ) and \
           ( self.point3 in selection ) and \
           ( self.point4 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), \
                                      selection.Position ( self.point3 ), \
                                      selection.Position ( self.point4 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class RestraintDistance ( Restraint ):
    """A distance restraint."""

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Distance"
    _attributable.update ( { "point1" : None ,
                             "point2" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintDistance, self )._CheckOptions ( )
        if not ( isinstance ( self.point1, int ) and \
                 isinstance ( self.point2, int ) ):
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        v12 = coordinates3.Displacement ( self.point1, self.point2 )
        r12 = v12.Norm2 ( )
        ( f, df ) = self.energyModel.Energy ( r12 )
        if gradients3 is not None:
            v12.Scale ( df / r12 )
            gradients3.AddScaledRow ( self.point1,  1.0, v12 )
            gradients3.AddScaledRow ( self.point2, -1.0, v12 )
        return ( f, r12 )

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class RestraintMultipleDistance ( Restraint ):
    """A multiple distance restraint."""

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Multiple Distance"
    _attributable.update ( { "distances" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintMultipleDistance, self )._CheckOptions ( )
        try:
            for ( point1, point2, weight ) in self.distances:
                if not ( isinstance ( point1, int ) and \
                         isinstance ( point2, int ) and \
                         isinstance ( weight, float ) ): raise
        except:
            raise TypeError ( "Invalid distances specification." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        value = 0.0
        temp  = []
        for ( point1, point2, weight ) in self.distances:
            v = coordinates3.Displacement ( point1, point2 )
            r = v.Norm2 ( )
            value += weight * r
            temp.append ( ( v, r ) )
        ( f, df ) = self.energyModel.Energy ( value )
        if gradients3 is not None:
            for ( ( point1, point2, weight ), ( v, r ) ) in zip ( self.distances, temp ):
                v.Scale ( df * weight / r )
                gradients3.AddScaledRow ( point1,  1.0, v )
                gradients3.AddScaledRow ( point2, -1.0, v )
        return ( f, value )

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = None
        # . Loop over distances.
        distances = []
        for ( point1, point2, weight ) in self.distances:
            distances.append ( ( point1 + indexIncrement, point2 + indexIncrement, weight ) )
        # . Create object.
        if len ( distances ) > 0:
            new = self.__class__ ( distances, Clone ( self.energyModel ) )
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        # . Loop over distances.
        reduced = []
        for ( point1, point2, weight ) in self.distances:
            if ( point1 in selection ) and ( point2 in selection ):
                reduced.append ( ( selection.Position ( point1 ), selection.Position ( point2 ), weight ) )
        # . Create object.
        if len ( reduced ) > 0:
            pruned = self.__class__ ( reduced, Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class RestraintMultipleTether ( Restraint ):
    """A multiple tether restraint."""

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Multiple Tether"
    _attributable.update ( { "reference" : None ,
                             "selection" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintMultipleTether, self )._CheckOptions ( )
        if not ( isinstance ( self.selection, Selection    ) and \
                 isinstance ( self.reference, Coordinates3 ) ):
            raise TypeError ( "Invalid reference coordinates or selection." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        fTotal = 0.0
        v12    = Vector3.Null ( )
        for i in self.selection:
            coordinates3[i].CopyTo ( v12 )
            v12.Add ( self.reference[i], scale = -1.0 )
            r12 = v12.Norm2 ( )
            ( f, dF ) = self.energyModel.Energy ( r12 )
            fTotal += f
            if ( gradients3 is not None ) and ( r12 != 0.0 ):
                v12.Scale ( dF / r12 )
                gradients3.AddScaledRow ( i, 1.0, v12 )
        return ( fTotal, None )

    def Merge ( self, indexIncrement ):
        """Merging."""
        indices = [ i + indexIncrement for i in self.selection ]
        return self.__class__ ( Selection.FromIterable ( indices ) ,
                                        Clone ( self.reference   ) ,
                                        Clone ( self.energyModel ) )
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        # . Get reduced selection.
        reduced = self.selection.Prune ( selection )
        if len ( reduced ) > 0:
            pruned = self.__class__ ( reduced, self.reference.Prune ( selection ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class RestraintOutOfPlaneBend ( Restraint ):
    """An out-of-plane bend restraint.

    Angle is between 1-2-3 plane and 1-4 vector for 1-4<23 bonding pattern.
    """

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Out-Of-Plane Bend"
    _attributable.update ( { "point1" : None ,
                             "point2" : None ,
                             "point3" : None ,
                             "point4" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintOutOfPlaneBend, self )._CheckOptions ( )
        if not ( isinstance ( self.point1, int ) and \
                 isinstance ( self.point2, int ) and \
                 isinstance ( self.point3, int ) and \
                 isinstance ( self.point4, int ) ):
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        # . Displacements.
        t12 = coordinates3.Displacement ( self.point1, self.point2 )
        t32 = coordinates3.Displacement ( self.point3, self.point2 )
        t34 = coordinates3.Displacement ( self.point3, self.point4 )
        t43 = coordinates3.Displacement ( self.point4, self.point3 )
        t41 = coordinates3.Displacement ( self.point4, self.point1 )
        t42 = coordinates3.Displacement ( self.point4, self.point2 )
        tm  = Clone ( t32 ) ; tm.Cross ( t12 )
        tn  = Clone ( tm  ) ; tn.Cross ( t41 )
        sm  = tm.Norm2  ( ) ; ism  = 1.0 / sm
        sn  = tn.Norm2  ( ) ; isn  = 1.0 / sn
        s32 = t32.Norm2 ( ) ; is32 = 1.0 / s32
        # . Factors from the paper.
        u41 =   t41.DotSelf (    ) ; is41 = 1.0 / math.sqrt ( u41 )
        E   = - t41.Dot     ( tm )
        C   =   tm.DotSelf  (    )
        D   =   E * E / C
        B   =   u41 - D
        # . Angle.
        cosine = math.sqrt ( B ) * is41
        cosine = min ( 1.0, max ( -1.0, cosine ) )
        thetaR = math.fabs ( math.acos ( cosine ) )
        theta  = thetaR * Units.Angle_Radians_To_Degrees
        if E < 0.0: sign = -1.0
        else:       sign =  1.0
        # . Energy.
        ( f, dF ) = self.energyModel.Energy ( sign * theta )
        # . Derivatives.
        if gradients3 is not None:
            dF   *= ( - Units.Angle_Radians_To_Degrees * is41 )
            dEdq1 = Clone ( t42 ) ; dEdq1.Cross ( t43 )
            dEdq2 = Clone ( t43 ) ; dEdq2.Cross ( t41 )
            dEdq3 = Clone ( t41 ) ; dEdq3.Cross ( t42 )
            if theta == 0.0:
                dF *= sign / math.sqrt ( C )
                g1  = dEdq1
                g2  = dEdq2
                g3  = dEdq3
            else:
                t42.CopyTo ( t12 ) ; t12.Add ( t41, scale = -1.0 ) # t42-t41
                t42.CopyTo ( t32 ) ; t32.Add ( t43, scale = -1.0 ) # t42-t43
                factor = t12.Dot ( t32 )
                t12.CopyTo ( tm ) ; tm.Scale ( 2.0 * s32 * s32 )       ; tm.Add ( t32, scale = -2.0 * factor )
                t32.CopyTo ( tn ) ; tn.Scale ( 2.0 * t12.DotSelf ( ) ) ; tn.Add ( t12, scale = -2.0 * factor )
                sine   = math.sin ( thetaR )
                factor = 2.0 * E / C
                g1     = dEdq1 ; g1.Scale ( factor )
                g2     = dEdq2 ; g2.Scale ( factor )
                g3     = dEdq3 ; g3.Scale ( factor )
                factor = - D / C
                g1.Add ( tm, scale = factor ) ; g2.Add ( tm, scale = -factor ) # tm = dCdq1
                g3.Add ( tn, scale = factor ) ; g2.Add ( tn, scale = -factor ) # tn = dCdq3
                g1.Add ( t41, scale = 2.0 )
                factor = - 0.5 / ( math.sqrt ( B ) * sine )
                g1.Scale ( factor ) ; g2.Scale ( factor ) ; g3.Scale ( factor )
                factor = cosine * is41 / sine
                g1.Add ( t41, scale = factor )
            g1.Scale ( dF ) ; g2.Scale ( dF ) ; g3.Scale ( dF )
            gradients3.AddScaledRow ( self.point1, 1.0, g1 ) ; gradients3.AddScaledRow ( self.point4, -1.0, g1 )
            gradients3.AddScaledRow ( self.point2, 1.0, g2 ) ; gradients3.AddScaledRow ( self.point4, -1.0, g2 )
            gradients3.AddScaledRow ( self.point3, 1.0, g3 ) ; gradients3.AddScaledRow ( self.point4, -1.0, g3 )
        # . Finish up.
#        print ( "SC>", f, theta, thetaR )
        return ( f, theta )

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        new.point3 += indexIncrement
        new.point4 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ) and \
           ( self.point3 in selection ) and \
           ( self.point4 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), \
                                      selection.Position ( self.point3 ), \
                                      selection.Position ( self.point4 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class RestraintTether ( Restraint ):
    """A tether restraint."""

    _attributable = dict ( Restraint._attributable )
    _classLabel   = "Tether"
    _attributable.update ( { "origin" : None ,
                             "point"  : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RestraintTether, self )._CheckOptions ( )
        if not ( isinstance ( self.point , int     ) and \
                 isinstance ( self.origin, Vector3 ) ):
            raise TypeError ( "Invalid point or origin." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        v12 = Vector3.Null ( )
        coordinates3[self.point].CopyTo ( v12 )
        v12.Add ( self.origin, scale = -1.0 )
        r12 = v12.Norm2 ( )
        ( f, df ) = self.energyModel.Energy ( r12 )
        if ( gradients3 is not None ) and ( r12 != 0.0 ):
            v12.Scale ( df / r12 )
            gradients3.AddScaledRow ( self.point, 1.0, v12 )
        return ( f, r12 )

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if self.point in selection:
            pruned = self.__class__ ( selection.Position ( self.point ), Clone ( self.origin ), Clone ( self.energyModel ) )
        return pruned

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
