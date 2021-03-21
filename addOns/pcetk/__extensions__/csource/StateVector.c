#include "Memory.h"
#include "StateVector.h"

/*
 * Allocate the state vector.
 */
StateVector *StateVector_Allocate ( const Integer nsites, Status *status )
{
    StateVector *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( StateVector ) ;
        if ( self != NULL )
        {
            self->sites         = NULL  ;
            self->substateSites = NULL  ;
            self->pairs         = NULL  ;
            self->nsites        = 0     ;
            self->nssites       = 0     ;
            self->npairs        = 0     ;
            if ( nsites > 0 )
            {
                self->sites = Memory_AllocateArrayOfTypes ( nsites, TitrSite );
                if ( self->sites == NULL ) StateVector_Deallocate ( &self ) ;
                else self->nsites = nsites ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self;
}

/*
 * Allocate a substate within the state vector.
 */
void StateVector_AllocateSubstate ( StateVector *self, const Integer nssites, Status *status )
{
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        Memory_Deallocate ( self->substateSites ) ;
        self->nssites = 0 ;
        if ( nssites > 0 )
        {
            self->substateSites = Memory_AllocateArrayOfReferences ( nssites, TitrSite ) ;
            if ( self->substateSites != NULL ) self->nssites = nssites ;
            else Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
}

/*
 * Allocate an array of pairs within the state vector.
 * The number of pairs and their contents are decided by the MCModelDefault module.
 */
void StateVector_AllocatePairs ( StateVector *self, const Integer npairs, Status *status )
{
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        Memory_Deallocate ( self->pairs ) ;
        self->npairs = 0 ;
        if ( npairs > 0 )
        {
            self->pairs = Memory_AllocateArrayOfTypes ( npairs, PairSite ) ;
            if ( self->pairs != NULL ) self->npairs = npairs ;
            else Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
}

/*
 * Deallocate the old pairs and allocate the new ones.
 */
void StateVector_ReallocatePairs ( StateVector *self, const Integer npairs, Status *status )
{
    StateVector_AllocatePairs ( self, npairs, status ) ;
}

/*
 * Deallocate the state vector, including the optional arrays of pairs and substate.
 */
void StateVector_Deallocate ( StateVector **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Memory_Deallocate ( (*self)->pairs         ) ;
        Memory_Deallocate ( (*self)->substateSites ) ;
        Memory_Deallocate ( (*self)->sites         ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*
 * Clone a state vector.
 */
StateVector *StateVector_Clone ( const StateVector *self, Status *status )
{
    StateVector *clone = NULL;
    if ( self != NULL )
    {
        clone = StateVector_Allocate ( self->nsites, status ) ;
        StateVector_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*
 * Copy a state vector to another vector.
 */
void StateVector_CopyTo ( const StateVector *self, StateVector *other, Status *status )
{
    if ( Status_IsOK ( status ) ) Memory_CopyTo ( self->sites, other->sites, self->nsites * sizeof ( TitrSite ) ) ;
}

/*
 * Set all sites of the vector to their initial instances.
 */
void StateVector_Reset ( const StateVector *self )
{
    if ( ( self        != NULL ) &&
         ( self->sites != NULL ) )
    {
        TitrSite *site = self->sites;
        Integer   i = self->nsites;
        for (; i > 0; i--, site++) { site->indexActive = site->indexFirst ; }
    }
}

/*
 * Set all sites of the substate to their initial instances.
 */
void StateVector_ResetSubstate ( const StateVector *self )
{
    if ( ( self != NULL ) && ( self->nssites > 0 ) )
    {
        TitrSite *site, **pointToSite = self->substateSites;
        Integer   i = self->nssites;
        for (; i > 0; i--, pointToSite++)
        {
            site = *pointToSite;
            site->indexActive = site->indexFirst;
        }
    }
}

/*
 * Set all sites of the vector to their final instances.
 */
void StateVector_ResetToMaximum (const StateVector *self)
{
    if ( ( self != NULL ) && ( self->sites != NULL ) )
    {
        TitrSite *site = self->sites;
        Integer   i = self->nsites;
        for (; i > 0; i--, site++) site->indexActive = site->indexLast ;
    }
}

/*
 * Set all sites of the vector to randomized instances.
 */
void StateVector_Randomize ( const StateVector *self, const RandomNumberGenerator *generator )
{
    if ( ( self        != NULL ) &&
         ( generator   != NULL ) &&
         ( self->sites != NULL ) )
    {
        TitrSite *site = self->sites;
        Integer   i = self->nsites;

        for (; i > 0; i--, site++) {
            site->indexActive = RandomNumberGenerator_NextCardinal (generator) % (site->indexLast - site->indexFirst + 1) + site->indexFirst;
        }
    }
}

/*
 * Set a state vector site.
 */
void StateVector_SetSite (const StateVector *self, const Integer indexSite, 
                          const Integer indexFirst, const Integer indexLast, Status *status) {
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( indexSite < 0 || indexSite >= self->nsites ) { Status_Set ( status, Status_IndexOutOfRange ); }
        else
        {
            TitrSite *site;
            site = &self->sites[indexSite];
            site->isSubstate   =  False       ;
            site->indexSite    =  indexSite   ;
            site->indexLast    =  indexLast   ;
            site->indexFirst   =  indexFirst  ;
            site->indexActive  =  indexFirst  ;
        }
    }
}

/*
 * Set a pair of strongly interacting sites.
 */
void StateVector_SetPair (const StateVector *self, const Integer indexPair, 
                          const Integer indexFirstSite, const Integer indexSecondSite, 
                          const Real Wmax, Status *status)
{
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( indexPair < 0 || indexPair >= self->npairs ) { Status_Set ( status, Status_IndexOutOfRange ); }
        else
        {
            PairSite *pair;
            pair       = &self->pairs[indexPair];
            pair->a    = &self->sites[indexFirstSite];
            pair->b    = &self->sites[indexSecondSite];
            pair->Wmax = Wmax;
        }
    }
}

/*
 * Get indices and maximum interaction energy of a pair of strongly interacting sites.
 */
void StateVector_GetPair (const StateVector *self, const Integer indexPair, 
                          Integer *indexFirstSite, Integer *indexSecondSite, Real *Wmax, 
                          Status *status)
{
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( indexPair < 0 || indexPair >= self->npairs ) { Status_Set ( status, Status_IndexOutOfRange ); }
        else
        {
            PairSite *pair ;
            TitrSite *site ;
            pair  = &self->pairs[indexPair];
            site  = pair->a;
            if ( indexFirstSite  != NULL ) *indexFirstSite  = site->indexSite;
            site  = pair->b;
            if ( indexSecondSite != NULL ) *indexSecondSite = site->indexSite ;
            if ( Wmax            != NULL ) *Wmax            = pair->Wmax ;
        }
    }
}

/*
 * Return true if the site belongs to a substate.
 */
Boolean StateVector_IsSubstate (const StateVector *self, const Integer siteIndex, Status *status)
{
    Boolean isSubstate = False ;
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( siteIndex < 0 || siteIndex >= self->nsites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            TitrSite *site;
            site = &self->sites[siteIndex];
            isSubstate = site->isSubstate ;
        }
    }
    return isSubstate ;
}

/*
 * Get the current protonation of a site, i.e. the local index of its currently "active" instance.
 */
Integer StateVector_GetItem ( const StateVector *self, const Integer siteIndex, Status *status )
{
    Integer instanceLocalIndex = -1;
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( siteIndex < 0 || siteIndex >= self->nsites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            TitrSite *site;
            site = &self->sites[siteIndex];
            instanceLocalIndex = site->indexActive - site->indexFirst;
        }
    }
    return instanceLocalIndex;
}

/*
 * Set the protonation of a site by defining a local index of its "active" instance.
 */
void StateVector_SetItem ( const StateVector *self, const Integer siteIndex, const Integer instanceLocalIndex, Status *status )
{
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( siteIndex < 0 || siteIndex >= self->nsites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            TitrSite *site;
            Integer instanceGlobalIndex;
            site = &self->sites[siteIndex];
            /* Translate local index to global index */
            instanceGlobalIndex = instanceLocalIndex + site->indexFirst;
            if ( instanceGlobalIndex < site->indexFirst || instanceGlobalIndex > site->indexLast ) Status_Set ( status, Status_InvalidArgument ) ;
            else site->indexActive = instanceGlobalIndex ;
        }
    }
}

/*
 * Get the current protonation of a site, i.e. the global index of its currently "active" instance.
 */
Integer StateVector_GetActualItem ( const StateVector *self, const Integer siteIndex, Status *status )
{
    Integer instanceGlobalIndex = -1;
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( siteIndex < 0 || siteIndex >= self->nsites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            TitrSite *site;
            site = &self->sites[siteIndex];
            instanceGlobalIndex = site->indexActive;
        }
    }
    return instanceGlobalIndex;
}

/*
 * Set the protonation of a site by defining a global index of its "active" instance.
 */
void StateVector_SetActualItem ( const StateVector *self, const Integer siteIndex, const Integer instanceGlobalIndex, 
                                Status *status )
{
    if ( ( self != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( siteIndex < 0 || siteIndex >= self->nsites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            TitrSite *site;
            site = &self->sites[siteIndex];
            if ( instanceGlobalIndex < site->indexFirst || instanceGlobalIndex > site->indexLast ) Status_Set ( status, Status_InvalidArgument ) ;
            else site->indexActive = instanceGlobalIndex ;
        }
    }
}

/*
 * Get the index of a site belonging to a substate.
 */
Integer StateVector_GetSubstateItem ( const StateVector *self, const Integer index, Status *status )
{
    Integer indexSite = -1 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {

        if ( index < 0 || index >= self->nssites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            TitrSite *site;
            site = self->substateSites[index];
            indexSite = site->indexSite ;
        }
    }
    return indexSite;
}

/*
 * Attach the selected site to a substate by passing its index.
 */
void StateVector_SetSubstateItem (const StateVector *self, const Integer selectedSiteIndex, const Integer index, 
                                  Status *status)
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( index < 0 || index >= self->nssites ) { Status_Set ( status, Status_IndexOutOfRange ) ; }
        else
        {
            if ( selectedSiteIndex < 0 || selectedSiteIndex >= self->nsites ) { Status_Set ( status, Status_InvalidArgument ) ; }
            else
            {
                TitrSite *site;
                site = &self->sites[selectedSiteIndex];
                site->isSubstate = True;
                self->substateSites[index] = site;
            }
        }
    }
}

/*
 * Increment the state vector.
 * After reaching the last vector, false is returned and the vector is back in its initial state.
 * True is returned as long as there are more vectors ahead.
 *
 * Incrementation algorithm by Timm Essigke.
 */
Boolean StateVector_Increment ( const StateVector *self )
{
    Boolean done = False ;
    if ( ( self != NULL ) && ( self->sites != NULL ) )
    {
        TitrSite *site = self->sites;
        Integer   i    = self->nsites;
        for (; i > 0; i--, site++) {
            if (site->indexActive < site->indexLast) {
                site->indexActive++;
                done = True ;
                break ;
            }
            else {
                site->indexActive = site->indexFirst;
            }
        }
    }
    return done ;
}

/*
 * Increment only within the substate of sites of the vector.
 */
Boolean StateVector_IncrementSubstate ( const StateVector *self )
{
    Boolean done = False ;
    if ( ( self != NULL ) && ( self->substateSites != NULL ) )
    {
        TitrSite *site, **pointToSite = self->substateSites;
        Integer   i = self->nssites;
        for (; i > 0; i--, pointToSite++) {
            site = *pointToSite;
            if (site->indexActive < site->indexLast) {
                site->indexActive++;
                done = True ;
                break ;
            }
            else {
                site->indexActive = site->indexFirst;
            }
        }
    }
    return done ;
}
