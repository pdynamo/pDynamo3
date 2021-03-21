# ifndef _MNDOPARAMETERS
# define _MNDOPARAMETERS

# include "Boolean.h"
# include "Integer.h"
# include "MachineTypes.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Types.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    /* . Counters. */
    Boolean     QDIATOMIC      ;
    Boolean    *QDIATOMICFLAGS ;
    Integer     atomicNumber   ;
    Integer     iii            ;
    Integer     iiid           ;
    Integer     ir016          ;
    Integer     ir066          ;
    Integer     ir244          ;
    Integer     ir266          ;
    Integer     ir466          ;
    Integer     nam1pm3g       ;
    Integer     ndiatomic      ;
    Integer     norbitals      ;
    Integer     npddg          ;
    Integer     qnd            ;
    Integer     qnp            ;
    Integer     qns            ;
    /* . Derived quantities. */
    Integer     nocteis        ;
    Cardinal16 *octeiindices   ;
    Real        hpp            ;
    Real       *octeivalues    ;
    /* . External (input) values. */
    Real        ad0            ;
    Real        alp0           ;
    Real        am0            ;
    Real        aq0            ;
    Real        betad0         ;
    Real        betap0         ;
    Real        betas0         ;
    Real        dd0            ;
    Real        eheat0         ;
    Real        eisol0         ;
    Real        f0sd0          ;
    Real        gphot0         ;
    Real        gpp0           ;
    Real        gp20           ;
    Real        gsp0           ;
    Real        gss0           ;
    Real        g2sd0          ;
    Real        hsp0           ;
    Real        pcore0         ;
    Real        qq0            ;
    Real        udd0           ;
    Real        upp0           ;
    Real        uss0           ;
    Real        zcore0         ;
    Real        zetad0         ;
    Real        zetap0         ;
    Real        zetas0         ;
    Real        zdn0           ;
    Real        zpn0           ;
    Real        zsn0           ;
    Real       *beta0          ;
    Real       *diatomica0     ;
    Real       *diatomicx0     ;
    Real       *fn10           ;
    Real       *fn20           ;
    Real       *fn30           ;
    Real       *pddgc0         ;
    Real       *pddge0         ;
    Real       *uspd0          ;
    /* . Internal values in atomic units. */
    Real        ad             ;
    Real        alp            ;
    Real        am             ;
    Real        aq             ;
    Real        betad          ;
    Real        betap          ;
    Real        betas          ;
    Real        dd             ;
    Real        eheat          ;
    Real        eisol          ;
    Real        f0sd           ;
    Real        gphot          ;
    Real        gpp            ;
    Real        gp2            ;
    Real        gsp            ;
    Real        gss            ;
    Real        g2sd           ;
    Real        hsp            ;
    Real        pcore          ;
    Real        qq             ;
    Real        udd            ;
    Real        upp            ;
    Real        uss            ;
    Real        zcore          ;
    Real        zetad          ;
    Real        zetap          ;
    Real        zetas          ;
    Real        zdn            ;
    Real        zpn            ;
    Real        zsn            ;
    Real       *beta           ;
    Real       *diatomica      ;
    Real       *diatomicx      ;
    Real       *fn1            ;
    Real       *fn2            ;
    Real       *fn3            ;
    Real       *pddgc          ;
    Real       *pddge          ;
    Real       *uspd           ;
    /* . Arrays. */
    Real        ddp[6]         ;
    Real        po[9]          ;
} MNDOParameters ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern MNDOParameters *MNDOParameters_Allocate                   ( void ) ;
extern void            MNDOParameters_AllocateDiatomic           (       MNDOParameters  *self, const Integer  n ) ;
extern void            MNDOParameters_AllocateFN123              (       MNDOParameters  *self, const Integer  n ) ;
extern void            MNDOParameters_AllocatePDDG               (       MNDOParameters  *self, const Integer  n ) ;
extern MNDOParameters *MNDOParameters_Clone                      ( const MNDOParameters  *self ) ;
extern void            MNDOParameters_CalculateOneCenterTEIs     (       MNDOParameters  *self ) ;
extern void            MNDOParameters_Deallocate                 (       MNDOParameters **self ) ;
extern void            MNDOParameters_DeallocateAtomicUnitArrays (       MNDOParameters  *self ) ;
extern void            MNDOParameters_FillBetaUspd               (       MNDOParameters  *self ) ;
extern void            MNDOParameters_ToAtomicUnits              (       MNDOParameters  *self ) ;

# endif
