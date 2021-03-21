/*==================================================================================================================================
! . This module defines terms for MNDO integral evaluation.
!=================================================================================================================================*/

/* . Remember that all "const" declarations can only be specified in one .c file! */

# include "MNDODefinitions.h"
# include "MNDOIntegralDefinitions.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Angular momentum of orbitals (s, p, d).
!---------------------------------------------------------------------------------------------------------------------------------*/
const Integer ORBITALAM[9] = { 0, 1, 1, 1, 2, 2, 2, 2, 2 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Core integrals in the local frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . sp orbitals removed from unique terms in CUNIQUE and NCUNIQUE. */
/* . Counters. */
const Integer  NCPOSITIVE [3] = { 0, 1,  4 } ;
const Integer  NCUNIQUE   [3] = { 0, 0,  6 } ;
const Integer _NCUniqueSPD[3] = { 1, 4, 10 } ;

/* . Integrals. */
/* . Unique indices must be such that i >= j. */
const Integer  CPOSITIVE [ 8] = { PYPY, PXPX,   DYZPY, DXZPX,   DYZDYZ, DXZDXZ,   DXYDXY, DX2Y2DX2Y2 } ;
const Integer  CUNIQUE   [12] = { DZ2, S,   DZ2, PZ,   DZ2, DZ2,   DXZ, PX,   DXZ, DXZ,   DX2Y2, DX2Y2 } ;
const Integer _CUniqueSPD[20] = { S, S,   PZ, S,   PZ, PZ,   PX, PX,   DZ2, S,   DZ2, PZ,   DZ2, DZ2,   DXZ, PX,   DXZ, DXZ,   DX2Y2, DX2Y2 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Two-electron integrals in the local frame.
! . sp orbitals removed from unique terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Counters. */
const Integer  NNEGATIVE[3][3] = { {  0,  0,   0 },
                                   {  0,  0,   7 },
                                   {  0,  7,  50 } } ;
const Integer  NPOSITIVE[3][3] = { {  0,  1,   4 },
                                   {  1, 12,  52 },
                                   {  4, 52, 237 } } ;
/*const Integer  NUNIQUE  [3][3] = { {  1,  4,  10 },
                                     {  4, 22,  64 },
                                     { 10, 64, 204 } } ;*/
const Integer  NUNIQUE  [3][3] = { {  0,  0,   6 },
                                   {  0,  0,  42 },
                                   {  6, 42, 182 } } ;

/* . Integrals. */
/* . Unique indices must be such that i >= j and k >= l. */
/* . Second atom is s. */
const Integer  SNEGATIVE[ 0] = {} ;
const Integer  SPOSITIVE[16] = {
     PYPY      , SS        , PXPX      , SS        ,       DYZPY     , SS        , DXZPX     , SS        ,       DYZDYZ    , SS        , DXZDXZ    , SS        ,
     DXYDXY    , SS        , DX2Y2DX2Y2, SS } ;
/*
const Integer  SUNIQUE[40] = {
     S         , S         , S         , S         ,       PZ        , S         , S         , S         ,       PZ        , PZ        , S         , S         ,       PX        , PX        , S         , S         ,       DZ2       , S         , S         , S ,
     DZ2       , PZ        , S         , S         ,       DZ2       , DZ2       , S         , S         ,       DXZ       , PX        , S         , S         ,       DXZ       , DXZ       , S         , S         ,       DX2Y2     , DX2Y2     , S         , S } ;
*/
const Integer  SUNIQUE[24] = {
     DZ2       , S         , S         , S         ,
     DZ2       , PZ        , S         , S         ,       DZ2       , DZ2       , S         , S         ,       DXZ       , PX        , S         , S         ,       DXZ       , DXZ       , S         , S         ,       DX2Y2     , DX2Y2     , S         , S } ;

/* . Second atom is sp. */
const Integer  SPNEGATIVE[28] = {
      DX2Y2S    , PYPY      , DX2Y2S    , PXPX      ,       DX2Y2PZ   , PYPY      , DXZPY     , PYPX      ,       DX2Y2PY   , PYS       , DX2Y2PX   , PXS       ,
      DX2Y2PY   , PYPZ      , DX2Y2PX   , PXPZ      ,       DX2Y2DZ2  , PYPY      , DX2Y2DZ2  , PXPX      ,       DX2Y2DYZ  , PYS       , DX2Y2DXZ  , PXS       ,
      DX2Y2DYZ  , PYPZ      , DX2Y2DXZ  , PXPZ } ;
const Integer  SPPOSITIVE[208] = {
      SS        , PYPY      , SS         , PXPX      ,       PZS       , PYPY      , PZS        , PXPX      ,       PZPZ       , PYPY      , PZPZ       , PXPX      ,
      PYS       , PYS       , PXS        , PXS       ,       PYS       , PYPZ      , PXS        , PXPZ      ,       PYPZ       , PYS       , PXPZ       , PXS       ,
      PYPZ      , PYPZ      , PXPZ       , PXPZ      ,       PYPY      , SS        , PXPX       , SS        ,       PYPY       , PZS       , PXPX       , PZS       ,
      PYPY      , PZPZ      , PXPX       , PZPZ      ,       PYPY      , PXPX      , PXPX       , PYPY      ,       PYPY       , PYPY      , PXPX       , PXPX      ,
      DZ2S      , PYPY      , DZ2S       , PXPX      ,       DZ2PZ     , PYPY      , DZ2PZ      , PXPX      ,       DZ2PY      , PYS       , DZ2PX      , PXS       ,
      DZ2PY     , PYPZ      , DZ2PX      , PXPZ      ,       DZ2DZ2    , PYPY      , DZ2DZ2     , PXPX      ,       DYZS       , PYS       , DXZS       , PXS       ,
      DYZS      , PYPZ      , DXZS       , PXPZ      ,       DYZPZ     , PYS       , DXZPZ      , PXS       ,       DYZPZ      , PYPZ      , DXZPZ      , PXPZ      ,
      DYZPX     , PYPX      , DXZPY      , PYPX      ,       DYZPY     , SS        , DXZPX      , SS        ,       DYZPY      , PZS       , DXZPX      , PZS       ,
      DYZPY     , PZPZ      , DXZPX      , PZPZ      ,       DYZPY     , PXPX      , DXZPX      , PYPY      ,       DYZPY      , PYPY      , DXZPX      , PXPX      ,
      DYZDZ2    , PYS       , DXZDZ2     , PXS       ,       DYZDZ2    , PYPZ      , DXZDZ2     , PXPZ      ,       DYZDYZ     , SS        , DXZDXZ     , SS        ,
      DYZDYZ    , PZS       , DXZDXZ     , PZS       ,       DYZDYZ    , PZPZ      , DXZDXZ     , PZPZ      ,       DYZDYZ     , PXPX      , DXZDXZ     , PYPY      ,
      DYZDYZ    , PYPY      , DXZDXZ     , PXPX      ,       DX2Y2PZ   , PXPX      , DXZPY      , PYPX      ,       DX2Y2DX2Y2 , PYPY      , DX2Y2DX2Y2 , PXPX      ,
      DXYS      , PYPX      , DX2Y2S     , PXPX      ,       DXYPZ     , PYPX      , DXZPY      , PYPX      ,       DXYPX      , PYS       , DX2Y2PX    , PXS       ,
      DXYPX     , PYPZ      , DX2Y2PX    , PXPZ      ,       DXYPY     , PXS       , DX2Y2PX    , PXS       ,       DXYPY      , PXPZ      , DX2Y2PX    , PXPZ      ,
      DXYDZ2    , PYPX      , DX2Y2DZ2   , PXPX      ,       DXYDXZ    , PYS       , DX2Y2DXZ   , PXS       ,       DXYDXZ     , PYPZ      , DX2Y2DXZ   , PXPZ      ,
      DXYDYZ    , PXS       , DX2Y2DXZ   , PXS       ,       DXYDYZ    , PXPZ      , DX2Y2DXZ   , PXPZ      ,       DXYDXY     , SS        , DX2Y2DX2Y2 , SS        ,
      DXYDXY    , PZS       , DX2Y2DX2Y2 , PZS       ,       DXYDXY    , PZPZ      , DX2Y2DX2Y2 , PZPZ      ,       DXYDXY     , PXPX      , DX2Y2DX2Y2 , PXPX      ,
      DXYDXY    , PYPY      , DX2Y2DX2Y2 , PXPX      } ;
/*
const Integer  SPUNIQUE[]   = {
      S         , S         , S         , S         ,       S         , S         , PZ        , S         ,       S         , S         , PZ        , PZ        ,       S         , S         , PX        , PX        ,       PZ        , S         , S         , S         ,
      PZ        , S         , PZ        , S         ,       PZ        , S         , PZ        , PZ        ,       PZ        , S         , PX        , PX        ,       PZ        , PZ        , S         , S         ,       PZ        , PZ        , PZ        , S         ,
      PZ        , PZ        , PZ        , PZ        ,       PZ        , PZ        , PX        , PX        ,       PX        , S         , PX        , S         ,       PX        , S         , PX        , PZ        ,       PX        , PZ        , PX        , S         ,
      PX        , PZ        , PX        , PZ        ,       PX        , PX        , S         , S         ,       PX        , PX        , PZ        , S         ,       PX        , PX        , PZ        , PZ        ,       PX        , PX        , PX        , PX        ,
      PX        , PX        , PY        , PY        ,       PY        , PX        , PY        , PX        ,       DZ2       , S         , S         , S         ,       DZ2       , S         , PZ        , S         ,       DZ2       , S         , PZ        , PZ        ,
      DZ2       , S         , PX        , PX        ,       DZ2       , PZ        , S         , S         ,       DZ2       , PZ        , PZ        , S         ,       DZ2       , PZ        , PZ        , PZ        ,       DZ2       , PZ        , PX        , PX        ,
      DZ2       , PX        , PX        , S         ,       DZ2       , PX        , PX        , PZ        ,       DZ2       , DZ2       , S         , S         ,       DZ2       , DZ2       , PZ        , S         ,       DZ2       , DZ2       , PZ        , PZ        ,
      DZ2       , DZ2       , PX        , PX        ,       DXZ       , S         , PX        , S         ,       DXZ       , S         , PX        , PZ        ,       DXZ       , PZ        , PX        , S         ,       DXZ       , PZ        , PX        , PZ        ,
      DXZ       , PX        , S         , S         ,       DXZ       , PX        , PZ        , S         ,       DXZ       , PX        , PZ        , PZ        ,       DXZ       , PX        , PX        , PX        ,       DXZ       , PX        , PY        , PY        ,
      DXZ       , PY        , PY        , PX        ,       DXZ       , DZ2       , PX        , S         ,       DXZ       , DZ2       , PX        , PZ        ,       DXZ       , DXZ       , S         , S         ,       DXZ       , DXZ       , PZ        , S         ,
      DXZ       , DXZ       , PZ        , PZ        ,       DXZ       , DXZ       , PX        , PX        ,       DXZ       , DXZ       , PY        , PY        ,       DYZ       , DXZ       , PY        , PX        ,       DX2Y2     , S         , PX        , PX        ,
      DX2Y2     , PX        , PX        , S         ,       DX2Y2     , PX        , PX        , PZ        ,       DX2Y2     , DZ2       , PX        , PX        ,       DX2Y2     , DXZ       , PX        , S         ,       DX2Y2     , DXZ       , PX        , PZ        ,
      DX2Y2     , DX2Y2     , S         , S         ,       DX2Y2     , DX2Y2     , PZ        , S         ,       DX2Y2     , DX2Y2     , PZ        , PZ        ,       DX2Y2     , DX2Y2     , PX        , PX } ;
*/
const Integer  SPUNIQUE[168] = {
      DZ2       , S         , S         , S         ,       DZ2       , S         , PZ        , S         ,       DZ2       , S         , PZ        , PZ        ,
      DZ2       , S         , PX        , PX        ,       DZ2       , PZ        , S         , S         ,       DZ2       , PZ        , PZ        , S         ,       DZ2       , PZ        , PZ        , PZ        ,       DZ2       , PZ        , PX        , PX        ,
      DZ2       , PX        , PX        , S         ,       DZ2       , PX        , PX        , PZ        ,       DZ2       , DZ2       , S         , S         ,       DZ2       , DZ2       , PZ        , S         ,       DZ2       , DZ2       , PZ        , PZ        ,
      DZ2       , DZ2       , PX        , PX        ,       DXZ       , S         , PX        , S         ,       DXZ       , S         , PX        , PZ        ,       DXZ       , PZ        , PX        , S         ,       DXZ       , PZ        , PX        , PZ        ,
      DXZ       , PX        , S         , S         ,       DXZ       , PX        , PZ        , S         ,       DXZ       , PX        , PZ        , PZ        ,       DXZ       , PX        , PX        , PX        ,       DXZ       , PX        , PY        , PY        ,
      DXZ       , PY        , PY        , PX        ,       DXZ       , DZ2       , PX        , S         ,       DXZ       , DZ2       , PX        , PZ        ,       DXZ       , DXZ       , S         , S         ,       DXZ       , DXZ       , PZ        , S         ,
      DXZ       , DXZ       , PZ        , PZ        ,       DXZ       , DXZ       , PX        , PX        ,       DXZ       , DXZ       , PY        , PY        ,       DYZ       , DXZ       , PY        , PX        ,       DX2Y2     , S         , PX        , PX        ,
      DX2Y2     , PX        , PX        , S         ,       DX2Y2     , PX        , PX        , PZ        ,       DX2Y2     , DZ2       , PX        , PX        ,       DX2Y2     , DXZ       , PX        , S         ,       DX2Y2     , DXZ       , PX        , PZ        ,
      DX2Y2     , DX2Y2     , S         , S         ,       DX2Y2     , DX2Y2     , PZ        , S         ,       DX2Y2     , DX2Y2     , PZ        , PZ        ,       DX2Y2     , DX2Y2     , PX        , PX } ;

/* . Second atom is spd. */
const Integer  SPDNEGATIVE[200] = {
      PYS       , DX2Y2PY   , PXS       , DX2Y2PX   ,       PYS       , DX2Y2DYZ  , PXS       , DX2Y2DXZ  ,       PYPZ      , DX2Y2PY   , PXPZ      , DX2Y2PX    ,
      PYPZ      , DX2Y2DYZ  , PXPZ      , DX2Y2DXZ  ,       PYPY      , DX2Y2S    , PXPX      , DX2Y2S    ,       PYPY      , DX2Y2PZ   , PXPX      , DX2Y2PZ    ,
      PYPY      , DX2Y2DZ2  , PXPX      , DX2Y2DZ2  ,       DZ2PY     , DX2Y2PY   , DZ2PX     , DX2Y2PX   ,       DZ2PY     , DX2Y2DYZ  , DZ2PX     , DX2Y2DXZ   ,
      DYZS      , DX2Y2PY   , DXZS      , DX2Y2PX   ,       DYZS      , DX2Y2DYZ  , DXZS      , DX2Y2DXZ  ,       DYZPZ     , DX2Y2PY   , DXZPZ     , DX2Y2PX    ,
      DYZPZ     , DX2Y2DYZ  , DXZPZ     , DX2Y2DXZ  ,       DYZPY     , DX2Y2S    , DXZPX     , DX2Y2S    ,       DYZPY     , DX2Y2PZ   , DXZPX     , DX2Y2PZ    ,
      DYZPY     , DX2Y2DZ2  , DXZPX     , DX2Y2DZ2  ,       DYZDZ2    , DX2Y2PY   , DXZDZ2    , DX2Y2PX   ,       DYZDZ2    , DX2Y2DYZ  , DXZDZ2    , DX2Y2DXZ   ,
      DYZDYZ    , DX2Y2S    , DXZDXZ    , DX2Y2S    ,       DYZDYZ    , DX2Y2PZ   , DXZDXZ    , DX2Y2PZ   ,       DYZDYZ    , DX2Y2DZ2  , DXZDXZ    , DX2Y2DZ2   ,
      DX2Y2S    , PYPY      , DX2Y2S    , PXPX      ,       DX2Y2S    , DYZPY     , DX2Y2S    , DXZPX     ,       DX2Y2S    , DYZDYZ    , DX2Y2S    , DXZDXZ     ,
      DX2Y2PZ   , PYPY      , DXZPY     , PYPX      ,       DX2Y2PZ   , DYZPY     , DXZPX     , DX2Y2PZ   ,       DX2Y2PZ   , DYZDYZ    , DXZPY     , DYZDXZ     ,
      DX2Y2PY   , PYS       , DX2Y2PX   , PXS       ,       DX2Y2PY   , PYPZ      , DX2Y2PX   , PXPZ      ,       DX2Y2PY   , DZ2PY     , DX2Y2PX   , DZ2PX      ,
      DX2Y2PY   , DYZS      , DX2Y2PX   , DXZS      ,       DX2Y2PY   , DYZPZ     , DX2Y2PX   , DXZPZ     ,       DX2Y2PY   , DYZDZ2    , DX2Y2PX   , DXZDZ2     ,
      DX2Y2PY   , DXYPX     , DX2Y2PX   , DXYPY     ,       DX2Y2PY   , DXYDXZ    , DX2Y2PX   , DXYDYZ    ,       DX2Y2DZ2  , PYPY      , DX2Y2DZ2  , PXPX       ,
      DX2Y2DZ2  , DYZPY     , DX2Y2DZ2  , DXZPX     ,       DX2Y2DZ2  , DYZDYZ    , DX2Y2DZ2  , DXZDXZ    ,       DX2Y2DYZ  , PYS       , DX2Y2DXZ  , PXS        ,
      DX2Y2DYZ  , PYPZ      , DX2Y2DXZ  , PXPZ      ,       DX2Y2DYZ  , DZ2PY     , DX2Y2DXZ  , DZ2PX     ,       DX2Y2DYZ  , DYZS      , DX2Y2DXZ  , DXZS       ,
      DX2Y2DYZ  , DYZPZ     , DX2Y2DXZ  , DXZPZ     ,       DX2Y2DYZ  , DYZDZ2    , DX2Y2DXZ  , DXZDZ2    ,       DX2Y2DYZ  , DXYPX     , DX2Y2DXZ  , DXYPY      ,
      DX2Y2DYZ  , DXYDXZ    , DX2Y2DXZ  , DXYDYZ    ,       DXYPX     , DX2Y2PY   , DX2Y2PX   , DXYPY     ,       DXYPX     , DX2Y2DYZ  , DX2Y2PX   , DXYDYZ     ,
      DXYDXZ    , DX2Y2PY   , DX2Y2DXZ  , DXYPY     ,       DXYDXZ    , DX2Y2DYZ  , DX2Y2DXZ  , DXYDYZ } ;
const Integer  SPDPOSITIVE[948] = {
      SS        , PYPY      , SS        , PXPX      ,       SS        , DYZPY     , SS        , DXZPX     ,       SS        , DYZDYZ    , SS        , DXZDXZ     ,
      SS        , DXYDXY    , SS        , DX2Y2DX2Y2 ,       PZS       , PYPY      , PZS       , PXPX      ,       PZS       , DYZPY     , PZS       , DXZPX      ,
      PZS       , DYZDYZ    , PZS       , DXZDXZ    ,       PZS       , DXYDXY    , PZS       , DX2Y2DX2Y2 ,       PZPZ      , PYPY      , PZPZ      , PXPX       ,
      PZPZ      , DYZPY     , PZPZ      , DXZPX     ,       PZPZ      , DYZDYZ    , PZPZ      , DXZDXZ    ,       PZPZ      , DXYDXY    , PZPZ      , DX2Y2DX2Y2 ,
      PXS       , DXYPY     , PXS       , DX2Y2PX   ,       PXS       , DXYDYZ    , PXS       , DX2Y2DXZ  ,       PXPZ      , DXYPY     , PXPZ      , DX2Y2PX    ,
      PXPZ      , DXYDYZ    , PXPZ      , DX2Y2DXZ  ,       PXPX      , DXYDXY    , PXPX      , DX2Y2DX2Y2 ,       PYS       , PYS       , PXS       , PXS        ,
      PYS       , PYPZ      , PXS       , PXPZ      ,       PYS       , DZ2PY     , PXS       , DZ2PX     ,       PYS       , DYZS      , PXS       , DXZS       ,
      PYS       , DYZPZ     , PXS       , DXZPZ     ,       PYS       , DYZDZ2    , PXS       , DXZDZ2    ,       PYS       , DXYPX     , PXS       , DX2Y2PX    ,
      PYS       , DXYDXZ    , PXS       , DX2Y2DXZ  ,       PYPZ      , PYS       , PXPZ      , PXS       ,       PYPZ      , PYPZ      , PXPZ      , PXPZ       ,
      PYPZ      , DZ2PY     , PXPZ      , DZ2PX     ,       PYPZ      , DYZS      , PXPZ      , DXZS      ,       PYPZ      , DYZPZ     , PXPZ      , DXZPZ      ,
      PYPZ      , DYZDZ2    , PXPZ      , DXZDZ2    ,       PYPZ      , DXYPX     , PXPZ      , DX2Y2PX   ,       PYPZ      , DXYDXZ    , PXPZ      , DX2Y2DXZ   ,
      PYPX      , DXZPY     , PXPX      , DX2Y2PZ   ,       PYPX      , DYZPX     , PXPX      , DX2Y2PZ   ,       PYPX      , DXYS      , PXPX      , DX2Y2S     ,
      PYPX      , DXYPZ     , PXPX      , DX2Y2PZ   ,       PYPX      , DXYDZ2    , PXPX      , DX2Y2DZ2  ,       PYPY      , SS        , PXPX      , SS         ,
      PYPY      , PZS       , PXPX      , PZS       ,       PYPY      , PZPZ      , PXPX      , PZPZ      ,       PYPY      , PXPX      , PXPX      , PYPY       ,
      PYPY      , PYPY      , PXPX      , PXPX      ,       PYPY      , DZ2S      , PXPX      , DZ2S      ,       PYPY      , DZ2PZ     , PXPX      , DZ2PZ      ,
      PYPY      , DZ2DZ2    , PXPX      , DZ2DZ2    ,       PYPY      , DXZPX     , PXPX      , DYZPY     ,       PYPY      , DXZDXZ    , PXPX      , DYZDYZ     ,
      PYPY      , DYZPY     , PXPX      , DXZPX     ,       PYPY      , DYZDYZ    , PXPX      , DXZDXZ    ,       PYPY      , DX2Y2DX2Y2 , PXPX      , DX2Y2DX2Y2 ,
      PYPY      , DXYDXY    , PXPX      , DX2Y2DX2Y2 ,       DZ2S      , PYPY      , DZ2S      , PXPX      ,       DZ2S      , DYZPY     , DZ2S      , DXZPX      ,
      DZ2S      , DYZDYZ    , DZ2S      , DXZDXZ    ,       DZ2S      , DXYDXY    , DZ2S      , DX2Y2DX2Y2 ,       DZ2PZ     , PYPY      , DZ2PZ     , PXPX       ,
      DZ2PZ     , DYZPY     , DZ2PZ     , DXZPX     ,       DZ2PZ     , DYZDYZ    , DZ2PZ     , DXZDXZ    ,       DZ2PZ     , DXYDXY    , DZ2PZ     , DX2Y2DX2Y2 ,
      DZ2PX     , DXYPY     , DZ2PX     , DX2Y2PX   ,       DZ2PX     , DXYDYZ    , DZ2PX     , DX2Y2DXZ  ,       DZ2PY     , PYS       , DZ2PX     , PXS        ,
      DZ2PY     , PYPZ      , DZ2PX     , PXPZ      ,       DZ2PY     , DZ2PY     , DZ2PX     , DZ2PX     ,       DZ2PY     , DYZS      , DZ2PX     , DXZS       ,
      DZ2PY     , DYZPZ     , DZ2PX     , DXZPZ     ,       DZ2PY     , DYZDZ2    , DZ2PX     , DXZDZ2    ,       DZ2PY     , DXYPX     , DZ2PX     , DX2Y2PX    ,
      DZ2PY     , DXYDXZ    , DZ2PX     , DX2Y2DXZ  ,       DZ2DZ2    , PYPY      , DZ2DZ2    , PXPX      ,       DZ2DZ2    , DYZPY     , DZ2DZ2    , DXZPX      ,
      DZ2DZ2    , DYZDYZ    , DZ2DZ2    , DXZDXZ    ,       DZ2DZ2    , DXYDXY    , DZ2DZ2    , DX2Y2DX2Y2 ,       DXZS      , DXYPY     , DXZS      , DX2Y2PX    ,
      DXZS      , DXYDYZ    , DXZS      , DX2Y2DXZ  ,       DXZPZ     , DXYPY     , DXZPZ     , DX2Y2PX   ,       DXZPZ     , DXYDYZ    , DXZPZ     , DX2Y2DXZ   ,
      DXZPX     , DXYDXY    , DXZPX     , DX2Y2DX2Y2 ,       DXZPY     , DXZPY     , DXZPX     , DX2Y2PZ   ,       DXZPY     , DYZPX     , DXZPX     , DX2Y2PZ    ,
      DXZPY     , DXYS      , DXZPX     , DX2Y2S    ,       DXZPY     , DXYPZ     , DXZPX     , DX2Y2PZ   ,       DXZPY     , DXYDZ2    , DXZPX     , DX2Y2DZ2   ,
      DXZDZ2    , DXYPY     , DXZDZ2    , DX2Y2PX   ,       DXZDZ2    , DXYDYZ    , DXZDZ2    , DX2Y2DXZ  ,       DXZDXZ    , DXYDXY    , DXZDXZ    , DX2Y2DX2Y2 ,
      DYZS      , PYS       , DXZS      , PXS       ,       DYZS      , PYPZ      , DXZS      , PXPZ      ,       DYZS      , DZ2PY     , DXZS      , DZ2PX      ,
      DYZS      , DYZS      , DXZS      , DXZS      ,       DYZS      , DYZPZ     , DXZS      , DXZPZ     ,       DYZS      , DYZDZ2    , DXZS      , DXZDZ2     ,
      DYZS      , DXYPX     , DXZS      , DX2Y2PX   ,       DYZS      , DXYDXZ    , DXZS      , DX2Y2DXZ  ,       DYZPZ     , PYS       , DXZPZ     , PXS        ,
      DYZPZ     , PYPZ      , DXZPZ     , PXPZ      ,       DYZPZ     , DZ2PY     , DXZPZ     , DZ2PX     ,       DYZPZ     , DYZS      , DXZPZ     , DXZS       ,
      DYZPZ     , DYZPZ     , DXZPZ     , DXZPZ     ,       DYZPZ     , DYZDZ2    , DXZPZ     , DXZDZ2    ,       DYZPZ     , DXYPX     , DXZPZ     , DX2Y2PX    ,
      DYZPZ     , DXYDXZ    , DXZPZ     , DX2Y2DXZ  ,       DYZPX     , PYPX      , DXZPY     , PYPX      ,       DYZPX     , DXZPY     , DXZPX     , DX2Y2PZ    ,
      DYZPX     , DYZPX     , DXZPX     , DX2Y2PZ   ,       DYZPX     , DYZDXZ    , DXZPY     , DYZDXZ    ,       DYZPX     , DXYS      , DXZPX     , DX2Y2S     ,
      DYZPX     , DXYPZ     , DXZPX     , DX2Y2PZ   ,       DYZPX     , DXYDZ2    , DXZPX     , DX2Y2DZ2  ,       DYZPY     , SS        , DXZPX     , SS         ,
      DYZPY     , PZS       , DXZPX     , PZS       ,       DYZPY     , PZPZ      , DXZPX     , PZPZ      ,       DYZPY     , PXPX      , DXZPX     , PYPY       ,
      DYZPY     , PYPY      , DXZPX     , PXPX      ,       DYZPY     , DZ2S      , DXZPX     , DZ2S      ,       DYZPY     , DZ2PZ     , DXZPX     , DZ2PZ      ,
      DYZPY     , DZ2DZ2    , DXZPX     , DZ2DZ2    ,       DYZPY     , DXZPX     , DXZPX     , DYZPY     ,       DYZPY     , DXZDXZ    , DXZPX     , DYZDYZ     ,
      DYZPY     , DYZPY     , DXZPX     , DXZPX     ,       DYZPY     , DYZDYZ    , DXZPX     , DXZDXZ    ,       DYZPY     , DX2Y2DX2Y2 , DXZPX     , DX2Y2DX2Y2 ,
      DYZPY     , DXYDXY    , DXZPX     , DX2Y2DX2Y2 ,       DYZDZ2    , PYS       , DXZDZ2    , PXS       ,       DYZDZ2    , PYPZ      , DXZDZ2    , PXPZ       ,
      DYZDZ2    , DZ2PY     , DXZDZ2    , DZ2PX     ,       DYZDZ2    , DYZS      , DXZDZ2    , DXZS      ,       DYZDZ2    , DYZPZ     , DXZDZ2    , DXZPZ      ,
      DYZDZ2    , DYZDZ2    , DXZDZ2    , DXZDZ2    ,       DYZDZ2    , DXYPX     , DXZDZ2    , DX2Y2PX   ,       DYZDZ2    , DXYDXZ    , DXZDZ2    , DX2Y2DXZ   ,
      DYZDXZ    , DXZPY     , DXZDXZ    , DX2Y2PZ   ,       DYZDXZ    , DYZPX     , DXZDXZ    , DX2Y2PZ   ,       DYZDXZ    , DXYS      , DXZDXZ    , DX2Y2S     ,
      DYZDXZ    , DXYPZ     , DXZDXZ    , DX2Y2PZ   ,       DYZDXZ    , DXYDZ2    , DXZDXZ    , DX2Y2DZ2  ,       DYZDYZ    , SS        , DXZDXZ    , SS         ,
      DYZDYZ    , PZS       , DXZDXZ    , PZS       ,       DYZDYZ    , PZPZ      , DXZDXZ    , PZPZ      ,       DYZDYZ    , PXPX      , DXZDXZ    , PYPY       ,
      DYZDYZ    , PYPY      , DXZDXZ    , PXPX      ,       DYZDYZ    , DZ2S      , DXZDXZ    , DZ2S      ,       DYZDYZ    , DZ2PZ     , DXZDXZ    , DZ2PZ      ,
      DYZDYZ    , DZ2DZ2    , DXZDXZ    , DZ2DZ2    ,       DYZDYZ    , DXZPX     , DXZDXZ    , DYZPY     ,       DYZDYZ    , DXZDXZ    , DXZDXZ    , DYZDYZ     ,
      DYZDYZ    , DYZPY     , DXZDXZ    , DXZPX     ,       DYZDYZ    , DYZDYZ    , DXZDXZ    , DXZDXZ    ,       DYZDYZ    , DX2Y2DX2Y2 , DXZDXZ    , DX2Y2DX2Y2 ,
      DYZDYZ    , DXYDXY    , DXZDXZ    , DX2Y2DX2Y2 ,       DX2Y2S    , DX2Y2PZ   , DX2Y2S    , DXZPX     ,       DX2Y2PZ   , PXPX      , DXZPY     , PYPX       ,
      DX2Y2PZ   , DXZPX     , DXZPX     , DX2Y2PZ   ,       DX2Y2PZ   , DXZDXZ    , DXZPY     , DYZDXZ    ,       DX2Y2PZ   , DX2Y2S    , DXZPX     , DX2Y2S     ,
      DX2Y2PZ   , DX2Y2PZ   , DXZPX     , DX2Y2PZ   ,       DX2Y2PZ   , DX2Y2DZ2  , DXZPX     , DX2Y2DZ2  ,       DX2Y2PY   , DX2Y2PY   , DX2Y2PX   , DX2Y2PX    ,
      DX2Y2PY   , DX2Y2DYZ  , DX2Y2PX   , DX2Y2DXZ  ,       DX2Y2DZ2  , DX2Y2PZ   , DX2Y2DZ2  , DXZPX     ,       DX2Y2DYZ  , DX2Y2PY   , DX2Y2DXZ  , DX2Y2PX    ,
      DX2Y2DYZ  , DX2Y2DYZ  , DX2Y2DXZ  , DX2Y2DXZ  ,       DX2Y2DX2Y2 , PYPY      , DX2Y2DX2Y2 , PXPX      ,       DX2Y2DX2Y2 , DYZPY     , DX2Y2DX2Y2 , DXZPX      ,
      DX2Y2DX2Y2 , DYZDYZ    , DX2Y2DX2Y2 , DXZDXZ    ,       DXYS      , PYPX      , DX2Y2S    , PXPX      ,       DXYS      , DXZPY     , DX2Y2S    , DXZPX      ,
      DXYS      , DYZPX     , DX2Y2S    , DXZPX     ,       DXYS      , DYZDXZ    , DX2Y2S    , DXZDXZ    ,       DXYS      , DXYS      , DX2Y2S    , DX2Y2S     ,
      DXYS      , DXYPZ     , DX2Y2S    , DXZPX     ,       DXYS      , DXYDZ2    , DX2Y2S    , DX2Y2DZ2  ,       DXYPZ     , PYPX      , DXZPY     , PYPX       ,
      DXYPZ     , DXZPY     , DXZPX     , DX2Y2PZ   ,       DXYPZ     , DYZPX     , DXZPX     , DX2Y2PZ   ,       DXYPZ     , DYZDXZ    , DXZPY     , DYZDXZ     ,
      DXYPZ     , DXYS      , DXZPX     , DX2Y2S    ,       DXYPZ     , DXYPZ     , DXZPX     , DX2Y2PZ   ,       DXYPZ     , DXYDZ2    , DXZPX     , DX2Y2DZ2   ,
      DXYPX     , PYS       , DX2Y2PX   , PXS       ,       DXYPX     , PYPZ      , DX2Y2PX   , PXPZ      ,       DXYPX     , DZ2PY     , DX2Y2PX   , DZ2PX      ,
      DXYPX     , DYZS      , DX2Y2PX   , DXZS      ,       DXYPX     , DYZPZ     , DX2Y2PX   , DXZPZ     ,       DXYPX     , DYZDZ2    , DX2Y2PX   , DXZDZ2     ,
      DXYPX     , DXYPX     , DX2Y2PX   , DX2Y2PX   ,       DXYPX     , DXYDXZ    , DX2Y2PX   , DX2Y2DXZ  ,       DXYPY     , PXS       , DX2Y2PX   , PXS        ,
      DXYPY     , PXPZ      , DX2Y2PX   , PXPZ      ,       DXYPY     , DZ2PX     , DX2Y2PX   , DZ2PX     ,       DXYPY     , DXZS      , DX2Y2PX   , DXZS       ,
      DXYPY     , DXZPZ     , DX2Y2PX   , DXZPZ     ,       DXYPY     , DXZDZ2    , DX2Y2PX   , DXZDZ2    ,       DXYPY     , DX2Y2PX   , DX2Y2PX   , DXYPY      ,
      DXYPY     , DX2Y2DXZ  , DX2Y2PX   , DXYDYZ    ,       DXYPY     , DXYPY     , DX2Y2PX   , DX2Y2PX   ,       DXYPY     , DXYDYZ    , DX2Y2PX   , DX2Y2DXZ   ,
      DXYDZ2    , PYPX      , DX2Y2DZ2  , PXPX      ,       DXYDZ2    , DXZPY     , DX2Y2DZ2  , DXZPX     ,       DXYDZ2    , DYZPX     , DX2Y2DZ2  , DXZPX      ,
      DXYDZ2    , DYZDXZ    , DX2Y2DZ2  , DXZDXZ    ,       DXYDZ2    , DXYS      , DX2Y2DZ2  , DX2Y2S    ,       DXYDZ2    , DXYPZ     , DX2Y2DZ2  , DXZPX      ,
      DXYDZ2    , DXYDZ2    , DX2Y2DZ2  , DX2Y2DZ2  ,       DXYDXZ    , PYS       , DX2Y2DXZ  , PXS       ,       DXYDXZ    , PYPZ      , DX2Y2DXZ  , PXPZ       ,
      DXYDXZ    , DZ2PY     , DX2Y2DXZ  , DZ2PX     ,       DXYDXZ    , DYZS      , DX2Y2DXZ  , DXZS      ,       DXYDXZ    , DYZPZ     , DX2Y2DXZ  , DXZPZ      ,
      DXYDXZ    , DYZDZ2    , DX2Y2DXZ  , DXZDZ2    ,       DXYDXZ    , DXYPX     , DX2Y2DXZ  , DX2Y2PX   ,       DXYDXZ    , DXYDXZ    , DX2Y2DXZ  , DX2Y2DXZ   ,
      DXYDYZ    , PXS       , DX2Y2DXZ  , PXS       ,       DXYDYZ    , PXPZ      , DX2Y2DXZ  , PXPZ      ,       DXYDYZ    , DZ2PX     , DX2Y2DXZ  , DZ2PX      ,
      DXYDYZ    , DXZS      , DX2Y2DXZ  , DXZS      ,       DXYDYZ    , DXZPZ     , DX2Y2DXZ  , DXZPZ     ,       DXYDYZ    , DXZDZ2    , DX2Y2DXZ  , DXZDZ2     ,
      DXYDYZ    , DX2Y2PX   , DX2Y2DXZ  , DXYPY     ,       DXYDYZ    , DX2Y2DXZ  , DX2Y2DXZ  , DXYDYZ    ,       DXYDYZ    , DXYPY     , DX2Y2DXZ  , DX2Y2PX    ,
      DXYDYZ    , DXYDYZ    , DX2Y2DXZ  , DX2Y2DXZ  ,       DXYDXY    , SS        , DX2Y2DX2Y2 , SS        ,       DXYDXY    , PZS       , DX2Y2DX2Y2 , PZS        ,
      DXYDXY    , PZPZ      , DX2Y2DX2Y2 , PZPZ      ,       DXYDXY    , PXPX      , DX2Y2DX2Y2 , PXPX      ,       DXYDXY    , PYPY      , DX2Y2DX2Y2 , PXPX       ,
      DXYDXY    , DZ2S      , DX2Y2DX2Y2 , DZ2S      ,       DXYDXY    , DZ2PZ     , DX2Y2DX2Y2 , DZ2PZ     ,       DXYDXY    , DZ2DZ2    , DX2Y2DX2Y2 , DZ2DZ2     ,
      DXYDXY    , DXZPX     , DX2Y2DX2Y2 , DXZPX     ,       DXYDXY    , DXZDXZ    , DX2Y2DX2Y2 , DXZDXZ    ,       DXYDXY    , DYZPY     , DX2Y2DX2Y2 , DXZPX      ,
      DXYDXY    , DYZDYZ    , DX2Y2DX2Y2 , DXZDXZ    ,       DXYDXY    , DX2Y2DX2Y2 , DX2Y2DX2Y2 , DXYDXY    ,       DXYDXY    , DXYDXY    , DX2Y2DX2Y2 , DX2Y2DX2Y2 } ;
/*
const Integer  SPDUNIQUE[]   = {
      S         , S         , S         , S         ,       S         , S         , PZ        , S         ,       S         , S         , PZ        , PZ        ,       S         , S         , PX        , PX        ,       S         , S         , DZ2       , S         ,
      S         , S         , DZ2       , PZ        ,       S         , S         , DZ2       , DZ2       ,       S         , S         , DXZ       , PX        ,       S         , S         , DXZ       , DXZ       ,       S         , S         , DX2Y2     , DX2Y2     ,
      PZ        , S         , S         , S         ,       PZ        , S         , PZ        , S         ,       PZ        , S         , PZ        , PZ        ,       PZ        , S         , PX        , PX        ,       PZ        , S         , DZ2       , S         ,
      PZ        , S         , DZ2       , PZ        ,       PZ        , S         , DZ2       , DZ2       ,       PZ        , S         , DXZ       , PX        ,       PZ        , S         , DXZ       , DXZ       ,       PZ        , S         , DX2Y2     , DX2Y2     ,
      PZ        , PZ        , S         , S         ,       PZ        , PZ        , PZ        , S         ,       PZ        , PZ        , PZ        , PZ        ,       PZ        , PZ        , PX        , PX        ,       PZ        , PZ        , DZ2       , S         ,
      PZ        , PZ        , DZ2       , PZ        ,       PZ        , PZ        , DZ2       , DZ2       ,       PZ        , PZ        , DXZ       , PX        ,       PZ        , PZ        , DXZ       , DXZ       ,       PZ        , PZ        , DX2Y2     , DX2Y2     ,
      PX        , S         , PX        , S         ,       PX        , S         , PX        , PZ        ,       PX        , S         , DZ2       , PX        ,       PX        , S         , DXZ       , S         ,       PX        , S         , DXZ       , PZ        ,
      PX        , S         , DXZ       , DZ2       ,       PX        , S         , DX2Y2     , PX        ,       PX        , S         , DX2Y2     , DXZ       ,       PX        , PZ        , PX        , S         ,       PX        , PZ        , PX        , PZ        ,
      PX        , PZ        , DZ2       , PX        ,       PX        , PZ        , DXZ       , S         ,       PX        , PZ        , DXZ       , PZ        ,       PX        , PZ        , DXZ       , DZ2       ,       PX        , PZ        , DX2Y2     , PX        ,
      PX        , PZ        , DX2Y2     , DXZ       ,       PX        , PX        , S         , S         ,       PX        , PX        , PZ        , S         ,       PX        , PX        , PZ        , PZ        ,       PX        , PX        , PX        , PX        ,
      PX        , PX        , PY        , PY        ,       PX        , PX        , DZ2       , S         ,       PX        , PX        , DZ2       , PZ        ,       PX        , PX        , DZ2       , DZ2       ,       PX        , PX        , DXZ       , PX        ,
      PX        , PX        , DXZ       , DXZ       ,       PX        , PX        , DYZ       , PY        ,       PX        , PX        , DYZ       , DYZ       ,       PX        , PX        , DX2Y2     , S         ,       PX        , PX        , DX2Y2     , PZ        ,
      PX        , PX        , DX2Y2     , DZ2       ,       PX        , PX        , DX2Y2     , DX2Y2     ,       PY        , PX        , PY        , PX        ,       PY        , PX        , DYZ       , DXZ       ,       DZ2       , S         , S         , S         ,
      DZ2       , S         , PZ        , S         ,       DZ2       , S         , PZ        , PZ        ,       DZ2       , S         , PX        , PX        ,       DZ2       , S         , DZ2       , S         ,       DZ2       , S         , DZ2       , PZ        ,
      DZ2       , S         , DZ2       , DZ2       ,       DZ2       , S         , DXZ       , PX        ,       DZ2       , S         , DXZ       , DXZ       ,       DZ2       , S         , DX2Y2     , DX2Y2     ,       DZ2       , PZ        , S         , S         ,
      DZ2       , PZ        , PZ        , S         ,       DZ2       , PZ        , PZ        , PZ        ,       DZ2       , PZ        , PX        , PX        ,       DZ2       , PZ        , DZ2       , S         ,       DZ2       , PZ        , DZ2       , PZ        ,
      DZ2       , PZ        , DZ2       , DZ2       ,       DZ2       , PZ        , DXZ       , PX        ,       DZ2       , PZ        , DXZ       , DXZ       ,       DZ2       , PZ        , DX2Y2     , DX2Y2     ,       DZ2       , PX        , PX        , S         ,
      DZ2       , PX        , PX        , PZ        ,       DZ2       , PX        , DZ2       , PX        ,       DZ2       , PX        , DXZ       , S         ,       DZ2       , PX        , DXZ       , PZ        ,       DZ2       , PX        , DXZ       , DZ2       ,
      DZ2       , PX        , DX2Y2     , PX        ,       DZ2       , PX        , DX2Y2     , DXZ       ,       DZ2       , DZ2       , S         , S         ,       DZ2       , DZ2       , PZ        , S         ,       DZ2       , DZ2       , PZ        , PZ        ,
      DZ2       , DZ2       , PX        , PX        ,       DZ2       , DZ2       , DZ2       , S         ,       DZ2       , DZ2       , DZ2       , PZ        ,       DZ2       , DZ2       , DZ2       , DZ2       ,       DZ2       , DZ2       , DXZ       , PX        ,
      DZ2       , DZ2       , DXZ       , DXZ       ,       DZ2       , DZ2       , DX2Y2     , DX2Y2     ,       DXZ       , S         , PX        , S         ,       DXZ       , S         , PX        , PZ        ,       DXZ       , S         , DZ2       , PX        ,
      DXZ       , S         , DXZ       , S         ,       DXZ       , S         , DXZ       , PZ        ,       DXZ       , S         , DXZ       , DZ2       ,       DXZ       , S         , DX2Y2     , PX        ,       DXZ       , S         , DX2Y2     , DXZ       ,
      DXZ       , PZ        , PX        , S         ,       DXZ       , PZ        , PX        , PZ        ,       DXZ       , PZ        , DZ2       , PX        ,       DXZ       , PZ        , DXZ       , S         ,       DXZ       , PZ        , DXZ       , PZ        ,
      DXZ       , PZ        , DXZ       , DZ2       ,       DXZ       , PZ        , DX2Y2     , PX        ,       DXZ       , PZ        , DX2Y2     , DXZ       ,       DXZ       , PX        , S         , S         ,       DXZ       , PX        , PZ        , S         ,
      DXZ       , PX        , PZ        , PZ        ,       DXZ       , PX        , PX        , PX        ,       DXZ       , PX        , PY        , PY        ,       DXZ       , PX        , DZ2       , S         ,       DXZ       , PX        , DZ2       , PZ        ,
      DXZ       , PX        , DZ2       , DZ2       ,       DXZ       , PX        , DXZ       , PX        ,       DXZ       , PX        , DXZ       , DXZ       ,       DXZ       , PX        , DYZ       , PY        ,       DXZ       , PX        , DYZ       , DYZ       ,
      DXZ       , PX        , DX2Y2     , S         ,       DXZ       , PX        , DX2Y2     , PZ        ,       DXZ       , PX        , DX2Y2     , DZ2       ,       DXZ       , PX        , DX2Y2     , DX2Y2     ,       DXZ       , PY        , PY        , PX        ,
      DXZ       , PY        , DYZ       , DXZ       ,       DXZ       , DZ2       , PX        , S         ,       DXZ       , DZ2       , PX        , PZ        ,       DXZ       , DZ2       , DZ2       , PX        ,       DXZ       , DZ2       , DXZ       , S         ,
      DXZ       , DZ2       , DXZ       , PZ        ,       DXZ       , DZ2       , DXZ       , DZ2       ,       DXZ       , DZ2       , DX2Y2     , PX        ,       DXZ       , DZ2       , DX2Y2     , DXZ       ,       DXZ       , DXZ       , S         , S         ,
      DXZ       , DXZ       , PZ        , S         ,       DXZ       , DXZ       , PZ        , PZ        ,       DXZ       , DXZ       , PX        , PX        ,       DXZ       , DXZ       , PY        , PY        ,       DXZ       , DXZ       , DZ2       , S         ,
      DXZ       , DXZ       , DZ2       , PZ        ,       DXZ       , DXZ       , DZ2       , DZ2       ,       DXZ       , DXZ       , DXZ       , PX        ,       DXZ       , DXZ       , DXZ       , DXZ       ,       DXZ       , DXZ       , DYZ       , PY        ,
      DXZ       , DXZ       , DYZ       , DYZ       ,       DXZ       , DXZ       , DX2Y2     , S         ,       DXZ       , DXZ       , DX2Y2     , PZ        ,       DXZ       , DXZ       , DX2Y2     , DZ2       ,       DXZ       , DXZ       , DX2Y2     , DX2Y2     ,
      DYZ       , DXZ       , PY        , PX        ,       DYZ       , DXZ       , DYZ       , DXZ       ,       DX2Y2     , S         , PX        , PX        ,       DX2Y2     , S         , DXZ       , PX        ,       DX2Y2     , S         , DXZ       , DXZ       ,
      DX2Y2     , S         , DX2Y2     , S         ,       DX2Y2     , S         , DX2Y2     , DZ2       ,       DX2Y2     , PX        , PX        , S         ,       DX2Y2     , PX        , PX        , PZ        ,       DX2Y2     , PX        , DZ2       , PX        ,
      DX2Y2     , PX        , DXZ       , S         ,       DX2Y2     , PX        , DXZ       , PZ        ,       DX2Y2     , PX        , DXZ       , DZ2       ,       DX2Y2     , PX        , DX2Y2     , PX        ,       DX2Y2     , PX        , DX2Y2     , DXZ       ,
      DX2Y2     , PX        , DXY       , PY        ,       DX2Y2     , PX        , DXY       , DYZ       ,       DX2Y2     , DZ2       , PX        , PX        ,       DX2Y2     , DZ2       , DXZ       , PX        ,       DX2Y2     , DZ2       , DXZ       , DXZ       ,
      DX2Y2     , DZ2       , DX2Y2     , S         ,       DX2Y2     , DZ2       , DX2Y2     , DZ2       ,       DX2Y2     , DXZ       , PX        , S         ,       DX2Y2     , DXZ       , PX        , PZ        ,       DX2Y2     , DXZ       , DZ2       , PX        ,
      DX2Y2     , DXZ       , DXZ       , S         ,       DX2Y2     , DXZ       , DXZ       , PZ        ,       DX2Y2     , DXZ       , DXZ       , DZ2       ,       DX2Y2     , DXZ       , DX2Y2     , PX        ,       DX2Y2     , DXZ       , DX2Y2     , DXZ       ,
      DX2Y2     , DXZ       , DXY       , PY        ,       DX2Y2     , DXZ       , DXY       , DYZ       ,       DX2Y2     , DX2Y2     , S         , S         ,       DX2Y2     , DX2Y2     , PZ        , S         ,       DX2Y2     , DX2Y2     , PZ        , PZ        ,
      DX2Y2     , DX2Y2     , PX        , PX        ,       DX2Y2     , DX2Y2     , DZ2       , S         ,       DX2Y2     , DX2Y2     , DZ2       , PZ        ,       DX2Y2     , DX2Y2     , DZ2       , DZ2       ,       DX2Y2     , DX2Y2     , DXZ       , PX        ,
      DX2Y2     , DX2Y2     , DXZ       , DXZ       ,       DX2Y2     , DX2Y2     , DX2Y2     , DX2Y2     ,       DX2Y2     , DX2Y2     , DXY       , DXY       ,       DXY       , DX2Y2     , DXY       , DX2Y2 } ;
*/
const Integer  SPDUNIQUE[728] = {
      S         , S         , DZ2       , S         ,
      S         , S         , DZ2       , PZ        ,       S         , S         , DZ2       , DZ2       ,       S         , S         , DXZ       , PX        ,       S         , S         , DXZ       , DXZ       ,       S         , S         , DX2Y2     , DX2Y2     ,
      PZ        , S         , DZ2       , S         ,
      PZ        , S         , DZ2       , PZ        ,       PZ        , S         , DZ2       , DZ2       ,       PZ        , S         , DXZ       , PX        ,       PZ        , S         , DXZ       , DXZ       ,       PZ        , S         , DX2Y2     , DX2Y2     ,
      PZ        , PZ        , DZ2       , S         ,
      PZ        , PZ        , DZ2       , PZ        ,       PZ        , PZ        , DZ2       , DZ2       ,       PZ        , PZ        , DXZ       , PX        ,       PZ        , PZ        , DXZ       , DXZ       ,       PZ        , PZ        , DX2Y2     , DX2Y2     ,
      PX        , S         , DZ2       , PX        ,       PX        , S         , DXZ       , S         ,       PX        , S         , DXZ       , PZ        ,
      PX        , S         , DXZ       , DZ2       ,       PX        , S         , DX2Y2     , PX        ,       PX        , S         , DX2Y2     , DXZ       ,
      PX        , PZ        , DZ2       , PX        ,       PX        , PZ        , DXZ       , S         ,       PX        , PZ        , DXZ       , PZ        ,       PX        , PZ        , DXZ       , DZ2       ,       PX        , PZ        , DX2Y2     , PX        ,
      PX        , PZ        , DX2Y2     , DXZ       ,
      PX        , PX        , DZ2       , S         ,       PX        , PX        , DZ2       , PZ        ,       PX        , PX        , DZ2       , DZ2       ,       PX        , PX        , DXZ       , PX        ,
      PX        , PX        , DXZ       , DXZ       ,       PX        , PX        , DYZ       , PY        ,       PX        , PX        , DYZ       , DYZ       ,       PX        , PX        , DX2Y2     , S         ,       PX        , PX        , DX2Y2     , PZ        ,
      PX        , PX        , DX2Y2     , DZ2       ,       PX        , PX        , DX2Y2     , DX2Y2     ,       PY        , PX        , DYZ       , DXZ       ,       DZ2       , S         , S         , S         ,
      DZ2       , S         , PZ        , S         ,       DZ2       , S         , PZ        , PZ        ,       DZ2       , S         , PX        , PX        ,       DZ2       , S         , DZ2       , S         ,       DZ2       , S         , DZ2       , PZ        ,
      DZ2       , S         , DZ2       , DZ2       ,       DZ2       , S         , DXZ       , PX        ,       DZ2       , S         , DXZ       , DXZ       ,       DZ2       , S         , DX2Y2     , DX2Y2     ,       DZ2       , PZ        , S         , S         ,
      DZ2       , PZ        , PZ        , S         ,       DZ2       , PZ        , PZ        , PZ        ,       DZ2       , PZ        , PX        , PX        ,       DZ2       , PZ        , DZ2       , S         ,       DZ2       , PZ        , DZ2       , PZ        ,
      DZ2       , PZ        , DZ2       , DZ2       ,       DZ2       , PZ        , DXZ       , PX        ,       DZ2       , PZ        , DXZ       , DXZ       ,       DZ2       , PZ        , DX2Y2     , DX2Y2     ,       DZ2       , PX        , PX        , S         ,
      DZ2       , PX        , PX        , PZ        ,       DZ2       , PX        , DZ2       , PX        ,       DZ2       , PX        , DXZ       , S         ,       DZ2       , PX        , DXZ       , PZ        ,       DZ2       , PX        , DXZ       , DZ2       ,
      DZ2       , PX        , DX2Y2     , PX        ,       DZ2       , PX        , DX2Y2     , DXZ       ,       DZ2       , DZ2       , S         , S         ,       DZ2       , DZ2       , PZ        , S         ,       DZ2       , DZ2       , PZ        , PZ        ,
      DZ2       , DZ2       , PX        , PX        ,       DZ2       , DZ2       , DZ2       , S         ,       DZ2       , DZ2       , DZ2       , PZ        ,       DZ2       , DZ2       , DZ2       , DZ2       ,       DZ2       , DZ2       , DXZ       , PX        ,
      DZ2       , DZ2       , DXZ       , DXZ       ,       DZ2       , DZ2       , DX2Y2     , DX2Y2     ,       DXZ       , S         , PX        , S         ,       DXZ       , S         , PX        , PZ        ,       DXZ       , S         , DZ2       , PX        ,
      DXZ       , S         , DXZ       , S         ,       DXZ       , S         , DXZ       , PZ        ,       DXZ       , S         , DXZ       , DZ2       ,       DXZ       , S         , DX2Y2     , PX        ,       DXZ       , S         , DX2Y2     , DXZ       ,
      DXZ       , PZ        , PX        , S         ,       DXZ       , PZ        , PX        , PZ        ,       DXZ       , PZ        , DZ2       , PX        ,       DXZ       , PZ        , DXZ       , S         ,       DXZ       , PZ        , DXZ       , PZ        ,
      DXZ       , PZ        , DXZ       , DZ2       ,       DXZ       , PZ        , DX2Y2     , PX        ,       DXZ       , PZ        , DX2Y2     , DXZ       ,       DXZ       , PX        , S         , S         ,       DXZ       , PX        , PZ        , S         ,
      DXZ       , PX        , PZ        , PZ        ,       DXZ       , PX        , PX        , PX        ,       DXZ       , PX        , PY        , PY        ,       DXZ       , PX        , DZ2       , S         ,       DXZ       , PX        , DZ2       , PZ        ,
      DXZ       , PX        , DZ2       , DZ2       ,       DXZ       , PX        , DXZ       , PX        ,       DXZ       , PX        , DXZ       , DXZ       ,       DXZ       , PX        , DYZ       , PY        ,       DXZ       , PX        , DYZ       , DYZ       ,
      DXZ       , PX        , DX2Y2     , S         ,       DXZ       , PX        , DX2Y2     , PZ        ,       DXZ       , PX        , DX2Y2     , DZ2       ,       DXZ       , PX        , DX2Y2     , DX2Y2     ,       DXZ       , PY        , PY        , PX        ,
      DXZ       , PY        , DYZ       , DXZ       ,       DXZ       , DZ2       , PX        , S         ,       DXZ       , DZ2       , PX        , PZ        ,       DXZ       , DZ2       , DZ2       , PX        ,       DXZ       , DZ2       , DXZ       , S         ,
      DXZ       , DZ2       , DXZ       , PZ        ,       DXZ       , DZ2       , DXZ       , DZ2       ,       DXZ       , DZ2       , DX2Y2     , PX        ,       DXZ       , DZ2       , DX2Y2     , DXZ       ,       DXZ       , DXZ       , S         , S         ,
      DXZ       , DXZ       , PZ        , S         ,       DXZ       , DXZ       , PZ        , PZ        ,       DXZ       , DXZ       , PX        , PX        ,       DXZ       , DXZ       , PY        , PY        ,       DXZ       , DXZ       , DZ2       , S         ,
      DXZ       , DXZ       , DZ2       , PZ        ,       DXZ       , DXZ       , DZ2       , DZ2       ,       DXZ       , DXZ       , DXZ       , PX        ,       DXZ       , DXZ       , DXZ       , DXZ       ,       DXZ       , DXZ       , DYZ       , PY        ,
      DXZ       , DXZ       , DYZ       , DYZ       ,       DXZ       , DXZ       , DX2Y2     , S         ,       DXZ       , DXZ       , DX2Y2     , PZ        ,       DXZ       , DXZ       , DX2Y2     , DZ2       ,       DXZ       , DXZ       , DX2Y2     , DX2Y2     ,
      DYZ       , DXZ       , PY        , PX        ,       DYZ       , DXZ       , DYZ       , DXZ       ,       DX2Y2     , S         , PX        , PX        ,       DX2Y2     , S         , DXZ       , PX        ,       DX2Y2     , S         , DXZ       , DXZ       ,
      DX2Y2     , S         , DX2Y2     , S         ,       DX2Y2     , S         , DX2Y2     , DZ2       ,       DX2Y2     , PX        , PX        , S         ,       DX2Y2     , PX        , PX        , PZ        ,       DX2Y2     , PX        , DZ2       , PX        ,
      DX2Y2     , PX        , DXZ       , S         ,       DX2Y2     , PX        , DXZ       , PZ        ,       DX2Y2     , PX        , DXZ       , DZ2       ,       DX2Y2     , PX        , DX2Y2     , PX        ,       DX2Y2     , PX        , DX2Y2     , DXZ       ,
      DX2Y2     , PX        , DXY       , PY        ,       DX2Y2     , PX        , DXY       , DYZ       ,       DX2Y2     , DZ2       , PX        , PX        ,       DX2Y2     , DZ2       , DXZ       , PX        ,       DX2Y2     , DZ2       , DXZ       , DXZ       ,
      DX2Y2     , DZ2       , DX2Y2     , S         ,       DX2Y2     , DZ2       , DX2Y2     , DZ2       ,       DX2Y2     , DXZ       , PX        , S         ,       DX2Y2     , DXZ       , PX        , PZ        ,       DX2Y2     , DXZ       , DZ2       , PX        ,
      DX2Y2     , DXZ       , DXZ       , S         ,       DX2Y2     , DXZ       , DXZ       , PZ        ,       DX2Y2     , DXZ       , DXZ       , DZ2       ,       DX2Y2     , DXZ       , DX2Y2     , PX        ,       DX2Y2     , DXZ       , DX2Y2     , DXZ       ,
      DX2Y2     , DXZ       , DXY       , PY        ,       DX2Y2     , DXZ       , DXY       , DYZ       ,       DX2Y2     , DX2Y2     , S         , S         ,       DX2Y2     , DX2Y2     , PZ        , S         ,       DX2Y2     , DX2Y2     , PZ        , PZ        ,
      DX2Y2     , DX2Y2     , PX        , PX        ,       DX2Y2     , DX2Y2     , DZ2       , S         ,       DX2Y2     , DX2Y2     , DZ2       , PZ        ,       DX2Y2     , DX2Y2     , DZ2       , DZ2       ,       DX2Y2     , DX2Y2     , DXZ       , PX        ,
      DX2Y2     , DX2Y2     , DXZ       , DXZ       ,       DX2Y2     , DX2Y2     , DX2Y2     , DX2Y2     ,       DX2Y2     , DX2Y2     , DXY       , DXY       ,       DXY       , DX2Y2     , DXY       , DX2Y2 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coefficients of charge interaction terms in the local frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Number of terms per charge distribution. */
const Integer NCHTERMS[45] = { 1, 1, 2, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 1, 2, 1, 1, 1, 0, 1, 3, 1, 1, 0, 1, 1, 1, 3, 1, 0, 1, 1, 1, 1, 1, 2, 1, 0, 1, 1, 1, 1, 1, 0, 2 } ;

/* . Indices of non-zero CH terms + 1. */
const Integer CHINDICES[675] = {
  0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,      0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   5,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,   0,      0,   0,   7,   0,   0,   0,   0,   0,   0,   0,   0,   0,   8,   0,   9,
  0,   0,   0,   0,   0,   0,  10,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  11,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  12,   0,   0,   0,   0,
  0,   0,  13,   0,   0,   0,   0,   0,   0,   0,   0,   0,  14,   0,  15,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,   0,      0,   0,   0,   0,   0,   0,   0,  17,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  18,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,  19,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,  20,   0,   0,   0,   0,   0,   0,   0,   0,   0,  21,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  22,   0,      0,   0,   0,   0,   0,   0,   0,   0,  23,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,  24,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  25,   0,      0,   0,  26,   0,   0,   0,   0,   0,   0,   0,   0,   0,  27,   0,  28,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  29,   0,   0,   0,      0,   0,   0,   0,   0,   0,  30,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,  31,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  32,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  33,   0,   0,   0,   0,
  0,   0,  34,   0,   0,   0,   0,   0,   0,   0,   0,   0,  35,   0,  36,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  37,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  38,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,  39,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  40,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  41,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  42,   0,   0,   0,      0,   0,  43,   0,   0,   0,   0,   0,   0,   0,   0,   0,  44,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  45,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,  46,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  47,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  48,   0,   0,   0,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  49,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  50,   0,      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,  51,   0,   0,   0,   0,   0,   0,   0,   0,   0,  52,   0,   0 } ;

/* . The non-zero CH terms indexed by CHINDICES (-1). */
const Real CHTERMS[52] = {
   1.000000000000000,    1.000000000000000,    1.000000000000000,    1.333333333333333,    1.000000000000000,
   1.000000000000000,    1.000000000000000,   -0.666666666666667,    1.000000000000000,    1.000000000000000,
   1.000000000000000,    1.000000000000000,    1.000000000000000,   -0.666666666666667,   -1.000000000000000,
   1.154700538379300,    1.154700538379300,   -0.577350269189630,   -0.577350269189630,    1.000000000000000,
   1.333333333333333,    1.000000000000000,    1.000000000000000,    1.000000000000000,    0.577350269189630,
   1.000000000000000,    0.666666666666667,    1.000000000000000,    1.000000000000000,    1.000000000000000,
   1.000000000000000,    0.577350269189630,    1.000000000000000,    1.000000000000000,    0.666666666666667,
  -1.000000000000000,    1.000000000000000,    1.000000000000000,   -1.000000000000000,   -1.154700538379300,
   1.000000000000000,   -1.000000000000000,    1.000000000000000,   -1.333333333333333,    1.000000000000000,
   1.000000000000000,    1.000000000000000,   -1.154700538379300,    1.000000000000000,    1.000000000000000,
   1.000000000000000,   -1.333333333333333 } ;
