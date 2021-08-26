/*
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for DUNEX Coupled Test, waves-ocean (SWAN/ROMS) two-way coupling.
** Based on Coupling Test
** Application flag:   DUNEX
** Input scripts:      coupling_dunex.in
**                     ocean_dunex.in
**                     swan_dunex.in
*/

#define ROMS_MODEL
#define SWAN_MODEL
#define MCT_LIB

#undef SWAN_IUNIFORM

#define UV_KIRBY
#define WDISS_WAVEMOD
#define WEC_VF
#define ROLLER_RENIERS

#define WET_DRY

#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#define MIX_S_TS
#define MASKING
#define UV_ADV
#undef  UV_C2ADVECTION 
#undef  TS_U3HADVECTION
#define TS_HSIMT
#define DJ_GRADPS
#define SOLVE3D
#undef  SPLINES_VVISC 
#undef  SPLINES_VDIFF

#undef  DIAGNOSTICS_UV
#undef  AVERAGES
#undef  AVERAGES_WEC

#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BPFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX
#define ANA_INITIAL 

#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# define SSW_LOGINT 
# define SSW_LOGINT_STOKES 
#endif

#define GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES 
# define ZOS_HSIG
# define TKE_WAVEDISS
# undef  LIMIT_VVISC
# define LIMIT_VDIFF
#endif




