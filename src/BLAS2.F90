! @@LICENSE@@ see Copyright notice in the top-directory

! Get default commands
#include "bud_utils.inc"
#include "BLAS.inc"


!> @defgroup BLAS-2 BLAS level 2 routines
!! @ingroup LA
!! BLAS level 2 routines.
!! @{

#define BUD_MOD_NAME BUD_CC3(BUD_MOD,_,BLAS2)
module BUD_MOD_NAME

  use BUD_CC3(BUD_MOD,_,rArray1D)
  use BUD_CC3(BUD_MOD,_,rArray2D)
  use BUD_CC3(BUD_MOD,_,dArray1D)
  use BUD_CC3(BUD_MOD,_,dArray2D)
  use BUD_CC3(BUD_MOD,_,cArray1D)
  use BUD_CC3(BUD_MOD,_,cArray2D)
  use BUD_CC3(BUD_MOD,_,zArray1D)
  use BUD_CC3(BUD_MOD,_,zArray2D)

#include "bud_common_declarations.inc"

  ! Define interfaces
  interface gemv
    module procedure BUD_BLAS(s,gemv,1)
    module procedure BUD_BLAS(s,gemv,2)
    module procedure BUD_BLAS(d,gemv,1)
    module procedure BUD_BLAS(d,gemv,2)
    module procedure BUD_BLAS(c,gemv,1)
    module procedure BUD_BLAS(c,gemv,2)
    module procedure BUD_BLAS(z,gemv,1)
    module procedure BUD_BLAS(z,gemv,2)
  end interface 

  interface hemv
    module procedure BUD_BLAS(c,hemv,1)
    module procedure BUD_BLAS(c,hemv,2)
    module procedure BUD_BLAS(z,hemv,1)
    module procedure BUD_BLAS(z,hemv,2)
  end interface

  interface hpmv
    module procedure BUD_BLAS(c,hpmv,1)
    module procedure BUD_BLAS(c,hpmv,2)
    module procedure BUD_BLAS(z,hpmv,1)
    module procedure BUD_BLAS(z,hpmv,2)
  end interface

  interface geru
    module procedure BUD_BLAS(c,geru,1)
    module procedure BUD_BLAS(c,geru,2)
    module procedure BUD_BLAS(z,geru,1)
    module procedure BUD_BLAS(z,geru,2)
  end interface

  interface gerc
    module procedure BUD_BLAS(c,gerc,1)
    module procedure BUD_BLAS(c,gerc,2)
    module procedure BUD_BLAS(z,gerc,1)
    module procedure BUD_BLAS(z,gerc,2)
  end interface

  interface her
    module procedure BUD_BLAS(c,her,1)
    module procedure BUD_BLAS(c,her,2)
    module procedure BUD_BLAS(z,her,1)
    module procedure BUD_BLAS(z,her,2)
  end interface

  interface hpr
    module procedure BUD_BLAS(c,hpr,1)
    module procedure BUD_BLAS(c,hpr,2)
    module procedure BUD_BLAS(z,hpr,1)
    module procedure BUD_BLAS(z,hpr,2)
  end interface

  interface symv
    module procedure BUD_BLAS(s,symv,1)
    module procedure BUD_BLAS(s,symv,2)
    module procedure BUD_BLAS(d,symv,1)
    module procedure BUD_BLAS(d,symv,2)
  end interface
  
  interface spmv
    module procedure BUD_BLAS(s,spmv,1)
    module procedure BUD_BLAS(s,spmv,2)
    module procedure BUD_BLAS(d,spmv,1)
    module procedure BUD_BLAS(d,spmv,2)
  end interface

  interface ger
    module procedure BUD_BLAS(s,ger,1)
    module procedure BUD_BLAS(s,ger,2)
    module procedure BUD_BLAS(d,ger,1)
    module procedure BUD_BLAS(d,ger,2)
  end interface

  interface syr
    module procedure BUD_BLAS(s,syr,1)
    module procedure BUD_BLAS(s,syr,2)
    module procedure BUD_BLAS(d,syr,1)
    module procedure BUD_BLAS(d,syr,2)
  end interface

contains

#define BUD_TYPE_VAR real
#define BUD_TYPE_VAR_PREC rr_
#define BUD_BLAS_VARP s
#define BUD_TYPE_VARP r
#define BUD_DIM 1
#include "BLAS2.inc"
#undef BUD_TYPE_VAR_PREC
#undef BUD_TYPE_VARP
#undef BUD_BLAS_VARP
#define BUD_TYPE_VAR_PREC rd_
#define BUD_TYPE_VARP d
#define BUD_DIM 1
#include "BLAS2.inc"

#undef BUD_TYPE_VAR
#define BUD_IS_COMPLEX
#define BUD_TYPE_VAR complex
#define BUD_TYPE_VAR_PREC rr_
#define BUD_TYPE_VARP c
#define BUD_DIM 1
#include "BLAS2.inc"
#undef BUD_TYPE_VAR_PREC
#undef BUD_TYPE_VARP
#define BUD_TYPE_VAR_PREC rd_
#define BUD_TYPE_VARP z
#define BUD_DIM 1
#include "BLAS2.inc"
  
end module
!> @}


! project-buds -- local file settings
!     Anything below this line may be overwritten by scripts
!     Below are non-editable settings

! Local Variables:
!  mode: f90
!  f90-if-indent: 2
!  f90-type-indent: 2
!  f90-associate-indent: 2
!  f90-continuation-indent: 2
!  f90-structure-indent: 2
!  f90-critical-indent: 2
!  f90-program-indent: 2
!  f90-do-indent: 2
! End:

