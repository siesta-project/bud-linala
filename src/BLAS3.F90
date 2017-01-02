! @@LICENSE@@ see Copyright notice in the top-directory

! Get default commands
#include "bud_utils.inc"
#include "BLAS.inc"


!> @defgroup BLAS-3 BLAS level 3 routines
!! @ingroup LA
!! BLAS level 3 routines.
!! @{

#define BUD_MOD_NAME BUD_CC3(BUD_MOD,_,BLAS3)
module BUD_MOD_NAME

#include "bud_common_declarations.inc"

  ! Define interfaces
  interface gemm
    module procedure BUD_BLAS(s,gemm,1)
    module procedure BUD_BLAS(s,gemm,2)
    module procedure BUD_BLAS(d,gemm,1)
    module procedure BUD_BLAS(d,gemm,2)
    module procedure BUD_BLAS(c,gemm,1)
    module procedure BUD_BLAS(c,gemm,2)
    module procedure BUD_BLAS(z,gemm,1)
    module procedure BUD_BLAS(z,gemm,2)
  end interface 

  interface symm
    module procedure BUD_BLAS(s,symm,1)
    module procedure BUD_BLAS(s,symm,2)
    module procedure BUD_BLAS(d,symm,1)
    module procedure BUD_BLAS(d,symm,2)
    module procedure BUD_BLAS(c,symm,1)
    module procedure BUD_BLAS(c,symm,2)
    module procedure BUD_BLAS(z,symm,1)
    module procedure BUD_BLAS(z,symm,2)
  end interface 

  interface syrk
    module procedure BUD_BLAS(s,syrk,1)
    module procedure BUD_BLAS(s,syrk,2)
    module procedure BUD_BLAS(d,syrk,1)
    module procedure BUD_BLAS(d,syrk,2)
    module procedure BUD_BLAS(c,syrk,1)
    module procedure BUD_BLAS(c,syrk,2)
    module procedure BUD_BLAS(z,syrk,1)
    module procedure BUD_BLAS(z,syrk,2)
  end interface

  interface hemm
    module procedure BUD_BLAS(c,hemm,1)
    module procedure BUD_BLAS(c,hemm,2)
    module procedure BUD_BLAS(z,hemm,1)
    module procedure BUD_BLAS(z,hemm,2)
  end interface

  interface herk
    module procedure BUD_BLAS(c,herk,1)
    module procedure BUD_BLAS(c,herk,2)
    module procedure BUD_BLAS(z,herk,1)
    module procedure BUD_BLAS(z,herk,2)
  end interface

contains

#define BUD_TYPE_VAR real
#define BUD_TYPE_VAR_PREC rr_
#define BUD_TYPE_VARP s
#define BUD_DIM 1
#include "BLAS3.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS3.inc"
#undef BUD_TYPE_VAR_PREC
#undef BUD_TYPE_VARP
#define BUD_TYPE_VAR_PREC rd_
#define BUD_TYPE_VARP d
#define BUD_DIM 1
#include "BLAS3.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS3.inc"

#undef BUD_TYPE_VAR
#define BUD_IS_COMPLEX
#define BUD_TYPE_VAR complex
#define BUD_TYPE_VAR_PREC rr_
#define BUD_TYPE_VARP c
#define BUD_DIM 1
#include "BLAS3.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS3.inc"
#undef BUD_TYPE_VAR_PREC
#undef BUD_TYPE_VARP
#define BUD_TYPE_VAR_PREC rd_
#define BUD_TYPE_VARP z
#define BUD_DIM 1
#include "BLAS3.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS3.inc"
  
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

