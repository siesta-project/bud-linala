! @@LICENSE@@ see Copyright notice in the top-directory

! Get default commands
#include "bud_utils.inc"
#include "BLAS.inc"


!> @defgroup BLAS-1 BLAS level 1 routines
!! @ingroup LA
!! BLAS level 1 routines.
!! @{

#define BUD_MOD_NAME BUD_CC3(BUD_MOD,_,BLAS1)
module BUD_MOD_NAME

#include "bud_common_declarations.inc"

  ! Define interfaces
  interface swap
    module procedure BUD_BLAS(s,swap,1)
    module procedure BUD_BLAS(s,swap_dir,1)
    module procedure BUD_BLAS(s,swap,2)
    module procedure BUD_BLAS(d,swap,1)
    module procedure BUD_BLAS(d,swap_dir,1)
    module procedure BUD_BLAS(d,swap,2)
    module procedure BUD_BLAS(c,swap,1)
    module procedure BUD_BLAS(c,swap_dir,2)
    module procedure BUD_BLAS(d,swap,2)
    module procedure BUD_BLAS(z,swap,1)
    module procedure BUD_BLAS(z,swap_dir,1)
    module procedure BUD_BLAS(z,swap,2)
  end interface 

  interface scal
    module procedure BUD_BLAS(s,scal,1)
    module procedure BUD_BLAS(s,scal_dir,1)
    module procedure BUD_BLAS(s,scal,2)
    module procedure BUD_BLAS(d,scal,1)
    module procedure BUD_BLAS(d,scal_dir,1)
    module procedure BUD_BLAS(d,scal,2)
    module procedure BUD_BLAS(c,scal,1)
    module procedure BUD_BLAS(c,scal_dir,1)
    module procedure BUD_BLAS(c,scal,2)
    module procedure BUD_BLAS(z,scal,1)
    module procedure BUD_BLAS(z,scal_dir,1)
    module procedure BUD_BLAS(z,scal,2)
  end interface

  interface copy
    module procedure BUD_BLAS(s,copy,1)
    module procedure BUD_BLAS(s,copy,2)
    module procedure BUD_BLAS(d,copy,1)
    module procedure BUD_BLAS(d,copy,2)
    module procedure BUD_BLAS(c,copy,1)
    module procedure BUD_BLAS(c,copy,2)
    module procedure BUD_BLAS(z,copy,1)
    module procedure BUD_BLAS(z,copy,2)
  end interface 

  interface axpy
    module procedure BUD_BLAS(s,axpy,1)
    module procedure BUD_BLAS(s,axpy_dir,1)
    module procedure BUD_BLAS(s,axpy,2)
    module procedure BUD_BLAS(d,axpy,1)
    module procedure BUD_BLAS(d,axpy_dir,1)
    module procedure BUD_BLAS(d,axpy,2)
    module procedure BUD_BLAS(c,axpy,1)
    module procedure BUD_BLAS(c,axpy_dir,1)
    module procedure BUD_BLAS(c,axpy,2)
    module procedure BUD_BLAS(z,axpy,1)
    module procedure BUD_BLAS(z,axpy_dir,1)
    module procedure BUD_BLAS(z,axpy,2)
  end interface 

  interface dot
    module procedure BUD_BLAS(s,dot,1)
    module procedure BUD_BLAS(s,dot_dir,1)
    module procedure BUD_BLAS(s,dot,2)
    module procedure BUD_BLAS(d,dot,1)
    module procedure BUD_BLAS(d,dot_dir,1)
    module procedure BUD_BLAS(d,dot,2)
  end interface 

  interface dotu
    module procedure BUD_BLAS(c,dotu,1)
    module procedure BUD_BLAS(c,dotu_dir,1)
    module procedure BUD_BLAS(c,dotu,2)
    module procedure BUD_BLAS(z,dotu,1)
    module procedure BUD_BLAS(z,dotu_dir,1)
    module procedure BUD_BLAS(z,dotu,2)
  end interface 

  interface dotc
    module procedure BUD_BLAS(c,dotc,1)
    module procedure BUD_BLAS(c,dotc_dir,1)
    module procedure BUD_BLAS(c,dotc,2)
    module procedure BUD_BLAS(z,dotc,1)
    module procedure BUD_BLAS(z,dotc_dir,1)
    module procedure BUD_BLAS(z,dotc,2)
  end interface 

  interface nrm2
    module procedure BUD_BLAS(s,nrm2,1)
    module procedure BUD_BLAS(s,nrm2_dir,1)
    module procedure BUD_BLAS(s,nrm2,2)
    module procedure BUD_BLAS(d,nrm2,1)
    module procedure BUD_BLAS(d,nrm2_dir,1)
    module procedure BUD_BLAS(d,nrm2,2)
    module procedure BUD_BLAS(c,nrm2,1)
    module procedure BUD_BLAS(c,nrm2_dir,1)
    module procedure BUD_BLAS(c,nrm2,2)
    module procedure BUD_BLAS(z,nrm2,1)
    module procedure BUD_BLAS(z,nrm2_dir,1)
    module procedure BUD_BLAS(z,nrm2,2)
  end interface 

  interface asum
    module procedure BUD_BLAS(s,asum,1)
    module procedure BUD_BLAS(s,asum_dir,1)
    module procedure BUD_BLAS(s,asum,2)
    module procedure BUD_BLAS(d,asum,1)
    module procedure BUD_BLAS(d,asum_dir,1)
    module procedure BUD_BLAS(d,asum,2)
    module procedure BUD_BLAS(c,asum,1)
    module procedure BUD_BLAS(c,asum_dir,1)
    module procedure BUD_BLAS(c,asum,2)
    module procedure BUD_BLAS(z,asum,1)
    module procedure BUD_BLAS(z,asum_dir,1)
    module procedure BUD_BLAS(z,asum,2)
  end interface 

contains

#define BUD_TYPE_VAR real
#define BUD_TYPE_VAR_PREC rr_
#define BUD_TYPE_VARP s
#define BUD_DIM 1
#include "BLAS1.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS1.inc"
#undef BUD_TYPE_VAR_PREC
#undef BUD_TYPE_VARP
#define BUD_TYPE_VAR_PREC rd_
#define BUD_TYPE_VARP d
#include "BLAS1.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS1.inc"

#undef BUD_TYPE_VAR
#define BUD_IS_COMPLEX
#define BUD_TYPE_VAR complex
#define BUD_TYPE_VAR_PREC rr_
#define BUD_TYPE_VARP c
#include "BLAS1.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS1.inc"
#undef BUD_TYPE_VAR_PREC
#undef BUD_TYPE_VARP
#define BUD_TYPE_VAR_PREC rd_
#define BUD_TYPE_VARP z
#include "BLAS1.inc"
#undef BUD_DIM
#define BUD_DIM 2
#include "BLAS1.inc"
  
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

