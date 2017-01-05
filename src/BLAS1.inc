! This is required for updating the BUD_DIM_SPEC
#include "bud_utils.inc"

#ifndef BUD_BLAS_VARP
# define BUD_BLAS_VARP BUD_TYPE_VARP
#endif
#ifndef BUD_BLAS_PRE
# define BUD_BLAS_PRE
#endif


!> @param[in] N size of array
!! @param[inout] X array to swap place with Y
!! @param[in] IX X increment
!! @param[inout] Y array to swap place with X
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,swap,BUD_DIM)(N, X, IX, Y, IY)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,swap)(N, AX BUD_DIM_SPEC(1), IX, AY BUD_DIM_SPEC(1), IY)
end subroutine


!> @param[in] N size of array
!! @param[in] A the scaling factor
!! @param[inout] X array to swap place with Y
!! @param[in] IX X increment
subroutine BUD_BLAS(BUD_BLAS_VARP,scal,BUD_DIM)(N, A, X, IX)
  integer, intent(in) :: N
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  integer, intent(in) :: IX
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)

  AX => array_p(X)
  call BUD_BLASX(BUD_BLAS_VARP,scal)(N, A, AX BUD_DIM_SPEC(1), IX)
end subroutine

!> @param[in] N size of array
!! @param[in] X array to copy into Y
!! @param[in] IX X increment
!! @param[out] Y copied array on output
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,copy,BUD_DIM)(N, X, IX, Y, IY)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,copy)(N, AX BUD_DIM_SPEC(1), IX, AY BUD_DIM_SPEC(1), IY)
end subroutine

!> @param[in] N size of array
!! @param[in] A factor on X
!! @param[in] X array to multiply by A and add to Y
!! @param[in] IX X increment
!! @param[inout] Y output array
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,axpy,BUD_DIM)(N, A, X, IX, Y, IY)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,axpy)(N, A, AX BUD_DIM_SPEC(1), IX, AY BUD_DIM_SPEC(1), IY)
end subroutine


#ifdef BUD_IS_COMPLEX

!> @param[in] N size of array
!! @param[in] X 1st array of dot-product
!! @param[in] IX X increment
!! @param[in] Y 2nd array of dot-product
!! @param[in] IY Y increment
function BUD_BLAS(BUD_BLAS_VARP,dotu,BUD_DIM)(N, X, IX, Y, IY) result(dot)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_BLAS_VARP,dotu)

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_BLAS_VARP,dotu)(N, AX BUD_DIM_SPEC(1), IX, AY BUD_DIM_SPEC(1), IY)
end function

!> @param[in] N size of array
!! @param[in] X 1st array of dot-product
!! @param[in] IX X increment
!! @param[in] Y 2nd array of dot-product
!! @param[in] IY Y increment
function BUD_BLAS(BUD_BLAS_VARP,dotc,BUD_DIM)(N, X, IX, Y, IY) result(dot)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_BLAS_VARP,dotc)

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_BLAS_VARP,dotc)(N, AX BUD_DIM_SPEC(1), IX, AY BUD_DIM_SPEC(1), IY)
end function


#else


!> @param[in] N size of array
!! @param[in] X 1st array of dot-product
!! @param[in] IX X increment
!! @param[in] Y 2nd array of dot-product
!! @param[in] IY Y increment
function BUD_BLAS(BUD_BLAS_VARP,dot,BUD_DIM)(N, X, IX, Y, IY) result(dot)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_BLAS_VARP,dot)

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_BLAS_VARP,dot)(N, AX BUD_DIM_SPEC(1), IX, AY BUD_DIM_SPEC(1), IY)
end function

#endif

!> @param[in] N size of array
!! @param[in] X array to take euclidean norm of
!! @param[in] IX X increment
function BUD_BLAS(BUD_BLAS_VARP,nrm2,BUD_DIM)(N, X, IX) result(nrm2)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  integer, intent(in) :: IX
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: nrm2
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),nrm2)

  AX => array_p(X)
  nrm2 = BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),nrm2)(N, AX BUD_DIM_SPEC(1), IX)
end function


!> @param[in] N size of array
!! @param[in] X array to take absolute sum of
!! @param[in] IX X increment
function BUD_BLAS(BUD_BLAS_VARP,asum,BUD_DIM)(N, X, IX) result(asum)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  integer, intent(in) :: IX
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: asum
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),asum)

  AX => array_p(X)
  asum = BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),asum)(N, AX BUD_DIM_SPEC(1), IX)
end function


#if BUD_DIM == 1


!> @param[inout] X array to swap place with Y
!! @param[inout] Y array to swap place with X
subroutine BUD_BLAS(BUD_BLAS_VARP,swap_dir,BUD_DIM)(X, Y)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,swap)(size(AX), AX BUD_DIM_SPEC(1), 1, AY BUD_DIM_SPEC(1), 1)
end subroutine

!> @param[in] A the scaling factor
!! @param[inout] X array to scale
subroutine BUD_BLAS(BUD_BLAS_VARP,scal_dir,BUD_DIM)(A, X)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)

  AX => array_p(X)
  call BUD_BLASX(BUD_BLAS_VARP,scal)(size(AX), A, AX BUD_DIM_SPEC(1), 1)
end subroutine

!> @param[in] A factor on X
!! @param[in] X array to multiply by A and add to Y
!! @param[inout] Y result array
subroutine BUD_BLAS(BUD_BLAS_VARP,axpy_dir,BUD_DIM)(A, X, Y)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,axpy)(size(AX), A, AX BUD_DIM_SPEC(1), 1, AY BUD_DIM_SPEC(1), 1)
end subroutine

# ifdef BUD_IS_COMPLEX

!> @param[in] X 1st array of dot-product
!! @param[in] Y 2nd array of dot-product
function BUD_BLAS(BUD_BLAS_VARP,dotu_dir,BUD_DIM)(X, Y) result(dot)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_BLAS_VARP,dotu)

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_BLAS_VARP,dotu)(size(AX), AX BUD_DIM_SPEC(1), 1, AY BUD_DIM_SPEC(1), 1)
end function

!> @param[in] X 1st array of dot-product
!! @param[in] Y 2nd array of dot-product
function BUD_BLAS(BUD_BLAS_VARP,dotc_dir,BUD_DIM)(X, Y) result(dot)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_BLAS_VARP,dotc)

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_BLAS_VARP,dotc)(size(AX), AX BUD_DIM_SPEC(1), 1, AY BUD_DIM_SPEC(1), 1)
end function
# else
!> @param[in] X 1st array of dot-product
!! @param[in] Y 2nd array of dot-product
function BUD_BLAS(BUD_BLAS_VARP,dot_dir,BUD_DIM)(X, Y) result(dot)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_BLAS_VARP,dot)

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_BLAS_VARP,dot)(size(AX), AX BUD_DIM_SPEC(1), 1, AY BUD_DIM_SPEC(1), 1)
end function
# endif

!> @param[in] X array to take euclidean norm of
function BUD_BLAS(BUD_BLAS_VARP,nrm2_dir,BUD_DIM)(X) result(nrm2)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: nrm2
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),nrm2)

  AX => array_p(X)
  nrm2 = BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),nrm2)(size(AX), AX BUD_DIM_SPEC(1), 1)
end function

!> @param[in] X array to take absolute sum of
function BUD_BLAS(BUD_BLAS_VARP,asum_dir,BUD_DIM)(X) result(asum)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(inout) :: X
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: asum
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), external :: BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),asum)

  AX => array_p(X)
  asum = BUD_BLASX(BUD_CC2(BUD_BLAS_PRE,BUD_BLAS_VARP),asum)(size(AX), AX BUD_DIM_SPEC(1), 1)
end function

#endif

#if BUD_DIM == 1
# undef BUD_DIM
# define BUD_DIM 2
# include "BLAS1.inc"
#endif

#if BUD_DIM == 2
# undef BUD_DIM
# undef BUD_BLAS_VARP
# undef BUD_BLAS_PRE
# undef BUD_TYPE_VAR_PREC
# undef BUD_TYPE_VARP
#endif


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

