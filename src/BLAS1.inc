!> @param[in] N size of array
!! @param[inout] X array to swap place with Y
!! @param[in] IX X increment 
!! @param[inout] Y array to swap place with X
!! @param[in] IY Y increment 
subroutine BUD_BLAS(BUD_TYPE_VARP,swap,BUD_DIM)(N, X, IX, Y, IY)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,swap)(N, AX, IX, AY, IY)
end subroutine


!> @param[in] N size of array
!! @param[in] A the scaling factor
!! @param[inout] X array to swap place with Y
!! @param[in] IX X increment 
subroutine BUD_BLAS(BUD_TYPE_VARP,scal,BUD_DIM)(N, A, X, IX)
  integer, intent(in) :: N
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  integer, intent(in) :: IX
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)

  AX => array_p(X)
  call BUD_BLASX(BUD_TYPE_VARP,scal)(N, A, AX, IX)
end subroutine

!> @param[in] N size of array
!! @param[in] X array to copy into Y
!! @param[in] IX X increment 
!! @param[out] Y copied array on output
!! @param[in] IY Y increment 
subroutine BUD_BLAS(BUD_TYPE_VARP,copy,BUD_DIM)(N, X, IX, Y, IY)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,copy)(N, AX, IX, AY, IY)
end subroutine

!> @param[in] N size of array
!! @param[in] A factor on X
!! @param[in] X array to multiply by A and add to Y
!! @param[in] IX X increment 
!! @param[inout] Y output array
!! @param[in] IY Y increment 
subroutine BUD_BLAS(BUD_TYPE_VARP,axpy,BUD_DIM)(N, A, X, IX, Y, IY)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,axpy)(N, A, AX, IX, AY, IY)
end subroutine


#ifdef BUD_IS_COMPLEX

!> @param[in] N size of array
!! @param[in] X 1st array of dot-product
!! @param[in] IX X increment 
!! @param[in] Y 2nd array of dot-product
!! @param[in] IY Y increment 
function BUD_BLAS(BUD_TYPE_VARP,dotu,BUD_DIM)(N, X, IX, Y, IY) result(dot)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_TYPE_VARP,dotu)(N, AX, IX, AY, IY)
end function

!> @param[in] N size of array
!! @param[in] X 1st array of dot-product
!! @param[in] IX X increment 
!! @param[in] Y 2nd array of dot-product
!! @param[in] IY Y increment 
function BUD_BLAS(BUD_TYPE_VARP,dotc,BUD_DIM)(N, X, IX, Y, IY) result(dot)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_TYPE_VARP,dotc)(N, AX, IX, AY, IY)
end function


#else


!> @param[in] N size of array
!! @param[in] X 1st array of dot-product
!! @param[in] IX X increment 
!! @param[in] Y 2nd array of dot-product
!! @param[in] IY Y increment 
function BUD_BLAS(BUD_TYPE_VARP,dot,BUD_DIM)(N, X, IX, Y, IY) result(dot)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_TYPE_VARP,dot)(N, AX, IX, AY, IY)
end function

#endif

!> @param[in] N size of array
!! @param[in] X array to take euclidean norm of
!! @param[in] IX X increment 
function BUD_BLAS(BUD_TYPE_VARP,nrm2,BUD_DIM)(N, X, IX) result(nrm2)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  integer, intent(in) :: IX
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: nrm2

  AX => array_p(X)
  nrm2 = BUD_BLASX(BUD_TYPE_VARP,nrm2)(N, AX, IX)
end function


!> @param[in] N size of array
!! @param[in] X array to take absolute sum of
!! @param[in] IX X increment 
function BUD_BLAS(BUD_TYPE_VARP,asum,BUD_DIM)(N, X, IX) result(asum)
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  integer, intent(in) :: IX
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: asum

  AX => array_p(X)
  asum = BUD_BLASX(BUD_TYPE_VARP,asum)(N, AX, IX)
end function


#if BUD_DIM == 1


!> @param[inout] X array to swap place with Y
!! @param[inout] Y array to swap place with X
subroutine BUD_BLAS(BUD_TYPE_VARP,swap_dir,BUD_DIM)(X, Y)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,swap)(size(AX), AX, 1, AY, 1)
end subroutine

!> @param[in] A the scaling factor
!! @param[inout] X array to scale
subroutine BUD_BLAS(BUD_TYPE_VARP,scal_dir,BUD_DIM)(A, X)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)

  AX => array_p(X)
  call BUD_BLASX(BUD_TYPE_VARP,scal)(size(AX), A, AX, 1)
end subroutine

!> @param[in] A factor on X
!! @param[in] X array to multiply by A and add to Y
!! @param[inout] Y result array
subroutine BUD_BLAS(BUD_TYPE_VARP,axpy_dir,BUD_DIM)(A, X, Y)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: A
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)

  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,axpy)(size(AX), A, AX, 1, AY, 1)
end subroutine

# ifdef BUD_IS_COMPLEX

!> @param[in] X 1st array of dot-product
!! @param[in] Y 2nd array of dot-product
function BUD_BLAS(BUD_TYPE_VARP,dotu_dir,BUD_DIM)(X, Y) result(dot)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_TYPE_VARP,dotu)(size(AX), AX, 1, AY, 1)
end function

!> @param[in] X 1st array of dot-product
!! @param[in] Y 2nd array of dot-product
function BUD_BLAS(BUD_TYPE_VARP,dotc_dir,BUD_DIM)(X, Y) result(dot)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_TYPE_VARP,dotc)(size(AX), AX, 1, AY, 1)
end function
# else
!> @param[in] X 1st array of dot-product
!! @param[in] Y 2nd array of dot-product
function BUD_BLAS(BUD_TYPE_VARP,dot_dir,BUD_DIM)(X, Y) result(dot)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: Y
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:), AY BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: dot

  AX => array_p(X)
  AY => array_p(Y)
  dot = BUD_BLASX(BUD_TYPE_VARP,dot)(size(AX), AX, 1, AY, 1)
end function
# endif

!> @param[in] X array to take euclidean norm of
function BUD_BLAS(BUD_TYPE_VARP,nrm2_dir,BUD_DIM)(X) result(nrm2)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: nrm2

  AX => array_p(X)
  nrm2 = BUD_BLASX(BUD_TYPE_VARP,nrm2)(size(AX), AX, 1)
end function

!> @param[in] X array to take absolute sum of
function BUD_BLAS(BUD_TYPE_VARP,asum_dir,BUD_DIM)(X) result(asum)
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(inout) :: X
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: asum

  AX => array_p(X)
  asum = BUD_BLASX(BUD_TYPE_VARP,asum)(size(AX), AX, 1)
end function


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

