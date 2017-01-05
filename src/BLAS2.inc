! This is required for updating the BUD_DIM_SPEC
#include "bud_utils.inc"

#ifndef BUD_BLAS_VARP
# define BUD_BLAS_VARP BUD_TYPE_VARP
#endif


!> @param[in] T transformation of matrix
!! @param[in] M rows of matrix A
!! @param[in] N columns of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[inout] Y matrix-vector result array
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,gemv,BUD_DIM)(T, M, N, alpha, A, LDA, X, IX, beta, Y, IY)
  character(len=1), intent(in) :: T
  integer, intent(in) :: M, N, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,gemv)(T, M, N, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AX(1), IX, beta, AY(1), IY)
end subroutine

#ifdef BUD_IS_COMPLEX
!> @param[in] UL upper-lower part of matrix used for multiplication
!! @param[in] N dimension of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[inout] Y matrix-vector result array
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,hemv,BUD_DIM)(UL, N, alpha, A, LDA, X, IX, beta, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N, LDA
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,hemv)(UL, N, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AX(1), IX, beta, AY(1), IY)
end subroutine

#else
!> @param[in] UL upper-lower part of matrix used for multiplication
!! @param[in] N dimension of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[inout] Y matrix-vector result array
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,symv,BUD_DIM)(UL, N, alpha, A, LDA, X, IX, beta, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N, LDA
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,symv)(UL, N, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AX(1), IX, beta, AY(1), IY)
end subroutine
#endif


#ifdef BUD_IS_COMPLEX
!> @param[in] UL upper-lower part of matrix used for multiplication
!! @param[in] N dimension of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A in packed format
!! @param[in] LDA leading dimension of A
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[inout] Y matrix-vector result array
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,hpmv,BUD_DIM)(UL, N, alpha, A, LDA, X, IX, beta, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N, LDA
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,hpmv)(UL, N, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AX(1), IX, beta, AY(1), IY)
end subroutine

#else
!> @param[in] UL upper-lower part of matrix used for multiplication
!! @param[in] N dimension of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A in packed format
!! @param[in] LDA leading dimension of A
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[inout] Y matrix-vector result array
!! @param[in] IY Y increment
subroutine BUD_BLAS(BUD_BLAS_VARP,spmv,BUD_DIM)(UL, N, alpha, A, LDA, X, IX, beta, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N, LDA
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,spmv)(UL, N, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AX(1), IX, beta, AY(1), IY)
end subroutine
#endif


#ifdef BUD_IS_COMPLEX
!> @param[in] M rows of matrix A
!! @param[in] N columns of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[in] Y matrix-vector result array
!! @param[in] IY Y increment
!! @param[inout] A matrix A
!! @param[in] LDA leading dimension of A
subroutine BUD_BLAS(BUD_BLAS_VARP,geru,BUD_DIM)(M, N, alpha, X, IX, Y, IY, A, LDA)
  integer, intent(in) :: M, N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,geru)(M, N, alpha, &
    AX(1), IX, AY(1), IY, AA BUD_DIM_SPEC(1), LDA)
end subroutine

!> @param[in] M rows of matrix A
!! @param[in] N columns of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[in] Y matrix-vector result array
!! @param[in] IY Y increment
!! @param[inout] A matrix A
!! @param[in] LDA leading dimension of A
subroutine BUD_BLAS(BUD_BLAS_VARP,gerc,BUD_DIM)(M, N, alpha, X, IX, Y, IY, A, LDA)
  integer, intent(in) :: M, N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,gerc)(M, N, alpha, &
    AX(1), IX, AY(1), IY, AA BUD_DIM_SPEC(1), LDA)
end subroutine

#else
!> @param[in] M rows of matrix A
!! @param[in] N columns of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[in] Y matrix-vector result array
!! @param[in] IY Y increment
!! @param[inout] A matrix A
!! @param[in] LDA leading dimension of A
subroutine BUD_BLAS(BUD_BLAS_VARP,ger,BUD_DIM)(M, N, alpha, X, IX, Y, IY, A, LDA)
  integer, intent(in) :: M, N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,ger)(M, N, alpha, &
    AX(1), IX, AY(1), IY, AA BUD_DIM_SPEC(1), LDA)
end subroutine
#endif


#ifdef BUD_IS_COMPLEX
!> @param[in] UL upper or lower part of matrix A
!! @param[in] N size of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[in] Y matrix-vector result array
!! @param[in] IY Y increment
!! @param[inout] A matrix A
!! @param[in] LDA leading dimension of A
subroutine BUD_BLAS(BUD_BLAS_VARP,her,BUD_DIM)(UL, N, alpha, X, IX, Y, IY, A, LDA)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,her)(UL, N, alpha, &
    AX(1), IX, AY(1), IY, AA BUD_DIM_SPEC(1), LDA)
end subroutine

!> @param[in] UL upper or lower part of matrix A
!! @param[in] N size of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[in] Y matrix-vector result array
!! @param[in] IY Y increment
!! @param[inout] A matrix A
!! @param[in] LDA leading dimension of A
subroutine BUD_BLAS(BUD_BLAS_VARP,hpr,BUD_DIM)(UL, N, alpha, X, IX, Y, IY, A, LDA)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,hpr)(UL, N, alpha, &
    AX(1), IX, AY(1), IY, AA BUD_DIM_SPEC(1), LDA)
end subroutine

#else
!> @param[in] UL upper or lower part of matrix A
!! @param[in] N size of matrix A
!! @param[in] alpha factor for multiplication
!! @param[in] X vector for matrix-vector multiplication
!! @param[in] IX X increment
!! @param[in] beta factor for output
!! @param[in] Y matrix-vector result array
!! @param[in] IY Y increment
!! @param[inout] A matrix A
!! @param[in] LDA leading dimension of A
subroutine BUD_BLAS(BUD_BLAS_VARP,syr,BUD_DIM)(UL, N, alpha, X, IX, Y, IY, A, LDA)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY, LDA
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_BLAS_VARP,syr)(UL, N, alpha, &
    AX(1), IX, AY(1), IY, AA BUD_DIM_SPEC(1), LDA)
end subroutine
#endif

#if BUD_DIM == 1
# undef BUD_DIM
# define BUD_DIM 2
# include "BLAS2.inc"
#endif

#if BUD_DIM == 2
# undef BUD_DIM
# undef BUD_BLAS_VARP
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
