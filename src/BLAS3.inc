! This is required for updating the BUD_DIM_SPEC
#include "bud_utils.inc"

#ifndef BUD_BLAS_VARP
# define BUD_BLAS_VARP BUD_TYPE_VARP
#endif


!> @param[in] TA transformation of matrix A
!! @param[in] TB transformation of matrix B
!! @param[in] M rows of matrix op(A)
!! @param[in] N columns of matrix op(B)
!! @param[in] K columns of matrix op(A)
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] B matrix B
!! @param[in] LDB leading dimension of B
!! @param[in] beta factor for multiplication
!! @param[inout] C output matrix
!! @param[in] LDC leading dimension of C
subroutine BUD_BLAS(BUD_BLAS_VARP,gemm,BUD_DIM)(TA, TB, M, N, K, alpha, A, LDA, B, LDB, &
     beta, C, LDC)
  character(len=1), intent(in) :: TA, TB
  integer, intent(in) :: M, N, K
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A, B, C
  integer, intent(in) :: LDA, LDB, LDC
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:), AB BUD_DIM_SPEC(:), AC BUD_DIM_SPEC(:)

  AA => array_p(A)
  AB => array_p(B)
  AC => array_p(C)
  call BUD_BLASX(BUD_BLAS_VARP,gemm)(TA, TB, M, N, K, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AB BUD_DIM_SPEC(1), LDB, beta, AC BUD_DIM_SPEC(1), LDC)
end subroutine


!> @param[in] S side
!! @param[in] UL upper-lower
!! @param[in] M rows of matrix op(A)
!! @param[in] N columns of matrix op(B)
!! @param[in] K columns of matrix op(A)
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] B matrix B
!! @param[in] LDB leading dimension of B
!! @param[in] beta factor for multiplication
!! @param[inout] C output matrix
!! @param[in] LDC leading dimension of C
subroutine BUD_BLAS(BUD_BLAS_VARP,symm,BUD_DIM)(S, UL, M, N, K, alpha, A, LDA, B, LDB, &
     beta, C, LDC)
  character(len=1), intent(in) :: S, UL
  integer, intent(in) :: M, N, K
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A, B, C
  integer, intent(in) :: LDA, LDB, LDC
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:), AB BUD_DIM_SPEC(:), AC BUD_DIM_SPEC(:)

  AA => array_p(A)
  AB => array_p(B)
  AC => array_p(C)
  call BUD_BLASX(BUD_BLAS_VARP,symm)(S, UL, M, N, K, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AB BUD_DIM_SPEC(1), LDB, beta, AC BUD_DIM_SPEC(1), LDC)
end subroutine

!> @param[in] UL upper-lower
!! @param[in] T transformation
!! @param[in] N order of matrix C
!! @param[in] K columns of matrix op(A)
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] beta factor for multiplication
!! @param[inout] C output matrix
!! @param[in] LDC leading dimension of C
subroutine BUD_BLAS(BUD_BLAS_VARP,syrk,BUD_DIM)(UL, T, N, K, alpha, A, LDA, beta, C, LDC)
  character(len=1), intent(in) :: UL, T
  integer, intent(in) :: N, K
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A, C
  integer, intent(in) :: LDA, LDC
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:), AC BUD_DIM_SPEC(:)

  AA => array_p(A)
  AC => array_p(C)
  call BUD_BLASX(BUD_BLAS_VARP,syrk)(UL, T, N, K, alpha, &
    AA BUD_DIM_SPEC(1), LDA, beta, AC BUD_DIM_SPEC(1), LDC)
end subroutine

#ifdef BUD_IS_COMPLEX

!> @param[in] S side
!! @param[in] UL upper-lower
!! @param[in] M rows of matrix C
!! @param[in] N columns of matrix C
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] B matrix B
!! @param[in] LDB leading dimension of B
!! @param[in] beta factor for multiplication
!! @param[inout] C output matrix
!! @param[in] LDC leading dimension of C
subroutine BUD_BLAS(BUD_BLAS_VARP,hemm,BUD_DIM)(S, UL, M, N, alpha, A, LDA, B, LDB, &
     beta, C, LDC)
  character(len=1), intent(in) :: S, UL
  integer, intent(in) :: M, N
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A, B, C
  integer, intent(in) :: LDA, LDB, LDC
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:), AB BUD_DIM_SPEC(:), AC BUD_DIM_SPEC(:)

  AA => array_p(A)
  AB => array_p(B)
  AC => array_p(C)
  call BUD_BLASX(BUD_BLAS_VARP,hemm)(S, UL, M, N, alpha, &
    AA BUD_DIM_SPEC(1), LDA, AB BUD_DIM_SPEC(1), LDB, beta, AC BUD_DIM_SPEC(1), LDC)
end subroutine

!> @param[in] UL upper-lower
!! @param[in] T transformation
!! @param[in] N order of matrix C
!! @param[in] K columns of matrix op(A)
!! @param[in] alpha factor for multiplication
!! @param[in] A matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] beta factor for multiplication
!! @param[inout] C output matrix
!! @param[in] LDC leading dimension of C
subroutine BUD_BLAS(BUD_BLAS_VARP,herk,BUD_DIM)(UL, T, N, K, alpha, A, LDA, beta, C, LDC)
  character(len=1), intent(in) :: UL, T
  integer, intent(in) :: N, K
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VARP,Array,BUD_DIM,D)), intent(in) :: A, C
  integer, intent(in) :: LDA, LDC
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:), AC BUD_DIM_SPEC(:)

  AA => array_p(A)
  AC => array_p(C)
  call BUD_BLASX(BUD_BLAS_VARP,herk)(UL, T, N, K, alpha, &
    AA BUD_DIM_SPEC(1), LDA, beta, AC BUD_DIM_SPEC(1), LDC)
end subroutine

#endif

#if BUD_DIM == 1
# undef BUD_DIM
# define BUD_DIM 2
# include "BLAS3.inc"
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
