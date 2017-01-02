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
subroutine BUD_BLAS(BUD_TYPE_VARP,gemv,DIM)(T, M, N, alpha, A, LDA, X, IX, BETA, Y, IY)
  character(len=1), intent(in) :: T
  integer, intent(in) :: M, N
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), intent(in) :: alpha, beta
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,gemv)(T, M, N, alpha, AA, LDA, AX, IX, beta, AY, IY)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,hemv,DIM)(UL, N, alpha, A, LDA, X, IX, BETA, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,hemv)(T, M, N, alpha, AA, LDA, AX, IX, beta, AY, IY)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,symv,DIM)(UL, N, alpha, A, LDA, X, IX, BETA, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,symv)(T, M, N, alpha, AA, LDA, AX, IX, beta, AY, IY)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,hpmv,DIM)(UL, N, alpha, A, LDA, X, IX, BETA, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,hpmv)(T, M, N, alpha, AA, LDA, AX, IX, beta, AY, IY)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,spmv,DIM)(UL, N, alpha, A, LDA, X, IX, BETA, Y, IY)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha, beta

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,spmv)(T, M, N, alpha, AA, LDA, AX, IX, beta, AY, IY)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,geru,DIM)(M, N, alpha, X, IX, Y, IY, A, LDA)
  integer, intent(in) :: M, N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,geru)(M, N, alpha, AX, IX, AY, IY, AA, LDA)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,gerc,DIM)(M, N, alpha, X, IX, Y, IY, A, LDA)
  integer, intent(in) :: M, N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,gerc)(M, N, alpha, AX, IX, AY, IY, AA, LDA)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,ger,DIM)(M, N, alpha, X, IX, Y, IY, A, LDA)
  integer, intent(in) :: M, N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,ger)(M, N, alpha, AX, IX, AY, IY, AA, LDA)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,her,DIM)(UL, N, alpha, X, IX, Y, IY, A, LDA)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,her)(UL, N, alpha, AX, IX, AY, IY, AA, LDA)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,hpr,DIM)(UL, N, alpha, X, IX, Y, IY, A, LDA)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,hpr)(UL, N, alpha, AX, IX, AY, IY, AA, LDA)
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
subroutine BUD_BLAS(BUD_TYPE_VARP,syr,DIM)(UL, N, alpha, X, IX, Y, IY, A, LDA)
  character(len=1), intent(in) :: UL
  integer, intent(in) :: N
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,BUD_DIM,D)), intent(in) :: A
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(in) :: X
  type(BUD_CC5(BUD_TYPE,BUD_TYPE_VAR_P,Array,1,D)), intent(inout) :: Y
  integer, intent(in) :: IX, IY
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AX(:), AY(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC), pointer :: AA BUD_DIM_SPEC(:)
  BUD_TYPE_VAR(BUD_TYPE_VAR_PREC) :: alpha

  AA => array_p(A)
  AX => array_p(X)
  AY => array_p(Y)
  call BUD_BLASX(BUD_TYPE_VARP,syr)(UL, N, alpha, AX, IX, AY, IY, AA, LDA)
end subroutine
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
