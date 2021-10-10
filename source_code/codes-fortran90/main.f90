program main
use plots


  !-------------------------------------------------------------------------!
  !                                                                         !
  !                     ALPHA VERSION - N/W(100) system                     !
  !                                                                         !
  !-------------------------------------------------------------------------!
  !  ==========================================================================
  !  This program uses the computes potential energy as a function of the N/W
  !  distance for three high-symmetry sites namely N atom located above the top,
  !  bridge and hollow sites.
  !
  !  It computes the Fourier coefficients for the Morse parameters and also
  !  contains a 3-D PES script for the computed potential energy surfaces
  !
  ! ===========================================================================

  implicit none (type, external)
  integer, parameter :: dp = selected_real_kind(p = 15, r = 307) ! Double precision
  external           :: DGESV

  !-------------------------------------------------------------------------!
  !                                                                         !
  !                     INITIALISE THE PARAMETERS                           !
  !                                                                         !
  !-------------------------------------------------------------------------!

  real(dp)   :: delta, Xtop, Ytop, Dtop, alphatop, reqtop
  real(dp)   :: Xhollow, Yhollow, Dhollow, alphahollow, reqhollow
  real(dp)   :: Xbridge, Ybridge, Dbridge, alphabridge, reqbridge
  real(dp)   :: a00,a01,a02,a10,a11,a12,a20,a21,a22
  integer    :: N, NRHS, LDA, LDB, INFO
  integer    :: IPIVOT(3)
  integer    :: N_points, Nbridge
  integer    :: i
  real(dp)   :: A_D(3,3), A_alpha(3,3), A_req(3,3), Coeff_D(3,1), Coeff_alpha(3,1), Coeff_req(3,1)
  real(dp)   :: pi, x, y, x_3d, y_3d
  real(dp)   :: x_start, x_end, y_start, y_end, z_start, z_end
  real(dp), dimension(:)  , allocatable :: Y_arr3d, Z_arr3d, X_arr3d, Z_arr1d

  pi    = 3.141592653589793238462643
  delta = 3.17*10._8**(-10) ! unit cell parameter

  !=========================================================================!
  !                     High Symmetry Site - TOP                            !
  !=========================================================================!
  Xtop        = 0.
  Ytop        = 0.
  Dtop        = 4.9984788
  alphatop    = 1.71646688
  reqtop      = 1.83275503

  !=========================================================================!
  !                     High Symmetry Site - HOLLOW                         !
  !=========================================================================!
  Xhollow      = 0.5
  Yhollow      = 0.5
  Dhollow      = 7.60346595
  alphahollow  = 1.2495198
  reqhollow    = 0.56501317

  !=========================================================================!
  !                     High Symmetry Site - BRIDGE                         !
  !=========================================================================!
  Xbridge      = 0.5
  Ybridge      = 0.
  Dbridge      = 6.55379066
  alphabridge  = 1.30169922
  reqbridge    = 0.95376732

  !=========================================================================!
  !                    COMPUTE FOURIER COEFFICIENTS                         !
  !=========================================================================!
  a00 = 1.
  a10 = 1.
  a20 = 1.

  a01 = cos(2*pi*Xtop/delta)    + cos(2*pi*Ytop/delta)
  a11 = cos(2*pi*Xbridge/delta) + cos(2*pi*Ybridge/delta)
  a21 = cos(2*pi*Xhollow/delta) + cos(2*pi*Yhollow/delta)

  a02 = cos(2*pi*(Xtop+Ytop)/delta)       + cos(2*pi*(Xtop-Ytop)/delta)
  a12 = cos(2*pi*(Xbridge+Ybridge)/delta) + cos(2*pi*(Xbridge-Ybridge)/delta)
  a22 = cos(2*pi*(Xhollow+Yhollow)/delta) + cos(2*pi*(Xhollow-Yhollow)/delta)

  A_D         = reshape([ a00, a10, a20, a01, a11, a21, a02, a12, a22], [ 3, 3 ])
  Coeff_D     = reshape([ Dtop, Dbridge, Dhollow], [ 3, 1 ])
  Coeff_alpha = reshape([ alphatop, alphabridge, alphahollow], [ 3, 1 ])
  Coeff_req   = reshape([ reqtop, reqbridge, reqhollow], [ 3, 1 ])

  ! ===========================================================================
  !     USING EXTERNAL SUBROUTINE DGESV TO SOLVE SYSTEMS OF LINEAR EQUATIONS
  ! ===========================================================================

  !  PURPOSE
  !  ==========================================================================
  !
  !  DGESV computes the solution to a real system of linear equations
  !     A * X = B,
  !  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
  !
  !  The LU decomposition with partial pivoting and row interchanges is
  !  used to factor A as
  !     A = P * L * U,
  !  where P is a permutation matrix, L is unit lower triangular, and U is
  !  upper triangular.  The factored form of A is then used to solve the
  !  system of equations A * X = B.
  ! ===========================================================================
  ! Arguments
  ! ===========
  !
  !  N       (input) INTEGER
  !          The number of linear equations, i.e., the order of the
  !          matrix A.  N >= 0.
  !
  !  NRHS    (input) INTEGER
  !          The number of right hand sides, i.e., the number of columns
  !          of the matrix B.  NRHS >= 0.
  !
  !  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  !          On entry, the N-by-N coefficient matrix A.
  !          On exit, the factors L and U from the factorization
  !          A = P*L*U; the unit diagonal elements of L are not stored.
  !
  !  LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= max(1,N).
  !
  !  IPIV    (output) INTEGER array, dimension (N)
  !          The pivot indices that define the permutation matrix P;
  !          row i of the matrix was interchanged with row IPIV(i).
  !
  !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
  !          On entry, the N-by-NRHS matrix of right hand side matrix B.
  !          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
  !
  !  LDB     (input) INTEGER
  !          The leading dimension of the array B.  LDB >= max(1,N).
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
  !                has been completed, but the factor U is exactly
  !                singular, so the solution could not be computed.
  !
  ! ===========================================================================

  N       = 3
  NRHS    = 1
  LDA     = 3
  LDB     = 3

  !=========================================================================!
  !                FOURIER COEFFICIENTS FOR MORSE PARAMETER D               !
  !=========================================================================!
  CALL DGESV(N, NRHS, A_D, LDA, IPIVOT, Coeff_D, LDB, INFO)

  if (info /= 0) then
    print '(a, i0)', 'Error: ', info
    stop
  end if

  print *, " Coefficients for morse parameter D (P0, P1, P2): "
  print *, Coeff_D

  !=========================================================================!
  !            FOURIER COEFFICIENTS FOR MORSE PARAMETER ALPHA               !
  !=========================================================================!
  A_alpha = reshape([ a00, a10, a20, a01, a11, a21, a02, a12, a22], [ 3, 3 ])

  CALL DGESV(N, NRHS, A_alpha, LDA, IPIVOT, Coeff_alpha, LDB, INFO)

    if (info /= 0) then
      print '(a, i0)', 'Error: ', info
      stop
    end if

  print *, " Coefficients for morse parameter alpha (P0, P1, P2): "
  print *, Coeff_alpha

  !=========================================================================!
  !                FOURIER COEFFICIENTS FOR MORSE PARAMETER R_EQ            !
  !=========================================================================!
  A_req = reshape([ a00, a10, a20, a01, a11, a21, a02, a12, a22], [ 3, 3 ])

  CALL DGESV(N, NRHS, A_req, LDA, IPIVOT, Coeff_req, LDB, INFO)

  if (info /= 0) then
    print '(a, i0)', 'Error: ', info
    stop
  end if

  print *, " Coefficients for morse parameter r_eq (P0, P1, P2): "
  print *, Coeff_req

  !-------------------------------------------------------------------------!
  !                                                                         !
  !                                3-D PES SCRIPT                           !
  !                                                                         !
  !-------------------------------------------------------------------------!


  N_points = 200 ! total number of points in the grid

  ! fixing a point on x-axis over which the potential surface is evaluated

  x_3d = 0.5

  ! Evaluating the potential energy surface over an array of points on y - axis
  y_start = 0.
  y_end   = 2.
  allocate(Y_arr3d(N_points))
  call linspace(y_start, y_end, Y_arr3d)


  ! Evaluating the potential energy surface over an array of points on z - axis
  z_start = 0.
  z_end   = 5.
  allocate(Z_arr3d(N_points))
  call linspace(z_start, z_end, Z_arr3d)

  ! Subroutine to compute the two dimensional potential
  call PES3D(x_3d, Y_arr3d, Z_arr3d, Coeff_D, Coeff_alpha, Coeff_req, N_points)


  !-------------------------------------------------------------------------!
  !                                                                         !
  !                     1-D CUTS FOR SYMMETRY SITES                         !
  !                                                                         !
  !-------------------------------------------------------------------------!
  Nbridge =300
  ! Evaluating the potential energy surface over an array of points on z - axis
  z_start  = 0.5
  z_end    = 10.

  ! Plotting 1D cuts for the symmetry sites:
  allocate(Z_arr1d(Nbridge))
  call linspace(z_start, z_end, Z_arr1d)


  !=========================================================================!
  !                           SYMMETRY SITES                                !
  !=========================================================================!

  ! COORDINATES FOR THE SYMMETRY - SITES
  x   = 0.5
  y   = 0

  ! Subroutine for plotting 1-D cut
  call plot1D( x, y, Nbridge, Coeff_D, Coeff_alpha, Coeff_req, Z_arr1d)


  deallocate(Y_arr3d)
  deallocate(Z_arr3d)
  deallocate(Z_arr1d)


end program main
