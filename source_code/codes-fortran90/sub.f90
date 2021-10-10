module plots

implicit none

contains
  subroutine linspace(from, to, array)

  !--------------------------------------------------------------------------!
  !                                                                          !
  !   SUBROUTINE TO MAKE AN ARRAY OF EQUALLY SPACED POINTS B/W TWO POINTS    !
  !  ======================================================================  !
  !                                                                          !
  !--------------------------------------------------------------------------!
  ! Arguments
  ! ===========
  !  from     (input) REAL - double precision
  !           Starting point for the grid
  !
  !  to       (input) REAL - double precision
  !           Ending point for the grid
  !
  !  array   (output) REAL - double precision
  !           Equally spaced array in the given range
  ! ===========================================================================

    integer , parameter   :: dp = selected_real_kind(p = 15, r = 307) ! Double precision
    integer               :: n, i
    real(dp), intent(in)  :: from, to
    real(dp), intent(out) :: array(:)
    real(dp)              :: range


    n     = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
      array(1) = from
      return
    end if


    do i=1, n
      array(i) = from + range * (i - 1) / (n - 1)
    end do

  end subroutine linspace

  subroutine PES3D(X, Y, Z, Dcoeff, Acoeff, Rcoeff, N)

  !---------------------------------------------------------------------------!
  !                                                                           !
  !                   SUBROUTINE TO MAKE THE 3-D PES SCRIPT                   !
  !              ==================================================           !                                                          !
  !                                                                           !
  !---------------------------------------------------------------------------!
  ! Arguments
  ! ===========
  !  N       (input) INTEGER
  !           Total number of points in the grid.
  !
  !  X       (input) REAL - double precision
  !          The fixed point on the x-axis over which the potential energy
  !          surfaces are evaluated
  !
  !  Y       (input) REAL - double precision - Array of dimension N
  !          Array of values on the y-axis over which the potential energy
  !          surfaces are evaluated
  !
  !  Z       (input) REAL - double precision - Array of dimension N
  !          Array of values on the z-axis over which the potential energy
  !          surfaces are evaluated
  !
  ! Dcoeff   (input) t
  !           Fourier Coefficients of the Morse parameter D
  !
  ! Acoeff   (input) REAl - double precesion - Dimension(3,1)
  !           Fourier Coefficients of the Morse parameter Alpha
  !
  ! Rcoeff   (input) REAl - double precesion - Dimension(3,1)
  !           Fourier Coefficients of the Morse parameter r_eq
  ! ===========================================================================

    integer , parameter                   :: dp = selected_real_kind(p = 15, r = 307) ! Double precision
    integer , intent(in)                  :: N
    real(dp), dimension(N,N)              :: potential3d
    real(dp), intent(in)                  :: Dcoeff(3,1), Acoeff(3,1), Rcoeff(3,1)
    real(dp), dimension(N)                :: D, alpha, req
    real(dp)                            :: pi, delta
    real(dp), intent(in)                  :: X
    real(dp), dimension(N), intent(in)    :: Y,Z
    integer                               :: i,j

    pi    = 3.141592653589793238462643
    delta = 3.17*10._8**(-10)


    open(unit = 50, file = 'potential3d.txt')
    open(unit = 60, file = 'y_arr3d.txt')
    open(unit = 70, file = 'z_arr3d.txt')

    Do i = 1,N
      Do j = 1,N

        D(i)     = Dcoeff(1,1) + Dcoeff(2,1)*(cos(2*pi*x/delta)  + cos(2*pi*y(i)/delta)) &
                         + Dcoeff(3,1)*(cos(2*pi*(x+y(i))/delta) + cos(2*pi*(x-y(i))/delta))

        alpha(i) = Acoeff(1,1) + Acoeff(2,1)*(cos(2*pi*x/delta)  + cos(2*pi*y(i)/delta)) &
                         + Acoeff(3,1)*(cos(2*pi*(x+y(i))/delta) + cos(2*pi*(x-y(i))/delta))

        req(i)   = Rcoeff(1,1) + Rcoeff(2,1)*(cos(2*pi*x/delta)   + cos(2*pi*y(i)/delta)) &
                          + Rcoeff(3,1)*(cos(2*pi*(x+y(i))/delta) + cos(2*pi*(x-y(i))/delta))

        potential3d(i,j) = D(i)*(exp(-2*alpha(i)*(Z(j) - req(i))) - (2 * exp(- alpha(i) *(Z(j) - req(i)))))

      end do

      write(60,*) Y(i)
      write(70,*) Z(i)

    end do

    Do i = 1,N
      Do j = 1,N
          write(50,*) potential3d(i,j)
      end do
      write(50,*) ' '
    end do



    close(50)
    close(60)
    close(70)

  end subroutine PES3D


  subroutine plot1D(x, y, N, Dcoeff, Acoeff, Rcoeff, Z)

  !---------------------------------------------------------------------------!
  !                                                                           !
  !            SUBROUTINE TO MAKE 1-D CUTS FOR SYMMETRY SITES                 !
  !                                                                           !
  !---------------------------------------------------------------------------!
  ! Arguments
  ! ===========
  !  N       (input) INTEGER
  !           Total number of points in the grid.
  !
  !  x       (input) REAL - double precision
  !          The fixed point on the x-axis over specified by the symmetry site
  !
  !  y       (input) REAL - double precision
  !          The fixed point on the y-axis over specified by the symmetry site
  !
  !  Z       (input) REAL - double precision - Array of dimension N
  !          Array of values on the z-axis over which the potential energy
  !          surfaces are evaluated
  !
  ! Dcoeff   (input) REAl - double precesion - Dimension(3,1)
  !           Fourier Coefficients of the Morse parameter D
  !
  ! Acoeff   (input) REAl - double precesion - Dimension(3,1)
  !           Fourier Coefficients of the Morse parameter Alpha
  !
  ! Rcoeff   (input) REAl - double precesion - Dimension(3,1)
  !           Fourier Coefficients of the Morse parameter r_eq
  ! ===========================================================================

    integer , parameter                 :: dp = selected_real_kind(p = 15, r = 307) ! Double precision
    integer , intent(in)                :: N
    real(dp), intent(in)                :: x,y
    real(dp), intent(in)                :: Dcoeff(3,1), Acoeff(3,1), Rcoeff(3,1)
    real(dp), dimension(N)              :: D, alpha, req, V
    real(dp)                            :: pi, delta
    real(dp), dimension(N), intent(in)  :: Z
    integer                             :: i

    pi    = 3.141592653589793238462643
    delta = 3.17*10._8**(-10)

    open(unit=650, file = 'potential.txt')
    open(unit=660, file = 'distance.txt')


    Do i = 1, N

      D(i)     = Dcoeff(1,1) + Dcoeff(2,1)*(cos(2*pi*x/delta) + cos(2*pi*y/delta)) &
                       + Dcoeff(3,1)*(cos(2*pi*(x+y)/delta) + cos(2*pi*(x-y)/delta))

      alpha(i) = Acoeff(1,1) + Acoeff(2,1)*(cos(2*pi*x/delta) + cos(2*pi*y/delta)) &
                       + Acoeff(3,1)*(cos(2*pi*(x+y)/delta) + cos(2*pi*(x-y)/delta))

      req(i)   = Rcoeff(1,1) + Rcoeff(2,1)*(cos(2*pi*x/delta) + cos(2*pi*y/delta)) &
                        + Rcoeff(3,1)*(cos(2*pi*(x+y)/delta) + cos(2*pi*(x-y)/delta))

      V(i)     = D(i)*(exp(-2*alpha(i)*(Z(i) - req(i))) - (2 * exp(- alpha(i) *(Z(i) - req(i)))))


      write(650,*) V(i)
      write(660,*) Z(i)

    End do

    close(650)
    close(660)

    end subroutine plot1D

end module plots
