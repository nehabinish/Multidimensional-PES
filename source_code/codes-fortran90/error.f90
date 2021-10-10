program error
  !--------------------------------------------------------------------------!
  !                                                                          !
  !                     ALPHA VERSION - N/W(100) system                      !
  !                                                                          !
  !--------------------------------------------------------------------------!
  !  ==========================================================================
  ! This program calculates an estimation of the mean error and mean deviation
  ! for the symmetry sites - top, bridge, hollow, top-bridge, top-hollow and
  ! bridge-hollow.
  ! ===========================================================================

  integer, parameter :: dp = selected_real_kind(p = 15, r = 307) ! Double precision

    !-------------------------------------------------------------------------!
    !                                                                         !
    !                     INITIALISE THE PARAMETERS                           !
    !                                                                         !
    !-------------------------------------------------------------------------!

  real(dp), dimension(:)  , allocatable :: CRP_hollow, CRP_top, CRP_bridge
  real(dp), dimension(:)  , allocatable :: CRP_topbridge, CRP_tophollow,CRP_bridgehollow

  real(dp), dimension(:)  , allocatable :: LEPS_hollow, LEPS_top, LEPS_bridge
  real(dp), dimension(:)  , allocatable :: LEPS_topbridge, LEPS_tophollow, LEPS_bridgehollow

  real(dp), dimension(:)  , allocatable :: Merror_h, Merror_b, Merror_t
  real(dp), dimension(:)  , allocatable :: Merror_tb, Merror_th, Merror_bh

  real(dp), dimension(:)  , allocatable :: RELerror_h, RELerror_b, RELerror_t
  real(dp), dimension(:)  , allocatable :: RELerror_tb, RELerror_th, RELerror_bh

  integer   :: N_hollowCRP, N_topCRP, N_bridgeCRP
  integer   :: N_topbridgeCRP, N_tophollowCRP, N_bridgehollowCRP

  integer   :: N_hollowLEPS, N_topLEPS, N_bridgeLEPS
  integer   :: N_topbridgeLEPS, N_tophollowLEPS, N_bridgehollowLEPS

  integer   :: I


  call read_file('pcrp_h.txt', CRP_hollow, N_hollowCRP)
  call read_file('leps_h.txt', LEPS_hollow, N_hollowLEPS)
  !print*, 'hollow'
  !print*, 'potential CRP',CRP_hollow, 'N CRP', N_hollowCRP
  !print*, 'potential LEPS',LEPS_hollow, 'N LEPS',N_hollowLEPS
  !print*, ' '

  call read_file('pcrp_b.txt', CRP_bridge, N_bridgeCRP)
  call read_file('leps_b.txt', LEPS_bridge, N_bridgeLEPS)
  !print*, 'bridge'
  !print*, 'potential CRP',CRP_bridge,'N CRP',  N_bridgeCRP
  !print*, 'potential LEPS',LEPS_bridge, 'N LEPS',N_bridgeLEPS
  !print*, ' '

  call read_file('pcrp_t.txt', CRP_top, N_topCRP)
  call read_file('leps_t.txt', LEPS_top, N_topLEPS)
  !print*, 'top'
  !print*,'potential CRP', CRP_top, 'N CRP', N_topCRP
  !print*, 'potential LEPS',LEPS_top,'N LEPS', N_topLEPS
  !print*, ' '

  call read_file('pcrp_tbridge.txt', CRP_topbridge, N_topbridgeCRP)
  call read_file('leps_tbridge.txt', LEPS_topbridge, N_topbridgeLEPS)
  !print*, 'top-bridge'
  !print*,'potential CRP', CRP_topbridge, 'N CRP', N_topbridgeCRP
  !print*, 'potential LEPS',LEPS_topbridge,'N LEPS', N_topbridgeLEPS
  !print*, ' '

  call read_file('pcrp_thollow.txt', CRP_tophollow, N_tophollowCRP)
  call read_file('leps_thollow.txt', LEPS_tophollow, N_tophollowLEPS)
  !print*, 'top-hollow'
  !print*, 'potential CRP',CRP_tophollow, 'N CRP', N_tophollowCRP
  !print*,'potential LEPS', LEPS_tophollow, 'N LEPS',N_tophollowLEPS
  !print*, ' '

  call read_file('pcrp_bridgehollow.txt', CRP_bridgehollow, N_bridgehollowCRP)
  call read_file('leps_bridgehollow.txt', LEPS_bridgehollow, N_bridgehollowLEPS)
  !print*, 'bridge hollow'
  !print*, 'potential CRP', CRP_bridgehollow, 'N CRP', N_bridgehollowCRP
  !print*, 'potential LEPS',LEPS_bridgehollow,'N LEPS', N_bridgehollowLEPS
  !print*, ' '

  !-------------------------------------------------------------------------!
  !                                                                         !
  !                    CALCULATE THE MEAN ABSOLOUTE ERROR                   !
  !                                                                         !
  !-------------------------------------------------------------------------!

  !=========================================================================!
  !                     High Symmetry Site - HOLLOW                         !
  !=========================================================================!

  allocate(Merror_h(N_hollowCRP))
  allocate(RELerror_h(N_hollowCRP))

  Do I = 1, N_hollowCRP

    Merror_h(I)   = abs(LEPS_hollow(I) - CRP_hollow(I))
    RELerror_h(I) = abs(((CRP_hollow(I) - LEPS_hollow(I))/ CRP_hollow(I))*100)

  End Do

  print *, ' '
  print *, ' Mean absolute error for High Symmetry Site: Hollow - ', (sum(Merror_h)/N_hollowCRP)
  print *, ' Mean relative error for High Symmetry Site: Hollow - ', (sum(RELerror_h)/N_hollowCRP),'%'

  deallocate(Merror_h)
  deallocate(RELerror_h)

  !=========================================================================!
  !                     High Symmetry Site - BRIDGE                         !
  !=========================================================================!

  allocate(Merror_b(N_bridgeCRP))
  allocate(RELerror_b(N_bridgeCRP))

  Do I = 1, N_bridgeCRP

    Merror_b(I)   = abs( LEPS_bridge(I) - CRP_bridge(I))
    RELerror_b(I) = abs(((CRP_bridge(I) - LEPS_bridge(I)) / CRP_bridge(I))*100)

  End Do

  print *, ' '
  print *, ' Mean absolute error for High Symmetry Site: Bridge - ', (sum(Merror_b)/N_bridgeCRP)
  print *, ' Mean relative error for High Symmetry Site: Bridge - ', (sum(RELerror_b)/N_bridgeCRP),'%'

  deallocate(Merror_b)
  deallocate(RELerror_b)

  !=========================================================================!
  !                     High Symmetry Site - TOP                            !
  !=========================================================================!

  allocate(Merror_t(N_topCRP))
  allocate(RELerror_t(N_topCRP))

  Do I = 1, N_topCRP

    Merror_t(I) = abs( LEPS_top(I) - CRP_top(I))
    RELerror_t(I) = abs(((CRP_top(I) - LEPS_top(I)) / CRP_top(I))*100)

  End Do

  print *, ' '
  print *, ' Mean absolute error for High Symmetry Site: Top - ', (sum(Merror_t)/N_topCRP)
  print *, ' Mean relative error for High Symmetry Site: Top - ', (sum(RELerror_t)/N_topCRP),'%'

  deallocate(Merror_t)
  deallocate(RELerror_t)

  !=========================================================================!
  !                  Low Symmetry Site - TOP-BRIDGE                         !
  !=========================================================================!

  allocate(Merror_tb(N_topbridgeCRP))
  allocate(RELerror_tb(N_topbridgeCRP))

  Do I = 1, N_topbridgeCRP

    Merror_tb(I) = abs( LEPS_topbridge(I) - CRP_topbridge(I))
    RELerror_tb(I) = abs(((CRP_topbridge(I) - LEPS_topbridge(I)) / CRP_topbridge(I))*100)

  End Do

  print *, ' '
  print *, ' Mean absolute error for Low Symmetry Site: Top Bridge - ', (sum(Merror_tb)/N_topbridgeCRP)
  print *, ' Mean relative error for Low Symmetry Site: Top Bridge - ', (sum(RELerror_tb)/N_topbridgeCRP),'%'

  deallocate(Merror_tb)
  deallocate(RELerror_tb)

  !=========================================================================!
  !                  Low Symmetry Site - TOP-HOLLOW                         !
  !=========================================================================!

  allocate(Merror_th(N_tophollowCRP))
  allocate(RELerror_th(N_tophollowCRP))

  Do I = 1, N_tophollowCRP

    Merror_th(I) = abs( LEPS_tophollow(I) - CRP_tophollow(I))
    RELerror_th(I) = abs(((CRP_tophollow(I) - LEPS_tophollow(I)) / CRP_tophollow(I))*100)

  End Do

  print *, ' '
  print *, ' Mean absolute error for Low Symmetry Site: Top Hollow - ', (sum(Merror_th)/N_tophollowCRP)
  print *, ' Mean relative error for Low Symmetry Site: Top Hollow- ', (sum(RELerror_th)/N_tophollowCRP),'%'

  deallocate(Merror_th)
  deallocate(RELerror_th)

  !=========================================================================!
  !                  Low Symmetry Site - BRIDGE-HOLLOW                         !
  !=========================================================================!

  allocate(Merror_bh(N_bridgehollowCRP))
  allocate(RELerror_bh(N_bridgehollowCRP))

  Do I = 1, N_bridgehollowCRP

    Merror_bh(I) = abs( LEPS_bridgehollow(I) - CRP_bridgehollow(I))
    RELerror_bh(I) = abs(((CRP_bridgehollow(I) - LEPS_bridgehollow(I)) / CRP_bridgehollow(I))*100)

  End Do

  print *, ' '
  print *, ' Mean absolute error for Low Symmetry Site: Bridge Hollow - ', (sum(Merror_bh)/N_bridgehollowCRP)
  print *, ' Mean relative error for Low Symmetry Site: Bridge Hollow - ', (sum(RELerror_bh)/N_bridgehollowCRP),'%'
  print *, ' '

  deallocate(Merror_bh)
  deallocate(RELerror_bh)

  contains

  subroutine read_file(filename, value, count)
    !-------------------------------------------------------------------------!
    !                                                                         !
    !                             SUBROUTINE TO READ FILES                    !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments                                                               !
    ! ===========                                                             !
    !  filename    (input) CHARACTER                                          !
    !              Input file to be read                                      !
    !                                                                         !
    !  value       (output) REAL - (kind = 8)                                 !
    !              Value to be read from the file                             !
    !                                                                         !
    !  count       (output) REAL - (kind = 8)                                 !
    !              Total number of lines in the file                          !
    ! ======================================================================= !

    implicit none

    character(len=*), intent(in) :: filename
    integer:: io,i
    integer, intent(out) :: count
    real(kind=8), intent(out), dimension(:), allocatable:: value

    open(unit=102, file = filename)
    count = 0
    Do
      READ(102,*,iostat=io)
      IF (io/=0) EXIT

      count = count + 1
    End Do

    rewind(102)

    allocate(value(count))

    read(102,*) value

    return
    deallocate(value)

    close(102)

  end subroutine

end program error
