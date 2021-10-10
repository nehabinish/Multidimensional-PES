program Global_min
  !--------------------------------------------------------------------------!
  !                                                                          !
  !                     ALPHA VERSION - N/W(100) system                      !
  !                                                                          !
  !--------------------------------------------------------------------------!
  !  ==========================================================================
  !  This program finds the global minimum of all the LEPS-PES and the CRP-PES
  ! ===========================================================================

  integer, parameter :: dp = selected_real_kind(p = 15, r = 307) ! Double precision

    !-------------------------------------------------------------------------!
    !                                                                         !
    !                     INITIALISE THE PARAMETERS                           !
    !                                                                         !
    !-------------------------------------------------------------------------!

  real(dp), dimension(:)  , allocatable :: CRP_hollow, CRP_top, CRP_bridge
  real(dp), dimension(:)  , allocatable :: CRP_topbridge, CRP_tophollow,CRP_bridgehollow
  real(dp), dimension(6)                :: CRPmin, LEPSmin
  real(dp), dimension(:)  , allocatable :: LEPS_hollow, LEPS_top, LEPS_bridge
  real(dp), dimension(:)  , allocatable :: LEPS_topbridge, LEPS_tophollow, LEPS_bridgehollow

  real(dp)  :: CRPmin_h, CRPmin_b, CRPmin_t
  real(dp)  :: CRPmin_th, CRPmin_tb, CRPmin_bh
  real(dp)  :: LEPSmin_h, LEPSmin_b, LEPSmin_t
  real(dp)  :: LEPSmin_th, LEPSmin_tb, LEPSmin_bh

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

  CRPmin_h   = minval(CRP_hollow)
  CRPmin(1)  = CRPmin_h
  LEPSmin_h  = minval(LEPS_hollow)
  LEPSmin(1) = LEPSmin_h

  CRPmin_b   = minval(CRP_bridge)
  CRPmin(2)  = CRPmin_b
  LEPSmin_b  = minval(LEPS_bridge)
  LEPSmin(2) = LEPSmin_b

  CRPmin_t   = minval(CRP_top)
  CRPmin(3)  = CRPmin_t
  LEPSmin_t  = minval(LEPS_top)
  LEPSmin(3) = LEPSmin_t

  CRPmin_th  = minval(CRP_tophollow)
  CRPmin(4)  = CRPmin_th
  LEPSmin_th = minval(LEPS_tophollow)
  LEPSmin(4) = LEPSmin_th

  CRPmin_tb  = minval(CRP_topbridge)
  CRPmin(5)  = CRPmin_tb
  LEPSmin_tb = minval(LEPS_topbridge)
  LEPSmin(5) = LEPSmin_tb

  CRPmin_bh  = minval(CRP_bridgehollow)
  CRPmin(6)  = CRPmin_bh
  LEPSmin_bh = minval(LEPS_bridgehollow)
  LEPSmin(6) = LEPSmin_bh

  Globalmin_CRP  = minval(CRPmin)
  Globalmin_LEPS = minval(LEPSmin)

  print *, ' '
  print *, 'The Global minimum for the CRP  - Potential Energy Surfaces : ', Globalmin_CRP, 'eV'
  print *, ' '
  print *, 'The Global minimum for the LEPS - Potential Energy Surfaces : ', Globalmin_LEPS, 'eV'
  print *, ' '

contains

    subroutine read_file(filename, value, count)

    !-----------------------------------------------------------------------!
    !                                                                       !
    !                             SUBROUTINE TO READ FILES                  !
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Arguments                                                             !
    ! ===========                                                           !
    !  filename    (input) CHARACTER                                        !
    !              Input file to be read                                    !
    !                                                                       !
    !  value       (output) REAL - (kind = 8)  - allocatable array          !
    !              Value to be read from the file                           !
    !                                                                       !
    !  count       (output) REAL - (kind = 8)                               !
    !              Total number of lines in the file                        !
    !                                                                       !
    ! ===================================================================== !

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


end program Global_min
