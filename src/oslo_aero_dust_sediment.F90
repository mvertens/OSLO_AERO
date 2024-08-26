module oslo_aero_dust_sediment

  !---------------------------------------------------------------------------------
  ! Routines to compute tendencies from sedimentation of dust
  ! Author: Phil Rasch
  !---------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp
  use physconst,      only: gravit, rair
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none
  private

  ! public routines
  public :: oslo_aero_dust_sediment_vel
  public :: oslo_aero_dust_sediment_tend

  ! private routines
  private :: getflx
  private :: cfint2
  private :: cfdotmc_pro

  real (r8), parameter :: vland  = 2.8_r8            ! dust fall velocity over land  (cm/s)
  real (r8), parameter :: vocean = 1.5_r8            ! dust fall velocity over ocean (cm/s)
  real (r8), parameter :: mxsedfac = 0.99_r8       ! maximum sedimentation flux factor

!===============================================================================
contains
!===============================================================================

  subroutine oslo_aero_dust_sediment_vel(ncol, icefrac, landfrac, ocnfrac, pmid, pdel, t, dustmr, pvdust)

    ! Compute gravitational sedimentation velocities for dust
    ! note that pvel is at the interfaces (loss from cell is based on pvel(ilev+1))

    ! Arguments
    integer,  intent(in)  :: ncol                 ! number of colums to process
    real(r8), intent(in)  :: icefrac (pcols)      ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)      ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)      ! ocean fraction (fraction)
    real(r8), intent(in)  :: pmid(pcols,pver)     ! pressure of midpoint levels (Pa)
    real(r8), intent(in)  :: pdel(pcols,pver)     ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t(pcols,pver)        ! temperature (K)
    real(r8), intent(in)  :: dustmr(pcols,pver)   ! dust (kg/kg)
    real(r8), intent(out) :: pvdust (pcols,pverp) ! vertical velocity of dust (Pa/s)

    ! Local variables
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real (r8) :: vfall(pcols)                       ! settling velocity of dust particles (m/s)
    integer   :: icol,ilev
    real (r8) :: lbound, ac, bc, cc

    ! dust fall velocity
    do ilev = 1,pver
       do icol = 1,ncol
          ! merge the dust fall velocities for land and ocean (cm/s) SHOULD ALSO ACCOUNT FOR ICEFRAC
          vfall(icol) = vland*landfrac(icol) + vocean*(1._r8-landfrac(icol))

          ! fall velocity (assume positive downward)
          pvdust(icol,ilev+1) = vfall(icol)
       end do
    end do
  end subroutine oslo_aero_dust_sediment_vel

  !===============================================================================
  subroutine oslo_aero_dust_sediment_tend ( ncol,   dtime,  pint, pmid, pdel, t, &
       dustmr, pvdust, dusttend, sfdust)

    !----------------------------------------------------------------------
    !  Apply Particle Gravitational Sedimentation
    ! -> note that pvel is at the interfaces (loss from cell is based on pvel(ilev+1))
    !----------------------------------------------------------------------

    ! Arguments
    integer,  intent(in)  :: ncol                      ! number of colums to process
    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), intent(in)  :: pint(pcols,pverp)         ! interfaces pressure (Pa)
    real(r8), intent(in)  :: pmid(pcols,pver)          ! midpoint pressures (Pa)
    real(r8), intent(in)  :: pdel(pcols,pver)          ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t(pcols,pver)             ! temperature (
    real(r8), intent(in)  :: dustmr(pcols,pver)        ! dust (kg/kg)
    real(r8), intent(in)  :: pvdust (pcols,pverp)      ! vertical velocity of dust drops  (Pa/s)
    real(r8), intent(out) :: dusttend(pcols,pver)      ! dust tend
    real(r8), intent(out) :: sfdust(pcols)             ! surface flux of dust (rain, kg/m/s)

    ! Local variables
    integer  :: icol,ilev
    real(r8) :: fxdust(pcols,pverp)                     ! fluxes at the interfaces, dust (positive = down)
    !----------------------------------------------------------------------

    ! initialize variables
    fxdust  (:ncol,:) = 0._r8 ! flux at interfaces (dust)
    dusttend(:ncol,:) = 0._r8 ! tend (dust)
    sfdust(:ncol)     = 0._r8 ! sedimentation flux out bot of column (dust)

    ! fluxes at interior points
    call getflx(ncol, pint, dustmr, pvdust, dtime, fxdust)

    ! calculate fluxes at boundaries
    do icol = 1,ncol
       fxdust(icol,1) = 0
       ! surface flux by upstream scheme
       fxdust(icol,pverp) = dustmr(icol,pver) * pvdust(icol,pverp) * dtime
    end do

    ! filter out any negative fluxes from the getflx routine
    do ilev = 2,pver
       fxdust(:ncol,ilev) = max(0._r8, fxdust(:ncol,ilev))
    end do

    ! Limit the flux out of the bottom of each cell to the water content in each phase.
    ! Apply mxsedfac to prevent generating very small negative cloud water/ice
    ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
    ! ***Should we include the flux in the top, to allow for thin surface layers?
    ! ***Requires simple treatment of cloud overlap, already included below.
    do ilev = 1,pver
       do icol = 1,ncol
          fxdust(icol,ilev+1) = min( fxdust(icol,ilev+1), mxsedfac * dustmr(icol,ilev) * pdel(icol,ilev) )
       end do
    end do

    ! Now calculate the tendencies
    do ilev = 1,pver
       do icol = 1,ncol
          ! net flux into cloud changes cloud dust/ice (all flux is out of cloud)
          dusttend(icol,ilev)  = (fxdust(icol,ilev) - fxdust(icol,ilev+1)) / (dtime * pdel(icol,ilev))
       end do
    end do

    ! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfdust(:ncol) = fxdust(:ncol,pverp) / (dtime*gravit)

  end subroutine oslo_aero_dust_sediment_tend

  !===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)

    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......

    ! arguments
    integer , intent(in)  :: ncol  ! number of colums to process
    real(r8), intent(out) :: flux(pcols,pverp)
    real(r8), intent(in)  :: xw(pcols,pverp)
    real(r8), intent(in)  :: vel(pcols,pverp)
    real(r8), intent(in)  :: deltat

    ! local variables
    integer   :: icol
    integer   :: ilev
    real (r8) :: psi(pcols,pverp)
    real (r8) :: phi(pcols,pverp-1)
    real (r8) :: fdot(pcols,pverp)
    real (r8) :: xx(pcols)
    real (r8) :: fxdot(pcols)
    real (r8) :: fxdd(pcols)
    real (r8) :: psistar(pcols)
    real (r8) :: xxk(pcols,pver)

    do icol = 1,ncol
       ! integral of phi
       psi(icol,1) = 0._r8
       ! fluxes at boundaries
       flux(icol,1) = 0
       flux(icol,pverp) = 0._r8
    end do

    ! integral function
    do ilev = 2,pverp
       do icol = 1,ncol
          psi(icol,ilev) = phi(icol,ilev-1)*(xw(icol,ilev)-xw(icol,ilev-1)) + psi(icol,ilev-1)
       end do
    end do

    ! calculate the derivatives for the interpolating polynomial
    call cfdotmc_pro (ncol, xw, psi, fdot)

    ! calculate fluxes at interior pts
    do ilev = 2,pver
       do icol = 1,ncol
          xxk(icol,ilev) = xw(icol,ilev)-vel(icol,ilev)*deltat
       end do
    end do
    do ilev = 2,pver
       call cfint2(ncol, xw, psi, fdot, xxk(1,ilev), fxdot, fxdd, psistar)
       do icol = 1,ncol
          flux(icol,ilev) = (psi(icol,ilev)-psistar(icol))
       end do
    end do

  end subroutine getflx

  !===============================================================================
  subroutine cfint2 (ncol, x, f, fdot, xin, fxdot, fxdd, psistar)

    ! arguments
    integer   , intent(in)  :: ncol                      ! number of colums to process
    real (r8) , intent(in)  :: x(pcols, pverp)
    real (r8) , intent(in)  :: f(pcols, pverp)
    real (r8) , intent(inout) :: fdot(pcols, pverp)
    real (r8) , intent(in)  :: xin(pcols)
    real (r8) , intent(out) :: fxdot(pcols)
    real (r8) , intent(out) :: fxdd(pcols)
    real (r8) , intent(out) :: psistar(pcols)

    ! local variables
    integer   :: icol
    integer   :: ilev
    integer   :: intz(pcols)
    real (r8) :: dx
    real (r8) :: s
    real (r8) :: c2
    real (r8) :: c3
    real (r8) :: xx
    real (r8) :: xinf
    real (r8) :: psi1, psi2, psi3, psim
    real (r8) :: cfint
    real (r8) :: cfnew
    real (r8) :: xins(pcols)
    real (r8) :: a, b, c ! the minmod function
    real (r8) :: minmod ! the minmod function
    real (r8) :: medan ! the minmod function

    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do icol = 1,ncol
       xins(icol) = medan(x(icol,1), xin(icol), x(icol,pverp))
       intz(icol) = 0
    end do

    ! first find the interval
    do ilev =  1,pverp-1
       do icol = 1,ncol
          if ((xins(icol)-x(icol,ilev))*(x(icol,ilev+1)-xins(icol)).ge.0._r8) then
             intz(icol) = ilev
          endif
       end do
    end do

    do icol = 1,ncol
       if (intz(icol).eq.0) then
          write(iulog,*) ' interval was not found for col icol ', icol
          call endrun('DUST_SEDIMENT_MOD:cfint2 -- interval was not found ')
       endif
    end do

    ! now interpolate
    do icol = 1,ncol
       ilev = intz(icol)
       dx = (x(icol,ilev+1)-x(icol,ilev))
       s = (f(icol,ilev+1)-f(icol,ilev))/dx
       c2 = (3*s-2*fdot(icol,ilev)-fdot(icol,ilev+1))/dx
       c3 = (fdot(icol,ilev)+fdot(icol,ilev+1)-2*s)/dx**2
       xx = (xins(icol)-x(icol,ilev))
       fxdot(icol) =  (3*c3*xx + 2*c2)*xx + fdot(icol,ilev)
       fxdd(icol) = 6*c3*xx + 2*c2
       cfint = ((c3*xx + c2)*xx + fdot(icol,ilev))*xx + f(icol,ilev)

       ! limit the interpolant
       psi1 = f(icol,ilev)+(f(icol,ilev+1)-f(icol,ilev))*xx/dx
       if (ilev.eq.1) then
          psi2 = f(icol,1)
       else
          psi2 = f(icol,ilev) + (f(icol,ilev)-f(icol,ilev-1))*xx/(x(icol,ilev)-x(icol,ilev-1))
       endif
       if (ilev+1.eq.pverp) then
          psi3 = f(icol,pverp)
       else
          psi3 = f(icol,ilev+1) - (f(icol,ilev+2)-f(icol,ilev+1))*(dx-xx)/(x(icol,ilev+2)-x(icol,ilev+1))
       endif
       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)
       if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
       endif
       psistar(icol) = cfnew
    end do

  end subroutine cfint2

  !===============================================================================
  subroutine cfdotmc_pro (ncol, x, f, fdot)

    ! prototype version; eventually replace with final SPITFIRE scheme
    ! calculate the derivative for the interpolating polynomial multi column version
    ! assumed variable distribution

    !     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
    !     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
    !     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
    !     .........d2.......d3.......d4.......d5.........     2,pver points
    !     .........s2.......s3.......s4.......s5.........     2,pver points
    !     .............dh2......dh3......dh4.............     2,pver-1 points
    !     .............eh2......eh3......eh4.............     2,pver-1 points
    !     ..................e3.......e4..................     3,pver-1 points
    !     .................ppl3......ppl4................     3,pver-1 points
    !     .................ppr3......ppr4................     3,pver-1 points
    !     .................t3........t4..................     3,pver-1 points
    !     ................fdot3.....fdot4................     3,pver-1 points


    ! arguments
    integer   , intent(in)  :: ncol               ! number of colums to process
    real (r8) , intent(in)  :: x(pcols, pverp)
    real (r8) , intent(in)  :: f(pcols, pverp)
    real (r8) , intent(out) :: fdot(pcols, pverp) ! derivative at nodes

    ! local variables
    integer  :: icol,ilev
    real(r8) :: a,b,c            ! work vars
    real(r8) :: s(pcols,pverp)   ! first divided differences at nodes
    real(r8) :: sh(pcols,pverp)  ! first divided differences between nodes
    real(r8) :: d(pcols,pverp)   ! second divided differences at nodes
    real(r8) :: dh(pcols,pverp)  ! second divided differences between nodes
    real(r8) :: e(pcols,pverp)   ! third divided differences at nodes
    real(r8) :: eh(pcols,pverp)  ! third divided differences between nodes
    real(r8) :: pp               ! p prime
    real(r8) :: ppl(pcols,pverp) ! p prime on left
    real(r8) :: ppr(pcols,pverp) ! p prime on right
    real(r8) :: qpl
    real(r8) :: qpr
    real(r8) :: ttt
    real(r8) :: t
    real(r8) :: tmin
    real(r8) :: tmax
    real(r8) :: delxh(pcols,pverp)
    real(r8) :: minmod           ! the minmod function
    real(r8) :: medan            ! the minmod function

    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do ilev = 1,pver
       ! first divided differences between nodes
       do icol = 1, ncol
          delxh(icol,ilev) = (x(icol,ilev+1)-x(icol,ilev))
          sh(icol,ilev) = (f(icol,ilev+1)-f(icol,ilev))/delxh(icol,ilev)
       end do

       ! first and second divided differences at nodes
       if (ilev.ge.2) then
          do icol = 1,ncol
             d(icol,ilev) = (sh(icol,ilev)-sh(icol,ilev-1))/(x(icol,ilev+1)-x(icol,ilev-1))
             s(icol,ilev) = minmod(sh(icol,ilev),sh(icol,ilev-1))
          end do
       endif
    end do

    ! second and third divided diffs between nodes
    do ilev = 2,pver-1
       do icol = 1, ncol
          eh(icol,ilev) = (d(icol,ilev+1)-d(icol,ilev))/(x(icol,ilev+2)-x(icol,ilev-1))
          dh(icol,ilev) = minmod(d(icol,ilev),d(icol,ilev+1))
       end do
    end do

    ! treat the boundaries
    do icol = 1,ncol
       e(icol,2) = eh(icol,2)
       e(icol,pver) = eh(icol,pver-1)
       !  outside level
       fdot(icol,1) = sh(icol,1) - d(icol,2)*delxh(icol,1) - eh(icol,2)*delxh(icol,1)*(x(icol,1)-x(icol,3))
       fdot(icol,1) = minmod(fdot(icol,1),3*sh(icol,1))
       fdot(icol,pverp) = sh(icol,pver) + d(icol,pver)*delxh(icol,pver) &
            + eh(icol,pver-1)*delxh(icol,pver)*(x(icol,pverp)-x(icol,pver-1))
       fdot(icol,pverp) = minmod(fdot(icol,pverp),3*sh(icol,pver))

       ! one in from boundary
       fdot(icol,2) = sh(icol,1) + d(icol,2)*delxh(icol,1) - eh(icol,2)*delxh(icol,1)*delxh(icol,2)
       fdot(icol,2) = minmod(fdot(icol,2),3*s(icol,2))
       fdot(icol,pver) = sh(icol,pver) - d(icol,pver)*delxh(icol,pver) &
            - eh(icol,pver-1)*delxh(icol,pver)*delxh(icol,pver-1)
       fdot(icol,pver) = minmod(fdot(icol,pver),3*s(icol,pver))
    end do

    do ilev = 3,pver-1
       do icol = 1,ncol
          e(icol,ilev) = minmod(eh(icol,ilev),eh(icol,ilev-1))
       end do
    end do

    do ilev = 3,pver-1
       do icol = 1,ncol
          ! p prime at ilev-0.5
          ppl(icol,ilev)=sh(icol,ilev-1) + dh(icol,ilev-1)*delxh(icol,ilev-1)

          ! p prime at ilev+0.5
          ppr(icol,ilev)=sh(icol,ilev)   - dh(icol,ilev)  *delxh(icol,ilev)
          t = minmod(ppl(icol,ilev),ppr(icol,ilev))

          ! derivate from parabola thru f(icol,ilev-1), f(icol,ilev), and f(icol,ilev+1)
          pp = sh(icol,ilev-1) + d(icol,ilev)*delxh(icol,ilev-1)

          ! quartic estimate of fdot
          fdot(icol,ilev) = pp - delxh(icol,ilev-1)*delxh(icol,ilev)*(eh(icol,ilev-1)*(x(icol,ilev+2)-x(icol,ilev)) &
               + eh(icol,ilev  )*(x(icol,ilev  )-x(icol,ilev-2)))/(x(icol,ilev+2)-x(icol,ilev-2))

          ! now limit it
          qpl = sh(icol,ilev-1) + delxh(icol,ilev-1)*minmod(d(icol,ilev-1)+ e(icol,ilev-1)*(x(icol,ilev)-x(icol,ilev-2)), &
               d(icol,ilev) - e(icol,ilev)*delxh(icol,ilev))
          qpr = sh(icol,ilev) + delxh(icol,ilev  )*minmod(d(icol,ilev) + e(icol,ilev)*delxh(icol,ilev-1), &
               d(icol,ilev+1)+e(icol,ilev+1)*(x(icol,ilev)-x(icol,ilev+2)))

          fdot(icol,ilev) = medan(fdot(icol,ilev), qpl, qpr)

          ttt = minmod(qpl, qpr)
          tmin = min(0._r8,3*s(icol,ilev),1.5_r8*t,ttt)
          tmax = max(0._r8,3*s(icol,ilev),1.5_r8*t,ttt)
          fdot(icol,ilev) = fdot(icol,ilev) + minmod(tmin-fdot(icol,ilev), tmax-fdot(icol,ilev))
       end do
    end do

  end subroutine cfdotmc_pro

end module oslo_aero_dust_sediment
