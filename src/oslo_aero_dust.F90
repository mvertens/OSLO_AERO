module oslo_aero_dust

  ! Calculate emission of all dusts.
  ! Note that the mobilization is calculated in the land model and
  ! the soil erodibility factor is applied here.

  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: mpicom, mstrid=>masterprocid, masterproc
  use spmd_utils,       only: mpi_logical, mpi_real8, mpi_character, mpi_integer,  mpi_success
  use namelist_utils,   only: find_group_name
  use ppgrid,           only: pcols
  use constituents,     only: cnst_name, pcnst
  use cam_logfile,      only: iulog
  use cam_abortutils,   only: endrun
  use shr_dust_emis_mod,only: is_dust_emis_zender, is_zender_soil_erod_from_atm, shr_dust_emis_readnl
  use soil_erod_mod,    only: soil_erod_init
  !
  use oslo_aero_share,  only: l_dst_a2, l_dst_a3

  implicit none
  private

  ! public routines
  public :: oslo_aero_dust_readnl
  public :: oslo_aero_dust_init
  public :: oslo_aero_dust_emis

  character(len=6), public :: dust_names(10)

  integer , parameter :: numberOfDustModes = 2
  real(r8), parameter :: emis_fraction_in_mode(numberOfDustModes) = (/0.13_r8, 0.87_r8 /)
  integer             :: tracerMap(numberOfDustModes) = (/-99, -99/) !index of dust tracers in the modes

  !Related to soil erodibility
  real(r8)          :: dust_emis_fact = 0._r8            ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

!===============================================================================
contains
!===============================================================================

  subroutine oslo_aero_dust_readnl(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file
    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(dust_emis_fact, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: dust_emis_fact")
    call mpi_bcast(soil_erod_file, len(soil_erod_file), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: soil_erod_file")

    call shr_dust_emis_readnl(mpicom, 'drv_flds_in')

    if ((soil_erod_file /= 'none') .and. (.not.is_zender_soil_erod_from_atm())) then
       call endrun(subname//': should not specify soil_erod_file if Zender soil erosion is not in CAM')
    end if

    if (masterproc) then
       if (is_dust_emis_zender()) then
          write(iulog,*) subname,': Zender_2003 dust emission method is being used.'
       end if
       if (is_zender_soil_erod_from_atm()) then
          write(iulog,*) subname,': Zender soil erod file is handled in atm'
          write(iulog,*) subname,': soil_erod_file = ',trim(soil_erod_file)
          write(iulog,*) subname,': dust_emis_fact = ',dust_emis_fact
       end if
    end if

  end subroutine oslo_aero_dust_readnl

  !===============================================================================
  subroutine oslo_aero_dust_init()

    ! local variables
    integer :: imode

    if (is_zender_soil_erod_from_atm()) then
       call soil_erod_init( dust_emis_fact, soil_erod_file )
    end if

    ! Set module variables
    tracerMap(1) = l_dst_a2
    tracerMap(2) = l_dst_a3

    dust_names(:)="      "
    do imode=1,numberOfDustModes
       dust_names(imode) = cnst_name(tracerMap(imode))
    end do

  end subroutine oslo_aero_dust_init

  !===============================================================================
  subroutine oslo_aero_dust_emis(lchnk, ncol, dstflx, cflx)

    !-----------------------------------------------------------------------
    ! Purpose: Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model and
    ! the soil erodibility factor is applied here if the zender scheme is used
    !-----------------------------------------------------------------------

    ! Arguments:
    integer  , intent(in)    :: lchnk
    integer  , intent(in)    :: ncol
    real(r8) , intent(in)    :: dstflx(pcols,4)
    real(r8) , intent(inout) :: cflx(pcols,pcnst) ! Surface fluxes

    ! Local variables
    integer  :: icol,imode
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: totalEmissionFlux(pcols)

    ! Note that following CESM use of "dust_emis_fact", the emissions are
    ! scaled by the INVERSE of the factor!!
    ! There is another random scale factor of 1.15 there. Adapting the exact
    ! same formulation as MAM now and tune later
    ! As of NE-380: Oslo dust emissions are 2/3 of CAM emissions
    ! gives better AOD close to dust sources

    if (is_zender_soil_erod_from_atm()) then
       ! Filter away unreasonable values for soil erodibility
       ! (using low values e.g. gives emissions in greenland..)
       where(soil_erodibility(:,lchnk) < 0.1_r8)
          soil_erod_tmp(:)=0.0_r8
       elsewhere
          soil_erod_tmp(:)=soil_erodibility(:,lchnk)
       end where

       totalEmissionFlux(:) = 0.0_r8
       do icol=1,ncol
          totalEmissionFlux(icol) = totalEmissionFlux(icol) + sum(dstflx(icol,:))
       end do

       do imode = 1,numberOfDustModes
          cflx(:ncol, tracerMap(imode)) = -1.0_r8*emis_fraction_in_mode(imode) &
               *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8
       end do

    else ! Leung emissions

       totalEmissionFlux(:) = 0.0_r8
       do icol=1,ncol
          totalEmissionFlux(icol) = totalEmissionFlux(icol) + sum(dstflx(icol,:))
       end do

       do imode = 1,numberOfDustModes
          cflx(:ncol, tracerMap(imode)) = -1.0_r8*emis_fraction_in_mode(imode) &
               *totalEmissionFlux(:ncol)/(dust_emis_fact)
       end do
    end if

  end subroutine oslo_aero_dust_emis

end module oslo_aero_dust
