module dust_model

  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl

  implicit none
  private

  integer, parameter :: numberOfDustModes = 2  !define in aerosoldef?
  character(len=6), public, dimension(10)        :: dust_names

  integer  :: tracerMap(numberOfDustModes) = (/-99, -99/) !index of dust tracers in the modes

  real(r8), parameter        :: emis_fraction_in_mode(numberOfDustModes) = (/0.13_r8, 0.87_r8 /)
  integer, parameter, public :: dust_nbin = numberOfDustModes

  !Related to soil erodibility
  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical, parameter, public :: dust_active = .TRUE.

  ! public routines
  public dust_readnl
  public dust_init
  public dust_emis

!===============================================================================
contains
!===============================================================================

  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: masterproc
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file
    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif
  end subroutine dust_readnl

  !===============================================================================
  subroutine dust_init()

    use soil_erod_mod, only: soil_erod_init
    use aerosoldef,    only: l_dst_a2, l_dst_a3
    use constituents,  only: cnst_name

    integer :: i

    call soil_erod_init( dust_emis_fact, soil_erod_file )

    ! Set module variables
    tracerMap(1) = l_dst_a2
    tracerMap(2) = l_dst_a3

    dust_names(:)="      "
    do i=1,numberOfDustModes
       dust_names(i) = cnst_name(tracerMap(i))
    end do

  end subroutine dust_init

  !===============================================================================
  subroutine dust_emis(state, cam_in)

    !----------------------------------------------------------------------- 
    ! Purpose: Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model and
    ! the soil erodibility factor is applied here.
    !-----------------------------------------------------------------------

    use ppgrid,        only: pcols
    use physics_types, only: physics_state
    use camsrfexch,    only: cam_in_t
    use soil_erod_mod, only: soil_erod_fact, soil_erodibility 

    ! Arguments:
    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t), target, intent(inout) :: cam_in  ! import state

    ! Local variables
    integer           :: lchnk
    integer           :: ncol
    integer           :: i,n
    real(r8)          :: soil_erod_tmp(pcols)
    real(r8)          :: totalEmissionFlux(pcols)
    real(r8), pointer :: cflx(:,:)

    lchnk = state%lchnk
    ncol = state%ncol

    !Filter away unreasonable values for soil erodibility
    !(using low values e.g. gives emissions in greenland..)
    where(soil_erodibility(:,lchnk) .lt. 0.1_r8)
       soil_erod_tmp(:)=0.0_r8
    elsewhere
       soil_erod_tmp(:)=soil_erodibility(:,lchnk)
    end where

    totalEmissionFlux(:)=0.0_r8
    do i=1,ncol
       totalEmissionFlux(i) = totalEmissionFlux(i) + sum(cam_in%dstflx(i,:))
    end do

    !Note that following CESM use of "dust_emis_fact", the emissions are 
    !scaled by the INVERSE of the factor!!
    !There is another random scale factor of 1.15 there. Adapting the exact
    !same formulation as MAM now and tune later
    !As of NE-380: Oslo dust emissions are 2/3 of CAM emissions

    cflx => cam_in%cflx
    do n=1, numberOfDustModes
       cflx(:ncol, tracerMap(n)) = -1.0_r8*emis_fraction_in_mode(n) &
            *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8  ! gives better AOD close to dust sources
    end do

  end subroutine dust_emis

end module dust_model
