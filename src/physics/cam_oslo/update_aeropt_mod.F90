module update_aeropt_mod

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use ppgrid            , only : pcols, pver
  use commondefinitions , only : nmodes, nbmodes
  use opttab            , only : cate, cat, fac, faq, fbc
  use aeropt_mod        , only : bep1, beg2to3, beb4, beg5to10
  use aeropt_mod        , only : bex440, bax440, bex500, bax500, bax550
  use aeropt_mod        , only : bex670, bax670, bex870, bax870
  use aeropt_mod        , only : bex550lt1, bex550gt1, backscx550

  implicit none

  type, public :: extinction_coeffs_type
    ! Modal total and absorption extiction coefficients (for AeroCom)
    ! for 440nm, 500nm, 550nm, 670nm and 870nm, and for d<1um (lt1) and d>1um (gt1).
    ! March 2009: + backscatter coefficient, backsc550 (km-1 sr-1).

     real(r8), allocatable :: bext440(:,:,:)
     real(r8), allocatable :: babs440(:,:,:)
     real(r8), allocatable :: bext500(:,:,:)
     real(r8), allocatable :: babs500(:,:,:)
     real(r8), allocatable :: bext550(:,:,:)
     real(r8), allocatable :: babs550(:,:,:)
     real(r8), allocatable :: bext670(:,:,:)
     real(r8), allocatable :: babs670(:,:,:)
     real(r8), allocatable :: bext870(:,:,:)
     real(r8), allocatable :: babs870(:,:,:)
     real(r8), allocatable :: bebg440(:,:,:)
     real(r8), allocatable :: bebg500(:,:,:)
     real(r8), allocatable :: bebg550(:,:,:)
     real(r8), allocatable :: babg550(:,:,:)
     real(r8), allocatable :: bebg670(:,:,:)
     real(r8), allocatable :: bebg870(:,:,:)
     real(r8), allocatable :: bebc440(:,:,:)
     real(r8), allocatable :: bebc500(:,:,:)
     real(r8), allocatable :: bebc550(:,:,:)
     real(r8), allocatable :: babc550(:,:,:)
     real(r8), allocatable :: bebc670(:,:,:)
     real(r8), allocatable :: bebc870(:,:,:)
     real(r8), allocatable :: beoc440(:,:,:)
     real(r8), allocatable :: beoc500(:,:,:)
     real(r8), allocatable :: beoc550(:,:,:)
     real(r8), allocatable :: baoc550(:,:,:)
     real(r8), allocatable :: beoc670(:,:,:)
     real(r8), allocatable :: beoc870(:,:,:)
     real(r8), allocatable :: besu440(:,:,:)
     real(r8), allocatable :: besu500(:,:,:)
     real(r8), allocatable :: besu550(:,:,:)
     real(r8), allocatable :: basu550(:,:,:)
     real(r8), allocatable :: besu670(:,:,:)
     real(r8), allocatable :: besu870(:,:,:)
     real(r8), allocatable :: bebg550lt1(:,:,:)
     real(r8), allocatable :: bebg550gt1(:,:,:)
     real(r8), allocatable :: bebc550lt1(:,:,:)
     real(r8), allocatable :: bebc550gt1(:,:,:)
     real(r8), allocatable :: beoc550lt1(:,:,:)
     real(r8), allocatable :: beoc550gt1(:,:,:)
     real(r8), allocatable :: besu550lt1(:,:,:)
     real(r8), allocatable :: besu550gt1(:,:,:)
     real(r8), allocatable :: backsc550(:,:,:)

   contains

     procedure :: allocate_coeffs
     procedure :: zero_coeffs
     procedure :: update_coeffs

  end type extinction_coeffs_type

  type(extinction_coeffs_type), public :: extinction_coeffs
  type(extinction_coeffs_type), public :: extinction_coeffsn

  public :: intaeropt0
  public :: intaeropt1
  public :: intaeropt2to3
  public :: intaeropt4
  public :: intaeropt5to10

! ==========================================================
contains
! ==========================================================

  subroutine allocate_coeffs(this)

    class(extinction_coeffs_type) :: this

    allocate(this_coeffs%bext440(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babs440(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bext500(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babs500(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bext550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babs550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bext670(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babs670(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bext870(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babs870(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg440(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg500(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babg550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg670(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg870(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc440(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc500(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%babc550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc670(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc870(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc440(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc500(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%baoc550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc670(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc870(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu440(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu500(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%basu550(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu670(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu870(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg550lt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebg550gt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc550lt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%bebc550gt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc550lt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%beoc550gt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu550lt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%besu550gt1(pcols,pver,0:nbmodes))
    allocate(this_coeffs%backsc550(pcols,pver,0:nbmodes))

  end subroutine allocate_coeffs

  ! ==========================================================
  subroutine zero_coeffs(this, kcomp, ncol)

    class(extinction_coeffs_type) :: this
    integer                      , intent(in)    :: kcomp
    integer                      , intent(in)    :: ncol

    integer :: k
    integer :: icol

    ! initialize all output fields to zero
    do k=1,pver
       do icol=1,ncol
          this%bext440(icol,k,kcomp) = 0.0_r8
          this%babs440(icol,k,kcomp) = 0.0_r8
          this%bext500(icol,k,kcomp) = 0.0_r8
          this%babs500(icol,k,kcomp) = 0.0_r8
          this%bext550(icol,k,kcomp) = 0.0_r8
          this%babs550(icol,k,kcomp) = 0.0_r8
          this%bext670(icol,k,kcomp) = 0.0_r8
          this%babs670(icol,k,kcomp) = 0.0_r8
          this%bext870(icol,k,kcomp) = 0.0_r8
          this%babs870(icol,k,kcomp) = 0.0_r8
          this%bebg440(icol,k,kcomp) = 0.0_r8
          this%bebg500(icol,k,kcomp) = 0.0_r8
          this%bebg550(icol,k,kcomp) = 0.0_r8
          this%babg550(icol,k,kcomp) = 0.0_r8
          this%bebg670(icol,k,kcomp) = 0.0_r8
          this%bebg870(icol,k,kcomp) = 0.0_r8
          this%bebc440(icol,k,kcomp) = 0.0_r8
          this%bebc500(icol,k,kcomp) = 0.0_r8
          this%bebc550(icol,k,kcomp) = 0.0_r8
          this%babc550(icol,k,kcomp) = 0.0_r8
          this%bebc670(icol,k,kcomp) = 0.0_r8
          this%bebc870(icol,k,kcomp) = 0.0_r8
          this%beoc440(icol,k,kcomp) = 0.0_r8
          this%beoc500(icol,k,kcomp) = 0.0_r8
          this%beoc550(icol,k,kcomp) = 0.0_r8
          this%baoc550(icol,k,kcomp) = 0.0_r8
          this%beoc670(icol,k,kcomp) = 0.0_r8
          this%beoc870(icol,k,kcomp) = 0.0_r8
          this%besu440(icol,k,kcomp) = 0.0_r8
          this%besu500(icol,k,kcomp) = 0.0_r8
          this%besu550(icol,k,kcomp) = 0.0_r8
          this%basu550(icol,k,kcomp) = 0.0_r8
          this%besu670(icol,k,kcomp) = 0.0_r8
          this%besu870(icol,k,kcomp) = 0.0_r8
          this%bebg550lt1(icol,k,kcomp) = 0.0_r8
          this%bebg550gt1(icol,k,kcomp) = 0.0_r8
          this%bebc550lt1(icol,k,kcomp) = 0.0_r8
          this%bebc550gt1(icol,k,kcomp) = 0.0_r8
          this%beoc550lt1(icol,k,kcomp) = 0.0_r8
          this%beoc550gt1(icol,k,kcomp) = 0.0_r8
          this%besu550lt1(icol,k,kcomp) = 0.0_r8
          this%besu550gt1(icol,k,kcomp) = 0.0_r8
          this%backsc550(icol,k,kcomp) = 0.0_r8
       end do
    end do

  end subroutine zero_coeffs

  ! ==========================================================
  subroutine update_coeffs(this, icol, k, kcomp)

    class(extinction_coeffs_type) :: this
    integer  , intent(in) :: icol
    integer  , intent(in) :: k
    integer  , intent(in) :: kcomp
    real(r8) , intent(in) :: opt(:)

    this%bext440(icol,k,kcomp)    = opt(1)
    this%bext500(icol,k,kcomp)    = opt(2)
    this%bext670(icol,k,kcomp)    = opt(3)
    this%bext870(icol,k,kcomp)    = opt(4)
    this%bebg440(icol,k,kcomp)    = opt(5)
    this%bebg500(icol,k,kcomp)    = opt(6)
    this%bebg670(icol,k,kcomp)    = opt(7)
    this%bebg870(icol,k,kcomp)    = opt(8)
    this%bebc440(icol,k,kcomp)    = opt(9)
    this%bebc500(icol,k,kcomp)    = opt(10)
    this%bebc670(icol,k,kcomp)    = opt(11)
    this%bebc870(icol,k,kcomp)    = opt(12)
    this%beoc440(icol,k,kcomp)    = opt(13)
    this%beoc500(icol,k,kcomp)    = opt(14)
    this%beoc670(icol,k,kcomp)    = opt(15)
    this%beoc870(icol,k,kcomp)    = opt(16)
    this%besu440(icol,k,kcomp)    = opt(17)
    this%besu500(icol,k,kcomp)    = opt(18)
    this%besu670(icol,k,kcomp)    = opt(19)
    this%besu870(icol,k,kcomp)    = opt(20)
    this%babs440(icol,k,kcomp)    = opt(21)
    this%babs500(icol,k,kcomp)    = opt(22)
    this%babs550(icol,k,kcomp)    = opt(23)
    this%babs670(icol,k,kcomp)    = opt(24)
    this%babs870(icol,k,kcomp)    = opt(25)
    this%bebg550lt1(icol,k,kcomp) = opt(26)
    this%bebg550gt1(icol,k,kcomp) = opt(27)
    this%bebc550lt1(icol,k,kcomp) = opt(28)
    this%bebc550gt1(icol,k,kcomp) = opt(29)
    this%beoc550lt1(icol,k,kcomp) = opt(30)
    this%beoc550gt1(icol,k,kcomp) = opt(31)
    this%besu550lt1(icol,k,kcomp) = opt(32)
    this%besu550gt1(icol,k,kcomp) = opt(33)
    this%backsc550(icol,k,kcomp)  = opt(34)
    this%babg550(icol,k,kcomp)    = opt(35)
    this%babc550(icol,k,kcomp)    = opt(36)
    this%baoc550(icol,k,kcomp)    = opt(37)
    this%basu550(icol,k,kcomp)    = opt(38)
    this%bebg550(icol,k,kcomp)    = opt(26)+opt(27)
    this%bebc550(icol,k,kcomp)    = opt(28)+opt(29)
    this%beoc550(icol,k,kcomp)    = opt(30)+opt(31)
    this%besu550(icol,k,kcomp)    = opt(32)+opt(33)
    this%bext550(icol,k,kcomp)    = bebg550(icol,k,kcomp)+bebc550(icol,k,kcomp) &
                                   +beoc550(icol,k,kcomp)+besu550(icol,k,kcomp)

  end subroutine update_coeffs

  ! ==========================================================
  subroutine intaeropt0 (lchnk, ncol, Nnatk, extinction_coeffs)

    ! Arguments
    integer                , intent(in)    :: lchnk                     ! chunk identifier
    integer                , intent(in)    :: ncol                      ! number of atmospheric columns
    real(r8)               , intent(in)    :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    type(extinction_coeffs), intent(inout) :: extinction_coeffs

    ! Local variables
    integer i, iv, ierr, k, kcomp, icol

    kcomp=0
    extinction_coeffs%zero_coeffs(kcomp, ncol)

    ! BC(ax) mode: update below to non-xero values
    do k = 1,pver
       do icol = 1,ncol
          if(Nnatk(icol,k,kcomp).gt.0) then
             bext440(icol,k,kcomp)=bex440
             babs440(icol,k,kcomp)=bax440
             bext500(icol,k,kcomp)=bex500
             babs500(icol,k,kcomp)=bax500
             bext550(icol,k,kcomp)=bex550lt1+bex550gt1
             babs550(icol,k,kcomp)=bax550
             bext670(icol,k,kcomp)=bex670
             babs670(icol,k,kcomp)=bax670
             bext870(icol,k,kcomp)=bex870
             babs870(icol,k,kcomp)=bax870
             bebg440(icol,k,kcomp)=bex440
             bebg500(icol,k,kcomp)=bex500
             bebg550(icol,k,kcomp)=bex550lt1+bex550gt1
             babg550(icol,k,kcomp)=bax550
             bebg670(icol,k,kcomp)=bex670
             bebg870(icol,k,kcomp)=bex870
             bebg550lt1(icol,k,kcomp)=bex550lt1
             bebg550gt1(icol,k,kcomp)=bex550gt1
             backsc550(icol,k,kcomp)=backscx550
          endif
       end do ! icol
    end do ! k

  end subroutine intaeropt0

  ! ==========================================================
  subroutine intaeropt1 (lchnk, ncol, xrh, irh1, mplus10, &
                         Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1, &
                         extinction_coeffs)

    ! arguments
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    integer  , intent(in) :: mplus10                    ! mode number (0) or number + 10 (1)
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xfombg(pcols,pver)         ! SOA/(SOA+H2SO4) for the background mode
    integer  , intent(in) :: ifombg1(pcols,pver)
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)
    type(extinction_coeffs) , intent(inout) :: extinction_coeffs

    ! local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, ierr, irelh, ifombg, ictot, ifac, kcomp, k, icol, kc10
    ! Temporary storage of often used array elements
    integer  :: t_irh1, t_irh2, t_ifo1, t_ifo2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xrh, t_rh1, t_rh2, t_fombg1, t_fombg2, t_xfombg
    real(r8) :: t_xct, t_cat1, t_cat2
    real(r8) :: d2mx(4), dxm1(4), invd(4)
    real(r8) :: opt4d(2,2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: opt1, opt2, opt(38)

    parameter :: (e=2.718281828_r8, eps=1.0e-60_r8)

    ! SO4/SOA(Ait) mode:
    kcomp = 1
    extinction_coeffs%zero_coeffs(kcomp, ncol)

    if(mplus10 == 0) then
       kc10 = kcomp
    else
       write(*,*) "mplus10=1 is no loger an option for kcomp=1."
       stop
    endif

     do k=1,pver 
        do icol=1,ncol

           if(Nnatk(icol,k,kc10).gt.0) then

              ! Collect all the vector elements into temporary storage
              ! to avoid cache conflicts and excessive cross-referencing

              t_irh1 = irh1(icol,k)
              t_irh2 = t_irh1+1
              t_ifo1 = ifombg1(icol,k)
              t_ifo2 = t_ifo1+1
              t_ict1 = ict1(icol,k,kcomp)
              t_ict2 = t_ict1+1
              t_ifc1 = ifac1(icol,k,kcomp)
              t_ifc2 = t_ifc1+1

              t_rh1  = rh(t_irh1)
              t_rh2  = rh(t_irh2)
              t_fombg1 = fombg(t_ifo1)
              t_fombg2 = fombg(t_ifo2)
              t_cat1 = cate(kcomp,t_ict1)
              t_cat2 = cate(kcomp,t_ict2)
              t_fac1 = fac(t_ifc1)
              t_fac2 = fac(t_ifc2)

              t_xrh  = xrh(icol,k)
              t_xct  = xct(icol,k,kc10)
              t_xfac = xfac(icol,k,kcomp)
              t_xfombg = xfombg(icol,k)

              ! partial lengths along each dimension (1-4) for interpolation 
              d2mx(1) = (t_rh2-t_xrh)
              dxm1(1) = (t_xrh-t_rh1)
              invd(1) = 1.0_r8/(t_rh2-t_rh1)
              d2mx(2) = (t_fombg2-t_xfombg)
              dxm1(2) = (t_xfombg-t_fombg1)
              invd(2) = 1.0_r8/(t_fombg2-t_fombg1)
              d2mx(3) = (t_cat2-t_xct)
              dxm1(3) = (t_xct-t_cat1)
              invd(3) = 1.0_r8/(t_cat2-t_cat1)
              d2mx(4) = (t_fac2-t_xfac)
              dxm1(4) = (t_xfac-t_fac1)
              invd(4) = 1.0_r8/(t_fac2-t_fac1)

              do iv=1,38  ! variable number
                 ! end points as basis for multidimentional linear interpolation  
                 opt4d(1,1,1,1) = bep1(iv,t_irh1,t_ifo1,t_ict1,t_ifc1)
                 opt4d(1,1,1,2) = bep1(iv,t_irh1,t_ifo1,t_ict1,t_ifc2)
                 opt4d(1,1,2,1) = bep1(iv,t_irh1,t_ifo1,t_ict2,t_ifc1)
                 opt4d(1,1,2,2) = bep1(iv,t_irh1,t_ifo1,t_ict2,t_ifc2)
                 opt4d(1,2,1,1) = bep1(iv,t_irh1,t_ifo2,t_ict1,t_ifc1)
                 opt4d(1,2,1,2) = bep1(iv,t_irh1,t_ifo2,t_ict1,t_ifc2)
                 opt4d(1,2,2,1) = bep1(iv,t_irh1,t_ifo2,t_ict2,t_ifc1)
                 opt4d(1,2,2,2) = bep1(iv,t_irh1,t_ifo2,t_ict2,t_ifc2)
                 opt4d(2,1,1,1) = bep1(iv,t_irh2,t_ifo1,t_ict1,t_ifc1)
                 opt4d(2,1,1,2) = bep1(iv,t_irh2,t_ifo1,t_ict1,t_ifc2)
                 opt4d(2,1,2,1) = bep1(iv,t_irh2,t_ifo1,t_ict2,t_ifc1)
                 opt4d(2,1,2,2) = bep1(iv,t_irh2,t_ifo1,t_ict2,t_ifc2)
                 opt4d(2,2,1,1) = bep1(iv,t_irh2,t_ifo2,t_ict1,t_ifc1)
                 opt4d(2,2,1,2) = bep1(iv,t_irh2,t_ifo2,t_ict1,t_ifc2)
                 opt4d(2,2,2,1) = bep1(iv,t_irh2,t_ifo2,t_ict2,t_ifc1)
                 opt4d(2,2,2,2) = bep1(iv,t_irh2,t_ifo2,t_ict2,t_ifc2)

                 ! interpolation in the fac, cat and fombg dimensions
                 call lininterpol4dim (d2mx, dxm1, invd, opt4d, opt1, opt2)

                 ! finally, interpolation in the rh dimension 
                 opt(iv)=((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) / (t_rh2-t_rh1)    
              end do ! iv=1,38 

              ! determin extinction coefficient
              extinction_coeffs%update_coeffs(icol, k, kcomp, opt)

           end if
        end do ! end of icol loop
     end do  ! end of k loop

  end subroutine intaeropt1

  ! ==========================================================
  subroutine intaeropt2to3 (lchnk, ncol, xrh, irh1, mplus10, &
       Nnatk, xct, ict1, xfac, ifac1, extinction_coeffs)

    !   Extended by Alf Kirkevaag to include SOA in September 2015

    ! Arguments
    integer  , intent(in) :: lchnk                       ! chunk identifier
    integer  , intent(in) :: ncol                        ! number of atmospheric columns
    integer  , intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)        
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)        

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, kcomp, k, icol, kc10
    ! Temporary storage of often used array elements
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) :: t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
    real(r8) :: d2mx(3), dxm1(3), invd(3)
    real(r8) :: opt3d(2,2,2)
    real(r8) :: opt1, opt2, opt(38)

    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! SO4(Ait), BC(Ait) and OC(Ait) modes:

    do kcomp=2,3        
       extinction_coeffs%zero_coeffs(kcomp, ncol)
    end do

    kcomp = 2 ! kcomp=3 is no longer used        
    do k=1,pver 
       do icol=1,ncol

          if(Nnatk(icol,k,kc10).gt.0) then

             !      Collect all the vector elements into temporary storage
             !      to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)

             !     partial lengths along each dimension (1-4) for interpolation 
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_cat2-t_xct)
             dxm1(2) = (t_xct-t_cat1)
             invd(2) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(3) = (t_fac2-t_xfac)
             dxm1(3) = (t_xfac-t_fac1)
             invd(3) = 1.0_r8/(t_fac2-t_fac1)

             do iv=1,38  ! variable number

                !  end points as basis for multidimentional linear interpolation  
                opt3d(1,1,1)=bep2to3(iv,t_irh1,t_ict1,t_ifc1,kcomp)
                opt3d(1,1,2)=bep2to3(iv,t_irh1,t_ict1,t_ifc2,kcomp)
                opt3d(1,2,1)=bep2to3(iv,t_irh1,t_ict2,t_ifc1,kcomp)
                opt3d(1,2,2)=bep2to3(iv,t_irh1,t_ict2,t_ifc2,kcomp)
                opt3d(2,1,1)=bep2to3(iv,t_irh2,t_ict1,t_ifc1,kcomp)
                opt3d(2,1,2)=bep2to3(iv,t_irh2,t_ict1,t_ifc2,kcomp)
                opt3d(2,2,1)=bep2to3(iv,t_irh2,t_ict2,t_ifc1,kcomp)
                opt3d(2,2,2)=bep2to3(iv,t_irh2,t_ict2,t_ifc2,kcomp)

                !     interpolation in the (fac and) cat dimension
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, opt1, opt2)

                !     finally, interpolation in the rh dimension
                opt(iv)=((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) /(t_rh2-t_rh1)    

             end do ! iv=1,38 

             ! determine extinction coefficient
             extinction_coeffs%update_coeffs(icol, k, kcomp, opt)

          end if ! Nnatk > 0
       end do ! icol
    end do ! k

  end subroutine intaeropt2to3

  ! ==========================================================
  subroutine intaeropt4 (lchnk, ncol, xrh, irh1, mplus10, Nnatk,   &
       xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, xfaq, ifaq1, &
       extinction_coeffs)

    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    integer  , intent(in) :: mplus10                    ! mode number (0) or number + 10 (1)
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
    real(r8) , intent(in) :: xfbcbg(pcols,pver)
    integer  , intent(in) :: ifbcbg1(pcols,pver)
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)        
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)        
    real(r8) , intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer  , intent(in) :: ifaq1(pcols,pver,nbmodes)
    type(extinction_coeffs), intent(inout) :: extinction_coeffs

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, kcomp, k, icol, kc10
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2,  t_ifa1, t_ifa2
    real(r8) :: t_fbcbg1, t_fbcbg2  
    integer  :: t_ifb1, t_ifb2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xrh, t_xct, t_rh1, t_rh2
    real(r8) :: t_cat1, t_cat2
    real(r8) :: t_xfbcbg
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: opt1, opt2, opt(38)
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! BC&OC(Ait) mode: 
    kcomp = 4
    extinction_coeffs%zero_coeffs(kcomp, ncol)

    if(mplus10==0) then
       kc10=kcomp
    else
       kc10=kcomp+10
    endif

    do k=1,pver 
       do icol=1,ncol
          if(Nnatk(icol,k,kc10).gt.0) then
             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ifb1 = ifbcbg1(icol,k)
             t_ifb2 = t_ifb1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifa1 = ifaq1(icol,k,kcomp)
             t_ifa2 = t_ifa1+1

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_fbcbg1 = fbcbg(t_ifb1)
             t_fbcbg2 = fbcbg(t_ifb2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_faq1 = faq(t_ifa1)
             t_faq2 = faq(t_ifa2)

             t_xrh  = xrh(icol,k)
             t_xfbcbg = xfbcbg(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)
             t_xfaq = xfaq(icol,k,kcomp)

             ! partial lengths along each dimension (1-5) for interpolation 
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_fbcbg2-t_xfbcbg)
             dxm1(2) = (t_xfbcbg-t_fbcbg1)
             invd(2) = 1.0_r8/(t_fbcbg2-t_fbcbg1)
             d2mx(3) = (t_cat2-t_xct)
             dxm1(3) = (t_xct-t_cat1)
             invd(3) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(4) = (t_fac2-t_xfac)
             dxm1(4) = (t_xfac-t_fac1)
             invd(4) = 1.0_r8/(t_fac2-t_fac1)
             d2mx(5) = (t_faq2-t_xfaq)
             dxm1(5) = (t_xfaq-t_faq1)
             invd(5) = 1.0_r8/(t_faq2-t_faq1)


             do iv=1,38  ! variable number

                opt5d(1,1,1,1,1)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,1,1,1,2)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,1,1,2,1)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,1,1,2,2)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,1,2,1,1)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,1,2,1,2)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,1,2,2,1)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,1,2,2,2)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(1,2,1,1,1)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,2,1,1,2)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,2,1,2,1)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,2,1,2,2)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,2,2,1,1)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,2,2,1,2)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,2,2,2,1)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,2,2,2,2)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,1,1,1,1)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,1,1,1,2)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,1,1,2,1)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,1,1,2,2)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,1,2,1,1)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,1,2,1,2)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,1,2,2,1)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,1,2,2,2)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,2,1,1,1)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,2,1,1,2)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,2,1,2,1)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,2,1,2,2)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,2,2,1,1)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,2,2,1,2)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,2,2,2,1)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,2,2,2,2)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac, cat and fbcbg dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, opt1, opt2)

                ! finally, interpolation in the rh dimension 
                opt(iv) = ((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) /(t_rh2-t_rh1)    

             end do ! iv=1,38 

             ! determine extinction coefficient
             extinction_coeffs%update_coeffs(icol, k, kcomp, opt)

          end if ! Nnatk > 0
       end do ! icol
    end do ! k

  end subroutine intaeropt4

  ! ==========================================================
  subroutine intaeropt5to10 (lchnk, ncol, xrh, irh1, Nnatk,    &
       xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
       extinction_coeffs)

    ! Arguments
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)        
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! modal (OC+BC)/(SO4+BC+OC)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfbc(pcols,pver,nbmodes)   ! modal BC/(OC+BC)
    integer  , intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer  , intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, kcomp, k, icol
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2
    integer  :: t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: t_fbc1, t_fbc2, t_xfbc
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xrh, t_xct, t_rh1, t_rh2
    real(r8) :: t_cat1, t_cat2
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: opt1, opt2, opt(38)
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.):

    do kcomp=5,10        
       ! zero extinction coefficients for this kcomp
       extinction_coeffs%zero_coeffs(kcomp, ncol)

       do k=1,pver 
          do icol=1,ncol
             if(Nnatk(icol,k,kcomp).gt.0) then
                ! Collect all the vector elements into temporary storage
                ! to avoid cache conflicts and excessive cross-referencing

                t_irh1 = irh1(icol,k)
                t_irh2 = t_irh1+1
                t_ict1 = ict1(icol,k,kcomp)
                t_ict2 = t_ict1+1
                t_ifc1 = ifac1(icol,k,kcomp)
                t_ifc2 = t_ifc1+1

                t_ifb1 = ifbc1(icol,k,kcomp)
                t_ifb2 = t_ifb1+1
                t_ifa1 = ifaq1(icol,k,kcomp)
                t_ifa2 = t_ifa1+1

                t_rh1  = rh(t_irh1)
                t_rh2  = rh(t_irh2)
                t_cat1 = cat(kcomp,t_ict1)
                t_cat2 = cat(kcomp,t_ict2)
                t_fac1 = fac(t_ifc1)
                t_fac2 = fac(t_ifc2)
                t_fbc1 = fbc(t_ifb1)
                t_fbc2 = fbc(t_ifb2)
                t_faq1 = faq(t_ifa1)
                t_faq2 = faq(t_ifa2)

                t_xrh  = xrh(icol,k)
                t_xct  = xct(icol,k,kcomp)
                t_xfac = xfac(icol,k,kcomp)
                t_xfbc = xfbc(icol,k,kcomp)
                t_xfaq = xfaq(icol,k,kcomp)
                
                ! partial lengths along each dimension (1-5) for interpolation 
                d2mx(1) = (t_rh2-t_xrh)
                dxm1(1) = (t_xrh-t_rh1)
                invd(1) = 1.0_r8/(t_rh2-t_rh1)
                d2mx(2) = (t_cat2-t_xct)
                dxm1(2) = (t_xct-t_cat1)
                invd(2) = 1.0_r8/(t_cat2-t_cat1)
                d2mx(3) = (t_fac2-t_xfac)
                dxm1(3) = (t_xfac-t_fac1)
                invd(3) = 1.0_r8/(t_fac2-t_fac1)
                d2mx(4) = (t_fbc2-t_xfbc)
                dxm1(4) = (t_xfbc-t_fbc1)
                invd(4) = 1.0_r8/(t_fbc2-t_fbc1)
                d2mx(5) = (t_faq2-t_xfaq)
                dxm1(5) = (t_xfaq-t_faq1)
                invd(5) = 1.0_r8/(t_faq2-t_faq1)


                do iv=1,38  ! variable number
                   opt5d(1,1,1,1,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, opt1, opt2)

                   ! finally, interpolation in the rh dimension 
                   opt(iv) = ((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) /(t_rh2-t_rh1)    

                end do ! iv=1,38 

                ! determine extinction coefficient
                extinction_coeffs%update_coeffs(icol, k, kcomp, opt)

             end if ! Nnatk > 0
          end do ! icol
       end do ! k
    end do ! kcomp

  end subroutine intaeropt5to10

end module update_aeropt_mod

