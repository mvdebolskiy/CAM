
subroutine opticsAtConstRh (lchnk, ncol, pint, rhoda, Nnatk, xrh, irh1, irf, &
     xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,           &
     xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,       &
     xfombg, ifombg1, vnbc, vaitbc, v_soana,            &
     extinction_coeffs, extinction_coeffsn)


  !     Extra AeroCom diagnostics requiring table look-ups with constant/fixed RH,
  !     i.e. for RH = (/"00","40","55","65","75","85" /) (see opttab.F90)

  use ppgrid
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_history,  only: outfld
  use constituents, only: pcnst
  use opttab
  use const
  use aerosoldef
  use commondefinitions
  use physics_types,   only: physics_state
  use interp_aeropt_mod, only : extinction_coeffs_type
  use interp_aeropt_mod, only : intaeropt0, intaeropt1


  implicit none
  !
  ! Input arguments
  !
  integer,  intent(in) :: lchnk                      ! chunk identifier
  integer,  intent(in) :: ncol                       ! number of atmospheric columns
  real(r8), intent(in) :: pint(pcols,pverp)          ! Model interface pressures (10*Pa)
  real(r8), intent(in) :: rhoda(pcols,pver)          ! Density of dry air (kg/m^3)
  real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
  integer,  intent(in) :: irh1(pcols,pver)
  integer,  intent(in) :: irf
  real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! aerosol mode number concentration  
  real(r8), intent(in) :: vnbc(pcols,pver)
  real(r8), intent(in) :: vaitbc(pcols,pver)
  real(r8), intent(in) :: v_soana(pcols,pver)
  real(r8), intent(in) :: xfombg(pcols,pver)
  integer,  intent(in) :: ifombg1(pcols,pver)
  real(r8), intent(in) :: xfbcbg(pcols,pver)
  integer,  intent(in) :: ifbcbg1(pcols,pver)
  real(r8), intent(in) :: xfbcbgn(pcols,pver)
  integer,  intent(in) :: ifbcbgn1(pcols,pver)
  real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
  integer,  intent(in) :: ict1(pcols,pver,nmodes)        
  real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! facm for use in the interpolations 
  integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
  real(r8), intent(in) :: xfbc(pcols,pver,nbmodes)   ! fbcm for use in the interpolations 
  integer,  intent(in) :: ifbc1(pcols,pver,nbmodes)
  real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! faqm for use in the interpolations 
  integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)
  type(extinction_coeffs_type) , intent(inout) :: extinction_coeffs
  type(extinction_coeffs_type) , intent(inout) :: extinction_coeffsn
  !
  !---------------------------Local variables-----------------------------
  !
  integer  :: i, k, icol, mplus10, irh
  integer  :: iloop
  real(r8) :: deltah
  real(r8) :: dod550rh(pcols), abs550rh(pcols)
  real(r8) :: ec550rh_aer(pcols,pver)
  real(r8) :: abs550rh_aer(pcols,pver)
  real(r8) :: bebglt1t(pcols,pver)
  real(r8) :: bebclt1t(pcols,pver)
  real(r8) :: beoclt1t(pcols,pver)
  real(r8) :: bes4lt1t(pcols,pver)
  real(r8) :: basu550tot(pcols,pver)
  real(r8) :: babc550tot(pcols,pver)
  real(r8) :: baoc550tot(pcols,pver)
  real(r8) :: babc550xt(pcols,pver)
  real(r8) :: baoc550xt(pcols,pver), &
  real(r8) :: ba550x(pcols,pver,nbmp1:nmodes)
  real(r8) :: belt1x(pcols,pver,nbmp1:nmodes)
  ! Additionl AeroCom Phase III output:   
  real(r8) :: ec440rh_aer(pcols,pver)
  real(r8) :: abs440rh_aer(pcols,pver)
  real(r8) :: ec870rh_aer(pcols,pver)
  real(r8) :: abs870rh_aer(pcols,pver)
  real(r8) :: be550lt1_aer(pcols,pver,0:nbmodes)
  real(r8) :: ec550rhlt1_aer(pcols,pver)
  real(r8) :: abs550rh_bc(pcols,pver)
  real(r8) :: abs550rh_oc(pcols,pver)
  real(r8) :: abs550rh_su(pcols,pver)
  real(r8) :: abs550rh_ss(pcols,pver)
  real(r8) :: abs550rh_du(pcols,pver)
  real(r8) :: ec550rhlt1_bc(pcols,pver)
  real(r8) :: ec550rhlt1_oc(pcols,pver)
  real(r8) :: ec550rhlt1_su(pcols,pver)
  real(r8) :: ec550rhlt1_ss(pcols,pver)
  real(r8) :: ec550rhlt1_du(pcols,pver) 
  real(r8) :: bedustlt1(pcols,pver) 
  real(r8) :: bedustgt1(pcols,pver)
  real(r8) :: besslt1(pcols,pver)
  real(r8) :: bessgt1(pcols,pver)
  real(r8) :: bbclt1xt(pcols,pver)
  real(r8) :: boclt1xt(pcols,pver)
  real(r8) :: bocgt1xt(pcols,pver)

  character(len=10) :: modeString
  character(len=20) :: varname
  !--------------------------------------------------

  belt1x(:,:,:) = 0._r8

  do iloop=1,1

     ! BC(ax) mode (hydrophobic, so no rhum needed here):
     call intaeropt0(lchnk, ncol, Nnatk, extinction_coeffs)

     ! SO4(Ait), BC(Ait) and OC(Ait) modes:
     mplus10=0
     call intaeropt1(lchnk, ncol, xrh, irh1, mplus10,  &
          Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1,&
          extinction_coeffs)

     mplus10=0
     call intaeropt2to3(lchnk, ncol, xrh, irh1, mplus10, &
          Nnatk, xct, ict1, xfac, ifac1,                   &
          extinction_coeffs)

     ! BC&OC(Ait) (4), OC&BC(Ait) mode
     mplus10=0
     call intaeropt4(lchnk, ncol, xrh, irh1, mplus10, Nnatk,  &
          xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, xfaq, ifaq1, &
          extinction_coeffs)

     ! SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
     call intaeropt5to10(lchnk, ncol, xrh, irh1, Nnatk,   &
          xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
          extinction_coeffs)

     ! then to the externally mixed SO4(n), BC(n) and OC(n) modes:
     mplus10=1
     call intaeropt2to3(lchnk, ncol, xrh, irh1, mplus10,  &
          Nnatk, xct, ict1, xfac, ifac1,                    &
          extinction_coeffsn)

     ! and finally the BC&OC(n) mode:
     mplus10=1
     call intaeropt4(lchnk, ncol, xrh, irh1, mplus10, Nnatk,    &
          xfbcbgn, ifbcbgn1, xct, ict1, xfac, ifac1, xfaq, ifaq1, &
          extinction_coeffsn)

  end do ! iloop


  ! Initialization
  do k=1,pver  
     do icol=1,ncol
        ec550rh_aer(icol,k)     = 0.0_r8 
        abs550rh_aer(icol,k)    = 0.0_r8 
        ec550rhlt1_aer(icol,k)  = 0.0_r8 
        abs550rh_bc(icol,k)     = 0.0_r8
        abs550rh_oc(icol,k)     = 0.0_r8
        abs550rh_su(icol,k)     = 0.0_r8
        abs550rh_ss(icol,k)     = 0.0_r8
        abs550rh_du(icol,k)     = 0.0_r8
        ec440rh_aer(icol,k)     = 0.0_r8 
        abs440rh_aer(icol,k)    = 0.0_r8 
        ec870rh_aer(icol,k)     = 0.0_r8 
        abs870rh_aer(icol,k)    = 0.0_r8 
        basu550tot(icol,k)      = 0.0_r8 
        babc550tot(icol,k)      = 0.0_r8 
        baoc550tot(icol,k)      = 0.0_r8 
        bebglt1t(icol,k)        = 0.0_r8
        bebclt1t(icol,k)        = 0.0_r8
        beoclt1t(icol,k)        = 0.0_r8
        bes4lt1t(icol,k)        = 0.0_r8
        bedustlt1(icol,k)       = 0.0_r8
        besslt1(icol,k)         = 0.0_r8
     end do
  end do
  do icol=1,ncol
     dod550rh(icol) = 0.0_r8 
     abs550rh(icol) = 0.0_r8 
  end do

  ! Calculation of extinction at given RH and absorption for all r and for r<0.5um
  do k=1,pver
     do icol=1,ncol

        do i=0,10
           ec550rh_aer(icol,k)  = ec550rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffs%bext550(icol,k,i)
           abs550rh_aer(icol,k) = abs550rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%babs550(icol,k,i)
           ec440rh_aer(icol,k)  = ec440rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffs%bext440(icol,k,i)
           abs440rh_aer(icol,k) = abs440rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%babs440(icol,k,i)
           ec870rh_aer(icol,k)  = ec870rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffs%bext870(icol,k,i)
           abs870rh_aer(icol,k) = abs870rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%babs870(icol,k,i)
           basu550tot(icol,k)   = basu550tot(icol,k)   + Nnatk(icol,k,i)*extinction_coeffs%basu550(icol,k,i)
           babc550tot(icol,k)   = babc550tot(icol,k)   + Nnatk(icol,k,i)*extinction_coeffs%babc550(icol,k,i)
           baoc550tot(icol,k)   = baoc550tot(icol,k)   + Nnatk(icol,k,i)*extinction_coeffs%baoc550(icol,k,i)
           bes4lt1t(icol,k)     = bes4lt1t(icol,k)     + Nnatk(icol,k,i)*extinction_coeffs%bes4lt1(icol,k,i)
           bebclt1t(icol,k)     = bebclt1t(icol,k)     + Nnatk(icol,k,i)*extinction_coeffs%bebclt1(icol,k,i)
           beoclt1t(icol,k)     = beoclt1t(icol,k)     + Nnatk(icol,k,i)*extinction_coeffs%beoclt1(icol,k,i)
        enddo
        do i=11,14
           ec550rh_aer(icol,k)  = ec550rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffsn%bext550(icol,k,i-10)
           abs550rh_aer(icol,k) = abs550rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffsn%babs550(icol,k,i-10)
           ec440rh_aer(icol,k)  = ec440rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffsn%bext440(icol,k,i-10)
           abs440rh_aer(icol,k) = abs440rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffsn%babs440(icol,k,i-10)
           ec870rh_aer(icol,k)  = ec870rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffsn%bext870(icol,k,i-10)
           abs870rh_aer(icol,k) = abs870rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffsn%babs870(icol,k,i-10)
           ba550x(icol,k,i)     = extinction_coeffsn%babs550(icol,k,i-10)
           belt1x(icol,k,i)     = bebglt1(icol,k,i-10) !???
        enddo

        do i=6,7
           bedustlt1(icol,k) = bedustlt1(icol,k) + Nnatk(icol,k,i)*bebglt1(icol,k,i)
        enddo
        do i=8,10
           besslt1(icol,k) = besslt1(icol,k) + Nnatk(icol,k,i)*bebglt1(icol,k,i)
        enddo
        ec550rhlt1_du(icol,k) = bedustlt1(icol,k)
        ec550rhlt1_ss(icol,k) = besslt1(icol,k)

        !soa: *(1-v_soan) for the sulfate volume fraction of mode 11
        bbclt1xt(icol,k) = Nnatk(icol,k,12)*belt1x(icol,k,12) &
                         + Nnatk(icol,k,14)*belt1x(icol,k,14)*vnbc(icol,k)
        !soa + v_soan part of mode 11 for the OC volume fraction of that mode
        boclt1xt(icol,k) = Nnatk(icol,k,13)*belt1x(icol,k,13) &
                         + Nnatk(icol,k,14)*belt1x(icol,k,14)*(1.0_r8-vnbc(icol,k)) 

        !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
        ec550rhlt1_su(icol,k) = bes4lt1t(icol,k)                         &  ! condensate
             + Nnatk(icol,k,1)*bebglt1(icol,k,1)*(1.0_r8-v_soana(icol,k))&  ! background, SO4(Ait) mode (1)
             + Nnatk(icol,k,5)*bebglt1(icol,k,5)                            ! background, SO4(Ait75) mode (5)
        ec550rhlt1_bc(icol,k) = bebclt1t(icol,k)+bbclt1xt(icol,k)        &  ! coagulated + n-mode BC (12)
             + Nnatk(icol,k,2)*bebglt1(icol,k,2)                        &  ! background, BC(Ait) mode (2)
             + Nnatk(icol,k,4)*bebglt1(icol,k,4)*vaitbc(icol,k)         &  ! background in OC&BC(Ait) mode (4)
             + Nnatk(icol,k,0)*bebglt1(icol,k,0)                           ! background, BC(ax) mode (0)
        !soa + v_soan part of mode 11 for the OC volume fraction of that mode
        ec550rhlt1_oc(icol,k) = beoclt1t(icol,k)+boclt1xt(icol,k)        &  ! coagulated + n-mode OC (13)
             + Nnatk(icol,k,3)*bebglt1(icol,k,3)                        &  ! background, OC(Ait) mode (3)
             + Nnatk(icol,k,4)*bebglt1(icol,k,4)*(1.0_r8-vaitbc(icol,k))&  ! background in OC&BC(Ait) mode (4)
             + Nnatk(icol,k,1)*bebglt1(icol,k,1)*v_soana(icol,k)

        ec550rhlt1_aer(icol,k) = ec550rhlt1_su(icol,k)+ec550rhlt1_bc(icol,k) &
                               + ec550rhlt1_oc(icol,k) + ec550rhlt1_ss(icol,k)+ec550rhlt1_du(icol,k)
        ec550rhlt1_aer(icol,k) = 1.e-3_r8*ec550rhlt1_aer(icol,k)

        abs550rh_du(icol,k) = Nnatk(icol,k,6)*babg550(icol,k,6) &
                            + Nnatk(icol,k,7)*babg550(icol,k,7)
        abs550rh_ss(icol,k) = Nnatk(icol,k,8)*babg550(icol,k,8) &
                            + Nnatk(icol,k,9)*babg550(icol,k,9) &
                            + Nnatk(icol,k,10)*babg550(icol,k,10)
        !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
        abs550rh_su(icol,k) = basu550tot(icol,k)                   &  ! condensate:w
             + (1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*babg550(icol,k,1) &  ! background, SO4(Ait) mode (1)
             + Nnatk(icol,k,5)*babg550(icol,k,5)    ! background, SO4(Ait75) mode (5)

        !soa: *(1-v_soan) for the sulfate volume fraction
        babc550xt(icol,k) = Nnatk(icol,k,12)*ba550x(icol,k,12)  &
                          + Nnatk(icol,k,14)*ba550x(icol,k,14)*vnbc(icol,k)

        baoc550xt(icol,k) = Nnatk(icol,k,13)*ba550x(icol,k,13) &
                          + Nnatk(icol,k,14)*ba550x(icol,k,14)*(1.0_r8-vnbc(icol,k)) 

        abs550rh_bc(icol,k) = babc550tot(icol,k)+babc550xt(icol,k) &     ! coagulated + n-mode BC (12)
             + Nnatk(icol,k,2)*babg550(icol,k,2) &        ! background, BC(Ait) mode (2)
             + vaitbc(icol,k)*Nnatk(icol,k,4)*babg550(icol,k,4) &        ! background in OC&BC(Ait) mode (4)
             + Nnatk(icol,k,0)*babg550(icol,k,0)          ! background, BC(ax) mode (0)

        abs550rh_oc(icol,k) = baoc550tot(icol,k)+baoc550xt(icol,k) &     ! coagulated + n-mode OC (13)
             + v_soana(icol,k)*Nnatk(icol,k,1)*babg550(icol,k,1) &       ! SOA fraction of mode 1
             + Nnatk(icol,k,3)*babg550(icol,k,3) &       ! background, OC(Ait) mode (3)
             + (1.0_r8-vaitbc(icol,k))*Nnatk(icol,k,4)*babg550(icol,k,4)         ! background in OC&BC(Ait) mode (4)

        deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
        dod550rh(icol) = dod550rh(icol)+ec550rh_aer(icol,k)*deltah
        abs550rh(icol) = abs550rh(icol)+abs550rh_aer(icol,k)*deltah

        ec550rh_aer(icol,k)  = 1.e-3_r8*ec550rh_aer(icol,k)
        abs550rh_aer(icol,k) = 1.e-3_r8*abs550rh_aer(icol,k)
        ec440rh_aer(icol,k)  = 1.e-3_r8*ec440rh_aer(icol,k)
        abs440rh_aer(icol,k) = 1.e-3_r8*abs440rh_aer(icol,k)
        ec870rh_aer(icol,k)  = 1.e-3_r8*ec870rh_aer(icol,k)
        abs870rh_aer(icol,k) = 1.e-3_r8*abs870rh_aer(icol,k)

        abs550rh_bc(icol,k)  = 1.e-3_r8*abs550rh_bc(icol,k)
        abs550rh_oc(icol,k)  = 1.e-3_r8*abs550rh_oc(icol,k)
        abs550rh_su(icol,k)  = 1.e-3_r8*abs550rh_su(icol,k)
        abs550rh_ss(icol,k)  = 1.e-3_r8*abs550rh_ss(icol,k)
        abs550rh_du(icol,k)  = 1.e-3_r8*abs550rh_du(icol,k)

     enddo
  enddo

  if(irf.eq.1) then

     call outfld('ECDRYAER',ec550rh_aer,pcols,lchnk)
     call outfld('ABSDRYAE',abs550rh_aer,pcols,lchnk)
     call outfld('OD550DRY',dod550rh,pcols,lchnk)       ! 2D variable
     call outfld('AB550DRY',abs550rh,pcols,lchnk)       ! 2D variable
     call outfld('ECDRY440',ec440rh_aer,pcols,lchnk)
     call outfld('ABSDR440',abs440rh_aer,pcols,lchnk)
     call outfld('ECDRY870',ec870rh_aer,pcols,lchnk)
     call outfld('ABSDR870',abs870rh_aer,pcols,lchnk)
     call outfld('ECDRYLT1',ec550rhlt1_aer,pcols,lchnk)
     !         Since we do not have enough look-up table info to take abs550rhlt1_aer,
     !         instead take out abs550rh for each constituent:
     call outfld('ABSDRYBC',abs550rh_bc,pcols,lchnk)
     call outfld('ABSDRYOC',abs550rh_oc,pcols,lchnk)
     call outfld('ABSDRYSU',abs550rh_su,pcols,lchnk)
     call outfld('ABSDRYSS',abs550rh_ss,pcols,lchnk)
     call outfld('ABSDRYDU',abs550rh_du,pcols,lchnk)

  elseif(irf.ge.2) then   ! only happens for AEROCOM_INSITU

     irh=RF(irf)

     modeString="  "
     write(modeString,"(I2)"),irh
     if(RF(irf).eq.0) modeString="00"

     !-          varName = "EC44RH"//trim(modeString)
     !-          call outfld(varName,ec440rh_aer(:,:),pcols,lchnk)
     varName = "EC55RH"//trim(modeString)
     call outfld(varName,ec550rh_aer(:,:),pcols,lchnk)
     !-          varName = "EC87RH"//trim(modeString)
     !-          call outfld(varName,ec870rh_aer(:,:),pcols,lchnk)

     !-          varName = "AB44RH"//trim(modeString)
     !-          call outfld(varName,abs440rh_aer(:,:),pcols,lchnk)
     varName = "AB55RH"//trim(modeString)
     call outfld(varName,abs550rh_aer(:,:),pcols,lchnk)
     !-          varName = "AB87RH"//trim(modeString)
     !-          call outfld(varName,abs870rh_aer(:,:),pcols,lchnk)

  end if ! irf

end subroutine opticsAtConstRh

