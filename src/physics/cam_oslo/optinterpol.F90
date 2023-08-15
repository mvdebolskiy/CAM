module optinterpol

  !--------------------------------------------------------------------------------
  ! Purpose: To interpolate between look-up table entries for SW optical aerosol properties.
  ! Optimized for speed by Arild Burud and Egil Storen (NoSerC), June-July 2002
  ! Updated for new kcomp1.out including condensed SOA - Alf Kirkevaag, May 2013.
  ! Extended for new SOA treatment for  kcomp1-4.out and treating SOA as coagulated OC
  ! for kcomp5-10 - Alf Kirkevaag, August 2015, and also rewritten to a more generalized
  ! for for interpolations using common subroutines interpol*dim.
  !--------------------------------------------------------------------------------

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use ppgrid            , only : pcols, pver
  use commondefinitions , only : nmodes, nbmodes
  use opttab_lw         , only : ka0, ka1, ka2to3, ka4, ka5to10, nlwbands
  use opttab            , only : nbands, e
  use opttab            , only : fombg, fbcbg, cate, cat, fac, faq, fbc, rh, eps
  use opttab            , only : om0, g0, be0, ke0
  use opttab            , only : om1, g1, be1, ke1
  use opttab            , only : om2to3, g2to3, be2to3, ke2to3
  use opttab            , only : om4, g4, be4, ke4
  use opttab            , only : om5to10, g5to10, be5to10, ke5to10
  use lininterpol_mod   , only : lininterpol3dim, lininterpol4dim, lininterpol5dim

  implicit none
  private

  public :: inputForInterpol
  public :: interpol0
  public :: interpol1
  public :: interpol2to3
  public :: interpol4
  public :: interpol5to10

!********************************************************************************************
contains
!********************************************************************************************

  subroutine inputForInterpol (lchnk, ncol, rhum, xrh, irh1, &
       f_soana, xfombg, ifombg1, faitbc, xfbcbg, ifbcbg1,    &
       fnbc, xfbcbgn, ifbcbgn1, Nnatk, Cam, xct, ict1,       &
       focm, fcm, xfac, ifac1, fbcm, xfbc, ifbc1, faqm, xfaq, ifaq1)

    !
    ! Input arguments
    integer, intent(in)  :: lchnk                      ! chunk identifier
    integer, intent(in)  :: ncol                       ! number of atmospheric columns
    real(r8), intent(in) :: rhum(pcols,pver)           ! level relative humidity (fraction)
    real(r8), intent(in) :: f_soana(pcols,pver)        ! SOA/(SOA+H2SO4) mass fraction for the background in mode 1
    real(r8), intent(in) :: faitbc(pcols,pver)         ! BC/(BC + OC) mass fraction for the background in mode 4
    real(r8), intent(in) :: fnbc(pcols,pver)           ! BC/(BC + OC) mass fraction for the background in mode 14
    real(r8), intent(in) :: focm(pcols,pver,4)         ! fraction of added mass which is either SOA condensate or OC coagulate
    real(r8), intent(in) :: Cam(pcols,pver,nbmodes)    ! added internally mixed SO4+BC+OC concentration for a normalized mode
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! aerosol mode number concentration
    real(r8), intent(in) :: fcm(pcols,pver,nbmodes)    ! fraction of added mass which is either BC or OC/SOA (carbonaceous)
    real(r8), intent(in) :: fbcm(pcols,pver,nbmodes)   ! fraction of added mass as BC/(BC+OC)
    real(r8), intent(in) :: faqm(pcols,pver,nbmodes)   ! fraction of added sulfate which is from aqueous phase (ammonium sulfate)
    !
    ! Output arguments
    real(r8), intent(out) :: xrh(pcols,pver)           ! rhum for use in the interpolations
    integer,  intent(out) :: irh1(pcols,pver)
    real(r8), intent(out) :: xfombg(pcols,pver)        ! f_soana for use in the interpolations (mode 1)
    integer,  intent(out) :: ifombg1(pcols,pver)
    real(r8), intent(out) :: xfbcbg(pcols,pver)        ! faitbc for use in the interpolations (mode 4)
    integer,  intent(out) :: ifbcbg1(pcols,pver)
    real(r8), intent(out) :: xfbcbgn(pcols,pver)       ! fnbc for use in the interpolations (mode 14)
    integer,  intent(out) :: ifbcbgn1(pcols,pver)
    real(r8), intent(out) :: xct(pcols,pver,nmodes)    ! Cam/Nnatk for use in the interpolations
    integer,  intent(out) :: ict1(pcols,pver,nmodes)
    real(r8), intent(out) :: xfac(pcols,pver,nbmodes)  ! focm (1-4) or fcm (5-10) for use in the interpolations
    integer,  intent(out) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(out) :: xfbc(pcols,pver,nbmodes)  ! fbcm for use in the interpolations
    integer,  intent(out) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(out) :: xfaq(pcols,pver,nbmodes)  ! faqm for use in the interpolations
    integer,  intent(out) :: ifaq1(pcols,pver,nbmodes)
    !
    ! Local variables
    integer k, icol, i, irelh
    real(r8) :: eps10 = 1.e-10_r8
    !------------------------------------------------------------------------
    !
    ! write(*,*) 'Before xrh-loop'
    do k=1,pver
       do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
       end do
    end do

    ! write(*,*) 'Before rh-loop'
    do irelh=1,9
       do k=1,pver
          do icol=1,ncol
             if(xrh(icol,k) >= rh(irelh) .and. xrh(icol,k)<=rh(irelh+1)) then
                irh1(icol,k)=irelh
             endif
          end do
       end do
    end do
    ! write(*,*) 'xrh, irh1, irh2 =', xrh(1,26), irh1(1,26), irh2(1,26)

    do k=1,pver
       do icol=1,ncol
          ! find common xfombg, ifombg1 and ifombg2 for use in the interpolation routines
          xfombg(icol,k) =min(max(f_soana(icol,k),fombg(1)),fombg(6))
          ifombg1(icol,k)=int(5.0_r8*xfombg(icol,k)-eps10)+1 ! Boer linkes til def. i opttab.F90
       end do
    enddo

    do k=1,pver
       do icol=1,ncol
          ! find common xfbcbg, ifbcbg1 and ifbcbg2 for use in the interpolation routines
          xfbcbg(icol,k) =min(max(faitbc(icol,k),fbcbg(1)),fbcbg(6)) ! Boer linkes til def. i opttab.F90
          ifbcbg1(icol,k)=min(max(int(4*log10(xfbcbg(icol,k))+6),1),5)

          ! find common xfbcbgn, ifbcbgn1 and ifbcbgn2 for use in the interpolation routines
          xfbcbgn(icol,k) =min(max(fnbc(icol,k),fbcbg(1)),fbcbg(6)) ! Boer linkes til def. i opttab.F90
          ifbcbgn1(icol,k)=min(max(int(4*log10(xfbcbgn(icol,k))+6),1),5)
       end do
    enddo

    do i=1,4
       do k=1,pver
          do icol=1,ncol
             ! find common xfac, ifac1 and ifac2 for use in the interpolation routines
             xfac(icol,k,i) =min(max(focm(icol,k,i),fac(1)),fac(6))
             ifac1(icol,k,i)=int(5.0_r8*xfac(icol,k,i)-eps10)+1 ! Boer linkes til def. i opttab.F90
          end do
       enddo
    enddo
    do i=5,nbmodes
       do k=1,pver
          do icol=1,ncol
             ! find common xfac, ifac1 and ifac2 for use in the interpolation routines
             xfac(icol,k,i) =min(max(fcm(icol,k,i),fac(1)),fac(6))
             ifac1(icol,k,i)=int(5.0_r8*xfac(icol,k,i)-eps10)+1 ! Boer linkes til def. i opttab.F90
          end do
       enddo
    enddo

    do i=1,nbmodes
       do k=1,pver
          do icol=1,ncol
             ! find common xfbc, ifbc1 and ifbc2 for use in the interpolation routines
             xfbc(icol,k,i) =min(max(fbcm(icol,k,i),fbc(1)),fbc(6)) ! Boer linkes til def. i opttab.F90
             ifbc1(icol,k,i)=min(max(int(4*log10(xfbc(icol,k,i))+6),1),5)
          end do
       enddo
    enddo

    do i=1,nbmodes
       do k=1,pver
          do icol=1,ncol
             ! find common xfaq, ifaq1 and ifaq2 for use in the interpolation routines
             xfaq(icol,k,i) =min(max(faqm(icol,k,i),faq(1)),faq(6))
             ifaq1(icol,k,i)=int(5.0_r8*xfaq(icol,k,i)-eps10)+1 ! Boer linkes til def. i opttab.F90
          end do
       enddo
    enddo

    ! find common xct, ict1 and ict2 for use in the interpolation routines ! Boer linkes til def. i opttab.F90
    do i=1,4
       do k=1,pver
          do icol=1,ncol
             xct(icol,k,i)=min(max(Cam(icol,k,i)/(Nnatk(icol,k,i)+eps),cate(i,1)),cate(i,16))
             if(i.le.2) then
                ict1(icol,k,i)=min(max(int(3*log10(xct(icol,k,i))+19.666_r8),1),15)
             elseif(i.eq.3) then ! mode not used
                xct(icol,k,i)=cate(i,1)
                ict1(icol,k,i)=1
             else
                ict1(icol,k,i)=min(max(int(3*log10(xct(icol,k,i))+13.903_r8),1),15)
             endif
          end do
       end do
    end do

    do i=5,10
       do k=1,pver
          do icol=1,ncol
             xct(icol,k,i)=min(max(Cam(icol,k,i)/(Nnatk(icol,k,i)+eps),cat(i,1)),cat(i,6))
             if(i.eq.5) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.824_r8),1),5)
             elseif(i.eq.6) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.523_r8),1),5)
             elseif(i.eq.7) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.699_r8),1),5)
             elseif(i.eq.8) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+5.921_r8),1),5)
             elseif(i.eq.9) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.301_r8),1),5)
             else
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.699_r8),1),5)
             endif
          end do
       end do
    end do

    do i=11,nmodes ! for the externally mixed modes 11-14 (now only 12 and 14)
       do k=1,pver
          do icol=1,ncol
             xct(icol,k,i)=cate(i-10,1)
             ict1(icol,k,i)=1
          end do
       end do
    end do

    return

  end subroutine inputForInterpol

 !********************************************************************************************
  subroutine interpol0 (lchnk, ncol, daylight, Nnatk, omega, gass, bex, ske, lw_on, kabs)
    !
    ! Arguments
    integer  , intent(in)  :: lchnk                              ! chunk identifier
    integer  , intent(in)  :: ncol                               ! number of atmospheric columns
    logical  , intent(in)  :: daylight(pcols)                    ! calculations also at (polar) night if daylight=.true.
    logical  , intent(in)  :: lw_on                              ! LW calculations are performed if true
    real(r8) , intent(in)  :: Nnatk(pcols,pver,0:nmodes)         ! modal aerosol number concentration
    real(r8) , intent(out) :: omega(pcols,pver,0:nmodes,nbands)  ! spectral modal single scattering albedo
    real(r8) , intent(out) :: gass(pcols,pver,0:nmodes,nbands)   ! spectral modal asymmetry factor
    real(r8) , intent(out) :: bex(pcols,pver,0:nmodes,nbands)    ! spectral modal extinction coefficient
    real(r8) , intent(out) :: ske(pcols,pver,0:nmodes,nbands)    ! spectral modal specific extinction coefficient
    real(r8) , intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands) ! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer i, kcomp, k, icol
    !---------------------------------------

    kcomp=0
    do i=1,nbands
       do icol=1,ncol
          do k=1,pver
             omega(icol,k,kcomp,i)=0.0_r8
             gass(icol,k,kcomp,i)=0.0_r8
             bex(icol,k,kcomp,i)=0.0_r8
             ske(icol,k,kcomp,i)=0.0_r8
          end do
       end do
    end do
    do i=1,nlwbands
       do icol=1,ncol
          do k=1,pver
             kabs(icol,k,kcomp,i)=0.0_r8
          end do
       end do
    end do

    ! SW optical parameters

    do k=1,pver
       do icol=1,ncol
          ! if(Nnatk(icol,k,kcomp)>0.0_r8) then
          if(daylight(icol)) then
             do i=1,nbands   ! i = wavelength index
                omega(icol,k,kcomp,i)=om0(i)
                gass(icol,k,kcomp,i)=g0(i)
                bex(icol,k,kcomp,i)=be0(i)
                ske(icol,k,kcomp,i)=ke0(i)
             end do          ! i
          else  ! daylight
             ! Need be and ke in   nband=4 for lw calculation
             bex(icol,k,kcomp,4)=be0(4)
             ske(icol,k,kcomp,4)=ke0(4)
          end if ! daylight
       end do ! icol
    end do ! k

    ! LW optical parameters

    if(lw_on) then
       do k=1,pver
          do icol=1,ncol
             do i=1,nlwbands   ! i = wavelength index
                kabs(icol,k,kcomp,i)=ka0(i)
             end do            ! i
          end do ! icol
       end do ! k

    endif ! lw_on

  end subroutine interpol0

  !********************************************************************************************
  subroutine interpol1 (lchnk, ncol, daylight, xrh, irh1, mplus10, Nnatk, xfombg, ifombg1, &
       xct, ict1, xfac, ifac1, omega, gass, bex, ske, lw_on, kabs)

    !
    ! Arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xfombg(pcols,pver)         ! SOA/(SOA+H2SO4) for the background mode
    integer,  intent(in) :: ifombg1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)

    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absoption coefficient
    !
    ! Local variables
    integer i, kcomp, k, icol, kc10
    real(r8) a, b
    integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifo1, t_ifo2
    real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2
    real(r8) t_cat1, t_cat2, t_fombg1, t_fombg2, t_xfombg
    real(r8) d2mx(4), dxm1(4), invd(4)
    real(r8) opt4d(2,2,2,2)
    real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) kabs1, kabs2
    !---------------------------------------

    ! write(*,*) 'Before kcomp-loop'
    do kcomp=1,1

       if(mplus10==0) then
          kc10=kcomp
       else
          kc10=kcomp+10
       endif

       ! write(*,*) 'Before init-loop', kc10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kc10,i)=0.0_r8
                gass(icol,k,kc10,i)=0.0_r8
                bex(icol,k,kc10,i)=0.0_r8
                ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kcomp)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifo1 = ifombg1(icol,k)
             t_ifo2 = t_ifo1+1

             t_rh1  = rh(t_irh1)
             !x      t_rh2  = t_rh1+1
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_fombg1 = fombg(t_ifo1)
             t_fombg2 = fombg(t_ifo2)

             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kcomp)
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


             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! single scattering albedo:

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=om1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=om1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=om1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=om1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=om1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=om1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=om1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=om1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=om1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=om1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=om1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=om1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=om1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=om1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=om1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=om1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) /(t_rh2-t_rh1)
                   !alt       omega(icol,k,kc10,i)=(d2mx(1)*ome1+dxm1(1)*ome2)*invd(1)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! asymmetry factor

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=g1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=g1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=g1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=g1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=g1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=g1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=g1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=g1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=g1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=g1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=g1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=g1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=g1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=g1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=g1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=g1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, ge1, ge2)

                   ! finally, interpolation in the rh dimension (dim. 1)
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) /(t_rh2-t_rh1)
                   !alt      gass(icol,k,kc10,i)=(d2mx(1)*ge1+dxm1(1)*ge2)*invd(1)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol extinction

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, bex1, bex2)

                   bex1=max(bex1,1.e-30_r8)
                   bex2=max(bex2,1.e-30_r8)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) /(t_rh2-t_rh1)
                      !alt        bex(icol,k,kc10,i)=(d2mx(1)*bex1+dxm1(1)*bex2)*invd(1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                      !alt        a=(log(bex2)-log(bex1))*invd(1)
                      !alt        b=(t_rh2*log(bex1)-t_rh1*log(bex2))*invd(1)
                      !alt        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

                ! if(bex(icol,k,kc10,8)<1.e-20_r8) then
                ! write(*,995) 'bex(8)=', kc10, t_xrh, t_xct, t_xfac, t_xfombg, bex(icol,k,kc10,8)
                ! endif
             else  ! daylight


                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol extinction used for size information in LW

                i=4

                ! end points as basis for multidimentional linear interpolation
                opt4d(1,1,1,1)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                opt4d(1,1,1,2)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                opt4d(1,1,2,1)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                opt4d(1,1,2,2)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                opt4d(1,2,1,1)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                opt4d(1,2,1,2)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                opt4d(1,2,2,1)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                opt4d(1,2,2,2)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                opt4d(2,1,1,1)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                opt4d(2,1,1,2)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                opt4d(2,1,2,1)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                opt4d(2,1,2,2)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                opt4d(2,2,1,1)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                opt4d(2,2,1,2)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                opt4d(2,2,2,1)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                opt4d(2,2,2,2)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                ! interpolation in the fac, cat and fombg dimensions
                call lininterpol4dim (d2mx, dxm1, invd, opt4d, bex1, bex2)

                bex1=max(bex1,1.e-30_r8)
                bex2=max(bex2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                if(t_xrh <= 0.37_r8) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol specific extinction

                ! end points as basis for multidimentional linear interpolation
                opt4d(1,1,1,1)=ke1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                opt4d(1,1,1,2)=ke1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                opt4d(1,1,2,1)=ke1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                opt4d(1,1,2,2)=ke1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                opt4d(1,2,1,1)=ke1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                opt4d(1,2,1,2)=ke1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                opt4d(1,2,2,1)=ke1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                opt4d(1,2,2,2)=ke1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                opt4d(2,1,1,1)=ke1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                opt4d(2,1,1,2)=ke1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                opt4d(2,1,2,1)=ke1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                opt4d(2,1,2,2)=ke1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                opt4d(2,2,1,1)=ke1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                opt4d(2,2,1,2)=ke1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                opt4d(2,2,2,1)=ke1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                opt4d(2,2,2,2)=ke1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                ! interpolation in the fac, cat and fombg dimensions
                call lininterpol4dim (d2mx, dxm1, invd, opt4d, ske1, ske2)

                ske1=max(ske1,1.e-30_r8)
                ske2=max(ske2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37_r8) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                   !alt        ske(icol,k,kc10,i)=(d2mx(1)*ske1+dxm1(1)*ske2)*invd(1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                   !alt        a=(log(ske2)-log(ske1))*invd(1)
                   !alt        b=(t_rh2*log(ske1)-t_rh1*log(ske2))*invd(1)
                   !alt        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i

             if (lw_on) then

                ! LW optical parameters
                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol specific absorption in LW

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=ka1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=ka1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=ka1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=ka1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=ka1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=ka1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=ka1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=ka1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=ka1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=ka1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=ka1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=ka1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=ka1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=ka1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=ka1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=ka1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30)
                   kabs2=max(kabs2,1.e-30)

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37) then
                      kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             endif ! lw_on

          end do ! icol
       end do ! k

       ! write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
       ! write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
       ! write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
       ! write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)

    end do  ! kcomp

    return
  end subroutine interpol1


  !********************************************************************************************
  subroutine interpol2to3 (lchnk, ncol, daylight, xrh, irh1, mplus10, Nnatk, &
       xct, ict1, xfac, ifac1, omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)

    ! Output arguments
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer i, kcomp, k, icol, kc10
    real(r8) a, b
    integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2,t_cat1, t_cat2
    real(r8) d2mx(3), dxm1(3), invd(3)
    real(r8) opt3d(2,2,2)
    real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) kabs1, kabs2
    !---------------------------------------

    ! write(*,*) 'Before kcomp-loop'
    ! do kcomp=2,3
    do kcomp=2,2

       if(mplus10==0) then
          kc10=kcomp
       else
          kc10=kcomp+10
       endif

       ! write(*,*) 'Before init-loop', kc10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kc10,i)=0.0_r8
                gass(icol,k,kc10,i)=0.0_r8
                bex(icol,k,kc10,i)=0.0_r8
                ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1

             ! write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
             ! write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
             ! write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
             ! write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)

             ! write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
             ! write(*,*) 't_fac1,t_fac2=',t_fac1,t_fac2

             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)

             ! partial lengths along each dimension (1-4) for interpolation
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_cat2-t_xct)
             dxm1(2) = (t_xct-t_cat1)
             invd(2) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(3) = (t_fac2-t_xfac)
             dxm1(3) = (t_xfac-t_fac1)
             invd(3) = 1.0_r8/(t_fac2-t_fac1)


             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! single scattering albedo:

                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=om2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=om2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=om2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=om2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=om2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=om2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=om2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=om2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                        /(t_rh2-t_rh1)
                   ! write(*,*) omega(icol,k,kc10,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! asymmetry factor

                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=g2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=g2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=g2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=g2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=g2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=g2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=g2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=g2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, ge1, ge2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                        /(t_rh2-t_rh1)
                   ! write(*,*) gass(icol,k,kc10,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol extinction

                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=be2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=be2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=be2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=be2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=be2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=be2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=be2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=be2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, bex1, bex2)

                   bex1=max(bex1,1.e-30)
                   bex2=max(bex2,1.e-30)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             else  ! daylight



                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol extinction used for LW size information

                i=4
                ! end points as basis for multidimentional linear interpolation
                opt3d(1,1,1)=be2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                opt3d(1,1,2)=be2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                opt3d(1,2,1)=be2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                opt3d(1,2,2)=be2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                opt3d(2,1,1)=be2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                opt3d(2,1,2)=be2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                opt3d(2,2,1)=be2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                opt3d(2,2,2)=be2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                ! interpolation in the (fac and) cat dimension
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, bex1, bex2)

                bex1=max(bex1,1.e-30)
                bex2=max(bex2,1.e-30)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before bex'
                if(t_xrh <= 0.37) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight



             do i=4,4            ! i = wavelength index

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol specific extinction

                ! end points as basis for multidimentional linear interpolation
                opt3d(1,1,1)=ke2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                opt3d(1,1,2)=ke2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                opt3d(1,2,1)=ke2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                opt3d(1,2,2)=ke2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                opt3d(2,1,1)=ke2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                opt3d(2,1,2)=ke2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                opt3d(2,2,1)=ke2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                opt3d(2,2,2)=ke2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                ! interpolation in the (fac and) cat dimension
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, ske1, ske2)

                ske1=max(ske1,1.e-30)
                ske2=max(ske2,1.e-30)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i



             if (lw_on) then

                ! LW optical parameters
                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol specific absorption in LW

                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=ka2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=ka2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=ka2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=ka2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=ka2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=ka2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=ka2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=ka2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30_r8)
                   kabs2=max(kabs2,1.e-30_r8)

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37_r8) then
                      kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             endif ! lw_on

          end do ! icol
       end do ! k

       ! write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
       ! write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
       ! write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
       ! write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)

    end do  ! kcomp

    return
  end subroutine interpol2to3

  !********************************************************************************************

  subroutine interpol4 (lchnk, ncol, daylight, xrh, irh1, mplus10, Nnatk, xfbcbg, ifbcbg1, &
       xct, ict1, xfac, ifac1, xfaq, ifaq1, &
       omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xfbcbg(pcols,pver)         ! mass fraction BC/(BC+OC) for the background mode
    integer,  intent(in) :: ifbcbg1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Output arguments
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer  :: i, kcomp, k, kc10, icol
    real(r8) :: a, b
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2, t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq, t_fbcbg1, t_fbcbg2, t_xfbcbg, t_fac1
    real(r8) :: t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: kabs1, kabs2
    !---------------------------------------

    ! write(*,*) 'Before kcomp-loop'
    do kcomp=4,4

       if(mplus10==0) then
          kc10=kcomp
       else
          kc10=kcomp+10
       endif

       ! write(*,*) 'Before init-loop', kc10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kc10,i)=0.0_r8
                gass(icol,k,kc10,i)=0.0_r8
                bex(icol,k,kc10,i)=0.0_r8
                ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifb1 = ifbcbg1(icol,k)
             t_ifb2 = t_ifb1+1
             t_ifa1 = ifaq1(icol,k,kcomp)
             t_ifa2 = t_ifa1+1

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_fbcbg1 = fbcbg(t_ifb1)
             t_fbcbg2 = fbcbg(t_ifb2)
             t_faq1 = faq(t_ifa1)
             t_faq2 = faq(t_ifa2)

             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)
             t_xfbcbg = xfbcbg(icol,k)
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

             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! single scattering albedo:

                   opt5d(1,1,1,1,1)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) /(t_rh2-t_rh1)
                   ! write(*,*) omega(icol,k,kc10,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! asymmetry factor

                   opt5d(1,1,1,1,1)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ge1, ge2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                        /(t_rh2-t_rh1)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol extinction

                   opt5d(1,1,1,1,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                   bex1=max(bex1,1.e-30_r8)
                   bex2=max(bex2,1.e-30_r8)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             else  ! daylight

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol extinction called for use in size estimate for use in LW
                i=4

                opt5d(1,1,1,1,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,1,1,1,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,1,1,2,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,1,1,2,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,1,2,1,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,1,2,1,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,1,2,2,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,1,2,2,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(1,2,1,1,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,2,1,1,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,2,1,2,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,2,1,2,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,2,2,1,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,2,2,1,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,2,2,2,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,2,2,2,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,1,1,1,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,1,1,1,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,1,1,2,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,1,1,2,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,1,2,1,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,1,2,1,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,1,2,2,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,1,2,2,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,2,1,1,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,2,1,1,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,2,1,2,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,2,1,2,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,2,2,1,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,2,2,1,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,2,2,2,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,2,2,2,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac, cat and fbcbg dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                bex1=max(bex1,1.e-30_r8)
                bex2=max(bex2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before bex'
                if(t_xrh <= 0.37_r8) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol specific extinction

                opt5d(1,1,1,1,1)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,1,1,1,2)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,1,1,2,1)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,1,1,2,2)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,1,2,1,1)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,1,2,1,2)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,1,2,2,1)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,1,2,2,2)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(1,2,1,1,1)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,2,1,1,2)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,2,1,2,1)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,2,1,2,2)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,2,2,1,1)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,2,2,1,2)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,2,2,2,1)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,2,2,2,2)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,1,1,1,1)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,1,1,1,2)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,1,1,2,1)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,1,1,2,2)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,1,2,1,1)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,1,2,1,2)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,1,2,2,1)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,1,2,2,2)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,2,1,1,1)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,2,1,1,2)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,2,1,2,1)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,2,1,2,2)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,2,2,1,1)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,2,2,1,2)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,2,2,2,1)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,2,2,2,2)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac, cat and fbcbg dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, ske1, ske2)

                ske1=max(ske1,1.e-30_r8)
                ske2=max(ske2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37_r8) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i



             if (lw_on) then

                ! LW optical parameters

                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol specific absorption

                   opt5d(1,1,1,1,1)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30_r8)
                   kabs2=max(kabs2,1.e-30_r8)

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37_r8) then
                      kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             endif ! lw_on

          end do ! icol
       end do ! k

       ! write(*,*) 'kcomp, omega(1,26,kc10,4)=', kcomp, omega(1,26,kc10,4)
       ! write(*,*) 'kcomp, gass(1,26,kc10,4)=', kcomp, gass(1,26,kc10,4)
       ! write(*,*) 'kcomp, bex(1,26,kc10,4)=', kcomp, bex(1,26,kc10,4)
       ! write(*,*) 'kcomp, ske(1,26,kc10,4)=', kcomp, ske(1,26,kc10,4)

    end do  ! kcomp

  end subroutine interpol4

  !********************************************************************************************
  subroutine interpol5to10 (lchnk, ncol, daylight, xrh, irh1, Nnatk, xct, ict1, &
       xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
       omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! modal (OC+BC)/(SO4+BC+OC)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbc(pcols,pver,nbmodes)   ! modal BC/(OC+BC)
    integer,  intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Output arguments
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient

    ! Local variables
    integer  :: i, kcomp, k, icol
    real(r8) :: a, b
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2
    integer  :: t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq, t_fbc1, t_fbc2, t_xfbc, t_fac1
    real(r8) :: t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: kabs1, kabs2
    !---------------------------------------

    ! write(*,*) 'Before kcomp-loop'
    do kcomp=5,10

       ! write(*,*) 'Before init-loop', kcomp
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kcomp,i)=0.0_r8
                gass(icol,k,kcomp,i)=0.0_r8
                bex(icol,k,kcomp,i)=0.0_r8
                ske(icol,k,kcomp,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kcomp,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

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


             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! single scattering albedo:

                   opt5d(1,1,1,1,1)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kcomp,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                        /(t_rh2-t_rh1)
                   ! write(*,*) omega(icol,k,kcomp,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! asymmetry factor

                   opt5d(1,1,1,1,1)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ge1, ge2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kcomp,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                        /(t_rh2-t_rh1)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol extinction

                   opt5d(1,1,1,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                   bex1=max(bex1,1.e-30_r8)
                   bex2=max(bex2,1.e-30_r8)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             else  ! daylight


                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol extinction  used for aerosol size estimate needed for LW calculations
                i=4
                opt5d(1,1,1,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,1,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,1,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,1,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,1,2,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,2,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,2,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,2,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,1,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,1,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,1,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,1,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,2,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,2,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,2,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,2,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,1,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,1,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,1,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,1,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,2,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,2,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,2,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,2,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,1,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,1,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,1,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,1,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,2,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,2,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,2,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,2,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                ! interpolation in the faq, fbc, fac and cat dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                bex1=max(bex1,1.e-30_r8)
                bex2=max(bex2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before bex'
                if(t_xrh <= 0.37_r8) then
                   bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight



             do i=4,4            ! i = wavelength index

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol specific extinction

                opt5d(1,1,1,1,1)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,1,1,2)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,1,2,1)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,1,2,2)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,1,2,1,1)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,2,1,2)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,2,2,1)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,2,2,2)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,1,1,1)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,1,1,2)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,1,2,1)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,1,2,2)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,2,1,1)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,2,1,2)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,2,2,1)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,2,2,2)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,1,1,1)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,1,1,2)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,1,2,1)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,1,2,2)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,2,1,1)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,2,1,2)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,2,2,1)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,2,2,2)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,1,1,1)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,1,1,2)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,1,2,1)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,1,2,2)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,2,1,1)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,2,1,2)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,2,2,1)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,2,2,2)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                ! interpolation in the faq, fbc, fac and cat dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, ske1, ske2)

                ske1=max(ske1,1.e-30_r8)
                ske2=max(ske2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37_r8) then
                   ske(icol,k,kcomp,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kcomp,i)=e**(a*t_xrh+b)
                endif

             end do ! i



             if (lw_on) then

                ! LW optical parameters

                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                   ! aerosol specific absorption

                   opt5d(1,1,1,1,1)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30_r8)
                   kabs2=max(kabs2,1.e-30_r8)

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37_r8) then
                      kabs(icol,k,kcomp,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kcomp,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             endif ! lw_on

          end do ! icol
       end do ! k
    end do  ! kcomp

  end subroutine interpol5to10

end module optinterpol
