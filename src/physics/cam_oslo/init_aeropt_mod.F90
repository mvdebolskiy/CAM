module aeropt_mod

  use shr_kind_mod      , only: r8 => shr_kind_r8
  use oslo_control      , only: oslo_getopts, dir_string_length
  use commondefinitions , only: nmodes, nbmodes
  use opttab            , only: cate, cat, fac, faq, fbc, rh, fombg, fbcbg
  use cam_logfile       , only: iulog

  implicit none
  private

  real(r8), public :: bep1    (38,10, 6, 16,   6       )
  real(r8), public :: bep2to3 (38,10,16,  6,   2:3     )
  real(r8), public :: bep4    (38,10, 6, 16,6, 6       )
  real(r8), public :: bep5to10(38,10, 6,  6,6, 6,  5:10)

  ! for initaeropt and intaeropt0:
  real(r8) , public:: bex440, bax440, bex500, bax500, bax550
  real(r8) , public:: bex670, bax670, bex870, bax870
  real(r8) , public:: bex550lt1, bex550gt1, backscx550

  public  :: initaeropt
  private :: set_data

contains

  subroutine initaeropt()

    !Purpose: To read in the AeroCom look-up tables for aerosol optical properties. 
    !  The grid for discrete input-values in the look-up tables is defined in opptab. 
    !  Tabulating the 'aerocomk'-files to save computing time.
    !  Updated for new kcomp1.out including condensed SOA - Alf Kirkevåg, May 2013
    !  Extended for new SOA treatment - Alf Kirkevaag, September 2015.
    !  Modified for optimized added masses and mass fractions for 
    !  concentrations from condensation, coagulation or cloud-processing 
    !  - Alf Kirkevaag, May 2016. 
    !  Modified for optimized added masses and mass fractions for concentrations from 
    !  condensation, coagulation or cloud-processing - Alf Kirkevaag, May 2016. 

    ! local variables
    integer  :: kcomp, irelh, ictot, ifac, ifbc, ifaq
    integer  :: ifombg, ifbcbg
    integer  :: ic, ifil, lin, iv
    character(len=dir_string_length) :: aerotab_table_dir
    !-----------------------------------------------------

    call oslo_getopts(aerotab_table_dir_out = aerotab_table_dir)

    open(11,file=trim(aerotab_table_dir)//'/aerocomk2.out' , form='formatted',status='old')
    open(12,file=trim(aerotab_table_dir)//'/aerocomk3.out' , form='formatted',status='old')
    open(13,file=trim(aerotab_table_dir)//'/aerocomk4.out' , form='formatted',status='old')
    open(14,file=trim(aerotab_table_dir)//'/aerocomk5.out' , form='formatted',status='old')
    open(15,file=trim(aerotab_table_dir)//'/aerocomk6.out' , form='formatted',status='old')
    open(16,file=trim(aerotab_table_dir)//'/aerocomk7.out' , form='formatted',status='old')
    open(17,file=trim(aerotab_table_dir)//'/aerocomk8.out' , form='formatted',status='old')
    open(18,file=trim(aerotab_table_dir)//'/aerocomk9.out' , form='formatted',status='old')
    open(19,file=trim(aerotab_table_dir)//'/aerocomk10.out', form='formatted',status='old')
    open(20,file=trim(aerotab_table_dir)//'/aerocomk0.out' , form='formatted',status='old')
    open(21,file=trim(aerotab_table_dir)//'/aerocomk1.out' , form='formatted',status='old')

    ! Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
    do ifil = 11,21
       call checkTableHeader (ifil)
    enddo
    !
    !-------------------------------------------
    ! Mode 0, BC(ax
    !-------------------------------------------
    !
    ifil = 11
    read (9+ifil,996) kcomp, relh,  &
         bex440, bax440, bex500, bax500, bax550, bex670, bax670, &
         bex870, bax870, bex550lt1, bex550gt1, backscx550
996 format(I2,f6.3,12e11.4)

    if (bex440<=0.0_r8) then
       write(*,*) 'bex440 =', bex440
       write(*,*) 'Error in initialization of bex1'
       stop
    endif
    write(iulog,*)'aerocom mode 0 ok'
    !
    !-------------------------------------------
    ! New mode 1 (H2SO4 and SOA + condensate from H2SO4 and SOA)
    !-------------------------------------------
    !
    ifil = 1
    do lin = 1,5760     ! 10x6x16x6
       call set_data(format_index=997, file_index=20+ifil, bep=bep1, mode1=.true.)
    end do  ! lin

    do irelh=1,10
       do ifombg=1,6
          do ictot=1,16
             do ifac=1,6
                if(bep1(1,irelh,ifombg,ictot,ifac)<=0.0_r8) then
                   write(*,*) 'bep1 =', irelh,ifombg, ictot, ifac, bep1(1,irelh,ifombg,ictot,ifac)
                   write(*,*) 'Error in initialization of bep1'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 1 ok' 
    ! 
    !-------------------------------------------
    ! Modes 2 to 3 (BC/OC + condesate from H2SO4 and SOA)
    !-------------------------------------------
    !
    ! do ifil = 2,3
    do ifil = 2,2
       do lin = 1,960     ! 10x16x6
          call set_data(format_index=994, file_index=9+ifil, bep=bep2to3, mode2to3=.true.)
       end do
    end do

    ! Prescribed dummy values for unused kcomp=3
    kcomp=3
    do irelh=1,10
       do ictot=1,16
          do ifac=1,6
             do iv=1,38
                bep2to3(iv,irelh,ictot,ifac,kcomp)=1.0_r8
             enddo
          enddo
       enddo
    enddo

    do kcomp=2,3
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                if(bep2to3(1,irelh,ictot,ifac,kcomp)<=0.0_r8) then
                   write(*,*) 'bep2to3 =', irelh, ictot, ifac, bep2to3(1,irelh,ictot,ifac,kcomp)
                   write(*,*) 'Error in initialization of bep2to3'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'aerocom mode 2-3 ok' 

    !
    !-------------------------------------------
    ! Mode 4 (BC&OC + condesate from H2SO4 and SOA + wetphase (NH4)2SO4)
    !-------------------------------------------
    !
    ifil = 4
    do lin = 1,34560     ! 10x16x6x6x6
       call set_data(format_index=995, file_index=9+ifil, bep=4, mode2to3=.true.)
    end do

    do irelh=1,10
       do ifbcbg=1,6
          do ictot=1,16
             do ifac=1,6
                do ifaq=1,6
                   if(bep4(1,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
                      write(*,*) 'bep4 =', irelh, ifbcbg, ictot, ifac, ifaq, bep4(1,irelh,ifbcbg,ictot,ifac,ifaq)
                      write(*,*) 'Error in initialization of bep4'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 4 ok'
    ! 
    !-------------------------------------------
    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
    !-------------------------------------------
    !
    do ifil = 5,10
       do lin = 1,12960     ! 10x6x6x6x6
          call set_data(format_index=993, file_index=9+ifil, bep=4, mode2to3=.true.)
       end do
    end do

    do kcomp=5,10
       do irelh=1,10
          do ictot=1,6
             do ifac=1,6
                do ifaq=1,6
                   if(bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
                      write(*,*) 'bep5to10 =', kcomp, irelh, ictot, ifac, ifbc, ifaq, &
                           bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)
                      write(*,*) 'Error in initialization of bep5to10'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 5-10 ok'

    do ifil=10,21
       close (ifil)
    end do

  end subroutine initaeropt

  subroutine set_data(this, format_index, file_index, bep, mode1, mode2to3, mode4, mode5to10)

    class(optical_properties_type) :: this
    integer  , intent(in)         :: format_index
    integer  , intent(in)         :: file_index
    real(r8) , intent(inout)      :: bep(:,:,:,:,:)
    logical, optional, intent(in) :: mode1
    logical, optional, intent(in) :: mode2to3
    logical, optional, intent(in) :: mode4
    logical, optional, intent(in) :: mode5to10

    ! local variables
    real(r8) :: catot, relh, frbcbg, frac, fabc, fraq   
    integer  :: kcomp, irelh, ictot, ifac, ifbc, ifaq
    integer  :: ifombg, ifbcbg
    real(r8) :: bext440, babs440, bext500, babs500, babs550    
    real(r8) :: bext670, babs670, bext870, babs870             
    real(r8) :: bebg440, babg440, bebg500, babg500, babg550    
    real(r8) :: bebg670, babg670, bebg870, babg870             
    real(r8) :: bebc440, babc440, bebc500, babc500, babc550    
    real(r8) :: bebc670, babc670, bebc870, babc870             
    real(r8) :: beoc440, baoc440, beoc500, baoc500, baoc550    
    real(r8) :: beoc670, baoc670, beoc870, baoc870             
    real(r8) :: besu440, basu440, besu500, basu500, basu550    
    real(r8) :: besu670, basu670, besu870, basu870             
    real(r8) :: bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1 
    real(r8) :: beoc550lt1, beoc550gt1, besu550lt1, besu550gt1 
    real(r8) :: backscat550
    real(r8) :: eps2 = 1.e-2_r8
    real(r8) :: eps4 = 1.e-4_r8
    real(r8) :: eps6 = 1.e-6_r8
    real(r8) :: eps7 = 1.e-7_r8
    !-----------------------------------------------------

    if  (format_index /= 993 .and. &
         format_index /= 994 .and. &
         format_index /= 995 .and. &
         format_index /= 996 .and. &
         format_index /= 997) then
       write(*,*) 'Error in format index'
       stop
    end if

    read (file_index,format_index) kcomp, relh, frombg, catot, frac, &
         bext440, bext500, bext670, bext870,             &
         bebg440, bebg500, bebg670, bebg870,             &
         bebc440, bebc500, bebc670, bebc870,             &
         beoc440, beoc500, beoc670, beoc870,             &
         besu440, besu500, besu670, besu870,             &
         babs440, babs500, babs550, babs670, babs870,    &
         bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
         beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
         backscat550, babg550, babc550, baoc550, basu550

    if (present(mode1)) then
       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             EXIT
          endif
       end do
       do ic=1,6
          if(abs(frombg-fombg(ic))<eps4) then
             ifombg=ic
             EXIT
          endif
       end do
       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             EXIT
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             EXIT
          endif
       end do
    else if (present(mode2to3)) then
       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             EXIT
          endif
       end do
       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             EXIT
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             EXIT
          endif
       end do
    else if (present(mode4)) then
       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
          endif
       end do
       do ic=1,6
          if(abs((frbcbg-fbcbg(ic))/fbcbg(ic))<eps2) then
             ifbcbg=ic
             exit
          endif
       end do
       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(fraq-faq(ic))<eps4) then
             ifaq=ic
             exit
          endif
       end do
    else if (present(mode5to10) then
       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs((catot-cat(kcomp,ic))/cat(kcomp,ic))<eps2) then
             ictot=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs((fabc-fbc(ic))/fbc(ic))<eps2) then
             ifbc=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(fraq-faq(ic))<eps4) then
             ifaq=ic
             exit
          endif
       end do
    end if

    bep(1,irelh,ifombg,ictot,ifac)  = this%bext440 ! unit km^-1
    bep(2,irelh,ifombg,ictot,ifac)  = this%bext500
    bep(3,irelh,ifombg,ictot,ifac)  = this%bext670 
    bep(4,irelh,ifombg,ictot,ifac)  = this%bext870
    bep(5,irelh,ifombg,ictot,ifac)  = this%bebg440
    bep(6,irelh,ifombg,ictot,ifac)  = this%bebg500
    bep(7,irelh,ifombg,ictot,ifac)  = this%bebg670
    bep(8,irelh,ifombg,ictot,ifac)  = this%bebg870
    bep(9,irelh,ifombg,ictot,ifac)  = this%bebc440  ! = 0
    bep(10,irelh,ifombg,ictot,ifac) = this%bebc500 ! = 0
    bep(11,irelh,ifombg,ictot,ifac) = this%bebc670 ! = 0
    bep(12,irelh,ifombg,ictot,ifac) = this%bebc870 ! = 0
    bep(13,irelh,ifombg,ictot,ifac) = this%beoc440
    bep(14,irelh,ifombg,ictot,ifac) = this%beoc500
    bep(15,irelh,ifombg,ictot,ifac) = this%beoc670
    bep(16,irelh,ifombg,ictot,ifac) = this%beoc870
    bep(17,irelh,ifombg,ictot,ifac) = this%besu440
    bep(18,irelh,ifombg,ictot,ifac) = this%besu500
    bep(19,irelh,ifombg,ictot,ifac) = this%besu670
    bep(20,irelh,ifombg,ictot,ifac) = this%besu870
    bep(21,irelh,ifombg,ictot,ifac) = this%babs440
    bep(22,irelh,ifombg,ictot,ifac) = this%babs500
    bep(23,irelh,ifombg,ictot,ifac) = this%babs550
    bep(24,irelh,ifombg,ictot,ifac) = this%babs670
    bep(25,irelh,ifombg,ictot,ifac) = this%babs870
    bep(26,irelh,ifombg,ictot,ifac) = this%bebg550lt1
    bep(27,irelh,ifombg,ictot,ifac) = this%bebg550gt1
    bep(28,irelh,ifombg,ictot,ifac) = this%bebc550lt1 ! = 0
    bep(29,irelh,ifombg,ictot,ifac) = this%bebc550gt1 ! = 0
    bep(30,irelh,ifombg,ictot,ifac) = this%beoc550lt1
    bep(31,irelh,ifombg,ictot,ifac) = this%beoc550gt1
    bep(32,irelh,ifombg,ictot,ifac) = this%besu550lt1
    bep(33,irelh,ifombg,ictot,ifac) = this%besu550gt1
    bep(34,irelh,ifombg,ictot,ifac) = this%backscat550
    bep(35,irelh,ifombg,ictot,ifac) = this%babg550
    bep(36,irelh,ifombg,ictot,ifac) = this%babc550 ! = this%0
    bep(37,irelh,ifombg,ictot,ifac) = this%baoc550
    bep(38,irelh,ifombg,ictot,ifac) = this%basu550

993 format(I2,f6.3,3e10.3,f5.2,38e10.3)
994 format(I2,f6.3,2e10.3,38e10.3)
995 format(I2,f6.3,3e10.3,f5.2,38e10.3)
996 format(I2,f6.3,12e11.4)
997 format(I2,f6.3,3e10.3,38e10.3)

  end subroutine set_data

end module aeropt_mod
