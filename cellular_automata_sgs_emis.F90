module cellular_automata_sgs_emis_mod

implicit none

contains

subroutine cellular_automata_emis_sgs(kstep,dtf,restart,first_time_step,domain, &
     ugrs,vgrs,qgrs,pgr,vvl,prsl,vfrac_cpl,fhour,vegtype_cpl,iopt_dveg, &
     ca_emis_anthro_cpl,ca_emis_dust_cpl,ca_emis_plume_cpl,ca_emis_seas_cpl, &
     ca_condition_diag, ca_plume_diag, ca_sgs_gbbepx_frp, &
     nblks,isc,iec,jsc,jec,npx,npy,nlev,nthresh,rcell, &
     nca,scells,tlives,nfracseed,nseed,ca_global,ca_sgs,iseed_ca, &
     ca_smooth,nspinup,ca_trigger,blocksize,mpiroot,mpicomm)

use kinddef,           only: kind_phys
use update_ca,         only: update_cells_sgs, update_cells_global, define_ca_domain
use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
use mpp_domains_mod,   only: domain2D
use block_control_mod, only: block_control_type, define_blocks_packed
use time_manager_mod, only: time_type
use mpi_wrapper,       only: mype,mp_reduce_sum,mp_bcst,mp_reduce_max,mp_reduce_min, &
                             mpi_wrapper_initialize
use mpp_domains_mod
use mpp_mod



implicit none
!L.Bengtsson, 2017-06

!L.Bengtsson, 2021-05
!Significant cleaning of old ideas
!Inclusion of restart capability
!Setting control variables as a function of dx and dt
!for scale adaptation. 

!This routine produces an output field CA_DEEP for coupling to convection (saSAS).
!CA_DEEP can be either number of plumes in a cluster (nca_plumes=true) or updraft 
!area fraction (nca_plumes=false)

integer,intent(in) :: kstep,scells,nca,tlives,nseed,iseed_ca,nspinup,mpiroot,mpicomm
real(kind=kind_phys), intent(in)    :: nfracseed,dtf,rcell,fhour
logical,intent(in) :: ca_global, ca_sgs, ca_smooth, restart,ca_trigger,first_time_step
integer, intent(in) :: nblks,isc,iec,jsc,jec,npx,npy,nlev,blocksize
real  , intent(out) :: nthresh
real(kind=kind_phys), intent(in)    :: ugrs(:,:,:)
real(kind=kind_phys), intent(in)    :: vgrs(:,:,:)
real(kind=kind_phys), intent(in)    :: qgrs(:,:,:)
real(kind=kind_phys), intent(in)    :: pgr(:,:)
real(kind=kind_phys), intent(in)    :: vvl(:,:,:)
real(kind=kind_phys), intent(in)    :: prsl(:,:,:)
real(kind=kind_phys), intent(inout) :: vfrac_cpl(:,:)
integer, intent(in) :: vegtype_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_emis_anthro_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_sgs_gbbepx_frp(:,:)
real(kind=kind_phys), intent(inout) :: ca_emis_dust_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_emis_plume_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_emis_seas_cpl(:,:)
real(kind=kind_phys), intent(out)   :: ca_condition_diag(:,:)
real(kind=kind_phys), intent(out)   :: ca_plume_diag(:,:)
type(domain2D),       intent(inout) :: domain

type(block_control_type)          :: Atm_block
type(random_stat) :: rstate
integer :: nlon, nlat, isize,jsize,nf,nn
integer :: inci, incj, nxc, nyc, nxch, nych, nx, ny
integer :: nxncells, nyncells
integer :: halo, k_in, i, j, k
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs
integer :: isdnx,iednx,jsdnx,jednx
integer :: iscnx,iecnx,jscnx,jecnx
integer :: ncells,nlives
integer, save :: initialize_ca
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000
real(kind=kind_phys), allocatable :: field_out(:,:,:),field_smooth(:,:]
real(kind=kind_phys), allocatable :: omega(:,:,:),pressure(:,:,:),humidity(:,:),uwind(:,:),vwind(:,:)
real(kind=kind_phys), allocatable :: vertvelsum(:,:),vertvelmean(:,:),dp(:,:,:),surfp(:,:)
real(kind=kind_phys), allocatable :: CA_EMIS_ANTHRO(:,:),CA_EMIS_DUST(:,:)
real(kind=kind_phys), allocatable :: CA_EMIS_PLUME(:,:),CA_EMIS_SEAS(:,:)
real(kind=kind_phys), allocatable :: vertvelhigh(:,:),cond_save(:,:)
integer, allocatable :: iini(:,:,:),ilives_in(:,:,:),ca_plumes(:,:)
real(kind=kind_phys), allocatable :: CA(:,:),condition(:,:),conditiongrid(:,:)
real(kind=kind_phys), allocatable :: noise1D(:),noise(:,:,:)
real(kind=kind_phys) :: condmax,livesmax,factor,dx,pi,re
type(domain2D)       :: domain_ncellx
logical,save         :: block_message=.true.
logical              :: nca_plumes
logical,save         :: first_flag

!nca         :: switch for number of cellular automata to be used.
!            :: for the moment only 1 CA can be used if ca_sgs = true
!ca_global   :: switch for global cellular automata
!ca_sgs      :: switch for cellular automata for deep convection
!nfracseed   :: switch for number of random cells initially seeded
!tlives      :: switch for time scale (s)
!nspinup     :: switch for number of itterations to spin up the ca
!scells      :: switch for CA cell size (m)
!ca_smooth   :: switch to smooth the cellular automata
!nca_plumes   :: compute number of CA-cells ("plumes") within a NWP gridbox.

! Initialize MPI and OpenMP
if (first_time_step) then
   call mpi_wrapper_initialize(mpiroot,mpicomm)
end if

halo=1
k_in=1

nca_plumes = .true.

if(first_time_step)then
   first_flag = .false.
   initialize_ca = 100000
endif

!----------------------------------------------------------------------------
! Get information about the compute domain, allocate fields on this
! domain

! Some security checks for namelist combinations:
 if(nca > 1)then
 write(0,*)'When ca_sgs=.True., nca has to be 1 - exiting'
 stop
 endif

 nlon=iec-isc+1
 nlat=jec-jsc+1
 isize=nlon+2*halo
 jsize=nlat+2*halo

 !Set time and length scales:
 call mpp_get_global_domain(domain,xsize=nx,ysize=ny,position=CENTER)
 pi=3.14159
 re=6371000.
 dx=0.5*pi*re/real(nx)
 ncells=int(dx/real(scells))
 nlives=int(real(tlives)/dtf)
 ncells = MIN(ncells,10)
 nlives = MAX(nlives,5)

 nthresh=rcell*real(ncells)*real(ncells)

 if(mype == 1)then
 write(*,*)'ncells=',ncells
 write(*,*)'nlives=',nlives
 write(*,*)'nthresh=',nthresh
 endif

 inci=ncells
 incj=ncells

!--- get params from domain_ncellx for building board and board_halo                                                                                

  !Get CA domain                                                                                                                                       
  call define_ca_domain(domain,domain_ncellx,ncells,nxncells,nyncells)
  call mpp_get_data_domain    (domain_ncellx,isdnx,iednx,jsdnx,jednx)
  call mpp_get_compute_domain (domain_ncellx,iscnx,iecnx,jscnx,jecnx)
  !write(1000+mpp_pe(),*) "nxncells,nyncells: ",nxncells,nyncells
  !write(1000+mpp_pe(),*) "iscnx,iecnx,jscnx,jecnx: ",iscnx,iecnx,jscnx,jecnx
  !write(1000+mpp_pe(),*) "isdnx,iednx,jsdnx,jednx: ",isdnx,iednx,jsdnx,jednx

  nxc = iecnx-iscnx+1
  nyc = jecnx-jscnx+1
  nxch = iednx-isdnx+1
  nych = jednx-jsdnx+1



 !Allocate fields:
 allocate(field_out(isize,jsize,1))
 allocate(field_smooth(nlon,nlat))
 allocate(omega(nlon,nlat,nlev))
 allocate(pressure(nlon,nlat,nlev))
 allocate(humidity(nlon,nlat))
 allocate(uwind(nlon,nlat))
 allocate(vwind(nlon,nlat))
 allocate(vertvelmean(nlon,nlat))
 allocate(vertvelsum(nlon,nlat))
 allocate(dp(nlon,nlat,nlev))
 allocate(surfp(nlon,nlat))
 allocate(CA_EMIS_ANTHRO(nlon,nlat))
 allocate(CA_EMIS_DUST(nlon,nlat))
 allocate(CA_EMIS_PLUME(nlon,nlat))
 allocate(CA_EMIS_SEAS(nlon,nlat))
 allocate(vertvelhigh(nxc,nyc))
 if(cond_scale==0) then
   allocate(cond_save(nlon,nlat))
 endif
 allocate(iini(nxc,nyc,nca))
 allocate(ilives_in(nxc,nyc,nca))
 allocate(ca_plumes(nlon,nlat))
 allocate(condition(nxc,nyc))
 allocate(conditiongrid(nlon,nlat))
 allocate(noise1D(nxc*nyc))
 allocate(noise(nxc,nyc,nca))

 !Initialize:
 humidity(:,:)=0.
 uwind(:,:) = 0.
 vwind(:,:) = 0.
 vertvelmean(:,:) =0.
 vertvelsum(:,:)=0.
 vertvelhigh(:,:)=0.
 condition(:,:)=0.
 conditiongrid(:,:)=0.
 ca_plumes(:,:) = 0
 noise(:,:,:) = 0.0
 noise1D(:) = 0.0
 iini(:,:,:) = 0
 ilives_in(:,:,:) = 0

 !Put the blocks of model fields into a 2d array - can't use nlev and blocksize directly,
 !because the arguments to define_blocks_packed are intent(inout) and not intent(in).
 levs=nlev
 blocksz=blocksize

 call define_blocks_packed('cellular_automata', Atm_block, isc, iec, jsc, jec, levs, &
                              blocksz, block_message)

    do blk = 1,Atm_block%nblks
      do ix = 1, Atm_block%blksz(blk)
        i = Atm_block%index(blk)%ii(ix) - isc + 1
        j = Atm_block%index(blk)%jj(ix) - jsc + 1
        uwind(i,j)         = ugrs(blk,ix,k350)
        vwind(i,j)         = vgrs(blk,ix,k350)
        conditiongrid(i,j) = max(0.0,vfrac_cpl(blk,ix))
        vegtype(i,j)       = vegtype_cpl(blk,ix)
        surfp(i,j)         = pgr(blk,ix)
        humidity(i,j)      = qgrs(blk,ix,k850) !about 850 hpa
        do k = 1,k350 !Lower troposphere
          omega(i,j,k)       = vvl(blk,ix,k) ! layer mean vertical velocity in pa/sec
          pressure(i,j,k)    = prsl(blk,ix,k) ! layer mean pressure in Pa
        enddo
      enddo
    enddo


    do blk = 1, Atm_block%nblks
      do ix = 1,Atm_block%blksz(blk)
        i = Atm_block%index(blk)%ii(ix) - isc + 1
        j = Atm_block%index(blk)%jj(ix) - jsc + 1
        CA_EMIS_ANTHRO(i,j)=ca_emis_anthro_cpl(blk,ix)*vfrac_cpl(blk,ix)
        CA_EMIS_DUST(i,j)=ca_emis_dust_cpl(blk,ix)*vfrac_cpl(blk,ix)
        CA_EMIS_PLUME(i,j)=ca_emis_plume_cpl(blk,ix)*vfrac_cpl(blk,ix)
        CA_EMIS_SEAS(i,j)=ca_emis_seas_cpl(blk,ix)*vfrac_cpl(blk,ix)
      enddo
    enddo


    do j=1,nlat
      do i =1,nlon
        dp(i,j,1)=(surfp(i,j)-pressure(i,j,1))
        do k=2,k350
          dp(i,j,k)=(pressure(i,j,k-1)-pressure(i,j,k))
        enddo
        count1=0.
        do k=1,k350
          count1=count1+1.
          vertvelsum(i,j)=vertvelsum(i,j)+(omega(i,j,k)*dp(i,j,k))
        enddo
      enddo
    enddo
    do j=1,nlat
      do i=1,nlon
        vertvelmean(i,j)=vertvelsum(i,j)/(surfp(i,j)-pressure(i,j,k350))
      enddo
    enddo

  condmax=maxval(condition)
  call mp_reduce_max(condmax)
  
if(kstep >=initialize_ca)then
  do nf=1,nca
     do j = 1,nyc
        do i = 1,nxc
           ilives_in(i,j,nf)=int(real(nlives)*(condition(i,j)/condmax))
        enddo
     enddo
  enddo

else

   do nf=1,nca
      do j = 1,nyc
         do i = 1,nxc
            ilives_in(i,j,nf)=0
         enddo
      enddo
   enddo

endif
                                                                                                                                        
!Generate random number, following stochastic physics code:
if(kstep == initialize_ca) then
   if (iseed_ca == 0) then
    ! generate a random seed from system clock and ens member number
    call system_clock(count, count_rate, count_max)
    ! iseed is elapsed time since unix epoch began (secs)
    ! truncate to 4 byte integer
    count_trunc = iscale*(count/iscale)
    count4 = count - count_trunc
  else
    ! don't rely on compiler to truncate integer(8) to integer(4) on
    ! overflow, do wrap around explicitly.
    count4 = mod(mype + iseed_ca + 2147483648, 4294967296) - 2147483648
  endif

  call random_setseed(count4)

  do nf=1,nca
    call random_number(noise1D)
    !Put on 2D:
    do j=1,nyc
      do i=1,nxc
        noise(i,j,nf)=noise1D(i+(j-1)*nxc)
      enddo
    enddo
   enddo

!Initiate the cellular automaton with random numbers larger than nfracseed
   do nf=1,nca
    do j = 1,nyc
      do i = 1,nxc
        if (noise(i,j,nf) > nfracseed ) then
          iini(i,j,nf)=1
        else
          iini(i,j,nf)=0
        endif
      enddo
    enddo
  enddo !nf

endif ! 

!Calculate neighbours and update the automata
 do nf=1,nca
   if(nf==1)then
     call set_condition(ca_emis_plume_cpl,.true.)
   elseif(nf==2)then
     call set_condition(ca_emis_dust_cpl,.false.)
   elseif(nf==3) then
     call set_condition(ca_emis_anthro_cpl,.false.)
   else
     call set_condition(ca_emis_seas_cpl,.false.)
   endif
   
   call update_cells_sgs(kstep,initialize_ca,first_flag,restart,first_time_step,iseed_ca,nca,nxc,nyc, &
                        nxch,nych,nlon,nlat,nxncells,nyncells,isc,iec,jsc,jec, &
                        npx,npy,isdnx,iednx,jsdnx,jednx,iscnx,iecnx,jscnx,jecnx,domain_ncellx,CA,ca_plumes,iini,ilives_in,        &
                        nlives,nfracseed,nseed,nspinup,nf,nca_plumes,ncells)

   livesmax=maxval(ilives_in)
   call mp_reduce_max(livesmax)
   
   if(nf==1)then
     CA_EMIS_PLUME(:,:)=CA(:,:)/livesmax
   elseif(nf==2)then
     CA_EMIS_DUST(:,:)=CA(:,:)/livesmax
   elseif(nf==3) then
     CA_EMIS_ANTHRO(:,:)=CA(:,:)/livesmax
   else
     CA_EMIS_SEAS(:,:)=CA(:,:)/livesmax
   endif
 enddo !nf (nca)

    !!Post-processesing - could be made into a separate sub-routine

    if(kstep == 1)then
      do j=1,nlat
        do i=1,nlon
          ca_plumes(i,j)=0.
        enddo
      enddo
    else
      do j=1,nlat
        do i=1,nlon
          if(conditiongrid(i,j) == 0)then
            ca_plumes(i,j)=0.
          endif
        enddo
      enddo
    endif

    !Put back into blocks 1D array to be passed to physics
    !or diagnostics output

    do blk = 1, Atm_block%nblks
      do ix = 1,Atm_block%blksz(blk)
        i = Atm_block%index(blk)%ii(ix) - isc + 1
        j = Atm_block%index(blk)%jj(ix) - jsc + 1

        ca_condition_diag(blk,ix)=conditiongrid(i,j)
        ca_plume_diag(blk,ix)=ca_plumes(i,j)

        ! ca_emis_anthro_cpl(blk,ix)=CA_EMIS_ANTHRO(i,j)/max(1.0,vfrac_cpl(blk,ix))
        ! ca_emis_dust_cpl(blk,ix)=CA_EMIS_DUST(i,j)/max(1.0,vfrac_cpl(blk,ix))
        ! ca_emis_plume_cpl(blk,ix)=CA_EMIS_PLUME(i,j)/max(1.0,vfrac_cpl(blk,ix))
        ! ca_emis_seas_cpl(blk,ix)=CA_EMIS_SEAS(i,j)/max(1.0,vfrac_cpl(blk,ix))
      enddo
    enddo

 deallocate(conditiongrid)
 deallocate(ssti)
 deallocate(lsmski)
 deallocate(lakei)
 deallocate(iini)
 deallocate(ilives_in)
 deallocate(condition)
 deallocate(CA)
 deallocate(ca_plumes)
 deallocate(CA_DEEP)
 deallocate(noise)
 deallocate(noise1D)


    subroutine normalize_output(ca_in,ca_out,save_condition)
      implicit none
      real(kind=kind_phys), intent(inout) :: ca_out(:,:)
      real(kind=kind_phys), intent(in) :: ca_in(:,:)
      logical, intent(in) :: save_condition
      integer :: blk,ix,i,j
      real(kind=kind_phys) :: minca,maxca,div,condmax,scale_at_condmax,scale
      real(kind=kind_phys) :: sendbuf(2)
      
      minca=1e20
      maxca=-1e20
      do blk = 1, Atm_block%nblks
        do ix = 1,Atm_block%blksz(blk)
          i = Atm_block%index(blk)%ii(ix) - isc + 1
          j = Atm_block%index(blk)%jj(ix) - jsc + 1
          ca_out(blk,ix)=ca_in(i,j) ! /max(1.0,vfrac_cpl(blk,ix))
          minca=min(minca,ca_out(blk,ix))
          maxca=max(maxca,ca_out(blk,ix))
        enddo
      enddo

      call mp_reduce_max(maxca)
      call mp_reduce_min(minca)

      div=1.0
      if(minca/=maxca) then
        div=maxca-minca
      endif

      if(cond_scale==0 .or. .not. save_condition .or. &
           (save_condition .and. allocated(cond_save))) then
        scale = 1.0
      else
        scale = cond_scale
      endif

      do blk = 1, Atm_block%nblks
        do ix = 1,Atm_block%blksz(blk)
          if(ca_out(blk,ix)/=0) then
            ca_out(blk,ix) = (ca_out(blk,ix)-minca)/div * scale
          endif
        enddo
      enddo

      if(save_condition .and. allocated(cond_save)) then
        ! Find conditiongrid/ca_out at maximum conditiongrid value,
        ! before conditiongrid was smoothed.
        condmax = 0
        scale_at_condmax = 0
        do blk = 1, Atm_block%nblks
          do ix = 1,Atm_block%blksz(blk)
            i = Atm_block%index(blk)%ii(ix) - isc + 1
            j = Atm_block%index(blk)%jj(ix) - jsc + 1
            if(cond_save(i,j)>condmax .and. ca_out(blk,ix)>0) then
              condmax = cond_save(i,j)
              scale_at_condmax = cond_save(i,j)/ca_out(blk,ix)
            endif
          enddo
        enddo

        ! Same as above, but for whole domain.
        sendbuf = (/ condmax, scale_at_condmax /)
        call mp_reduce_maxloc(sendbuf)
        cond_scale = max(0.0,sendbuf(2))

        ! Now that we have a scale, re-normalize
        do blk = 1, Atm_block%nblks
          do ix = 1,Atm_block%blksz(blk)
            if(ca_out(blk,ix)/=0) then
              ca_out(blk,ix) = ca_out(blk,ix) * cond_scale
            endif
          enddo
        enddo
      endif
    end subroutine normalize_output
    
    subroutine set_condition(ca_in,save_condition)
      implicit none
      real(kind=kind_phys), intent(in) :: ca_in(:,:)
      integer :: blk,ix,i,j,ih,jh,inci,incj
      logical, intent(in) :: save_condition
      real(kind=kind_phys) :: condmax
      
      init_weight=max(0.0,min(1.0,fhour))
      conditiongrid = 0

      ! if(init_weight>0.0) then
      !   do blk = 1,Atm_block%nblks
      !     do ix = 1, Atm_block%blksz(blk)
      !       i = Atm_block%index(blk)%ii(ix) - isc + 1
      !       j = Atm_block%index(blk)%jj(ix) - jsc + 1
      !       field_in(i+(j-1)*nlon,1)=ca_sgs_gbbepx_frp(blk,ix)*(1.0-init_weight)
      !     enddo
      !   enddo
      ! else
      !   field_in=0.0
      ! endif

      do blk = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(blk)
          i = Atm_block%index(blk)%ii(ix) - isc + 1
          j = Atm_block%index(blk)%jj(ix) - jsc + 1
          ih=i+halo
          jh=j+halo
          field_out(ih,jh,1) = ca_in(blk,ix)
        enddo
      enddo

      call atmosphere_scalar_field_halo(field_out,halo,isize,jsize,k_in,isc,iec,jsc,jec,npx,npy,domain_for_coupler)

      condmax=0
      do blk = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(blk)
          i = Atm_block%index(blk)%ii(ix) - isc + 1
          j = Atm_block%index(blk)%jj(ix) - jsc + 1
          ih=i+halo
          jh=j+halo
          field_smooth(i,j)=(8.0*field_out(ih,jh,1)+4.0*field_out(ih-1,jh,1)+ &
               4.0*field_out(ih,jh-1,1)+4.0*field_out(ih+1,jh,1)+&
               4.0*field_out(ih,jh+1,1)+2.0*field_out(ih-1,jh-1,1)+&
               2.0*field_out(ih-1,jh+1,1)+2.0*field_out(ih+1,jh+1,1)+&
               2.0*field_out(ih+1,jh-1,1))/32.
          conditiongrid(i,j) = max(0.0,vfrac_cpl(blk,ix)*(field_smooth(i,j)*init_weight*1000.0 + ca_sgs_gbbepx_frp(blk,ix)*(1.0-init_weight)))
          condmax = max(condmax,conditiongrid(i,j))
        enddo
      enddo

      if(save_condition .and. allocated(cond_save)) then
        cond_save = conditiongrid
      endif

      call mp_reduce_max(condmax)

      if(condmax>0) then
        do j=1,nlat
          do i=1,nlon
            conditiongrid(i,j) = conditiongrid(i,j)/condmax
          enddo
        enddo
      endif
      
      if(save_condition) then
        do blk = 1, Atm_block%nblks
          do ix = 1,Atm_block%blksz(blk)
            i = Atm_block%index(blk)%ii(ix) - isc + 1
            j = Atm_block%index(blk)%jj(ix) - jsc + 1

            ca_condition_diag(blk,ix)=conditiongrid(i,j)
          enddo
        enddo
      endif
      
      inci=ncells
      incj=ncells
      do j=1,nyc
        do i=1,nxc
          ilives(i,j,nf)=real(nlives)*conditiongrid(inci/ncells,incj/ncells)
          if(i.eq.inci)then
            inci=inci+ncells
          endif
        enddo
        inci=ncells
        if(j.eq.incj)then
          incj=incj+ncells
        endif
      enddo

      !Vertical velocity has its own variable in order to condition on combination
      !of "condition" and vertical velocity.

      inci=ncells
      incj=ncells
      do j=1,nyc
        do i=1,nxc
          vertvelhigh(i,j)=vertvelmean(inci/ncells,incj/ncells)
          if(i.eq.inci)then
            inci=inci+ncells
          endif
        enddo
        inci=ncells
        if(j.eq.incj)then
          incj=incj+ncells
        endif
      enddo
    end subroutine set_condition

  end subroutine cellular_automata_sgs_emis

end subroutine cellular_automata_sgs_emis

end module cellular_automata_sgs_emis_mod
