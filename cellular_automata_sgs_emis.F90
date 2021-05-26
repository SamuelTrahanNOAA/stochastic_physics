module cellular_automata_sgs_emis_mod

  implicit none

contains

  subroutine cellular_automata_sgs_emis(kstep,ugrs,vgrs,qgrs,pgr,vvl,prsl,vfrac_cpl, &
       ca_emis_anthro_cpl,ca_emis_dust_cpl,ca_emis_plume_cpl,ca_emis_seas_cpl, &
       ca_condition_diag, ca_plume_diag, ca_sgs_gbbepx_frp, domain_for_coupler, &
       nblks,isc,iec,jsc,jec,npx,npy,nlev,fhour,vegtype_cpl,iopt_dveg, &
       nca,ncells,nlives,nfracseed,nseed,nthresh,ca_global,ca_sgs,iseed_ca, &
       ca_smooth,nspinup,blocksize,cond_scale,mpiroot, mpicomm)

    use kinddef,           only: kind_phys
    use halo_exchange,     only: atmosphere_scalar_field_halo
    use update_ca,         only: update_cells_sgs, update_cells_global
    use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
    use mpp_domains_mod,   only: domain2D
    use block_control_mod, only: block_control_type, define_blocks_packed
    use mpi_wrapper,       only: mype,mp_reduce_sum,mp_bcst,mp_reduce_max,mp_reduce_min, &
         mpi_wrapper_initialize,is_master,mp_reduce_maxloc


    implicit none

    !L.Bengtsson, 2017-06

    !This program evolves a cellular automaton uniform over the globe given
    !the flag ca_global, if instead ca_sgs is .true. it evolves a cellular automata conditioned on
    !perturbed grid-box mean field. The perturbations to the mean field are given by a
    !stochastic gaussian skewed (SGS) distribution.

    !If ca_global is .true. it weighs the number of ca (nca) together to produce 1 output pattern
    !If instead ca_sgs is given, it produces nca ca:
    ! 1 CA_DEEP = deep convection
    ! 2 CA_SHAL = shallow convection
    ! 3 CA_TURB = turbulence

    !PLEASE NOTE: This is considered to be version 0 of the cellular automata code for FV3GFS, some functionally
    !is missing/limited.

    integer,intent(in) :: kstep,ncells,nca,nlives,nseed,iseed_ca,nspinup,mpiroot,mpicomm,iopt_dveg
    real(kind=kind_phys), intent(in)    :: nfracseed,nthresh,fhour
    logical,intent(in) :: ca_global, ca_sgs, ca_smooth
    integer, intent(in) :: nblks,isc,iec,jsc,jec,npx,npy,nlev,blocksize
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
    type(domain2D), intent(inout)       :: domain_for_coupler
    real(kind=kind_phys), intent(inout) :: cond_scale

    type(block_control_type)          :: Atm_block
    type(random_stat) :: rstate
    integer :: nlon, nlat, isize,jsize,nf,nn
    integer :: inci, incj, nxc, nyc, nxch, nych
    integer :: halo, k_in, i, j, k
    integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
    integer :: blocksz,levs,k350,k850
    integer(8) :: count, count_rate, count_max, count_trunc
    integer(8) :: iscale = 10000000000
    integer, allocatable :: iini(:,:,:),ilives(:,:,:),iini_g(:,:,:),ilives_g(:,:),ca_plumes(:,:), vegtype(:,:)
    real(kind=kind_phys), allocatable :: field_out(:,:,:), field_in(:,:),field_smooth(:,:),Detfield(:,:,:)
    real(kind=kind_phys), allocatable :: omega(:,:,:),pressure(:,:,:),cloud(:,:),humidity(:,:),uwind(:,:),vwind(:,:)
    real(kind=kind_phys), allocatable :: vertvelsum(:,:),vertvelmean(:,:),dp(:,:,:),surfp(:,:),shalp(:,:),gamt(:,:)
    real(kind=kind_phys), allocatable :: CA(:,:),condition(:,:),rho(:,:),conditiongrid(:,:)
    real(kind=kind_phys), allocatable :: CA_EMIS_ANTHRO(:,:),CA_EMIS_DUST(:,:)
    real(kind=kind_phys), allocatable :: CA_EMIS_PLUME(:,:),CA_EMIS_SEAS(:,:)
    real(kind=kind_phys), allocatable :: noise1D(:),vertvelhigh(:,:),noise(:,:,:),cond_save(:,:)
    real(kind=kind_phys) :: psum,csum,CAmean,sq_diff,CAstdv,count1,lambda
    real(kind=kind_phys) :: Detmax(nca),Detmin(nca),Detmean(nca),phi,stdev,delt
    logical,save         :: block_message=.true.
    logical              :: nca_plumes, init_condmax
    real(kind=kind_phys) :: init_weight

    !nca         :: switch for number of cellular automata to be used.
    !ca_global   :: switch for global cellular automata
    !ca_sgs      :: switch for cellular automata conditioned on SGS perturbed vertvel.
    !nfracseed   :: switch for number of random cells initially seeded
    !nlives      :: switch for maximum number of lives a cell can have
    !nspinup     :: switch for number of itterations to spin up the ca
    !ncells      :: switch for higher resolution grid e.g ncells=4
    !               gives 4x4 times the FV3 model grid resolution.
    !ca_smooth   :: switch to smooth the cellular automata
    !nthresh     :: threshold of perturbed vertical velocity used in case of sgs
    !nca_plumes   :: compute number of CA-cells ("plumes") within a NWP gridbox.

    ! Initialize MPI and OpenMP
    if (kstep==0) then
      call mpi_wrapper_initialize(mpiroot,mpicomm)
      return ! fields are not available yet.
    end if
    
    halo=1
    k_in=1

    if (nlev .EQ. 64) then
      k350=29
      k850=13
    elseif (nlev .EQ. 127) then
      k350=61
      k850=28
    else ! make a guess
      k350=int(nlev/2)
      k850=int(nlev/5)
      print*,'this level selection is not supported, making an approximation for k350 and k850'
    endif
    nca_plumes = .true.
    !----------------------------------------------------------------------------
    ! Get information about the compute domain, allocate fields on this
    ! domain

    ! Some security checks for namelist combinations:
    if(nca > 5)then
      write(0,*)'Namelist option nca cannot be larger than 5 - exiting'
      stop
    endif

    nlon=iec-isc+1
    nlat=jec-jsc+1
    isize=nlon+2*halo
    jsize=nlat+2*halo

    inci=ncells
    incj=ncells

    nxc=nlon*ncells
    nyc=nlat*ncells

    nxch=nxc+2*halo
    nych=nyc+2*halo

    !Allocate fields:
    allocate(vegtype(nlon,nlat))
    allocate(cloud(nlon,nlat))
    allocate(omega(nlon,nlat,nlev))
    allocate(pressure(nlon,nlat,nlev))
    allocate(humidity(nlon,nlat))
    allocate(uwind(nlon,nlat))
    allocate(vwind(nlon,nlat))
    allocate(dp(nlon,nlat,nlev))
    allocate(rho(nlon,nlat))
    allocate(surfp(nlon,nlat))
    allocate(vertvelmean(nlon,nlat))
    allocate(vertvelsum(nlon,nlat))
    allocate(field_in(nlon*nlat,1))
    allocate(field_out(isize,jsize,1))
    allocate(field_smooth(nlon,nlat))
    allocate(iini(nxc,nyc,nca))
    allocate(ilives(nxc,nyc,nca))
    allocate(iini_g(nxc,nyc,nca))
    allocate(ilives_g(nxc,nyc))
    allocate(vertvelhigh(nxc,nyc))
    allocate(condition(nxc,nyc))
    allocate(conditiongrid(nlon,nlat))
    allocate(shalp(nlon,nlat))
    allocate(gamt(nlon,nlat))
    allocate(Detfield(nlon,nlat,nca))
    allocate(CA(nlon,nlat))
    allocate(ca_plumes(nlon,nlat))
    allocate(CA_EMIS_ANTHRO(nlon,nlat))
    allocate(CA_EMIS_DUST(nlon,nlat))
    allocate(CA_EMIS_PLUME(nlon,nlat))
    allocate(CA_EMIS_SEAS(nlon,nlat))
    allocate(noise(nxc,nyc,nca))
    allocate(noise1D(nxc*nyc))
    if(cond_scale==0) then
      allocate(cond_save(nlon,nlat))
    endif

    !Initialize:
    Detfield(:,:,:)=0.
    vertvelmean(:,:) =0.
    vertvelsum(:,:)=0.
    cloud(:,:)=0.
    humidity(:,:)=0.
    uwind(:,:) = 0.
    condition(:,:)=0.
    conditiongrid(:,:)=0.
    vertvelhigh(:,:)=0.
    ca_plumes(:,:) = 0
    noise(:,:,:) = 0.0
    noise1D(:) = 0.0
    iini(:,:,:) = 0
    ilives(:,:,:) = 0
    iini_g(:,:,:) = 0
    ilives_g(:,:) = 0
    Detmax(:)=0.
    Detmin(:)=0.

    field_in=0
    field_out=0
    field_smooth=0

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

    field_out=0.
    
!    CA_EMIS_ANTHRO(:,:) = 0.0
!    CA_EMIS_DUST(:,:) = 0.0
!    CA_EMIS_PLUME(:,:) = 0.0
!    CA_EMIS_SEAS(:,:) = 0.0

    !Compute layer averaged vertical velocity (Pa/s)
    vertvelsum=0.
    vertvelmean=0.

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

    !Generate random number, following stochastic physics code:
    !if(kstep==2) then
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

      !Initiate the cellular automaton with random numbers larger than nfracseed

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
    !endif ! kstep=0

    !In case we want to condition the cellular automaton on a large scale field
    !we here set the "condition" variable to a different model field depending
    !on nf. (this is not used if ca_global = .true.)


    do nf=1,nca !update each ca
      
      if(nf==1)then
        call set_condition(ca_emis_plume_cpl,.true.)
      elseif(nf==2)then
        call set_condition(ca_emis_dust_cpl,.false.)
      elseif(nf==3) then
        call set_condition(ca_emis_anthro_cpl,.false.)
      else
        call set_condition(ca_emis_seas_cpl,.false.)
      endif

      !Calculate neighbours and update the automata
      !If ca-global is used, then nca independent CAs are called and weighted together to create one field; CA

      call update_cells_sgs(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,isc,iec,jsc,jec, &
           npx,npy,domain_for_coupler,CA,ca_plumes,iini,ilives,        &
           nlives,ncells,nfracseed,nseed,nthresh,nspinup,nf,nca_plumes)

      if(nf==1)then
        CA_EMIS_PLUME(:,:)=CA(:,:)/nlives
      elseif(nf==2)then
        CA_EMIS_DUST(:,:)=CA(:,:)/nlives
      elseif(nf==3) then
        CA_EMIS_ANTHRO(:,:)=CA(:,:)/nlives
      else
        CA_EMIS_SEAS(:,:)=CA(:,:)/nlives
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

    call normalize_output(CA_EMIS_PLUME,ca_emis_plume_cpl,.true.)
    call normalize_output(CA_EMIS_DUST,ca_emis_dust_cpl,.false.)
    call normalize_output(CA_EMIS_ANTHRO,ca_emis_anthro_cpl,.false.)
    call normalize_output(CA_EMIS_SEAS,ca_emis_seas_cpl,.false.)

    deallocate(omega)
    deallocate(pressure)
    deallocate(humidity)
    deallocate(dp)
    deallocate(conditiongrid)
    deallocate(shalp)
    deallocate(gamt)
    deallocate(rho)
    deallocate(surfp)
    deallocate(vertvelmean)
    deallocate(vertvelsum)
    deallocate(field_in)
    deallocate(field_out)
    deallocate(field_smooth)
    deallocate(iini)
    deallocate(ilives)
    deallocate(condition)
    deallocate(Detfield)
    deallocate(CA)
    deallocate(ca_plumes)
    deallocate(CA_EMIS_ANTHRO)
    deallocate(CA_EMIS_DUST)
    deallocate(CA_EMIS_PLUME)
    deallocate(CA_EMIS_SEAS)
    deallocate(noise)
    deallocate(noise1D)
    if(allocated(cond_save)) then
      deallocate(cond_save)
    endif
  contains

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
          field_in(i+(j-1)*nlon,1)=ca_in(blk,ix)
        enddo
      enddo

      call atmosphere_scalar_field_halo(field_out,halo,isize,jsize,k_in,field_in,isc,iec,jsc,jec,npx,npy,domain_for_coupler)

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

end module cellular_automata_sgs_emis_mod
