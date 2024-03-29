subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid
  ! files
  character(LEN=5)                    :: nsort, nocpu
  character(LEN = 80)                 :: filenamex,filenamey,filenamez
  integer::uleidx,uleidy,uleidz

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !  if(myid .EQ. 1) write(*,*) 'STATUS: starting cooling_fine, ilevel=', ilevel
  ! Valeska
  !--------------------------------------------------------------------
  ! The write option permits us to construct a column density map projected 
  ! in x, y, and z by calculating the column densities along positive and 
  ! negative directions with respect to a slice that passes through the 
  ! respective mid planes.
  !--------------------------------------------------------------------
  if(writing) then
     call title(ifout-1, nsort)
     call title(myid, nocpu)
     filenamex= TRIM(nsort)//'_test_densX_'//TRIM(nocpu)//'.dat'        !ex.:00001_test_densX_00010.dat
     filenamey= TRIM(nsort)//'_test_densY_'//TRIM(nocpu)//'.dat'
     filenamez= TRIM(nsort)//'_test_densZ_'//TRIM(nocpu)//'.dat'
     
     uleidx = myid + 100                                                !integer
     uleidy = myid + 200
     uleidz = myid + 300
     
     open(unit=uleidx, file=filenamex, form='formatted', status='unknown',position='append')
     open(unit=uleidy, file=filenamey, form='formatted', status='unknown',position='append')
     open(unit=uleidz, file=filenamez, form='formatted', status='unknown',position='append')
  end if
  !--------------------------------------------------------------------
  ! Valeska

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
#ifdef grackle
     if(use_grackle==0)then
        if(myid==1)write(*,*)'Computing new cooling table'
        call set_table(dble(aexp))
     endif
#else
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
#endif
  endif

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  use radiation_parameters, ONLY: sinks_opt_thin
#ifdef grackle
  use grackle_parameters
#endif
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
  use rt_parameters, only: nGroups, iGroups
  use rt_hydro_commons
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,kappaAbs, a_r,is_kIR_T,rt_vc
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
#if defined(grackle) && !defined(WITHOUTMPI)
  integer::info
#endif
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant=1.
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  logical,dimension(1:nvector),save::insink
#ifdef RT
  integer::ig,iNp,il
  real(kind=8),dimension(1:nvector),save:: ekk_new
  logical,dimension(1:nvector),save::cooling_on=.true.
  real(dp)::scale_Np,scale_Fp,work,Npc,Npnew, kIR, E_rad, TR
  real(dp),dimension(1:ndim)::Fpnew
  real(dp),dimension(nIons, 1:nvector),save:: xion
  real(dp),dimension(nGroups, 1:nvector),save:: Np, Np_boost=0d0, dNpdt=0d0
  real(dp),dimension(ndim, nGroups, 1:nvector),save:: Fp, Fp_boost=0d0, dFpdt=0d0
  real(dp),dimension(ndim, 1:nvector),save:: p_gas, u_gas
  real(kind=8)::f_trap, NIRtot, EIR_trapped, unit_tau, tau, Np2Ep, aexp_loc
  real(dp),dimension(nDim, nDim):: tEdd ! Eddington tensor
  real(dp),dimension(nDim):: flux
  real(kind=8),dimension(1:nvector),save::T2_new
#endif
#ifdef grackle
  real(kind=8),dimension(1:nvector),save::T2_new
#endif
#if NENER>0
  integer::irad
#endif

  real(dp) :: barotrop1D,mincolumn_dens
  real(dp)                                   :: x0, y0, z0,coeff_chi,cst2, coef
  double precision                           :: v_extinction
  integer::uleidx,uleidy,uleidz,uleidh,igrid,ii,indc2,iskip2,ind_ll


  !-------------- SPHERICAL DIRECTIONS ------------------------------------------------!
!  real(dp),dimension(1:nvector,1:ndir)                 :: col_dens                     !
  real(dp),dimension(1:nvector,1:NdirExt_m,1:NdirExt_n):: column_dens,column_dens_loc  !
  real(dp),dimension(1:nvector,1:NdirExt_m,1:NdirExt_n):: H2column_dens,H2column_dens_loc  ! H2 column density
!  real(dp),dimension(ndir)                             :: vcol_dens                    !
  real(dp),dimension(1:NdirExt_m,1:NdirExt_n)          :: vcolumn_dens
  real(dp),dimension(1:3)                              :: xpos
  real(dp)                                             :: dx_cross_int, dx_cross_loc
  integer                                              :: index_m,index_n, mmmm,nnnn
  integer                                              :: m, n, mloop, nloop, nl   
  integer, dimension(1:NdirExt_n)                      :: deltan1, deltan2
  integer                                              :: deltam
  !-----   simple_chem   --------------------------------------------------------------!

  !Valeska
!  vcol_dens(:) = 0.
  column_dens(:,:,:) = 0.
  H2column_dens(:,:,:) = 0.
  vcolumn_dens(:,:) = 0.
  
  if(writing) then
     !---position of reference to calculate the column density maps---
     ! we add 1.0D-09 in order to avoid the exact center (there are not cells centered in 0.5L)
     x0 = 0.5D0 + 1.0D-09
     y0 = 0.5D0 + 1.0D-09
     z0 = 0.5D0 + 1.0D-09
     
     !---Units uleidx, uleidy, uleidz---
     uleidx = myid + 100
     uleidy = myid + 200
     uleidz = myid + 300
  end if
  
294 FORMAT(I10,4ES14.5) 
295 FORMAT(5ES14.5)   
296 FORMAT(I10,5ES14.5)                
  
  if(numbtot(1,ilevel)==0)return
  
  !get the column density within the box from the grid faces

  !-----   EXTERNAL CONTRIBUTION   -----
  if(extinction)  call column_density(ind_grid,ngrid,ilevel,column_dens,H2column_dens) 
  !Valeska

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif

#ifdef RT
#if NGROUPS>0
  if(rt_isIRtrap) then
     ! For conversion from photon number density to photon energy density:
     Np2Ep = scale_Np * group_egy(iIR) * ev_to_erg                       &
           * rt_pressBoost / scale_d / scale_v**2
  endif
#endif
  aexp_loc=aexp
  ! Allow for high-z UV background in noncosmo sims:
  if(.not. cosmo .and. haardt_madau .and. aexp_ini .le. 1.)              &
       aexp_loc = aexp_ini
#endif

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf.eq.0)cycle

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do

     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

#ifdef RT
     ! Floor density (prone to go negative with strong rad. pressure):
     do i=1,nleaf
        uold(ind_leaf(i),1) = max(uold(ind_leaf(i),1),smallr)
     end do
     ! Initialise gas momentum and velocity for photon momentum abs.:
     do i=1,nleaf
        p_gas(:,i) = uold(ind_leaf(i),2:ndim+1) * scale_d * scale_v
        u_gas(:,i) = uold(ind_leaf(i),2:ndim+1) &
                     /uold(ind_leaf(i),1) * scale_v
     end do

#if NGROUPS>0
     if(rt_isIRtrap) then  ! Gather also trapped photons for solve_cooling
        iNp=iGroups(iIR)
        do i=1,nleaf
           il=ind_leaf(i)
           rtuold(il,iNp) = rtuold(il,iNp) + uold(il,iIRtrapVar)/Np2Ep
           if(rt_smooth) &
                rtunew(il,iNp)= rtunew(il,iNp) + uold(il,iIRtrapVar)/Np2Ep
        end do
     endif

     if(rt_vc) then      ! Add/remove radiation work on gas. Eq A6 in RT15
        iNp=iGroups(iIR)
        do i=1,nleaf
           il=ind_leaf(i)
           NIRtot = rtuold(il,iNp)
           kIR  = kappaSc(iIR)
           if(is_kIR_T) then                        ! kIR depends on T_rad
              ! For rad. temperature,  weigh the energy in each group by
              ! its opacity over IR opacity (derived from IR temperature)
              E_rad = group_egy(iIR) * ev_to_erg * NIRtot * scale_Np
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)     ! IR temp.
              kIR = kappaAbs(iIR) * (TR/10d0)**2
              do ig=1,nGroups
                 if(i .ne. iIR)                                          &
                      E_rad = E_rad + kappaAbs(ig) / kIR                 &
                            * max(rtuold(il,iGroups(ig)),smallNp)        &
                            * ev_to_erg * scale_Np
              end do
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)   ! Rad. temp.
              ! Set the IR opacity according to the rad. temperature:
              kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
           endif
           kIR = kIR*scale_d*scale_l           !  Convert to code units
           flux = rtuold(il,iNp+1:iNp+ndim)
           work = scale_v/c_cgs * kIR * sum(uold(il,2:ndim+1)*flux) &
                * Zsolar(i) * dtnew(ilevel)       ! Eq A6

           uold(il,ndim+2) = uold(il,ndim+2) &    ! Add work to gas energy
                + work * group_egy(iIR) &
                * ev_to_erg / scale_d / scale_v**2 / scale_l**3

           rtuold(il,iNp) = rtuold(il,iNp) - work !Remove from rad density
           rtuold(il,iNp) = max(rtuold(il,iNp),smallnp)
           call reduce_flux(rtuold(il,iNp+1:iNp+ndim),rtuold(il,iNp)*rt_c)
        enddo
     endif
#endif
#endif

     ! Compute thermal pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        err(i)=0.0d0
     end do
#if NENER>0
     do irad=0,nener-1
        do i=1,nleaf
           err(i)=err(i)+uold(ind_leaf(i),inener+irad)
        end do
     end do
#endif
     do i=1,nleaf
        emag(i)=0.0d0
     end do
#ifdef SOLVERmhd
     do idim=1,ndim
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+ndim+2)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
#endif
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i)-err(i)-emag(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2

!        if(ind .eq. 1 .and. abs(xg(ind_grid(i),1)-0.5) .le. dx .and. abs(xg(ind_grid(i),2)-0.5) .le. dx) then
!           write(*,*) 'T2, nH, scale_T2',T2(i),nH(i),scale_T2
!        do ii=0,nIons-1
!            write(*,*) 'iions', ii, uold(ind_leaf(i),iIons+ii)/uold(ind_leaf(i),1)
!        end do
!        endif

     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=MAX(exp(-nH(i)/0.01),1.0D-20)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = 0. !T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     if(cooling)then
        ! Compute thermal temperature by substracting polytrope
        do i=1,nleaf
           T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),T2max)
        end do
     endif

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

#ifdef RT
     if(neq_chem) then
        ! Get the ionization fractions
        do ii=0,nIons-1
           do i=1,nleaf
              xion(1+ii,i) = uold(ind_leaf(i),iIons+ii)/uold(ind_leaf(i),1)
           end do
        end do

        ! Get photon densities and flux magnitudes
        do ig=1,nGroups
           iNp=iGroups(ig)
           do i=1,nleaf
              il=ind_leaf(i)
              Np(ig,i)        = scale_Np * rtuold(il,iNp)
              Fp(1:ndim, ig, i) = scale_Fp * rtuold(il,iNp+1:iNp+ndim)
           enddo
           if(rt_smooth) then                           ! Smooth RT update
              do i=1,nleaf !Calc addition per sec to Np, Fp for current dt
                 il=ind_leaf(i)
                 Npnew = scale_Np * rtunew(il,iNp)
                 Fpnew = scale_Fp * rtunew(il,iNp+1:iNp+ndim)
                 dNpdt(ig,i)   = (Npnew - Np(ig,i)) / dtcool
                 dFpdt(:,ig,i) = (Fpnew - Fp(:,ig,i)) / dtcool
              end do
           end if
        end do

        if(cooling .and. delayed_cooling) then
           cooling_on(1:nleaf)=.true.
           do i=1,nleaf
              if(uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1) .gt. 1d-3) &
                   cooling_on(i)=.false.
           end do
        end if
        if(isothermal)cooling_on(1:nleaf)=.false.
     endif

     if(rt_vc) then ! Do the Lorentz boost. Eqs A4 and A5. in RT15
        do i=1,nleaf
           do ig=1,nGroups
              Npc=Np(ig,i)*rt_c_cgs
              call cmp_Eddington_tensor(Npc,Fp(:,ig,i),tEdd)
              Np_boost(ig,i) = - 2d0/c_cgs/rt_c_cgs * sum(u_gas(:,i)*Fp(:,ig,i))
              do idim=1,ndim
                 Fp_boost(idim,ig,i) =  &
                      -u_gas(idim,i)*Np(ig,i) * rt_c_cgs/c_cgs &
                      -sum(u_gas(:,i)*tEdd(idim,:))*Np(ig,i)*rt_c_cgs/c_cgs
              end do
           end do
           Np(:,i)   = Np(:,i) + Np_boost(:,i)
           Fp(:,:,i) = Fp(:,:,i) + Fp_boost(:,:,i)
        end do
     endif
#endif

     ! grackle tabular cooling
#ifdef grackle
     if(use_grackle==1)then
        gr_rank = 3
        do i = 1, gr_rank
           gr_dimension(i) = 1
           gr_start(i) = 0
           gr_end(i) = 0
        enddo
        gr_dimension(1) = nvector
        gr_end(1) = nleaf - 1

        if(cosmo)then
           my_grackle_units%a_value = aexp
           my_grackle_units%density_units = scale_d
           my_grackle_units%length_units = scale_l
           my_grackle_units%time_units = scale_t
           my_grackle_units%velocity_units = scale_v
        endif

        do i = 1, nleaf
           gr_density(i) = uold(ind_leaf(i),1)
           if(metal)then
              gr_metal_density(i) = uold(ind_leaf(i),imetal)
           else
              gr_metal_density(i) = uold(ind_leaf(i),1)*0.02*z_ave
           endif
           gr_energy(i) = T2(i)/(scale_T2*(gamma-1.0))
           gr_HI_density(i) = X*gr_density(i)
           gr_HeI_density(i) = (1.0-X)*gr_density(i)
           gr_DI_density(i) = 2.0*3.4e-5*gr_density(i)
        enddo
        ! Update grid properties
        my_grackle_fields%grid_rank = gr_rank
        my_grackle_fields%grid_dx = dx_loc

        iresult = solve_chemistry(my_grackle_units, my_grackle_fields, dtnew(ilevel))
        if(iresult.eq.0)then
            write(*,*) 'Grackle: error in solve_chemistry'
#ifndef WITHOUTMPI
            call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
            stop
#endif
        endif

        do i = 1, nleaf
           T2_new(i) = gr_energy(i)*scale_T2*(gamma-1.0)
        end do
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     else
        ! Compute net cooling at constant nH
        if(cooling.and..not.neq_chem)then
           call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)

        endif
     endif
#else
!     ! Compute net cooling at constant nH

!commented by PH to make it compatible with "FRIG" version 19/02/2017
!     if(cooling.and..not.neq_chem)then
!        call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
!     endif
!#endif
!#ifdef RT
!     if(neq_chem) then
!        T2_new(1:nleaf) = T2(1:nleaf)
!        call rt_solve_cooling(T2_new, xion, Np, Fp, p_gas, dNpdt, dFpdt  &
!                         , nH, cooling_on, Zsolar, dtcool, aexp_loc,nleaf)
!        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
!     endif
!#endif
!end of commented by PH

     if(cooling.and..not.neq_chem)then
!        call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
! USE Audit & Hennebelle cooling function
        call solve_cooling_frig(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)

!           write(*,*) 'solve cooling no RT'
     endif
#endif
#ifdef RT
     if(neq_chem) then
        T2_new(1:nleaf) = T2(1:nleaf)
        if(sinks_opt_thin)then
           do i=1,nleaf!!!raph 28/04/2020 for sink opt thin
              insink(i)=in_sink(ind_leaf(i))
           enddo
        endif
        call rt_solve_cooling(T2_new, xion, Np, Fp, p_gas, dNpdt, dFpdt  &
                         , nH, cooling_on, Zsolar, dtcool, aexp_loc,nleaf,insink)
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)

!     do i=1,nleaf
!        if(ind .eq. 1 .and. abs(xg(ind_grid(i),1)-0.5) .le. dx .and. abs(xg(ind_grid(i),2)-0.5) .le. dx) then
!           write(*,*) 'a T2, nH, scale_T2',T2(i),nH(i),scale_T2
!           write(*,*) 'a T2new, nH, scale_T2',T2_new(i),nH(i),scale_T2
!        endif
!     end do

     endif
#endif





#ifdef RT
     if(.not. static) then
        ! Update gas momentum and kinetic energy:
        do i=1,nleaf
           uold(ind_leaf(i),2:1+ndim) = p_gas(:,i) /scale_d /scale_v
        end do
        ! Energy update ==================================================
        ! Calculate NEW pressure from updated momentum
        ekk_new(1:nleaf) = 0d0
        do i=1,nleaf
           do idim=1,ndim
              ekk_new(i) = ekk_new(i) &
                   + 0.5*uold(ind_leaf(i),idim+1)**2 / uold(ind_leaf(i),1)
           end do
        end do
        do i=1,nleaf
           ! Update the pressure variable with the new kinetic energy:
           uold(ind_leaf(i),ndim+2) = uold(ind_leaf(i),ndim+2)           &
                                    - ekk(i) + ekk_new(i)
        end do
        do i=1,nleaf
           ekk(i)=ekk_new(i)
        end do

#if NGROUPS>0
        if(rt_vc) then ! Photon work: subtract from the IR ONLY radiation
           do i=1,nleaf
              Np(iIR,i) = Np(iIR,i) + (ekk(i) - ekk_new(i))              &
                   /scale_d/scale_v**2 / group_egy(iIR) / ev_to_erg
           end do
        endif
#endif
        ! End energy update ==============================================
     endif ! if(.not. static)
#endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Deal with cooling
     if(cooling.or.neq_chem)then
        ! Compute net energy sink
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Compute initial fluid internal energy
        do i=1,nleaf
           T2(i) = T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/max(uold(ind_leaf(i),1),smallr)
              if(cooling_switch > 1d-3)then
                 delta_T2(i) = MAX(delta_T2(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Compute polytrope internal energy
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0)
     end do

     ! Update fluid internal energy
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           T2(i) = T2(i) + delta_T2(i)
        end do
     endif

     ! Update total fluid energy
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     else if((cooling .or. neq_chem) .and. .not.fld)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=t_diss*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=max(uold(ind_leaf(i),idelay)*damp_factor,0d0)
        end do
     endif

#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
        do ii=0,nIons-1
           do i=1,nleaf
              uold(ind_leaf(i),iIons+ii) = xion(1+ii,i)*nH(i)
           end do
        end do
     endif
#if NGROUPS>0
     if(rt) then
        ! Update photon densities and flux magnitudes
        do ig=1,nGroups
           do i=1,nleaf
              rtuold(ind_leaf(i),iGroups(ig)) = (Np(ig,i)-Np_boost(ig,i)) /scale_Np
              rtuold(ind_leaf(i),iGroups(ig)) = &
                   max(rtuold(ind_leaf(i),iGroups(ig)),smallNp)
              rtuold(ind_leaf(i),iGroups(ig)+1:iGroups(ig)+ndim)         &
                               = (Fp(1:ndim,ig,i)-Fp_boost(1:ndim,ig,i)) /scale_Fp           
           enddo
        end do
     endif

     ! Split IR photons into trapped and freeflowing
     if(rt_isIRtrap) then
        if(nener .le. 0) then
           print*,'Trying to store E_trapped pressure, but NERAD too small!!'
           STOP
        endif
        iNp=iGroups(iIR)
        unit_tau = 1.5d0 * dx_loc * scale_d * scale_l
        do i=1,nleaf
           il=ind_leaf(i)
           NIRtot =max(rtuold(il,iNp),smallNp)      ! Total photon density
           kIR  = kappaSc(iIR)
           if(is_kIR_T) then                        ! kIR depends on T_rad
              ! For rad. temperature,  weigh the energy in each group by
              ! its opacity over IR opacity (derived from IR temperature)
              E_rad = group_egy(iIR) * ev_to_erg * NIRtot * scale_Np
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)     ! IR temp.
              kIR = kappaAbs(iIR) * (TR/10d0)**2
              do ig=1,nGroups
                 if(i .ne. iIR)                                          &
                      E_rad = E_rad + kappaAbs(ig) / kIR                 &
                            * max(rtuold(il,iGroups(ig)),smallNp)        &
                            * ev_to_erg * scale_Np
              end do
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)   ! Rad. temp.
              ! Set the IR opacity according to the rad. temperature:
              kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
           endif
           tau = nH(i) * Zsolar(i) * unit_tau * kIR
           f_trap = 0d0             ! Fraction IR photons that are trapped
           if(tau .gt. 0d0) f_trap = min(max(exp(-1d0/tau), 0d0), 1d0)
           ! Update streaming photons, trapped photons, and tot energy:
           rtuold(il,iNp) = max(smallnp,(1d0-f_trap) * NIRtot) ! Streaming
           rtuold(il,iNp+1:iNp+ndim) = &            ! Limit streaming flux
                                  rtuold(il,iNp+1:iNp+ndim) * (1d0-f_trap)
           EIR_trapped = max(0d0, NIRtot-rtuold(il,iNp)) * Np2Ep ! Trapped
           ! Update tot energy due to change in trapped radiation energy:
           uold(il,ndim+2)=uold(il,ndim+2)-uold(il,iIRtrapVar)+EIR_trapped
           ! Update the trapped photon energy:
           uold(il,iIRtrapVar) = EIR_trapped

           call reduce_flux(rtuold(il,iNp+1:iNp+ndim),rtuold(il,iNp)*rt_c)
        end do ! i=1,nleaf

     endif  !rt_isIRtrap
#endif
#endif

  end do
  ! End loop over cells

end subroutine coolfine1

#ifdef RT
!************************************************************************
subroutine cmp_Eddington_tensor(Npc,Fp,T_Edd)

! Compute Eddington tensor for given radiation variables
! Npc     => Photon number density times light speed
! Fp     => Photon number flux
! T_Edd  <= Returned Eddington tensor
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::Npc
  real(dp),dimension(1:ndim)::Fp ,u
  real(dp),dimension(1:ndim,1:ndim)::T_Edd
  real(dp)::iterm,oterm,Np_c_sq,Fp_sq,fred_sq,chi
  integer::p,q
!------------------------------------------------------------------------
  if(Npc .le. 0.d0) then
     write(*,*)'negative photon density in cmp_Eddington_tensor. -EXITING-'
     call clean_stop
  endif
  T_Edd(:,:) = 0.d0
  Np_c_sq = Npc**2
  Fp_sq = sum(Fp**2)              !  Sq. photon flux magnitude
  u(:) = 0.d0                           !           Flux unit vector
  if(Fp_sq .gt. 0.d0) u(:) = Fp/sqrt(Fp_sq)
  fred_sq = Fp_sq/Np_c_sq           !      Reduced flux, squared
  chi = max(4.d0-3.d0*fred_sq, 0.d0)   !           Eddington factor
  chi = (3.d0+ 4.d0*fred_sq)/(5.d0 + 2.d0*sqrt(chi))
  iterm = (1.d0-chi)/2.d0               !    Identity term in tensor
  oterm = (3.d0*chi-1.d0)/2.d0          !         Outer product term
  do p = 1, ndim
     do q = 1, ndim
        T_Edd(p,q) = oterm * u(p) * u(q)
     enddo
     T_Edd(p,p) = T_Edd(p,p) + iterm
  enddo

end subroutine cmp_Eddington_tensor
#endif

!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
subroutine pressure_eos(rho_temp,Enint_temp,Peos)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the pressure from the density and 
  ! internal volumic energy. Inputs/output are in code units
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Peos

  Peos = (gamma-1.d0)*Enint_temp

  return

end subroutine pressure_eos
!===========================================================================================
!===========================================================================================
!===========================================================================================
!===========================================================================================
subroutine temperature_eos(rho_temp,Enint_temp,Teos,ht)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  use cooling_module      ,only:kB,mH
  use radiation_parameters,only:mu_gas
  use units_commons
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the temperature from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  integer , intent(out):: ht 
  real(dp), intent(out):: Teos
  real(dp)::rho,Enint

  rho   = rho_temp*scale_d
  Enint = Enint_temp*scale_d*scale_v**2 

  Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))

  ht=1

  return

end subroutine temperature_eos
!===========================================================================================
!===========================================================================================
!===========================================================================================
!===========================================================================================
subroutine enerint_eos(rho_temp,temp_temp,Eeos)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  use cooling_module      ,only:kB,mH
  use radiation_parameters,only:mu_gas
  use units_commons
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the internal volumic energy from  
  ! the density and the temperature. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: temp_temp,rho_temp
  real(dp), intent(out):: Eeos
  real(dp)::rho,temp

  rho  = rho_temp * scale_d
  temp = temp_temp

  Eeos = rho*kB/(mu_gas*mH*(gamma-1.0))*temp/(scale_d*scale_v**2)

  return

end subroutine enerint_eos
!==================================================================================
!==================================================================================
!==================================================================================
!==================================================================================
subroutine soundspeed_eos(rho_temp,Enint_temp,Cseos)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the sound speed from the internal volumic energy 
  ! and the temperature. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Cseos

  Cseos = sqrt(gamma*(gamma-1.d0)*Enint_temp/rho_temp)

  return

end subroutine soundspeed_eos
!==================================================================================
!==================================================================================
!==================================================================================
!==================================================================================
function cmp_Cv_eos(rho,Enint)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  use cooling_module      ,only:kB,mH
  use radiation_parameters,only:mu_gas
  use units_commons
  implicit none
  !--------------------------------------------------------------
  ! This function computes the Cv from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp)   :: rho,Enint
  real(dp)   :: cmp_Cv_eos

  cmp_Cv_eos = rho*kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2

end function cmp_Cv_eos
!==================================================================================
!==================================================================================
!==================================================================================
!==================================================================================
double precision function barotrop1D(rhon)
  use hydro_commons
  use amr_parameters, only : n_star
  use radiation_parameters, only : Tr_floor
  implicit none

  real(dp)::inp,ll,rhon
  integer :: j

  if(analytical_barotrop)then
     barotrop1D = Tr_floor * ( 1.0d0 + (rhon/n_star)**(gamma-1.0d0) )
  else
     inp=rhon ! in g.cc
     ll=(1.d0+(log10(inp)-rhomin_barotrop)/drho_barotrop)
     j=dble(floor(ll))
     barotrop1D=(ll-j)*(temp_barotrop(j+1))+(1.d0-(ll-j))*(temp_barotrop(j))
     barotrop1D=10.0d0**barotrop1D     ! temperature in K
  endif

end function barotrop1D
!!!
!! all the rest has been moved
