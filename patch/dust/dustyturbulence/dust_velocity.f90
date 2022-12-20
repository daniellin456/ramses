subroutine set_vdust(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the dust diffusion scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust,idim
  integer,dimension(1:nvector),save::ind_grid
  logical:: d_cycle_ok
  integer :: icycle, ncycle

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel


  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call vdustfine1(ind_grid,ngrid,ilevel)
  end do

  

  !if(simple_boundary)call make_boundary_dust(ilevel)
  if(simple_boundary)call make_boundary_hydro(ilevel)

111 format('   Entering dust_velocity for level ',i2)


end subroutine set_vdust
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine vdustfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use cooling_module,ONLY:kB,mH
  use cloud_module
  use radiation_parameters
  use units_commons, only : scale_m


  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !---------------------------------------------------------------------!
  ! This routine gathers first hydro variables from neighboring grids  -!
  ! to set initial conditions in a 6x6x6 grid. It interpolate from     -!
  ! coarser level missing grid variables. It then calls the            -!
  ! dust diffusion solver that compute the flux. This flux is zeroed at-!
  ! coarse-fine boundaries, since contribution from finer levels has   -!
  ! already been taken into account. Conservative variables are updated-!
  ! and stored in array unew(:), both at the current level and at the  -!
  ! coarser level if necessary.                                        -!
  !---------------------------------------------------------------------!

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,0:twondim  ,1:ndust*ndim),save::u1dust
  real(dp),dimension(1:nvector,1:twotondim,1:ndust*ndim),save::u2dust
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim),save::vloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim),save::vdloc

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer ,dimension(1:nvector                   ),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::idust,ht
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::  ncycle,icycle
  real(dp):: dt_dustcycle
  logical :: d_cycle_ok
  real(dp)::dx,scale,oneontwotondim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::sum_dust,sum_dust_new,sum_dust_old
  real(dp)::d,u,v,w,A,B,C,enint,e_kin,e_mag,pressure,cs, temp
  real(dp)::rho_gas, pi, t_stop,t_stop_floor,dens_floor,d0,r0

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
    
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do idim=1,ndim
              do idust=1,ndust
                 do i=1,nbuffer
                    u1dust(i,j,ndim*(idust-1)+idim)=v_dust(ibuffer_father(i,j),idust,idim)
                 end do
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
        call interpol_hydro_dust(u1dust,ind1,u2dust,nbuffer)

     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather dust velocity
        do idust=1,ndust
           do idim=1,ndim
              do i=1,nexist
                 vdloc(ind_exist(i),i3,j3,k3,idust,idim)=v_dust(ind_cell(i),idust,idim)
              end do
              do i=1,nbuffer
                 vdloc(ind_nexist(i),i3,j3,k3,idust,idim)=u2dust(i,ind_son,ndim*(idust-1)+idim)
              end do
           end do
     end do
 
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids


  call cmpvdust(uloc,vloc,vdloc,dx,dx,dx,dtnew(ilevel),ncache)

   !--------------------------------------------------------
   !Udate at level ilevel for the dust velocity
   !--------------------------------------------------------
  do idim=1,ndim
     do k2=k2min,k2max
        do j2=j2min,j2max
           do i2=i2min,i2max
              ind_son=1+i2+2*j2+4*k2
              iskip=ncoarse+(ind_son-1)*ngridmax
              do i=1,ncache
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              i3=1+i2
              j3=1+j2
              k3=1+k2
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    do idust=1,ndust
                       v_dust(ind_cell(i),idust,idim)=vloc(i,i3,j3,k3,idust,idim)
                    enddo
                 end if
              end do
           end do
        end do
     end do
end do


 
end subroutine vdustfine1


!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine cmpvdust(uin,vout,vdin,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use hydro_commons
  use units_commons
  use const
  use cloud_module
  use cooling_module,ONLY:kB,mH
  use radiation_parameters

  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim)::vout
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim)::vdin
  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::qin
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 
  ! declare local variables
  integer ::i, j, k, l, ht
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::idust,idim,irad
  real(dp) :: d,u,v,w,e_mag,e_kin, sum_dust, enint
  real(dp) ::pressure, cs,A,B,C,wnorm,vmax, Mach_dv,impl,bcell , tcell
  real(dp) ::  dAy, dAz,dBx,dBz,dCx,dCy,divB
  real(dp),dimension(1:ndim) :: fpress,fmag,grad_P,curl_BcrossB
  real(dp),dimension(1:ndust)  :: t_stop,t_gyr
  real(dp)  ::pi,tstop_tot,t_stop_floor, epsilon_0
  real(dp), dimension(1:ndust) ::d_grain,l_grain
  real(dp),dimension(1:ndim):: ew
  real(dp),dimension(1:ndust):: dustMRN
  real(dp),dimension(1:ndust,1:ndim):: wdust_Hydro
  
#if NIMHD==1
  real(dp):: betaadbricolo,etaohmdissbricolo, betaad,etaohmdiss,eta_hall_chimie, etahall,ionisrate,e_el,B_gauss
  real(dp),dimension(1:ndust):: Gamma_dust
  real(dp),dimension(1:ndim) :: curl_B,curl_BcrossBcrossB,Enimhd,b_dir,Enimhd_cross_B,Enimhd_par,Enimhd_ortho
  real(dp),dimension(1:ndust,1:ndim):: wdust_Hydro_cross_B,wdust_Hydro_par,wdust_Hydro_ortho,w_mh,w_em

#endif  
  epsilon_0 = dust_ratio(1)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  
  vmax=vdust_max/scale_v
  pi =3.14159265358979323846_dp
#if NIMHD==1  
  e_el = 1.602e-20
#endif
#if NDUST>0
     do idust =1,ndust
        dustMRN(idust) = dust_ratio(idust)/(1.0d0+dust_ratio(idust))
     end do     
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
     sum_dust=0.0d0
     do idust =1,ndust
        sum_dust = sum_dust + dustMRN(idust)
     end do
#endif
  if(mrn.eqv..true.) then
     call size_dust(l_grain)
     do idust=1,ndust
       l_grain(idust) = l_grain(idust)/scale_l
       d_grain(idust)=grain_dens(idust)/scale_d
    end do
  else
     do idust=1,ndust
       d_grain(idust)=grain_dens(idust)/scale_d
       l_grain(idust)=grain_size(idust)/scale_l

    end do
 endif
 call ctoprimdust(uin,qin,dt,ngrid)

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
                 d= max(qin(l,i,j,k,1),smallr)
                 tstop_tot=0.0d0
                 t_stop=0.0d0
                 sum_dust=0.0d0
                 do idust = 1,ndust
                    sum_dust= sum_dust + qin(l,i,j,k,firstindex_ndust+idust)
                 end do
                 !magnetic field and required derivatives to get rotB
                 A=0.5d0*(uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))
                 B=0.5d0*(uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))                 
                 C=0.5d0*(uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))
                 e_mag=0.5d0*(A**2+B**2+C**2)

                 u=0.0d0; v=0.0d0; w=0.0d0

                 u = qin(l,i,j,k,2)
                 if(ndim>1)v = qin(l,i,j,k,3)
                 if(ndim>2)w = qin(l,i,j,k,4)
                 e_kin=0.5d0*d*(u**2+v**2+w**2)

#if NENER>0
                 do irad=1,nener
                    e_mag=e_mag+uin(l,i,j,k,8+irad)
                 end do
#endif
                 if(energy_fix)then
                    enint=uin(l,i,j,k,nvar)
                 else
                    enint=uin(l,i,j,k,5)-e_kin- e_mag
                 end if
                 
                 call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)

                 do idust = 1,ndust
                    if(.not.charged_dust_dynamics) then
                       t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/max((d- uin(l,i,j,k,firstindex_ndust+idust)),smallr)
#if NDUSTPSCAL>0
                       if(stepinski)  t_stop(idust) =  d_grain(idust)*(uin(l,i,j,k,firstindex_dustpscal+idust)/uin(l,i,j,k,firstindex_ndust+idust))*SQRT(pi*gamma/8.0_dp)/cs/max((d- uin(l,i,j,k,firstindex_ndust+idust)),smallr)
                       
#endif
                    else
                       !No back reaction here
                       t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d
#if NDUSTPSCAL>0
                       if(stepinski)  t_stop(idust) =  d_grain(idust)*(uin(l,i,j,k,firstindex_dustpscal+idust)/uin(l,i,j,k,firstindex_ndust+idust))*SQRT(pi*gamma/8.0_dp)/cs/d
                       
#endif
                    endif

                    !! Carefull need to change with crossing time
!!$                    if(reduce_tstop)  t_stop(idust) = min(t_stop(idust),Stokes_max*sqrt(3.0*pi/32.0d0/max(d,smallr))) ! Stopping time limited with maximum stokes /!\ this is for the collapse patch

                    
                    if(K_drag)  t_stop(idust) = uin(l,i,j,k,firstindex_ndust+idust)/K_dust(idust)                     ! Linear drag
                    if(dust_diffusion_test) t_stop(idust) = K_dust(idust)/(1.d0-qin(l,i,j,k,firstindex_ndust+idust))  ! For Barrenblatt test

                    tstop_tot= tstop_tot-t_stop(idust)*qin(l,i,j,k,firstindex_ndust+idust)
                    
                    !Kwok correction
                    if(kwok_correction)then
                       wnorm= sqrt(sum(vdin(l,i,j,k,idust,1:ndim)**2.0))
                       Mach_dv = wnorm / cs/max((1.0d0-uin(l,i,j,k,firstindex_ndust+idust)/d),smallr/d)
                       t_stop(idust) =  t_stop(idust) /(1.0d0+(9.0d0*pi*Mach_dv**2/128.0d0))**0.5
                    end if
                 end do
                 
                 if(charged_dust_dynamics) tstop_tot=0.0d0
                 
                 dAy=(0.5d0*(uin(l,i,j+1,k,6)+uin(l,i,j+1,k,nvar+1))-0.5d0*(uin(l,i,j-1,k,6)+uin(l,i,j-1,k,nvar+1)))*0.5d0/dy
                 dAz=(0.5d0*(uin(l,i,j,k+1,6)+uin(l,i,j,k+1,nvar+1))-0.5d0*(uin(l,i,j,k-1,6)+uin(l,i,j,k-1,nvar+1)))*0.5d0/dz
                 dBx=(0.5d0*(uin(l,i+1,j,k,7)+uin(l,i+1,j,k,nvar+2))-0.5d0*(uin(l,i-1,j,k,7)+uin(l,i-1,j,k,nvar+2)))*0.5d0/dx
                 dBz=(0.5d0*(uin(l,i,j,k+1,7)+uin(l,i,j,k+1,nvar+2))-0.5d0*(uin(l,i,j,k-1,7)+uin(l,i,j,k-1,nvar+2)))*0.5d0/dz
                 dCy=(0.5d0*(uin(l,i,j+1,k,8)+uin(l,i,j+1,k,nvar+3))-0.5d0*(uin(l,i,j-1,k,8)+uin(l,i,j-1,k,nvar+3)))*0.5d0/dy
                 dCx=(0.5d0*(uin(l,i+1,j,k,8)+uin(l,i+1,j,k,nvar+3))-0.5d0*(uin(l,i-1,j,k,8)+uin(l,i-1,j,k,nvar+3)))*0.5d0/dx
  
                 !pressure gradient and magnetic curl
                 do idim=1,ndim
                    if(idim==1) then
                       grad_P(idim)=(qin(l,i+1,j,k,5)-qin(l,i-1,j,k,5))/(2.0d0*dx)
                       Curl_BcrossB(idim) = ((dAz-dCx)*C-(dBx-dAy)*B)
                    endif
                    if(idim==2) then
                       grad_P(idim)=(qin(l,i,j+1,k,5)-qin(l,i,j-1,k,5))/(2.0d0*dx)
                       Curl_BcrossB(idim)= ((dBx-dAy)*A-(dCy-dBz)*C)
                    end if
                    if(idim==3) then
                       grad_P(idim)=(qin(l,i,j,k+1,5)-qin(l,i,j,k-1,5))/(2.0d0*dx)
                       Curl_BcrossB(idim) = ((dCy-dBz)*B-(dAz-dCx)*A)
                    endif 
                 end do
                 do idust=1,ndust
                    t_stop(idust) = t_stop(idust)+tstop_tot
                 end do

                 do idust = 1,ndust
                    do idim =1,ndim
                       wdust_Hydro(idust,idim) =  t_stop(idust)*(grad_P(idim)-curl_BcrossB(idim))/d
                    end do            
                 end do
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !Charged grains modif
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
#if NIMHD==1
                 if(charged_dust_dynamics) then
                    bcell = (A**2+B**2+C**2)
                    
                    call temperature_eos(d*(1.0d0-sum_dust),enint,tcell,ht)
                    !If charged grains one must determine Enimhd
                    !Curlbbb for Ambipolar diffusion and CurlB for ohmic dissipation
                    Curl_B(1) = dCy-dBz
                    Curl_B(2) = dAz-dCx
                    Curl_B(3) = dBx-dAy
                    
                    curl_BcrossBcrossB(1)= curl_BcrossB(2)*C-curl_BcrossB(3)*B
                    curl_BcrossBcrossB(2)= curl_BcrossB(3)*A-curl_BcrossB(1)*C
                    curl_BcrossBcrossB(3)= curl_BcrossB(1)*B-curl_BcrossB(2)*A

                    b_dir(1) = A/sqrt(bcell)
                    b_dir(2) = B/sqrt(bcell)
                    b_dir(3) = C/sqrt(bcell)
                    
                    ionisrate=default_ionisrate
                    Enimhd=0.0d0 
                    !Ohmic dissipation
                    if(nmagdiffu .eq.1 .or.nmagdiffu2 .eq.1) then
                       etaohmdiss= etaohmdissbricolo(d,bcell,tcell,dt,dx,ionisrate)
                       Enimhd=Enimhd+etaohmdiss*curl_B
                    endif
                    !Ambipolar resistivity
                    if(nambipolar .eq. 1.or.nambipolar2 .eq. 1) then
                       betaad = betaadbricolo(d,d,dt,bcell,bcell,dx,0,tcell,ionisrate)
                       Enimhd=Enimhd+betaad*curl_BcrossBcrossB
                    endif
#if HALL==1
                    !Hall effect
                    if(nhall.eq.1) then
                       !TODO  
                    end if
#endif
                    
                    Enimhd_cross_B(1)=Enimhd(2)*b_dir(3)-Enimhd(3)*b_dir(2)
                    Enimhd_cross_B(2)=Enimhd(3)*b_dir(1)-Enimhd(1)*b_dir(3)
                    Enimhd_cross_B(3)=Enimhd(1)*b_dir(2)-Enimhd(2)*b_dir(1)
                    
                    Enimhd_par(1)=(Enimhd(1)*b_dir(1)+Enimhd(2)*b_dir(2)+Enimhd(3)*b_dir(3))*b_dir(1)
                    Enimhd_par(2)=(Enimhd(1)*b_dir(1)+Enimhd(2)*b_dir(2)+Enimhd(3)*b_dir(3))*b_dir(2)
                    Enimhd_par(3)=(Enimhd(1)*b_dir(1)+Enimhd(2)*b_dir(2)+Enimhd(3)*b_dir(3))*b_dir(3)
                    
                    do idust=1,ndust
                       
                       wdust_Hydro_cross_B(idust,1)=wdust_hydro(idust,2)*b_dir(3)-wdust_hydro(idust,3)*b_dir(2)
                       wdust_Hydro_cross_B(idust,2)=wdust_hydro(idust,3)*b_dir(1)-wdust_hydro(idust,1)*b_dir(3)
                       wdust_Hydro_cross_B(idust,3)=wdust_hydro(idust,1)*b_dir(2)-wdust_hydro(idust,2)*b_dir(1)
                       
                       wdust_Hydro_par(idust,1)=(wdust_hydro(idust,1)*b_dir(1)+wdust_hydro(idust,2)*b_dir(2)+wdust_hydro(idust,3)*b_dir(3))*b_dir(1)
                       wdust_Hydro_par(idust,2)=(wdust_hydro(idust,1)*b_dir(1)+wdust_hydro(idust,2)*b_dir(2)+wdust_hydro(idust,3)*b_dir(3))*b_dir(2)
                       wdust_Hydro_par(idust,3)=(wdust_hydro(idust,1)*b_dir(1)+wdust_hydro(idust,2)*b_dir(2)+wdust_hydro(idust,3)*b_dir(3))*b_dir(3)
                    end do
                    Enimhd_ortho= Enimhd-Enimhd_par
                    wdust_Hydro_ortho=wdust_Hydro-wdust_Hydro_par

                    ! We compute tgyr in cgs to avoid making units errors and rescale it
                    B_gauss= 4.0*pi*sqrt(Bcell*scale_m*scale_l**2/scale_t**2.)
                    do idust=1,ndust
                       t_gyr(idust)=min(Z_dust(idust)*e_el*B_gauss/(4./3.*pi*d_grain(idust)*l_grain(idust)**3.*scale_m)/scale_t,1e15*t_stop(idust))
#if NDUSTPSCAL>0
                       if(stepinski) t_gyr(idust)=min(Z_dust(idust)*e_el*B_gauss/(4./3.*pi*d_grain(idust)*(uin(l,i,j,k,firstindex_dustpscal+idust)/uin(l,i,j,k,firstindex_ndust+idust))**3.*scale_m)/scale_t,1e15*t_stop(idust))
#endif                      
                    end do
                    gamma_dust=t_stop/t_gyr

                    do idust=1,ndust
                       do idim=1,ndim
                          w_em(idust,idim)=1.0d0/(sqrt(Bcell))*(gamma_dust(idust)**2./(1.0d0+gamma_dust(idust)**2.)*Enimhd_cross_B(idim)+gamma_dust(idust)/(1.0d0+gamma_dust(idust)**2.)*Enimhd_par(idim))

                          w_mh(idust,idim) = gamma_dust(idust)**2./(1.0d0+gamma_dust(idust)**2.)*wdust_Hydro_cross_B(idust,idim)+ gamma_dust(idust)/(1.0d0+gamma_dust(idust)**2.)*wdust_Hydro_ortho(idust,idim) + wdust_Hydro_par(idust,idim)

                       end do
                    end do
                 endif
#endif

                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !End of charged grains modif
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               
              if(.not.charged_dust_dynamics) then
                 do idust = 1,ndust
                    do idim =1,ndim
                       vout(l,i,j,k,idust,idim)= wdust_Hydro(idust,idim) 
                    end do
                 end do
              else
                  do idust = 1,ndust
                     do idim =1,ndim
                        vout(l,i,j,k,idust,idim)= w_em(idust,idim)+w_mh(idust,idim)
                     end do
                  end do
               end if
               
               !Regularisation to avoid unphysical velocities (or annoyingly small dt :()
               if (reduce_wdust) then   
                  do idust=1,ndust
                     wnorm= sqrt(sum(vout(l,i,j,k,idust,1:ndim)**2.0))
                       do idim=1,ndim
                          if(vmax_barycenter)then
                             vmax = f_vmax*sqrt(u*u+v*v+w*w)
                          end if
                          if(vmax_cs)then
                             vmax = f_vmax*cs
                          end if
                          if(vmax_dust_lim)then
                             vmax = vdust_max/scale_l*scale_t
                          end if
                          if(wnorm>vmax) ew(idim)= vout(l,i,j,k,idust,idim)/wnorm        
                          if(wnorm>vmax) vout(l,i,j,k,idust,idim)=  ew(idim)*min(wnorm,vmax)
                       end do
                    end do
                 end if
              end do
           end do
        end do
     end do
     
     
end subroutine cmpvdust

subroutine ctoprimdust(uin,q,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  use radiation_parameters,only:small_er
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  

  integer ::i, j, k, l, idim
  real(dp)::eint, smalle, smallp, etot
  real(dp),dimension(1:nvector),save::eken,emag,erad

  ! EOS
  real(dp)  :: pp_eos

#if NENER>0
  integer::irad
#endif
#if NVAR>8+NENER
  integer::n
#endif
  real(dp):: sum_dust
#if NDUST>0
  integer:: idust,n_d
#endif  
  
  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma

 

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do
           ! Debug
           if(debug)then
              do l = 1, ngrid
                 if(uin(l,i,j,k,1).le.smallr)then
                    write(*,*)'negative density'
                    write(*,*)uin(l,i,j,k,1)
                    stop
                 end if
              end do
           end if

           ! Compute velocities
           do l = 1, ngrid
              q(l,i,j,k,2) = uin(l,i,j,k,2)/q(l,i,j,k,1)
              q(l,i,j,k,3) = uin(l,i,j,k,3)/q(l,i,j,k,1)
              q(l,i,j,k,4) = uin(l,i,j,k,4)/q(l,i,j,k,1)
           end do

           ! Compute cell centered magnetic field
           DO l = 1, ngrid
              q(l,i,j,k,6) = (uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))*half
              q(l,i,j,k,7) = (uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))*half
              q(l,i,j,k,8) = (uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))*half
           END DO

           ! Compute specific kinetic energy and magnetic energy
           do l = 1, ngrid
              eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
              emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
           end do

           ! Compute non-thermal pressure
           erad = zero
#if NENER>0
           do irad = 1,nent
              do l = 1, ngrid
                 q(l,i,j,k,8+irad) = (gamma_rad(irad)-one)*uin(l,i,j,k,8+irad)
                 erad(l) = erad(l)+uin(l,i,j,k,8+irad)
              end do
           enddo
           do irad = 1,ngrp
              do l = 1, ngrid
                 q(l,i,j,k,firstindex_er+irad) = uin(l,i,j,k,firstindex_er+irad)
                 erad(l) = erad(l)+uin(l,i,j,k,firstindex_er+irad)
              end do
           enddo
#endif

           
           ! Compute thermal pressure through EOS
           do l = 1, ngrid
              sum_dust = 0.0d0
#if NDUST>0              
              do idust = 1, ndust
                 sum_dust=sum_dust + uin(l,i,j,k,firstindex_ndust+idust)/max(uin(l,i,j,k,1),smallr)
              end do
#endif  
              etot = uin(l,i,j,k,5) - emag(l) -erad(l)
              eint = etot-eken(l)*q(l,i,j,k,1)
              if(energy_fix)eint=uin(l,i,j,k,nvar)
            
              call pressure_eos((1.0d0-sum_dust)*q(l,i,j,k,1),eint,pp_eos)
              if(dust_diffusion_test) q(l,i,j,k,5)=uin(l,i,j,k,5)*(gamma-1.0d0)
              q(l,i,j,k,5)= MAX(pp_eos,smallp)
           end do

        end do
     end do
  end do
#if NVAR>8+NENER
  do n = firstindex_pscal+1, firstindex_pscal+npscal
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/max(uin(l,i,j,k,1),smallr)
              end do
           end do
        end do
     end do
  end do
#endif

#if NDUST>0
#if NDUSTPSCAL>0
  do n_d=1,ndust
     do n=1,ndust_pscal
        do k = ku1, ku2
           do j = ju1, ju2
              do i = iu1, iu2
                 do l = 1, ngrid
                    q(l,i,j,k,firstindex_dustpscal+(n-1)*ndust+n_d) = uin(l,i,j,k,firstindex_dustpscal+(n-1)*ndust+n_d)/q(l,i,j,k,firstindex_ndust+n_d)
                 end do
              end do
           end do
     end do        
     end do
  end do
#endif
#endif

  
end subroutine ctoprimdust


!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_boundary_dust(ilevel)
  use amr_commons
  use hydro_commons
  use hydro_parameters
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor=1
  integer::i,idust,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::switch,dx,dx_loc,scale
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ndust,1:ndim),save::ff

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over boundaries
  do ibound=1,nboundary

     ! Compute direction of reference neighbors
     boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
     if(boundary_dir==1)inbor=2
     if(boundary_dir==2)inbor=1
     if(boundary_dir==3)inbor=4
     if(boundary_dir==4)inbor=3
     if(boundary_dir==5)inbor=6
     if(boundary_dir==6)inbor=5

     ! Compute index of reference cells
     ! Reflexive boundary
     if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     ! Free boundary
     if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! Imposed boundary (used only for flag1)
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)

     ! Vector sign switch for reflexive boundary conditions
     gs=(/1,1,1/)
     if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
     if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
     if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do

        ! Gather neighboring reference grid
        do i=1,ngrid
           ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather neighboring reference cell
           iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
           do i=1,ngrid
              ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
           end do

           ! Wall and free boundary conditions
           if((boundary_type(ibound)/10).ne.2)then

              ! Gather reference hydro variables
              do ivar=1,ndim
                 do idust=1,ndust
                 do i=1,ngrid
                    ff(i,idust,ivar)=v_dust(ind_cell_ref(i),idust,ivar)
                 end do
                 end do
              end do
              ! Scatter to boundary region
              do ivar=1,ndim
                 switch=gs(ivar)
                 do idust=1,ndust
                 do i=1,ngrid
                    v_dust(ind_cell(i),idust,ivar)=ff(i,idust,ivar)*switch
                 end do
                 end do
              end do

              ! Imposed boundary conditions
          
           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_force for level ',I2)

end subroutine make_boundary_dust


