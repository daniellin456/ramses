subroutine dust_diffusion_fine(ilevel)
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
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust,idim,nstepinski
  integer,dimension(1:nvector),save::ind_grid
 
  if(numbtot(1,ilevel)==0)return

  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call dustdifffine1(ind_grid,ngrid,ilevel)
  end do
111 format('   Entering dust_diffusion_fine for level ',i2)


end subroutine dust_diffusion_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dustdifffine1(ind_grid,ncache,ilevel)
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
  
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvardust),save::uloc  

  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust+1+ndust_pscal*ndust,1:ndim),save::flux_dust

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer ,dimension(1:nvector                   ),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer, dimension(1:ndust,1:ndim):: dv_ind1
  integer, dimension(1:ndust+1+ndust_pscal*ndust,1:ndim):: dv_ind2

  integer::idust,ht,idust_pscal
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
  real(dp), dimension(1:ndust) ::d_grain,l_grain
  real(dp):: epsilon_0
  real(dp),dimension(1:ndust):: dustMRN
  real(dp),dimension(1:ndim):: dvgtemp
  real(dp)::dust_temp,flxmax
  epsilon_0 = dust_ratio(1)  
  pi =3.14159265358979323846_dp
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0
  u1dust=0.0d0
  u2dust=0.0d0
  flux_dust = 0.0d0  
  uloc = 0.0d0

  ok   = .false.
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
  
 do idim=1,ndim
    do idust=1,ndust
       !Index for u1 u2
       dv_ind1(idust,idim)=ndim*(idust-1)+idim
       !Index for uloc
       dv_ind2(idust,idim)=ndust+1+ndust_pscal*ndust+ndim*(idust-1)+idim
    end do
    !Internal energy
    
    dv_ind2(ndust+ndust*ndust_pscal+1,idim)=ndust+1+ndust_pscal*ndust+ndim*ndust+idim
 end do
 
 !Associate the pscal to their dust velocity
 do idim=1,ndim
    do idust=1,ndust
       do idust_pscal=1,ndust_pscal
          !print *,(idust_pscal-1)*ndust+idust,dv_ind2(idust,idim)
          dv_ind2(ndust+(idust_pscal-1)*ndust+idust,idim)=dv_ind2(idust,idim)
       end do
    end do
 end do
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

     !Check if neighboring grid exists
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
     !If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do idim=1,ndim
              do idust=1,ndust
                 do i=1,nbuffer
                    u1dust(i,j,dv_ind1(idust,idim))=v_dust(ibuffer_father(i,j),idust,idim)
                 end do
              end do
           end do
           do ivar=1,nvar+3 
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
        call interpol_hydro_dust(u1dust,ind1,u2dust,nbuffer)

     end if
     !Loop over 2x2x2 cells
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
        !Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        !Gather dust variables

        do i=1,nexist
           do idust=1,ndust
              uloc(ind_exist(i),i3,j3,k3,idust)=uold(ind_cell(i),firstindex_ndust+idust)
           end do
           
           !Dust pscals
           do idust=1,ndust
              do idust_pscal=1,ndust_pscal
                 uloc(ind_exist(i),i3,j3,k3,ndust+(idust_pscal-1)*ndust+idust)=uold(ind_cell(i),firstindex_dustpscal+(idust_pscal-1)*ndust+idust)
              end do
           end do
           
           !We retrieve enint
           d= max(uold(ind_cell(i),1),smallr)
           A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
           e_mag=0.5d0*(A**2.0+B**2.0+C**2.0)
           u=0.0d0; v=0.0d0; w=0.0d0
           u = uold(ind_cell(i),2)/d
#if NDIM>1
           v = uold(ind_cell(i),3)/d
#endif
#if NDIM>2                 
           w = uold(ind_cell(i),4)/d
#endif
           e_kin=0.5d0*d*(u**2+v**2+w**2)
           
#if NENER>0
           do irad=1,nener
              e_mag=e_mag+uold(ind_cell(i),8+irad)
           end do
#endif
           if(energy_fix)then
              enint=uold(ind_cell(i),nvar)
           else
              enint=uold(ind_cell(i),5)-e_kin-e_mag
              
           end if
           !Internal Energy
           uloc(ind_exist(i),i3,j3,k3,ndust+1+ndust_pscal*ndust)=enint

           sum_dust =0.0d0
           do idust=1,ndust
              sum_dust = sum_dust + uold(ind_cell(i),firstindex_ndust+idust)/d 
           end do
 
          !Drift velocity 
           do idim= 1,ndim
              dvgtemp(idim)=0.0d0
                 do idust=1,ndust
                    uloc(ind_exist(i),i3,j3,k3,dv_ind2(idust,idim))= v_dust(ind_cell(i),idust,idim)
                    !dvgas
                    dvgtemp(idim)=dvgtemp(idim)- v_dust(ind_cell(i),idust,idim)*uold(ind_cell(i),firstindex_ndust+idust)/d/(1.0-sum_dust)

                 end  do
                 uloc(ind_exist(i),i3,j3,k3,dv_ind2(ndust+1+ndust_pscal*ndust,idim))=dvgtemp(idim)
              end do
           end do

           do i=1,nbuffer
              d= max(u2(i,ind_son,1),smallr)
              A=0.5d0*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1))
              B=0.5d0*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2))
              C=0.5d0*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3))
              e_mag=0.5d0*(A**2.0+B**2.0+C**2.0)
              
              u=0.0d0; v=0.0d0; w=0.0d0
              u = u2(i,ind_son,2)/d
#if NDIM>1
              v =u2(i,ind_son,3)/d
#endif
#if NDIM>2                 
              w = u2(i,ind_son,4)/d
#endif
                 
              e_kin=0.5d0*d*(u**2+v**2+w**2)

#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+u2(i,ind_son,8+irad)
              end do
#endif
              if(energy_fix)then
                 enint=u2(i,ind_son,nvar)
              else
                 enint=u2(i,ind_son,5)-e_kin- e_mag
              end if
              uloc(ind_nexist(i),i3,j3,k3,ndust+1+ndust_pscal*ndust)=enint
              do idust=1,ndust
                 uloc(ind_nexist(i),i3,j3,k3,idust)=u2(i,ind_son,firstindex_ndust+idust)
              end do
              !Dust pscals
              do idust=1,ndust
                 do idust_pscal=1,ndust_pscal
                    uloc(ind_nexist(i),i3,j3,k3,ndust+(idust_pscal-1)*ndust+idust)=u2(i,ind_son,firstindex_dustpscal+(idust_pscal-1)*ndust+idust)
                 end do
              end do
              sum_dust =0.0d0
              do idust =1 ,ndust
                 sum_dust= sum_dust + u2(i,ind_son,firstindex_ndust+idust)/d
              end do
             !Drift velocity
              do idim= 1,ndim
                 dvgtemp(idim)=0.0d0

                 do idust=1,ndust
                    uloc(ind_nexist(i),i3,j3,k3,dv_ind2(idust,idim))= u2dust(i,ind_son,dv_ind1(idust,idim))

                    dvgtemp(idim)=dvgtemp(idim)-u2(i,ind_son,firstindex_ndust+idust)/d/(1.0d0-sum_dust)*u2dust(i,ind_son,dv_ind1(idust,idim))

                 end do
                 uloc(ind_nexist(i),i3,j3,k3,dv_ind2(ndust+1+ndust_pscal*ndust,idim))= dvgtemp(idim)

              end do

           end do
 
        end do

     end do
  end do
end do
     !End loop over cells
end do
end do
!end do
!End loop over neighboring grids
  
  !-----------------------------------------------
  ! Compute flux_dust due to dust diffusion
  !-----------------------------------------------
   call dustdiff_predict(uloc,flux_dust,dx,dx,dx,dtnew(ilevel),ncache,dv_ind2) 

  !Reset flux_dustes at refined interfaces
  do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
      do i3=i3min,i3max+i0
         do idust=1,(ndust+1)+ndust_pscal*ndust
            do i=1,ncache
               if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                  flux_dust(i,i3,j3,k3,idust,idim)=0.0d0
               end if
            end do
         end do      
      end do
      end do
      end do
   end do
  !--------------------------------------------------------
  !Conservative  update at level ilevel for the dust flux_dustes
   !--------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1   
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
              if(update_eint.and..not.dust_diffusion_test) unew(ind_cell(i),5)=unew(ind_cell(i),5) +(flux_dust(i,i3,j3,k3,ndust+1+ndust_pscal*ndust,idim)&
                   &-flux_dust(i,i3+i0,j3+j0,k3+k0,ndust+1+ndust_pscal*ndust,idim))
              do idust=1,ndust+ndust_pscal*ndust
                 unew(ind_cell(i),firstindex_ndust+idust)=unew(ind_cell(i),firstindex_ndust+idust) +flux_dust(i,i3,j3,k3,idust,idim)&
                      &-flux_dust(i,i3+i0,j3+j0,k3+k0,idust,idim)
              enddo                     
           end if
     end do
  end do
  end do
  end do
end do

  if(ilevel>levelmin)then
  !-----------------------------------------------------
  ! update at level ilevel-1
  !-----------------------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     !----------------------
     ! Left flux_dust at boundary
      !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh)   = i
        end if
     end do
     ! Conservative update of the flux_dust
     do idust=1,ndust+ndust_pscal*ndust
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
              do i3=i3min,i3max-i0
                 do i=1,nb_noneigh
                    unew(ind_buffer(i),firstindex_ndust+idust)=unew(ind_buffer(i),firstindex_ndust+idust) &
                         &-flux_dust(ind_cell(i),i3,j3,k3,idust,idim)*oneontwotondim
           end do
        end do
     end do
  end do
end do
if(update_eint.and..not.dust_diffusion_test) then 
   ! Loop over boundary cells
   do k3=k3min,k3max-k0
      do j3=j3min,j3max-j0
         do i3=i3min,i3max-i0
            do i=1,nb_noneigh
               unew(ind_buffer(i),5)=unew(ind_buffer(i),firstindex_ndust+ndust+1+ndust_pscal*ndust) &
                         &-flux_dust(ind_cell(i),i3,j3,k3,idust,idim)*oneontwotondim
            end do
         end do
      end do
   end do
end if
!-----------------------
     ! Right flux_dust at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh)   = i
        end if
     end do
     ! Conservative update of the flux_dust
     do idust=1,ndust+ndust_pscal*ndust
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
             unew(ind_buffer(i),firstindex_ndust+idust)=unew(ind_buffer(i),firstindex_ndust+idust) &
                   &+flux_dust(ind_cell(i),i3+i0,j3+j0,k3+k0,idust,idim)*oneontwotondim
       
           end do
        end do
        end do
     end do
  end do
  if(update_eint.and..not.dust_diffusion_test) then 
     ! Loop over boundary cells
   do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
         do i3=i3min,i3max+i0
            do i=1,nb_noneigh
               unew(ind_buffer(i),5)=unew(ind_buffer(i),firstindex_ndust+ndust+1+ndust_pscal*ndust) &
                    &-flux_dust(ind_cell(i),i3+i0,j3+j0,k3+k0,idust,idim)*oneontwotondim
            end do
         end do
      end do
   end do
end if
end do
  ! End loop over dimensions
end if
 
end subroutine dustdifffine1

