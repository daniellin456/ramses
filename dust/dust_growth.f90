subroutine dust_growth_fine(ilevel)
  use amr_commons
  use hydro_commons
  use units_commons
  use poisson_commons
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
!  use smol2other
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif  
  integer :: info_growth
  integer::ilevel,nstp
  !--------------------------------------------------------------------------
  ! Enforcing relaxation damping to velocity components at early times,     !
  ! of in geometric regions. This now happens after set_uold, so do to uold !
  !--------------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::ix,iy,iz,id1,ig1,ih1,id2,ig2,ih2
  integer::ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp),dimension(ndim)::vel

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp):: epsilon_0,sum_dust,vrel,s_dust,vfrag_adim,x0,dx,pi,scale,dx_loc,v_brownian,d,u,v,w,enint,cs,ekin,emag,e,alpha_disk,rossby_disk,stokes_number,sc_number,vrel1,a,b,c,t_stop
  integer::idust,jdust,irad
  integer:: maxiter_growth,maxiter_growth_loc,maxiter_growth_all
  real(dp),dimension(1:ndust):: dustMRN,dust0,epsilondust
  real(dp), dimension(1:ndust)::d_grain,l_grain,m_grain
  real(dp),dimension(1:ndust+1)::l_grain_interface,m_grain_interface
  real(dp), dimension(1:ndust,1:ndust):: vdrift
  real(dp),dimension(1:ndust)::rhodust_loc
  integer::nstep_growth,isubcycle_growth

  integer,dimension(1:3,1:2,1:8)::iii,jjj
  epsilon_0 = dust_ratio(1)
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ! Mesh spacing at that level

  !Put to 1 if you want infos on growth
  info_growth=0
  maxiter_growth=0
  maxiter_growth_loc=0
  maxiter_growth_all=0

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  x0 = boxlen/2.0d0
  pi=ACOS(-1.0d0)

  dx=0.5d0**ilevel
  dx_loc=dx*scale

  !Sous-cyclage smol
  isubcycle_growth=0
  alpha_disk=1d-2
  rossby_disk=3.0d0
  vfrag_adim=v_frag/scale_l*scale_t

  
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
  

  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    end do

    ! Loop over cells
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
         ind_cell(i)=ind_grid(i)+iskip
      end do
          ! Gather cell centre positions
      do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
         end do
      end do
      ! Rescale position from code units to user units
      do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
         end do
      end do
#if NDUSTPSCAL>0      
      if (stepinski) then      
         do i=1,ngrid
            d=max(uold(ind_cell(i),1),smallr)
            u=uold(ind_cell(i),2)/d
            v=uold(ind_cell(i),3)/d
            w=uold(ind_cell(i),4)/d
            A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
            B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
            C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
            e=0.5d0*d*(u**2+v**2+w**2)+0.5d0*(A**2+B**2+C**2) 
#if NENER>0
            do irad=1,nener
               e=e+uold(ind_cell(i),8+irad)
            end do
#endif
            if(energy_fix)then
               enint=uold(ind_cell(i),nvar)
            else
               enint=uold(ind_cell(i),5)-e
            end if
            sum_dust=0.0d0
            do idust=1,ndust
               sum_dust=sum_dust+uold(ind_cell(i),firstindex_ndust+idust)/d
            end do
            call soundspeed_eos(d*(1.0d0-sum_dust),enint,cs)

            do idust=1,ndust
               vrel1= sqrt(v_dust(ind_cell(i),idust,1)**2.0+v_dust(ind_cell(i),idust,2)**2.0+v_dust(ind_cell(i),idust,3)**2.0)
               s_dust=uold(ind_cell(i),firstindex_dustpscal+idust)/uold(ind_cell(i),firstindex_ndust+idust)
               t_stop =  d_grain(idust)*s_dust*SQRT(pi*gamma/8.0d0)/cs/d

               Stokes_number=(t_stop/sqrt(3.0*pi/32.0d0/d))
               sc_number=(1.0-Stokes_number)*sqrt(1.0+vrel1**2./(sqrt(sqrt(2.0)*Rossby_disk*alpha_disk)*cs)**2.)
               vrel=sqrt(2.0**(3./2.)*rossby_disk*alpha_disk)*sqrt(1.0-sc_number)/sc_number*cs
               if(vrel.le.vfrag_adim) then
               !print *, s_dust,dtnew(ilevel),dtnew(ilevel)*uold(ind_cell(i),firstindex_dustpscal+idust)/d_grain(idust)*vrel
                  s_dust=s_dust+dtnew(ilevel)*uold(ind_cell(i),firstindex_dustpscal+idust)/d_grain(idust)*vrel
               else
                  s_dust=max(s_dust-dtnew(ilevel)*uold(ind_cell(i),firstindex_dustpscal+idust)/d_grain(idust)*vrel,l_grain(1)) ! We impose a minimum grain size 
               endif
               uold(ind_cell(i),firstindex_dustpscal+idust)=s_dust*uold(ind_cell(i),firstindex_ndust+idust)
            end do
         
      enddo
   end if
#endif
   
   enddo
enddo
111 format('   dust_growth_fine for level ',I2)

end subroutine dust_growth_fine
