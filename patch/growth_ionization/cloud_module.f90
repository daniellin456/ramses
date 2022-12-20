module cloud_module
  use amr_parameters
  use hydro_parameters,only:Msun
  
  real(dp)::bl_fac=1.   !multiply calculated boxlen by this factor

  !Initial conditions parameters for the dense core
  logical ::bb_test=.false. ! Activate Boss & Bodenheimer inital conditions instead of 1/R^2 density profile
  logical ::uniform_bmag=.false. ! Activate uniform magnetic field initial conditions for BE-like initial density profile
  real(dp)::mass_c=1.         !cloud mass in solar mass
  real(dp)::contrast=100.d0   !density contrast (used when bb_test=.true.)
  real(dp)::cont=1.           !density contrast (used when bb_test=.false.)
  real(dp)::rap=1.            !axis ratio
  real(dp)::ff_sct=1.         !freefall time / sound crossing time
  real(dp)::ff_rt=1.          !freefall time / rotation time
  real(dp)::ff_act=1.         !freefall time / Alfven crossing time
  real(dp)::ff_vct=1.         !freefall time / Vrms crossing time
  real(dp)::theta_mag=0.      !angle between magnetic field and rotation axis

  real(dp):: C2_vis=0.0d0 !Von Neumann & Richtmeyer artificial viscosity coefficient 3 en principe
  real(dp):: alpha_dense_core=0.5d0
  real(dp):: beta_dense_core=0.0d0
  real(dp):: crit=0.0d0
  real(dp):: delta_rho=0.0d0
  real(dp):: Mach=0.0d0
  real(dp):: r0_box=4.0d0

  ! PMS evolution related stuff
  logical :: rt_feedback=.false.       ! take into account RT feedback
 ! logical :: rt_protostar_m1=.false.   ! take into account RT feedback with M1
  logical :: rt_protostar_fld=.false.  ! take into account RT feedback with FLD
  logical :: PMS_evol=.false.          ! Take into account PMS evolution subgrid model
  logical :: Hosokawa_track=.false.    ! Take into account PMS evolution subgrid model
  real(dp):: dt_lsink_update=50        ! frequency of the sink luminosity update with PMS evolution (in yr)
  real(dp):: epsilonlib=0.0            ! Fraction of energy absorbed by the prostostars at the accretion shock
  real(dp):: mprotostar=0.0009546*Msun ! initial mass of the protostar (1 Mjup)
  real(dp):: rstar_init=2.5            ! Initial radius of the protostar in Rsun
  integer :: modell=0
  integer :: modrestart=0              ! name of model you want to restart from, this is an input
  real(dp):: facc_star_lum=0.75d0      ! fraction of the accretion luminosity radiated by the sinks
  real(dp):: facc_star=0.5d0           ! fraction of the sink accreted mass actually accreted by the star
  real(dp):: facc_star_mom=1.0d0       ! fraction of the angular momentum accreted by the sinks
  integer :: lum_injection=0           ! sink luminosity injection distribution (0=uniform, 1=peaked)
  integer::nmdot_PMS,nm_PMS,ndata_PMS
  integer ,allocatable,dimension(:)::nb_ligne_PMS
  real(dp),allocatable,dimension(:,:,:)::data_PMS

end module cloud_module

subroutine read_cloud_params(nml_ok)

  use amr_parameters
  use feedback_module
  use clfind_commons
  use cloud_module

  implicit none
  logical::nml_ok
  real(dp)::cellsize
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),parameter::pcincm=3.086d18

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/cloud_params/bl_fac
!  namelist/cloud_params/alpha_dense_core,beta_dense_core,crit,delta_rho &
!       & ,mass_c,rap,cont,ff_sct,ff_rt,ff_act,ff_vct,theta_mag,bb_test &
!       & ,contrast,Mach,uniform_bmag

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist

  ! Get some units out there
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)



end subroutine read_cloud_params

!#########################################################
!#########################################################
!#########################################################
subroutine boundary_frig(ilevel)
  Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
!  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine set up open boundary conditions which deals properly with div B 
  ! it uses the 2 last cells of the domain
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,j
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::pi,time
  integer ::ivar,jgrid,ind_cell_vois
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,Cwnm
  real(dp)::dx_min, fact, Emag,Emag0

! STG HACK - ignore if not MHD
! TODO: Take boundary cleaner and use for non-MHD solver
#ifndef SOLVERmhd
  return
#endif 



  !if you want to activate zero gradient boundary conditions with proper div B 
  ! remove the return
  return



  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  Cwnm = sqrt(8000./scale_T2)

  pi=ACOS(-1.0d0)

  time = t * Cwnm / boxlen



  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  dx_min = (0.5D0**levelmin)/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
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
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        

       do i=1,ngrid



        ! STG HACK CHECK FOR BORKED
        if (uold(ind_cell(i),1) .lt. 0d0) then
           write(*,*) "DENSITY < 0 BEFORE VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        if (uold(ind_cell(i),5) .lt. 0d0) then
           write(*,*) "TOTAL CELL ENERGY < 0 BEFORE VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        do j=5,8
           if (isnan(uold(ind_cell(i),j))) then
              write(*,*) "VARIABLE IS NAN BEFORE VELOCITY_FINE, OH NO", ind_cell(i),j,uold(ind_cell(i),1)
              call clean_stop
           end if
        end do


        !impose vanishing gradient conditions at the x  faces
        if(  xx(i,1) .lt. 2.*dx_min ) then 

             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),2))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! Subtract non-thermal pressure terms
!#if NENER>0
!           do irad=1,nener
!              uold(ind_cell(i),5) = uold(ind_cell(i),5) - uold(ind_cell(i),8+irad)
!           end do
!#endif

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+1) = uold(ind_cell_vois,6)
 

              uold(ind_cell(i),6)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),7) - uold(ind_cell(i),8) 
           else
              !should be equal to uold(ind_cell(i),7) of the preceeding case 
              uold(ind_cell(i),nvar+1) =  uold(ind_cell_vois,6) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),7)  - uold(ind_cell(i),8) 

              !ensure div B
              uold(ind_cell(i),6) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),7) - uold(ind_cell(i),8) 
           endif




           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

           ! Add back the non-thermal pressure
!#if NENER>0
!           do irad=1,nener
!              uold(ind_cell(i),5) = uold(ind_cell(i),5) + uold(ind_cell(i),8+irad)
!           end do
!#endif

        endif



        !impose vanishing gradient conditions at the x  faces
        if(  xx(i,1) .gt. boxlen-2.*dx_min ) then 

             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),1))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 
             ind_cell_vois = ind_cell_vois + ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 
              uold(ind_cell(i),6) = uold(ind_cell_vois,nvar+1)
 
              uold(ind_cell(i),nvar+1) = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 
           else
              !should be equal to uold(ind_cell(i),9) of the preceeding case 
              uold(ind_cell(i),6) =  uold(ind_cell(i),7) + uold(ind_cell(i),8)  + uold(ind_cell_vois,nvar+1) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 

              !ensure div B
              uold(ind_cell(i),nvar+1) = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 
           endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif





        !impose vanishing gradient conditions at the y  faces
        if(  xx(i,2) .lt. 2.*dx_min ) then 


             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),4))


           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 3 (1) is 4 (2)
           if(ind .eq. 3 .or. ind .eq. 4 .or. ind .eq. 7 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - 2*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 3 .or. ind .eq. 4 .or. ind .eq. 7 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+2) = uold(ind_cell_vois,7)
 
              uold(ind_cell(i),7)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6) - uold(ind_cell(i),8) 
           else
              !should be equal to uold(ind_cell(i),7) of the preceeding case 
              uold(ind_cell(i),nvar+2) =  uold(ind_cell(i),nvar+1 ) + uold(ind_cell_vois,7) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6)  - uold(ind_cell(i),8) 

              !ensure div B
              uold(ind_cell(i),7) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),6) - uold(ind_cell(i),8) 
           endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif

        if(  xx(i,2) .gt. boxlen-2.*dx_min ) then 
             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),3))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 3 (4) is 1 (2)
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 5 .or. ind .eq. 6) then 
             ind_cell_vois = ind_cell_vois + 2*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 5 .or. ind .eq. 6) then 
              uold(ind_cell(i),7) = uold(ind_cell_vois,nvar+2)
 
              uold(ind_cell(i),nvar+2)  = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+3) 
           else
              !should be equal to uold(ind_cell(i),10) of the preceeding case 
              uold(ind_cell(i),7) =  uold(ind_cell(i),6 ) + uold(ind_cell_vois,nvar+2) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1)  - uold(ind_cell(i),nvar+3) 

              !ensure div B
              uold(ind_cell(i),nvar+2) =  uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8)  -uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+3) 
           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 


        endif



        if(  xx(i,3) .lt. 2.*dx_min ) then 


             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),6))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 5 (6) is 1 (2)
           if(ind .eq. 5 .or. ind .eq. 6 .or. ind .eq. 7 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - 4*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 5 .or. ind .eq. 6 .or. ind .eq. 7 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+3) = uold(ind_cell_vois,8)
 
              uold(ind_cell(i),8)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6) - uold(ind_cell(i),7) 
           else
              !should be equal to uold(ind_cell(i),8) of the preceeding case 
              uold(ind_cell(i),nvar+3) =  uold(ind_cell(i), nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell_vois,8) - uold(ind_cell(i),6)  - uold(ind_cell(i),7) 

              !ensure div B
              uold(ind_cell(i),8) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),6) - uold(ind_cell(i),7) 

           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))


           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 


        endif


        if(  xx(i,3) .gt. boxlen-2.*dx_min ) then 


             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),5))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 1 (2) is 5 (6)
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 3 .or. ind .eq. 4) then 
             ind_cell_vois = ind_cell_vois + 4*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 3 .or. ind .eq. 4) then 
              uold(ind_cell(i),8) = uold(ind_cell_vois,nvar+3)
 
              uold(ind_cell(i),nvar+3)  = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+2) 
           else
              !should be equal to uold(ind_cell(i),nvar+3) of the preceeding case 
              uold(ind_cell(i),8) =  uold(ind_cell(i), 6) + uold(ind_cell(i),7) + uold(ind_cell_vois,nvar+3) - uold(ind_cell(i),nvar+1)  - uold(ind_cell(i),nvar+2) 

              !ensure div B
              uold(ind_cell(i),nvar+3) =  uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8)  -uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+2) 

           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif

        uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

        ! STG HACK CHECK FOR BORKED
        if (uold(ind_cell(i),5) .lt. 0d0) then
           write(*,*) "TOTAL ENERGY < 0 AFTER VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        do j=5,8
           if (isnan(uold(ind_cell(i),j))) then
              write(*,*) "VARIABLE IS NAN AFTER VELOCITY_FINE, OH NO", ind_cell(i),j
              call clean_stop
           end if
        end do


       enddo



       
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine boundary_frig
!#########################################################
!#########################################################
!#########################################################
!#########################################################
#if NIMHD==1
! We compute,eta_a,eta_o eta_h and psi and dump them in uold
subroutine nmhd_resist(ilevel)
  Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
  use cooling_module,only: mh,kb
  use radiation_parameters,only:mu_gas
  use hydro_commons
  use hydro_parameters,only:psi0
  use units_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine set up the magnetic resistivities in a self
  ! consistent way. It was put here to keep the makefile as
  ! it is.. this should be moved if this ever gets out of a
  ! patch.
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,j,ht
  integer::info,ibound,nx_loc,idim,neul=5,irad
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C,B_adim,dtlim
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all,ekin,emag,enint
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::pi,time
  integer ::ivar,jgrid,ind_cell_vois
  real(dp)::dx_min, fact
  real(dp):: etaobricolo,etaabricolo,gammaad_loc,eta_h

  real(dp) :: Temp_gas
  real(dp) :: dust_to_gas_ratio,sum_dust
  integer  :: idust
  real(dp) :: psi,chi
#if THERMAL_IONIZ>0
  real(dp) :: eps,ni,ns
#endif
  real(dp) :: nh,zeta_adim
  real(dp) :: d_grain,dustMRN
#if NDUST>0
  real(dp),dimension(1:ndust) :: n_dust
  real(dp), dimension(1:ndust) :: l_grain,l_grain_loc,mgrain
#endif
  integer :: niter_ionis_max,niter_ionis

#if NEXTINCT>0
    real(dp) :: ion_rate_fit
#endif
  dust_to_gas_ratio=0.006
  if(verbose)write(*,111)ilevel
  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale = dble(nx_loc)/boxlen
  dx_loc=dx/scale
  dx_min = (0.5D0**levelmin)/scale
  pi = acos(-1.0d0)
#if NDUST>0
  call size_dust(l_grain,ndust)
  d_grain=grain_dens(1)/scale_d
  l_grain=l_grain/scale_l  ! Adimensionning
  if(icy_grains)mgrain=4.0d0/3.d0*pi*d_grain*(l_grain-ice_mantle/scale_l)**3.
  if(.not.icy_grains)mgrain=4.0d0/3.d0*pi*d_grain*l_grain**3.
#endif

  !! d_grain == densite des grains adim
  !! l_grain == taille adim
  !! m_grain == masse adimensionne de 1 grain !
  zeta_adim=zeta_ionis*scale_t
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  dtlim=dtnew(ilevel)

  !We now resolve eq 31. from Marchand et al., 2021 : 1-psi=theta eps exp psi
  !!!!!!!!!!!!!!!!!!!!!Iterations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
  !print*, psi0
  !!!!!!!!!!!!!!!!!!!!!!End of iterations!!!!!!!!!!!!!!!!!!!!!!!
  ncache=active(ilevel)%ngrid
  niter_ionis_max=0
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid

           !Gather variables
           
           d=max(uold(ind_cell(i),1),smallr)
           nH=max(d,smallr)/((mu_gas*mH)/scale_m)!,1d5*scale_l**3)
           
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           
           ekin=0.5d0*d*(u**2.+v**2.+w**2.)
           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
           emag=0.5d0*(A**2+B**2+C**2)
#if NENER>0
           do irad=1,nener
              emag=emag+uold(ind_cell(i),8+irad)
           end do
#endif
           if(energy_fix) then
              enint=uold(ind_cell(i),nvar)
           else
              enint=uold(ind_cell(i),5)-emag-ekin
           endif
           sum_dust=0.0d0
#if NDUST>0
           do idust=1,ndust
              sum_dust=sum_dust+uold(ind_cell(i),firstindex_ndust+idust)/d
           end do
#endif           
           call temperature_eos(d*(1.0d0-sum_dust),enint,temp_gas,ht)

#if NEXTINCT>0
           zeta_adim=ion_rate_fit(uold(ind_cell(i),firstindex_extinct+1),lmh_rate)*scale_t
           if(.not.extinction) zeta_adim=default_ionisrate*scale_t
#else
           zeta_adim=default_ionisrate*scale_t
#endif

           !Dust abundances
#if NDUST>0
           do idust=1,ndust
              n_dust(idust) = uold(ind_cell(i),firstindex_ndust+idust)/(mgrain(idust)) !We cheat so under nH=1d5
           end do
           
#endif           
#if ICOAG>0
           chi = uold(ind_cell(i),firstindex_pscal+icoag)/uold(ind_cell(i),1)*scale_t
#endif

           !Get B in gauss
           B_adim= sqrt((A**2+B**2+C**2))

#if THERMAL_IONIZ >0
           psi=uold(ind_cell(i),firstindex_pscal+ipsi_therm)
           eps=uold(ind_cell(i),firstindex_pscal+ieps_therm)
           ni=uold(ind_cell(i),firstindex_pscal+ini_therm)
           ns=uold(ind_cell(i),firstindex_pscal+ins_therm)
#else
           psi=uold(ind_cell(i),firstindex_pscal+ipsi)
#endif



! Call ionization algo
           call calc_ioniz(nh,temp_gas,zeta_adim,B_adim,psi,&
#if THERMAL_IONIZ >0
                           &eps,ni,ns,&
#endif
#if ICOAG>0
                           &chi,&
#endif
#if NDUST>0
                           &n_dust,&
#endif
                           &etaobricolo,eta_h,etaabricolo,dx,dtlim,niter_ionis)

           niter_ionis_max=max(niter_ionis_max,niter_ionis)

           
#if THERMAL_IONIZ >0
           uold(ind_cell(i),firstindex_pscal+ipsi_therm)=psi
           uold(ind_cell(i),firstindex_pscal+ieps_therm)=eps
           uold(ind_cell(i),firstindex_pscal+ini_therm)=ni
           uold(ind_cell(i),firstindex_pscal+ins_therm)=ns
#else
           uold(ind_cell(i),firstindex_pscal+ipsi)=psi
#endif
              
           uold(ind_cell(i),firstindex_pscal+i_ohm)  = etaobricolo
           uold(ind_cell(i),firstindex_pscal+i_ad)   = etaabricolo !/!\ Ramses uses betaad and not etaad, this is eta_a/bsquare
           uold(ind_cell(i),firstindex_pscal+i_hall) = eta_h

           !endif 
           !print *, 'toto', eta_o,eta_a
        enddo
     end do
     ! End loop over cells
     
   end do
  ! End loop over grids
  if(myid.eq.1) print*," Max number of iteration for ionisation = ", niter_ionis_max
  111 format('   Entering nmhd_resist for level ',i2)







end subroutine nmhd_resist


subroutine calc_ioniz(nh,temp_gas,zeta_adim,B_adim,psiini,&
#if THERMAL_IONIZ >0
                     &eps,ni,ns,&
#endif
#if ICOAG>0
                     &chi,&
#endif
#if NDUST>0
                     &n_dust,&
#endif
                     &eta_o,eta_h,eta_a,dx,dtlim,niter_ionis)
  Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
  use cooling_module,only: mh,kb
  use radiation_parameters,only:mu_gas
  use hydro_commons
  use hydro_parameters,only:psi0,thetai
  use units_commons
  implicit none
  integer :: j
  real(dp) :: nh,temp_gas,B_adim,dx
  real(dp) :: eta_o,eta_h,eta_a,etaobricolo,etaabricolo
  real(dp) :: psiini,zeta_adim
  real(dp) :: B_gauss
#if THERMAL_IONIZ >0
  real(dp) :: eps,ni,ns
#endif
#if ICOAG>0
  real(dp),intent(in) :: chi
#endif
#if NDUST>0
  real(dp),intent(in),dimension(1:ndust) :: n_dust
#endif
#if THERMAL_IONIZ>0
  real(dp) :: sigmav_se,vs,n_s,eps0
  real(dp), dimension(4) :: fonct
  real(dp), dimension(4,4) :: jacob,invjacob
  real(dp) ::omegas_ions_s,sigmas_ions_s,sigmav_ions_s,t_sions_s,mu_s,vrms_s
#endif
#if NDUST>0
  real(dp), dimension(1:ndust) :: l_grain,mgrain,tau_k,alpha_k,J_tauk,z_k,n_k
  real(dp), dimension(1:ndust) :: t_sdust,vrms_dust,sigmav_dust,sigmas_dust,omegas_dust
  real(dp):: d_grain
#if THERMAL_IONIZ>0
  real(dp),dimension(1:ndust) :: dzkdeps,djkdeps
#endif
#else
#if ICOAG>0
  real(dp), dimension(1:nbincoag) :: l_grain,mgrain,tau_k,alpha_k,J_tauk,z_k,n_k
  real(dp), dimension(1:nbincoag) :: t_sdust,vrms_dust,sigmav_dust,sigmas_dust,omegas_dust
  integer  :: ichi
#if THERMAL_IONIZ>0
  real(dp),dimension(1:nbincoag) :: dzkdeps,djkdeps
#endif
#endif  
#endif
  real(dp) :: t_sions,t_sel,vrms_i,vrms_el,sigmav_ions,sigmav_el,sigmas_ions,omegas_ions,sigmas_el,omegas_el
  real(dp) :: as_He_dust,as_He_ions,as_He_el,mu_i,mu_e,clight
  real(dp) :: sigma_par,sigma_perp,sigma_H
  real(dp) :: psi_iter, psi_0,psi
  real(dp) :: fpsi,dfdpsi,epsone
  real(dp) :: convergence_ionis,eps_psi,n_i,cross_sectJ_k,n_e
  real(dp) :: depsilon_dpsi,dni_dpsi,dni_dpsi_1,dni_dpsi_2,dfdpsi0,dfdpsi1,dfdpsi2,dfdpsi3
  real(dp) :: eps_theta
  integer  :: niter_ionis,niter_ionis_max
  real(dp) :: xres,dtlim,dtt
  real(dp) :: pi
  real(dp) :: sigmav_ie,vi

  pi = acos(-1.0d0)

#if NDUST>0
  call size_dust(l_grain,ndust)
  d_grain=grain_dens(1)/scale_d
  l_grain=l_grain/scale_l  ! Adimensionning
  if(icy_grains)mgrain=4.0d0/3.d0*pi*d_grain*(l_grain-ice_mantle/scale_l)**3.
  if(.not.icy_grains)mgrain=4.0d0/3.d0*pi*d_grain*l_grain**3.

#endif

  clight=3d10
  as_He_dust=1.28
  as_He_ions=1.14
  as_He_el=1.16
  B_gauss= B_adim*sqrt((4.d0*pi*scale_d*(scale_v)**2))

#if NDUST>0
     n_k(:)=n_dust(:)
#endif

           !n_k adimensionne
#if ICOAG>0
           !idem remplir n_k
    if(chi<chi_min) then
      ichi=1
    else if(chi>chi_max) then
      ichi=ntotchi
    else
      ichi=floor(2+dlog(chi/chi_min)/dlog(chi_max/chi_min)*(ntotchi-1))
    end if
    do j=1,nbincoag
      l_grain(j)=table_coag(ichi,1+4*(j-1)+1)/scale_l
      mgrain(j) =table_coag(ichi,1+4*(j-1)+3)/scale_m
      n_k(j)    =table_coag(ichi,1+4*(j-1)+4)*nh
    end do
#endif           
    !Loop independant quantities
    !if(debug_ionis) temp_gas=10.0d0
    sigmav_ie=(2d-7)/sqrt(temp_gas/300.0d0)/(scale_v*scale_l**2.)
    vi=sqrt(8.0d0*kB*temp_gas/(pi*mu_ions*mH))/scale_v
    tau_k= (l_grain*scale_l)*kB*temp_gas/e_el_stat**2. ! l_grain in cgs to keep tau_k dimensionless
    alpha_k=sqrt(8.0d0/(pi*tau_k))


    !print *, n_k
    !Psi0 determined, serves as a guess for psi
    convergence_ionis=1d4
    niter_ionis=0
    !!!!!!!!!!!!!!!!!!!!!Iterations!!!!!!!!!!!!!!!!!!!!!!!

    psi=psiini


#if THERMAL_IONIZ>0



     sigmav_se=sigcoef_a_ioniz*(temp_gas/300d0)**(sigcoef_b_ioniz)
     vs=sqrt(8.0d0*kB*temp_gas/(pi*mu_ions_s*mH))!/scale_v
     eps0=1d0/thetai
     vi=vi*scale_v
     zeta_adim=zeta_adim/scale_t
     l_grain=l_grain*scale_l
     n_k=n_k/scale_l**3d0
     nh=nh/scale_l**3d0
     sigmav_ie=sigmav_ie*scale_v*scale_l**2d0
     eps_psi=eps
     n_i=ni
     n_s=ns


     ! WARNING: eps is actually eps-eps0 here !!! (but not for .not. thermal_ioniz)

     eps_theta=(eps_psi+eps0)*thetai 

     fonct(1)=ffpsi()
     fonct(2)=ffcharge()
     fonct(3)=ffreco()
     fonct(4)=ffrecos()


    do while(convergence_ionis>epsilon_ionis .and.niter_ionis<nitermax_ionis)

        !Jacobian matrix
       call calc_jacob()

        !Inverse of jacobian
       call inversejacob()

       psi=psi-sum(invjacob(1,1:4)*fonct(1:4))
       eps_psi=eps_psi-sum(invjacob(2,1:4)*fonct(1:4))
       n_i=n_i-sum(invjacob(3,1:4)*fonct(1:4))
       n_s=n_s-sum(invjacob(4,1:4)*fonct(1:4))
       eps_theta=(eps_psi+eps0)*thetai 

       fonct(1)=ffpsi()
       fonct(2)=ffcharge()
       fonct(3)=ffreco()
       fonct(4)=ffrecos()
       convergence_ionis=maxval(abs(fonct))
       niter_ionis=niter_ionis+1

    end do
    !We recompute 
    n_e=(eps_psi+eps0)*(n_i+n_s*qis)
    z_k(:)=psi*tau_k(:)+(-eps_psi**2.0*thetai**2.0-2d0*eps_psi*thetai)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0)
!#if NDUST>0
    !if(debug_ionis) then
      !! print*,firstindex_ndust,ndust,firstindex_ndust+ndust,firstindex_pscal
       !do idust=1,ndust
          !uold(ind_cell(i),firstindex_ndust+ndust+idust)=-z_k(idust)
       !end do
      !uold(ind_cell(i),firstindex_ndust+ndust+ndust+1)=n_i/scale_l**3.
      !uold(ind_cell(i),firstindex_ndust+ndust+ndust+2)=n_i*eps_psi/scale_l**3.
    !endif
!#endif 

    psiini=psi
    eps=eps_psi
    ni=n_i
    ns=n_s

    l_grain=l_grain/scale_l
    n_k=n_k*scale_l**3d0
    n_i=n_i*scale_l**3d0
    n_s=n_s*scale_l**3d0
    n_e=n_e*scale_l**3d0
    nh=nh*scale_l**3d0

#else
! Else if not thermal_ioniz


    


    do while(convergence_ionis>epsilon_ionis .and.niter_ionis<nitermax_ionis)
       
       eps_psi=(1.0d0-psi)/thetai*exp(-psi)
       eps_theta=eps_psi*thetai
       z_k(:)=psi*tau_k(:)+(1.0d0-eps_theta**2.0)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0)
       J_tauk(:)=(1.0d0-psi)+(2.0d0/tau_k(:))*(eps_theta**2.0+eps_theta)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0)
       n_i=-sum(z_k(:)*n_k(:))/(1.0d0-eps_psi)
       cross_sectJ_k=sum(n_k(:)*pi*l_grain(:)**2.0*J_tauk(:))

       !f(psi)
       fpsi=sigmav_ie*eps_psi*n_i**2.0+n_i*vi*cross_sectJ_k
       fpsi=fpsi/(zeta_adim*nh)-1d0

       !df/dpsi
       depsilon_dpsi=-eps_psi*(2.0d0-psi)/(1.0d0-psi)


       !dni_dpsi part : 1
       !dni_dpsi_1 = sum(n_k(:)*(-eps_theta**4.0 + (4.0d0*thetai**2.-thetai**3.*alpha_k(:))*eps_psi**2.+(2.0d0*thetai*alpha_K(:)-4.0d0*thetai**2.0)*eps_psi +1.0d0 - thetai*alpha_k(:))&
       dni_dpsi_1 = sum(n_k(:)*(-eps_theta**4.0 + (4.0d0-thetai*alpha_k(:))*eps_theta**2.+(2.0d0*alpha_K(:)-4.0d0*thetai)*eps_theta +1.0d0 - thetai*alpha_k(:))&
            &/(1.0d0 + eps_theta*alpha_k(:) + eps_theta**2.0d0)**2.0)
       !dni_dpsi part : 2
       dni_dpsi_2 = -(psi/(1.0d0-eps_psi)*depsilon_dpsi+1.0d0)*sum(n_k(:)*tau_k(:))/(1.0d0-eps_psi)

       !dni_dpsi
       dni_dpsi   = -1.0d0/(1.0d0-eps_psi)**2.*depsilon_dpsi*dni_dpsi_1+dni_dpsi_2

       dfdpsi0=dni_dpsi/n_i+depsilon_dpsi*(2.0d0*eps_theta+1.0d0)/(eps_psi*(1.0+eps_theta))
       dfdpsi1=sum(n_k(:)*pi*l_grain(:)**2.0/(tau_k(:)*(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0d0)))
       dfdpsi2=depsilon_dpsi*sum(n_k(:)*pi*l_grain(:)**2.0*(2.0d0*eps_psi*thetai**2.+thetai*alpha_k(:))/(tau_k(:)*(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0d0)**2.0d0))
       dfdpsi3=sum(n_k(:)*pi*l_grain(:)**2.0)


       dfdpsi=sigmav_ie*n_i*(n_i*depsilon_dpsi+2.0d0*eps_psi*dni_dpsi)
       dfdpsi=dfdpsi+2.0d0*n_i*vi*(eps_theta**2.+eps_theta)*(dfdpsi0*dfdpsi1+dfdpsi2)
       dfdpsi=dfdpsi+n_i*vi*(dni_dpsi/n_i*(1.0d0-psi)-1.0d0)*dfdpsi3
       dfdpsi=dfdpsi/(zeta_adim*nh)

       !Psi^n+1=Psi^n-f(psi)/dfdpsi(psi)
       psi=psi-fpsi/dfdpsi
       if(psi<psi0) then
          psi = psi0*0.9999999d0
          if(convergence_ionis==1d3) then
            convergence_ionis=epsilon_ionis*0.1  ! To not loop indefinitely on psi<psi0
          else
            convergence_ionis=1d3 !We cheated so we shouldn't converge at this step
          end if
       else if (psi>0.0d0) then
          psi=-1d-5
          convergence_ionis=1d4 !We cheated so we shouldn't converge at this step
       else
          convergence_ionis=abs(fpsi)
       end if
       niter_ionis=niter_ionis+1

    end do
    !!!!!!!!!!!!!!!!!!!!!!End of iterations!!!!!!!!!!!!!!!!!!!!!!!

    !We recompute 
    eps_psi=(1.0d0-psi)/thetai*exp(-psi)
    z_k(:)=psi*tau_k(:)+(1.0d0-eps_psi**2.0*thetai**2.0)/(1.0d0+eps_psi*thetai*alpha_k(:)+eps_psi**2*thetai**2.0)
    n_i=sum(-z_k(:)*n_k(:)/(1.0d0-eps_psi))
    n_e=eps_psi*n_i
!#if NDUST>0
    !if(debug_ionis) then
      !! print*,firstindex_ndust,ndust,firstindex_ndust+ndust,firstindex_pscal
       !do idust=1,ndust
          !uold(ind_cell(i),firstindex_ndust+ndust+idust)=-z_k(idust)
       !end do
      !uold(ind_cell(i),firstindex_ndust+ndust+ndust+1)=n_i/scale_l**3.
      !uold(ind_cell(i),firstindex_ndust+ndust+ndust+2)=n_i*eps_psi/scale_l**3.
    !endif
!#endif           
    psiini=psi

#endif
! End not. thermal ioniz






           !Resistivities
    sigmav_dust(:)=pi*(l_grain(:)*scale_l)**2.*sqrt(8.0*kB*temp_gas/(pi*2.0d0*mH))*(1.0d0+sqrt(pi/(2*tau_k(:))))

    t_sdust(:)=1.0d0/as_He_dust*((mgrain(:)*scale_m+2.0d0*mH)/(2.0d0*mH))/sigmav_dust(:)/(nH/scale_l**3.)
    
    mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
    mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)
#if THERMAL_IONIZ>0
    mu_s=2.0d0*mH*mu_ions_s*mH/(2.0d0*mH+mu_ions_s*mH)
#endif
           
    vrms_i=sqrt(8.0d0*kB*temp_gas/(pi*mu_i))*1d-5  ! These velocities need to be in km/s for the Pinto & Galli 2008 fit
    vrms_el=sqrt(8.0d0*kB*temp_gas/(pi*mu_e))*1d-5
#if THERMAL_IONIZ>0
    vrms_s=sqrt(8.0d0*kB*temp_gas/(pi*mu_s))*1d-5  ! These velocities need to be in km/s for the Pinto & Galli 2008 fit
#endif

    sigmav_el=3.16d-11*vrms_el**1.3
    sigmav_ions=2.4d-9*vrms_i**0.6
#if THERMAL_IONIZ>0
    sigmav_ions_s=2.4d-9*vrms_s**0.6
#endif

    t_sel=1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/(nH/scale_l**3.)
    t_sions=1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/(nH/scale_l**3.)
#if THERMAL_IONIZ>0
    t_sions_s=1.0d0/as_He_ions*((mu_ions_s*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions_s/(nH/scale_l**3.)
#endif

    sigmas_el=(n_e/scale_l**3.)*e_el_stat**2.*t_sel/m_el
    sigmas_ions=(n_i/scale_l**3.)*e_el_stat**2*t_sions/(mu_ions*mH)
    sigmas_dust=n_k/scale_l**3.*(z_k*e_el_stat)**2.*t_sdust/(mgrain*scale_m)
#if THERMAL_IONIZ>0
    sigmas_ions_s=(n_s/scale_l**3.)*e_el_stat**2*t_sions_s/(mu_ions_s*mH)
#endif

    omegas_el=-e_el_stat*B_gauss/clight/m_el
    omegas_ions=e_el_stat*B_gauss/clight/(mu_ions*mH)
    omegas_dust=z_k*e_el_stat*B_gauss/clight/(mgrain*scale_m)
#if THERMAL_IONIZ>0
    omegas_ions_s=e_el_stat*B_gauss/clight/(mu_ions_s*mH)
#endif

    sigma_par=sum(sigmas_dust(:))+sigmas_ions+sigmas_el
    sigma_perp=sigmas_ions/(1.0d0+(omegas_ions*t_sions)**2.0)+sigmas_el/(1.0d0+(omegas_el*t_sel)**2.0)+sum(sigmas_dust(:)/(1.0d0+(omegas_dust(:)*t_sdust(:))**2.0))
    sigma_H=-sigmas_ions*(omegas_ions*t_sions)/(1.0d0+(omegas_ions*t_sions)**2.0)-sigmas_el*(omegas_el*t_sel)/(1.0d0+(omegas_el*t_sel)**2.0)-sum(sigmas_dust(:)*(omegas_dust(:)*t_sdust(:))/(1.0d0+(omegas_dust(:)*t_sdust(:))**2.0))
#if THERMAL_IONIZ>0
    sigma_par=sigma_par+sigmas_ions_s
    sigma_perp=sigma_perp+sigmas_ions_s/(1.0d0+(omegas_ions_s*t_sions_s)**2.0)
    sigma_H=sigma_H-sigmas_ions_s*(omegas_ions_s*t_sions_s)/(1.0d0+(omegas_ions_s*t_sions_s)**2.0)
#endif

    eta_o=1.0d0/sigma_par*(clight**2./(4.*acos(-1.0d0)))/scale_l**2*scale_t
    eta_H=sigma_H/(sigma_perp**2.+sigma_H**2.)*(clight**2./(4.*acos(-1.0d0)))/scale_l**2*scale_t
    eta_a=(sigma_perp/(sigma_perp**2.+sigma_H**2.)-1.0d0/sigma_par)*(clight**2./(4.*acos(-1.0d0)))/scale_l**2*scale_t
    etaobricolo=eta_o

    if(nminitimestep.eq.1) then

       if(dtlim.ne.0.d0) then
          xres=eta_o
          if(xres.ne.0.d0) then
             dtt=coefohm*dx*dx/xres   !dtohm pour la cellule
          else
             dtt=1.d39
          endif
          if (dtt.le.dtlim) then
             etaobricolo=coefohm*dx*dx/(dtlim)
          endif
       endif
    endif
    
    etaabricolo=eta_a/B_adim**2d0

    if(nminitimestep.eq.1) then
       
       if(dtlim.ne.0.d0) then
          xres=eta_a
          if(xres.ne.0.d0) then
             dtt=coefad*dx*dx/xres   !dtAD pour la cellule
          else
             dtt=1.d39
          endif
          if (dtt.le.dtlim) then   ! on compare bien dtAD calcule pour la cellule (rhocelln) avec le temps de la simu
             etaabricolo=coefad*dx*dx/(dtlim*B_adim**2d0)
             !      write(*,*) 'la où ça seuille rho et B valent : ', rhocelln, bsquare
             !ici dtlim est le dt le plus petit : normal, ou seuillé si besoin est.
          else
             etaabricolo=eta_a/B_adim ! le betaad normal calcule avec rho a l'interface
          endif
       endif
       
    endif

    eta_o=etaobricolo
    eta_h=eta_h
    eta_a=etaabricolo

    !print*, nh/scale_l**3d0,psi,n_i,eps_psi,eta_o,eta_h,eta_a




#if THERMAL_IONIZ>0
contains

! Mean charge sum
double precision function nkzk()
  nkzk=sum(n_k(:)*(psi*tau_k(:)+(-eps_psi**2d0*thetai**2d0-2d0*eps_psi*thetai)/(1d0+eps_theta*alpha_k(:)+eps_theta**2d0)))
end function nkzk

! Mean Jk sum
double precision function Jk()
  Jk=sum(n_k(:)*pi*l_grain(:)**2d0*( (1d0-psi) + &
     & 2d0/tau_k(:)*(eps_theta**2d0 + eps_theta)/(eps_theta**2d0+alpha_k(:)*eps_theta+1d0)))
end function Jk

! Equation of psi
double precision function ffpsi()
  ffpsi=(1d0-psi)/(eps_theta*dexp(psi))-1d0
end function ffpsi

! Equation of charge neutrality
double precision function ffcharge()
  ffcharge= 1d0-(eps_psi+eps0)*(n_i+n_s*qis)/(n_i+n_s)+1d0/(n_i+n_s)*nkzk()
end function ffcharge

! Equation on recombination
double precision function ffreco()
  real(dp) :: adimreco
  adimreco=((zeta_adim+ksi_ioniz*n_s)*nh)
  ffreco=1d0-(nh*abund_ioniz+n_i)/nh-sigmav_ie/adimreco*(eps_psi+eps0)*n_i*(n_i+qis*n_s)&
        &-n_i*vi/adimreco*Jk()-kis_ioniz*n_i*(nh*abund_ioniz-n_s)/adimreco
end function ffreco

! Equation on species s recombination
double precision function ffrecos()
  real(dp) :: adimrecos
  if(abund_ioniz<1d-30 .or. n_s==0d0 .or. abs(n_s-abund_ioniz*nh)/(abund_ioniz*nh)<1d-6) then
    ffrecos=0d0
  else
    adimrecos=(zeta_adim+alpt_ioniz*temp_gas**0.5d0*dexp(-Temp_ioniz/temp_gas)*nh)*nh*abund_ioniz
    ffrecos=(1d0+kis_ioniz*n_i*nh*abund_ioniz/adimrecos)*(1d0-n_s/(abund_ioniz*nh))&  ! Have to multiply first term by ns0
           &-sigmav_se/adimrecos*(eps_psi+eps0)*n_s*(n_i+n_s*qis)&
           &-n_s*vs/adimrecos*Jk()&
           &-ksi_ioniz*n_s*(nh-nh*abund_ioniz-n_i)/adimrecos
  endif
end function ffrecos



subroutine calc_jacob()
  implicit none
  double precision :: adim4,adim3
    ! Jacobian matrix

    dzkdeps(:)=(-eps_theta**2d0*thetai*alpha_k(:)-4d0*eps_theta*thetai-alpha_k(:)*thetai)&
              &/(1d0+alpha_k(:)*eps_theta+eps_theta**2d0)**2d0
    djkdeps(:)=2d0/tau_k(:)*((alpha_k(:)-1d0)*eps_theta**2d0*thetai+thetai+2d0*eps_theta*thetai)&
              &/(1d0+alpha_k(:)*eps_theta+eps_theta**2d0)**2d0


    jacob(1,1)=(psi-2d0)/(eps_theta*dexp(psi))
    jacob(1,2)=-(1d0-psi)/((eps_psi+eps0)*eps_theta*dexp(psi))
    jacob(1,3)=0d0
    jacob(1,4)=0d0

    jacob(2,1)=1d0/(n_i+n_s)*sum(n_k(:)*tau_k(:))
    jacob(2,2)=-(n_i+qis*n_s)/(n_i+n_s)+1d0/(n_i+n_s)*sum(n_k(:)*dzkdeps(:))
    jacob(2,3)=(eps_psi+eps0)*n_s*(1d0-qis)/(n_i+n_s)**2d0-1d0/(n_i+n_s)**2d0*nkzk()
    jacob(2,4)=(eps_psi+eps0)*n_i*(qis-1d0)/(n_i+n_s)**2d0-1d0/(n_i+n_s)**2d0*nkzk()


    adim3=((zeta_adim+ksi_ioniz*n_s)*nh)
    jacob(3,1)=n_i/adim3*vi*sum(n_k(:)*pi*l_grain(:)**2d0)
    jacob(3,2)=-sigmav_ie/adim3*n_i*(n_i+qis*n_s)-n_i/adim3*vi*sum(n_k(:)*pi*l_grain(:)**2d0*djkdeps(:))
    jacob(3,3)=-1d0/nh-sigmav_ie*(eps_psi+eps0)*(2d0*n_i+qis*n_s)/adim3-vi/adim3*Jk()-kis_ioniz*(nh*abund_ioniz-n_s)/adim3
    jacob(3,4)=-sigmav_ie*(eps_psi+eps0)*n_i*(zeta_adim*qis-ksi_ioniz*n_i)/((zeta_adim+ksi_ioniz*n_s)**2d0*nh)+ksi_ioniz*n_i*vi/((zeta_adim+ksi_ioniz*n_s)**2d0*nh)*Jk()+kis_ioniz*n_i*(zeta_adim+ksi_ioniz*nh*abund_ioniz)/(nh*(zeta_adim+ksi_ioniz*n_s)**2d0)

    adim4=((zeta_adim+alpt_ioniz*temp_gas**0.5d0*dexp(-Temp_ioniz/temp_gas)*nh)*nh*abund_ioniz)
    jacob(4,1)=n_s*vs/adim4*sum(n_k(:)*pi*l_grain(:)**2d0)
    jacob(4,2)=-sigmav_se*n_s*(n_i+qis*n_s)/adim4-n_s*vs/adim4*sum(n_k(:)*pi*l_grain(:)**2d0*djkdeps(:))
    jacob(4,3)=kis_ioniz*abund_ioniz*nh/adim4*(1d0-n_s/(abund_ioniz*nh))-sigmav_se*(eps_psi+eps0)*n_s/adim4+ksi_ioniz*n_s/adim4  ! Have to multiply the first term by ns0 because of adim4
    jacob(4,4)=-1d0/(abund_ioniz*nh)-kis_ioniz*n_i/adim4-sigmav_se*(eps_psi+eps0)*(n_i+2d0*qis*n_s)/adim4-vs/adim4*Jk()-ksi_ioniz*(nh-abund_ioniz*nh-n_i)/adim4


end subroutine calc_jacob


subroutine inversejacob()
  implicit none
  Double precision :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv= &
      & jacob(1,1)*(jacob(2,2)*( jacob(3,3)*jacob(4,4)-jacob(3,4)*jacob(4,3))&
                  &+jacob(2,3)*(-jacob(3,2)*jacob(4,4)+jacob(3,4)*jacob(4,2))&
                  &+jacob(2,4)*( jacob(3,2)*jacob(4,3)-jacob(3,3)*jacob(4,2)))&
      &+jacob(1,2)*(jacob(2,1)*(-jacob(3,3)*jacob(4,4)+jacob(3,4)*jacob(4,3))&
                  &+jacob(2,3)*( jacob(3,1)*jacob(4,4)-jacob(3,4)*jacob(4,1))&
                  &+jacob(2,4)*(-jacob(3,1)*jacob(4,3)+jacob(3,3)*jacob(4,1)))&
      &+jacob(1,3)*(jacob(2,1)*( jacob(3,2)*jacob(4,4)-jacob(3,4)*jacob(4,2))&
                  &+jacob(2,2)*(-jacob(3,1)*jacob(4,4)+jacob(3,4)*jacob(4,1))&
                  &+jacob(2,4)*( jacob(3,1)*jacob(4,2)-jacob(3,2)*jacob(4,1)))&
      &+jacob(1,4)*(jacob(2,1)*(-jacob(3,2)*jacob(4,3)+jacob(3,3)*jacob(4,2))&
                  &+jacob(2,2)*( jacob(3,1)*jacob(4,3)-jacob(3,3)*jacob(4,1))&
                  &+jacob(2,3)*(-jacob(3,1)*jacob(4,2)+jacob(3,2)*jacob(4,1)))
     detinv=1d0/detinv

    ! Calculate the inverse of the matrix
    invjacob(1,1) = detinv*(jacob(2,2)*(jacob(3,3)*jacob(4,4)-jacob(3,4)*jacob(4,3))+jacob(2,3)*(jacob(3,4)*jacob(4,2)&
      &-jacob(3,2)*jacob(4,4))+jacob(2,4)*(jacob(3,2)*jacob(4,3)-jacob(3,3)*jacob(4,2)))
    invjacob(2,1) = detinv*(jacob(2,1)*(jacob(3,4)*jacob(4,3)-jacob(3,3)*jacob(4,4))+jacob(2,3)*(jacob(3,1)*jacob(4,4)&
      &-jacob(3,4)*jacob(4,1))+jacob(2,4)*(jacob(3,3)*jacob(4,1)-jacob(3,1)*jacob(4,3)))
    invjacob(3,1) = detinv*(jacob(2,1)*(jacob(3,2)*jacob(4,4)-jacob(3,4)*jacob(4,2))+jacob(2,2)*(jacob(3,4)*jacob(4,1)&
      &-jacob(3,1)*jacob(4,4))+jacob(2,4)*(jacob(3,1)*jacob(4,2)-jacob(3,2)*jacob(4,1)))
    invjacob(4,1) = detinv*(jacob(2,1)*(jacob(3,3)*jacob(4,2)-jacob(3,2)*jacob(4,3))+jacob(2,2)*(jacob(3,1)*jacob(4,3)&
      &-jacob(3,3)*jacob(4,1))+jacob(2,3)*(jacob(3,2)*jacob(4,1)-jacob(3,1)*jacob(4,2)))
    invjacob(1,2) = detinv*(jacob(1,2)*(jacob(3,4)*jacob(4,3)-jacob(3,3)*jacob(4,4))+jacob(1,3)*(jacob(3,2)*jacob(4,4)&
      &-jacob(3,4)*jacob(4,2))+jacob(1,4)*(jacob(3,3)*jacob(4,2)-jacob(3,2)*jacob(4,3)))
    invjacob(2,2) = detinv*(jacob(1,1)*(jacob(3,3)*jacob(4,4)-jacob(3,4)*jacob(4,3))+jacob(1,3)*(jacob(3,4)*jacob(4,1)&
      &-jacob(3,1)*jacob(4,4))+jacob(1,4)*(jacob(3,1)*jacob(4,3)-jacob(3,3)*jacob(4,1)))
    invjacob(3,2) = detinv*(jacob(1,1)*(jacob(3,4)*jacob(4,2)-jacob(3,2)*jacob(4,4))+jacob(1,2)*(jacob(3,1)*jacob(4,4)&
      &-jacob(3,4)*jacob(4,1))+jacob(1,4)*(jacob(3,2)*jacob(4,1)-jacob(3,1)*jacob(4,2)))
    invjacob(4,2) = detinv*(jacob(1,1)*(jacob(3,2)*jacob(4,3)-jacob(3,3)*jacob(4,2))+jacob(1,2)*(jacob(3,3)*jacob(4,1)&
      &-jacob(3,1)*jacob(4,3))+jacob(1,3)*(jacob(3,1)*jacob(4,2)-jacob(3,2)*jacob(4,1)))
    invjacob(1,3) = detinv*(jacob(1,2)*(jacob(2,3)*jacob(4,4)-jacob(2,4)*jacob(4,3))+jacob(1,3)*(jacob(2,4)*jacob(4,2)&
      &-jacob(2,2)*jacob(4,4))+jacob(1,4)*(jacob(2,2)*jacob(4,3)-jacob(2,3)*jacob(4,2)))
    invjacob(2,3) = detinv*(jacob(1,1)*(jacob(2,4)*jacob(4,3)-jacob(2,3)*jacob(4,4))+jacob(1,3)*(jacob(2,1)*jacob(4,4)&
      &-jacob(2,4)*jacob(4,1))+jacob(1,4)*(jacob(2,3)*jacob(4,1)-jacob(2,1)*jacob(4,3)))
    invjacob(3,3) = detinv*(jacob(1,1)*(jacob(2,2)*jacob(4,4)-jacob(2,4)*jacob(4,2))+jacob(1,2)*(jacob(2,4)*jacob(4,1)&
      &-jacob(2,1)*jacob(4,4))+jacob(1,4)*(jacob(2,1)*jacob(4,2)-jacob(2,2)*jacob(4,1)))
    invjacob(4,3) = detinv*(jacob(1,1)*(jacob(2,3)*jacob(4,2)-jacob(2,2)*jacob(4,3))+jacob(1,2)*(jacob(2,1)*jacob(4,3)&
      &-jacob(2,3)*jacob(4,1))+jacob(1,3)*(jacob(2,2)*jacob(4,1)-jacob(2,1)*jacob(4,2)))
    invjacob(1,4) = detinv*(jacob(1,2)*(jacob(2,4)*jacob(3,3)-jacob(2,3)*jacob(3,4))+jacob(1,3)*(jacob(2,2)*jacob(3,4)&
      &-jacob(2,4)*jacob(3,2))+jacob(1,4)*(jacob(2,3)*jacob(3,2)-jacob(2,2)*jacob(3,3)))
    invjacob(2,4) = detinv*(jacob(1,1)*(jacob(2,3)*jacob(3,4)-jacob(2,4)*jacob(3,3))+jacob(1,3)*(jacob(2,4)*jacob(3,1)&
      &-jacob(2,1)*jacob(3,4))+jacob(1,4)*(jacob(2,1)*jacob(3,3)-jacob(2,3)*jacob(3,1)))
    invjacob(3,4) = detinv*(jacob(1,1)*(jacob(2,4)*jacob(3,2)-jacob(2,2)*jacob(3,4))+jacob(1,2)*(jacob(2,1)*jacob(3,4)&
      &-jacob(2,4)*jacob(3,1))+jacob(1,4)*(jacob(2,2)*jacob(3,1)-jacob(2,1)*jacob(3,2)))
    invjacob(4,4) = detinv*(jacob(1,1)*(jacob(2,2)*jacob(3,3)-jacob(2,3)*jacob(3,2))+jacob(1,2)*(jacob(2,3)*jacob(3,1)&
      &-jacob(2,1)*jacob(3,3))+jacob(1,3)*(jacob(2,1)*jacob(3,2)-jacob(2,2)*jacob(3,1)))
end subroutine inversejacob

#endif








end subroutine calc_ioniz






#endif
!#########################################################
!#########################################################
!#########################################################
!#########################################################

! Added just to test the algo (same function as in the dust folder)
subroutine init_dust_ratio_loc(dustratio,epsilondust,ndustbins)
  use amr_commons
  use hydro_commons
  implicit none
  
  real(dp) :: dustratio,size_max_loc,size_min_loc
  real(dp), dimension(1:ndustbins):: epsilondust
  real(dp), dimension(1:ndustbins+1):: sdust
  real(dp) :: epsilon_0,Anorm
  integer  :: idust,ndustbins
  
  size_min_loc  = 5.0d-7
  size_max_loc  = 2.5d-5
  do idust =1,ndustbins+1
     sdust(idust) = 10**(log10(size_max_loc/size_min_loc)*dble(idust-1)/dble(ndustbins)+log10(size_min_loc))
  enddo
  Anorm=0.0d0
  do idust=1,ndustbins
     Anorm = Anorm + (sdust(idust+1)**(0.5)-sdust(idust)**(0.5))
  end do 
  do idust=1,ndustbins
     epsilondust(idust)= dustratio*(sdust(idust+1)**(0.5)-sdust(idust)**(0.5))/Anorm
  enddo
end subroutine init_dust_ratio_loc

subroutine size_dust_loc(sdust,ndustbins)
  use amr_commons
  use hydro_commons
  implicit none
  
  real(dp), dimension(1:ndustbins):: sdust
  real(dp), dimension(1:ndustbins+1):: sdust_interval
  integer  :: idust,ndustbins
  real(dp) :: size_max_loc,size_min_loc

  size_min_loc  = 5.0d-7
  size_max_loc  = 2.5d-5
  do idust =1,ndustbins+1
     sdust_interval(idust) = 10**(log10(size_max_loc/size_min_loc)*dble(idust-1)/dble(ndustbins)+log10(size_min_loc))
  enddo
  !We compute the average dust size in the bin to get the correct stopping time 
  do idust =1,ndustbins
     sdust(idust) = sqrt(sdust_interval(idust)*sdust_interval(idust+1))
     if(icy_grains) sdust(idust) = sdust(idust) + ice_mantle
  enddo
end subroutine size_dust_loc


