!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar, idust, i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xn,x0,sum_dust, t_stop, rho_0, epsilon_0, delta_rho0,v0, pi,P_0,rho_gas, dusttogas, phidust
      pi =3.14159265358979323846_dp
      dusttogas =  1.0d0!0.10
      rho_gas =1.0_dp
      rho_0 = rho_gas + dusttogas*rho_gas
      epsilon_0 =  dusttogas*rho_gas/rho_0
      phidust=0.8
      delta_rho0= 0.0001
      v0  = delta_rho0
      q(1:nn,2)=0.0d0
      q(1:nn,3)=0.0d0
      q(1:nn,4)=0.0d0
      q(1:nn,5)=0.0d0
      q(1:nn,6)=0.0d0
      q(1:nn,7)=0.0d0
      q(1:nn,8)=0.0d0
      q(1:nn,nvar+1)=0.0d0
      q(1:nn,nvar+2)=0.0d0
      q(1:nn,nvar+3)=0.0d0
  
      do ivar=9,nvar
         q(1:nn,ivar)=0.0d0
      end do
   
      do i=1, nn
         ! Compute position in normalized coordinates
         xn=abs(x(i,1))
        
         sum_dust =0.0d0
         q(i,1)=rho_0
         q(i,2)=v0*sin(2.0_dp*pi*xn)
         q(i,firstindex_ndust+1) = epsilon_0*phidust
         q(i,firstindex_ndust+2) = epsilon_0*(1.0d0-phidust)
         q(i,5)=(1.0d0-epsilon_0)*q(i,1)
      end do
      
     ! Convert primitive to conservative variables
     ! density -> density
     u(1:nn,1)=q(1:nn,1)
     ! velocity -> momentum
     u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
     u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
     u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
     ! kinetic energy
     u(1:nn,5)=0.0d0
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
     ! pressure -> total fluid energy
     u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
     ! magnetic energy -> total fluid energy
     u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
     u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
     u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
     u(1:nn,6:8)=q(1:nn,6:8)
     u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
     ! non-thermal pressure -> non-thermal energy
     ! non-thermal energy   -> total fluid energy
     do ivar=1,nener-ngrp
        u(1:nn,8+ivar)=q(1:nn,8+ivar)/(gamma_rad(ivar)-1.0d0)
        u(1:nn,5)=u(1:nn,5)+u(1:nn,8+ivar)
     enddo
    ! Radiative transfer
#if NGRP>0
     ! radiative energy   -> total fluid energy
     do ivar=1,ngrp
        u(1:nn,firstindex_er+ivar)= q(1:nn,firstindex_er+ivar)
        u(1:nn,5)=u(1:nn,5)+ u(1:nn,firstindex_er+ivar)
     enddo
#if USE_M_1==1
     ! radiative flux
     do ivar=1,ndim*ngrp
        do i=1,ncache
           u(1:nn,fisrtindex_fr+ivar)=q(1:nn,firstindex+ivar)
        end do
        write(ilun)xdp
     end do
#endif
#endif
#endif
#if NEXTINCT>0
     ! Extinction
     if(extinction)u(1:nn,firstindex_extinct+nextinct)=1.0D0
#endif
#if NPSCAL>0
     ! passive scalars
     do ivar=1,npscal
        u(1:nn,firstindex_pscal+ivar)=q(1:nn,1)*q(1:nn,firstindex_pscal+ivar)
     end do
     ! Internal energy
     u(1:nn,nvar)=q(1:nn,5)/(gamma-1.0d0)
#endif
#if NDUST>0
     ! dust
     do ivar=1,ndust
        u(1:nn,firstindex_ndust+ivar)=q(1:nn,1)*q(1:nn,firstindex_ndust+ivar)
     end do
#endif
 
end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

!!$     xx=x(i,1)
!!$#if NDIM > 1
!!$     yy=x(i,2)
!!$#endif
!!$#if NDIM > 2
!!$     zz=x(i,3)
!!$#endif
!!$     ! ABC
!!$     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
!!$     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
!!$     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega
     
!!$     v(i,1)=vx
!!$#if NDIM > 1
!!$     v(i,2)=vy
!!$#endif
!!$#if NDIM > 2
!!$     v(i,3)=vz
!!$#endif

  end do


end subroutine velana
