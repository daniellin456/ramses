!================================================================
!================================================================
!================================================================
! This routine reinforces initial conditions for RAMSES for some(x,y,z).
! Positions are in user units:
!   x(i,1:3) are in [0,boxlen]**ndim.
! U is the conservative variable vector. Conventions are here:
!   U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,>5): d.scalar
! Q is the primitive variable vector. Conventions are here:
!   Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,>5): scalar
! If nvar > 5, remaining variables (6:nvar) are treated as passive
! scalars in the hydro solver.
!   U(:,:) and Q(:,:) are in user units.
!
!
!  IMPORTANT: This DOES NOT reset the magnetic field or passive variables.
!   The former is because I don't know how to do that while preserving Div B = 0
!
!   The latter is because at this point we are not using any.
!
!
!================================================================
!================================================================
!================================================================
subroutine relax_condinit(xc,yc,zc,dx,rho,vel,pres)

  use amr_parameters
  use amr_commons
  use hydro_parameters
  use poisson_parameters
  implicit none

  real(dp)::xc,yc,zc,dx
  real(dp)::rho,vel(ndim),pres

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  integer::ivar,i,ii,jj,kk,i_rad
  real(dp)::x0,y0,z0,xc0,yc0,zc0
  real(dp)::xl,yl,zl,xr,yr,zr,R0,Rin,Rout,Rout2,Hout
  real(dp)::V0,P0,P1,pi,f,rr,abar,alpha,beta
  real(dp)::m,power,pm1,pp1,thetamax,GM,GM_code
  real(dp)::theta,g,D0,D1,dens,T0,r_cyl,v_cyl,r_sph,exp_arg,Rin0,Rout0,hR_clip
  real(dp)::pmin,R_integral,R_left,R_right,tmini,temp_zero,temp_minus1
  real(dp)::dR,dz,r_loop,z_loop,rho_l,rho_r,temp,temp_l,temp_r,&
             f_cent_out,f_grav_inR,f_grav_inZ,R200,R1,H1,c,Reff,p_tot,dens_tot
  real(dp)::Rmin_grav,Rmax_grav,dx_max,pres_l,pres_r,dR_scale,dZ_scale,dPdR,dPdz
  real(dp)::r_cyl_l,r_cyl_r,r_loop_l,r_loop_r,f_grav_inZ_l,f_grav_inZ_r,dPdz_R
  real(dp)::dPdz_R_gal,sech,sech_term, dPdR2_z0,dPdR2_z,dPdz2_R,pres_int,pres_intm1
  real(dp)::pres_p1,pres_l_p1,pres_r_p1,P_fact_cen0,P_fact_cen,soft
  real(dp)::v_turb(3),mom_tot(2),ztrans1,ztrans2,rtrans1,dztrans,dz_here,ftrans1,max_z, &
            maxz_fact,maxz_fact2,maxz_fact0,alpha2,aspect_ratio
  real(dp)::exp_arg_l,exp_arg_r,r_bound,eps,dPdR_z0,rom_int,f_table,fgrav_R
  real(dp)::k5,k4,k3,delta_trans,theta_trans,eta_trans,dx_ss
  real(dp)::rcyl_ratio,rsph_ratio,z1,rsph1,pres1,pres0,rho1,rho0,rho10,vcyl2,delta_pres
  real(dp),dimension(20)::func_params
  integer::nR_loop,iR,nZ_loop,iZ,iiz,iter,iter_tot,max_iter,i_table
  real(dp)::C0,C1,C2,C3,Q1,Q2,Q0,pres_zero,pres_minus1,xscl,r_sph_l,r_sph_r


  real(dp)::Mass_of_R,MPrime_of_R,Rho_of_R,RhoPrime_of_R,Pres_of_R,PPrime_of_R,Z_of_R,ZPrime_of_R
  external::Mass_of_R,MPrime_of_R,Rho_of_R,RhoPrime_of_R,Pres_of_R,PPrime_of_R,Z_of_R,ZPrime_of_R
  external :: dPdR_z0, dPdz_R, fgrav_R, dPdz_R_gal, sech, dPdR2_z0,dPdR2_z,dPdz2_R

  pi=ACOS(-1.0d0)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  dx_max = boxlen/2.d0**nlevelmax  ! max here means for maximum refinement level
  dR_scale = 1.d0
  dZ_scale = 40.d0
  eps = 1.d-9

! Gressel values ... simpler disk.  --> type ####1
! gamma = 7./5.
! R0 = 5 AU
! T0 = 108 K
! abar = 2.024 --> Chosen to Make T0 = Pres / ( rho R_idealGas abar ) = 108 K @ R0
! Sigma = 150 g / cm^2 @ R0
!  This gives D0 = Sigma / (2 * R0 * .0427878) ! I have double checked this!
!     .0427878 = integral of (exp(400(1/sqrt(1+x^2) - 1))) from 0 to .05.
!       --> D0 = 2.3434e-11 (only for unflared model, q = -1)
! theta_max = 0.4 radians above / below midplane. (only for unflared model, q = -1)
! R on = 0.5 - 5.5 AU
! res is 1024 x 384
!  v_theta = sqrt(G Msun / R ) * sqrt( (p+q)/400 + (1+q) - qR/r)
!    p = -3/2   q = -1 , gives: (note in this code I use alpha and beta instead of p and q)
!     v_theta = sqrt(G Msun / R) * r/R * sqrt( max(R/r - 1/240, 0) )
!  aspect ratio = H/R = 0.05.
!  beta_p0 = 2p/B^2 = 1e5 in midplane. I add this.
!   Bz = sqrt(beta_p0/2/p)
!

! Bai values for below: --> type 2
!   alpha = 2.0
!   D0 = 1
!   T0 = 1
!   GM = 1
!   m = 0.5
!   beta0 = 1e4 = q5/(B^2/8 pi) --> B0^2 + Bz^2 = ?

! Setting up an ideal MHD wedge version of Bai & Stone 2017
! Note, best way to get B field is use vector potential on
! corners and take first order finite difference. This well
! have left / right states be the same on neighbouring cells
!

! No B, Keplarian disk --> type 3

! Saving parameters to local variables to tidy up their use below
  Rin  =disk_innerR
  Rin0 =disk_innerR0
  Rout =disk_outerR
  Rout0=disk_outerR0
  Rout2=disk_outerR2 ! This is a boundary condition relevant R, see HTL2017, where beyond R>0.9 d and v are set.

  R0=disk_R0
  R1=disk_R1 ! scale radius
  H1=disk_H1 ! scale height

  Hout=Rout*H1/R1
  R200=disk_R200 ! virial radius, applicable for virial halo sims
  c=disk_conc ! halo concentration

  ! alpha = rho power law index ; beta = T power law index
  alpha = disk_rho_pow
  beta  = disk_cs_pow
  abar  = disk_abar
  power = (1.d0 - alpha)/2.d0 ! this is the flaring index.
  pm1 = power - 1.d0
  pp1 = power + 1.d0
  m = disk_m ! for Bai

  ! disk_xc etc in unis of dx. How many / what fraction of dx offset to I want to put the center. Usually 0 but I've looked at dx/2.
  x0 = boxlen/2.d0 + disk_xc*dx_max
  y0 = boxlen/2.d0 + disk_yc*dx_max
  z0 = boxlen/2.d0 + disk_zc*dx_max

  ! Values
  V0 = disk_V0
  D0 = disk_D0
  D1 = disk_D1
  T0 = disk_T0
  P0 = 8.3144598d7 * T0 * (D0/scale_d) / abar / (scale_l**2 / scale_t**2)
  P1 = 8.3144598d7 * T0 * (D1/scale_d) / abar / (scale_l**2 / scale_t**2)
  pmin = P1/1000.d3 ! arbitrarily setting lower pressure, factor 1000

  Tmini = pmin/(D1/scale_d) ! Here and below I call T = P/rho; it's really cs^2

  thetaMax=disk_thetaMax
  hR_clip = 4.6054 ! means at boundary you are 1./100. that of clipping point

  GM = disk_Mass ! Msun
  GM = GM*6.674d-8*1.98855d33 ! cgs
  GM_code = GM/scale_l**3*scale_t**2 ! GM in code units

  ! General Statement on 'disk_setup'
  !  The One's digit indicates the setup model--> 1 = Gressel, 2 = Bai
  !  Subsequent digits indicate a different method of doing the setup, for new methods I do a bitwise shift.
  !  I increment that digit when it's a similar method, but intrinsically different. To be cleaned up in future
  if ( disk_setup == 11 ) then ! Gressel, stepping cell by cell. Naive method that doesn't really work.
    stop "disk_setup = 11 not recommended. Try 25001 or 26001"
    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(14) = grav_angle
    func_params(15) = grav_width

    ! grav_params(6:7) regions without gravity
    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    nR_loop = 1+int(abs(r_cyl - R0)/dx*dR_scale)
    nz_loop = 1+int(abs(zc)/dx*dZ_scale)
    max_iter = max(nR_loop,nZ_loop)*16
    dR = sign(dx_max,r_cyl-R0)/dR_scale
    dz = sign(dx_max,zc)/dZ_scale

    pres = P0
    ! Integrate along R
    r_loop = R0
    iter = 0
    iter_tot = 0
    do while (r_loop .ne. r_cyl)
      r_loop = r_loop + dR
      if (dR .gt. 0.d0) then
        if (r_loop .gt. r_cyl) then
          dR = dR - (r_loop - r_cyl)
          r_loop = r_cyl
        endif
      else
        if (r_loop .lt. r_cyl) then
          dR = dR - (r_loop - r_cyl)
          r_loop = r_cyl
        endif
      endif

      if (r_loop .lt. Rin) then ! take value at Rin and have it linearly go to zero at Rin0
        v_cyl = max((r_loop - Rin0)/(Rin - Rin0),0.d0) * &
                    sqrt(GM/(max(Rin,gravity_params(2))*scale_l))/(scale_l/scale_t)
!          rho = D0 * (Rin/R0)**(alpha)*exp(-hR_clip*(1.d0 - r_loop/Rin))/scale_d
      else if (r_loop .lt. Rout) then ! V = Vcirc * sqrt ( (p+q)/400.d0 + (1+q) - qR/r )
        v_cyl = sqrt(GM/(max(r_loop,gravity_params(2))*scale_l)) * &
                  sqrt( max(0.d0,(alpha + beta)/400.d0 + 1))/(scale_l/scale_t) ! note, the first term assumes beta = -1
      else
        v_cyl = max((Rout0 - r_loop)/(Rout0 - Rout),0.d0) * &
                    sqrt(GM/(max(Rout,gravity_params(2))*scale_l)) * &
                    sqrt( max(0.d0,(alpha + beta)/400.d0 + 1))/(scale_l/scale_t)
!          rho = D0 * (Rout/R0)**(alpha)*exp(-hR_clip*(1.d0 - r_loop/Rout)/ &
!                       (1.d0 - boxlen/1.4142/Rout))/scale_d ! exp(-hR_clip) in corner
      endif
      rho = D0 * (r_loop/R0)**(alpha)/scale_d
      rho = max(rho,D1/scale_d)

      if ( (Rmin_grav < r_loop ) .and. (r_loop < Rmax_grav) ) then
        if (r_loop .gt. Rin) then
          f_grav_inR = GM_code/(max(r_loop,gravity_params(2)))**2 ! >0 means inwards
          f_cent_out = f_grav_inR
        else
          f_grav_inR = GM_code*r_loop/Rin**3 ! >0 means inwards
          f_cent_out = v_cyl**2/max(r_loop,gravity_params(2))
        endif
      else
        f_grav_inR = 0.d0
        v_cyl = 0.d0
        f_cent_out = 0.d0
      endif

      dPdR = rho*(f_cent_out - f_grav_inR)
      if (pres + dPdR*dR .lt. 0.5d0*pmin) then
        r_loop = r_loop - dR
        dR = 1.25*(pmin - Pres)/dPdR
        iter =iter + 1
!        if (iter .gt. max_iter) then
!          pres = pmin
!          r_loop = r_cyl
!        endif
      else
        if (pres + dPdR*dR .gt. 10.d0*pmin) then
          pres = pres + dPdR*dR
          dR = 2*dR
          if (abs(dR) .gt. dx_max/dR_scale) dR = sign(dx_max,dR)/dR_scale
        else
          pres = pres + dPdR*dR
        endif
      endif
      if (pres .lt. pmin) then
        pres = pmin
        dPdR = 0.d0
      endif

      iter = 0
      iter_tot = iter_tot + 1
!        if (iter_tot .gt. max_iter) then
!          pres = pmin
!          r_loop = r_cyl
!        endif
    enddo

    ! Now integrate along z
    z_loop = 0.d0
    iter = 0
    iter_tot = 0
    do while (pres .gt. pmin .and. z_loop .ne. zc)
      z_loop = z_loop + dz
      if (dz .gt. 0.d0) then
        if (z_loop .gt. zc) then
          dz = dz - (z_loop - zc)
          z_loop = zc
        endif
      else
        if (z_loop .lt. zc) then
          dz = dz - (z_loop - zc)
          z_loop = zc
        endif
      endif


!        if (r_cyl .lt. Rin) then ! take value at Rin and have it linearly go to zero at Rin0
!          r_loop = sqrt(Rin**2 + z_loop**2)
!          exp_arg = 400.d0*(1.d0 - Rin/r_loop)*(Rin/R0)**(-beta-1.d0)
!          rho = D0 * (Rin/R0)**(alpha)*exp(-exp_arg)*exp(-hR_clip*(1.d0 - r_cyl/Rin))/scale_d
!        else if (r_cyl .lt. Rout) then ! V = Vcirc * sqrt ( (p+q)/400.d0 + (1+q) - qR/r )
!        else
!          r_loop = sqrt(Rout**2 + z_loop**2)
!          exp_arg = 400.d0*(1.d0 - Rout/r_loop)*(Rout/R0)**(-beta-1.d0)
!          rho = D0 * (Rout/R0)**(alpha)*exp(-exp_arg)*exp(-hR_clip*(1.d0 - r_cyl/Rout)/ &
!                       (1.d0 - boxlen/1.4142/Rout))/scale_d ! exp(-hR_clip) in corner
!        endif
      r_loop = sqrt(r_cyl**2 + z_loop**2)
      exp_arg = 400.d0*(1.d0 - r_cyl/r_loop)*(r_cyl/R0)**(-beta-1.d0)
      rho = D0 * (r_cyl/R0)**(alpha)*exp(-exp_arg)/scale_d
      rho = max(rho,D1/scale_d)

      if ( (Rmin_grav < r_cyl ) .and. (r_cyl < Rmax_grav) ) then
        if (r_loop .gt. Rin) then
          f_grav_inZ = GM_code*z_loop/ &
                             (max(r_loop,gravity_params(2)))**3 ! >0 means inwards
        else
          f_grav_inZ = GM_code*z_loop/Rin**3 ! >0 means inwards
        endif
      else
        f_grav_inZ = 0.d0
      endif

      dPdz = -rho*f_grav_inZ
      if (pres + dPdz*dz .le. 0.5d0*pmin) then
        z_loop = z_loop - dz
        dz = 1.25d0*(pmin - Pres)/dPdz
        iter =iter + 1
!          if (iter .gt. max_iter) then
!            pres = pmin
!            z_loop = zc
!          endif
      else
        if (pres + dPdz*dz .gt. 10.d0*pmin) then
          pres = pres + dPdz*dz
          dz = 2*dz
          if (abs(dz) .gt. dx_max/dZ_scale) dz = sign(dx_max,dz)/dZ_scale
        else
          pres = pres + dPdz*dz
        endif
      endif
      if (pres .lt. pmin) then
        pres = pmin
        dPdZ = 0.d0
      endif
      iter = 0
      iter_tot = iter_tot + 1
!        if (iter_tot .gt. max_iter) then
!          pres = pmin
!          z_loop = zc
!        endif
    enddo

    ! Recalculate azimuthal velocity
    if (r_sph .lt. Rin) then ! take value at Rin and have it linearly go to zero at Rin0
      v_cyl = max((r_cyl - Rin0)/(Rin - Rin0),0.d0) * &
                  sqrt(GM/(max(Rin,gravity_params(2))*scale_l))/(scale_l/scale_t)
!        exp_arg = 400.d0*(1.d0 - Rin/sqrt(Rin**2+zc**2))*(Rin/R0)**(-beta-1.d0)
!        rho = D0 * (Rin/R0)**(alpha)*exp(-exp_arg)*exp(-hR_clip*(1.d0 - r_cyl/Rin))/scale_d

    else if (r_cyl .lt. Rin) then ! V = Vcirc * sqrt ( (p+q)/400.d0 + (1+q) - qR/r )
      v_cyl = sqrt(GM/(r_sph*scale_l)) * &
            sqrt( max(0.d0,-beta))/(scale_l/scale_t) ! note, the first term assumes beta = -1
    else if (r_cyl .lt. Rout) then ! V = Vcirc * sqrt ( (p+q)/400.d0 + (1+q) - qR/r )
      v_cyl = sqrt(GM/(r_cyl*scale_l)) * &
                  sqrt( max(0.d0,(alpha + beta)/400.d0 + (1+beta) - &
                  beta*r_cyl/r_sph))/(scale_l/scale_t) ! note, the first term assumes beta = -1
    else
      v_cyl = max((Rout0 - r_cyl)/(Rout0 - Rout),0.d0) * &
                  sqrt(GM/(max(Rout,gravity_params(2))*scale_l)) * &
                  sqrt( max(0.d0,(alpha + beta)/400.d0 + (1+beta) - &
                  beta*Rout/max(sqrt(zc**2 + Rout**2),gravity_params(2))))/(scale_l/scale_t)
!       exp_arg = 400.d0*(1.d0 - Rout/sqrt(Rout**2+zc**2))*(Rout/R0)**(-beta-1.d0)
!       rho = D0 * (Rout/R0)**(alpha)*exp(-exp_arg)*exp(-hR_clip*(1.d0 - r_loop/Rout)/ &
!                    (1.d0 - boxlen/1.4142/Rout))/scale_d ! exp(-hR_clip) in corner
    endif
    exp_arg = 400.d0*(1.d0 - r_cyl/r_sph)*(r_cyl/R0)**(-beta-1.d0)
    rho = D0 * (r_cyl/R0)**(alpha)*exp(-exp_arg)/scale_d
    rho = max(rho,D1/scale_d)

    if ( (Rmin_grav > r_cyl ) .or. (r_cyl > Rmax_grav) ) then
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif

  elseif ( disk_setup == 101 ) then ! Gressel, Romberg integration, cylindrical boundary,
    !  transition calculation: integrate to (Rin, Rin, Rin), (constant pressure within Rin)
    !  velocity is circular within boundary, scaling linearly to zero.
    ! Romberg integration, but it's overkill, since the P is mostly analytically. Variations
    !  from analytic can be lumped into v_cyl.
    stop "disk_setup = 101 not recommended. Try 25001 or 26001"
    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! buffer space for tracing R for dPdz

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale

    r_loop = R0
    call Romberg(R0,max(r_cyl,Rin),eps,r_loop,rom_int,dPdR_z0,func_params)
    pres = P0 + rom_int
    rho = D0 * (max(r_cyl,Rin)/R0)**(alpha)/scale_d
    rho = max(rho,D1/scale_d)
    temp = pres/rho

    ! Inside is a boundary condition, so just use the same P. In the dPdz, the calculation knows this
    ! to the left
    call Romberg(max(r_cyl,Rin),max(r_cyl-dr/2.d0,Rin),eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_l = pres + rom_int
    rho_l = D0 * (max(r_cyl-dr/2.d0,Rin)/R0)**(alpha)/scale_d
    rho_l = max(rho_l,D1/scale_d)
    temp_l = pres_l/rho_l

    call Romberg(max(r_cyl,Rin),max(r_cyl+dr/2.d0,Rin),eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_r = pres + rom_int
    rho_r = D0 * (max(r_cyl+dr/2.d0,Rin)/R0)**(alpha)/scale_d
    rho_r = max(rho_r,D1/scale_d)
    temp_r = pres_r/rho_r

    ! Recalculate full density
    if (r_cyl .gt. Rin) then ! key here is R < Rin + dr/2 use the value in this cell.
      exp_arg = min(400.d0*(1.d0 - r_cyl/r_sph)*&
                  (r_cyl/R0)**(-beta-1.d0),disk_exp_limit)
    else
      exp_arg = min(400.d0*(1.d0 - (Rin)/sqrt((Rin)**2+zc**2))*&
                   ((Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    endif
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl,Rin)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
      pres = temp*rho
    endif

    ! Recalculate full density_l
    if (r_cyl-dr/2.d0 .gt. Rin) then
      exp_arg = min(400.d0*(1.d0 - (r_cyl-dr/2.d0)/sqrt((r_cyl-dr/2.d0)**2 + zc**2))*&
                  ((r_cyl-dr/2.d0)/R0)**(-beta-1.d0),disk_exp_limit)
    else
      exp_arg = min(400.d0*(1.d0 - (Rin)/sqrt((Rin)**2+zc**2))*&
                   ((Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    endif
    rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl-dr/2.d0,Rin)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres_l) then
      pres_l =  pres_l + rom_int
    else ! set to midplane temp
      pres_l = temp_l*rho_l
    endif

    ! Recalculate full density_l
    if (r_cyl+dr/2.d0 .gt. Rin) then
      exp_arg = min(400.d0*(1.d0 - (r_cyl+dr/2.d0)/sqrt((r_cyl+dr/2.d0)**2 + zc**2))*&
                  ((r_cyl+dr/2.d0)/R0)**(-beta-1.d0),disk_exp_limit)
    else
      exp_arg = min(400.d0*(1.d0 - (Rin)/sqrt((Rin)**2+zc**2))*&
                   ((Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    endif
    rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl+dr/2.d0,Rin)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres_r) then
      pres_r = max(pres_r + rom_int,pmin)
    else ! set to midplane temp
      pres_r = temp_r*rho_r
    endif

    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    if (r_cyl .lt. Rin) then
      ! want v_cyl = (r_cyl/Rin)*v_cyl(Rin)
      ! dPdR is set up for the boundary, so use it to get v_cyl(Rin)

      !  Just putting as circular for Rin,z that goes to zero at centre.
      v_cyl = r_cyl*sqrt(GM_code/(sqrt(Rin**2 + zc**2))**3)
    else
      dPdR = (pres_r - pres_l) / dR !! Double check sign is right for pos and neg dR --> Pretty sure
      v_cyl = sqrt(max(GM_code*r_cyl**2/r_sph**3 + r_cyl*dPdR/rho,0.d0)) ! do I need to change this for inside Rin? --> Should be self consistent

    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif

  elseif ( disk_setup == 201 ) then ! Gressel, Romberg integration, cylindrical boundary,
    !  transition calculation: integrate to (Rin, Rin, Rin+dr/2), (constant pressure within Rin)
    !  velocity is circular within boundary, scaling linearly to zero.
    !  Save velocity,fluid information in a table that is used to set relaxation and gravity within boundary.
    ! Romberg integration, but it's overkill, since the P is mostly analytically. Variations
    !  from analytic can be lumped into v_cyl.
    stop "disk_setup = 101 not recommended. Try 25001 or 26001"
    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    r_loop = R0
    call Romberg(R0,max(r_cyl,Rin),eps,r_loop,rom_int,dPdR_z0,func_params)
    pres = P0 + rom_int
    rho = D0 * (max(r_cyl,Rin)/R0)**(alpha)/scale_d
    rho = max(rho,D1/scale_d)
    temp = pres/rho

    ! Inside is a boundary condition, so just use the same P. In the dPdz, the calculation knows this
    ! to the left
    if (r_cyl .gt. Rin) then
      call Romberg(r_cyl,max(r_cyl-dr/2.d0,Rin),eps,r_loop,rom_int,dPdR_z0,func_params)
      pres_l = pres + rom_int
      rho_l = D0 * (max(r_cyl-dr/2.d0,Rin)/R0)**(alpha)/scale_d
      rho_l = max(rho_l,D1/scale_d)
      temp_l = pres_l/rho_l
    else
      pres_l = pres
      rho_l = rho
      temp_l = temp
    endif

    call Romberg(max(r_cyl,Rin),max(r_cyl+dr/2.d0,Rin+dr/2.d0),eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_r = pres + rom_int
    rho_r = D0 * (max(r_cyl+dr/2.d0,Rin+dr/2.d0)/R0)**(alpha)/scale_d
    rho_r = max(rho_r,D1/scale_d)
    temp_r = pres_r/rho_r

    ! Recalculate full density
    if (r_cyl .gt. Rin) then
      exp_arg = min(400.d0*(1.d0 - r_cyl/r_sph)*&
                  (r_cyl/R0)**(-beta-1.d0),disk_exp_limit)
    else
      exp_arg = min(400.d0*(1.d0 - Rin/sqrt(Rin**2+zc**2))*&
                   (Rin/R0)**(-beta-1.d0),disk_exp_limit)
    endif
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl,Rin)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
      pres = temp*rho
    endif

    ! Recalculate full density_l
    if (r_cyl .gt. Rin) then
      exp_arg = min(400.d0*(1.d0 - max(r_cyl-dr/2.d0,Rin)/sqrt((max(r_cyl-dr/2.d0,Rin))**2 + zc**2))*&
                  (max(r_cyl-dr/2.d0,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      func_params(9) = max(r_cyl-dr/2.d0,Rin)

      call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
      if (rom_int .gt. pmin-pres_l) then
        pres_l =  pres_l + rom_int
      else ! set to midplane temp
        pres_l = temp_l*rho_l
      endif
    else
      rho_l = rho
      pres_l = pres
    endif

    ! Recalculate full density_r
    exp_arg = min(400.d0*(1.d0 - max(r_cyl+dr/2.d0,Rin+dr/2.d0)/ &
                 sqrt((max(r_cyl+dr/2.d0,Rin+dr/2.d0))**2 + zc**2))*&
                 (max(r_cyl+dr/2.d0,Rin+dr/2.d0)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl+dr/2.d0,Rin+dr/2.d0)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres_r) then
      pres_r = pres_r + rom_int
    else ! set to midplane temp
      pres_r = temp_r*rho_r
    endif

    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    if (r_cyl .lt. Rin) then
      ! want v_cyl = (r_cyl/Rin)*v_cyl(Rin)
      ! dPdR is set up for the boundary, so use it to get v_cyl(Rin)

      !  Actually extending the outer value. Also need to save the value at this heigh in Boundary_table
      dPdR = (pres_r - pres_l) / dr
      v_cyl = sqrt(max(GM_code*Rin**2/(sqrt(Rin**2+zc**2))**3 + &
                      Rin*dPdR/rho,0.d0)) ! do I need to change this for inside Rin? --> Should be self consistent

      v_cyl = min(r_cyl/(disk_vramp*Rin),1.d0)*v_cyl
    else
      dPdR = (pres_r - pres_l) / dr !! (max(r_cyl+dr/2.d0,Rin) - max(r_cyl-dR/2.d0,Rin-dr/2.d0)) !! Double check sign is right for pos and neg dR --> Pretty sure
      v_cyl = sqrt(max(GM_code*r_cyl**2/r_sph**3 + r_cyl*dPdR/rho,0.d0)) ! do I need to change this for inside Rin? --> Should be self consistent

    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif

  elseif ( disk_setup == 1001) then ! Gressel, 2nd derivative extension to central pseudo boundary
    !  Still uses Romberg integration, cylindrical boundary,
    !  transition calculation: integrate to (Rin, Rin, Rin+dr/2), (constant pressure within Rin)
    !  velocity is circular within boundary, scaling linearly to zero.
    !  Save Pressure interpolant coefficients in a table that is used to set gravity and relaxation within boundary.

    ! This method really doesn't solve the boundary issues.
    stop "disk_setup = 1001 not recommended. Try 25001 or 26001"

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive
    r_loop = R0

    ! integrate from reference Radius, R0, to R of cell. Note, if cell is within the interpolation region,
    !   R < Rin+1.5dR, then I need to know what P is doing at the boundary.
    R_integral = max(r_cyl,Rin+1.5d0*dR)
    call Romberg(R0,R_integral,eps,r_loop,rom_int,dPdR_z0,func_params)
    pres = P0 + rom_int
    rho = D0 * (R_integral/R0)**(alpha)/scale_d
    rho = max(rho,D1/scale_d)
    temp = pres/rho

    ! Inside is a boundary condition, so just use the same P. In the dPdz, the calculation knows this
    ! to the left
    R_left = max(r_cyl-dR,Rin+0.5d0*dR)
    call Romberg(R_integral,R_left,eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_l = pres + rom_int
    rho_l = D0 * (R_left/R0)**(alpha)/scale_d
    rho_l = max(rho_l,D1/scale_d)
    temp_l = pres_l/rho_l

    R_right = max(r_cyl+dR,Rin+2.5d0*dr)
    call Romberg(R_integral,R_right,eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_r = pres + rom_int
    rho_r = D0 * (R_right/R0)**(alpha)/scale_d
    rho_r = max(rho_r,D1/scale_d)
    temp_r = pres_r/rho_r

    ! Recalculate full density
    if (r_cyl .eq. R_integral) then
      exp_arg = min(400.d0*(1.d0 - r_cyl/r_sph)*&
                  (r_cyl/R0)**(-beta-1.d0),disk_exp_limit)
    else
      exp_arg = min(400.d0*(1.d0 - (Rin+1.5d0*dR)/sqrt((Rin+1.5d0*dR)**2+zc**2))*&
                   ((Rin+1.5d0*dR)/R0)**(-beta-1.d0),disk_exp_limit)
    endif
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl,Rin+1.5d0*dR)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
      pres = temp*rho
    endif

    ! Recalculate full density_l
    exp_arg = min(400.d0*(1.d0 - max(r_cyl-dR,Rin+0.5d0*dR)/sqrt((max(r_cyl-dR,Rin+0.5d0*dR))**2 + zc**2))*&
                (max(r_cyl-dR,Rin+0.5d0)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl-dr/2.d0,Rin)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres_l) then
      pres_l =  pres_l + rom_int
    else ! set to midplane temp
      pres_l = temp_l*rho_l
    endif

    ! Recalculate full density_r
    exp_arg = min(400.d0*(1.d0 - max(r_cyl+dR,Rin+2.5d0*dR)/ &
                sqrt((max(r_cyl+dR,Rin+2.5d0*dr))**2 + zc**2))*&
                (max(r_cyl+dR,Rin+2.5d0*dR)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    func_params(9) = max(r_cyl+dR,Rin+2.5d0*dR)

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz_R,func_params)
    if (rom_int .gt. pmin-pres_r) then
      pres_r = pres_r + rom_int
    else ! set to midplane temp
      pres_r = temp_r*rho_r
    endif

    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    if (r_cyl .lt. Rin+1.5d0*dR) then
      ! P is an analytic interpolation / extrapolation of pres, pres_l, pres_r. Solve dPdR analytically
      !   Note: In region where dP/dR is zero, need to have v_cyl go to zero so central v is zero:
      !     This gives v_cyl = (r_cyl/Rin)*v_cyl(Rin)


      !
      !  Quadratic fit to the l, cen, r states at x=1,2,3
      Q2 =  (     pres_l - 2.d0*pres +      pres_r)/2.d0
      Q1 = -(5.d0*pres_l - 8.d0*pres + 3.d0*pres_r)/2.d0
      Q0 =  (3.d0*pres_l - 3.d0*pres +      pres_r)

      ! Keep the predicted value at x=0, corresponding to the first boundary cell
      pres_zero   =  max(Q0, pmin)
      ! Artificially move the x=-1 point to x=-2. This leads to a smooth transition to the constant state.
      pres_minus1 =  max(Q2 - Q1 + Q0, pmin)

      !  Refit for cubic, with point y_minus1, y_zero, y(0), y(1), but corresponding to x=-2,0,1,2
      !  Key here is I'm dragging the -1 point over to -2, this makes the 1st derivative approach zero in boundary
      C3 =-(    pres_minus1 -  6.*pres_zero +  8.*pres_l - 3.*pres)/24.
      C2 = (    pres_minus1 -  2.*pres_zero +                 pres)/8.
      C1 =-(    pres_minus1 + 12.*pres_zero - 16.*pres_l + 3.*pres)/12.
      C0 =                        pres_zero

      xscl = ( max(r_cyl,Rin - 2.5d0*dR) - (Rin - 0.5d0*dR) )/dR
      pres = C3*xscl**3 + C2*xscl**2 + C1*xscl + C0
      if (r_cyl .gt. Rin-2.5d0*dR) then
        dPdR = 3.d0*C3*xscl**2 + 2.0d0*C2*xscl + C1
      else
        dPdR = 0.d0
      endif

      rho = D0 * (max(r_cyl,Rin-2.5d0*dR)/R0)**(alpha)/scale_d
      exp_arg = min(400.d0*(1.d0 - max(r_cyl,Rin-2.5d0*dR) / sqrt( (max(r_cyl,Rin-2.5d0*dR))**2 + zc**2 ))*&
                  (max(r_cyl,Rin-2.5d0*dR)/R0)**(-beta-1.d0),disk_exp_limit)
      rho  = max(rho*exp(-exp_arg),D1/scale_d)

      v_cyl = sqrt(max(  GM_code*max(r_cyl,Rin-2.5d0*dR)**2/(sqrt( (max(r_cyl,Rin-2.5d0*dR))**2 + zc**2))**3 + &
                            max(r_cyl,Rin-2.5d0*dR)*dPdR/rho,  0.d0  )) ! do I need to change this for inside Rin? --> Should be self consistent
      if ( r_cyl .lt. Rin - 2.5d0*dR ) then
        v_cyl = v_cyl * r_cyl / (Rin - 2.5d0*dR)
      endif

    else
      dPdR = (pres_r - pres_l) / (2.d0*dR)
      v_cyl = sqrt(max(GM_code*r_cyl**2/r_sph**3 + r_cyl*dPdR/rho,0.d0))
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 2001) then ! Gressel, same as 1001 but with Ttilde = P/rho is interpolated instead.
    ! This method still doesn't solve the boundary issues.
    stop "disk_setup = 1001 not recommended. Try 25001 or 26001"

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive
    r_loop = R0

    ! integrate from reference Radius, R0, to R of cell. Note, if cell is within the interpolation region,
    !   R < Rin+1.5dR, then I need to know what P is doing at the boundary.
    R_integral = max(r_cyl,Rin+1.5d0*dR)
    call Romberg(R0,R_integral,eps,r_loop,rom_int,dPdR_z0,func_params)
    pres = P0 + rom_int
    rho = D0 * (R_integral/R0)**(alpha)/scale_d
    rho = max(rho,D1/scale_d)
    temp = pres/rho

    ! Inside is a boundary condition, so just use the same P. In the dPdz, the calculation knows this
    ! to the left
    R_left = max(r_cyl-dR,Rin+0.5d0*dR)
    call Romberg(R_integral,R_left,eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_l = pres + rom_int
    rho_l = D0 * (R_left/R0)**(alpha)/scale_d
    rho_l = max(rho_l,D1/scale_d)
    temp_l = pres_l/rho_l

    R_right = max(r_cyl+dR,Rin+2.5d0*dr)
    call Romberg(R_integral,R_right,eps,r_loop,rom_int,dPdR_z0,func_params)
    pres_r = pres + rom_int
    rho_r = D0 * (R_right/R0)**(alpha)/scale_d
    rho_r = max(rho_r,D1/scale_d)
    temp_r = pres_r/rho_r


    if (r_cyl .lt. Rin+1.5d0*dR) then
      ! T is an analytic interpolation / extrapolation of pres, pres_l, pres_r. Solve dPdR analytically
      !   Note: In region where dP/dR is zero, need to have v_cyl go to zero so central v is zero:
      !     This gives v_cyl = (r_cyl/Rin)*v_cyl(Rin)

      !
      !  Quadratic fit to the l, cen, r states at x=1,2,3
      Q2 =  (     temp_l - 2.d0*temp +      temp_r)/2.d0
      Q1 = -(5.d0*temp_l - 8.d0*temp + 3.d0*temp_r)/2.d0
      Q0 =  (3.d0*temp_l - 3.d0*temp +      temp_r)

      ! val predicted at x=0, saved for later
      temp_zero   =  max(Q0, Tmini)
      ! val predicted at x=-1, saved for later to be moved to x=-2.
      temp_minus1 =  max(Q2 - Q1 + Q0, Tmini)

      !  Refit for cubic, with point y_minu1, y_zero, y(0), y(1), but corresponding to x=-2,0,1,2
      !  Key here is I'm dragging the -1 point over to -2, this makes the 1st derivative approach zero in boundary
      C3 =-(    temp_minus1 -  6.*temp_zero +  8.*temp_l - 3.*temp)/24.
      C2 = (    temp_minus1 -  2.*temp_zero +                 temp)/8.
      C1 =-(    temp_minus1 + 12.*temp_zero - 16.*temp_l + 3.*temp)/12.
      C0 =                        temp_zero

    endif

    ! Recalculate full density
    if (r_cyl .eq. R_integral) then
      exp_arg = min(400.d0*(1.d0 - r_cyl/r_sph)*&
                  (r_cyl/R0)**(-beta-1.d0),disk_exp_limit)
    else
      exp_arg = min(400.d0*(1.d0 - (Rin+1.5d0*dR)/sqrt((Rin+1.5d0*dR)**2+zc**2))*&
                   ((Rin+1.5d0*dR)/R0)**(-beta-1.d0),disk_exp_limit)
    endif
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    pres = temp*rho

    ! Recalculate full density_l
    exp_arg = min(400.d0*(1.d0 - max(r_cyl-dR,Rin+0.5d0*dR)/sqrt((max(r_cyl-dR,Rin+0.5d0*dR))**2 + zc**2))*&
                (max(r_cyl-dR,Rin+0.5d0)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
    pres_l = temp_l*rho_l

    ! Recalculate full density_r
    exp_arg = min(400.d0*(1.d0 - max(r_cyl+dR,Rin+2.5d0*dR)/ &
                sqrt((max(r_cyl+dR,Rin+2.5d0*dr))**2 + zc**2))*&
                (max(r_cyl+dR,Rin+2.5d0*dR)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    pres_r = temp_r*rho_r

    if (r_cyl .lt. Rin+1.5d0*dR) then
      ! P is an analytic interpolation / extrapolation of pres, pres_l, pres_r. Solve dPdR analytically
      !   Note: In region where dP/dR is zero, need to have v_cyl go to zero so central v is zero:
      !     This gives v_cyl = (r_cyl/Rin)*v_cyl(Rin)

      xscl = ( max(r_cyl,Rin - 2.5d0*dR) - (Rin - 0.5d0*dR) )/dR
      temp = C3*xscl**3 + C2*xscl**2 + C1*xscl + C0
      rho = D0 * (max(r_cyl,Rin-2.5d0*dR)/R0)**(alpha)/scale_d
      exp_arg = min(400.d0*(1.d0 - max(r_cyl,Rin-2.5d0*dR) / sqrt( (max(r_cyl,Rin-2.5d0*dR))**2 + zc**2 ))*&
                  (max(r_cyl,Rin-2.5d0*dR)/R0)**(-beta-1.d0),disk_exp_limit)
      rho  = max(rho*exp(-exp_arg),D1/scale_d)
      pres = temp*rho

      if ( r_cyl .gt. Rin - 2.5d0*dR) then
        dPdR = (3.d0*C3*xscl**2 + 2.0d0*C2*xscl + C1)*rho + &
               temp*rho*(alpha/r_cyl + exp_arg*( (-beta-1.d0)/r_cyl - zc**2/r_sph**3/(1. - r_cyl/r_sph) ) )
        v_cyl = sqrt(max( GM_code*r_cyl**2/r_sph**3 + &
                           r_cyl*dPdR/rho,  0.d0  )) ! do I need to change this for inside Rin? --> Should be self consistent
      else
        dPdR = 0.d0
        R_integral = Rin - 2.5d0*dR
        v_cyl = sqrt(max( GM_code*R_integral**2/(sqrt(R_integral**2+zc**2))**3 + &
                           R_integral*dPdR/rho,  0.d0  )) ! do I need to change this for inside Rin? --> Should be self consistent
        v_cyl = v_cyl * r_cyl / R_integral
      endif

    else
      dPdR = (pres_r - pres_l) / (2.d0*dR)
      v_cyl = sqrt(max(GM_code*r_cyl**2/r_sph**3 + r_cyl*dPdR/rho,0.d0))
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 10001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! Closer, but this method still doesn't solve the boundary issues.
    stop "disk_setup = 10001 not recommended. Try 25001 or 26001"

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(12) = disk_testR
    func_params(13) = disk_testz
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl=sqrt(xc**2+yc**2)
    endif
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    r_loop = R0
    call Romberg(R0,r_cyl,eps,r_loop,rom_int,dPdR2_z0,func_params)
    pres = P0 + rom_int
    rho = D0 * (max(r_cyl,Rin)/R0)**(alpha)/scale_d
    rho = max(rho,D1/scale_d)
    temp = pres/rho

    ! Inside is a boundary condition, so just use the same P. In the dPdz, the calculation knows this
    ! to the left
    r_cyl_l = max(r_cyl-dr,0.d0)
    call Romberg(r_cyl,r_cyl_l,eps,r_loop,rom_int,dPdR2_z0,func_params)
    pres_l = pres + rom_int
    rho_l = D0 * (max(r_cyl_l,Rin)/R0)**(alpha)/scale_d
    rho_l = max(rho_l,D1/scale_d)
    temp_l = pres_l/rho_l

    ! to the right
    r_cyl_r = r_cyl+dr
    call Romberg(r_cyl,r_cyl_r,eps,r_loop,rom_int,dPdR2_z0,func_params)
    pres_r = pres + rom_int
    rho_r = D0 * (max(r_cyl_r,Rin)/R0)**(alpha)/scale_d
    rho_r = max(rho_r,D1/scale_d)
    temp_r = pres_r/rho_r

    ! Recalculate full density
    ! Point of Caution. The cylindrical component above puts Density as if it's at Rin. I think this is fine,
    !  the only points where this matters is high above the plane of the disk, where the exponent sets it Dmin.
    Reff = max(r_cyl, sqrt(Rin**2-zc**2))
    if (r_sph .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph,Rin))*&  ! previously r_cyl/r_sph
                  (Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho = max(rho*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl
    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
      pres = temp*rho
    endif

    ! Recalculate full density_l
    Reff = max(r_cyl_l, sqrt(Rin**2-zc**2))
    r_sph_l = sqrt( r_cyl_l**2 + zc**2)
    if (r_sph_l .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_l,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl_l

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
    if (rom_int .gt. pmin-pres_l) then
      pres_l =  pres_l + rom_int
    else ! set to midplane temp
      pres_l = temp_l*rho_l
    endif

    ! Recalculate full density_r
    Reff = max(r_cyl_r, sqrt(Rin**2-zc**2))
    r_sph_r = sqrt( r_cyl_r**2 + zc**2)
    if (r_sph_r .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_r,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl_r

    call Romberg(0.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
    if (rom_int .gt. pmin-pres_r) then
      pres_r =  pres_r + rom_int
    else ! set to midplane temp
      pres_r = temp_r*rho_r
    endif

    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    dPdR = (pres_r - pres_l) / (r_cyl_r - r_cyl_l)
    if (r_sph .gt. 0.d0) then
      v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + r_cyl*dPdR/rho,0.d0))
    else
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 20001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! I'm a bit worried about doing the first integral in the midplane. Here I want to now integrate:
    !   (R0,0) --> (R0,dx/2) --> (R,dx/2) --> (R,z)

    ! Note, with radial softening, the inner M(<r), and aspect_ratio should change. This is not accounted for here.
    !
    stop "disk_setup = 20001 not quite right. Try 25001 or 26001"

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(12) = disk_testR
    func_params(13) = disk_testz
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl=sqrt(xc**2+yc**2)
    endif
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    ! First: Assume Isothermal at R0:
    exp_arg = 400.d0*(1.d0 - R0/sqrt(R0**2 + (dx/2.d0)**2))
    pres = P0*exp(-exp_arg)

    r_loop = R0
    func_params(9) = dx/2.d0
    if (disk_testR .lt. 0.d0) then
      r_cyl_l = max(r_cyl-dr,0.d0)
      r_cyl_r =     r_cyl+dr
      if (r_cyl .gt. R0) then
        call Romberg(R0,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_l = pres + rom_int
        call Romberg(r_cyl_l,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres = pres_l + rom_int
        call Romberg(r_cyl,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_r = pres + rom_int
      else
        call Romberg(R0,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_r = max(pres + rom_int, pmin)
        if (pres_r .gt. pmin) then
          call Romberg(r_cyl_r,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
          pres = max(pres_r + rom_int,pmin)
          if (pres_r .gt. pmin) then
            call Romberg(r_cyl,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
            pres_l = max(pres + rom_int,pmin)
          else
            pres_l = pmin
          endif
        else
          pres = pmin
          pres_l = pmin
        endif
      endif
    else
      call Romberg(R0,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres = pres + rom_int
      pres_l = pres
      pres_r = pres
    endif

    ! I only use Reff in the exponent of rho, because it doesn't diverge, it goes to zero.
    !  Here I'm just getting the rho at dx/2 above the midplane.
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-(dx/2.d0)**2))
    exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                  (max(r_cyl,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    temp = pres/rho

    if (disk_testR .lt. 0.d0) then
      rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_l, sqrt(Rin**2-(dx/2.d0)**2))
      exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                   (max(r_cyl_l,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      temp_l = pres_l/rho_l

      rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_r, sqrt(Rin**2-(dx/2.d0)**2))
      exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                   (max(r_cyl_r,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
      temp_r = pres_r/rho_r

    else
      rho_l = rho
      temp_l = temp
      rho_r = rho
      temp_r = temp
    endif

    ! Recalculate full density
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-zc**2))
    if (r_sph .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph,Rin))*&  ! previously r_cyl/r_sph
                  (Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho = max(rho*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl
    call Romberg(dx/2.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
      pres = temp*rho
    endif

    ! Recalculate full density_l
    if (disk_testR .lt. 0.d0) then
      rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_l, sqrt(Rin**2-zc**2))
      r_sph_l = sqrt(r_cyl_l**2 + zc**2)
      if (r_sph_l .gt. 0.d0) then
        exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_l,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
        rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      endif
      func_params(9) = r_cyl_l

      call Romberg(dx/2.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
      if (rom_int .gt. pmin-pres_l) then
        pres_l =  pres_l + rom_int
      else ! set to midplane temp
        pres_l = temp_l*rho_l
      endif

      ! Recalculate full density_r
      rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_r, sqrt(Rin**2-zc**2))
      r_sph_r = sqrt(r_cyl_r**2 + zc**2)
      if (r_sph_r .gt. 0.d0) then
        exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_r,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
        rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
      endif
      func_params(9) = r_cyl_r

      call Romberg(dx/2.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
      if (rom_int .gt. pmin-pres_r) then
        pres_r =  pres_r + rom_int
      else ! set to midplane temp
        pres_r = temp_r*rho_r
      endif
    else
      rho_l = rho
      rho_r = rho
      pres_l = pres
      pres_r = pres
    endif
    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    dPdR = (pres_r - pres_l) / (r_cyl_r - r_cyl_l)
    if (disk_testR .lt. 0.d0) then
      if (r_sph .gt. 0.d0) then
        v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + r_cyl*dPdR/rho,0.d0))
      else
        v_cyl = 0.d0
      endif
    else
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 21001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! I'm a bit worried about doing the first integral in the midplane. Here I want to now integrate:
    !   (R0,0) --> (R0,dx/2) --> (R,dx/2) --> (R,z)

    ! Same as 20001 except in the radial direction using P analytic solution and only integreating vertically.

    ! Note, with radial softening, the inner M(<r), and aspect_ratio should change. This is not accounted for here.
    !
    stop "disk_setup = 21001 not quite right. Try 25001, or 26001"

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(12) = disk_testR
    func_params(13) = disk_testz
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen = func_params(2)/func_params(5)

    ! Note, self-grav off
    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl=sqrt(xc**2+yc**2)
    endif
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    ! I only use Reff in the exponent of rho, because it doesn't diverge, it goes to zero.
    !  Here I'm just getting the rho at dx/2 above the midplane.
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-(dx/2.d0)**2))
    exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                  (max(r_cyl,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    pres = (rho/max(R_cyl,Rin))*P_fact_cen
    temp = pres/rho

    if (disk_testR .lt. 0.d0) then
      rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_l, sqrt(Rin**2-(dx/2.d0)**2))
      exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                   (max(r_cyl_l,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      pres_l = (rho_l/Reff)*P_fact_cen
      temp_l = pres_l/rho_l

      rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_r, sqrt(Rin**2-(dx/2.d0)**2))
      exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                   (max(r_cyl_r,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
      pres_r = (rho_r/Reff)*P_fact_cen
      temp_r = pres_r/rho_r

    else
      rho_l = rho
      pres_l = pres
      temp_l = temp
      rho_r = rho
      pres_r = pres
      temp_r = temp
    endif

    ! Recalculate full density
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-zc**2))
    if (r_sph .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph,Rin))*&  ! previously r_cyl/r_sph
                  (Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho = max(rho*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl
    call Romberg(dx/2.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
      pres = temp*rho
    endif

    ! Recalculate full density_l
    if (disk_testR .lt. 0.d0) then
      rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_l, sqrt(Rin**2-zc**2))
      r_sph_l = sqrt(r_cyl_l**2 + zc**2)
      if (r_sph_l .gt. 0.d0) then
        exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_l,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
        rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      endif
      func_params(9) = r_cyl_l

      call Romberg(dx/2.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
      if (rom_int .gt. pmin-pres_l) then
        pres_l =  pres_l + rom_int
      else ! set to midplane temp
        pres_l = temp_l*rho_l
      endif

      ! Recalculate full density_r
      rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_r, sqrt(Rin**2-zc**2))
      r_sph_r = sqrt(r_cyl_r**2 + zc**2)
      if (r_sph_r .gt. 0.d0) then
        exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_r,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
        rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
      endif
      func_params(9) = r_cyl_r

      call Romberg(dx/2.d0,zc,eps,z_loop,rom_int,dPdz2_R,func_params)
      if (rom_int .gt. pmin-pres_r) then
        pres_r =  pres_r + rom_int
      else ! set to midplane temp
        pres_r = temp_r*rho_r
      endif
    else
      rho_l = rho
      rho_r = rho
      pres_l = pres
      pres_r = pres
    endif
    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    dPdR = (pres_r - pres_l) / (r_cyl_r - r_cyl_l)
    if (disk_testR .lt. 0.d0) then
      if (r_sph .gt. 0.d0) then
        v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + r_cyl*dPdR/rho,0.d0))
      else
        v_cyl = 0.d0
      endif
    else
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 22001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! I'm a bit worried about doing the first integral in the midplane. Here I want to now integrate:
    !   (R0,0) --> (R0,dx/2) --> (R,dx/2) --> (R,z)

    ! Whole thing is analytic now. This is basically correct, except the inner M(<r),
    !  and aspect_ratio should change. This is not accounted for here.
    !
    stop "disk_setup = 22001 not quite right. Try 25001, or 26001"

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(12) = disk_testR
    func_params(13) = disk_testz
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen = func_params(2)/func_params(5)
    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    else
      D1 = D0*exp(-disk_exp_limit)
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    endif

    ! Note, self-grav off
    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl=sqrt(xc**2+yc**2)
    endif
    r_sph=sqrt(r_cyl**2+zc**2)

    Reff = max(r_cyl, sqrt(max(Rin**2-zc**2,0.d0)))
    ! Recalculate full density - first midplane
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d

    max_z = maxz_fact * Reff
    ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
    ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
    dztrans = ztrans2 - ztrans1
    if (abs(zc) .le. ztrans1) then
      if (r_sph .gt. 0.d0) then
        exp_arg = 400.d0*(1.d0 - Reff/max(r_sph,Rin))*(Reff/R0)**(-beta-1.d0)
        rho = rho*exp(-exp_arg)
      endif
    else if (abs(zc) .ge. ztrans2) then
      rho = rho*exp(-disk_exp_limit)
    else
      rtrans1 = sqrt(Reff**2 + ztrans1**2)
      ftrans1 = -400.d0*(1.d0 - Reff/rtrans1)
      delta_trans = ftrans1 + disk_exp_limit

      theta_trans = -400.d0 * Reff * ztrans1 / rtrans1 ** 3 * dztrans
        eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

      k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
      k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
      k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

      dz_here = ztrans2 - abs(zc)

      exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
      rho = rho*exp(-exp_arg)
    endif

    pres = rho/max(R_cyl,Rin) * P_fact_cen


    if (disk_testR .lt. 0.d0) then
      Reff = max(r_cyl - dr, sqrt(max(Rin**2-zc**2,0.d0)))
      ! Recalculate full density - first midplane
      rho_l = D0 * (max(R_cyl-dr,Rin)/R0)**(alpha)/scale_d

      max_z = maxz_fact * Reff
      ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
      ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
      dztrans = ztrans2 - ztrans1
      if (abs(zc) .le. ztrans1) then
        if (r_sph .gt. 0.d0) then
          exp_arg = 400.d0*(1.d0 - Reff/max(r_sph,Rin))*(Reff/R0)**(-beta-1.d0)
          rho_l = rho_l*exp(-exp_arg)
        endif
      else if (abs(zc) .ge. ztrans1) then
        rho_l = rho_l*exp(-disk_exp_limit)
      else
        rtrans1 = sqrt(Reff**2 + ztrans1**2)
        ftrans1 = -400.d0*(1.d0 - Reff/rtrans1)
        delta_trans = ftrans1 + disk_exp_limit

        theta_trans = -400.d0 * Reff * ztrans1 / rtrans1 ** 3 * dztrans
          eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

        k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
        k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
        k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

        dz_here = ztrans2 - abs(zc)

        exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
        rho_l = rho_l*exp(-exp_arg)
      endif

      pres_l = rho_l/max(R_cyl-dR,Rin) * P_fact_cen

      Reff = max(r_cyl+dR, sqrt(max(Rin**2-zc**2,0.d0)))
      ! Recalculate full density - first midplane
      rho_r = D0 * (max(R_cyl+dR,Rin)/R0)**(alpha)/scale_d

      max_z = maxz_fact * Reff
      ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
      ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
      dztrans = ztrans2 - ztrans1
      if (abs(zc) .le. ztrans1) then
        if (r_sph .gt. 0.d0) then
          exp_arg = 400.d0*(1.d0 - Reff/max(r_sph,Rin))*(Reff/R0)**(-beta-1.d0)
          rho_r = rho_r*exp(-exp_arg)
        endif
      else if (abs(zc) .ge. ztrans1) then
        rho_r = rho_r*exp(-disk_exp_limit)
      else
        rtrans1 = sqrt(Reff**2 + ztrans1**2)
        ftrans1 = -400.d0*(1.d0 - Reff/rtrans1)
        delta_trans = ftrans1 + disk_exp_limit

        theta_trans = -400.d0 * Reff * ztrans1 / rtrans1 ** 3 * dztrans
          eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

        k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
        k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
        k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

        dz_here = ztrans2 - abs(zc)

        exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
        rho_r = rho_r*exp(-exp_arg)
      endif

      pres_r = rho_r/max(R_cyl+dR,Rin) * P_fact_cen
    else
      pres_l = pres
      pres_r = pres
    endif


    dPdR = (pres_r - pres_l) / (2*dr)
    if (disk_testR .lt. 0.d0) then
      if (r_sph .gt. 0.d0) then
        v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + r_cyl*dPdR/rho,0.d0))
      else
        v_cyl = 0.d0
      endif
    else
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 23001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! I'm a bit worried about doing the first integral in the midplane. Here I want to now integrate:
    !   (R0,0) --> (R0,dx/2) --> (R,dx/2) --> (R,z)

    ! 23001 --> Still setting P = rho/R GM/alpha2, but subsampling within cell to get more average value.
    !  --> question on whether this would be in HSE though.
    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(12) = disk_testR
    func_params(13) = disk_testz
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen = func_params(2)/func_params(5)
    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    else
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    endif

    ! Note, self-grav off
    dx_ss = dx/2.d0**disk_nsubsample
    dR = dx_ss

    dens_tot = 0.d0
    mom_tot = 0.d0
    p_tot = 0.d0

    do kk=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
      zc = zc0 + kk*dx_ss/2.d0
      do jj=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
        yc = yc0 + jj*dx_ss/2.d0
        do ii=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
          xc = xc0 + ii*dx_ss/2.d0

          if (disk_testR .gt. 0.d0) then
            xc = disk_testR
            yc = 0.d0
            r_cyl = disk_testR
          elseif (disk_testZ .gt. 0.d0) then
            r_cyl=sqrt(xc**2+yc**2)
            zc = disk_testz
          else
            r_cyl=sqrt(xc**2+yc**2)
          endif
          r_sph=sqrt(r_cyl**2+zc**2)

          Reff = max(r_cyl, sqrt(max(Rin**2-zc**2,0.d0)))
          ! Recalculate full density - first midplane
          rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d

          max_z = maxz_fact * Reff
          ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
          ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
          dztrans = ztrans2 - ztrans1
          if (abs(zc) .le. ztrans1) then
            if (r_sph .gt. 0.d0) then
              exp_arg = 400.d0*(1.d0 - Reff/max(r_sph,Rin))*(Reff/R0)**(-beta-1.d0)
              rho = rho*exp(-exp_arg)
            endif
          else if (abs(zc) .ge. ztrans2) then
            rho = rho*exp(-disk_exp_limit)
          else
            rtrans1 = sqrt(Reff**2 + ztrans1**2)
            ftrans1 = -400.d0*(1.d0 - Reff/rtrans1)
            delta_trans = ftrans1 + disk_exp_limit

            theta_trans = -400.d0 * Reff * ztrans1 / rtrans1 ** 3 * dztrans
              eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

            k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
            k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
            k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

            dz_here = ztrans2 - abs(zc)

            exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
            rho = rho*exp(-exp_arg)
          endif

          pres = rho/max(R_cyl,Rin) * P_fact_cen


          if (disk_testR .lt. 0.d0) then
            Reff = max(r_cyl - dr, sqrt(max(Rin**2-zc**2,0.d0)))
            ! Recalculate full density - first midplane
            rho_l = D0 * (max(R_cyl-dr,Rin)/R0)**(alpha)/scale_d

            max_z = maxz_fact * Reff
            ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
            ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
            dztrans = ztrans2 - ztrans1
            if (abs(zc) .le. ztrans1) then
              if (r_sph .gt. 0.d0) then
                exp_arg = 400.d0*(1.d0 - Reff/max(r_sph,Rin))*(Reff/R0)**(-beta-1.d0)
                rho_l = rho_l*exp(-exp_arg)
              endif
            else if (abs(zc) .ge. ztrans2) then
              rho_l = rho_l*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(Reff**2 + ztrans1**2)
              ftrans1 = -400.d0*(1.d0 - Reff/rtrans1)
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -400.d0 * Reff * ztrans1 / rtrans1 ** 3 * dztrans
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(zc)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho_l = rho_l*exp(-exp_arg)
            endif

            pres_l = rho_l/max(R_cyl-dR,Rin) * P_fact_cen

            Reff = max(r_cyl+dR, sqrt(max(Rin**2-zc**2,0.d0)))
            ! Recalculate full density - first midplane
            rho_r = D0 * (max(R_cyl+dR,Rin)/R0)**(alpha)/scale_d

            max_z = maxz_fact * Reff
            ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
            ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
            dztrans = ztrans2 - ztrans1
            if (abs(zc) .le. ztrans1) then
              if (r_sph .gt. 0.d0) then
                exp_arg = 400.d0*(1.d0 - Reff/max(r_sph,Rin))*(Reff/R0)**(-beta-1.d0)
                rho_r = rho_r*exp(-exp_arg)
              endif
            else if (abs(zc) .ge. ztrans2) then
              rho_r = rho_r*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(Reff**2 + ztrans1**2)
              ftrans1 = -400.d0*(1.d0 - Reff/rtrans1)
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -400.d0 * Reff * ztrans1 / rtrans1 ** 3 * dztrans
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(zc)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho_r = rho_r*exp(-exp_arg)
            endif

            pres_r = rho_r/max(R_cyl+dR,Rin) * P_fact_cen
          else
            pres_l = pres
            pres_r = pres
          endif


          dPdR = (pres_r - pres_l) / (2.d0*dr)
          if (disk_testR .lt. 0.d0) then
            if (r_sph .gt. 0.d0) then
              v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + disk_shear_forcing*r_cyl*dPdR/rho,0.d0))
            else
              v_cyl = 0.d0
            endif
          else
            v_cyl = 0.d0
          endif


          ! Add to totals
          dens_tot = dens_tot + rho
          p_tot = p_tot + pres
          if (r_cyl .gt. 0.d0) then
            mom_tot(1) = mom_tot(1) - v_cyl*yc/r_cyl*rho
            mom_tot(2) = mom_tot(2) + v_cyl*xc/r_cyl*rho
          endif

          ! Finish subsample loops
        enddo ! ii
      enddo ! jj
    enddo ! kk

    ! average totals
    mom_tot = mom_tot / dens_tot
    rho = dens_tot/8**disk_nsubsample
    pres = p_tot/8**disk_nsubsample

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    vel(1)= mom_tot(1) + v_turb(1)
    vel(2)= mom_tot(2) + v_turb(2)
    vel(3)= v_turb(3)

  elseif ( disk_setup == 23101 ) then ! Gressel, but with the same inner/outer boundary behavior as Hennebelle, Lesur, & Teyssier 2017 (HLT17)
    ! Note, inner hydro boundary has Rin set T constant, not rho. Gravity softening, params(2), sets that softening, which need not be the same
    ! Specifically for HLT17, alpha = -2.5, beta = -1, disk_aspect=0.158114

    ! 23101 --> Still setting P = rho/R GM/alpha2, but subsampling within cell to get more average value.

    soft = gravity_params(2)      ! inside soft mass will be assumed constant: F(r) = -GMtot/Rsoft^3 * r
    Rmin_grav = gravity_params(6) ! inside Rmin_grav, grav is set to zero. Note, Rmin_grav < 0 means cylindrical r
    Rmax_grav = gravity_params(7) ! outside Rmax_grav, grav is set to zero. Note, Rmax_grav < 0 means cylindrical r

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen0 = GM_code*disk_aspect**2  ! P = Pfact * rho / R, so Pfact = cs_0^2 * Rd = disk_aspect**2*G*Mstar, see eqn 2 in HLT17.
    !! Note, P_fact_cen is only constant if soft < Rin. Otherwise need the factor (Rin/soft)**3
    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact0 = sqrt( 1.d0/(1.d0 - disk_exp_limit*disk_aspect**2)**2 - 1.d0)
    else
      maxz_fact0 = sqrt( 1.d0/(1.d0 - disk_exp_limit*disk_aspect**2)**2 - 1.d0)
    endif

    ! Note, self-grav off
    dx_ss = dx/2.d0**disk_nsubsample
    dR = dx_ss

    dens_tot = 0.d0
    mom_tot = 0.d0
    p_tot = 0.d0

    do kk=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
      zc = zc0 + kk*dx_ss/2.d0
      do jj=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
        yc = yc0 + jj*dx_ss/2.d0
        do ii=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
          xc = xc0 + ii*dx_ss/2.d0

          if (disk_testR .gt. 0.d0) then
            xc = disk_testR
            yc = 0.d0
            r_cyl = disk_testR
          elseif (disk_testZ .gt. 0.d0) then
            r_cyl=sqrt(xc**2+yc**2)
            zc = disk_testz
          else
            r_cyl=sqrt(xc**2+yc**2)
          endif
          r_sph=sqrt(r_cyl**2+zc**2)

          ! Recalculate full density - first midplane
          Reff = max(R_cyl,dR*2.d0)
          rho = D0 * (Reff/R0)**(alpha)/scale_d ! setting minimum R to dr/2

          if (Reff .lt. Rin) then
            aspect_ratio = disk_aspect*sqrt(Reff/Rin)
            maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*aspect_ratio**2)**2 - 1.d0)
            P_fact_cen = P_fact_cen0*(Reff/Rin)
          else
            aspect_ratio = disk_aspect
            maxz_fact = maxz_fact0
            P_fact_cen = P_fact_cen0
          endif

          if (Reff .lt. Rout2) then
            max_z = maxz_fact * Reff
            if (disk_smooth_vertRho) then
              ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
              ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
              dztrans = ztrans2 - ztrans1
            else
               ztrans1 = max_z
               ztrans2 = max_z
            endif

            if (abs(zc) .le. ztrans1) then
              if (r_sph .gt. 0.d0) then
                exp_arg = (1.d0 - Reff/max(r_sph,Reff))*(Reff/R0)**(-beta-1.d0)/aspect_ratio**2
                rho = rho*exp(-exp_arg)
              endif
            else if (abs(zc) .ge. ztrans2) then
              rho = rho*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(Reff**2 + ztrans1**2)
              ftrans1 = -(1.d0 - Reff/rtrans1)/aspect_ratio**2
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -Reff * ztrans1 / rtrans1 ** 3 * dztrans / aspect_ratio**2
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(zc)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho = rho*exp(-exp_arg)
            endif
            rho = max(rho,D1/scale_D)
            if (Reff .gt. Rout) rho = max(rho / 100.d0, D1/scale_D)
          else
            rho = D1/scale_D
          endif

          pres = rho/Reff * P_fact_cen
          if (soft .gt. Rin) then
            pres = pres * (Rin/soft)**3
          endif

          ! Get left/right states to make approximation of P gradient
          if (disk_testR .lt. 0.d0 .and. Reff .lt. Rout2) then
            Reff = max(r_cyl - dR, dR*2.d0)

            if (Reff .lt. Rin) then
              aspect_ratio = disk_aspect*sqrt(Reff/Rin)
              maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*aspect_ratio**2)**2 - 1.d0)
              P_fact_cen = P_fact_cen0*(Reff/Rin)
            else
              aspect_ratio = disk_aspect
              maxz_fact = maxz_fact0
              P_fact_cen = P_fact_cen0
            endif
            ! Recalculate full density - first midplane
            rho_l = D0 * (Reff/R0)**(alpha)/scale_d

            max_z = maxz_fact * Reff
            if (disk_smooth_vertRho) then
              ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
              ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
              dztrans = ztrans2 - ztrans1
            else
               ztrans1 = max_z
               ztrans2 = max_z
            endif

            if (abs(zc) .le. ztrans1) then
              if (r_sph .gt. 0.d0) then
                exp_arg = (1.d0 - Reff/max(r_sph,Reff))*(Reff/R0)**(-beta-1.d0)/aspect_ratio**2
                rho_l = rho_l*exp(-exp_arg)
              endif
            else if (abs(zc) .ge. ztrans2) then
              rho_l = rho_l*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(Reff**2 + ztrans1**2)
              ftrans1 = -(1.d0 - Reff/rtrans1)/aspect_ratio**2
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -Reff * ztrans1 / rtrans1 ** 3 * dztrans/aspect_ratio**2
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(zc)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho_l = rho_l*exp(-exp_arg)
            endif
            rho_l = max(rho_l,D1/scale_D)
            if (Reff .gt. Rout) rho_l = max(rho_l / 100.d0, D1/scale_D)

            pres_l = rho_l/Reff * P_fact_cen

            Reff = max(r_cyl+dR, 2.d0*dR)
            if (Reff .lt. Rin) then
              aspect_ratio = disk_aspect*sqrt(Reff/Rin)
              maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*aspect_ratio**2)**2 - 1.d0)
              P_fact_cen = P_fact_cen0*(Reff/Rin)
            else
              aspect_ratio = disk_aspect
              maxz_fact = maxz_fact0
              P_fact_cen = P_fact_cen0
            endif

            ! Recalculate full density - first midplane
            rho_r = D0 * (Reff/R0)**(alpha)/scale_d

            max_z = maxz_fact * Reff
            if (disk_smooth_vertRho) then
              ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
              ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
              dztrans = ztrans2 - ztrans1
            else
              ztrans1 = max_z
              ztrans2 = max_z
            endif
            if (abs(zc) .le. ztrans1) then
              if (r_sph .gt. 0.d0) then
                exp_arg = (1.d0 - Reff/max(r_sph,Reff))*(Reff/R0)**(-beta-1.d0)/aspect_ratio**2
                rho_r = rho_r*exp(-exp_arg)
              endif
            else if (abs(zc) .ge. ztrans2) then
              rho_r = rho_r*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(Reff**2 + ztrans1**2)
              ftrans1 = -(1.d0 - Reff/rtrans1)/aspect_ratio**2
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -Reff * ztrans1 / rtrans1 ** 3 * dztrans/aspect_ratio**2
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(zc)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho_r = rho_r*exp(-exp_arg)
            endif
            rho_r = max(rho_r,D1/scale_D)
            if (Reff .gt. Rout) rho_r = max(rho_r / 100.d0, D1/scale_D)

            pres_r = rho_r/Reff * P_fact_cen
          else
            pres_l = pres
            pres_r = pres
          endif


          dPdR = (pres_r - pres_l) / (2.d0*dr)
          if (disk_testR .lt. 0.d0) then
            if (r_sph .gt. 0.d0) then
              v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + disk_shear_forcing*r_cyl*dPdR/rho,0.d0))
            else
              v_cyl = 0.d0
            endif
          else
            v_cyl = 0.d0
          endif


          ! Add to totals
          dens_tot = dens_tot + rho
          p_tot = p_tot + pres
          if (r_cyl .gt. 0.d0) then
            mom_tot(1) = mom_tot(1) - v_cyl*yc/r_cyl*rho
            mom_tot(2) = mom_tot(2) + v_cyl*xc/r_cyl*rho
          endif

    ! Finish subsample loops
        enddo ! ii
      enddo ! jj
    enddo ! kk

    ! average totals
    mom_tot = mom_tot / dens_tot
    rho = dens_tot/8**disk_nsubsample
    pres = p_tot/8**disk_nsubsample

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    vel(1)= mom_tot(1) + v_turb(1)
    vel(2)= mom_tot(2) + v_turb(2)
    vel(3)= v_turb(3)


  elseif ( disk_setup == 24001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed so v has continuous 2nd derivatives.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun and going to zero at centre.

    ! Similar, I interpolate P and rho to be constant third derivatives. This puts limit on how much the dens and pres can increase.

    ! Polynomials for r <= R are:
    !   v(r) = v(R) * [  (99/8)(r/R)^3  -  (77/4)(r/R)^4  + (63/8)(r/R)^5  ] == v(R) * P5(r/R)
    !     --> same as if M(r<R) = Mtot * (r/R) * P5(r/R)^2
    !   rho & pres = rho(R) & pres(R) * [ Delta - P7(r/R;Delta,alpha) ]
    !          --> rho : Delta = 2 ; alpha = -1.5
    !               --> P7(r/R;2,-1.5) = (31/16)(r/R)^4 + (69/16)(r/R)^5 - (143/16)(r/R)^6 + (59/16)(r/R)^7
    !          --> Pres : Delta = 3 ; alpha = -2.5
    !               --> P7(r/R;3,-2.5) = (275/48)(r/R)^4 + (87/16)(r/R)^5 - (265/16)(r/R)^6 + (355/48)(r/R)^7

    ! 24001 --> Still setting P = rho/R GM/alpha2, but subsampling within cell to get more average value.
    !  --> Crux here is as soon as you soften, the P ~ rho/R * const is wrong within Rin, because M(<r) is a function of z.
    !  ---> But... it's not locally a strong function of z.
    !  --> OK, but if I use Mtot to integrate up for r<R then the grav is right for z > Reff, but too high below.
    !              if I use Mint(r=R), then grav is too low at high z ... BUT THIS DOESN'T MATTER BECAUSE I HIT CONST P ANYWAY!
    ! What this means is P_fact should be a function of r, following Mtot given above: P(x) * P5(x)^2.
    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen0 = ( GM_code )/400.d0

    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    else
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    endif

    ! Note, self-grav off
    dx_ss = dx/2.d0**disk_nsubsample
    dR = dx_ss

    dens_tot = 0.d0
    mom_tot = 0.d0
    p_tot = 0.d0

    do kk=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
      zc = zc0 + kk*dx_ss/2.d0
      do jj=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
        yc = yc0 + jj*dx_ss/2.d0
        do ii=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
          xc = xc0 + ii*dx_ss/2.d0

    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl=sqrt(xc**2+yc**2)
    endif
    r_sph=sqrt(r_cyl**2+zc**2)

    rcyl_ratio = r_cyl/Rin

    if (rcyl_ratio .ge. 1.d0) then
      ! Recalculate full density - first midplane
      rho = max(D0 * (R_cyl/R0)**(alpha),D1)/scale_d

      max_z = maxz_fact * R_cyl
      ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
      ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
      dztrans = ztrans2 - ztrans1
      if (abs(zc) .le. ztrans1) then
        if (r_sph .gt. 0.d0) then
          exp_arg = 400.d0*(1.d0 - R_cyl/r_sph)*(R_cyl/R0)**(-beta-1.d0)
          rho = rho*exp(-exp_arg)
        endif
      else if (abs(zc) .ge. ztrans2) then
        rho = rho*exp(-disk_exp_limit)
      else
        rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
        ftrans1 = -400.d0*(1.d0 - R_cyl/rtrans1)
        delta_trans = ftrans1 + disk_exp_limit

        theta_trans = -400.d0 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
          eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

        k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
        k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
        k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

        dz_here = ztrans2 - abs(zc)

        exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
        rho = rho*exp(-exp_arg)
      endif
      rho = max(rho,D1/scale_D)
      pres = rho/R_cyl * P_fact_cen0
      v_cyl = sqrt(GM/(r_cyl*scale_l)) * &
                sqrt( max(0.d0,(alpha + beta)/400.d0 + (1+beta) - &
                beta*r_cyl/r_sph))/(scale_l/scale_t)
    else
      rsph_ratio = r_sph/Rin
      P_fact_cen = P_fact_cen0 * Mass_of_r(rcyl_ratio) * rcyl_ratio

      z1 = zc/Z_of_r(rcyl_ratio)
      rsph1 = sqrt(Rin**2 + z1**2)
      ! Recalculate full density - first midplane
      rho = max(D0 * (Rin/R0)**(alpha),D1)/scale_d ! rho(R=Rin,0)
      pres1 = rho/Rin * P_fact_cen0 ! P(Rin,0)

      rho10 = rho ! rho(Rin,0) saved for later
      rho0  = max(rho*Rho_of_r(rcyl_ratio),D1/scale_D) ! rho(R<Rin,0)
      pres0 = pres1*Pres_of_r(rcyl_ratio) ! pres(R<Rin,0)
      max_z = maxz_fact * Rin

      ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
      ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
      dztrans = ztrans2 - ztrans1

      if (abs(z1) .le. ztrans1) then
        if (rsph1 .gt. 0.d0) then
          exp_arg = 400.d0*(1.d0 - Rin/rsph1)*(Rin/R0)**(-beta-1.d0)
          rho = rho*exp(-exp_arg)
        endif
      else if (abs(z1) .ge. ztrans2) then
        rho = rho*exp(-disk_exp_limit)
      else
        rtrans1 = sqrt(Rin**2 + ztrans1**2)
        ftrans1 = -400.d0*(1.d0 - Rin/rtrans1)
        delta_trans = ftrans1 + disk_exp_limit

        theta_trans = -400.d0 * Rin * ztrans1 / rtrans1 ** 3 * dztrans
          eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

        k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
        k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
        k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

        dz_here = ztrans2 - abs(z1)

        exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
        rho = rho*exp(-exp_arg)
      endif

      rho1 = max(rho,D1/scale_D) ! rho(Rin,z_tilde)
      if (abs(z1) .ge. ztrans2) then
        rho = rho1
      else
        rho  = max(rho*Rho_of_r(rcyl_ratio),D1/scale_D)
      endif
!  Can't do the originally smooth P version, since it sets P above a background, so high T above disk.
!   However, can raise it a bit to try and soften central shear. Basically, pres0 and rho0 terms are scaled by disk_raiseP_factor
      delta_pres = (disk_raiseP_factor*rho0 - rho)*P_fact_cen/R_cyl ! Note P_fact_cen has a R_cyl/Rin factor
      pres = disk_raiseP_factor*pres0 - delta_pres

      rho = rho * Mass_of_r(rcyl_ratio)/Mass_of_r(rsph_ratio)
!          pres = rho*P_fact_cen/R_cyl

      ! Now to get velocity. This is TRICKY!

      ! Something below isn't right ... it gives a strong kink in v_cyl at Rin ... not sure why yet.
      stop "24001 gives a weird v_cyl solution, so I recommend using 26001"
      if (r_sph .gt. 0.d0) then
        vcyl2 = GM*rsph_ratio*Mass_of_r(rsph_ratio)*r_cyl**2/r_sph**3/scale_l
        vcyl2 = vcyl2/(scale_l/scale_t)**2
      else
        vcyl2 = 0.d0
      endif

      ! Gradient in midplane P
!  NO LONGER APPLIES
!  Unless disk_raiseP_factor is > 0
      vcyl2 = vcyl2 - disk_raiseP_factor*(r_cyl/rho)*(pres1*PPrime_of_r(rcyl_ratio)/Rin)

      ! Radial Gradient in rho(R,0) - rho(R,z)
! NO LONGER TAKE OFF MIDPLANE VALUE
! Unless disk_raiseP_factor > 0
      vcyl2 = vcyl2 - P_fact_cen/rho * (rho1 - disk_raiseP_factor*rho10)*RhoPrime_of_R(rcyl_ratio)/Rin
!          vcyl2 = vcyl2 - P_fact_cen/rho * rho1*RhoPrime_of_R(rcyl_ratio)/Rin

      ! Gradient in the exponent term since z = f(R)
! This term remains
      vcyl2 = vcyl2 + P_fact_cen * 400.d0*z1**3/zc/rsph1**3 * ZPrime_of_r(rcyl_ratio) / (Mass_of_r(rcyl_ratio)/Mass_of_r(rsph_ratio))

      ! Finally gradient in M(<r)
! Remove midplane value
! Unless disk_raiseP_factor > 0
      vcyl2 = vcyl2 - (r_cyl/rho)*delta_pres/Rin*MPrime_of_r(rcyl_ratio)/Mass_of_r(rcyl_ratio)
!          vcyl2 = vcyl2 - (r_cyl/rho)*pres/Rin*MPrime_of_r(rcyl_ratio)/Mass_of_r(rcyl_ratio)

      ! Get vcyl
      v_cyl = sqrt(max(0.d0,vcyl2)) ! Double check gradient in rho(R,0) term and units

    endif

    ! Add to totals
    dens_tot = dens_tot + rho
    p_tot = p_tot + pres
    if (r_cyl .gt. 0.d0) then
      mom_tot(1) = mom_tot(1) - v_cyl*yc/r_cyl*rho
      mom_tot(2) = mom_tot(2) + v_cyl*xc/r_cyl*rho
    endif

    ! Finish subsample loops
        enddo ! ii
      enddo ! jj
    enddo ! kk

    ! average totals
    mom_tot = mom_tot / dens_tot
    rho = dens_tot/8**disk_nsubsample
    pres = p_tot/8**disk_nsubsample

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    vel(1)= mom_tot(1) + v_turb(1)
    vel(2)= mom_tot(2) + v_turb(2)
    vel(3)= v_turb(3)

    pmin = min(pres,pmin)

  elseif ( disk_setup == 25001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed so v has continuous 2nd derivatives.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun and going to zero at centre.

    ! Similar, I interpolate P and rho to be constant third derivatives. This puts limit on how much the dens and pres can increase.

    ! Polynomials for r <= R are:
    !   v(r) = v(R) * [  (99/8)(r/R)^3  -  (77/4)(r/R)^4  + (63/8)(r/R)^5  ] == v(R) * P5(r/R)
    !     --> same as if M(r<R) = Mtot * (r/R) * P5(r/R)^2
    !   rho & pres = rho(R) & pres(R) * [ Delta - P7(r/R;Delta,alpha) ]
    !          --> rho : Delta = 2 ; alpha = -1.5
    !               --> P7(r/R;2,-1.5) = (31/16)(r/R)^4 + (69/16)(r/R)^5 - (143/16)(r/R)^6 + (59/16)(r/R)^7
    !          --> Pres : Delta = 3 ; alpha = -2.5
    !               --> P7(r/R;3,-2.5) = (275/48)(r/R)^4 + (87/16)(r/R)^5 - (265/16)(r/R)^6 + (355/48)(r/R)^7

    ! 25001 --> Still setting P = rho/R GM/alpha2, but subsampling within cell to get more average value.
    !  --> Crux here is as soon as you soften, the P ~ rho/R * const is wrong within Rin, because M(<r) is a function of z.
    !  ---> But... it's not locally a strong function of z.
    !  --> OK, but if I use Mtot to integrate up for r<R then the grav is right for z > Reff, but too high below.
    !              if I use Mint(r=R), then grav is too low at high z ... BUT THIS DOESN'T MATTER BECAUSE I HIT CONST P ANYWAY!
    ! What this means is P_fact should be a function of r, following Mtot given above: P(x) * P5(x)^2.
    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen0 = ( GM_code )/400.d0

    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    else
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit/400.d0)**2 - 1.d0)
    endif

    ! Note, self-grav off
    dx_ss = dx/2.d0**disk_nsubsample
    dR = dx_ss

    dens_tot = 0.d0
    mom_tot = 0.d0
    p_tot = 0.d0

    do kk=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
      zc = zc0 + kk*dx_ss/2.d0
      do jj=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
        yc = yc0 + jj*dx_ss/2.d0
        do ii=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
          xc = xc0 + ii*dx_ss/2.d0

    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl_l = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl_l=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl_l=max(sqrt(xc**2+yc**2),dR)
    endif

    do i_Rad=-1,3,2 ! want to do i_rad=0 last so values are saved at end of loop.
      r_cyl = max(r_cyl_l + i_Rad*dR,dR)
      if (i_Rad .eq. 3) r_cyl = r_cyl_l

      r_sph=sqrt(r_cyl**2+zc**2)

      rcyl_ratio = r_cyl/Rin

    if (rcyl_ratio .ge. 1.d0) then
      ! Recalculate full density - first midplane
      rho = max(D0 * (R_cyl/R0)**(alpha),D1)/scale_d

      max_z = maxz_fact * R_cyl
      ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
      ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
      dztrans = ztrans2 - ztrans1
      if (abs(zc) .le. ztrans1) then
        if (r_sph .gt. 0.d0) then
          exp_arg = 400.d0*(1.d0 - R_cyl/r_sph)*(R_cyl/R0)**(-beta-1.d0)
          rho = rho*exp(-exp_arg)
        endif
      else if (abs(zc) .ge. ztrans2) then
        rho = rho*exp(-disk_exp_limit)
      else
        rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
        ftrans1 = -400.d0*(1.d0 - R_cyl/rtrans1)
        delta_trans = ftrans1 + disk_exp_limit

        theta_trans = -400.d0 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
          eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

        k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
        k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
        k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

        dz_here = ztrans2 - abs(zc)

        exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
        rho = rho*exp(-exp_arg)
      endif

      rho = max(rho,D1/scale_D)
      pres = rho/R_cyl * P_fact_cen0

    else
      rsph_ratio = r_sph/Rin
      P_fact_cen = P_fact_cen0 * Mass_of_r(rcyl_ratio) * rcyl_ratio

!          z1 = zc/Z_of_r(rcyl_ratio)
!          rsph1 = sqrt(Rin**2 + z1**2)
      ! Note pinching flaring anymore ...
      z1 = max(abs(zc),dR)
      rsph1 = max(r_sph,sqrt(r_cyl**2 + z1**2))

      ! Recalculate full density - first midplane
      rho = max(D0 * (Rin/R0)**(alpha)/scale_d, D1/scale_D) ! rho(R=Rin,0)
      pres1 = rho/Rin * P_fact_cen0 ! P(Rin,0)

      rho10 = rho ! rho(Rin,0) saved for later
      rho0  = max(rho*Rho_of_r(rcyl_ratio), D1/scale_D) ! rho(R<Rin,0)
      pres0 = pres1*Pres_of_r(rcyl_ratio) ! pres(R<Rin,0)
!          max_z = maxz_fact * Rin
      max_z = maxz_fact * R_cyl

      ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
      ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
      dztrans = ztrans2 - ztrans1

      if (abs(z1) .le. ztrans1) then
        if (rsph1 .gt. 0.d0) then
          exp_arg = 400.d0*(1.d0 - R_cyl/rsph1)*(Rin/R0)**(-beta-1.d0)
          rho = rho*exp(-exp_arg)
        endif
      else if (abs(z1) .ge. ztrans2) then
        rho = rho*exp(-disk_exp_limit)
      else
        rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
        ftrans1 = -400.d0*(1.d0 - R_cyl/rtrans1)
        delta_trans = ftrans1 + disk_exp_limit

        theta_trans = -400.d0 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
          eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

        k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
        k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
        k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

        dz_here = ztrans2 - abs(z1)

        exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
        rho = rho*exp(-exp_arg)
      endif

      rho1 = max(rho,D1/scale_D) ! rho(Rin,z_tilde)
      rho  = max(rho*Rho_of_r(rcyl_ratio),D1/scale_d)

!  Can't do the originally smooth P version, since it sets P above a background, so high T above disk.
!   However, can raise it a bit to try and soften central shear. Basically, pres0 and rho0 terms are scaled by disk_raiseP_factor
      delta_pres = (disk_raiseP_factor*rho0 - rho)*P_fact_cen/max(R_cyl,dR) ! Note P_fact_cen has a R_cyl/Rin factor
      pres = disk_raiseP_factor*pres0 - delta_pres
      rho = rho*Mass_of_r(rcyl_ratio)/Mass_of_r(rsph_ratio) ! This actually corrects it for analytic HSE.

!          pres = rho*P_fact_cen/R_cyl

    endif

    if (i_Rad .eq. -1) then
      pres_l = pres
    endif
    if (i_Rad .eq. 1) then
      pres_r = pres
    endif

! Done loop
    enddo

    ! Now to get Pressure gradient
    dPdR = (pres_R - pres_L)/(r_cyl + dR - max(r_cyl - dR,dR) )

    if (disk_testR .lt. 0.d0) then
      if (r_sph .gt. 0.d0) then
        v_cyl = sqrt(max(GM_code*Mass_of_R(r_sph/Rin)*r_cyl**2/r_sph**2/max(Rin,r_sph) + disk_shear_forcing*r_cyl*dPdR/rho,0.d0))
      else
        v_cyl = 0.d0
      endif
    else
      v_cyl = 0.d0
    endif



    ! Add to totals
    dens_tot = dens_tot + rho
    p_tot = p_tot + pres
    if (r_cyl .gt. 0.d0) then
      mom_tot(1) = mom_tot(1) - v_cyl*yc/r_cyl*rho
      mom_tot(2) = mom_tot(2) + v_cyl*xc/r_cyl*rho
    endif

    ! Finish subsample loops
        enddo ! ii
      enddo ! jj
    enddo ! kk

    ! average totals
    mom_tot = mom_tot / dens_tot
    rho = dens_tot/8**disk_nsubsample
    pres = p_tot/8**disk_nsubsample

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    vel(1)= mom_tot(1) + v_turb(1)
    vel(2)= mom_tot(2) + v_turb(2)
    vel(3)= v_turb(3)

    pmin = min(pres,pmin)

  elseif ( disk_setup == 26001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed so v has continuous 2nd derivatives.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun and going to zero at centre.

    ! Similar, I interpolate P and rho to be constant third derivatives. This puts limit on how much the dens and pres can increase.

    ! Polynomials for r <= R are:
    !   v(r) = v(R) * [  (99/8)(r/R)^3  -  (77/4)(r/R)^4  + (63/8)(r/R)^5  ] == v(R) * P5(r/R)
    !     --> same as if M(r<R) = Mtot * (r/R) * P5(r/R)^2
    !   rho & pres = rho(R) & pres(R) * [ Delta - P7(r/R;Delta,alpha) ]
    !          --> rho : Delta = 2 ; alpha = -1.5
    !               --> P7(r/R;2,-1.5) = (31/16)(r/R)^4 + (69/16)(r/R)^5 - (143/16)(r/R)^6 + (59/16)(r/R)^7
    !          --> Pres : Delta = 3 ; alpha = -2.5
    !               --> P7(r/R;3,-2.5) = (275/48)(r/R)^4 + (87/16)(r/R)^5 - (265/16)(r/R)^6 + (355/48)(r/R)^7

    ! 26001 --> Still setting P = rho/R GM/alpha2, but subsampling within cell to get more average value.

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen0 = ( GM_code )*disk_aspect**2

    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*disk_aspect**2)**2 - 1.d0)
    else
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*disk_aspect**2)**2 - 1.d0)
    endif

    ! Note, self-grav off
    dx_ss = dx/2.d0**disk_nsubsample
    dR = dx_ss

    dens_tot = 0.d0
    mom_tot = 0.d0
    p_tot = 0.d0

    do kk=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
      zc = zc0 + kk*dx_ss/2.d0
      do jj=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
        yc = yc0 + jj*dx_ss/2.d0
        do ii=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
          xc = xc0 + ii*dx_ss/2.d0

          if (disk_testR .gt. 0.d0) then
            xc = disk_testR
            yc = 0.d0
            r_cyl_l = disk_testR
          elseif (disk_testZ .gt. 0.d0) then
            r_cyl_l=sqrt(xc**2+yc**2)
            zc = disk_testz
          else
            r_cyl_l=max(sqrt(xc**2+yc**2),dR)
          endif

          do i_Rad=-1,3,2 ! want to do i_rad=0 last so values are saved at end of loop.
            r_cyl = max(r_cyl_l + i_Rad*dR,dR)
            if (i_Rad .eq. 3) r_cyl = r_cyl_l

            r_sph=sqrt(r_cyl**2+zc**2)

            rcyl_ratio = r_cyl/Rin

          if (rcyl_ratio .ge. 1.d0) then
            ! Recalculate full density - first midplane
            rho = max(D0 * (R_cyl/R0)**(alpha),D1)/scale_d

            max_z = maxz_fact * R_cyl
            ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
            ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
            dztrans = ztrans2 - ztrans1
            if (abs(zc) .le. ztrans1) then
              if (r_sph .gt. 0.d0) then
                exp_arg = 400.d0*(1.d0 - R_cyl/r_sph)*(R_cyl/R0)**(-beta-1.d0)
                rho = rho*exp(-exp_arg)
              endif
            else if (abs(zc) .ge. ztrans2) then
              rho = rho*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
              ftrans1 = -400.d0*(1.d0 - R_cyl/rtrans1)
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -400.d0 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(zc)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho = rho*exp(-exp_arg)
            endif

            rho = max(rho,D1/scale_D)
            pres = rho/R_cyl * P_fact_cen0

          else
            rsph_ratio = r_sph/Rin
            P_fact_cen = P_fact_cen0 * Mass_of_r(rcyl_ratio) * rcyl_ratio

  !          z1 = zc/Z_of_r(rcyl_ratio)
  !          rsph1 = sqrt(Rin**2 + z1**2)
            ! Note pinching flaring anymore ...
            z1 = max(abs(zc),dR)
            rsph1 = max(r_sph,sqrt(r_cyl**2 + z1**2))

            ! Recalculate full density - first midplane
            rho = max(D0 * (Rin/R0)**(alpha)/scale_d, D1/scale_D) ! rho(R=Rin,0)
            pres1 = rho/Rin * P_fact_cen0 ! P(Rin,0)

            rho10 = rho ! rho(Rin,0) saved for later
            rho0  = max(rho*Rho_of_r(rcyl_ratio), D1/scale_D) ! rho(R<Rin,0)
            pres0 = pres1*Pres_of_r(rcyl_ratio) ! pres(R<Rin,0)
  !          max_z = maxz_fact * Rin

            alpha2 = 400.d0*(1.d0 - 0.75d0*(1 - Rcyl_ratio)**3)**2
            maxz_fact2 = sqrt( 1.d0/(1.d0 - disk_exp_limit/alpha2)**2 - 1.d0)

            max_z = maxz_fact2 * R_cyl

            ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
            ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
            dztrans = ztrans2 - ztrans1

            if (abs(z1) .le. ztrans1) then
              if (rsph1 .gt. 0.d0) then
                exp_arg = alpha2*(1.d0 - R_cyl/rsph1)*(R_cyl/R0)**(-beta-1.d0)
                rho = rho*exp(-exp_arg)
              endif
            else if (abs(z1) .ge. ztrans2) then
              rho = rho*exp(-disk_exp_limit)
            else
              rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
              ftrans1 = -alpha2*(1.d0 - R_cyl/rtrans1)
              delta_trans = ftrans1 + disk_exp_limit

              theta_trans = -alpha2 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
                eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

              k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
              k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
              k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

              dz_here = ztrans2 - abs(z1)

              exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
              rho = rho*exp(-exp_arg)
            endif

            rho1 = max(rho,D1/scale_D) ! rho(Rin,z_tilde)
            rho  = max(rho*Rho_of_r(rcyl_ratio),D1/scale_d)

  !  Can't do the originally smooth P version, since it sets P above a background, so high T above disk.
  !   However, can raise it a bit to try and soften central shear. Basically, pres0 and rho0 terms are scaled by disk_raiseP_factor
            delta_pres = (disk_raiseP_factor*rho0 - rho)*P_fact_cen*400.d0/alpha2/max(R_cyl,dR) ! Note P_fact_cen has a R_cyl/Rin factor
            pres = disk_raiseP_factor*pres0 - delta_pres
            rho = rho*Mass_of_r(rcyl_ratio)/Mass_of_r(rsph_ratio)

  !          pres = rho*P_fact_cen/R_cyl

          endif

          if (i_Rad .eq. -1) then
            pres_l = pres
          endif
          if (i_Rad .eq. 1) then
            pres_r = pres
          endif

  ! Done loop
          enddo

          ! Now to get Pressure gradient
          dPdR = (pres_R - pres_L)/(r_cyl + dR - max(r_cyl - dR,dR) )

          if (disk_testR .lt. 0.d0) then
            if (r_sph .gt. 0.d0) then
              v_cyl = sqrt(max(GM_code*Mass_of_R(r_sph/Rin)*r_cyl**2/r_sph**2/max(Rin,r_sph) + disk_shear_forcing*r_cyl*dPdR/rho,0.d0))
            else
              v_cyl = 0.d0
            endif
          else
            v_cyl = 0.d0
          endif



          ! Add to totals
          dens_tot = dens_tot + rho
          p_tot = p_tot + pres
          if (r_cyl .gt. 0.d0) then
            mom_tot(1) = mom_tot(1) - v_cyl*yc/r_cyl*rho
            mom_tot(2) = mom_tot(2) + v_cyl*xc/r_cyl*rho
          endif

          ! Finish subsample loops
        enddo ! ii
      enddo ! jj
    enddo ! kk

    ! average totals
    mom_tot = mom_tot / dens_tot
    rho = dens_tot/8**disk_nsubsample
    pres = p_tot/8**disk_nsubsample

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    vel(1)= mom_tot(1) + v_turb(1)
    vel(2)= mom_tot(2) + v_turb(2)
    vel(3)= v_turb(3)

    pmin = min(pres,pmin)

  elseif ( disk_setup == 30001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! Here I'm going to do vertical HSE directly, similar to disk_setup 11.
    !  In midplane this is hard, Going to first go from (R0,0) --> (R0,dx/2) assuming constant T.
    !   Then Romberg integrate to (R,dx/2)
    !   Then:
    !     Let P(z=dx/2)    = P1           dPdz(z=dx/2)    = dP1
    !         P(z=dx+dx/2) = P2 etc.      dPdz(z=dx+dx/2) = dP2 etc.
    !         P(-z) = P(z)
    !
    !         e.g.,  dP2 = (P3 - P1)/2dz
    !                dP_j = (P_j+1 - P_j-1)/2dz
    !         Except: dP1 = (P2 - P1)/2dz
    !            where P1 is already known.

    !         So P2 = (2dz)*dP1 + P1 (recall dP_j is also known from HSE alone).
    !            P3 = (2dz)*dP2 + P1
    !            P4 = (2dz)*dP3 + P2 = (2dz)*(dP3+dP1) + P1
    !            P5 = (2dz)*dP4 + P3 = (2dz)*(dP4+dP2) + P1
    !            P6 = (2dz)*dP5 + P4 = (2dz)*(dP5+dP3+dP1) + P1
    !            P_j = (2dz)*dP_j-1 + P_j-2 = (2dz)*(dP_j-1 + dP_j-3) + P_j-4

    !
    ! so iiz = {1,2},iz-1,2

    !         So P_j, j even     = (2dz)*sum(1 <= k odd  < j) [dP_k] + P1
    !            P_j, j odd > 1  = (2dz)*sum(2 <= k even < j) [dP_k] + P1

    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(12) = disk_testR
    func_params(13) = disk_testz
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    if (disk_testR .gt. 0.d0) then
      xc = disk_testR
      yc = 0.d0
      r_cyl = disk_testR
    elseif (disk_testZ .gt. 0.d0) then
      r_cyl=sqrt(xc**2+yc**2)
      zc = disk_testz
    else
      r_cyl=sqrt(xc**2+yc**2)
    endif
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    ! First: Assume Isothermal at R0:
    exp_arg = 400.d0*(1.d0 - R0/sqrt(R0**2 + (dx/2.d0)**2))
    pres = P0*exp(-exp_arg)

    r_loop = R0
    func_params(9) = dx/2.d0
    if (disk_testR .lt. 0.d0) then
      r_cyl_l = max(r_cyl-dr,0.d0)
      r_cyl_r =     r_cyl+dr
      if (r_cyl .gt. R0) then
        call Romberg(R0,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_l = pres + rom_int
        call Romberg(r_cyl_l,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres = pres_l + rom_int
        call Romberg(r_cyl,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_r = pres + rom_int
      else
        call Romberg(R0,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_r = max(pres + rom_int, pmin)
        if (pres_r .gt. pmin) then
          call Romberg(r_cyl_r,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
          pres = max(pres_r + rom_int,pmin)
          if (pres_r .gt. pmin) then
            call Romberg(r_cyl,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
            pres_l = max(pres + rom_int,pmin)
          else
            pres_l = pmin
          endif
        else
          pres = pmin
          pres_l = pmin
        endif
      endif
    else
      call Romberg(R0,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres = pres + rom_int
      pres_l = pres
      pres_r = pres
    endif

    ! I only use Reff in the exponent of rho, because it doesn't diverge, it goes to zero.
    !  Here I'm just getting the rho at dx/2 above the midplane.
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-(dx/2.d0)**2))
    exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                  (max(r_cyl,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    temp = pres/rho

    if (disk_testR .lt. 0.d0) then
      rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_l, sqrt(Rin**2-(dx/2.d0)**2))
      exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                    (max(r_cyl_l,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      temp_l = pres_l/rho_l

      rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_r, sqrt(Rin**2-(dx/2.d0)**2))
      exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                    (max(r_cyl_r,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
      rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
      temp_r = pres_r/rho_r
    else
      rho_l = rho
      rho_r = rho
      temp_l = temp
      temp_r = temp
    endif

    ! Recalculate full density
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-zc**2))
    if (r_sph .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph,Rin))*&  ! previously r_cyl/r_sph
                  (Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho = max(rho*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl

    iz = int((abs(zc) + dx/2.d0 + 0.1d0*dx)/dx)

!         So P_j, j even     = (2dz)*sum(1 <= k odd  < j) [dP_k] + P1
!            P_j, j odd > 1  = (2dz)*sum(2 <= k even < j) [dP_k] + P1
!   z_j = (j-1/2)dx = (2j-1)dx/2 --> j = int(z_j/dx + 1/2) = int( (z_j + dx/2)/dx )
    if (mod(iz,2) == 0) then
      pres_int = pres/(2.d0*dx) ! P @ 0
      pres_intm1 = pres_int     ! P @ 1

      ! P @ cell 2
      pres_int = pres_int + dPdz2_R(dx/2.d0, func_params)

      do iiz=3,iz-1,2
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz-3)*dx/2.d0, func_params), pres_int  ) ! Podd  += dPdz @ iiz-1
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1) ! Peven += dPdz @ iiz
      enddo
      pres = pres_int*2.d0*dx

    else
      pres_intm1 = pres/(2.d0*dx) ! P @ 0
      pres_int = pres_intm1       ! P @ 1

      ! P @ cell 2
      pres_intm1 = pres_intm1 + dPdz2_R(dx/2.d0, func_params)

      do iiz=2,iz-1,2
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1) ! Podd  += dPdz @ iiz
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz+1)*dx/2.d0, func_params), pres_int  ) ! Peven += dPdz @ iiz+1
      enddo
      pres = 2.d0*dx*pres_int
    endif

    if (pres .lt. pmin) then
      pres = pmin
    endif

    if (disk_testR .lt. 0.d0) then
      ! Recalculate full density_l
      rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_l, sqrt(Rin**2-zc**2))
      r_sph_l = sqrt(r_cyl_l**2 + zc**2)
      if (r_sph_l .gt. 0.d0) then
        exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_l,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
        rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
      endif
      func_params(9) = r_cyl_l

      if (mod(iz,2) == 0) then
        pres_int = pres_l/(2.d0*dx)
        pres_intm1 = pres_int

        ! P @ cell 2
        pres_int = pres_int + dPdz2_R(dx/2.d0, func_params)

        do iiz=3,iz-1,2
          pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz-3)*dx/2.d0, func_params), pres_int  )
          pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
        enddo

        pres_l = 2.d0*dx*pres_int
      else
        pres_int = pres_l/(2.d0*dx)
        pres_intm1 = pres_int

        ! P @ cell 2
        pres_intm1 = pres_intm1 + dPdz2_R(dx/2.d0, func_params)

        do iiz=2,iz-1,2
          pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
          pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz+1)*dx/2.d0, func_params), pres_int  )
        enddo

        pres_l = 2.d0*dx*pres_int
      endif

      if (pres_l .lt. pmin) then
        pres_l = pmin
      endif

      ! Recalculate full density_r
      rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
      Reff = max(r_cyl_r, sqrt(Rin**2-zc**2))
      r_sph_r = sqrt(r_cyl_r**2 + zc**2)
      if (r_sph_r .gt. 0.d0) then
        exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_r,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
        rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
      endif
      func_params(9) = r_cyl_r

      if (mod(iz,2) == 0) then
        pres_int = pres_r/(2.d0*dx)
        pres_intm1 = pres_int

        ! P @ cell 2
        pres_int = pres_int + dPdz2_R(dx/2.d0, func_params)

        do iiz=3,iz-1,2
          pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz-3)*dx/2.d0, func_params), pres_int  )
          pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
        enddo

        pres_r = 2.d0*dx*pres_int
      else
        pres_int = pres_r/(2.d0*dx)
        pres_intm1 = pres_int

        ! P @ cell 2
        pres_intm1 = pres_intm1 + dPdz2_R(dx/2.d0, func_params)

        do iiz=2,iz-1,2
          pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
          pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz+1)*dx/2.d0, func_params), pres_int  )
        enddo

        pres_r = 2.d0*dx*pres_int
      endif

      if (pres_r .lt. pmin) then
        pres_r = pmin
      endif
    else
      pres_l = pres
      pres_r = pres
    endif

    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    if (disk_testR .lt. 0.d0) then
      dPdR = (pres_r - pres_l) / (r_cyl_r - r_cyl_l)
      if (r_sph .gt. 0.d0) then
        v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + r_cyl*dPdR/rho,0.d0))
      else
        v_cyl = 0.d0
      endif
    else
      dPdR = 0.d0
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif

  elseif ( disk_setup == 40001 ) then ! Gressel
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed spherically with uniform density.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun linearly going to zero.

    ! Here I'm going to do vertical HSE directly, similar to disk_setup 11.
    !  In midplane this is hard, Going to first go from (R0,0) --> (R0,dx/2) and (R0,3dx2/) assuming constant T.
    !   Then Romberg integrate to (R,dx/2) and (R,3dx/2)
    !   Then:
    !     Let P(z=dx/2)    = P1           dPdz(z=dx/2)    = dP1
    !         P(z=3dx/2)   = P2 etc.      dPdz(z=3dx/2)   = dP2 etc.
    !         P(-z) = P(z)
    !
    !         e.g.,  dP2 = (P3 - P1)/2dz
    !                dP_j = (P_j+1 - P_j-1)/2dz
    !         Except: dP1 = (P2 - P1)/2dz
    !            where P1 is already known.

    !            P3 = (2dz)*dP2 + P1
    !            P4 = (2dz)*dP3 + P2  --> But ensure monotonic
    !            P5 = (2dz)*dP4 + P3 = (2dz)*(dP4+dP2) + P1
    !            P6 = (2dz)*dP5 + P4 = (2dz)*(dP5+dP3+dP1) + P1
    !            P_j = (2dz)*dP_j-1 + P_j-2 = (2dz)*(dP_j-1 + dP_j-3) + P_j-4

    !
    ! so iiz = {1,2},iz-1,2

    !         So P_j, j even     = (2dz)*sum(1 <= k odd  < j) [dP_k] + P1
    !            P_j, j odd > 1  = (2dz)*sum(2 <= k even < j) [dP_k] + P1

    ! Question, what happens once one P_j hits Pmin?
    !
    func_params(1) = D0/scale_d
    func_params(2) = GM_code
    func_params(3) = disk_exp_limit
    func_params(4) = R0
    func_params(5) = 400.d0
    func_params(6) = alpha
    func_params(7) = beta
    func_params(8) = abar
    func_params(9) = R0
    func_params(10) = D1/scale_d
    func_params(11) = Rin
    func_params(14) = grav_angle
    func_params(15) = grav_width

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    ! First: Assume Isothermal at R0:
    exp_arg = 400.d0*(1.d0 - R0/sqrt(R0**2 + (dx/2.d0)**2))
    pres = P0*exp(-exp_arg)

    exp_arg = 400.d0*(1.d0 - R0/sqrt(R0**2 + (3.d0*dx/2.d0)**2))
    pres_p1 = P0*exp(-exp_arg)

    ! Pres +- dR
    r_loop = R0
    func_params(9) = dx/2.d0
    r_cyl_l = max(r_cyl-dr,0.d0)
    r_cyl_r =     r_cyl+dr
    if (r_cyl .gt. R0) then
      call Romberg(R0,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_l = pres + rom_int
      call Romberg(r_cyl_l,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres = pres_l + rom_int
      call Romberg(r_cyl,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_r = pres + rom_int
    else
      call Romberg(R0,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_r = max(pres + rom_int, pmin)
      if (pres_r .gt. pmin) then
        call Romberg(r_cyl_r,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres = max(pres_r + rom_int,pmin)
        if (pres .gt. pmin) then
          call Romberg(r_cyl,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
          pres_l = max(pres + rom_int,pmin)
        else
          pres_l = pmin
        endif
      else
        pres = pmin
        pres_l = pmin
      endif
    endif

    ! Pres_p1 +- dR
    r_loop = R0
    func_params(9) = 3.d0*dx/2.d0
    r_cyl_l = max(r_cyl-dr,0.d0)
    r_cyl_r =     r_cyl+dr
    if (r_cyl .gt. R0) then
      call Romberg(R0,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_l_p1 = pres_p1 + rom_int
      call Romberg(r_cyl_l,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_p1 = pres_l_p1 + rom_int
      call Romberg(r_cyl,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_r_p1 = pres_p1 + rom_int
    else
      call Romberg(R0,r_cyl_r,eps,r_loop,rom_int,dPdR2_z,func_params)
      pres_r_p1 = max(pres_p1 + rom_int, pmin)
      if (pres_r_p1 .gt. pmin) then
        call Romberg(r_cyl_r,r_cyl,eps,r_loop,rom_int,dPdR2_z,func_params)
        pres_p1 = max(pres_r_p1 + rom_int,pmin)
        if (pres_r_p1 .gt. pmin) then
          call Romberg(r_cyl,r_cyl_l,eps,r_loop,rom_int,dPdR2_z,func_params)
          pres_l_p1 = max(pres_p1 + rom_int,pmin)
        else
          pres_l_p1 = pmin
        endif
      else
        pres_p1 = pmin
        pres_l_p1 = pmin
      endif
    endif


    ! I only use Reff in the exponent of rho, because it doesn't diverge, it goes to zero.
    !  Here I'm just getting the rho at dx/2 above the midplane.
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-(dx/2.d0)**2))
    exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                  (max(r_cyl,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    rho = max(rho*exp(-exp_arg),D1/scale_d)
    temp = pres/rho


    rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl_l, sqrt(Rin**2-(dx/2.d0)**2))
    exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                  (max(r_cyl_l,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
    temp_l = pres_l/rho_l


    rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl_r, sqrt(Rin**2-(dx/2.d0)**2))
    exp_arg = min(400.d0*(1.d0 - Reff/sqrt(Reff**2 + (dx/2.d0)**2))*&  ! previously r_cyl/r_sph
                  (max(r_cyl_r,Rin)/R0)**(-beta-1.d0),disk_exp_limit)
    rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    temp_r = pres_r/rho_r


    ! Recalculate full density
    rho = D0 * (max(R_cyl,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl, sqrt(Rin**2-zc**2))
    if (r_sph .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph,Rin))*&  ! previously r_cyl/r_sph
                  (Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho = max(rho*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl

    iz = int((abs(zc) + dx/2.d0 + 0.1d0*dx)/dx)

!         So P_j, j even     = (2dz)*sum(1 <= k odd  < j) [dP_k] + P1
!            P_j, j odd > 1  = (2dz)*sum(2 <= k even < j) [dP_k] + P1
!   z_j = (j-1/2)dx = (2j-1)dx/2 --> j = int(z_j/dx + 1/2) = int( (z_j + dx/2)/dx )
    if (mod(iz,2) == 0) then
      pres_intm1 = pres/(2.d0*dx)  ! P1
      pres_int = pres_p1/(2.d0*dx) ! P2

      do iiz=3,iz-1,2
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz-3)*dx/2.d0, func_params), pres_int  )
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
      enddo
      pres = pres_int*2.d0*dx

    else
      pres_int = pres/(2.d0*dx) ! P1
      pres_intm1 = pres_p1/(2.d0*dx) ! P2

      do iiz=2,iz-1,2
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz+1)*dx/2.d0, func_params), pres_int  )
      enddo
      pres = 2.d0*dx*pres_int
    endif

    if (pres .lt. pmin) then
      pres = pmin
    endif

    ! Recalculate full density_l
    rho_l = D0 * (max(R_cyl_l,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl_l, sqrt(Rin**2-zc**2))
    r_sph_l = sqrt(r_cyl_l**2 + zc**2)
    if (r_sph_l .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_l,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho_l = max(rho_l*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl_l

    if (mod(iz,2) == 0) then
      pres_int = pres_l_p1/(2.d0*dx) ! P2
      pres_intm1 = pres_l/(2.d0*dx)  ! P1

      do iiz=3,iz-1,2
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz-3)*dx/2.d0, func_params), pres_int  )
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
      enddo
      pres_l = pres_int*2.d0*dx

    else
      pres_int = pres_l/(2.d0*dx) ! P1
      pres_intm1 = pres_l_p1/(2.d0*dx) ! P2

      do iiz=2,iz-1,2
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz+1)*dx/2.d0, func_params), pres_int  )
      enddo
      pres_l = 2.d0*dx*pres_int
    endif

    if (pres_l .lt. pmin) then
      pres_l = pmin
    endif

    ! Recalculate full density_r
    rho_r = D0 * (max(R_cyl_r,Rin)/R0)**(alpha)/scale_d
    Reff = max(r_cyl_r, sqrt(Rin**2-zc**2))
    r_sph_r = sqrt(r_cyl_r**2 + zc**2)
    if (r_sph_r .gt. 0.d0) then
      exp_arg = min(400.d0*(1.d0 - Reff/max(r_sph_r,Rin))*(Reff/R0)**(-beta-1.d0),disk_exp_limit)
      rho_r = max(rho_r*exp(-exp_arg),D1/scale_d)
    endif
    func_params(9) = r_cyl_r

    if (mod(iz,2) == 0) then
      pres_int = pres_r_p1/(2.d0*dx) ! P2
      pres_intm1 = pres_r/(2.d0*dx)  ! P1

      do iiz=3,iz-1,2
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz-3)*dx/2.d0, func_params), pres_int  )
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
      enddo
      pres_r = pres_int*2.d0*dx

    else
      pres_int = pres_r/(2.d0*dx) ! P1
      pres_intm1 = pres_r_p1/(2.d0*dx) ! P2

      do iiz=2,iz-1,2
        pres_int   = min(pres_int   + dPdz2_R((2*iiz-1)*dx/2.d0, func_params), pres_intm1)
        pres_intm1 = min(pres_intm1 + dPdz2_R((2*iiz+1)*dx/2.d0, func_params), pres_int  )
      enddo
      pres_r = 2.d0*dx*pres_int
    endif

    if (pres_r .lt. pmin) then
      pres_r = pmin
    endif

    pres   =min(pres  ,disk_Textra*temp  *rho  )
    pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
    pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    dPdR = (pres_r - pres_l) / (r_cyl_r - r_cyl_l)
    if (r_sph .gt. 0.d0) then
      v_cyl = sqrt(max(GM_code*r_cyl**2/max(r_sph,Rin)**3 + r_cyl*dPdR/rho,0.d0))
    else
      v_cyl = 0.d0
    endif

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif


  elseif ( disk_setup == 2 ) then ! Bai

    rr=sqrt(xc**2+yc**2)

    rho= D0 * (rr/R0)**(-alpha) * f(theta)
    vel = 0.d0
    vel(3) = ( GM - (alpha + 1) * T0 * R0 * g(theta) ) / rr
    pres = T0 * (R0/rr) * g(theta) * rho

  else if ( disk_setup == 3 ) then ! NoB Kepler ring
    gravity_params(1) = GM*scale_t**2/scale_l**3 ! GM in code units
    gravity_params(2) = 7e10/scale_l ! softening
    ! grav_params(3:5) = position of sun = (0,0,0) = default grav_params value
    rr=sqrt((xc-R0)**2+yc**2)

    if (rr < Rout) then
      rho= D0/scale_d
    else
      rho= 1d-16/scale_d
    endif
    vel= 0.d0
    ! Some question about whether v_z should be R * Omega , or r * Omega. I'm using the first below, but multiply by rr/xc if you change
    vel(3)= sqrt(GM/(xc*scale_l))*scale_t/scale_l
    pres= 8.314472d7 * T0 * D0 / abar /(scale_l**2 / scale_t**2)

  else if ( disk_setup == 4 ) then ! Weird ramp for tesdting
    gravity_params(1) = GM*scale_t**2/scale_l**3 ! GM in code units
    gravity_params(2) = 7e10/scale_l ! softening
    ! grav_params(3:5) = position of sun = (0,0,0) = default grav_params value
    rho= D0*(1.d0 + xc)/scale_d
    vel= 0.d0
    vel(3)= V0
    pres= 8.314472d7 * T0 * D0 / abar /(scale_l**2 / scale_t**2)

  elseif (disk_setup == 5) then ! Galaxy gas is sech z, sech R, background DM, star particles are only bulge
    ! Need to duplicate gravana. For now doing analytic bulge too.

    ! At large R where sech is small, sech ~ 2/exp(x)
    !  --> D0 * 2 / exp(r_cyl/R1) = D1
    !  --> r_cyl ~ R1 * log(2D0/D1)
    !    Use half this to be sure
      R0 = 0.5d0 * R1 * log(2.d0*D0/D1)

      func_params(1) = D0/scale_d
      func_params(2) = R1
      func_params(3) = H1
      func_params(5) = D1/scale_d
      func_params(6) = GM_code
      func_params(7) = R200
      func_params(8) = c
      func_params(9) = Rout
      func_params(10) = alpha
      func_params(11) = scale_l
      func_params(12) = scale_t

      Rmin_grav = gravity_params(6)
      Rmax_grav = gravity_params(7)

      gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

      if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

      ! Note, self-grav off
    r_cyl=sqrt(xc**2+yc**2)
    r_sph=sqrt(r_cyl**2+zc**2)

    ! Need to integrate from P0 at (R0,0) to P at (r_cyl,zc)
    dR = dx_max/dR_scale ! so dR is positive

    if ( r_cyl .gt. Rout ) then ! mid plane
      rho = D1/scale_d*(r_cyl/Rout)**alpha
      rho_l = D1/scale_d*((r_cyl-0.5d0*dR)/Rout)**alpha
      rho_r = D1/scale_d*((r_cyl+0.5d0*dR)/Rout)**alpha
      pres = P0 * sech(Rout/R1) * (r_cyl/Rout)**beta
      pres_l = P0 * sech(Rout/R1) * ((r_cyl-dR/2.d0)/Rout)**beta
      pres_r = P0 * sech(Rout/R1) * ((r_cyl+dR/2.d0)/Rout)**beta
    else
      rho = max(D0 * sech(r_cyl/R1)/scale_d,D1/scale_d)
      pres = max(P0 * sech(r_cyl/R1), P1)
      temp = pres/rho ! --> I assume constant midplane temperature though

      rho_l = max(D0 * sech(max(r_cyl-dr/2.d0,0.d0)/R1)/scale_d,D1/scale_d)
      pres_l = max(P0 * sech(max(r_cyl-dr/2.d0,0.d0)/R1),P1)
      temp_l = pres_l/rho_l

      rho_r = max(D0 * sech((r_cyl+dr/2.d0)/R1)/scale_d,D1/scale_d)
      pres_r = max(P0 * sech((r_cyl+dr/2.d0)/R1),P1)
      temp_r = pres_r/rho_r
    endif

    ! z direction
    if ( sqrt( (r_cyl/Rout)**2 + (zc/Hout)**2 ) .gt. 1.0d0) then
      sech_term = (D1/scale_d)/rho * (r_sph/Rout)**alpha
    else
      sech_term = sech(zc/H1)
    endif

    ! Recalculate full density
    rho = rho * sech_term
    rho = max(rho,D1/scale_d*(r_sph/Rout)**alpha)
    func_params(9) = r_cyl
    call Romberg(0.d0,abs(zc),eps,z_loop,rom_int,dPdz_R_gal,func_params)
    if (rom_int .gt. pmin-pres) then
      pres = pres + rom_int
    else ! set to midplane temp
!          pres = temp*rho
      pres = pmin
    endif

    ! Recalculate full density_l
    r_sph_l = sqrt( (r_cyl-0.5d0*dR)**2 + zc**2)
    if ( sqrt( ((r_cyl-0.5d0*dR)/Rout)**2 + (zc/Hout)**2 ) .gt. 1.0d0) then
      sech_term = (D1/scale_d)/rho * (r_sph_l/Rout)**alpha
    else
      sech_term = sech(zc/H1)
    endif
    rho_l = rho_l * sech_term
    rho_l = max(rho_l,D1/scale_d*(r_sph_l/Rout)**alpha)
    func_params(9) = max(r_cyl-dr/2.d0,Rin)
    call Romberg(0.d0,abs(zc),eps,z_loop,rom_int,dPdz_R_gal,func_params)
    if (rom_int .gt. pmin-pres_l) then
      pres_l =  pres_l + rom_int
    else ! set to midplane temp
!          pres_l = temp_l*rho_l
      pres_l = pmin
    endif

    ! Recalculate full density_r
    r_sph_r = sqrt( (r_cyl+0.5d0*dR)**2 + zc**2)
    if ( sqrt( ((r_cyl+0.5d0*dR)/Rout)**2 + (zc/Hout)**2 ) .gt. 1.0d0) then
      sech_term = (D1/scale_d)/rho * (r_sph_r/Rout)**alpha
    else
      sech_term = sech(zc/H1)
    endif
    rho_r = rho_r * sech_term
    rho_r = max(rho_r,D1/scale_d*(r_sph_r/Rout)**alpha)
    func_params(9) = r_cyl+dr/2.d0
    call Romberg(0.d0,abs(zc),eps,z_loop,rom_int,dPdz_R_gal,func_params)
    if (rom_int .gt. pmin-pres_r) then
      pres_r = pres_r + rom_int
    else ! set to midplane temp
!          pres_r = temp_r*rho_r
      pres_r = pmin
    endif

!        pres   =min(pres  ,disk_Textra*temp  *rho  )
!        pres_l =min(pres_l,disk_Textra*temp_l*rho_l)
!        pres_r =min(pres_r,disk_Textra*temp_r*rho_r)

    dPdR = (pres_r - pres_l) / ( r_cyl+dr/2.d0 - max(r_cyl-dR/2.d0,0.d0) )
    ! v_cyl**2 / r_cyl = -fgrav_R + dPdR/rho
    v_cyl = sqrt(max( -fgrav_R(zc,r_cyl,r_sph,func_params)*r_cyl + r_cyl*dPdR/rho,0.d0))

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    if (r_cyl .gt. 0.d0) then
      vel(1)= -v_cyl*yc/r_cyl + v_turb(1)
      if (ndim .gt. 1) then
        vel(2)=  v_cyl*xc/r_cyl + v_turb(2)
        if (ndim .gt. 2) then
          vel(3)= v_turb(3)
        endif
      endif
    else
      vel(1)= v_turb(1)
      vel(2)= v_turb(2)
      vel(3)= v_turb(3)
    endif
  else
    stop "Unknown disk_setup"
  endif

end subroutine relax_condinit
!================================================================
!================================================================
!================================================================
