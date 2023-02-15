!! comes from frig version and then all the specific development done by Valeska
!! to take into account extinction have been moved there
!! PH 19/01/2017
!=======================================================================
subroutine solve_cooling_parametric(nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  ! BRIDGE FUNCTION WITH SAME INTERFACE AS 
  ! Input/output variables to this function
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! zsolar - Metallicity in solar units (Zphys / 0.02)
  ! boost - raditation boost - exp(-100*nH) if self_shielding=True
  ! dt - cooling timestep in seconds
  ! deltaT2 - temperature change in K/mu (??)
  ! ncell - number of elements in the vector
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost
  ! Input/output variables to analytic function calc_temp 
  real(kind=8)::NN,TT, dt_tot_unicode
  ! Temporary variables
  integer::i
  real(kind=8)::TT_ini, mu
  ! Units
  real(kind=8) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
  mu = 1.4
  scale_T2 = scale_T2 * mu
  do i=1,ncell
     NN = nH(i) ! NOTE!! THE CODE BELOW ASSUMES scale_nH=1 !!
                ! SO WE LEAVE THIS AS IT IS TO KEEP UNITS CONSISTENCY
     TT = T2(i) / scale_T2
     TT_ini = TT
     dt_tot_unicode = dt / scale_t
     call calc_by_parametric_func(NN,TT,dt_tot_unicode)
     deltaT2(i) = (TT - TT_ini) * scale_T2
  end do
end subroutine solve_cooling_parametric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  calc_by_parametric_func(NN,TT,dt_tot_unicode)
    use amr_parameters
    use hydro_commons

    implicit none

    integer :: n,i,j,k,idim, iter, itermax,ii

    real(dp) :: dt, dt_tot, temps, dt_max, itermoy
    real(dp) :: rho,temp,dt_tot_unicode

    !alpha replaced by alpha_ct because of conflict with another alpha by PH 19/01/2017
    real(dp) :: mm,uma, kb, alpha_ct,mu,kb_mm
    real(dp) :: NN,TT, TTold, ref,ref2,dRefdT, eps, vardt,varrel, dTemp
    real(dp) :: rhoutot2
    real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
    ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
    mu = 1.4
    !
    ! Cette routine fonctionne en cgs
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    kb  =  1.38062d-16   ! erg/degre
    !  uma =  1.660531e-24  ! gramme
    !  mu  =  1.4
    !  mm = mu*uma
    !  kb_mm = kb / mm
    !  TT = TT  / kb  !/ kb_mm

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    if( TT .le. 0.) then
        TT = 50. / scale_T2
        return
    endif

    !if( TT*scale_T2 .gt. 50.) then
    !TT = 50. / scale_T2
    !return
    !endif

    vardt = 10.**(1./10.); varrel = 0.2

    dt_tot = dt_tot_unicode * scale_t ! * 3.08d18 / sqrt(kb_mm)
    TT     = TT * scale_T2

    !  nn = (rho/(gramme/cm3)) /mm

    itermax = 0 ; itermoy = 0.



    if (NN .le. smallr) then
        if( NN .le. 0)  write(*,*) 'prob dens',NN
        NN = smallr  !max(NN,smallr)
    endif


    ! alpha_ct = NN*kb_mm/(gamma-1.)
    alpha_ct = NN*kb/(gamma-1.)

    ! eps - a small offset of T to find gradient in T
    eps = 1d-5
    ! TODO : Parametric cooling function
    TT = (-lambda_0 / alpha_ct * NN ** power_m * dt_tot * (1-power_n) + TT ** (1-power_n)) ** (1/(1-power_n))

   !  iter  = 0 ; temps = 0.
   !  do while ( temps < dt_tot)
   !      if (TT .lt.0) then
   !          write(*,*) 'prob Temp',TT, NN
   !          !         write(*,*) 'repair assuming isobariticity'
   !          NN = max(NN,smallr)
   !          TT = min(4000./NN,8000.)  !2.*4000. / NN
   !      endif


   !      TTold = TT

   !      ! Calculate cooling rate
   !      !NN is assumed to be in cc and TT in Kelvin
   !      if (TT < 10035.d0) then
   !          call cooling_low(TT,NN,ref)
   !          call cooling_low(TT*(1d0+eps),NN,ref2)
   !      else
   !          call cooling_high(TT,NN,ref)
   !          call cooling_high(TT*(1d0+eps),NN,ref2)
   !      end if
        
   !      ! dT = T*(1+eps)-T = eps*T
   !      dRefdT = (ref2-ref)/(TT*eps)

   !      ! TODO STG - COPY THIS FUNCTION UP TO HERE, USE ref, drefdT TO 
   !      !            REPLACE rt_cmp_metals SOMEHOW


   !      !       write(*,*) 'check',TTold, TT,NN,ref,dRefdT,iter


   !      if (iter == 0) then
   !          if (dRefDT .ne. 0.) then
   !              dt = abs(1.0E-1 * alpha_ct/dRefDT)
   !          else
   !              dt = 1.0E-1 * dt_tot
   !          endif
   !          dt_max = dt_tot - temps
   !          if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
   !      endif

   !      dTemp = ref/(alpha_ct/dt - dRefdT)

   !      eps = abs(dTemp/TT)
   !      if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

   !      TT = TTold + dTemp
   !      if (TT < 0.) then
   !          write(*,*) 'Temperature negative !!!'
   !          write(*,*) 'TTold,TT   = ',TTold,TT
   !          write(*,*) 'rho   = ',rho
   !          TT = 100.  !*kelvin
   !      endif


   !      iter = iter + 1

   !      temps = temps + dt

   !      dt = vardt*varrel*dt/Max(vardt*eps, varrel)

   !      dt_max = dt_tot - temps
   !      if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
   !      !        write(*,987) temps, TT
   !      !987     format(E10.3,2x,E10.3)
   !      !        read (*,*)
   !  enddo


    !  if (TT .ge. 50.)  TT=50.

    !!now convert temperature in code units
    TT = TT / scale_T2

    return
end subroutine calc_by_parametric_func
