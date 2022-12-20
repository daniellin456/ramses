subroutine file_descriptor_hydro(filename)
  use amr_commons
  use hydro_commons
  use rt_hydro_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ivar,ilun,ivar_bef
#if NDUST>0
  integer :: idust
#endif  
  if(verbose)write(*,*)'Entering file_descriptor_hydro'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
#ifdef RT
  if(output_rtvar_in_hydro)then
     write(ilun,'("nvar        =",I11)')nvar+4+NGroups*(ndim+1)+ndust*ndim
  else
     write(ilun,'("nvar        =",I11)')nvar+4
  end if
#else
  write(ilun,'("nvar        =",I11)')nvar+4+ndust*ndim
#endif
  ivar=1
  write(ilun,'("variable #",I2,": density")')ivar
  if(write_conservative) then
     ivar=2
     write(ilun,'("variable #",I2,": momentum_x")')ivar
     ivar=3
     write(ilun,'("variable #",I2,": momentum_y")')ivar
     ivar=4
     write(ilun,'("variable #",I2,": momentum_z")')ivar
  else
     ivar=2
     write(ilun,'("variable #",I2,": velocity_x")')ivar
     ivar=3
     write(ilun,'("variable #",I2,": velocity_y")')ivar
     ivar=4
     write(ilun,'("variable #",I2,": velocity_z")')ivar
  endif
  ivar=5
  write(ilun,'("variable #",I2,": B_left_x")')ivar
  ivar=6
  write(ilun,'("variable #",I2,": B_left_y")')ivar
  ivar=7
  write(ilun,'("variable #",I2,": B_left_z")')ivar
  ivar=8
  write(ilun,'("variable #",I2,": B_right_x")')ivar
  ivar=9
  write(ilun,'("variable #",I2,": B_right_y")')ivar
  ivar=10
  write(ilun,'("variable #",I2,": B_right_z")')ivar
#if NENER>NGRP
  if(write_conservative) then
#if NCR>0  
     ! CR energies
     do ivar=1,ncr
        write(ilun,'("variable #",I2,": cosmic_rays_energy_",I1)')10+ivar,ivar
     end do
#endif
     ! Non-thermal energies
     do ivar=1+ncr,nent
        write(ilun,'("variable #",I2,": non_thermal_energy_",I1)')10+ivar,ivar
     end do
  else
#if NCR>0
     ! CR pressures
     do ivar=1,ncr
        write(ilun,'("variable #",I2,": cosmic_rays_pressure_",I1)')10+ivar,ivar
     end do
#endif
     ! Non-thermal pressures
     do ivar=1+ncr,nent
        write(ilun,'("variable #",I2,": non_thermal_pressure_",I1)')10+ivar,ivar
     end do
  end if
#endif
  if(write_conservative) then
     ivar=11+nent
     write(ilun,'("variable #",I2,": total_energy")')ivar
  else
     ivar=11+nent
     write(ilun,'("variable #",I2,": thermal_pressure")')ivar
  endif
#if NGRP>0
  ! Radiative energies
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_energy_",I1)')firstindex_er+3+ivar,ivar
  end do
#endif
#if USE_M_1==1
  ! Radiative fluxes
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_flux_x",I1)')firstindex_fr+3       +ivar,ivar
  end do
if(ndim>1) then
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_flux_y",I1)')firstindex_fr+3+  ngrp+ivar,ivar
  end do
endif
if(ndim>2) then
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_flux_z",I1)')firstindex_fr+3+2*ngrp+ivar,ivar
  end do
endif
#endif
#if NEXTINCT>0
  ! Extinction
  do ivar=1,nextinct
     write(ilun,'("variable #",I2,": extinction",I1)')firstindex_extinct+3+ivar,ivar
  end do
#endif
#if NPSCAL>0
#if NIMHD==1
  ! Passive scalars excluding current and internal energy
  if(write_conservative) then
#if NDUST>0
     ivar_bef=1
     if(firstindex_ndust-firstindex_pscal.ge.1) then
        do ivar=1,firstindex_ndust-firstindex_pscal
           if(ivar<10)write(ilun,'("variable #",I2,": passive_scalar_cons_",I1)')firstindex_pscal+3+ivar,ivar
           if(ivar.ge.10)write(ilun,'("variable #",I2,": passive_scalar_cons_",I2)')firstindex_pscal+3+ivar,ivar
        end do
        ivar= firstindex_ndust-firstindex_pscal
     end if
     do ivar=1,ndust
        if(ivar<10)write(ilun,'("variable #",I2,": density_dust_",I1)')firstindex_ndust+3+ivar,ivar
        if(ivar.ge.10)write(ilun,'("variable #",I2,": density_dust_",I2)')firstindex_ndust+3+ivar,ivar
     end do
     do ivar=1,ndust_pscal
        if(ivar<10)write(ilun,'("variable #",I2,": size_density_dust__",I1)')firstindex_dustpscal+3+ivar,ivar
        if(ivar.ge.10)write(ilun,'("variable #",I2,": size_density_dust_",I2)')firstindex_dustpscal+3+ivar,ivar
     end do
     do ivar=1, npscal-4-ndust-ndust_pscal
        if(ivar_bef<10)write(ilun,'("variable #",I2,": passive_scalar_cons_",I1)')firstindex_pscal+3+ivar_bef,ivar_bef
        if(ivar_bef.ge.10)write(ilun,'("variable #",I2,": passive_scalar_cons_",I2)')firstindex_pscal+3+ivar_bef,ivar_bef
        ivar_bef=ivar_bef+1
     end do
#else
     do ivar=1,npscal-4
        write(ilun,'("variable #",I2,": passive_scalar_cons_",I1)')firstindex_pscal+3+ivar,ivar
     end do
#endif       
  else
#if NDUST>0
     ivar_bef=1
     if(firstindex_ndust-firstindex_pscal.ge.1) then
        do ivar=1,firstindex_ndust-firstindex_pscal
           if(ivar<10)write(ilun,'("variable #",I2,": passive_scalar_",I1)')firstindex_pscal+3+ivar,ivar
           if(ivar.ge.10)write(ilun,'("variable #",I2,": passive_scalar_",I2)')firstindex_pscal+3+ivar,ivar
        end do
        ivar= firstindex_ndust-firstindex_pscal
     end if
     do ivar=1,ndust
        if(ivar<10)write(ilun,'("variable #",I2,": dust_ratio_",I1)')firstindex_ndust+3+ivar,ivar
        if(ivar.ge.10)write(ilun,'("variable #",I2,": dust_ratio_",I2)')firstindex_ndust+3+ivar,ivar
     end do
     do ivar=1,ndust_pscal
        if(ivar<10)write(ilun,'("variable #",I2,": dust_size_",I1)')firstindex_dustpscal+3+ivar,ivar
        if(ivar.ge.10)write(ilun,'("variable #",I2,": dust_size_",I2)')firstindex_dustpscal+3+ivar,ivar
     end do
     do ivar=1, npscal-4-ndust-ndust_pscal
        if(ivar_bef<10)write(ilun,'("variable #",I2,": passive_scalar_",I1)')firstindex_pscal+3+ivar_bef,ivar_bef
        if(ivar_bef.ge.10)write(ilun,'("variable #",I2,": passive_scalar_",I2)')firstindex_pscal+3+ivar_bef,ivar_bef
        ivar_bef=ivar_bef+1
     end do
#else
     do ivar=1,npscal-4
        write(ilun,'("variable #",I2,": passive_scalar_",I1)')firstindex_pscal+3+ivar,ivar
     end do
#endif 
  endif
  ivar=npscal-3
  write(ilun,'("variable #",I2,": current_x")')firstindex_pscal+3+ivar
  ivar=npscal-2
  write(ilun,'("variable #",I2,": current_y")')firstindex_pscal+3+ivar
  ivar=npscal-1
  write(ilun,'("variable #",I2,": current_z")')firstindex_pscal+3+ivar
#else
  ! Passive scalars excluding internal energy
  if(write_conservative) then
#if NDUST>0
     do ivar=1,npscal-1
        if(firstindex_pscal+ivar.le.firstindex_ndust .or. firstindex_pscal+ivar > firstindex_ndust+ndust ) then
           if(ivar<10)write(ilun,'("variable #",I2,": passive_scalar_cons_",I1)')firstindex_pscal+3+ivar,ivar
           if(ivar.ge.10)write(ilun,'("variable #",I2,": passive_scalar_cons_",I2)')firstindex_pscal+3+ivar,ivar
        else if(firstindex_pscal+ivar>firstindex_ndust .and. firstindex_pscal+ivar .le. firstindex_ndust+ndust) then
           idust=firstindex_pscal+ivar-firstindex_ndust
           if(idust<10)write(ilun,'("variable #",I2,": density_dust_",I1)')firstindex_pscal+3+ivar,idust
           if(idust.ge.10)write(ilun,'("variable #",I2,": density_dust_",I2)')firstindex_pscal+3+ivar,idust
        endif
     end do
#else
     do ivar=1,npscal-1
        write(ilun,'("variable #",I2,": passive_scalar_cons_",I1)')firstindex_pscal+3+ivar,ivar
     end do
#endif          
  else
#if NDUST>0
     do ivar=1,npscal-1
        if(firstindex_pscal+ivar .le. firstindex_ndust .or. firstindex_pscal+ivar > firstindex_ndust+ndust ) then
           if(ivar-npscal+ndust<10)write(ilun,'("variable #",I2,": passive_scalar_",I1)')firstindex_pscal+3+ivar,ivar-npscal+ndust
           if(ivar-npscal+ndust.ge.10)write(ilun,'("variable #",I2,": passive_scalar_",I2)')firstindex_pscal+3+ivar,ivar-npscal+ndust
        else if(firstindex_pscal+ivar.ge. firstindex_ndust .and. firstindex_pscal+ivar .le. firstindex_ndust+ndust) then
           idust=firstindex_pscal+ivar-firstindex_ndust

           if(idust<10) write(ilun,'("variable #",I2,": dust_ratio_",I1)')firstindex_pscal+3+ivar,idust
           if(idust.ge.10) write(ilun,'("variable #",I2,": dust_ratio_",I2)')firstindex_pscal+3+ivar,idust
        endif
     end do
#else
     do ivar=1,npscal-1
        write(ilun,'("variable #",I2,": passive_scalar_",I1)')firstindex_pscal+3+ivar,ivar
     end do
#endif       
  endif
#endif
  ivar=npscal
  write(ilun,'("variable #",I2,": internal_energy")')firstindex_pscal+3+ivar
#endif
  ! Temperature
  ivar=firstindex_pscal+3+npscal+1
  write(ilun,'("variable #",I2,": temperature")')ivar
#if NDUST>0
  do idust=1,ndust
     if(idust<10)write(ilun,'("variable #",I2,": velocity_drift_",I1,"_x")')ivar+idust,idust
     if(idust.ge.10)write(ilun,'("variable #",I2,": velocity_drift_",I2,"_x")')ivar+idust,idust

#if NDIM>1
     if(idust<10)write(ilun,'("variable #",I2,": velocity_drift_",I1,"_y")')ivar+idust+1,idust
     if(idust.ge.10)write(ilun,'("variable #",I2,": velocity_drift_",I2,"_y")')ivar+idust+1,idust

#endif
#if NDIM>2
     if(idust<10)write(ilun,'("variable #",I2,": velocity_drift_",I1,"_z")')ivar+idust+2,idust
     if(idust.ge.10)write(ilun,'("variable #",I2,": velocity_drift_",I2,"_z")')ivar+idust+2,idust

#endif          
     ivar = ivar +(ndim-1)
  end do
#endif  
#ifdef RT

  if(output_rtvar_in_hydro)then
     do ivar=1,nGroups
        write(ilun,'("variable #",I2,": photon_number_flux_",I1)')nvar+4+ndust*ndim+ivar,ivar
        write(ilun,'("variable #",I2,": photon_number_xflux_",I1)')nvar+4+ndust*ndim+ivar+1,ivar
        write(ilun,'("variable #",I2,": photon_number_yflux_",I1)')nvar+4+ndust*ndim+ivar+2,ivar
        write(ilun,'("variable #",I2,": photon_number_zflux_",I1)')nvar+4+ndust*ndim+ivar+3,ivar
     end do
  end if
#endif

  close(ilun)

end subroutine file_descriptor_hydro

subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  use rt_hydro_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::dummy_io,info2
#endif
  character(LEN=80)::filename
  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,ht,idim
  real(dp)::d,u,v,w,A,B,C,e
  integer,allocatable,dimension(:)::ind_grid
  real(dp)::cmp_temp,p
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1121
#if NENER>0
  integer::irad
#endif

  real(dp) :: sum_dust
#if NDUST>0
  integer :: idust,idust_pscal
#endif  
  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
!   if(eos) then 
!      write(ilun)nvar+3+1
!   else
!      write(ilun)nvar+3
!   endif
#ifdef RT
  if(output_rtvar_in_hydro)then
     write(ilun)nvar+4+NGroups*(ndim+1)
  else
     write(ilun)nvar+4
  end if
#else
  write(ilun)nvar+4
#endif
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do ivar=1,4
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else ! Write velocity field
                    if(write_conservative) then
                       do i=1,ncache
                          xdp(i)=uold(ind_grid(i)+iskip,ivar)
                       end do
                    else
                       do i=1,ncache
                          xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 endif
                 write(ilun)xdp
              end do
              do ivar=6,8 ! Write left B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
              do ivar=nvar+1,nvar+3 ! Write right B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
#if NENER>NGRP
              ! Write non-thermal pressures
              if(write_conservative) then
                 do ivar=1,nent
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,8+ivar)
                    end do
                    write(ilun)xdp
                 end do
              else
                 do ivar=1,nent
                    do i=1,ncache
                       xdp(i)=(gamma_rad(ivar)-1d0)*uold(ind_grid(i)+iskip,8+ivar)
                    end do
                    write(ilun)xdp
                 end do
              endif
#endif
              if(write_conservative) then
                 do i=1,ncache ! Write total energy
                    xdp(i)=uold(ind_grid(i)+iskip,5)
                 enddo
                 write(ilun)xdp
              else
                 do i=1,ncache ! Write thermal pressure
                    d=max(uold(ind_grid(i)+iskip,1),smallr)
                    if(energy_fix) then
                       e=uold(ind_grid(i)+iskip,nvar)
                    else
                       u=uold(ind_grid(i)+iskip,2)/d
                       v=uold(ind_grid(i)+iskip,3)/d
                       w=uold(ind_grid(i)+iskip,4)/d
                       A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                       B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                       C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                       e=uold(ind_grid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                       do irad=1,nener
                          e=e-uold(ind_grid(i)+iskip,8+irad)
                       end do
#endif
                    endif
                    sum_dust=0.0d0
#if NDUST>0
                    do idust=1,ndust
                       sum_dust=sum_dust+uold(ind_grid(i)+iskip,firstindex_ndust+idust)/d
                    end do
#endif                            
                    call pressure_eos((1.0d0-sum_dust)*d,e,p)
                    xdp(i)=p
                 end do
                 write(ilun)xdp
              endif

#if NGRP>0
              do ivar=1,ngrp ! Write radiative energy if any
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,firstindex_er+ivar)
                 end do
                 write(ilun)xdp
              end do
#if USE_M_1==1
              do ivar=1,nfr ! Write radiative flux if any
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,firstindex_fr+ivar)
                 end do
                 write(ilun)xdp
              end do
#endif
#endif
#if NEXTINCT>0
              ! Write extinction if activated
              do i=1,ncache
                 xdp(i)=uold(ind_grid(i)+iskip,firstindex_extinct+1)
              end do
              write(ilun)xdp
#endif

#if NPSCAL>0
#if NIMHD==1
              if(write_conservative) then
                 do ivar=1,npscal-4 ! Write conservative passive scalars if any
                    do i=1,ncache                      
                       xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)                       
                    end do
                    write(ilun)xdp
                 end do
              else
                 do ivar=1,npscal-4 ! Write passive scalars if any
                    do i=1,ncache
                       if(firstindex_pscal+ivar.ge.firstindex_dustpscal+1.and.firstindex_pscal+ivar.le.firstindex_dustpscal+ndust)then ! /!\ only one dust pscal is allowed for now
                          xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)/uold(ind_grid(i)+iskip,firstindex_pscal+(ivar-ndust))
                       else  
                          xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                       endif 
                       
                    end do
                    write(ilun)xdp
                 end do
              endif

              do ivar=npscal-3,npscal-1 ! Write current
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)
                 end do
                 write(ilun)xdp
              end do

#else
              if(write_conservative) then
                 do ivar=1,npscal-1 ! Write conservative passive scalars if any
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)
                    end do
                    write(ilun)xdp
                 end do
              else
                 do ivar=1,npscal-1 ! Write passive scalars if any
                    do i=1,ncache                     
                       xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                    write(ilun)xdp
                 end do
              endif
#endif
              
              ! Write internal energy
              do i=1,ncache
                 xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+npscal)
              end do
              write(ilun)xdp
              
#endif
              
              ! Write temperature
              do i=1,ncache
                 d=max(uold(ind_grid(i)+iskip,1),smallr)
                 if(energy_fix) then
                    e=uold(ind_grid(i)+iskip,nvar)
                 else
                    u=uold(ind_grid(i)+iskip,2)/d
                    v=uold(ind_grid(i)+iskip,3)/d
                    w=uold(ind_grid(i)+iskip,4)/d
                    A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                    B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                    C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                    e=uold(ind_grid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                    do irad=1,nener
                       e=e-uold(ind_grid(i)+iskip,8+irad)
                    end do
#endif
                 endif
                    sum_dust=0.0d0
#if NDUST>0
                    do idust=1,ndust
                       sum_dust=sum_dust+uold(ind_grid(i)+iskip,firstindex_ndust+idust)/d
                    end do
#endif       
                 call temperature_eos((1.0d0-sum_dust)*d,e,cmp_temp,ht)                 
                 xdp(i)=cmp_temp
              end do
              write(ilun)xdp
#if NDUST>0
           do idust=1,ndust
              do idim=1,ndim           
                 do i=1,ncache
                    xdp(i)=v_dust(ind_grid(i)+iskip,idust,idim)
                 end do
                 write(ilun)xdp
              end do
           end do
#endif
           
#ifdef RT
              if(output_rtvar_in_hydro)then
                 do ivar=1,nGroups
                    do i=1,ncache
                       xdp(i)=rtuold(ind_grid(i)+iskip,iGroups(ivar))!*rt_c
                    end do
                    write(ilun)xdp
                    do idim=1,ndim
                       ! Store photon flux
                       do i=1,ncache
                          xdp(i)=rtuold(ind_grid(i)+iskip,iGroups(ivar)+idim)
                       end do
                       write(ilun)xdp
                    enddo
                 end do
              end if
#endif
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

end subroutine backup_hydro





