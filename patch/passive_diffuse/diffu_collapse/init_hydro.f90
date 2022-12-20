subroutine init_hydro
  use amr_commons
  use hydro_commons
  use radiation_parameters
  use rt_parameters
  use rt_hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::dummy_io,info,info2
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
#if NENER>0
  integer::irad
#endif
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,idim
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  real(dp)::d,u,v,w,A,B,C,e
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1108
  real(dp)::sum_dust
  integer::nvar_expected = nvar + 4 
#ifdef RT
  integer::nvar_rt = 0 ! Number of rt variables to read from hydro files
#endif
#if NDUST>0
  integer::idust
#endif
  integer::irecorder=0
  if(verbose)write(*,*)'Entering init_hydro'

  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar+3))
  allocate(unew(1:ncell,1:nvar+3))
  uold=0.0d0; unew=0.0d0
  if(fld)then
     allocate(rad_flux(1:ncell,1:nvar_bicg))
     allocate(urad(1:ncell,1:nvar_bicg))
     allocate(frad(1:ncell,1:ndim))
     allocate(in_sink(1:ncell))
     rad_flux=0.0d0; urad=0.0d0; frad=0.0d0; in_sink=.false.
  endif
  if(momentum_feedback)then
     allocate(pstarold(1:ncell))
     allocate(pstarnew(1:ncell))
     pstarold=0.0d0; pstarnew=0.0d0
  endif
  
#if NDUST>0  
  allocate(v_dust(1:ncell,1:ndust,1:ndim))
  v_dust=0.0d0
#endif  
#if NIMHD==1
  if(pressure_fix .or. nambipolar2.eq.1 .or.nmagdiffu2.eq.1)then
#else
  if(pressure_fix)then
#endif     
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  ! Variables for BICG scheme
  ! 1 : r
  ! 2 : p
  ! 3 : r*
  ! 4 : M-1
  ! 5 : 
  ! 6 : z and Ap
  ! 7 : p*
  ! 8 : p*A
  ! 9 : z*
  allocate(kappaR_bicg(1:ncell,1:ngrp))
  ! if FLD: matrix of size ngrpxngrp (because matrix only on Eg)
  ! if  M1: matrix of size (1+nrad)x(1+nrad) (on T,Eg,Fg)
  allocate(var_bicg(1:ncell,1:nvar_bicg,1:10+2*ndim))
  allocate(precond_bicg(1:ncell,1:nvar_bicg,1:nvar_bicg))
  if(store_matrix) then
     allocate(mat_residual_glob(1:ncell,1:nvar_bicg,1:nvar_bicg),residual_glob(1:ncell,1:nvar_bicg))
     allocate(coeff_glob_left(1:ncell,1:nvar_bicg,1:nvar_bicg,1:ndim),coeff_glob_right(1:ncell,1:nvar_bicg,1:nvar_bicg,1:ndim))
  else
     allocate(mat_residual_glob(1,1:nvar_bicg,1:nvar_bicg),residual_glob(1,1:nvar_bicg))
     allocate(coeff_glob_left(1,1:nvar_bicg,1:nvar_bicg,1:ndim),coeff_glob_right(1,1:nvar_bicg,1:nvar_bicg,1:ndim))
  endif
  kappar_bicg=0.0d0;var_bicg=0.0d0;precond_bicg=0.0d0
  mat_residual_glob=0.0d0;residual_glob=0.0d0
  coeff_glob_left=0.0d0;coeff_glob_right=0.0d0
  
  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     endif
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
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
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
!      if( (eos .and. nvar2.ne.nvar+3+1) .or. (.not.eos .and. nvar2.ne.nvar+3) )then
!     if(nvar2.ne.nvar+4)then
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar 
!     if(.not.(neq_chem.or.rt) .and. nvar2 .ne. nvar_expected)then
     if(.not.(neq_chem.or.rt.or.seed_pscal.or.seed_high_T) .and. nvar2.ne.nvar_expected)then
!----------------------------------------------!
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar_expected
        call clean_stop
     end if
#ifdef RT

     if(output_rtvar_in_hydro) then
        nvar_rt = NGroups * (ndim + 1)
        nvar_expected = nvar + 4 + nvar_rt
     end if

     if((neq_chem.or.rt).and.nvar2.lt.nvar_expected)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar_expected
        if(myid==1) write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to zero'
     end if
     if((neq_chem.or.rt).and.nvar2.gt.nvar_expected)then ! Not OK to drop variables
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found   =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar_expected
        call clean_stop
     end if
#endif
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
     if(seed_pscal.or.seed_high_T)then
       if(myid==1) then
        write(*,*)'Found nvar2  =',nvar2
        write(*,*)'Expected=',nvar_expected
        write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to zero'
        write(*,*)'Ignore old passive scalars from hydro.tmp'
        write(*,*)'Seed new passive scalars'
       end if
     end if
!-----------------------------------------------
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
              irecorder=0
                 iskip=ncoarse+(ind-1)*ngridmax
                 ! Loop over conservative variables
                 do ivar=1,4
                    read(ilun)xx
                    if(ivar==1)then ! Read density
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else
                       if(write_conservative) then ! Read momentum field
                          do i=1,ncache
                             uold(ind_grid(i)+iskip,ivar)=xx(i)
                          end do
                       else ! Read velocity field
                          do i=1,ncache
                             uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                          end do
                       endif
                    end if
                 end do
                 do ivar=6,8 ! Read left B field
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 do ivar=nvar+1,nvar+3 ! Read right B field
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
#if NENER>NGRP
                 if(write_conservative) then
                    ! Read non-thermal energies
                    do ivar=9,8+nent
                       read(ilun)xx
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)
                       end do
                    end do
                 else
                    ! Read non-thermal pressures --> non-thermal energies
                    do ivar=9,8+nent
                       read(ilun)xx
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)/(gamma_rad(ivar-8)-1.0d0)
                       end do
                    end do
                 endif
#endif

                 if(write_conservative) then
                    read(ilun)xx ! Read total energy
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,5)=xx(i)
                    enddo
                 else
                    read(ilun)xx ! Read pressure
                    if(.not.eos) then
                       do i=1,ncache
                          e=xx(i)/(gamma-1d0)
                          d=max(uold(ind_grid(i)+iskip,1),smallr)
                          u=uold(ind_grid(i)+iskip,2)/d
                          v=uold(ind_grid(i)+iskip,3)/d
                          w=uold(ind_grid(i)+iskip,4)/d
                          A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                          B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                          C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                          uold(ind_grid(i)+iskip,5)=e+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
                       end do
                    endif
                 endif

#if USE_FLD==1
                 do ivar=1,ngrp
                    read(ilun)xx ! Read radiative energy if any
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,firstindex_er+ivar) = xx(i)
                    end do
                 end do
#endif
#if USE_M_1==1
                 do ivar=1,nfr
                    read(ilun)xx ! Read radiative flux if any
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,firstindex_fr+ivar) = xx(i)
                    end do
                 end do
#endif

#if NEXTINCT>0
                 !Read extinction parameter
                 do ivar=1,nextinct
                    read(ilun)xx ! Read extinction if activated
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,firstindex_extinct+ivar) = xx(i)
                    end do
                 end do
#endif

#if NPSCAL>0
#if NIMHD==1
                 if(write_conservative) then
#ifdef RT
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4-nvar_rt)-4 ! Read conservative passive scalars if any
#else
                    !do ivar=1,npscal-4 ! Read conservative passive scalars if any
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4)-4 ! Read conservative passive scalars if any
#endif
                       read(ilun)xx
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
! discard old passive scalars, write only if no passive scalar seeding)
                       if(.not.(seed_pscal.or.(seed_high_T.and.T_seed(ivar-firstindex_pscal).gt.0)))then
!-----------------------------------------------
                       do i=1,ncache
                          !uold(ind_grid(i)+iskip,firstindex_pscal+ivar)=xx(i)
                          uold(ind_grid(i)+iskip,ivar)=xx(i)
                       end do
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
                       end if
!-----------------------------------------------
                    end do
                 else
#ifdef RT
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4-nvar_rt)-4 ! Read passive scalars if any
#else
                    !do ivar=1,npscal-4 ! Read passive scalars if any
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4)-4 ! Read passive scalars if any
#endif
                       read(ilun)xx
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar (discard old passive scalars)
                       if(.not.(seed_pscal.or.(seed_high_T.and.T_seed(ivar-firstindex_pscal).gt.0)))then
!-----------------------------------------------
                       do i=1,ncache
                          !uold(ind_grid(i)+iskip,firstindex_pscal+ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
                       end if
!-----------------------------------------------
                    end do
                 endif
#ifdef RT
                 do ivar=min(nvar,nvar2-4)-3,min(nvar,nvar2-4-nvar_rt)-1 ! Read current
#else
                 !do ivar=npscal-3,npscal-1 ! Read current
                 do ivar=min(nvar,nvar2-4)-3,min(nvar,nvar2-4)-1 ! Read current
#endif
                    read(ilun)xx
                    do i=1,ncache
                       !uold(ind_grid(i)+iskip,firstindex_pscal+ivar)=xx(i)
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do                 
#else
                 if(write_conservative) then
#ifdef RT
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4-nvar_rt)-1 ! Read conservative passive scalars if any
#else
                    !do ivar=1,npscal-1 ! Read conservative passive scalars if any
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4)-1 ! Read conservative passive scalars if any
#endif
                       read(ilun)xx
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar (discard old passive scalars)
                       if(.not.(seed_pscal.or.(seed_high_T.and.T_seed(ivar-firstindex_pscal).gt.0)))then
!-----------------------------------------------
                       do i=1,ncache
                          !uold(ind_grid(i)+iskip,firstindex_pscal+ivar)=xx(i)
                          uold(ind_grid(i)+iskip,ivar)=xx(i)
                       end do
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
                       end if
!-----------------------------------------------
                    end do
                 else
#ifdef RT
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4-nvar_rt)-1 ! Read passive scalars if any
#else
                    !do ivar=1,npscal-1 ! Read passive scalars if any
                    do ivar=firstindex_pscal+1,min(nvar,nvar2-4)-1 ! Read passive scalars if any
#endif
                       read(ilun)xx
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar (discard old passive scalars)
                       if(.not.(seed_pscal.or.(seed_high_T.and.T_seed(ivar-firstindex_pscal).gt.0)))then
!-----------------------------------------------
                       do i=1,ncache
                          !uold(ind_grid(i)+iskip,firstindex_pscal+ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
!-----------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
                       end if
!-----------------------------------------------
                    end do
                 endif
#endif

                 ! Read internal energy
                 read(ilun)xx
                 do i=1,ncache
                    uold(ind_grid(i)+iskip,firstindex_pscal+npscal)=xx(i)
                 end do

#endif

                 ! Read in the temperature
                 read(ilun)xx
!-------------------------------------------------
!! 04/2022 added by ynl for seeded high temperature recorder
!#if NPSCAL>0
!                 if(seed_high_T)then
!#if NIMHD==1
!                   do ivar=1,npscal-4
!#else
!                   do ivar=1,npscal-1
!#endif
!                      if(T_seed(ivar)>0.)then
!                         do i=1,ncache
!                            if(xx(i)>T_seed(ivar))then
!                               uold(ind_grid(i)+iskip,firstindex_pscal+ivar) = uold(ind_grid(i)+iskip,1)
!                            end if
!                         end do
!                      end if
!                   end do
!                 end if
!#endif
!-------------------------------------------------
                 if(.not.write_conservative) then
                    if(eos) then
                       !if eos, update the total energy
                       do i=1,ncache
                          d=max(uold(ind_grid(i)+iskip,1),smallr)
                          if(energy_fix) then
                             e=uold(ind_grid(i)+iskip,nvar)
                          else
                             sum_dust= 0.0d0
#if NDUST>0
                             do idust=1,ndust
                                sum_dust= sum_dust+ uold(ind_grid(i)+iskip,firstindex_ndust+idust)/d
                             enddo

#endif
                             call enerint_eos((1.0d0-sum_dust)*d,xx(i),e)                             
                          endif
                          u=uold(ind_grid(i)+iskip,2)/d
                          v=uold(ind_grid(i)+iskip,3)/d
                          w=uold(ind_grid(i)+iskip,4)/d
                          A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                          B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                          C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                          uold(ind_grid(i)+iskip,5)=e+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
                       end do
                    endif

#if NENER>0
                    do i=1,ncache
                       do irad=1,nener
                          uold(ind_grid(i)+iskip,5)=uold(ind_grid(i)+iskip,5)+uold(ind_grid(i)+iskip,8+irad)
                       end do
                    end do
#endif
                 endif

#ifdef RT
                 if(output_rtvar_in_hydro) then
                    ! Read-only
                    do ivar=1,nGroups
                       read(ilun)xx
                       do idim=1,ndim
                          read(ilun)xx
                       enddo
                    end do
                 end if
#endif
#if NDUST>0
                 !Read the dust velocity (needed for Kwok correction) and practical to have it in the outputs : TODO add outout v_dust
                 do idust=1,ndust
                    do idim=1,ndim   
                       read(ilun)xx        
                       do i=1,ncache
                          v_dust(ind_grid(i)+iskip,idust,idim)=xx(i)
                       end do
                    end do
                 end do
#endif
              end do
              deallocate(ind_grid,xx)
           end if
        end do
!-------------------------------------------------
! 04/2022 added by ynl for seeded passive scalar
     if(seed_pscal) call init_passive(ilevel)
!-------------------------------------------------
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

#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'

  end if

end subroutine init_hydro
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_passive(ilevel)
! adapted from mhd/init_flow_fine/init_flow_fine
  use amr_commons
  use hydro_commons
!  use units_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::i,igrid,ncache,iskip,ngrid
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar+3),save::uu

  if(numbtot(1,ilevel)==0)return
!  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  uu(1:nvector,1:nvar+3) = 0.0d0 
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
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid
  if(ncache.gt.0)then
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
          uu(i,1) = uold(ind_cell(i),1)
        end do
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              ! Gather cell centre positions
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              ! Rescale position from code units to user units
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))!*scale
           end do
        end do
             ! Call initial condition routine
        call seed_passive(xx,uu,dx,ngrid)
        ! Scatter variables
#if NIMHD==1
        do ivar=1,npscal-4
#else
        do ivar=1,npscal-1
#endif
           do i=1,ngrid
              uold(ind_cell(i),firstindex_pscal+ivar)=uu(i,firstindex_pscal+ivar)
           end do
        end do
      end do
    end do
  end if
  ! End loop over cells
111 format('   Entering init_passive for level ',I2)
end subroutine init_passive
!################################################################
!################################################################
!################################################################
!################################################################
subroutine seed_passive(x,u,dx,nn)
! adapted from mhd/condinit
! adapted from mhd/init_flow_fine/region_condinit
  use amr_commons
  use hydro_commons

  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u  ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position.
  integer::i,k
  real(dp)::vol,r,xn,yn,zn,en,rn,zn2
#if NVAR > NDIM + 2
  integer::ivar
#endif
  !================================================================
  ! This routine seeds passive scalars for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  !================================================================
  u(1:nvector,firstindex_pscal+1:nvar+3) = 0.0d0
  do k=1,nseeds
     ! Volume elements
     vol=dx**ndim
     ! For "point" regions only:
     if(seed_type(k) .eq. 0)then ! 0 for 'point' seeds
        do i=1,nn
!            xn=1.0; yn=1.0; zn=1.0
            xn=max(1.0-2.0*abs(x(i,1)-x_seed(k,1))/dx,0.0_dp)
            yn=max(1.0-2.0*abs(x(i,2)-x_seed(k,2))/dx,0.0_dp)
            zn=max(1.0-2.0*abs(x(i,3)-x_seed(k,3))/dx,0.0_dp)
            if(xn.gt.0..and.yn.gt.0..and.zn.gt.0.) then
            ! If point lies within cell, set pscal to density
!                write(*,*)"seeded no", k, xn, yn, zn, log(dx)/log(0.5D0)
                write(*,*)"seeded no", k, x(i,1),x(i,2),x(i,3), log(dx)/log(0.5D0)
                u(i,firstindex_pscal+k)=u(i,1)
            end if
        end do
     end if
!     ! For "ring" regions only:
!     if(seed_type(k) .eq. 1)then ! 1 for 'ring' seeds
!        do i=1,nn
!            r =  (ax_seed(k,1)**2+ax_seed(k,2)**2+ax_seed(k,3)**2)**(0.5)
!            ax_seed(k,:) = ax_seed(k,:)/r
!            xn=x(i,1)-x_seed(k,1)
!            yn=x(i,2)-x_seed(k,2)
!            zn=x(i,3)-x_seed(k,3)
!            r = (xn**2+yn**2+zn**2)**(0.5)
!            en = xn*ax_seed(k,1) + yn*ax_seed(k,2) + zn*ax_seed(k,3)
!            rn = (r**2-en**2)**(0.5)
!            xn = max(1.0-2.0*abs(xn*(1.0_dp-r_seed(k)/rn))/dx,0.0_dp)
!            yn = max(1.0-2.0*abs(yn*(1.0_dp-r_seed(k)/rn))/dx,0.0_dp)
!            zn = max(1.0-2.0*abs(zn*(1.0_dp-r_seed(k)/rn))/dx,0.0_dp)
!            if(xn.gt.0..and.yn.gt.0..and.zn.gt.0.) then
!            ! If point lies within cell, set pscal to density
!                u(i,firstindex_pscal+k)=u(i,1)
!            end if
!        end do
!     end if
!     ! For "square" regions only:
!     if(seed_type(k) .eq. 2)then ! 2 for 'square' seeds
!        ! Exponent of choosen norm
!        en=2
!        do i=1,nn
!           ! Compute position in normalized coordinates
!           xn=0.0d0; yn=0.0d0; zn=0.0d0
!           xn=2.0d0*abs(x(i,1)-x_seed(k,1))/r_seed(k)
!           yn=2.0d0*abs(x(i,2)-x_seed(k,2))/r_seed(k)
!           zn=2.0d0*abs(x(i,3)-x_seed(k,3))/r_seed(k)
!           ! Compute cell "radius" relative to region center
!           if(exp_region(k)<10)then
!              r=(xn**en+yn**en+zn**en)**(1.0/en)
!           else
!              r=max(xn,yn,zn)
!           end if
!           ! If cell lies within region,
!           ! Set pscal to density
!           if(r<1.0)then
!               u(i,firstindex_pscal+k)=u(i,1)
!           end if
!        end do
!     end if

  end do
  return
end subroutine seed_passive
!################################################################
!################################################################
!################################################################
!################################################################
