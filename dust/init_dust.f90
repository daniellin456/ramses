subroutine init_dust_ratio(dustratio,epsilondust)
  use amr_commons
  use hydro_commons
  use hydro_parameters, ONLY: ndust_max_mrn,smalleps
  implicit none
  
  real(dp) :: dustratio
  real(dp), dimension(1:ndust):: epsilondust
  real(dp), dimension(1:ndust+1):: sdust
  real(dp) :: epsilon_0,Anorm
  integer  :: idust

  
  !You can replace these lines to use a custom dust distribution
  !Default is MRN
  
  epsilon_0 = dustratio/(1.0d0+dustratio)
  do idust =1,ndust+1
     sdust(idust) = 10**(log10(size_max/size_min)*dble(idust-1)/dble(ndust)+log10(size_min))
  enddo
  Anorm=0.0d0
  do idust=1,ndust_max_mrn
     Anorm = Anorm + (sdust(idust+1)**(4.0d0-mrn_index)-sdust(idust)**(4.0d0-mrn_index))
  end do
  do idust=1,ndust
     epsilondust(idust)=smalleps
  end do
  do idust=1,ndust_max_mrn
     epsilondust(idust)= epsilon_0*(sdust(idust+1)**(4.0d0-mrn_index)-sdust(idust)**(4.0d0-mrn_index))/Anorm
  enddo
  if(ndust_max_mrn<ndust) then
     epsilondust(ndust_max_mrn+1)=(epsilondust(ndust_max_mrn)+smalleps)/2.0d0
   endif 
end subroutine init_dust_ratio

subroutine size_dust(sdust)
  use amr_commons
  use hydro_commons
  implicit none
  
  real(dp), dimension(1:ndust):: sdust
  real(dp), dimension(1:ndust+1):: sdust_interval
  integer  :: idust
  
  !You can replace these lines to use a custom dust bining, /!\ must be consistent with what you put in init_dust_ratio
  !Default is log bining
  
  do idust =1,ndust+1
     sdust_interval(idust) = 10**(log10(size_max/size_min)*dble(idust-1)/dble(ndust)+log10(size_min))
  enddo
  !We compute the average dust size in the bin to get the correct stopping time 
  do idust =1,ndust
     sdust(idust) = sqrt(sdust_interval(idust)*sdust_interval(idust+1))
  enddo

  !print*, sdust

  !print *, 4./3.*acos(-1.0d0)*sdust**3.
  
end subroutine size_dust


subroutine size_dust_interface(sdust)
  use amr_commons
  use hydro_commons
  implicit none
  
  real(dp), dimension(1:ndust+1):: sdust
  integer  :: idust
 
  
  do idust =1,ndust+1
     sdust(idust) = 10**(log10(size_max/size_min)*dble(idust-1)/dble(ndust)+log10(size_min))
  enddo

end subroutine size_dust_interface


