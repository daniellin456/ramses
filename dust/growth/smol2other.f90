!--------------------------------------
! MODULE: interface with other codes
! physical input parameters:
! ndusttypes     -> number of dust species
! rhogas         -> gas mass density
! rhodust        -> dust mass density array
! rhograin       -> intrinsic density of grains
! dv             -> 2D array of relative velocities
! grainssizegrid -> sizegrid linked to rhodust
! dthydro        -> hydrodynamic timestep
! timetosec      -> coefficient to convert unit time in second
!
! physical output parameter:
! rhodust        -> dust mass density array evolved after 1 hydro timestep
!
!
! algorithm parameters:
! k                -> order or polynomials to approximate rhodust by polynomial piecewise function
! grainsmassgrid   -> massgrid linked to rhodust
! grainsmassgridDL -> massgrid dimensionless linked to rhodust
! grainsmassbinsDL -> array of mass at middle of bins linked to rhodust
! tabflux          -> 6D array with dim (ndusttypes,k+1,k+1,k+1,ndusttypes,ndusttypes) to evaluate coagulation flux in DG scheme
! tabPhi           -> 3D array with dim (ndusttypes+1,k+1,2) for polynomials Legendre basis in DG scheme
! tabhc            -> 2D array with dim (ndusttypes,k+1) for normalisation coefficient for polynomials Legendre basis in DG scheme
! dvDL             -> 2D array with dim (ndusttypes,ndusttypes) of relative velocities dimensionless
!
!
! First run grain_growth_init to initialise arrays for DG scheme
! Then run grain_growth to evolve rhodust
!
! link between rhodust and dust mass distribution g in Smoluchowski equation (cgs)
! g unit is cm^(-3)
! then in each bin gi = rhodust_i/(dm_i) -> cm^(-3)
! with dm_i size of the bin i
!
!
!--------------------------------------
module smol2other


   implicit none

   interface
      module subroutine grain_growth_init(ndusttypes,ulength,umass,utime,dv,&
                                          grainsmassgrid,grainsmassbins,grainsmassgridDL,grainsmassbinsDL,&
                                          tabflux,tabhc,dvDL)
         use precision        

         implicit none
         integer,  intent(in)                                                :: ndusttypes
         real(wp), intent(in)                                                :: ulength,umass,utime
         real(wp), intent(in),  dimension(ndusttypes)                        :: grainsmassbins
         real(wp), intent(in),  dimension(ndusttypes+1)                      :: grainsmassgrid
         real(wp), intent(in),  dimension(ndusttypes,ndusttypes)             :: dv
         real(wp), intent(out), dimension(ndusttypes,ndusttypes,ndusttypes)  :: tabflux
         real(wp), intent(out), dimension(ndusttypes)                        :: tabhc
         real(wp), intent(out), dimension(ndusttypes,ndusttypes)             :: dvDL
         real(wp), intent(out), dimension(ndusttypes+1)                      :: grainsmassgridDL
         real(wp), intent(out), dimension(ndusttypes)                        :: grainsmassbinsDL

      end subroutine grain_growth_init
   end interface


   interface
      module subroutine grain_growth(ndusttypes,ulength,grainsmassgrid,grainsmassgridDL,grainsmassbinsDL,&
                                       tabhc,tabflux,rhodust,dthydro)
         use precision

         implicit none
         integer,  intent(in)                                                  :: ndusttypes
         real(wp), intent(in)                                                  :: dthydro,ulength
         real(wp), intent(in),    dimension(ndusttypes+1)                      :: grainsmassgrid,grainsmassgridDL
         real(wp), intent(in),    dimension(ndusttypes)                        :: grainsmassbinsDL
         real(wp), intent(in),    dimension(ndusttypes,ndusttypes,ndusttypes)  :: tabflux
         real(wp), intent(in),    dimension(ndusttypes)                        :: tabhc
         real(wp), intent(inout), dimension(ndusttypes)                        :: rhodust
      end subroutine grain_growth
   end interface


end module smol2other



