program cap2d_ecm
  use global_com 
  use global_geom,only:nglunk,nsuinf,ncond
  implicit none

  integer::idx,jdx
  real(kind=dp)::cap

  ! Input file storing all the geometry information
  open(unit=16,file='geo_mesh.inp',status='old')  
  ! Output file storing output info
  open(unit=17,file='info.out',status='unknown')  

  ! read geometry info for conductor and dielectric
  call ingeo

  !// LEI: Moved to interfaces.f90
  !// ! Total number of unknown
  !// nglunk=nsuinf(1)+nsuinf(2)
  !// print*,'nglunk',nglunk
  
  ! Solve for capacitance
  call ky_ecm_init(.true., 1.0d-4, 4)
  
  call ky_simulate
  do idx=0,ncond-1
     do jdx=0,ncond-1
        ! Get capacitance between conductor idx and jdx
        call ky_getC(idx,jdx,cap)
     end do
  end do
  ! deallocate array
  call ky_clear_edge
     
  print*,'Closing files-safely here.'
  close(16,status='keep')
  close(17,status='keep')
  stop
end program cap2d_ecm
