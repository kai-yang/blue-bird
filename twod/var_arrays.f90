subroutine var_arrays
  use global_com,only:dp
  use global_dim,only:rj
  use global_geom,only:nglunk
  implicit none
  integer::alloc_err

  allocate(rj(1:nglunk), stat = alloc_err)
  if (alloc_err/=0) print*,'Memory fault in RJ'
  rj(1:nglunk)=cmplx(0.0_dp,0.0_dp,dp)

  return
end subroutine var_arrays


