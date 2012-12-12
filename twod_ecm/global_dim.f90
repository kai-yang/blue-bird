module global_dim
  use global_com,only:dp
  implicit none
  save
  real(kind=dp),dimension(:,:),allocatable::pmatrix
  real(kind=dp),dimension(:),allocatable::rj
  real(kind=dp),dimension(:),allocatable::zpast
  real(kind=dp),dimension(:),allocatable::prcdin
end module global_dim
