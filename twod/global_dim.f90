module global_dim
  use global_com,only:dp
  implicit none
  save
  real(kind=dp),dimension(:,:),allocatable::pmatrix
  complex(kind=dp),dimension(:),allocatable::rj
  complex(kind=dp),dimension(:),allocatable::zpast
  complex(kind=dp),dimension(:,:,:),allocatable::prcdin
  integer,dimension(:,:),allocatable::IPIV
  real(kind=dp),dimension(:),allocatable::frequencies
end module global_dim
