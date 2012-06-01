module global_geom 
  use global_com,only:dp
  implicit none
  save
  integer::nglunk,nsuunk,add_point_ptr, add_edge_ptr
! pec body parameters
  real(kind=dp),dimension(:,:),allocatable::sunod
  real(kind=dp),dimension(:,:),allocatable::unormal
  integer,allocatable::npat_cond(:)
  real(kind=dp)::edge_av
  real(kind=dp),allocatable::rho_max(:,:),z_min(:),z_max(:)

  integer,dimension(3)::nsuinf
! nsuinf(3) is not equal to nsuunk!
  integer,dimension(:,:),allocatable::nsupan,nsuedgn,nsupae
  integer,dimension(:,:),allocatable::nsuedn,nsuedp,basis_nodes
! juction edge paramaters:
  integer::njunction_edges
  integer,dimension(:),allocatable::junction_store
  integer,dimension(:,:),allocatable::junction_patches
  integer::max_pae
! ******* WARNING: JUNCTION EDGES ARE INVISIBLE TO NSUPAE
!         the second patch of a junction edge does not list the junction edge
!         as one of its edges
  integer::disconnected_body_count
  type body_patch
     integer,pointer::id(:)
  end type body_patch
  type(body_patch)::patches_of_body(256) !at most 100 disconnected bodies

! Resistor/Impedance loading commons!
  integer::nnloaded
  integer,allocatable,dimension(:):: nload
  complex(kind=dp),allocatable,dimension(:)::rrload
  real(kind=dp)::estimated_edge_av
contains
  function sunod1(k)
    implicit none
    integer,intent(in)::k
    real(kind=dp)::sunod1(1:2)
!    close(unit=99,status='keep')
    sunod1(1:2)=sunod(1:2,k)
    return
  end function sunod1

  function unormal1(k)
    implicit none
    integer,intent(in)::k
    real(kind=dp)::unormal1(1:3)
!    read(unit=95,rec=k) unormal1(1:3)
    unormal1(1:3)=unormal(1:3,k)
    return
  end function unormal1
end module global_geom 
