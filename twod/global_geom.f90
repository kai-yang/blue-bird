module global_geom 
  use global_com,only:dp
  implicit none
  save
  integer::nglunk,nsuunk,add_point_ptr, add_edge_ptr
! pec body parameters
  integer,dimension(:),allocatable::edg_cond,cond_sta,edg_dmg,dmg_sta
  real(kind=dp),allocatable::edg_dmg_epsr(:,:),cond_epsr(:)
  real(kind=dp),dimension(:,:,:),allocatable::edg_coord,edg_dmg_coord
  real(kind=dp),dimension(:,:),allocatable::sunod
  real(kind=dp),dimension(:,:),allocatable::unormal
  integer,allocatable::npat_cond(:)
  real(kind=dp)::edge_av,edge_av_dmg
  real(kind=dp),allocatable::rho_max(:,:),z_min(:),z_max(:)

  integer,dimension(3)::nsuinf

  integer,dimension(:,:),allocatable::nsupan,nsuedgn,nsupae
  integer,dimension(:,:),allocatable::nsuedn,nsuedp,basis_nodes

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
