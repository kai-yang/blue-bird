module global_geom 
  use global_com,only:dp
  implicit none
  save
  integer::nglunk,nsuunk,ncond,ncond_node,ncond_edg,ndiel_edg,add_cond_node_ptr,add_cond_edge_ptr,add_diel_edge_ptr
! Geometry array
  integer,dimension(2)::nsuinf
  integer,dimension(:),allocatable::edg_cond,cond_sta
  real(kind=dp),allocatable::edg_diel_epsr(:,:),cond_epsr(:)
  real(kind=dp),allocatable::node_coord(:,:)
  real(kind=dp),dimension(:,:,:),allocatable::edg_coord,edg_diel_coord
  real(kind=dp)::edge_av_cond,edge_av_diel,edge_av
  real(kind=dp)::scale
end module global_geom 
