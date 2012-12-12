subroutine ingeo
  ! Read geometry info
  use global_com,only:dp
  use global_geom,only:nsuinf,ncond,ncond_node,ncond_edg,node_coord,&
       ndiel_edg,cond_epsr,edg_coord,edg_cond,edg_diel_coord,&
       edg_diel_epsr,edge_av_cond,edge_av_diel,edge_av,scale
  implicit none

  integer::ii,ids(2)
  !// real(kind=dp),allocatable::node_coord(:,:)

  scale=1.d0 ! Scale of the geometry

  ! Read the following info:
  ! (1) number of nodes for conductor: ncond_node
  ! (2) number of edges for conductor: ncond_edg
  ! (3) number of conductor: ncond
  ! (4) number of edges for dielectric (conformal and background dielectric): ndiel_edg
  read(16,*)ncond_node,ncond_edg,ncond,ndiel_edg
  if (ncond==0) then
     print*,'# of conductors can not be 0'
     stop
  else
     print*,'The # of conductors is: ',ncond
     print*,'The # of conductor edges is: ',ncond_edg
     nsuinf(1)=ncond_edg ! # of conductor edge
     ! store the permittivity of dielectric enclosing each conductor edge
     allocate(cond_epsr(1:ncond_edg))
     cond_epsr(1:ncond_edg)=1.d0 ! default value is 1 (free space)
  end if

  if (ndiel_edg==0) then
     nsuinf(2)=0
     print*,'No dielectric'
  else
     nsuinf(2)=ndiel_edg ! # of dielectric edge
     print*,'The # of dieletric edges is: ',ndiel_edg
  end if

  ! store coordinates of node for conductor
  allocate(node_coord(1:2,1:ncond_node))
  ! store coordinates of edge node for conductor
  allocate(edg_coord(1:2,1:2,1:ncond_edg))
  ! store relation between edge and conductor
  allocate(edg_cond(1:ncond_edg))
  node_coord(:,:)=0.d0
  edg_coord(:,:,:)=0.d0
  edg_cond(:)=0

  ! Read node of conductor
  do ii=1,ncond_node
     read(16,*)node_coord(1,ii),node_coord(2,ii)
     !// print*,ii,node_coord(1,ii),node_coord(2,ii)
  end do

  ! Read connectivity of node
  ! The edge for conductor is counter-clockwise
  do ii=1,ncond_edg
     read(16,*)ids(1),ids(2),edg_cond(ii),cond_epsr(ii)
     edg_coord(1:2,1,ii)=node_coord(1:2,ids(1))*scale
     edg_coord(1:2,2,ii)=node_coord(1:2,ids(2))*scale
!     print*,ii,edg_coord(1:2,1,ii),edg_coord(1:2,2,ii)
  end do

  !// Moved to "ky_simulate"
  !// deallocate(node_coord)
  !// ! Calculate the maximum and average edge length for conductors
  !// ! Find the first edge id for each conductor
  !// call ky_end_cond_edge

  ! store coordinates of edge node for dielectric
  allocate(edg_diel_coord(1:2,1:2,1:ndiel_edg))
  ! store outer/inner permittivity of edge node for dielectric
  allocate(edg_diel_epsr(1:2,1:ndiel_edg))
  edg_diel_epsr(:,:)=0.d0
  edg_diel_coord(:,:,:)=0.d0

  ! Read edge info for dielectric edge
  do ii=1,ndiel_edg
     read(16,*)edg_diel_coord(1,1,ii),edg_diel_coord(2,1,ii),edg_diel_coord(1,2,ii),&
          edg_diel_coord(2,2,ii),edg_diel_epsr(1,ii),edg_diel_epsr(2,ii)
     edg_diel_coord(1:2,1:2,ii)=edg_diel_coord(1:2,1:2,ii)*scale
!     print*,ii,edg_diel_coord(1:2,1,ii),edg_diel_coord(1:2,2,ii)
!     print*,ii,edg_diel_epsr(1:2,ii)
  end do

  !// Moved to "ky_simulate"
  !// ! Calculate the maximum and average edge length for dielectric
  !// call ky_end_diel_edge
  !// edge_av=(edge_av_cond*nsuinf(1)+edge_av_diel*nsuinf(2))/(nsuinf(1)+nsuinf(2))
  !// print*,'Average Edge length (m): ',edge_av
  !// print*,'Geometry for conductor and dielectric is generated'
  return
end subroutine ingeo

