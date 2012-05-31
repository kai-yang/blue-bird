subroutine insu
!
!  read the number of surfaces and for every surface the number of
!  nodes and the coordinates of those nodes and the number of patches
!  and the nodes of every patch.
!
  use global_geom,only:nsuinf,nsuedgn,&
       nsuunk,sunod1,sunod,nglunk,edge_av
  use global_com,only:closed_body,dp,real_mem,int_mem,ncond
  use misc_dbl,only:cened_dbl
  implicit none
  integer::nnod,nedg,nbe,n1,n2,n3,k,j,i,ik
  integer::nt1,nt2,nt3
  logical::edge_found
  integer,dimension(2)::no
  integer,allocatable::temp_junction(:)
  logical,allocatable::is_not_a_junction_edge(:)
  integer::patch_no,edge_no,node,node1,node2,node_no
  integer,allocatable::count_nodeedg(:)
  integer,allocatable::dummy_nodeedg(:),dummy_nodeedg2(:)
  type edges_of_nodes
     integer,pointer::p(:),ed(:)
  end type edges_of_nodes
  type(edges_of_nodes),allocatable::nsuned(:)
  ! pec body parameters
  real(kind=dp),allocatable::super_edp1(:),super_edp2(:)
  real(kind=dp)::dummy(3)
  real(kind=dp)::max_leng,av_leng,leng


!     read the number of nodes and patches 
!  read(11,*) nnod,nedg
!  nsuinf(1)=nnod                                                  
!  nsuinf(2)=nedg
  nnod = nsuinf(1)
  nedg = nsuinf(2)
  !assert(nnod > 0)
  !assert(nedg > 0)

! allocate Arrays
!  print*,'yaba',nnod,npat,nedg
!  allocate(sunod(2,nnod)) 
!  allocate(nsuedgn(2,nedg))
! initialize 
!  sunod(:,:)=0.0_dp
!  nsuedgn(:,:)=0
!                                                                       
!  read the coordinates of all the nodes                                
!                                                                       
 ! do j=1,nnod                                                  
 !    read(11,*) sunod(1,j),sunod(2,j)
 !    sunod(:,j)=sunod(:,j)*2.54d-5
 ! end do
!                                                                       
!     read the nodes for every patch                                    
!                                                                       
!  do j=1,nedg
!     read(11,*) nsuedgn(1,j),nsuedgn(2,j)
!  end do



  if (ncond/=1) then
     call reorder_edge
  end if

  av_leng=0.d0;max_leng=0.d0
  do j=1,nedg
     call cened_dbl(j,leng)
     av_leng=av_leng+leng
     max_leng=max(max_leng,leng)
  end do
  av_leng=av_leng/nedg
  edge_av=av_leng
  print*,'# of edges, Maximum and Average Edge length (m):',&
       nedg,max_leng,av_leng
  write(17,*) '# of edges,Maximum and Average Edge length (m):',&
       nedg,max_leng,av_leng
  return
end subroutine insu

subroutine reorder_edge
  use global_com,only:dp,ncond
  use global_geom,only:nglunk,nsuinf,nsuedgn,npat_cond,sunod

  integer::dummy,count,idx,cond_id
  real(kind=dp)::rm(2)
  real(kind=dp),allocatable::nsuedgn_tmp(:,:)
  type cap
     integer,pointer::p(:) !patch ids
  end type cap
  type(cap),allocatable::pat_cond(:)

  allocate(npat_cond(1:ncond))
  npat_cond(1:ncond)=0
  allocate(pat_cond(1:ncond))

  if (ncond==1) then
     npat_cond(1)=nsuinf(2)
  else 
     ! count number of patches belonging to different conductor
     do dummy=1,nsuinf(2)
        cond_id=nsuedgn(3,dummy)
        npat_cond(cond_id)=npat_cond(cond_id)+1
     end do
  end if
  
  do dummy=1,ncond
     allocate(pat_cond(dummy)%p(1:npat_cond(dummy)))
     print*,npat_cond(dummy),'lines belong to Conductor is',dummy
  end do
  print*,'The total number of lines is ',nsuinf(2)

  ! Re-initialize to 0
  npat_cond(1:ncond)=0
  if (ncond==1) then
     do dummy=1,nsuinf(2)
        pat_cond(1)%p(dummy)=dummy
     end do
  else 
     do dummy=1,nsuinf(2)
        cond_id=nsuedgn(3,dummy)
        npat_cond(cond_id)=npat_cond(cond_id)+1
        pat_cond(cond_id)%p(npat_cond(cond_id))=dummy
     end do
  end if

  ! convert different points to a 1D array based on the conductor ID
  allocate(nsuedgn_tmp(2,nsuinf(2)))
  idx=0
  do dummy=1,ncond
     do count=1,npat_cond(dummy)
        idx=idx+1
        nsuedgn_tmp(:,idx)=nsuedgn(:,pat_cond(dummy)%p(count))
     end do
  end do
  nsuedgn=nsuedgn_tmp
  deallocate(nsuedgn_tmp)
  return
end subroutine reorder_edge
