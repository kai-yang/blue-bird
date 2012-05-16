subroutine ky_init_ncond(num_cond)
  use global_com,only:ncond
  implicit none

  integer,intent(in)::num_cond

  ncond=num_cond
  if (ncond==0) then
     print*,'# of conductors can not be 0'
     stop
  else
     print*,'The # of conductors is: ',ncond
  end if

  return
end subroutine ky_init_ncond

subroutine ky_init_ndmg(num_dmg)
  use global_com,only:ndmg,ncond
  use global_geom,only:nsuinf,cond_epsr
  implicit none

  integer,intent(in)::num_dmg

  ndmg=num_dmg
  if (ndmg==0) then
     print*,'No conformal dielectric'
     nsuinf(2)=0
  else
     print*,'The # of conformal dieletrics is: ',ndmg
     allocate(cond_epsr(1:ncond))
     cond_epsr(1:ncond)=1.d0 ! default value is 1
  end if

  return
end subroutine ky_init_ndmg

subroutine ky_init_edge(num_edge)
  use global_geom,only:nsuinf,edg_cond,edg_coord
  use global_com,only:dp,ncond
  implicit none

  integer,intent(in)::num_edge

  integer::nedg

  nedg=num_edge
  nsuinf(1)=nedg ! # of conductor edge
  allocate(edg_coord(1:2,1:2,1:nedg))
  allocate(edg_cond(1:nedg))
  edg_cond(:)=0
  edg_coord(:,:,:)=0.d0
  return
end subroutine ky_init_edge

subroutine ky_init_edge_dmg(num_edge_dmg)
  use global_geom,only:nsuinf,edg_dmg,edg_dmg_coord,edg_dmg_epsr
  use global_com,only:dp,ndmg
  implicit none

  integer,intent(in)::num_edge_dmg

  integer::nedg_dmg

  nedg_dmg=num_edge_dmg
  nsuinf(2)=nedg_dmg ! # of conformal dielectric edge
  allocate(edg_dmg_coord(1:2,1:2,1:nedg_dmg))
  allocate(edg_dmg(1:nedg_dmg))
  allocate(edg_dmg_epsr(1:nedg_dmg))
  edg_dmg(:)=0
  edg_dmg_epsr(:)=0.d0
  edg_dmg_coord(:,:,:)=0.d0
  return
end subroutine ky_init_edge_dmg

subroutine ky_add_edge(idx,x1,y1,x2,y2,cond_id)
  use global_geom,only:nsuinf,edg_cond,&
       edg_coord
  use global_com,only:dp,ncond
  implicit none

  integer,intent(in)::idx,cond_id
  real(kind=dp),intent(in)::x1,y1,x2,y2

  edg_cond(idx)=cond_id
  edg_coord(1:2,1,idx)=(/x1,y1/)
  edg_coord(1:2,2,idx)=(/x2,y2/)
  if (cond_id/=1) then ! check validity of input
     if (idx > 1 .and. edg_cond(idx)<edg_cond(idx-1)) then
        print*,'The cond_id is not in order. Please reorder them', edg_cond(idx), edg_cond(idx-1)
        stop
     end if
     if (cond_id>ncond) then
        print*,'cond_id is larger than # of conductor',ncond
        stop
     end if
  end if

  return
end subroutine ky_add_edge

subroutine ky_add_edge_dmg(idx,x1,y1,x2,y2,dmg_id,epsr)
  use global_geom,only:nsuinf,edg_dmg,&
       edg_dmg_coord,edg_dmg_epsr
  use global_com,only:dp,ndmg
  implicit none

  integer,intent(in)::idx,dmg_id
  real(kind=dp),intent(in)::x1,y1,x2,y2,epsr

  edg_dmg(idx)=dmg_id
  edg_dmg_epsr(idx)=epsr
  edg_dmg_coord(1:2,1,idx)=(/x1,y1/)
  edg_dmg_coord(1:2,2,idx)=(/x2,y2/)
  if (dmg_id/=1) then ! check validity of input
     if (idx > 1 .and. edg_dmg(idx)<edg_dmg(idx-1)) then
        print*,'The diel_id is not in order. Please reorder them', &
             edg_dmg(idx), edg_dmg(idx-1)
        stop
     end if
     if (dmg_id>ndmg) then
        print*,'dmg_id is larger than # of conformal dielectric',ndmg
        stop
     end if
  end if

  return
end subroutine ky_add_edge_dmg

subroutine ky_end_edge
  use global_geom,only:nsuinf,edg_cond,&
       edg_coord,cond_sta,edge_av
  use global_com,only:dp,ncond
  use misc_dbl,only:cened_dbl
  implicit none

  integer::j
  integer,allocatable::nedg_cond(:)
  real(kind=dp)::av_leng,max_leng,leng

  allocate(cond_sta(1:ncond+1)) ! There are ncond ranging from 0 to ncond-1
  allocate(nedg_cond(1:ncond))
  cond_sta(:)=0
  nedg_cond(:)=0

  av_leng=0.d0;max_leng=0.d0
  do j=1,nsuinf(1)
     nedg_cond(edg_cond(j))=nedg_cond(edg_cond(j))+1
     call cened_dbl(j,leng)
     av_leng=av_leng+leng
     max_leng=max(max_leng,leng)
  end do
  av_leng=av_leng/nsuinf(1)
  edge_av=av_leng
  print*,'# of conductor edges, Maximum and Average Edge length (m):',&
       nsuinf(1),max_leng,av_leng
  write(17,*) '# of conductor edges,Maximum and Average Edge length (m):',&
       nsuinf(1),max_leng,av_leng
  cond_sta(1)=1
  do j=1,ncond
     cond_sta(j+1)=cond_sta(j)+nedg_cond(j)
  end do
  deallocate(nedg_cond)
  return
end subroutine ky_end_edge

subroutine ky_end_edge_dmg
  use global_geom,only:nsuinf,edge_av_dmg,edg_dmg,&
       edg_dmg_epsr,cond_epsr
  use global_com,only:dp
  use misc_dbl,only:cened_dmg_dbl
  implicit none

  integer::j
  real(kind=dp)::av_leng,max_leng,leng

  av_leng=0.d0;max_leng=0.d0
  do j=1,nsuinf(2)
     call cened_dmg_dbl(j,leng)
     av_leng=av_leng+leng
     max_leng=max(max_leng,leng)
     cond_epsr(edg_dmg(j))=edg_dmg_epsr(j)
  end do
  av_leng=av_leng/nsuinf(2)
  edge_av_dmg=av_leng
  print*,'# of conformal dielectrc edges, Maximum and Average Edge length (m):',&
       nsuinf(2),max_leng,av_leng
  write(17,*) '# of conductor edges,Maximum and Average Edge length (m):',&
       nsuinf(2),max_leng,av_leng
  return
end subroutine ky_end_edge_dmg
