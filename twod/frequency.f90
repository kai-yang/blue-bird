program cap2d_layered
  use global_com 
  use global_geom,only:nsuunk,nglunk,nsuinf,edg_dmg_epsr
  use layers,only:nlayers,ky_init_diel,ky_add_diel,ky_end_diel

  implicit none

  integer::num_diel,num_cond,num_edge,cond_id,num_dmg,num_edge_dmg
  integer::i,idx,jdx,num_x,num_z,num_x_dmg,num_z_dmg,start
  integer,allocatable::id_edg(:),id_edg_dmg(:),epsr_edg_dmg(:)
  real(kind=dp)::cap
  real(kind=dp),allocatable::coord_edg(:,:,:),xL(:),xH(:),zL(:),zH(:),dx(:),dz(:)
  real(kind=dp),allocatable::coord_edg_dmg(:,:,:),xL_dmg(:),xH_dmg(:),zL_dmg(:),zH_dmg(:),dx_dmg(:),dz_dmg(:)
  real(kind=dp),allocatable::zcoord(:),eps_r(:),epsr_dmg(:)

  open(unit=17,file='info.out',status='unknown')  

  ! Construct layered media info
  num_diel=3
  allocate(zcoord(1:num_diel+1),eps_r(1:num_diel))
  
  ! Cao 1
!!$  zcoord(1:num_diel-1)=(/0.d0,1.d0/)*2.54d-5 ! it doesn't include PEC coord, start from (zlow_of_layer(2))
!!$  zcoord(num_diel)=-1.d0
!!$  eps_r(:)=(/-1.d0,2.d0,1.d0/)

  ! Cao 2
!!$  zcoord(1:num_diel-1)=(/0.d0,0.5d0/)*2.54d-5 ! it doesn't include PEC coord
!!$  zcoord(num_diel)=-1.d0
!!$  eps_r(:)=(/-1.d0,6.8d0,1.d0/)

  ! Zutter
!!$  zcoord(1:num_diel-1)=(/0.d0,2.d2,3.d2/)*1.d-6 ! it doesn't include PEC coord
!!$  zcoord(num_diel)=-1.d0
!!$  eps_r(:)=(/-1.d0,4.3d0,3.2d0,1.d0/)

  ! Oh
!!$  zcoord(1:num_diel-1)=(/0.d0,6.d2,6.015d2/)*1.d-6 ! it doesn't include PEC coord
!!$  zcoord(num_diel)=-1.d0
!!$  eps_r(:)=(/-1.d0,11.d0,5.d0,1.d0/)

  ! Huang
  zcoord(1:num_diel-1)=(/0.d0,5.d0/)*2.54d-5 ! it doesn't include PEC coord
  zcoord(num_diel)=-1.d0
  eps_r(:)=(/-1.d0,5.d0,1.d0/)

  ! Initiate layered media info
  call ky_init_diel(num_diel)
  do i=1,num_diel
     call ky_add_diel(i,zcoord(i),eps_r(i))
  end do
  call ky_end_diel

  ! Construct conductor info
  num_cond=2
  call ky_init_ncond(num_cond)
  allocate(xL(1:num_cond),xH(1:num_cond),zL(1:num_cond),zH(1:num_cond))
  allocate(dx(1:num_cond),dz(1:num_cond))

  ! Cao 1
!!$  xL(:)=(/0.d0,5*2.54d-5/)
!!$  xH(:)=(/3*2.54d-5,8*2.54d-5/)
!!$  zL(:)=(/2.54d-5,2.54d-5/)
!!$  zH(:)=(/2*2.54d-5,2*2.54d-5/)

  ! Cao 2
!!$  xL(:)=(/-0.3d0,0.1d0/)*2.54d-5
!!$  xH(:)=(/-0.1d0,0.3d0/)*2.54d-5
!!$  zL(:)=(/0.6d0,0.2d0/)*2.54d-5
!!$  zH(:)=(/0.7d0,0.4d0/)*2.54d-5

  ! Zutter1
!!$  xL(:)=(/0.d0,5.d2,1.d3/)*1.d-6
!!$  xH(:)=(/3.5d2,8.5d2,1.35d3/)*1.d-6
!!$  zL(:)=(/2.d2,3.d2,3.d2/)*1.d-6
!!$  zH(:)=(/2.7d2,3.7d2,3.7d2/)*1.d-6

  ! Zutter2
!!$  xL(:)=(/0.d0,8.d0,1.9d1,3.3d1,5.d1,0.d0,8.d0,1.9d1,3.3d1,5.d1/)*1.d-6
!!$  xH(:)=(/2.d0,1.d1,2.1d1,3.5d1,5.2d1,2.d0,1.d1,2.1d1,3.5d1,5.2d1/)*1.d-6
!!$  zL(:)=(/6.d2,6.d2,6.d2,6.d2,6.d2,6.015d2,6.015d2,6.015d2,6.015d2,6.015d2/)*1.d-6
!!$  zH(:)=(/6.01d2,6.01d2,6.01d2,6.01d2,6.01d2,6.025d2,6.025d2,6.025d2,6.025d2,6.025d2/)*1.d-6

  ! Huang
  xL(:)=(/-4.d0,1.d0/)*2.54d-5
  xH(:)=(/-1.d0,4.d0/)*2.54d-5
  zL(:)=(/1.d0,1.d0/)*2.54d-5
  zH(:)=(/2.d0,2.d0/)*2.54d-5

  num_x=10; num_z=5
  dx(:)=(xH(:)-xL(:))/num_x
  dz(:)=(zH(:)-zL(:))/num_z

  num_edge=2*(num_x+num_z)*num_cond
  print*,'# of conductor edge is:',num_edge
  call ky_init_edge(num_edge)
  allocate(coord_edg(2,2,num_edge),id_edg(num_edge))
  ! cond
  do i=1,num_cond
     start=2*(num_x+num_z)*(i-1)
     do idx=0,num_z-1
        coord_edg(1,1,start+idx+1)=xL(i)
        coord_edg(2,1,start+idx+1)=zL(i)+idx*dz(i)
        coord_edg(1,2,start+idx+1)=xL(i)
        coord_edg(2,2,start+idx+1)=zL(i)+(idx+1)*dz(i)
        id_edg(start+idx+1)=i
     end do
     do idx=0,num_x-1
        coord_edg(1,1,start+num_z+idx+1)=xL(i)+idx*dx(i)
        coord_edg(2,1,start+num_z+idx+1)=zH(i)
        coord_edg(1,2,start+num_z+idx+1)=xL(i)+(idx+1)*dx(i)
        coord_edg(2,2,start+num_z+idx+1)=zH(i)
        id_edg(start+num_z+idx+1)=i
     end do
     do idx=0,num_z-1
        coord_edg(1,1,start+num_x+num_z+idx+1)=xH(i)
        coord_edg(2,1,start+num_x+num_z+idx+1)=zH(i)-idx*dz(i)
        coord_edg(1,2,start+num_x+num_z+idx+1)=xH(i)
        coord_edg(2,2,start+num_x+num_z+idx+1)=zH(i)-(idx+1)*dz(i)
        id_edg(start+num_x+num_z+idx+1)=i
     end do
     do idx=0,num_x-1
        coord_edg(1,1,start+num_x+2*num_z+idx+1)=xH(i)-idx*dx(i)
        coord_edg(2,1,start+num_x+2*num_z+idx+1)=zL(i)
        coord_edg(1,2,start+num_x+2*num_z+idx+1)=xH(i)-(idx+1)*dx(i)
        coord_edg(2,2,start+num_x+2*num_z+idx+1)=zL(i)
        id_edg(start+num_x+2*num_z+idx+1)=i
     end do
  end do

  do i=1,num_edge
!     print*,i,coord_edg(:,1,i),coord_edg(:,2,i)
!     print*,i,id_edg(i)
     call ky_add_edge(i,coord_edg(1,1,i),coord_edg(2,1,i),&
          coord_edg(1,2,i),coord_edg(2,2,i),id_edg(i))
  end do
  call ky_end_edge

  ! Construct damage info
  num_dmg=0
  call ky_init_ndmg(num_dmg)
  if (num_dmg/=0) then
     allocate(xL_dmg(1:num_dmg),xH_dmg(1:num_dmg),zL_dmg(1:num_dmg),zH_dmg(1:num_dmg))
     allocate(dx_dmg(1:num_dmg),dz_dmg(1:num_dmg))
     allocate(epsr_dmg(1:num_dmg))

     ! Huang
     xL_dmg(:)=(/-4.1d0,0.9d0/)*2.54d-5
     xH_dmg(:)=(/-0.9d0,4.1d0/)*2.54d-5
     zL_dmg(:)=(/0.9d0,0.9d0/)*2.54d-5
     zH_dmg(:)=(/2.1d0,2.1d0/)*2.54d-5
     epsr_dmg(1:num_dmg)=(/1.d0,1.d0/)

     num_x_dmg=10; num_z_dmg=5
     dx_dmg(:)=(xH_dmg(:)-xL_dmg(:))/num_x_dmg
     dz_dmg(:)=(zH_dmg(:)-zL_dmg(:))/num_z_dmg

     num_edge_dmg=2*(num_x_dmg+num_z_dmg)*num_dmg
     print*,'# of conformal dielectric edge is:',num_edge_dmg
     call ky_init_edge_dmg(num_edge_dmg)
     allocate(coord_edg_dmg(2,2,num_edge_dmg),id_edg_dmg(num_edge_dmg),epsr_edg_dmg(num_edge_dmg))

     ! conformal dielectric (the edge is clockwise)
     do i=1,num_dmg
        start=2*(num_x_dmg+num_z_dmg)*(i-1)
        do idx=0,num_z_dmg-1
           coord_edg_dmg(1,1,start+idx+1)=xL_dmg(i)
           coord_edg_dmg(2,1,start+idx+1)=zL_dmg(i)+idx*dz_dmg(i)
           coord_edg_dmg(1,2,start+idx+1)=xL_dmg(i)
           coord_edg_dmg(2,2,start+idx+1)=zL_dmg(i)+(idx+1)*dz_dmg(i)
           id_edg_dmg(start+idx+1)=i
           epsr_edg_dmg(start+idx+1)=epsr_dmg(i)
        end do
        do idx=0,num_x_dmg-1
           coord_edg_dmg(1,1,start+num_z_dmg+idx+1)=xL_dmg(i)+idx*dx_dmg(i)
           coord_edg_dmg(2,1,start+num_z_dmg+idx+1)=zH_dmg(i)
           coord_edg_dmg(1,2,start+num_z_dmg+idx+1)=xL_dmg(i)+(idx+1)*dx_dmg(i)
           coord_edg_dmg(2,2,start+num_z_dmg+idx+1)=zH_dmg(i)
           id_edg_dmg(start+num_z_dmg+idx+1)=i
           epsr_edg_dmg(start+num_z_dmg+idx+1)=epsr_dmg(i)
        end do
        do idx=0,num_z_dmg-1
           coord_edg_dmg(1,1,start+num_x_dmg+num_z_dmg+idx+1)=xH_dmg(i)
           coord_edg_dmg(2,1,start+num_x_dmg+num_z_dmg+idx+1)=zH_dmg(i)-idx*dz_dmg(i)
           coord_edg_dmg(1,2,start+num_x_dmg+num_z_dmg+idx+1)=xH_dmg(i)
           coord_edg_dmg(2,2,start+num_x_dmg+num_z_dmg+idx+1)=zH_dmg(i)-(idx+1)*dz_dmg(i)
           id_edg_dmg(start+num_x_dmg+num_z_dmg+idx+1)=i
           epsr_edg_dmg(start+num_x_dmg+num_z_dmg+idx+1)=epsr_dmg(i)
        end do
        do idx=0,num_x_dmg-1
           coord_edg_dmg(1,1,start+num_x_dmg+2*num_z_dmg+idx+1)=xH_dmg(i)-idx*dx_dmg(i)
           coord_edg_dmg(2,1,start+num_x_dmg+2*num_z_dmg+idx+1)=zL_dmg(i)
           coord_edg_dmg(1,2,start+num_x_dmg+2*num_z_dmg+idx+1)=xH_dmg(i)-(idx+1)*dx_dmg(i)
           coord_edg_dmg(2,2,start+num_x_dmg+2*num_z_dmg+idx+1)=zL_dmg(i)
           id_edg_dmg(start+num_x_dmg+2*num_z_dmg+idx+1)=i
           epsr_edg_dmg(start+num_x_dmg+2*num_z_dmg+idx+1)=epsr_dmg(i)
        end do
     end do

     do i=1,num_edge_dmg
!        print*,i,coord_edg_dmg(:,1,i),coord_edg_dmg(:,2,i)
!        print*,i,id_edg_dmg(i),epsr_edg_dmg(i)
        call ky_add_edge_dmg(i,coord_edg_dmg(1,1,i),coord_edg_dmg(2,1,i),&
             coord_edg_dmg(1,2,i),coord_edg_dmg(2,2,i),id_edg_dmg(i),epsr_edg_dmg(i))
     end do
     call ky_end_edge_dmg
  end if

  call ky_simulate
  do idx=0,num_cond-1
     do jdx=idx,num_cond-1
        call ky_getC(idx,jdx,cap)
     end do
  end do
  call ky_clear_edge
  call ky_clear_all
     
  print*,'Closing files-safely here.'
  close(17,status='keep')                                        
  stop
end program cap2d_layered


