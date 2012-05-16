subroutine ky_simulate
  use global_com 
  use global_geom,only:nsuunk,nglunk,nsuinf,edge_av
  use global_dim,only:pmatrix,rj
  use quadratures,only:determine_quadrature,nqp_t
  use layers,only:init_layers,find_rho_z,&
       init_interpolation,fill_Green_stored_array
  use mat_vec_mult,only:initialize_r0
  implicit none

  real(kind=dp)::tim_start,tim_gf,tim_all,tim_dirfield,tim_extra
  integer::Itim_start,Itim_gf,Itim_all,Itim_dirfield,Itim_extra

  call inmom   ! read in the input parameters for the MOM algorithm

  call system_clock(COUNT=Itim_start,COUNT_RATE=Itim_rate,COUNT_MAX=Itim_max)
  tim_start=real(Itim_start)/real(Itim_rate)
  print*,'As long as time is less than',real(Itim_max)/real(Itim_rate),'secs, timing is OK'
  print*,'---------------------------------------------------------- '
  print*,'Serial Capacitance extraction by K. Yang'
  print*,'---------------------------------------------------------- '
  write(17,*) '-------------------------------------------------------- '
  write(17,*) 'Serial Capacitance extraction by K. Yang'
  write(17,*) '-------------------------------------------------------- '

!  call insu    ! surface discretization info.
  R_a=factor2*edge_av
  print *,'R_a=',R_a
  write(17,*) 'Distance for singularity extraction (R_a)=',R_a
  nsuunk=nsuinf(1)
  nglunk=nsuinf(1)

  if (nsuunk>0) then
     call determine_quadrature
     call find_rho_z
  end if

  if (is_iter) then
     call precon_init
  end if

  call init_layers
  call init_interpolation
  call system_clock(Itim_extra,Itim_rate)
  tim_extra=real(Itim_extra)/real(Itim_rate)
  print*,'TIMING::::::::::Extra',tim_extra-tim_start
  call fill_Green_stored_array

  call system_clock(Itim_gf,Itim_rate)
  tim_gf=real(Itim_gf)/real(Itim_rate)
  print*,'TIMING::::::::::GF table',tim_gf-tim_extra

  print*,'Entering to dirfield'
  call dir_field_mom_only
  print*,'out of dirfield'

call system_clock(Itim_dirfield,Itim_rate);tim_dirfield=real(Itim_dirfield)/real(Itim_rate)
  print*,'TIMING::::::::::Dirfield',tim_dirfield-tim_gf
  print*,'TIMING::::::::::ALL-BEFORE-SOLVE=',tim_dirfield-tim_start
  write(17,*) 'TIMING::::::::::Dirfield',tim_dirfield-tim_gf
  write(17,*) 'TIMING::::::::::ALL-BEFORE-SOLVE=',tim_dirfield-tim_start
!*************************************************************
!****  First synchronization point
!*************************************************************
  if (is_iter) then
     call var_arrays     ! allocates variable arrays 
     call initialize_r0
  else
     call invert_pmatrix ! direct solver
  end if
  call cal_C
  call system_clock(Itim_all,Itim_rate)
  tim_all=real(Itim_all)/real(Itim_rate)

  print*,'TIMING::::::::::TOTAL Solve time',tim_all-tim_dirfield
  print*,'TIMING::::::::::TOTAL Cap',tim_all-tim_start
  
  write(17,*) 'TIMING::::::::::TOTAL Solve time',tim_all-tim_dirfield
  write(17,*) 'TIMING::::::::::TOTAL Cap',tim_all-tim_start
  return
end subroutine ky_simulate


subroutine ky_num_node_num_edge(num_node, num_edge)
  use global_com,only:dp
  use global_geom,only:nsuinf, sunod, nsuedgn, add_point_ptr, add_edge_ptr
  implicit none

  integer,intent(in)::num_node, num_edge

  nsuinf(1)=num_node
  nsuinf(2)=num_edge
  
  allocate(sunod(2,num_node)) 
  allocate(nsuedgn(2,num_edge))
  sunod(:,:)=0.0_dp
  nsuedgn(:,:)=0
  add_point_ptr = 1
  add_edge_ptr = 1
  return
end subroutine ky_num_node_num_edge

subroutine ky_add_point(xx,yy)
  use global_com,only:dp
  use global_geom,only:add_point_ptr, sunod
  implicit none

  real(kind=dp),intent(in)::xx,yy

  !assert(this_pnt < nnod+1)
  sunod(1,add_point_ptr)=xx*2.54d-5
  sunod(2,add_point_ptr)=yy*2.54d-5
  add_point_ptr = add_point_ptr + 1
  return
end subroutine ky_add_point

subroutine ky_add_edge(from,to)
  use global_geom,only:add_edge_ptr, nsuedgn
  implicit none

  integer,intent(in)::from,to

  !assert(add_edge_ptr < nedg)
  nsuedgn(1,add_edge_ptr) = from
  nsuedgn(2,add_edge_ptr) = to
  add_edge_ptr = add_edge_ptr + 1
  return
end subroutine ky_add_edge

subroutine ky_is_layers(is_layers)
  use layers,only:is_multilayer
  implicit none

  logical,intent(in)::is_layers

  is_multilayer=is_layers
  return
end subroutine ky_is_layers

subroutine ky_num_layers(num_layers)
  use global_com,only:dp,real_mem,complex_mem
  use layers,only:nlayers,h_of_layer,zlow_of_layer,eps_t,&
       Z_0,GammaL_mn,GammaR_mn,kz_wave,k_prop2

  integer,intent(in)::num_layers

  real(kind=dp)::mem_est

  nlayers = num_layers
  if (nlayers==1) then
     print*,'The layered medium is only free space'
  else if (nlayers<1) then
     print*,'ERROR::NLAYERS MUST BE A POSITIVE INTEGER GREAT THAN 1!'
     stop
  else
     print*,'The number of layered medium is: ',nlayers
  end if
  
  mem_est=(4.d0*nlayers+1)*real_mem+(5.d0*nlayers)*complex_mem
  print*,'The layered medium requires memory (MB): ',mem_est/1024.d0/1024.d0
  ! Essential array
  allocate(h_of_layer(1:nlayers),zlow_of_layer(1:nlayers+1),eps_t(1:nlayers),&
       Z_0(1:nlayers),GammaL_mn(1:nlayers),GammaR_mn(1:nlayers),&
       kz_wave(1:nlayers),k_prop2(1:nlayers))
  return
end subroutine ky_num_layers

subroutine ky_set_zref(zref)
  use global_com,only:dp
  use layers,only:zlow_of_layer

  zlow_of_layer(2)=zref
  return
end subroutine ky_set_zref

subroutine ky_set_layer(ilayer,eps,height)
  use global_com,only:dp
  use layers,only:nlayers,h_of_layer,zlow_of_layer,eps_t

  integer,intent(in)::ilayer
  real(kind=dp),intent(in)::eps,height

  integer::i

  i=ilayer
  eps_t(i)=eps
  h_of_layer(i)=height
  if (i==1) then
     if (h_of_layer(i)==-1.d0) then           
        zlow_of_layer(i)=-1.d0 ! the height of the bottom layer is infinite (half space)
     else
        zlow_of_layer(i)=zlow_of_layer(2)  ! the height of the bottom layer is 0 (PEC or PMC)
     end if
  else if (i==nlayers) then
     if (h_of_layer(i)==-1.d0) then           
        zlow_of_layer(i+1)=-1.d0 ! the height of the top layer is infinite
     else
        zlow_of_layer(i+1)=zlow_of_layer(i-1)+h_of_layer(i-1)  ! the height of the top layer is 0 (PEC or PMC)
     end if
  else
     ! when layer_s==nlayers, zlow_of_layer(nlayers+1) should be handled carefully
     zlow_of_layer(i+1)=zlow_of_layer(i)+h_of_layer(i) 
  end if

  return
end subroutine ky_set_layer

subroutine ky_set_tol(tol)
  use global_com,only:dp
  use layers,only:threshold

  real(kind=dp),intent(in)::tol

  threshold=tol
  return
end subroutine ky_set_tol

subroutine ky_set_misc
  use global_com,only:dp
  use layers,only:eps_t_max,eps_t,C_s

  eps_t_max=maxval(eps_t(:)) 
  C_s(1:4)=(/1.d0,1.d0,1.d0,1.d0/)
  return
end subroutine ky_set_misc

subroutine parse_layers(file_no)
  use global_com,only:dp
  use layers,only:nlayers,h_of_layer,zlow_of_layer,eps_t,&
       Z_0,GammaL_mn,GammaR_mn,kz_wave,k_prop2
  implicit none

  integer,intent(in)::file_no

  logical::is_layers
  integer::num_layers,i
  real(kind=dp)::eps,height,zref,tol

  read(file_no,*) is_layers
  call ky_is_layers(is_layers)

  if (.not. is_layers) then  
     num_layers=1
  else
     read(file_no,*) num_layers
  end if
  call ky_num_layers(num_layers)

  if (.not. is_layers) then           
     eps_t(1)=1.d0; h_of_layer(1)=-1.d0
     zlow_of_layer(1)=-1.d0; zlow_of_layer(2)=-1.d0
  else
     read(file_no,*) zref
     call ky_set_zref(zref)     
     ! In order to handle the infinity case, we use the following notations.
     ! h_of_layer(i)==-1.0 means the height is infinite
     ! In order to handle PEC case, we use the following notation: eps_t=-1 (Inf), Z_0=0
     do i=1,nlayers
        read(file_no,*) eps,height
        call ky_set_layer(i,eps,height)
     end do
     read(file_no,*) tol
     call ky_set_tol(tol)
  end if
  call ky_set_misc
  return
end subroutine parse_layers

subroutine parse_geom(file_no)
  use global_com,only:dp
  implicit none

  integer,intent(in)::file_no

  real(kind=dp)::xx,yy
  integer::from,to,j,nnod,nedg

  read(file_no,*) nnod,nedg
  call ky_num_node_num_edge(nnod, nedg)
  do j=1,nnod                                                  
     read(file_no,*) xx,yy
     call ky_add_point(xx,yy)
  end do
  do j=1,nedg
     read(file_no,*) from,to
     call ky_add_edge(from,to)
  end do  
  return
end subroutine parse_geom


