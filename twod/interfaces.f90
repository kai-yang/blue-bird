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
  
  call ky_set_zref(0.0)
  call ky_set_misc
  call ky_set_tol(1.d-4)

  !call inmom   ! read in the input parameters for the MOM algorithm

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
  print *, num_node, num_edge
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
  print *, sunod(:, add_point_ptr)
  add_point_ptr = add_point_ptr + 1
  return
end subroutine ky_add_point

subroutine ky_add_edge(from,to,cid)
  use global_geom,only:add_edge_ptr, nsuedgn
  implicit none

  integer,intent(in)::from,to, cid

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
  real(kind=dp),intent(in)::zref

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
  integer::from,to,j,nnod,nedg,cid

  read(file_no,*) nnod,nedg
  call ky_num_node_num_edge(nnod, nedg)
  do j=1,nnod                                                  
     read(file_no,*) xx,yy
     call ky_add_point(xx,yy)
  end do
  do j=1,nedg
     read(file_no,*) from,to,cid
     call ky_add_edge(from,to,cid)
  end do  
  return
end subroutine parse_geom



subroutine invert_pmatrix
  use global_com 
  use global_dim,only:pmatrix
  use global_geom,only:nglunk,nsuinf
  use mat_vec_mult,only:tot_near_time,tot_reduction_time,&
       invert_preconditioner,matvec,initialize_r0
!                         
  implicit none
  real(kind=dp)::mem_est(2)
  ! matrix inversion using LU decomposition
  integer,allocatable::ipiv1(:)
  real(kind=dp),allocatable::work(:)
  integer::info(2)

!*****************************************************************************
!* START OF INITIALIZATION                                                   *
!*****************************************************************************
!***************************************
! For the first incidence, allocate arrays
mem_est(2)=mem_est(1)
  allocate(ipiv1(1:nsuinf(2)),work(1:nsuinf(2)))
!* END OF INITIALIZATION                                                     *
!*****************************************************************************
!*****************************************************************************
!* START OF Iterative Solver
!*****************************************************************************
  print*,'Start to directly solving'
  call DGETRF(nsuinf(2),nsuinf(2),pmatrix,nsuinf(2),ipiv1,info(1))
  call DGETRI(nsuinf(2),pmatrix,nsuinf(2),ipiv1,work,nsuinf(2),info(2))
  if (info(1)/=0 .or. info(2)/=0) then
     print*,'inversion of matrix info',info(1),info(2)
     stop
  else
     print*,'potential matrix has been inverted'
  end if
!*****************************************************************************
!* END OF Iterative Solver
!*****************************************************************************
  return
end subroutine invert_pmatrix

subroutine solve(rhsd)
  use global_com 
  use global_dim,only:rj
  use global_geom,only:nglunk
  use mat_vec_mult,only:tot_near_time,tot_reduction_time,&
       invert_preconditioner,matvec,initialize_r0
!                         
  implicit none
  complex(kind=dp),intent(in)::rhsd(1:nglunk)
  real(kind=dp)::err,rerr_local,rerr_sum,bmag_local
  integer::it,me,iter,dummy,i
  integer::k,j,dumb(1)
  real(kind=dp)::mem_est(2)

!*****************************************************************************
!* START OF INITIALIZATION                                                   *
!*****************************************************************************
!***************************************
! For the first incidence, allocate arrays, form and pre-FFT AIM matrices
mem_est(2)=mem_est(1)

call invert_preconditioner
maxit=max(nglunk,1000) ! do at least 1000 iterations before giving up
!     maxit=10000
call initialize_r0
  !***************************************
  ! Initialize the arrays
  rj(1:nglunk)=cmplx(0.0_dp,0.0_dp,dp)
  ! else the initial guess is the solution at the previous frequency
 !***************************************
!*****************************************************************************
!* END OF INITIALIZATION                                                     *
!*****************************************************************************
!*****************************************************************************
!* START OF Iterative Solver
!*****************************************************************************
  print*,'Start to iteratively solving'
  iter=maxit; err=dtol
  call ztfqmr_serial(nglunk,rhsd(1:nglunk),rj(1:nglunk),err,iter) 
!*****************************************************************************
!* END OF Iterative Solver
!*****************************************************************************
  return
end subroutine solve

subroutine cal_C
  use global_com,only:dp,ncond,is_iter,eps0d
  use global_dim,only:pmatrix,rj
  use global_geom,only:nglunk,nsuinf,npat_cond
  use misc_dbl,only:cened_dbl
!                         
  implicit none
  integer::dummy,count,idx
  integer,allocatable::cond_sta(:)
  real(kind=dp)::rm(3),leng
  real(kind=dp),allocatable::tot_Q(:,:)
  complex(kind=dp)::rhsd(1:nglunk)

  allocate(cond_sta(1:ncond+1))
  if (ncond==1) then
     cond_sta(1:2)=1
  else
     cond_sta(1)=1
     do dummy=1,ncond
        cond_sta(dummy+1)=cond_sta(dummy)+npat_cond(dummy)
     end do
  end if

  allocate(tot_Q(1:ncond,1:ncond))
  tot_Q(:,:)=0.d0
  if (is_iter) then
    if (ncond==1) then
       rhsd(:)=cmplx(1.0_dp,0.0_dp,dp)
       call solve(rhsd)   ! iterative method-of-moments solver
       do count=1,nsuinf(2)
          call cened_dbl(count,leng)
          tot_Q=tot_Q+rj(count)*leng
       end do
    else
       do dummy=1,ncond
          do count=1,ncond
             rhsd(:)=cmplx(0.0_dp,0.0_dp,dp)
             rhsd(cond_sta(count):cond_sta(count+1)-1)=cmplx(1.0_dp,0.0_dp,dp)
             call solve(rhsd)   ! iterative method-of-moments solver
             do idx=cond_sta(dummy),cond_sta(dummy+1)-1
                call cened_dbl(idx,leng)
                tot_Q(dummy,count)=tot_Q(dummy,count)+rj(idx)*leng
             end do
          end do
       end do
    end if
  else
     ! Multiply with proper voltage to obtain charge density on each conductor
     ! For 1 conductor, charge density is the summation of cmatrix
     if (ncond==1) then
        do dummy=1,nsuinf(2)
           call cened_dbl(dummy,leng)
           tot_Q=tot_Q+sum(pmatrix(dummy,1:nsuinf(2)))*leng
        end do
     else
        do dummy=1,ncond
           do idx=cond_sta(dummy),cond_sta(dummy+1)-1
              call cened_dbl(idx,leng)
              do count=1,ncond
                 tot_Q(dummy,count)=tot_Q(dummy,count)+&
                      sum(pmatrix(idx,cond_sta(count):cond_sta(count+1)-1))*leng
              end do
           end do
        end do
     end if
  end if
  do dummy=1,ncond
     do count=1,ncond
        print*,'Cap',dummy,count,'is',tot_Q(dummy,count)
     end do
  end do
end subroutine cal_C
!                                                                       
!*********************************************************************
!                                                                       
subroutine inmom                                                   
  ! This subroutine inputs algorithmic and electromagnetic information from
  ! mom.inp, and does some of the preprocesseing necessary to construct the 
  ! incident pulse. 
  ! Continually modified as the code expands 
  use global_com 
  use global_geom,only:diel_const,loss_factor
  use misc_dbl,only:rvpcr_dbl
  use quadratures,only:nqp_s,nqp_t,total_maxqp_t
  implicit none
  real(kind=dp)::epsilon_r,mu_r,dum1,dum2
  integer::i,su_order(2),noqp(2)

  pid=4.0_dp*atan(1.0_dp) 

  read(12,*) factor2     !factor2
  read(12,*) ncond       !ncond
  read(12,*) is_iter     !direct solver or iterative solve
  
  if (is_iter) then
     read(12,*) dtol     !tolerance of iterative solve
  end if
  read (12,*) diel_const,loss_factor  ! relative dielectic constant and loss factor

  read(12,*) istore_cur !if istore_cur == 1 then store all cur. 
  read(12,*) epsilon_r,mu_r,loss_sigma
  if (loss_sigma/=0.0_dp) then
     print*,'lossy part is incomplete... correct the fielde_dbl and aim.f90'
     stop
  end if
  rmu0d=pid*4.d-7*mu_r; cd=299792458._dp;eps0d=epsilon_r/(rmu0d*cd**2)
  cd=cd/sqrt(epsilon_r*mu_r); vlite=cd; eta0d=sqrt(rmu0d/eps0d)
  wod=2.d0*pid*1.d1

  ! Orders of integrations in the order: surface,wire,volume
  read(12,*) su_order(1),su_order(2) 
  read(12,*) i, prcd_blocksize, precon_dist
  if (i==1) then
     diag_precon=.true.;block_diag_precon=.false.;near_field_precon=.false.
     group_precon=.false.
     print*,'prcd_blocksize and precon_dist are meaningless for diag_precon'
  elseif(i==2) then
     diag_precon=.false.;block_diag_precon=.true.; near_field_precon=.false.
     group_precon=.false.
     print*,'prcd_blocksize is # of unknowns in each block, &
          & precon_dist is meaningless for block_diag_precon'
  elseif(i==3) then
     print*,'This pre-conditioner (near-field-precon) is not implemented here'
     stop
!     diag_precon=.false.;block_diag_precon=.false.; near_field_precon=.true.
!     group_precon=.false.
!     precon_dist=lambda_fmax*precon_dist
!     print*,'prcd_blocksize is meaningless for near_field_precon &
!          & neighboorhood radius is precon_dist*lambda_fmax'
  elseif(i==4) then
     print*,'This pre-conditioner (block_diag_precon) is not implemented here'
     stop
!     diag_precon=.false.;block_diag_precon=.false.; near_field_precon=.false.
!     group_precon=.true. 
!     precon_dist=lambda_fmax*precon_dist
!     print*,'prcd_blocksize is meaningless for group_precon &
!          & group radius is precon_dist*lambda_fmax'
     ! prcd_blocksize measures the # of cells used for determining the
     ! blocks, when grouping unknowns, prcd_dist, is the length-scale     
     ! used for constructing the blocks, i.e.,lambda*precon_dist
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! All the inputs in unit 12 are read here (except for multiple excitation data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,2
     if (su_order(i)<=10) then
        noqp(i)=su_order(i)
     else
        noqp(i)=3
     end if
  end do
  nqp_s=noqp(1); nqp_t=noqp(2); 
  total_maxqp_t=nqp_t

  write(17,*) 'factor2(rnear=factor2*lambda_min):',factor2
  
  if (is_iter) then
     write(17,*) 'Iterative solver tolerance', dtol
  end if
  write(17,*) 'Epsilon_r,mu_r,loss_sigma:',epsilon_r,mu_r,loss_sigma
  write(17,*) 'Order and # of quadrature points(source,test):',su_order(1),nqp_s,su_order(2),nqp_t
  write(17,*) 'The blocksize of the preconditioner is:',prcd_blocksize
  if (loss_sigma>0.0_dp) then
     print*,'skin depth is(m):',vlite*sqrt(2.0_DP)/wod/&
          sqrt(sqrt(1.0_DP+(loss_sigma/wod/eps0d)**2)-1.0_DP)
  end if
  return
end subroutine inmom

subroutine precon_init
  use global_com!,only:nprcd_blocks,prcd_blocksize,last_block_size,myid,&
!       numprocs,my_precon_ids,nprecon_ids,me_to_precon,dp,diag_precon,&
!       block_diag_precon,near_field_precon,group_precon,max_INF
  use global_dim,only:ipiv,prcdin
  use global_geom,only:nglunk,nsuunk
  implicit none

!************************************************
!**************** PRE-CONDITIONER ***************
!************************************************ 
  if (diag_precon) then
     allocate(prcdin(nglunk,1,1))
  elseif(block_diag_precon) then
     nprecon_ids=nglunk
     print*,'preconditions this many basis:',nprecon_ids
     if (prcd_blocksize>nprecon_ids) then
        prcd_blocksize=nprecon_ids
        ! if the blocksize is greater than the # of unknowns on this proc.
        ! reduce the blocksize... then there is only one block on this 
        ! processor
     end if
     nprcd_blocks=floor(real(nprecon_ids,dp)/real(prcd_blocksize,dp))
     if (nprcd_blocks*prcd_blocksize/=nprecon_ids) then
        last_block_size=nprecon_ids-nprcd_blocks*prcd_blocksize
        ! last-block-is-small
        nprcd_blocks=nprcd_blocks+1
     else
        last_block_size=prcd_blocksize
        ! last-block-is-not-small
     end if
     allocate(prcdin(prcd_blocksize,prcd_blocksize,nprcd_blocks))
     allocate(IPIV(prcd_blocksize,nprcd_blocks))
  elseif(near_field_precon) then
  elseif(group_precon) then
  end if
!************************************************
  return
end subroutine precon_init


