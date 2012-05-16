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

subroutine ky_getC(idx,jdx,cap)
  use global_com,only:dp,tot_Q
  implicit none

  integer,intent(in)::idx,jdx
  real(kind=dp),intent(out)::cap

  cap=abs(tot_Q(idx+1,jdx+1)*1.d12)

  print*,'The capacitance between conductor',idx+1,'and',jdx+1,'is:',cap
  return
end subroutine ky_getC

subroutine ky_clear_edge
  use global_geom,only:edg_cond,edg_coord,cond_sta
  use global_com,only:tot_Q,is_iter
  use quadratures,only:qp_s,wght_s,qp_t,wght_t
  use layers,only:z_min,z_max,rho_max,layers_eff,map_layer,num_rho,num_z,&
       drho,dz,gf_table_same,gf_table_diff,nlayers_eff,drho_ij,dz_ij
  use global_dim,only:pmatrix,zpast
  implicit none

  integer::i,j

  deallocate(edg_cond,edg_coord,tot_Q)
  deallocate(cond_sta)
  deallocate(qp_s,wght_s,qp_t,wght_t)
  deallocate(layers_eff,map_layer,z_min,z_max,rho_max,&
       num_rho,num_z,drho,dz,drho_ij,dz_ij)

  do i=1,nlayers_eff
     deallocate(gf_table_same(i)%Gf_grid_array_t,gf_table_same(i)%Gf_grid_array_h)
     do j=1,nlayers_eff
        if (i/=j) then
           deallocate(gf_table_diff(j,i)%Gf_grid_array)
        end if
     end do
  end do
  deallocate(gf_table_same,gf_table_diff)

  if (is_iter) then
     deallocate(zpast)
  else
     deallocate(pmatrix)
  end if
  return
end subroutine ky_clear_edge

subroutine ky_clear_all
  use layers,only:h_of_layer,zlow_of_layer,eps_t,&
       Z_0,GammaL_mn,GammaR_mn,kz_wave,k_prop2
  implicit none

  deallocate(h_of_layer,zlow_of_layer,eps_t,&
       Z_0,GammaL_mn,GammaR_mn,&
       kz_wave,k_prop2)

  return
end subroutine ky_clear_all

subroutine invert_pmatrix
  use global_com 
  use global_dim,only:pmatrix
  use global_geom,only:nglunk,nsuinf
  use mat_vec_mult,only:tot_near_time,tot_reduction_time,&
       invert_preconditioner,matvec,initialize_r0
!                         
  implicit none
  real(kind=dp)::tot_solve_time,mem_est(2)
  ! matrix inversion using LU decomposition
  integer,allocatable::ipiv1(:)
  real(kind=dp),allocatable::work(:)
  integer::info(2),ntotunk

!*****************************************************************************
!* START OF INITIALIZATION                                                   *
!*****************************************************************************
!***************************************
! For the first incidence, allocate arrays
mem_est(2)=mem_est(1)
  ntotunk=nsuinf(1)+nsuinf(2)
  allocate(ipiv1(1:ntotunk),work(1:ntotunk))
!* END OF INITIALIZATION                                                     *
!*****************************************************************************
!*****************************************************************************
!* START OF Iterative Solver
!*****************************************************************************
  print*,'Start to directly solve the matrix'
  call DGETRF(ntotunk,ntotunk,pmatrix,ntotunk,ipiv1,info(1))
  call DGETRI(ntotunk,pmatrix,ntotunk,ipiv1,work,ntotunk,info(2))
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
  use global_geom,only:nglunk,nsuinf
  use mat_vec_mult,only:tot_near_time,tot_reduction_time,&
       invert_preconditioner,matvec,initialize_r0
!                         
  implicit none
  complex(kind=dp),intent(in)::rhsd(1:nsuinf(1)+nsuinf(2))
  real(kind=dp)::err,rerr_local,rerr_sum,bmag_local
  integer::it,me,iter,dummy,i
  integer::k,j,dumb(1),ntotunk
  real(kind=dp)::mem_est(2)

  ntotunk=nsuinf(1)+nsuinf(2)
!*****************************************************************************
!* START OF INITIALIZATION                                                   *
!*****************************************************************************
!***************************************
! For the first incidence, allocate arrays, form and pre-FFT AIM matrices
mem_est(2)=mem_est(1)

call invert_preconditioner
maxit=max(ntotunk,1000) ! do at least 1000 iterations before giving up
!     maxit=10000
  !***************************************
  ! Initialize the arrays
  rj(1:ntotunk)=cmplx(0.0_dp,0.0_dp,dp)
  ! else the initial guess is the solution at the previous frequency
 !***************************************
!*****************************************************************************
!* END OF INITIALIZATION                                                     *
!*****************************************************************************
!*****************************************************************************
!* START OF Iterative Solver
!*****************************************************************************
  print*,'Start to iteratively solve the matrix'
  iter=maxit; err=dtol
  call ztfqmr_serial(ntotunk,rhsd(1:ntotunk),rj(1:ntotunk),err,iter) 
!*****************************************************************************
!* END OF Iterative Solver
!*****************************************************************************
  return
end subroutine solve

subroutine cal_C
  use global_com,only:dp,ncond,is_iter,eps0d,tot_Q
  use global_dim,only:pmatrix,rj
  use global_geom,only:nglunk,nsuinf,npat_cond,cond_sta,cond_epsr
  use misc_dbl,only:cened_dbl
!                         
  implicit none
  integer::dummy,count,idx,ntotunk
  real(kind=dp)::rm(3),leng
  complex(kind=dp)::rhsd(1:nsuinf(1)+nsuinf(2))

  ntotunk=nsuinf(1)+nsuinf(2)
  allocate(tot_Q(1:ncond,1:ncond))
  tot_Q(:,:)=0.d0
  if (is_iter) then
    if (ncond==1) then
       rhsd(1:nsuinf(1))=cmplx(1.0_dp,0.0_dp,dp)
       rhsd(nsuinf(1)+1:ntotunk)=cmplx(0.0_dp,0.0_dp,dp)
       call solve(rhsd)   ! iterative method-of-moments solver
       do count=1,nsuinf(1)
          call cened_dbl(count,leng)
          tot_Q=tot_Q+rj(count)*leng*cond_epsr(count)
       end do
    else
       do dummy=1,ncond
          do count=1,ncond
             rhsd(:)=cmplx(0.0_dp,0.0_dp,dp)
             rhsd(cond_sta(count):cond_sta(count+1)-1)=cmplx(1.0_dp,0.0_dp,dp)
             call solve(rhsd)   ! iterative method-of-moments solver
             do idx=cond_sta(dummy),cond_sta(dummy+1)-1
                call cened_dbl(idx,leng)
                tot_Q(dummy,count)=tot_Q(dummy,count)+rj(idx)*leng*cond_epsr(count)
             end do
          end do
       end do
    end if
  else
     ! Multiply with proper voltage to obtain charge density on each conductor
     ! For 1 conductor, charge density is the summation of cmatrix
     if (ncond==1) then
        do dummy=1,nsuinf(1)
           call cened_dbl(dummy,leng)
           tot_Q=tot_Q+sum(pmatrix(dummy,1:nsuinf(1)))*leng*cond_epsr(count)
        end do
     else
        do dummy=1,ncond
           do idx=cond_sta(dummy),cond_sta(dummy+1)-1
              call cened_dbl(idx,leng)
              do count=1,ncond
                 tot_Q(dummy,count)=tot_Q(dummy,count)+&
                      sum(pmatrix(idx,cond_sta(count):cond_sta(count+1)-1))*leng*cond_epsr(count)
              end do
           end do
        end do
     end if
  end if
  return
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
  real(kind=dp)::dum1,dum2,freq
  integer::i,su_order(2),noqp(2)

  pid=4.0_dp*atan(1.0_dp) 

  factor2=2.5d0     !factor2
  is_iter=.false.     !direct solver or iterative solve
  if (is_iter) then
     dtol=1.d-4     !tolerance of iterative solve
  end if

  rmu0d=pid*4.d-7; cd=299792458._dp;eps0d=1.d0/(rmu0d*cd**2)
  vlite=cd; eta0d=sqrt(rmu0d/eps0d)
  freq=1.d1
  wod=2.d0*pid*freq

  ! Orders of integrations in the order: surface,wire,volume
  su_order(1:2)=(/10,1/)
  i=1; prcd_blocksize=75; precon_dist=0.5d0
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
  write(17,*) 'Order and # of quadrature points(source,test):',su_order(1),nqp_s,su_order(2),nqp_t
  write(17,*) 'The blocksize of the preconditioner is:',prcd_blocksize
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
