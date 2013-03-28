subroutine ky_simulate
  !// LEI: Different with original implementation, I moved the "ky_end_cond_edge" and "ky_end_diel_edge"
  !// function calls into this subroutine to make the subroutine suitable for both stand-alone program
  !// as well as library interface.
  !// For same reason, the "setmom" has been changed to "ky_init" and moved out from this subroutine.
  use global_com
  use global_geom, only:nsuunk, nglunk, nsuinf, edge_av, node_coord, edge_av_cond, edge_av_diel
  use global_dim, only:pmatrix, rj
  use quadratures, only:determine_quadrature, nqp_t
  use mat_vec_mult, only:initialize_r0
  implicit none

#ifdef ENABLED_FOR_MAIN_PROGRAM
  real(kind=dp)::tim_start, tim_gf, tim_all, tim_dirfield, tim_extra
  integer::Itim_start, Itim_gf, Itim_all, Itim_dirfield, Itim_extra
#endif

  !// LEI: calling "ky_end_cond_edge" and "ky_end_diel_edge"
  deallocate(node_coord)
  !// Calculate the maximum and average edge length for conductors.
  !// Find the first edge id for each conductor
  call ky_end_cond_edge

  !// Calculate the maximum and average edge length for dielectric
  call ky_end_diel_edge
  edge_av = (edge_av_cond * nsuinf(1) + edge_av_diel * nsuinf(2)) / (nsuinf(1) + nsuinf(2))
#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'Average Edge length (m): ',edge_av
  print*,'Geometry for conductor and dielectric is generated'
#endif

  !// LEI: Moved from "frequency.f90"
  ! Total number of unknown
  nglunk = nsuinf(1) + nsuinf(2)
#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*, 'nglunk', nglunk
#endif

#ifdef ENABLED_FOR_MAIN_PROGRAM
  call system_clock(COUNT=Itim_start, COUNT_RATE=Itim_rate, COUNT_MAX=Itim_max)
  tim_start = real(Itim_start) / real(Itim_rate)
  print*,'---------------------------------------------------------- '
  print*,'Capacitance extraction by K. Yang'
  print*,'---------------------------------------------------------- '
  write(17,*) '-------------------------------------------------------- '
  write(17,*) 'Capacitance extraction by K. Yang'
  write(17,*) '-------------------------------------------------------- '
#endif

  ! The maximum distance used for singularity extraction
  R_a = max(1.d-7, factor2 * edge_av)
#ifdef ENABLED_FOR_MAIN_PROGRAM
  print *,'R_a=',R_a
  write(17,*) 'Distance for singularity extraction (R_a)=',R_a
#endif

  if (nglunk > 0) then
     ! numerical integration rule (weight and quadrature)
     call determine_quadrature
  end if

  ! is_iter: Iteratively solve matrix equation or invert matrix directly
  if (is_iter) then
     ! Initialize precondition for interative solver
     call precon_init
  end if

#ifdef ENABLED_FOR_MAIN_PROGRAM
  call system_clock(Itim_extra, Itim_rate)
  tim_extra = real(Itim_extra) / real(Itim_rate)
  print*, 'TIMING::::::::::Extra',tim_extra-tim_start
#endif
  
  call dir_field_mom_only ! Fill the interaction matrix

#ifdef ENABLED_FOR_MAIN_PROGRAM
  call system_clock(Itim_dirfield, Itim_rate)
  tim_dirfield = real(Itim_dirfield) / real(Itim_rate)
  print*,'TIMING::::::::::Dirfield',tim_dirfield-tim_extra
  print*,'TIMING::::::::::ALL-BEFORE-SOLVE=',tim_dirfield-tim_start
  write(17,*) 'TIMING::::::::::Dirfield',tim_dirfield-tim_extra
  write(17,*) 'TIMING::::::::::ALL-BEFORE-SOLVE=',tim_dirfield-tim_start
  if (is_iter .NEQV. .TRUE.) then
     open(unit=18, file='pmatrix_orig.txt', status='unknown')
     call debug_print_pmatrix
     close(18, status='keep')
  end if
#endif
  
  !*************************************************************
  !****  First synchronization point
  !*************************************************************
  if (is_iter) then
     call var_arrays		! allocates variable arrays 
     call initialize_r0		!
  else
     call invert_pmatrix	! direct solver
  end if
  call cal_C			! Calculate cap matrix

#ifdef ENABLED_FOR_MAIN_PROGRAM
  call system_clock(Itim_all, Itim_rate)
  tim_all = real(Itim_all) / real(Itim_rate)
  print*,'TIMING::::::::::TOTAL Solve time',tim_all-tim_dirfield
  print*,'TIMING::::::::::TOTAL Cap',tim_all-tim_start
  write(17,*) 'TIMING::::::::::TOTAL Solve time',tim_all-tim_dirfield
  write(17,*) 'TIMING::::::::::TOTAL Cap',tim_all-tim_start
  if (is_iter .NEQV. .TRUE.) then
     open(unit=18, file='pmatrix_invert.txt', status='unknown')
     call debug_print_pmatrix
     close(18, status='keep')
  end if
#endif
  
  return
end subroutine ky_simulate


subroutine ky_getC(idx,jdx,cap)
  use global_com,only:dp,tot_Q
  implicit none

  integer,intent(in)::idx,jdx
  real(kind=dp),intent(out)::cap

  ! return capacitance between two conductors
  cap=tot_Q(idx+1,jdx+1)
#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'The capacitance between conductor',idx+1,'and',jdx+1,'is:',cap
#endif
  return
end subroutine ky_getC


subroutine ky_clear_edge
  use global_com, only:tot_Q, is_iter
  use global_geom, only:edg_cond, cond_sta,edg_diel_epsr, cond_epsr, edg_coord, edg_diel_coord
  use quadratures, only:qp_s, wght_s, qp_t, wght_t
  use global_dim, only:pmatrix, zpast, rj, prcdin
  use mat_vec_mult,only:r0_initial
  implicit none

  ! Clean all the array
  deallocate(edg_cond,cond_sta,edg_diel_epsr,tot_Q)
  deallocate(cond_epsr,edg_coord,edg_diel_coord)
  deallocate(qp_s,wght_s,qp_t,wght_t)
  deallocate(r0_initial)

  if (is_iter) then
     deallocate(zpast,rj,prcdin)
  else
     deallocate(pmatrix)
  end if
  return
end subroutine ky_clear_edge


subroutine invert_pmatrix
  use global_com,only:dp
  use global_dim,only:pmatrix
  use global_geom,only:nglunk

  implicit none
  ! matrix inversion using LU decomposition
  integer,allocatable::ipiv1(:)
  real(kind=dp),allocatable::work(:)
  integer::info(2)

  !*****************************************************************************
  !* START OF INITIALIZATION                                                   *
  !*****************************************************************************
  ! Allocate arrays
  allocate(ipiv1(1:nglunk),work(1:nglunk))
  !* END OF INITIALIZATION                                                     *
  
  !*****************************************************************************
  !* START OF Direct Solver
  !*****************************************************************************
#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'Start to directly solve the matrix'
#endif
  call DGETRF(nglunk,nglunk,pmatrix,nglunk,ipiv1,info(1))
  call DGETRI(nglunk,pmatrix,nglunk,ipiv1,work,nglunk,info(2))
  if (info(1)/=0 .or. info(2)/=0) then
     print*,'inversion of matrix failed',info(1),info(2)
     stop
  else
#ifdef ENABLED_FOR_MAIN_PROGRAM
     print*,'potential matrix has been inverted'
#endif
     deallocate(ipiv1,work)
  end if
  !*****************************************************************************
  !* END OF Direct Solver
  !*****************************************************************************
  return
end subroutine invert_pmatrix


subroutine solve(count,rhsd)
  use global_com,only:dp,maxit,dtol
  use global_dim,only:rj
  use global_geom,only:nglunk
  use mat_vec_mult,only:invert_preconditioner,r0_initial

  implicit none
  integer,intent(in)::count
  real(kind=dp),intent(in)::rhsd(1:nglunk)
  real(kind=dp)::err
  integer::iter

  !*****************************************************************************
  !* START OF INITIALIZATION                                                   *
  !*****************************************************************************
  ! Invert preconditioner only once
  if (count==1) then
     call invert_preconditioner
  end if
  maxit=max(nglunk,1000) ! do at least 1000 iterations before giving up
  !***************************************
  ! Initialize the arrays
  rj(1:nglunk)=0.d0
  !*****************************************************************************
  !* END OF INITIALIZATION                                                     *
  !*****************************************************************************

  !*****************************************************************************
  !* START OF Iterative Solver
  !*****************************************************************************
#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'Start to iteratively solve the matrix for excitation: ',count
#endif
  iter=maxit; err=dtol
  call dtfqmr_serial(nglunk,rhsd(1:nglunk),rj(1:nglunk),err,iter) 
  !*****************************************************************************
  !* END OF Iterative Solver
  !*****************************************************************************
  return
end subroutine solve


subroutine cal_C
  use global_com,only:dp,is_iter,eps0d,tot_Q
  use global_dim,only:pmatrix,rj
  use global_geom,only:nglunk,nsuinf,ncond,cond_sta,cond_epsr
  use misc_dbl,only:cened_dbl
  use mat_vec_mult,only:r0_initial

  implicit none
  integer::dummy,count,idx
  real(kind=dp)::leng
  real(kind=dp)::rhsd(1:nglunk)

  allocate(tot_Q(1:ncond,1:ncond))
  tot_Q(:,:)=0.d0
  if (is_iter) then
    if (ncond==1) then
       rhsd(1:nsuinf(1))=1.d0
       rhsd(nsuinf(1)+1:nglunk)=0.d0
       call solve(1,rhsd)   ! iterative method-of-moments solver
       do count=1,nsuinf(1)
          call cened_dbl(count,leng)
          tot_Q=tot_Q+rj(count)*leng*cond_epsr(count)
       end do
    else
       do count=1,ncond
          rhsd(:)=0.d0
          rhsd(cond_sta(count):cond_sta(count+1)-1)=1.d0
          call solve(count,rhsd)   ! iterative method-of-moments solver
          do dummy=1,ncond
             do idx=cond_sta(dummy),cond_sta(dummy+1)-1
                call cened_dbl(idx,leng)
                tot_Q(dummy,count)=tot_Q(dummy,count)+rj(idx)*leng*cond_epsr(idx)
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
           tot_Q=tot_Q+sum(pmatrix(dummy,1:nsuinf(1)))*leng*cond_epsr(dummy)
        end do
     else
        do dummy=1,ncond
           do idx=cond_sta(dummy),cond_sta(dummy+1)-1
              call cened_dbl(idx,leng)
              do count=1,ncond
                 tot_Q(dummy,count)=tot_Q(dummy,count)+&
                      sum(pmatrix(idx,cond_sta(count):cond_sta(count+1)-1))*leng*cond_epsr(idx)
              end do
           end do
        end do
     end if
  end if
  return
end subroutine cal_C


subroutine ky_ecm_init(use_iter_solver, tolerance, num_quadratures)
  !// LEI: Modified from "setmom" subroutine.
  use global_com
  use quadratures, only:nqp_s, nqp_t, total_maxqp_t
  implicit none
  
  logical,intent(in)::use_iter_solver
  real(kind=dp),intent(in)::tolerance
  integer,intent(in)::num_quadratures
  integer::i, su_order(2), noqp(2)

  pid = 4.0_dp * atan(1.0_dp)	! constant PI
  factor2 = 2.5d0		! factor2
  is_iter = use_iter_solver	! direct or iterative solver
  if (is_iter) then		! tolerance of iterative solver
     if (tolerance < 1.0e-4) then
        dtol = tolerance
     else
        dtol = 1.0d-4
     end if
  end if
  rmu0d=pid*4.d-7; cd=299792458._dp;eps0d=1.d0/(rmu0d*cd**2)

  ! Orders of integrations in the order: source,test
  ! Here always use order (1,1)
  su_order(1:2)=(/4,1/)
  su_order(1) = num_quadratures
  do i=1,2
     if (su_order(i)<=10) then
        noqp(i)=su_order(i)
     else
        noqp(i)=3
     end if
  end do
  nqp_s=noqp(1); nqp_t=noqp(2); 
  total_maxqp_t=nqp_t

#ifdef ENABLED_FOR_MAIN_PROGRAM
  write(17,*) 'factor2:',factor2  
  if (is_iter) then
     write(17,*) 'Iterative solver tolerance', dtol
  end if
  write(17,*) 'Order and # of quadrature points(source,test):',su_order(1),nqp_s,su_order(2),nqp_t
#endif
  return
end subroutine ky_ecm_init


subroutine setmom                                                   
  ! This subroutine inputs algorithmic and electromagnetic information from
  ! mom.inp, and does some of the preprocesseing necessary to construct the 
  ! incident pulse. 
  ! Continually modified as the code expands 
  use global_com 
  use quadratures,only:nqp_s,nqp_t,total_maxqp_t
  implicit none
  integer::i,su_order(2),noqp(2)

  pid=4.0_dp*atan(1.0_dp) ! constant PI
  factor2=2.5d0     ! factor2
  is_iter=.true.   ! direct solver or iterative solve
  if (is_iter) then
     dtol=1.d-4     ! tolerance of iterative solve
  end if
  rmu0d=pid*4.d-7; cd=299792458._dp;eps0d=1.d0/(rmu0d*cd**2)

  ! Orders of integrations in the order: source,test
  ! Here always use order (1,1)
  su_order(1:2)=(/4,1/)
  do i=1,2
     if (su_order(i)<=10) then
        noqp(i)=su_order(i)
     else
        noqp(i)=3
     end if
  end do
  nqp_s=noqp(1); nqp_t=noqp(2); 
  total_maxqp_t=nqp_t
  write(17,*) 'factor2:',factor2  
  if (is_iter) then
     write(17,*) 'Iterative solver tolerance', dtol
  end if
  write(17,*) 'Order and # of quadrature points(source,test):',su_order(1),nqp_s,su_order(2),nqp_t
  return
end subroutine setmom


subroutine precon_init
  use global_com,only:dp
  use global_dim,only:prcdin
  use global_geom,only:nglunk
  implicit none
  !************************************************
  !**************** PRE-CONDITIONER ***************
  !************************************************ 
  ! diagonal precondition
  allocate(prcdin(nglunk))
  return
end subroutine precon_init


subroutine ky_end_cond_edge
  use global_com,only:dp
  use global_geom,only:nsuinf,ncond,cond_sta,edge_av_cond,edg_cond
  use misc_dbl,only:cened_dbl
  implicit none

  integer::j
  integer,allocatable::nedg_cond(:)
  real(kind=dp)::av_leng,max_leng,leng

  ! store the first edge id for each conductor
  allocate(cond_sta(1:ncond+1))
  ! count number of edges for each conductor
  allocate(nedg_cond(1:ncond))
  cond_sta(:)=0
  nedg_cond(:)=0

  ! average and maximum conductor edge length
  av_leng=0.d0;max_leng=0.d0
  do j=1,nsuinf(1)
     nedg_cond(edg_cond(j))=nedg_cond(edg_cond(j))+1
     call cened_dbl(j,leng) ! find edge length for conductor
     av_leng=av_leng+leng
     max_leng=max(max_leng,leng)
  end do
  av_leng=av_leng/nsuinf(1)
  edge_av_cond=av_leng

#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'Maximum and Average Conductor Edge length (m): ',max_leng,av_leng
  write(17,*) 'Maximum and Average Conductor Edge length (m):',max_leng,av_leng
#endif
  
  cond_sta(1)=1
  do j=1,ncond
     cond_sta(j+1)=cond_sta(j)+nedg_cond(j)
  end do
  deallocate(nedg_cond)
  return
end subroutine ky_end_cond_edge


subroutine ky_end_diel_edge
  use global_geom,only:nsuinf,edge_av_diel
  use global_com,only:dp
  use misc_dbl,only:cened_diel_dbl
  implicit none

  integer::j
  real(kind=dp)::av_leng,max_leng,leng

  av_leng=0.d0;max_leng=0.d0
  do j=1,nsuinf(2)
     call cened_diel_dbl(j,leng) ! find edge length for dielectric
     av_leng=av_leng+leng
     max_leng=max(max_leng,leng)
  end do
  av_leng=av_leng/nsuinf(2)
  edge_av_diel=av_leng

#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'Maximum and Average Dielectric Edge length (m):',max_leng,av_leng
  write(17,*) 'Maximum and Average Dielectric Edge length (m):',max_leng,av_leng
#endif
  
  return
end subroutine ky_end_diel_edge


subroutine ky_ecm_num_cond(num_cond)
  use global_com,only:dp
  use global_geom,only:ncond
  implicit none

  integer,intent(in)::num_cond
  ncond = num_cond
  return
end subroutine ky_ecm_num_cond


subroutine ky_ecm_init_cond_edges(num_nodes, num_edges)
  use global_com,only:dp
  use global_geom,only:nsuinf, ncond_node, ncond_edg, cond_epsr, node_coord, edg_coord, edg_cond, &
       add_cond_node_ptr, add_cond_edge_ptr, scale
  implicit none

  integer,intent(in)::num_nodes, num_edges
  ncond_node = num_nodes
  ncond_edg = num_edges
  nsuinf(1) = num_edges
  add_cond_node_ptr = 1
  add_cond_edge_ptr = 1
  allocate(cond_epsr(1:ncond_edg))
  cond_epsr(1:ncond_edg) = 1.d0	! default value is 1 (free space)

  allocate(node_coord(1:2, 1:ncond_node))
  allocate(edg_coord(1:2, 1:2, 1:ncond_edg))
  allocate(edg_cond(1:ncond_edg))
  node_coord(:,:) = 0.d0
  edg_coord(:,:,:) = 0.d0
  edg_cond(:) = 0
  scale = 1.d0 ! Scale of the geometry
  return
end subroutine ky_ecm_init_cond_edges
  

subroutine ky_ecm_init_diel_edges(num_edges)
  use global_com,only:dp
  use global_geom,only:nsuinf, ndiel_edg, edg_diel_coord, edg_diel_epsr, add_diel_edge_ptr
  implicit none

  integer,intent(in)::num_edges
  ndiel_edg = num_edges
  nsuinf(2) = num_edges

  allocate(edg_diel_coord(1:2, 1:2, ndiel_edg))
  allocate(edg_diel_epsr(1:2, 1:ndiel_edg))
  edg_diel_epsr(:,:) = 0.d0
  edg_diel_coord(:,:,:) = 0.d0
  add_diel_edge_ptr = 1

  return
end subroutine ky_ecm_init_diel_edges


subroutine ky_ecm_add_cond_node(xx, yy)
  use global_com,only:dp
  use global_geom,only:node_coord, add_cond_node_ptr
  implicit none
  
  real(kind=dp),intent(in)::xx, yy
  
  node_coord(1, add_cond_node_ptr) = xx
  node_coord(2, add_cond_node_ptr) = yy
  add_cond_node_ptr = add_cond_node_ptr + 1
  return
end subroutine ky_ecm_add_cond_node


subroutine ky_ecm_add_cond_edge(from, to, cid, outer_eps)
  use global_com, only:dp
  use global_geom, only:edg_cond, edg_coord, cond_epsr, node_coord, add_cond_edge_ptr, scale
  implicit none

  integer,intent(in)::from, to, cid
  real(kind=dp),intent(in)::outer_eps

  cond_epsr(add_cond_edge_ptr) = outer_eps
  edg_cond(add_cond_edge_ptr) = cid
  edg_coord(1:2, 1, add_cond_edge_ptr) = node_coord(1:2, from) * scale
  edg_coord(1:2, 2, add_cond_edge_ptr) = node_coord(1:2, to) * scale
  add_cond_edge_ptr = add_cond_edge_ptr + 1

  return
end subroutine ky_ecm_add_cond_edge


subroutine ky_ecm_add_diel_edge(x1, y1, x2, y2, epsr1, epsr2)
  use global_com, only:dp
  use global_geom, only:edg_diel_coord, edg_diel_epsr, scale, add_diel_edge_ptr
  implicit none

  real(kind=dp),intent(in)::x1, y1, x2, y2, epsr1, epsr2

  edg_diel_coord(1, 1, add_diel_edge_ptr) = x1 * scale
  edg_diel_coord(2, 1, add_diel_edge_ptr) = y1 * scale
  edg_diel_coord(1, 2, add_diel_edge_ptr) = x2 * scale
  edg_diel_coord(2, 2, add_diel_edge_ptr) = y2 * scale
  edg_diel_epsr(1, add_diel_edge_ptr) = epsr1
  edg_diel_epsr(2, add_diel_edge_ptr) = epsr2
  add_diel_edge_ptr = add_diel_edge_ptr + 1

  return
end subroutine ky_ecm_add_diel_edge


subroutine debug_print_pmatrix
  use global_dim, only:pmatrix
  use global_geom,only:nglunk
  implicit none

  integer::i, j
  do i=1,nglunk
     do j=1,nglunk
        write(18,*), i, j, pmatrix(i, j), j, i, pmatrix(j, i)
     end do
  end do
end subroutine debug_print_pmatrix
