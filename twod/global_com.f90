module global_com 

!
!  This module contains all the constants which are common to the whole code
!
  implicit none
  save
  integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
! MOM commons 
  integer::ncond ! number of conductor
  real(kind=dp)::wod,f0d,R_a,factor2
  complex(kind=dp),parameter::c1=(0.0_dp,1.0_dp)

  real(kind=dp)::eta0d
  real(kind=dp)::pid,cd,vlite,eps0d,rmu0d,loss_sigma
  real(kind=dp)::fmax_all,no_freqs
! memory counters
  real(kind=dp)::real_mem=8.0_dp,complex_mem=16.0_dp, int_mem=4.0_dp
  integer::maxit

  real(kind=dp)::lambda_fmax,epss
  logical::closed_body,is_iter,is_pec
  ! FOR MONOSTATIC SCATTERING!
  real(kind=dp)::time_far_allinc=0.0_dp,time_near_allinc=0.0_dp,time_mom_allinc=0.0_dp,&
       time_it_allinc=0.0_dp  
! iterative solver
  real(kind=dp)::dtol

  real(kind=dp)::tot_near_time=0.0_dp,tot_reduction_time=0.0_dp
  integer::istore_cur
  integer::max_INF=1e9
! PRE-CONDITIONER
  real(kind=dp)::precon_dist
  integer::prcd_blocksize,nprcd_blocks,last_block_size,nprecon_ids
  integer,allocatable::my_precon_ids(:),me_to_precon(:)
  logical::diag_precon,block_diag_precon,group_precon,near_field_precon!,precon_all_near,&
! TIMING ROUTINES
  integer::Itim_rate,Itim_max,Itim_global
end module global_com
