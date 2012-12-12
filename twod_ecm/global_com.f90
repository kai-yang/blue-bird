module global_com 

!
!  This module contains all the constants which are common to the whole code
!
  implicit none
  save
  integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
! MOM commons 
  real(kind=dp)::R_a,factor2
  real(kind=dp)::pid,cd,eps0d,rmu0d
  real(kind=dp),allocatable::tot_Q(:,:)
! memory counters
  real(kind=dp)::real_mem=8.0_dp,int_mem=4.0_dp
  integer::maxit

  logical::is_iter
  ! FOR MONOSTATIC SCATTERING!
  real(kind=dp)::time_far_allinc=0.0_dp,time_near_allinc=0.0_dp,time_mom_allinc=0.0_dp,&
       time_it_allinc=0.0_dp  
! iterative solver
  real(kind=dp)::dtol

  real(kind=dp)::tot_near_time=0.0_dp,tot_reduction_time=0.0_dp
  integer::max_INF=1e9
! TIMING ROUTINES
  integer::Itim_rate,Itim_max,Itim_global
end module global_com
