program cap2d_layered
!*********************************************************************  
!  A simple & clean serial MOM+EFIE code for surfaces
!  Author: A. E. Yilmaz
!
!  This code is a simplified (teaching) version of the FDAIM_CFIE_SWJV_CKT code
!  that was created by A. E. Yilmaz and H. Bagci in the period 2002-2007,
!  which in turn was based on the codes by A.E.Yilmaz, K.Aygun, B.Shanker,
!  A.Ergin, N.Gres, D.Weile, B.Fisher, E.Michielssen
!
!  Latest Update: October 24, 2007
!  Computational Electromagnetics Group,University of Texas at Austin (CEMG-UT)
!*********************************************************************
  use global_com 
  use global_geom,only:nsuunk,nglunk,nsuinf,edge_av
  use global_dim,only:frequencies,pmatrix,rj
  use global_fast,only:do_MOM_only
  use quadratures,only:determine_quadrature,nqp_t
  use layers,only:is_multilayer,inlayers,init_layers,rs,ro,find_layer,&
  find_height,layer_s,layer_o,fill_Layered_Green,find_rho_z,&
  init_interpolation,green_mode, green_index
  use mat_vec_mult,only:initialize_r0
  implicit none

  integer::i,j
  real(kind=dp)::tim_start,tim_gf,tim_all,tim_dirfield,tim_extra,tim_solve
  integer::Itim_start,Itim_gf,Itim_all,Itim_dirfield,Itim_extra,Itim_solve
  real(kind=dp)::coords(3,2),shift(3,2)
  complex(kind=dp)::go
  integer::resulti,color
  character,allocatable::nameall(:)

  call inmom   ! read in the input parameters for the MOM algorithm
  call determine_quadrature

  ! mesh file for pec surfaces+wires+swjs
  open(unit=12,file='mom.inp',status='old')                 
  ! main input file
  open(unit=66,file='layers.inp',status='old')                 
  ! multilayered media input file
  call parse_layers(66)
  close(66,status='keep')                                    

  call system_clock(COUNT=Itim_start,COUNT_RATE=Itim_rate,COUNT_MAX=Itim_max)
  tim_start=real(Itim_start)/real(Itim_rate)

  call calculate_green_table ! green table is calculated here and stored

  call system_clock(Itim_gf,Itim_rate)
  tim_gf=real(Itim_gf)/real(Itim_rate)
  print*,'TIMING::::::::::GF table',tim_gf-tim_start

  do i=1,2
     open(unit=11,file='geo_pec.inp',status='old')           
     call parse_geom(11)
     close(11,status='keep')                                        

!     call inlayers   ! read in the input parameters for multilayered media

     !print*,'---------------------------------------------------------- '
     !print*,'Serial Capacitance extraction by K. Yang'
     !print*,'---------------------------------------------------------- '
     
     call insu    ! surface discretization info.
     R_a=factor2*edge_av
     print *,'R_a=',R_a
     nsuunk=nsuinf(2)
     nglunk=nsuinf(2)     
     
     if (is_iter) then
        call precon_init
     end if

     call system_clock(Itim_extra,Itim_rate)
     tim_extra=real(Itim_extra)/real(Itim_rate)
     
     print*,'Entering to dirfield'
     call dir_field_mom_only
     print*,'out of dirfield'
     
     call system_clock(Itim_dirfield);tim_dirfield=real(Itim_dirfield)/real(Itim_rate)
     !print*,'TIMING::::::::::Dirfield',tim_dirfield-tim_extra
     !print*,'TIMING::::::::::ALL-BEFORE-SOLVE=',tim_dirfield-tim_start
!*************************************************************
!****  First synchronization point
!*************************************************************
     if (is_iter) then
        call var_arrays     ! allocates variable arrays 
     else
        call invert_pmatrix ! direct solver
     end if
     call system_clock(Itim_solve,Itim_rate)
     tim_solve=real(Itim_solve)/real(Itim_rate)
     !print*,'TIMING::::::::::TOTAL Solve time',tim_solve-tim_dirfield
     
     call cal_C
     call system_clock(Itim_all,Itim_rate)
     tim_all=real(Itim_all)/real(Itim_rate)
     
     call ky_clear_local

     !print*,'TIMING::::::::::TOTAL Solve time',tim_all-tim_solve
     !print*,'TIMING::::::::::TOTAL Cap',tim_all-tim_start         
  end do
  !print*,'Closing files-safely here.'
  stop     
end program cap2d_layered

