subroutine field (me,ne,wghts_phi)
  !     this routine computes the scalar potential and 
  !     compresses them in the wghts array.  this is then used to 
  !     assemble Z_MOM.
  !     last modified : 2011 K. Yang
  use global_com,only:dp,pid
  use global_geom,only:nsuinf
  use misc_dbl,only:cenedg_dbl,cenedg_dmg_dbl
  use quadratures,only:nqp_t,wght_t,total_maxqp_t
  use layers,only:find_layer,layer_me,layer_ne,is_multilayer,k_prop,eps_t,k0

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),intent(out)::wghts_phi
  real(kind=dp),dimension(2)::rm,rn
  real(kind=dp),dimension(2)::rm_g
  integer::i,patch_pos,ne_dmg
  complex(kind=dp)::wghts_phi_direct,wghts_phi_ps,wghts_phi_ns,wghts_phi_tmp

  wghts_phi=0.d0; wghts_phi_tmp=cmplx(0.d0,0.d0,dp)
  !------------------------------------------------------------------------
  !     testing on a conductor surface element:
  !------------------------------------------------------------------------
  call cenedg_dbl(me,rm)     
  call find_layer(rm(:),layer_me)
  rm_g(1:2)=rm
  !       source 
  !       if surface ++++++++++++++++++++++++++++++++++++++++
  source_s:if (ne<=nsuinf(1)) then
     call cenedg_dbl(ne,rn)
  else
     ne_dmg=ne-nsuinf(1)
     call cenedg_dmg_dbl(ne_dmg,rn)
  end if source_s

  call find_layer(rn(:),layer_ne)
  k_prop=k0*sqrt(eps_t(layer_ne))

  if (layer_ne==layer_me) then
     patch_pos=1
  else if (layer_ne<layer_me) then
     patch_pos=2
  else
     patch_pos=3
  end if
  
  select case (patch_pos)
  case (1) ! 5 terms to be seperated to calculate analytically or numerically
     call find_direct_intg(ne,me,rn,rm,rm_g,wghts_phi_direct)
     wghts_phi_tmp=wghts_phi_tmp+wghts_phi_direct
!     print*,'dir',wghts_phi_direct
           
     if (is_multilayer) then
        call find_ps_intg(ne,me,rn,rm,rm_g,wghts_phi_ps) ! potentially singular term 
        wghts_phi_tmp=wghts_phi_tmp+wghts_phi_ps
!        print*,'ps',wghts_phi_ps
!        stop
              
        call find_ns_intg(ne,me,rn,rm,rm_g,wghts_phi_ns) ! non-singular term
        wghts_phi_tmp=wghts_phi_tmp+wghts_phi_ns
!        print*,'ns',wghts_phi_ns
!        stop
     end if
  case (2)
     if (is_multilayer) then
        call find_ps_intg2(ne,me,rn,rm,rm_g,wghts_phi_ps) ! potentially singular term 
        wghts_phi_tmp=wghts_phi_tmp+wghts_phi_ps
!        print*,'ps',wghts_phi_ps
!        stop
              
        call find_ns_intg2(ne,me,rn,rm,rm_g,wghts_phi_ns) ! non-singular term
        wghts_phi_tmp=wghts_phi_tmp+wghts_phi_ns
!        print*,'ns',wghts_phi_ns
!        stop
     end if
  case (3)
     if (is_multilayer) then
        call find_ps_intg3(ne,me,rn,rm,rm_g,wghts_phi_ps) ! potentially singular term 
        wghts_phi_tmp=wghts_phi_tmp+wghts_phi_ps
!        print*,'ps',wghts_phi_ps
!        stop
              
        call find_ns_intg3(ne,me,rn,rm,rm_g,wghts_phi_ns) ! non-singular term
        wghts_phi_tmp=wghts_phi_tmp+wghts_phi_ns
!        print*,'ns',wghts_phi_ns
!        stop
     end if
  end select
!  print*,patch_pos!,wghts_phi_tmp
  wghts_phi=real(wghts_phi_tmp,dp)

  return 
end subroutine field

subroutine find_direct_intg(ne,me,rn,rm,rm_g,wghts_phi_direct)
  use global_com,only:dp,r_a
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_direct
  real(kind=dp)::dist,rnear

  wghts_phi_direct=0.d0
  rnear=r_a  ! distance for singularity extractions
  dist=sqrt(sum((rm(:)-rn(:))**2))
  ! near far is based on the center-to-center distances in all cases
  if (dist<=rnear) then
!     print*,'dir near'
     call source_surf_near(ne,me,rm_g(1:2),wghts_phi_direct)
  else
!     print*,'dir far'
     call source_surf_far(ne,me,rm_g(1:2),wghts_phi_direct)
  end if
  return
end subroutine find_direct_intg

subroutine find_ps_intg(ne,me,rn,rm,rm_g,wghts_phi_ps)
  use global_com,only:dp,r_a
  use layers,only:layer_me,layer_ne,zlow_of_layer,h_of_layer,&
       GammaR_mn,GammaL_mn
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_ps
  integer::i,ps
  real(kind=dp)::dist,rnear,dist_ps(4)
  real(kind=dp)::coeff(4),rm_g_tmp(2,4),rm_tmp(2,4)
  complex(kind=dp)::phimn

  wghts_phi_ps=cmplx(0.d0,0.d0,dp)
  rnear=r_a  ! distance for singularity extractions

  ! For m=n, 4 terms to be singular, otherwise, 9 terms to be evaluated
  ! Whether it's necessary to evaluate these terms
  dist_ps(1)=zlow_of_layer(layer_ne+1)
  dist_ps(2)=zlow_of_layer(layer_ne)
  dist_ps(3)=h_of_layer(layer_ne)
  dist_ps(4)=dist_ps(3)

  ! testing qp transformation
  do i=1,4
     rm_g_tmp(1:2,i)=rm_g(1:2) ! Do the transformation wrt z
  end do
  rm_g_tmp(2,1)=2.d0*zlow_of_layer(layer_ne+1)-rm_g(2)
  rm_g_tmp(2,2)=2.d0*zlow_of_layer(layer_ne)-rm_g(2)
  rm_g_tmp(2,3)=rm_g(2)+2.d0*h_of_layer(layer_ne)
  rm_g_tmp(2,4)=rm_g(2)-2.d0*h_of_layer(layer_ne)

  ! center point of testing patch transformation
  do i=1,4
     rm_tmp(1:2,i)=rm(1:2) ! Do the transformation wrt z
  end do
  rm_tmp(2,1)=2.d0*zlow_of_layer(layer_ne+1)-rm(2)
  rm_tmp(2,2)=2.d0*zlow_of_layer(layer_ne)-rm(2)
  rm_tmp(2,3)=rm(2)+2.d0*h_of_layer(layer_ne)
  rm_tmp(2,4)=rm(2)-2.d0*h_of_layer(layer_ne)

  ! For Gq,ps
  coeff(1)=GammaR_mn(layer_ne)
  coeff(2)=GammaL_mn(layer_ne)
  coeff(3)=GammaR_mn(layer_ne)*GammaL_mn(layer_ne)
  coeff(4)=coeff(3)

  do ps=1,4
     if (dist_ps(ps)/=-1.d0) then          
        dist=sqrt(sum((rm_tmp(:,ps)-rn(:))**2))
        ! near far is based on the center-to-center distances in all cases
        if (dist<=rnear) then
!           print*,'ps near'
           call source_surf_near(ne,me,rm_g_tmp(1:2,ps),phimn)
        else
!           print*,'ps far'
           call source_surf_far(ne,me,rm_g_tmp(1:2,ps),phimn)
        end if
        wghts_phi_ps=wghts_phi_ps+coeff(ps)*phimn
     end if
  end do
  return
end subroutine find_ps_intg

subroutine find_ps_intg2(ne,me,rn,rm,rm_g,wghts_phi_ps) ! o>s
  use global_com,only:dp,r_a
  use layers,only:layer_me,layer_ne,zlow_of_layer,h_of_layer,&
       GammaR_mn,GammaL_mn
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_ps
  integer::i,ps
  real(kind=dp)::dist,rnear,dist_ps(10)
  real(kind=dp)::coeff(10),rm_g_tmp(2,10),rm_tmp(2,10)
  real(kind=dp)::h_tot,TauR_vv_prod_coef
  complex(kind=dp)::phimn

  wghts_phi_ps=cmplx(0.d0,0.d0,dp)
  rnear=r_a  ! distance for singularity extractions

  TauR_vv_prod_coef=1.d0
  h_tot=0.d0
  do i=layer_ne+1,layer_me-1
     h_tot=h_tot+h_of_layer(i)
     TauR_vv_prod_coef=TauR_vv_prod_coef*(1.d0+GammaR_mn(i))
  end do

  ! For m/=n, 10 terms might be singular
  ! dist_ps(1),(2) are always finite
  dist_ps(1)=0.d0
  dist_ps(2)=0.d0
  dist_ps(3)=zlow_of_layer(layer_ne)
  dist_ps(4)=h_of_layer(layer_ne)
  dist_ps(5)=dist_ps(4)
  dist_ps(6)=zlow_of_layer(layer_me+1)
  dist_ps(7)=zlow_of_layer(layer_me+1)
  if (h_of_layer(layer_ne)==-1.d0 .or. zlow_of_layer(layer_me+1)==-1.d0) then
     dist_ps(8)=-1.d0
  else
     dist_ps(8)=0.d0 ! this term need to be calculated
  end if
  dist_ps(9)=dist_ps(8)
  dist_ps(10)=dist_ps(8)

  ! testing qp transformation
  do i=1,10
     rm_g_tmp(1:2,i)=rm_g(1:2) ! Do the transformation wrt z
  end do
  rm_g_tmp(2,1)=zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm_g(2)
  rm_g_tmp(2,2)=zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm_g(2)
  rm_g_tmp(2,3)=-(zlow_of_layer(layer_ne+1)-2.d0*zlow_of_layer(layer_ne)-zlow_of_layer(layer_me)+h_tot+rm_g(2))
  rm_g_tmp(2,4)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm_g(2)
  rm_g_tmp(2,5)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm_g(2))
  rm_g_tmp(2,6)=zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)-zlow_of_layer(layer_me)+h_tot-rm_g(2)
  rm_g_tmp(2,7)=zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)-zlow_of_layer(layer_me)+h_tot-rm_g(2)
  rm_g_tmp(2,8)=-(zlow_of_layer(layer_ne+1)-2.d0*zlow_of_layer(layer_ne)+2.d0*zlow_of_layer(layer_me+1)&
       -zlow_of_layer(layer_me)+h_tot-rm_g(2))
  rm_g_tmp(2,9)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)&
       -zlow_of_layer(layer_me)+h_tot-rm_g(2)
  rm_g_tmp(2,10)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)&
       -zlow_of_layer(layer_me)+h_tot-rm_g(2))
  
  ! center point of testing patch transformation
  do i=1,10
     rm_tmp(1:2,i)=rm(1:2) ! Do the transformation wrt z
  end do
  rm_tmp(2,1)=zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm(2)
  rm_tmp(2,2)=zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm(2)
  rm_tmp(2,3)=-(zlow_of_layer(layer_ne+1)-2.d0*zlow_of_layer(layer_ne)-zlow_of_layer(layer_me)+h_tot+rm(2))
  rm_tmp(2,4)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm(2)
  rm_tmp(2,5)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_me)+h_tot+rm(2))
  rm_tmp(2,6)=zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)-zlow_of_layer(layer_me)+h_tot-rm(2)
  rm_tmp(2,7)=zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)-zlow_of_layer(layer_me)+h_tot-rm(2)
  rm_tmp(2,8)=-(zlow_of_layer(layer_ne+1)-2.d0*zlow_of_layer(layer_ne)+2.d0*zlow_of_layer(layer_me+1)&
       -zlow_of_layer(layer_me)+h_tot-rm(2))
  rm_tmp(2,9)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)&
       -zlow_of_layer(layer_me)+h_tot-rm(2)
  rm_tmp(2,10)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne+1)+2.d0*zlow_of_layer(layer_me+1)&
       -zlow_of_layer(layer_me)+h_tot-rm(2))
  
  ! For Gq,ps
  coeff(1)=1.d0
  coeff(2)=GammaR_mn(layer_ne)
  coeff(3)=GammaL_mn(layer_ne)
  coeff(4)=GammaR_mn(layer_ne)*GammaL_mn(layer_ne)
  coeff(5)=coeff(4)
  coeff(6)=GammaR_mn(layer_me)
  coeff(7)=GammaR_mn(layer_ne)*GammaR_mn(layer_me)
  coeff(8)=GammaL_mn(layer_ne)*GammaR_mn(layer_me)
  coeff(9)=GammaR_mn(layer_ne)*GammaL_mn(layer_ne)*GammaR_mn(layer_me)
  coeff(10)=coeff(9)
  coeff(1:10)=coeff(1:10)*TauR_vv_prod_coef

  do ps=1,10
     if (dist_ps(ps)/=-1.d0) then          
        dist=sqrt(sum((rm_tmp(:,ps)-rn(:))**2))
        ! near far is based on the center-to-center distances in all cases
        if (dist<=rnear) then
!           print*,'ps near'
           call source_surf_near(ne,me,rm_g_tmp(1:2,ps),phimn)
        else
!           print*,'ps far'
           call source_surf_far(ne,me,rm_g_tmp(1:2,ps),phimn)
        end if
        wghts_phi_ps=wghts_phi_ps+coeff(ps)*phimn
     end if
  end do
  return
end subroutine find_ps_intg2

subroutine find_ps_intg3(ne,me,rn,rm,rm_g,wghts_phi_ps) ! o<s
  use global_com,only:dp,r_a
  use layers,only:layer_me,layer_ne,zlow_of_layer,h_of_layer,&
       GammaR_mn,GammaL_mn
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_ps
  integer::i,ps
  real(kind=dp)::dist,rnear,dist_ps(10)
  real(kind=dp)::coeff(10),rm_g_tmp(2,10),rm_tmp(2,10)
  real(kind=dp)::h_tot,TauL_vv_prod_coef
  complex(kind=dp)::phimn

  wghts_phi_ps=cmplx(0.d0,0.d0,dp)
  rnear=r_a  ! distance for singularity extractions

  TauL_vv_prod_coef=1.d0
  h_tot=0.d0
  do i=layer_me+1,layer_ne-1
     h_tot=h_tot+h_of_layer(i)
     TauL_vv_prod_coef=TauL_vv_prod_coef*(1.d0+GammaL_mn(i))
  end do

  ! For m/=n, 10 terms might be singular
  dist_ps(1)=0.d0
  dist_ps(2)=zlow_of_layer(layer_ne+1)
  dist_ps(3)=0.d0
  dist_ps(4)=h_of_layer(layer_ne)
  dist_ps(5)=dist_ps(4)
  dist_ps(6)=zlow_of_layer(layer_me)
  if (h_of_layer(layer_ne)==-1.d0 .or. zlow_of_layer(layer_me)==-1.d0) then
     dist_ps(7)=-1.d0
  else
     dist_ps(7)=0.d0 ! this term need to be calculated
  end if
  dist_ps(8)=zlow_of_layer(layer_me)
  dist_ps(9)=dist_ps(7)
  dist_ps(10)=dist_ps(7)

  ! testing qp transformation
  do i=1,10
     rm_g_tmp(1:2,i)=rm_g(1:2) ! Do the transformation wrt z
  end do
  rm_g_tmp(2,1)=-(-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm_g(2))
  rm_g_tmp(2,2)=2.d0*zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm_g(2)
  rm_g_tmp(2,3)=-(zlow_of_layer(layer_ne)-2.d0*zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm_g(2))
  rm_g_tmp(2,4)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm_g(2)
  rm_g_tmp(2,5)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm_g(2))
  rm_g_tmp(2,6)=-(-zlow_of_layer(layer_ne)-2.d0*zlow_of_layer(layer_me)+zlow_of_layer(layer_me+1)+h_tot+rm_g(2))
  rm_g_tmp(2,7)=2.d0*zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm_g(2)
  rm_g_tmp(2,8)=-(zlow_of_layer(layer_ne)-2.d0*zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm_g(2))
  rm_g_tmp(2,9)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm_g(2)
  rm_g_tmp(2,10)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm_g(2))
  
  ! center point of testing patch transformation
  do i=1,10
     rm_tmp(1:2,i)=rm(1:2) ! Do the transformation wrt z
  end do
  rm_tmp(2,1)=-(-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm(2))
  rm_tmp(2,2)=2.d0*zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm(2)
  rm_tmp(2,3)=-(zlow_of_layer(layer_ne)-2.d0*zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm(2))
  rm_tmp(2,4)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm(2)
  rm_tmp(2,5)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)+h_tot-rm(2))
  rm_tmp(2,6)=-(-zlow_of_layer(layer_ne)-2.d0*zlow_of_layer(layer_me)+zlow_of_layer(layer_me+1)+h_tot+rm(2))
  rm_tmp(2,7)=2.d0*zlow_of_layer(layer_ne+1)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm(2)
  rm_tmp(2,8)=-(zlow_of_layer(layer_ne)-2.d0*zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm(2))
  rm_tmp(2,9)=2.d0*h_of_layer(layer_ne)+zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm(2)
  rm_tmp(2,10)=-(2.d0*h_of_layer(layer_ne)-zlow_of_layer(layer_ne)+zlow_of_layer(layer_me+1)&
       -2.d0*zlow_of_layer(layer_me)+h_tot+rm(2))
  
  ! For Gq,ps
  coeff(1)=1.d0
  coeff(2)=GammaR_mn(layer_ne)
  coeff(3)=GammaL_mn(layer_ne)
  coeff(4)=GammaR_mn(layer_ne)*GammaL_mn(layer_ne)
  coeff(5)=coeff(4)
  coeff(6)=GammaL_mn(layer_me)
  coeff(7)=GammaR_mn(layer_ne)*GammaL_mn(layer_me)
  coeff(8)=GammaL_mn(layer_ne)*GammaL_mn(layer_me)
  coeff(9)=GammaR_mn(layer_ne)*GammaL_mn(layer_ne)*GammaL_mn(layer_me)
  coeff(10)=coeff(9)
  coeff(1:10)=coeff(1:10)*TauL_vv_prod_coef

  do ps=1,10
     if (dist_ps(ps)/=-1.d0) then          
        dist=sqrt(sum((rm_tmp(:,ps)-rn(:))**2))
        ! near far is based on the center-to-center distances in all cases
        if (dist<=rnear) then
!           print*,'ps near'
           call source_surf_near(ne,me,rm_g_tmp(1:2,ps),phimn)
        else
!           print*,'ps far'
           call source_surf_far(ne,me,rm_g_tmp(1:2,ps),phimn)
        end if
        wghts_phi_ps=wghts_phi_ps+coeff(ps)*phimn
     end if
  end do
  return
end subroutine find_ps_intg3

subroutine find_ns_intg(ne,me,rn,rm,rm_g,wghts_phi_ns)
  use global_com,only:dp,r_a
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_ns
  real(kind=dp)::dist,rnear

  wghts_phi_ns=cmplx(0.d0,0.d0,dp)
  rnear=r_a  ! distance for singularity extractions
  dist=sqrt(sum((rm(:)-rn(:))**2))
  ! near far is based on the center-to-center distances in all cases
  if (dist<=rnear) then
!     print*,'ns near'
     call source_surf_ns_near(ne,me,rm_g(1:2),wghts_phi_ns)
  else
!     print*,'ns far'
     call source_surf_ns_far(ne,me,rm_g(1:2),wghts_phi_ns)
  end if
  return
end subroutine find_ns_intg

subroutine find_ns_intg2(ne,me,rn,rm,rm_g,wghts_phi_ns)
  use global_com,only:dp,r_a
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_ns
  real(kind=dp)::dist,rnear

  wghts_phi_ns=cmplx(0.d0,0.d0,dp)
  rnear=r_a  ! distance for singularity extractions
  dist=sqrt(sum((rm(:)-rn(:))**2))
  ! near far is based on the center-to-center distances in all cases
  if (dist<=rnear) then
!     print*,'ns near'
     call source_surf_ns_near2(ne,me,rm_g(1:2),wghts_phi_ns)
  else
!     print*,'ns far'
     call source_surf_ns_far2(ne,me,rm_g(1:2),wghts_phi_ns)
  end if
  return
end subroutine find_ns_intg2

subroutine find_ns_intg3(ne,me,rn,rm,rm_g,wghts_phi_ns)
  use global_com,only:dp,r_a
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  complex(kind=dp),intent(out)::wghts_phi_ns
  real(kind=dp)::dist,rnear

  wghts_phi_ns=cmplx(0.d0,0.d0,dp)
  rnear=r_a  ! distance for singularity extractions
  dist=sqrt(sum((rm(:)-rn(:))**2))
  ! near far is based on the center-to-center distances in all cases
  if (dist<=rnear) then
!     print*,'ns near'
     call source_surf_ns_near3(ne,me,rm_g(1:2),wghts_phi_ns)
  else
!     print*,'ns far'
     call source_surf_ns_far3(ne,me,rm_g(1:2),wghts_phi_ns)
  end if
  return
end subroutine find_ns_intg3
!==========================================================================
!==========================================================================
! SOURCE SURFACE (FAR)
!==========================================================================
!==========================================================================
subroutine source_surf_far(ne,me,rm,wghts_phi)
  ! ne: basis ID
  ! me: testing ID
  ! rm: q.p. coordinates
  ! wghts: vector potential,scalar potential, and H field contributions to Z-matrix
  ! source is a surface function
  ! testing is a surface function (can be a wire function if EFIE_only)
  use global_com,only:dp,eps0d,pid,c1
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers,only:layer_ne,eps_t,k_prop,euler

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  real(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_dmg
  real(kind=dp)::aaint,distances(nqp_s),R(2,nqp_s)
  complex(kind=dp)::const

  wghts_phi=cmplx(0.d0,0.d0,dp)
  const=-0.25d0*c1-(euler+dlog(0.5d0*k_prop))/(2.d0*pid)

  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if

  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
  end do
    
  do source=1,nqp_s
     !EFIE
     R(:,source)=rm(1:2)-rsrc(1:2,source)
     distances(source)=sqrt(dot_product(R(:,source),R(:,source)))!R
     aaint=wght_s(source)*log(distances(source))
     phimn(source)=aaint
  end do
  wghts_phi=wghts_phi-sum(phimn(1:nqp_s))
  wghts_phi=wghts_phi*leng_s/(2.d0*pid*eps0d*eps_t(layer_ne))
  wghts_phi=wghts_phi+const*leng_s/(eps0d*eps_t(layer_ne))
  return
end subroutine source_surf_far
!==========================================================================
!==========================================================================
! SOURCE SURFACE (NEAR)
!==========================================================================
!==========================================================================
subroutine source_surf_near(ne,me,rm,wghts_phi)
  use global_com,only:dp,eps0d,pid,c1
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers,only:layer_ne,eps_t,euler,k_prop

  implicit none                      
  !
  !depending on the testing scheme one_or_seven is 1 or 7
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  real(kind=dp)::v1(2),v2(2)
  real(kind=dp)::aaint1,aaint2(2)
  integer::test,ne_dmg
  real(kind=dp)::phimn1,leng_s
  complex(kind=dp)::const

  wghts_phi=cmplx(0.d0,0.d0,dp); aaint1=0.d0; aaint2=0.d0
  const=-0.25d0*c1-(euler+dlog(0.5d0*k_prop))/(2.d0*pid)

  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if
  ! analytical e- and h-field integrals
  call sintg_2d(rm(1:2),v1,v2,aaint1,aaint2)
  phimn1=aaint1
  wghts_phi=wghts_phi-phimn1
  wghts_phi=wghts_phi/(2.d0*pid*eps0d*eps_t(layer_ne))
  wghts_phi=wghts_phi+const*leng_s/(eps0d*eps_t(layer_ne))
  return
end subroutine source_surf_near

subroutine source_surf_ns_far(ne,me,rm,wghts_phi)
  ! ne: basis ID
  ! me: testing ID
  ! rm: q.p. coordinates
  ! wghts: vector potential,scalar potential, and H field contributions to Z-matrix
  ! source is a surface function
  ! testing is a surface function (can be a wire function if EFIE_only)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  complex(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_dmg
  complex(kind=dp)::Gf_nsigu,Gf_t_nsigu(3),Gf_h_nsigu(3)
  complex(kind=dp)::Gf,Gf_tmp,Gf_tmp1,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5
  complex(kind=dp)::Gf_ref_t,Gf_ref_h

  wghts_phi=cmplx(0.d0,0.d0,dp)
  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if

  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
  end do
    
  ro(:)=rm(1:2)    
  call find_layer(ro(:),layer_o)

  do source=1,nqp_s
     !EFIE
     rs(:)=rsrc(1:2,source)
     call find_layer(rs(:),layer_s)
     call find_height

     ! Direct calculation
!     call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!     Gf_t_nsigu=Gf_tmp2
!     Gf_h_nsigu=Gf_tmp3
!     Gf_nsigu=Gf_t_nsigu(1)+Gf_h_nsigu(1)
!     print*,'fdir',ne,me,Gf_nsigu!,Gf_t_nsigu(1),Gf_h_nsigu(1)
!!$
     ! Interpolation
     call Gf_interpolation_2d(rs,ro,Gf_t_nsigu(1),Gf_h_nsigu(1))
     Gf_nsigu=Gf_t_nsigu(1)+Gf_h_nsigu(1)
!     print*,'inter',ne,me,Gf_nsigu!,Gf_t_nsigu(1),Gf_h_nsigu(1)
     phimn(source)=wght_s(source)*Gf_nsigu
  end do
  wghts_phi=wghts_phi+sum(phimn(1:nqp_s)) ! '-' sign has been included in Gf
  wghts_phi=wghts_phi*leng_s/eps0d
  return
end subroutine source_surf_ns_far

subroutine source_surf_ns_near(ne,me,rm,wghts_phi)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers

  implicit none                      
  !
  !depending on the testing scheme one_or_seven is 1 or 7
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  complex(kind=dp),dimension(3*nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,3*nqp_s),rhon(2,3*nqp_s),sub_area(3*nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,n_src,ne_dmg
  complex(kind=dp)::Gf_nsigu,Gf_t_nsigu(3),Gf_h_nsigu(3)
  complex(kind=dp)::Gf,Gf_tmp,Gf_tmp1,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5
  complex(kind=dp)::Gf_ref_t,Gf_ref_h
  logical::divide_surf

  wghts_phi=0.d0
  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if
  
  n_src=nqp_s
  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
     sub_area(source)=wght_s(source)
  end do

  ro(:)=rm(1:2)    
  call find_layer(ro(:),layer_o)
  call how2divide_edge

  do source=1,n_src
     !EFIE
     rs(:)=rsrc(1:2,source)
     call find_layer(rs(:),layer_s)
     call find_height

     ! Direct calculation
!     call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!     Gf_t_nsigu(1)=Gf_tmp2
!     Gf_h_nsigu(1)=Gf_tmp3
!     Gf_nsigu=Gf_t_nsigu(1)+Gf_h_nsigu(1)
!     print*,'ndir',ne,me,Gf_nsigu!,Gf_t_nsigu,Gf_h_nsigu

     ! Interpolation
     call Gf_interpolation_2d(rs,ro,Gf_t_nsigu(1),Gf_h_nsigu(1))
     Gf_nsigu=Gf_t_nsigu(1)+Gf_h_nsigu(1)
!     print*,'inter',ne,me,Gf_nsigu!,Gf_t_nsigu(1),Gf_h_nsigu(1)

     phimn(source)=sub_area(source)*Gf_nsigu
  end do

  wghts_phi=wghts_phi+sum(phimn(1:n_src))
  wghts_phi=wghts_phi*leng_s/eps0d
  return
contains
  subroutine how2divide_edge
    ! decides if we will divide the source triangle into 3
    ! and if so,which point to use as the 4th point in the dielectric case
    ! i.e.,the divider point is set
    real(kind=dp)::d,vo1(2),pl(2)
    real(kind=dp),dimension(2)::divider
    integer::inside,test
    if (me<=nsuinf(1)) then !divide the self term for surface testing
       if(ne==me) then
          divide_surf=.true.
          divider(1:2)=rm(1:2)
       else
          divide_surf=.false.
       end if
    end if
    if (divide_surf) then
       call divide_edge_due2pointinside(divider(1:2))
       n_src=2*nqp_s
    end if
    return
  end subroutine how2divide_edge

  subroutine divide_edge_due2pointinside(divider)
    real(kind=dp),intent(in)::divider(2)
    real(kind=dp)::leng_sub,coeff

    leng_sub=sqrt(sum((divider(:)-v1(:))**2))
    coeff=leng_sub/leng_s
    do source=1,nqp_s
       rsrc(:,source)=(divider(:)-v1(:))*qp_s(source)+v1(:)
       sub_area(source)=coeff*wght_s(source)
    end do
    
    leng_sub=sqrt(sum((v2(:)-divider(:))**2))
    coeff=leng_sub/leng_s
    do source=1,nqp_s
       rsrc(:,source+nqp_s)=(v2(:)-divider(:))*qp_s(source)+divider(:)
       sub_area(source+nqp_s)=coeff*wght_s(source)
    end do
    return
  end subroutine divide_edge_due2pointinside
end subroutine source_surf_ns_near

subroutine source_surf_ns_far2(ne,me,rm,wghts_phi)
  ! ne: basis ID
  ! me: testing ID
  ! rm: q.p. coordinates
  ! wghts: vector potential,scalar potential, and H field contributions to Z-matrix
  ! source is a surface function
  ! testing is a surface function (can be a wire function if EFIE_only)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  complex(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_dmg
  complex(kind=dp)::Gf_nsigu(3)
  complex(kind=dp)::Gf,Gf_tmp,Gf_tmp1,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5

  wghts_phi=cmplx(0.d0,0.d0,dp)

  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if

  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
  end do
    
  ro(:)=rm(1:2)    
  call find_layer(ro(:),layer_o)

  do source=1,nqp_s
     !EFIE
     rs(:)=rsrc(1:2,source)
     call find_layer(rs(:),layer_s)
     call find_height

     ! Direct calculation
!     call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!     Gf_nsigu(1)=Gf_tmp
!     print*,'dir',ne,me,Gf_nsigu(1)

     ! Interpolation
     call Gf_interpolation_3d(rs,ro,Gf_nsigu(1))
!     print*,'inter',ne,me,Gf_nsigu(1)
     phimn(source)=wght_s(source)*Gf_nsigu(1)
  end do
  wghts_phi=wghts_phi+sum(phimn(1:nqp_s)) ! '-' sign has been included in Gf
  wghts_phi=wghts_phi*leng_s/eps0d
  return
end subroutine source_surf_ns_far2

subroutine source_surf_ns_near2(ne,me,rm,wghts_phi)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers

  implicit none                      
  !
  !depending on the testing scheme one_or_seven is 1 or 7
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  complex(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s),sub_area(nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_dmg
  complex(kind=dp)::Gf_nsigu(3)
  complex(kind=dp)::Gf,Gf_tmp,Gf_tmp1,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5
  logical::divide_surf

  wghts_phi=cmplx(0.d0,0.d0,dp)
  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if
  
  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
     sub_area(source)=wght_s(source)
  end do

  ro(:)=rm(1:2)    
  call find_layer(ro(:),layer_o)

  do source=1,nqp_s
     !EFIE
     rs(:)=rsrc(1:2,source)
     call find_layer(rs(:),layer_s)
     call find_height

     ! Direct calculation
!     call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!     Gf_nsigu(1)=Gf_tmp
!     print*,'dir',ne,me,Gf_nsigu(1)

     ! Interpolation
     call Gf_interpolation_3d(rs,ro,Gf_nsigu(1))
!     print*,'inter',ne,me,Gf_nsigu(1)
     phimn(source)=sub_area(source)*Gf_nsigu(1)
  end do

  wghts_phi=wghts_phi+sum(phimn(1:nqp_s))
  wghts_phi=wghts_phi*leng_s/eps0d
  return
end subroutine source_surf_ns_near2

subroutine source_surf_ns_far3(ne,me,rm,wghts_phi)
  ! ne: basis ID
  ! me: testing ID
  ! rm: q.p. coordinates
  ! wghts: vector potential,scalar potential, and H field contributions to Z-matrix
  ! source is a surface function
  ! testing is a surface function (can be a wire function if EFIE_only)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  complex(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_dmg
  complex(kind=dp)::Gf_nsigu(3)
  complex(kind=dp)::Gf,Gf_tmp,Gf_tmp1,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5

  wghts_phi=cmplx(0.d0,0.d0,dp)
  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if

  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
  end do
    
  ro(:)=rm(1:2)    
  call find_layer(ro(:),layer_o)

  do source=1,nqp_s
     !EFIE
     rs(:)=rsrc(1:2,source)
     call find_layer(rs(:),layer_s)
     call find_height

     ! Direct calculation
!     call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!     Gf_nsigu(1)=Gf_tmp
!     print*,'dir',ne,me,Gf_nsigu(1)

     ! Interpolation
     call Gf_interpolation_3d(rs,ro,Gf_nsigu(1))
!     print*,'inter',ne,me,Gf_nsigu(1)
     phimn(source)=wght_s(source)*Gf_nsigu(1)
  end do
  wghts_phi=wghts_phi+sum(phimn(1:nqp_s)) ! '-' sign has been included in Gf
  wghts_phi=wghts_phi*leng_s/eps0d
  return
end subroutine source_surf_ns_far3

subroutine source_surf_ns_near3(ne,me,rm,wghts_phi)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_dmg_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_dmg_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  use layers

  implicit none                      
  !
  !depending on the testing scheme one_or_seven is 1 or 7
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  complex(kind=dp),intent(out)::wghts_phi

  complex(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s),sub_area(nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_dmg
  complex(kind=dp)::Gf_nsigu(3)
  complex(kind=dp)::Gf,Gf_tmp,Gf_tmp1,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5
  logical::divide_surf

  wghts_phi=cmplx(0.d0,0.d0,dp)
  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_dmg=ne-nsuinf(1)
     v1(:)=edg_dmg_coord(:,1,ne_dmg)
     v2(:)=edg_dmg_coord(:,2,ne_dmg)
     call cened_dmg_dbl(ne_dmg,leng_s) 
  end if
  
  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
     sub_area(source)=wght_s(source)
  end do

  ro(:)=rm(1:2)    
  call find_layer(ro(:),layer_o)

  do source=1,nqp_s
     !EFIE
     rs(:)=rsrc(1:2,source)
     call find_layer(rs(:),layer_s)
     call find_height

     ! Direct calculation
!     call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!     Gf_nsigu(1)=Gf_tmp
!     print*,'dir',ne,me,Gf_nsigu(1)

     ! Interpolation
     call Gf_interpolation_3d(rs,ro,Gf_nsigu(1))
!     print*,'inter',ne,me,Gf_nsigu(1)
     phimn(source)=sub_area(source)*Gf_nsigu(1)
  end do

  wghts_phi=wghts_phi+sum(phimn(1:nqp_s))
  wghts_phi=wghts_phi*leng_s/eps0d
  return
end subroutine source_surf_ns_near3
