subroutine field_diel (me,ne,wghts_phi)
  !     this routine computes the scalar potential and 
  !     compresses them in the wghts array.  this is then used to 
  !     assemble Z_MOM.
  !     last modified : 2011 K. Yang
  use global_com,only:dp,pid,eps0d
  use global_geom,only:nsuinf,edg_diel_coord,edg_diel_epsr
  use misc_dbl,only:cenedg_dbl,cenedg_diel_dbl
  use quadratures,only:nqp_t,wght_t,total_maxqp_t

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),intent(out)::wghts_phi
  real(kind=dp),dimension(2)::rm,rn,vm1,vm2,vm12
  real(kind=dp),dimension(2)::rm_g,um1
  integer::i,ne_diel,me_diel
  real(kind=dp)::self_E
  logical::is_image
  real(kind=dp)::wghts_phi_direct,wghts_phi_image,wghts_phi_tmp
  
  wghts_phi=0.d0; wghts_phi_tmp=0.d0
  !------------------------------------------------------------------------
  !     testing on a conductor surface element:
  !------------------------------------------------------------------------
  me_diel=me-nsuinf(1)
  call cenedg_diel_dbl(me_diel,rm)     
  rm_g(1:2)=rm
  vm1(:)=edg_diel_coord(:,1,me_diel)
  vm2(:)=edg_diel_coord(:,2,me_diel)
  vm12(:)=vm2(:)-vm1(:)
  um1(:)=(/vm12(2),-vm12(1)/) ! counter-clockwise
  um1(:)=um1(:)/sqrt(dot_product(um1(:),um1(:)))
  !       source 
  !       if surface ++++++++++++++++++++++++++++++++++++++++
  source_s:if (ne<=nsuinf(1)) then
     call cenedg_dbl(ne,rn)
  else
     ne_diel=ne-nsuinf(1)
     call cenedg_diel_dbl(ne_diel,rn)
  end if source_s

  ! Image contribution
  is_image=.false.
  call find_direct_intg_diel(ne,me,rn,rm,rm_g,um1,is_image,wghts_phi_direct)
  wghts_phi_tmp=wghts_phi_tmp+wghts_phi_direct
!  print*,'dir',wghts_phi_direct
!  stop

  is_image=.true.
  rn(2)=-rn(2)
  call find_direct_intg_diel(ne,me,rn,rm,rm_g,um1,is_image,wghts_phi_image)
  wghts_phi_tmp=wghts_phi_tmp-wghts_phi_image
!  print*,'image',wghts_phi_image
  if (me==ne) then
     ! self term (dielectric-dielectric surface)
     self_E=(edg_diel_epsr(1,me_diel)+edg_diel_epsr(2,me_diel))/(edg_diel_epsr(1,me_diel)-edg_diel_epsr(2,me_diel))/(2.d0*eps0d)
     wghts_phi_tmp=wghts_phi_tmp-self_E
!     print*,self_E
!     print*,edg_diel_epsr(1:2,me_diel)
  end if
!  print*,wghts_phi_tmp
  wghts_phi=wghts_phi_tmp
  return 
end subroutine field_diel

subroutine find_direct_intg_diel(ne,me,rn,rm,rm_g,um1,is_image,wghts_phi)
  use global_com,only:dp,r_a
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2),um1(2)
  logical,intent(in)::is_image
  real(kind=dp),intent(out)::wghts_phi

  real(kind=dp)::dist,rnear

  wghts_phi=0.d0
  rnear=r_a  ! distance for singularity extractions
  dist=sqrt(sum((rm(:)-rn(:))**2))

  ! near far is based on the center-to-center distances in all cases
  if (dist<=rnear) then
!     print*,ne,'dir near'
     call source_surf_near_diel(ne,me,rm_g(1:2),um1,is_image,wghts_phi)
!     print*,ne,me,wghts_phi_direct
  else
!     print*,ne,'dir far'
     call source_surf_far_diel(ne,me,rm_g(1:2),um1,is_image,wghts_phi)
!     print*,ne,me,wghts_phi_direct
  end if
  return
end subroutine find_direct_intg_diel
!==========================================================================
!==========================================================================
! SOURCE SURFACE (FAR)
!==========================================================================
!==========================================================================
subroutine source_surf_far_diel(ne,me,rm,um1,is_image,wghts_phi)
  ! ne: basis ID
  ! me: testing ID
  ! rm: q.p. coordinates
  ! wghts: vector potential,scalar potential, and H field contributions to Z-matrix
  ! source is a surface function
  ! testing is a surface function (can be a wire function if EFIE_only)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_diel_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_diel_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm,um1
  logical,intent(in)::is_image
  real(kind=dp),intent(out)::wghts_phi

  real(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_diel
  real(kind=dp)::distances2(nqp_s),R(2,nqp_s)
  real(kind=dp)::aaint(2)

  wghts_phi=0.d0

  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_diel=ne-nsuinf(1)
     v1(:)=edg_diel_coord(:,1,ne_diel)
     v2(:)=edg_diel_coord(:,2,ne_diel)
     call cened_diel_dbl(ne_diel,leng_s) 
  end if

  if (is_image) then
     v1(2)=-v1(2)
     v2(2)=-v2(2)
  end if

  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
  end do
    
  do source=1,nqp_s
     !EFIE
     R(:,source)=rm(1:2)-rsrc(1:2,source)
     distances2(source)=dot_product(R(:,source),R(:,source))!R^2
     aaint(:)=wght_s(source)*R(:,source)/distances2(source)
     phimn(source)=um1(1)*aaint(1)+um1(2)*aaint(2)
  end do
  wghts_phi=sum(phimn(1:nqp_s))*leng_s
  wghts_phi=-wghts_phi/(2.d0*pid*eps0d)
  return
end subroutine source_surf_far_diel
!==========================================================================
!==========================================================================
! SOURCE SURFACE (NEAR)
!==========================================================================
!==========================================================================
subroutine source_surf_near_diel(ne,me,rm,um1,is_image,wghts_phi)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_diel_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_diel_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
  implicit none                      
  !
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm,um1
  logical,intent(in)::is_image
  real(kind=dp),intent(out)::wghts_phi

  real(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::aaint1,aaint2(2)
  integer::source,test,ne_diel
  real(kind=dp)::R(2,nqp_s)
  real(kind=dp)::leng_s
  real(kind=dp)::aaint(2),phimn1

  wghts_phi=0.d0
  aaint1=0.d0; aaint2=0.d0

  if (ne<=nsuinf(1)) then
     v1(:)=edg_coord(:,1,ne)
     v2(:)=edg_coord(:,2,ne)
     call cened_dbl(ne,leng_s) 
  else
     ne_diel=ne-nsuinf(1)
     v1(:)=edg_diel_coord(:,1,ne_diel)
     v2(:)=edg_diel_coord(:,2,ne_diel)
     call cened_diel_dbl(ne_diel,leng_s) 
  end if

  if (is_image) then
     v1(2)=-v1(2)
     v2(2)=-v2(2)
  end if

  do source=1,nqp_s
     rhon(:,source)=(v2(:)-v1(:))*qp_s(source)
     rsrc(:,source)=v1(:)+rhon(:,source)
  end do
  ! analytical integrals
  call sintg_2d(rm(1:2),v1,v2,aaint1,aaint2(:))
  phimn1=dot_product(um1,aaint2)
  wghts_phi=-phimn1/(2.d0*pid*eps0d)
  return
end subroutine source_surf_near_diel
