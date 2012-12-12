subroutine field (me,ne,wghts_phi)
  !     this routine computes the scalar potential and 
  !     compresses them in the wghts array.  this is then used to 
  !     assemble Z_MOM.
  !     last modified : 2011 K. Yang
  use global_com,only:dp,pid
  use global_geom,only:nsuinf
  use misc_dbl,only:cenedg_dbl,cenedg_diel_dbl
  use quadratures,only:nqp_t,wght_t,total_maxqp_t

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),intent(out)::wghts_phi 
  real(kind=dp),dimension(2)::rm,rn
  real(kind=dp),dimension(2)::rm_g
  integer::i,ne_diel
  logical::is_image
  real(kind=dp)::wghts_phi_direct,wghts_phi_image,wghts_phi_tmp

  wghts_phi=0.d0; wghts_phi_tmp=0.d0
  !------------------------------------------------------------------------
  !     testing on a conductor surface element:
  !------------------------------------------------------------------------
  call cenedg_dbl(me,rm)     
  rm_g(1:2)=rm
  !       source 
  !       if surface ++++++++++++++++++++++++++++++++++++++++
  source_s:if (ne<=nsuinf(1)) then
     call cenedg_dbl(ne,rn)
  else
     ne_diel=ne-nsuinf(1)
     call cenedg_diel_dbl(ne_diel,rn)
  end if source_s

  is_image=.false.
  call find_direct_intg(ne,me,rn,rm,rm_g,is_image,wghts_phi_direct)
  wghts_phi_tmp=wghts_phi_tmp+wghts_phi_direct
!  print*,'dir',wghts_phi_direct

  ! image the source
  is_image=.true. 
  rn(2)=-rn(2)
  call find_direct_intg(ne,me,rn,rm,rm_g,is_image,wghts_phi_image)
  wghts_phi_tmp=wghts_phi_tmp-wghts_phi_image
!  print*,'image',wghts_phi_image
!  print*,wghts_phi_tmp
!  stop
  wghts_phi=wghts_phi_tmp
  return 
end subroutine field

subroutine find_direct_intg(ne,me,rn,rm,rm_g,is_image,wghts_phi)
  use global_com,only:dp,r_a
  implicit none

  integer,intent(in)::ne,me
  real(kind=dp),intent(in)::rn(2),rm(2),rm_g(2)
  logical,intent(in)::is_image
  real(kind=dp),intent(out)::wghts_phi

  real(kind=dp)::dist,rnear

  wghts_phi=0.d0
  rnear=r_a  ! distance for singularity extractions
  dist=sqrt(sum((rm(:)-rn(:))**2))
     
  ! near far is based on the center-to-center distances in all cases
  if (dist<=rnear) then
!     print*,'dir near'
     call source_surf_near(ne,me,rm_g(1:2),is_image,wghts_phi)
  else
!     print*,'dir far'
     call source_surf_far(ne,me,rm_g(1:2),is_image,wghts_phi)
  end if
  return
end subroutine find_direct_intg
!==========================================================================
!==========================================================================
! SOURCE SURFACE (FAR)
!==========================================================================
!==========================================================================
subroutine source_surf_far(ne,me,rm,is_image,wghts_phi)
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
!  use layers,only:layer_ne,eps_t,k_prop,euler

  implicit none                      
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  logical,intent(in)::is_image
  real(kind=dp),intent(out)::wghts_phi

  real(kind=dp),dimension(nqp_s)::phimn
  real(kind=dp)::v1(2),v2(2),rsrc(2,nqp_s),rhon(2,nqp_s)
  real(kind=dp)::leng_s
  integer::test,source,ne_diel
  real(kind=dp)::aaint,distances(nqp_s),R(2,nqp_s)

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
     distances(source)=sqrt(dot_product(R(:,source),R(:,source)))!R
     aaint=wght_s(source)*log(distances(source))
     phimn(source)=aaint
  end do
  wghts_phi=-sum(phimn(1:nqp_s))*leng_s
  wghts_phi=wghts_phi/(2.d0*pid*eps0d)
  return
end subroutine source_surf_far
!==========================================================================
!==========================================================================
! SOURCE SURFACE (NEAR)
!==========================================================================
!==========================================================================
subroutine source_surf_near(ne,me,rm,is_image,wghts_phi)
  use global_com,only:dp,eps0d,pid
  use global_geom,only:edg_coord,edg_diel_coord,nsuinf
  use misc_dbl,only:cened_dbl,cened_diel_dbl
  use quadratures,only:qp_s,wght_s,total_maxqp_t,nqp_s
!  use layers,only:layer_ne,eps_t,euler,k_prop

  implicit none                      
  !
  !depending on the testing scheme one_or_seven is 1 or 7
  integer,intent(in)::me,ne
  real(kind=dp),dimension(2),intent(in)::rm
  logical,intent(in)::is_image
  real(kind=dp),intent(out)::wghts_phi

  real(kind=dp)::v1(2),v2(2)
  real(kind=dp)::aaint1,aaint2(2)
  integer::test,ne_diel
  real(kind=dp)::phimn1,leng_s

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
  
  ! analytical e- and h-field integrals
  call sintg_2d(rm(1:2),v1,v2,aaint1,aaint2)
  phimn1=aaint1
  wghts_phi=-phimn1
  wghts_phi=wghts_phi/(2.d0*pid*eps0d)
  return
end subroutine source_surf_near
