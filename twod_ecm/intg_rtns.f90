!*********************************************************************  
subroutine sintg_2d(vo,v1,v2,lnp,rhat_ov_r)
! computes analytical integration over line:
! ln(P) dL', (rho-rho')/abs(rho-rho')^2 dL'
! combines them and returns analytical integration over line of:
! lnP=log(rho-rho') dL'
! this subroutine evaluates the analytical integrations associated
! with line segaments
! vo : observation point
! v1,v2 : vertices of the line
! one_ov_r : scalar potential contribution
  use global_com,only:dp
  use global_geom,only:edge_av
  use misc_dbl,only:rvpcr_dbl
  implicit none                      
  real(kind=dp),dimension(2),intent(in)::vo,v1,v2
  real(kind=dp),intent(out)::lnp,rhat_ov_r(2)
  
  real(kind=dp)::vt1(2),vt2(2),rm,rp,vo1(2),vo2(2),ul(2),up0(2)
  real(kind=dp)::un1_tmp(3),up0_tmp(3),vo2_tmp(3),ul_tmp(3)
  real(kind=dp)::pp,pm,rlp,rlm,uu(2),p0,t1,betai,one_ov_p

  lnp=0.0_dp; one_ov_p=0.d0; rhat_ov_r(:)=0.0_dp
  vt1=v1; vt2=v2

  vo1=vt1-vo
  pm=sqrt(dot_product(vo1,vo1)) ! P-
  vo2=vt2-vo
  pp=sqrt(dot_product(vo2,vo2)) ! P+

  ul=vt2-vt1!pp-pm
  ul=ul/sqrt(dot_product(ul,ul))
  ! Compute unit vector tangent to edge (ul)
  rlp=dot_product(vo2,ul)                                       
  rlm=dot_product(vo1,ul)                                            
  ! Compute length from P0 (projection of pl onto line of source edge)to vertices
  ! l+ (rlp) and l- (rlm) 
  uu=(/-ul(2),ul(1)/) ! uu is perpendicular to ul (normal)
  uu=uu/sqrt(dot_product(uu,uu))

  p0=dot_product(vo2,uu)!P0 
  p0=abs(p0)

  vo2_tmp(1)=vo2(1); vo2_tmp(2)=0.d0; vo2_tmp(3)=vo2(2)
  ul_tmp(1)=ul(1); ul_tmp(2)=0.d0; ul_tmp(3)=ul(2)
  call rvpcr_dbl(vo2_tmp,ul_tmp,un1_tmp) 
  call rvpcr_dbl(ul_tmp,un1_tmp,up0_tmp) 
  up0(1)=up0_tmp(1); up0(2)=up0_tmp(3)
  t1 = dot_product(up0, up0)
  if (t1 .NE. 0) then
     up0=up0/sqrt(dot_product(up0,up0))
  endif

  ! this is absolute value in wilton, ..., butler, TAP, Mar.1984
  if ((abs(p0)<1.d-8*edge_av)) then
     ! anint=0 since p0,t2p,t2m ->0
     t1=rlp*log(pp)-rlm*log(pm)-(rlp-rlm)
     lnp=t1
     rhat_ov_r(:)=0.d0
  else
     betai=p0*(atan(rlp/p0)-atan(rlm/p0))
     t1=rlp*log(pp)-rlm*log(pm)-(rlp-rlm)
     lnp=t1+betai
     rhat_ov_r(:)=up0(:)*(atan(rlp/p0)-atan(rlm/p0))
  endif
  rhat_ov_r(:)=rhat_ov_r(:)+ul(:)*(log(pp)-log(pm))
  rhat_ov_r(:)=-rhat_ov_r(:) ! In TAP, the figure shows rs-ro rather than ro-rs
  return
end subroutine sintg_2d
