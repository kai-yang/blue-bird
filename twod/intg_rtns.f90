!*********************************************************************  
subroutine sintg_3d(un1,vo,v1,v2,v3,v12,v13,one_ov_r,one_ov_r3)
! computes analytical integration over triangle:
! 1/R dS',1/R^3 dS'
! combines them and returns analytical integration over triangle of:
! one_ov_r=1/R dS', one_ov_r3=1/R dS'
!  this subroutine evaluates the analytical integrations associated
!  with planar triangular pec surface patches
!  vo : observation point
!  v1,v2,v3 : vertices of the triangle
!  one_ov_r : scalar potential contribution
!  one_ov_r3 : derivative of scalar potential contribution
! assumes (v12=v2-v1)x(v13=v3-v1) is the normal of this patch given in un1
  use global_com,only:dp
  use global_geom,only:edge_av
  use misc_dbl,only:rvpcr_dbl
  implicit none                      
  real(kind=dp),intent(in)::un1(3)
  real(kind=dp),dimension(3),intent(in)::vo,v1,v2,v3,v12,v13
  real(kind=dp),intent(out)::one_ov_r,one_ov_r3
  
  real(kind=dp)::vo1(3),abs_d,pl(3),rhov(3)
  integer::isignd,i
  real(kind=dp)::vt1(3),vt2(3),vt3(3),vtt(3),rm,vo2(3),rp,pp(3),pm(3),ul(3)
  real(kind=dp)::rlp,rlm,uu(3),p0,r0_squared,t1,betai
! find distance d and normal vector n                                         
  vo1=v1-vo
  abs_d=dot_product(vo1,un1)
  pl=vo+abs_d*un1
  rhov=pl-v1! vector from source vertex to obs point projection 
  if (abs_d<0) then  
     abs_d=-abs_d                                                         
     isignd=-1
  else
     isignd=1
  endif
! determine contribution for every edge                             
  one_ov_r=0.0_dp
  one_ov_r3=0.0_dp
  vt1=v1; vt2=v2; vt3=v3 
  do i=1,3 ! loop over source patch edges
     vo1=vt1-vo 
     rm=sqrt(dot_product(vo1,vo1)) ! R-
     vo2=vt2-vo 
     rp=sqrt(dot_product(vo2,vo2)) ! R+
     pp=vt2-pl! P+ 
     pm=vt1-pl! P-
     ul=vt2-vt1!pp-pm
     ul=ul/sqrt(dot_product(ul,ul))
     ! Compute unit vector tangent to edge (ul)
     rlp=dot_product(pp,ul)                                            
     rlm=dot_product(pm,ul)                                            
! Compute length from P0 (projection of pl onto line of source edge)to vertices
! l+ (rlp) and l- (rlm) 
     call rvpcr_dbl(ul,un1,uu)  
     p0=dot_product(pp,uu)!P0 
     ! this is absolute value in wilton, ..., butler, TAP, Mar.1984
     ! but is without absolute value in graglia, TAP, Oct. 1993 
     ! (this is t_i^0 in graglia's paper)
     ! graglias notation helps get rid of t3, and hence is preferred here.
     r0_squared=p0**2+abs_d**2   !R0  
     if ((abs(p0)<1.d-8*edge_av)) then
        ! anint=0 since p0,t2p,t2m ->0
        if (abs_d<1.d-8*edge_av) then 
           ! i.e., we are almost in the same line as this edge
           if (rlp<0) then !rlm is also <0
              t1=-log((rp-rlp)/(rm-rlm)) ! this avoids a 0/0 division
              ! and gives the correct integral when rp->-rlp and rm->-rlm
           elseif (rlm>0) then ! rlp is also >0
              t1=log((rp+rlp)/(rm+rlm))
           else !rlp>0 and rlm<0
              if ((rm+rlm)<1.d-12*edge_av) then
                 print*,'MFIE-You are testing (almost) on the edge of the source!!!'
                 t1=0.0_dp
                 ! we set it to zero although it is really infinity
                 ! but this will no be multiplied by r0_squared
              else
                 t1=log((rp+rlp)/(rm+rlm))
              end if
           end if
           ! because r0_squared is very small we don't need to do the above
           ! we could just set t1=0
        else                                                         
           t1=log((rp+rlp)/(rm+rlm))   
        endif
     else
        betai=atan(p0*rlp/(r0_squared+abs_d*rp))-&
             atan(p0*rlm/(r0_squared+abs_d*rm))
        t1=log((rp+rlp)/(rm+rlm))               
        one_ov_r=one_ov_r+p0*t1-abs_d*betai
        one_ov_r3=one_ov_r3+isignd*betai
     endif
     vtt=vt1; vt1=vt2; vt2=vt3; vt3=vtt        
  end do
  ! one_ov_r3 holds: 1/R^3 dS'*d
  return
end subroutine sintg_3d

!*********************************************************************  
subroutine sintg_2d(vo,v1,v2,lnp,rhat_ov_r)
! computes analytical integration over line:
! ln(P) dL', (rho-rho')/abs(rho-rho')^2 dL'
! combines them and returns analytical integration over line of:
! lnP=log(rho-rho') dL'
!  this subroutine evaluates the analytical integrations associated
!  with pec line segaments
!  vo : observation point
!  v1,v2 : vertices of the line
!  one_ov_r : scalar potential contribution
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
  up0=up0/sqrt(dot_product(up0,up0))

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
