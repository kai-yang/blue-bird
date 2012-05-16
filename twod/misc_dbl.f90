module misc_dbl
  use global_com,only:dp
  implicit none
  save
contains
  !*********************************************************************  
  subroutine rvpcr_dbl(v1,v2,vr)                                      
    ! computes the cross product of v1,v2                        
    implicit none
    real(kind=dp),intent(in)::v1(3),v2(3)
    real(kind=dp),intent(out)::vr(3)              
    !                                                                       
    vr(1)=v1(2)*v2(3)-v2(2)*v1(3)                         
    vr(2)=v1(3)*v2(1)-v2(3)*v1(1)                         
    vr(3)=v1(1)*v2(2)-v2(1)*v1(2)                         
    !                                                                       
    return
  end subroutine rvpcr_dbl
  !*********************************************************************  
  subroutine cened_dbl(m,rl)                                       
    ! Finds the length(rl) of the edge m
    use global_geom,only:edg_coord
    implicit none
    integer,intent(in)::m
    real(kind=dp),intent(out)::rl!,rc(3)
    real(kind=dp),dimension(2)::r1,r2                            
    !                                                                       
    r1(:)=edg_coord(:,1,m)
    r2(:)=edg_coord(:,2,m)
    !  rc=0.5_dp*(r1+r2)
    !                                                                       
    rl=sqrt(dot_product(r1-r2,r1-r2))
    return
  end subroutine cened_dbl
  !*********************************************************************  
  subroutine cened_dmg_dbl(m,rl)                                       
    ! Finds the length(rl) of the edge m
    use global_geom,only:edg_dmg_coord
    implicit none
    integer,intent(in)::m
    real(kind=dp),intent(out)::rl!,rc(3)
    real(kind=dp),dimension(2)::r1,r2                            
    !                                                                       
    r1(:)=edg_dmg_coord(:,1,m)
    r2(:)=edg_dmg_coord(:,2,m)
    !  rc=0.5_dp*(r1+r2)
    !                                                                       
    rl=sqrt(dot_product(r1-r2,r1-r2))
    return
  end subroutine cened_dmg_dbl
  !*********************************************************************  
  subroutine cenedg_dbl(m,rc)
    ! Find the coordinate of the center of patch m
    use global_geom,only:edg_coord
    implicit none
    integer,intent(in)::m
    real(kind=dp),dimension(2),intent(out)::rc

    rc(1:2)=(edg_coord(:,1,m)+edg_coord(:,2,m))/2.0_dp
    
    return
  end subroutine cenedg_dbl
  !*********************************************************************  
  subroutine cenedg_dmg_dbl(m,rc)
    ! Find the coordinate of the center of patch m
    use global_geom,only:edg_dmg_coord
    implicit none
    integer,intent(in)::m
    real(kind=dp),dimension(2),intent(out)::rc

    rc(1:2)=(edg_dmg_coord(:,1,m)+edg_dmg_coord(:,2,m))/2.0_dp
    
    return
  end subroutine cenedg_dmg_dbl
  !*********************************************************************
  ! TESTING FUNCTIONS 
  !*********************************************************************
  subroutine tstvec_sur(me,tstv,rm)
    ! me is the global testing function number
    ! rm holds the coordinates of the quadrature points
    ! tstv is  the weighted values of the RWG basis at q. points
    ! updated by Hakan Bagci(Feb 2004)
    use global_geom,only:edg_coord
    use quadratures,only:qp_t,wght_t,nqp_t
    implicit none
    integer,intent(in)::me
    real(kind=dp),intent(out)::tstv
    real(kind=dp),dimension(2),intent(out)::rm
    real(kind=dp)::v1(2),v2(2),rho(2)
    integer::i
    
    call cenedg_dbl(me,rm)
    tstv=wght_t(1)

    return
  end subroutine tstvec_sur
end module misc_dbl
