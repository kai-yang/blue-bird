!*********************************************************************
!  this subroutine solves the system ax=b using transpose free quasi-
!  minimal residual (tfqmr) algorithm. reference:
!    siam j. sci. compt. vol.14, no.2, pp. 470-482, march 93
!       by roland. w. freund
!
!  the program terminates when the required precision is reached.
!  if the required precision is not established only n iterations
!  are carried out.     a.a.ergin may 1995
!
! Modified by a.e. yilmaz 2003
! needs mat_vec_mult and initialize_r0 routine
! that computes an initial random vector r0_initial
subroutine ztfqmr_serial(ntotal,b,x,err,iter)
  use global_com
  use global_geom,only:nglunk
  use mat_vec_mult,only:matvec,r0_initial,precon
  implicit none
  integer,intent(in)::ntotal
  complex(kind=dp),dimension(1:nglunk)::x,bb,b
  real(kind=dp)::err,rerr
  integer::iter,itmax,it
  complex(kind=dp),dimension(1:nglunk)::w,yo,ayo,ye,aye,r,d,v
  real(kind=dp)::ta,we,cm
  complex(kind=dp)::we_local,we_sum,rerr_local,rerr_sum,err_local,err_sum
  complex(kind=dp)::ta_local,ta_sum,bmag_local,bmag_sum1,dumb_ali(6)
  complex(kind=dp)::etha,rho,rho_local,amgis,amgis_local,ahpla,dum,dum_local,beta
  real(kind=dp)::bmag

  itmax=iter
  if (iter.eq.0) itmax=ntotal
  call precon(b(1:nglunk),bb(1:nglunk))
  !
  !  set initial values
  !
  d(1:nglunk)=cmplx(0.0_dp,0.0_dp,dp)
  call matvec(x(1:nglunk),r(1:nglunk))
  r(1:nglunk)=bb(1:nglunk)-r(1:nglunk) !residual from the initial guess
  w(1:nglunk)=r(1:nglunk)
  yo(1:nglunk)=r(1:nglunk)
  call matvec(yo(1:nglunk),ayo(1:nglunk))
  v(1:nglunk)=ayo(1:nglunk)
  we =0.0_dp
  etha=cmplx(0.0_dp,0.0_dp,dp)
  
  ta_local=dot_product(r(1:nglunk),r(1:nglunk))
  rho_local=dot_product(r0_initial(1:nglunk),&
       r(1:nglunk))
  bmag_local=dot_product(bb(1:nglunk),bb(1:nglunk))
  
  dumb_ali(1:3)=(/ta_local,rho_local,bmag_local/)
  dumb_ali(4:6)=dumb_ali(1:3)
  ta_sum=dumb_ali(4);rho=dumb_ali(5);bmag_sum1=dumb_ali(6)

  ta=sqrt(abs(ta_sum))
  bmag=sqrt(abs(bmag_sum1))
  rerr=ta/bmag
  iters: do it=1,itmax
     amgis_local=dot_product(r0_initial(1:nglunk),v(1:nglunk))
     amgis=amgis_local

     ahpla=rho/amgis
     ye(1:nglunk)=yo(1:nglunk)-ahpla*v(1:nglunk)
     call matvec(ye(1:nglunk),aye(1:nglunk))
   !  start odd (2n-1) m loop
     d(1:nglunk)=yo(1:nglunk)+(we*we*etha/ahpla)*d(1:nglunk)
     w(1:nglunk)=w(1:nglunk)-ahpla*ayo(1:nglunk)
     we_local=dot_product(w(1:nglunk),w(1:nglunk))
     we_sum=we_local

     we=sqrt(abs(we_sum))/ta
     cm=1.0d0/sqrt(1.0d0+we*we)
     ta=ta*we*cm
     etha=ahpla*cm*cm
     x(1:nglunk)=x(1:nglunk)+etha*d(1:nglunk)
!  check if the result has converged.
!a        if (err*bmag .gt. ta*sqrt(2.*it)) then
!
!  start even (2n)  m loop
     d(1:nglunk)=ye(1:nglunk)+(we*we*etha/ahpla)*d(1:nglunk)
     w(1:nglunk)=w(1:nglunk)-ahpla*aye(1:nglunk)
     we_local=dot_product(w(1:nglunk),w(1:nglunk))
     we_sum=we_local

     we=sqrt(abs(we_sum))/ta
     cm=1.0d0/sqrt(1.0d0+we*we)
     ta=ta*we*cm
     etha=ahpla*cm*cm
     x(1:nglunk)=x(1:nglunk)+etha*d(1:nglunk)
   !  check if the result has converged.
     if (mod(it-1,5)==0 .or. rerr<5.0_dp*err) then
        call matvec(x(1:nglunk),r(1:nglunk))
        r(1:nglunk)=bb(1:nglunk) -r(1:nglunk)
        rerr_local=dot_product(r(1:nglunk),r(1:nglunk))
        rerr_sum=rerr_local
        rerr=sqrt(abs(rerr_sum))/bmag
        print*,'#ofiter,error:',it,rerr
        if (err > rerr) then
           err=rerr
           iter=it
           return
        endif
     end if
     !  make preparations for next iteration
     dum_local=dot_product( r0_initial(1:nglunk),w(1:nglunk))
     dum=dum_local

     beta=dum/rho
     rho=dum
     yo(1:nglunk)=w(1:nglunk)+beta*ye(1:nglunk)
     call matvec(yo(1:nglunk),ayo(1:nglunk))
     !MAGIC
     v(1:nglunk)=ayo(1:nglunk)+beta*( aye(1:nglunk)+beta*v(1:nglunk) )
  enddo iters
  !
  call matvec(x(1:nglunk),r(1:nglunk))
  !MAGIC
  r(1:nglunk)=bb(1:nglunk)-r(1:nglunk)
  err_local=dot_product(r(1:nglunk),r(1:nglunk))
  err_sum=err_local

  err=sqrt(abs(err_sum))/bmag
  iter=itmax
  return
end subroutine ztfqmr_serial
