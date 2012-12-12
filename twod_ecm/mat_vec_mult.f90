module mat_vec_mult
  use global_com
  use global_geom,only:nglunk,nsuunk
  use global_dim,only:prcdin
  implicit none
  save
  real(kind=dp),allocatable::r0_initial(:)
contains
  subroutine matvec(xin,xout1)
    use global_dim,only:zpast
    implicit none
    real(kind=dp),intent(in)::xin(1:nglunk)
    real(kind=dp),intent(out)::xout1(1:nglunk)
    real(kind=dp)::xout(1:nglunk)
    real(kind=dp)::time_near_it,time_dummy(2)
    integer::Itime_near_it,Itime_dummy(2)

call system_clock(Itime_dummy(2)); time_dummy(2)=real(Itime_dummy(2))/real(Itim_rate)
    xout(:)=0.d0
    ! matrix-vector multiplication
    ! C*Q
    call multiply_near
call system_clock(Itime_dummy(1)); time_dummy(1)=real(Itime_dummy(1))/real(Itim_rate)
    time_near_it=time_dummy(1)-time_dummy(2)
    tot_near_time=tot_near_time+time_near_it

    ! only for iterative solver
    ! left multiply preconditioner matrix
    ! P*C*Q
    call precon(xout,xout1)
    return
  contains
    subroutine multiply_near
      integer::dummy,zinteract_counter
      zinteract_counter=0
      do dummy=1,nglunk
         xout(1:nglunk)=&
              xout(1:nglunk)+&
              xin(dummy)*zpast(zinteract_counter:zinteract_counter+nglunk-1)
         zinteract_counter=zinteract_counter+nglunk
      end do
      return
    end subroutine multiply_near
  end subroutine matvec

  subroutine precon(u,v)
    ! is a left pre-conditioner!
    implicit none
    real(kind=dp),intent(in)::u(1:nglunk)
    real(kind=dp),intent(out)::v(1:nglunk)

    ! diagonal preconditioner
    v(1:nglunk)=&
         u(1:nglunk)*prcdin(1:nglunk)
    return
  end subroutine precon

  subroutine invert_preconditioner
    implicit none

    ! diagonal preconditioner
    prcdin(1:nglunk)=&
         1.0_dp/prcdin(1:nglunk)
    return
  end subroutine invert_preconditioner

  subroutine initialize_r0
    implicit none
    integer::jran,i
    allocate(r0_initial(1:nglunk))

    jran=-3
    do i=1,nglunk
       r0_initial(i)=2.0_dp*real(ran1(jran),dp)-1.0_dp
    end do
    return
  end subroutine initialize_r0

  function ran1(idum)
    implicit none
    integer,intent(inout)::idum
    integer,PARAMETER::IM=2147483647,IQ=127773,IR=2836,IA=16807,NTAB=32
    integer::NDIV
    real,parameter::EPS=1.2E-7,RNMX=1.-EPS
    real::am,ran1
    integer::j,k,iv(NTAB)=0,iy=0
    SAVE iv,iy

    NDIV=1+(IM-1)/NTAB
    AM=1./IM

    if (idum<= 0 .or. iy== 0) then
       idum=max(-idum,1)
       do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum .lt. 0) idum=idum+IM
          if (j .le. NTAB) iv(j)=idum
       end do
       iy=iv(1)
    end if
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum .lt. 0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
    ran1=min(AM*iy,RNMX)
    return
  end function ran1
end module mat_vec_mult





