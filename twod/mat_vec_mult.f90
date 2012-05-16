module mat_vec_mult
  use global_com ! nb_local used to be here
  use global_geom,only:nglunk,nsuunk
  use global_dim,only:IPIV,prcdin
  implicit none
  save
!  real(kind=dp)::bmag_sum !used in the computation of the true relative error
!  complex(kind=dp),allocatable::rhsd_store(:)
  complex(kind=dp),allocatable::r0_initial(:)
contains
  subroutine matvec(xin,xout1)!,IPAR)
    use global_dim,only:zpast
    implicit none
    complex(kind=dp),intent(in)::xin(1:nglunk)
    complex(kind=dp),intent(out)::xout1(1:nglunk)
    complex(kind=dp)::xout(1:nglunk)
    real(kind=dp)::time_near_it,time_reduction_it,time_dummy(2)
    integer::Itime_near_it,Itime_reduction_it,Itime_dummy(2)
    complex(kind=dp)::x_local(1:nglunk)

call system_clock(Itime_dummy(2)); time_dummy(2)=real(Itime_dummy(2))/real(Itim_rate)
    x_local(:)=cmplx(0.0_dp,0.0_dp,dp)
    call multiply_near
call system_clock(Itime_dummy(1)); time_dummy(1)=real(Itime_dummy(1))/real(Itim_rate)
    time_near_it=time_dummy(1)-time_dummy(2)
    xout(1:nglunk)=x_local(1:nglunk)

    time_reduction_it=0.0_dp
    tot_near_time=tot_near_time+time_near_it
    tot_reduction_time=tot_reduction_time+time_reduction_it
    call precon(xout,xout1)
    return
  contains
    subroutine multiply_near
      integer::dummy,zinteract_counter
      complex(kind=dp)::wrk(3) 
      zinteract_counter=0
      do dummy=1,nglunk
         x_local(1:nglunk)=&
              x_local(1:nglunk)+&
              xin(dummy)*zpast(zinteract_counter:zinteract_counter+nglunk-1)
         zinteract_counter=zinteract_counter+nglunk
      end do
      return
    end subroutine multiply_near
  end subroutine matvec

  subroutine precon(u,v)
    ! is a left pre-conditioner!
!    use global_com,only:nprecon_ids,my_precon_ids
    implicit none
    complex(kind=dp),intent(in)::u(1:nglunk)
    complex(kind=dp),intent(out)::v(1:nglunk)
    complex(kind=dp),allocatable::precon_all(:),precon_loc(:),&
         invert_for_precon(:,:)
    integer,allocatable::IPIV_precon(:)
    integer::k,shift,info,j
    integer::x,y,z,total
    if (diag_precon) then
       ! diagonal preconditioner
       v(1:nglunk)=&
            u(1:nglunk)*prcdin(1:nglunk,1,1)
    elseif(block_diag_precon) then
       ! 1. map from local unknowns to global unknowns
       ! 2. do an all-to-all communication :(
       ! 3. multiply with the inverse of the block
       ! 4. do an all-to-all communication
       ! 5. map from global unknowns to local unknowns
       ! step 1: 
       allocate(precon_all(1:nglunk),precon_loc(1:nprecon_ids))

       precon_all(:)=cmplx(0.0_dp,0.0_dp,dp)
       precon_all(1:nglunk)=u(1:nglunk)
       
       precon_loc(1:nprecon_ids)=precon_all(1:nprecon_ids)
       ! prcdin is really the transpose of M^-1 ! 
       ! done for ease of multiplication.
       shift=0
       do k=1,nprcd_blocks-1
          call zgetrs('T',prcd_blocksize,1,&
               prcdin(1:prcd_blocksize,1:prcd_blocksize,k),prcd_blocksize,&
               IPIV(1:prcd_blocksize,k),&
               precon_loc(shift+1:shift+prcd_blocksize),prcd_blocksize,info)
!          do j=1,prcd_blocksize
!             v(shift+j)=&
!                  dot_product(prcdin(1:prcd_blocksize,j,k),&
!                  u(shift+1:shift+prcd_blocksize))
!          end do
          shift=shift+prcd_blocksize
       end do
       ! treat the last block specially because it might not be equal sized
       ! with the others
       k=nprcd_blocks
!       do j=1,last_block_size
!          v(shift+j)=&
!               dot_product(prcdin(1:last_block_size,j,k),&
!               u(shift+1:shift+last_block_size))
!       end do
       call zgetrs('T',last_block_size,1,&
            prcdin(1:last_block_size,1:last_block_size,k),last_block_size,&
            IPIV(1:last_block_size,k),precon_loc(shift+1:shift+last_block_size),&
            last_block_size,info)
       precon_all(:)=cmplx(0.0_dp,0.0_dp,dp)
       precon_all(my_precon_ids(1:nprecon_ids))=precon_loc(1:nprecon_ids)
       v(1:nglunk)=precon_all(1:nglunk)
       deallocate(precon_all,precon_loc)
    elseif(near_field_precon) then
    elseif(group_precon) then
    end if
    return
  end subroutine precon

  subroutine invert_preconditioner
    integer::k
    integer::info!,lwork
!    complex(kind=dp)::work(prcd_blocksize)
    if (diag_precon) then
       ! diagonal preconditioner
       prcdin(1:nglunk,1,1)=&
            1.0_dp/prcdin(1:nglunk,1,1)
    elseif(block_diag_precon) then
       ! simple-block-diagonal preconditioner
       do k=1,nprcd_blocks-1
!          if (myid==0 .and. k==1) then
!             print*,'before:',prcdin(:,:,k)
!          end if
          call zgetrf(prcd_blocksize,prcd_blocksize,prcdin(:,:,k),prcd_blocksize,&
               IPIV(1:prcd_blocksize,k),info)
!          lwork=prcd_blocksize
!          call zgetri(prcd_blocksize,prcdin(:,:,k),prcd_blocksize,&
!               IPIV(1:prcd_blocksize,k),work(1:prcd_blocksize),lwork,info)
!          if (myid==0 .and. k==1) then
!             print*,'after:',prcdin(:,:,k)
!          end if
       end do
       k=nprcd_blocks
       call zgetrf(last_block_size,last_block_size,&
            prcdin(1:last_block_size,1:last_block_size,k),last_block_size,&
            IPIV(1:last_block_size,k),info)
!       lwork=last_block_size
!       call zgetri(last_block_size,&
!            prcdin(1:last_block_size,1:last_block_size,k),last_block_size,&
!            IPIV(1:last_block_size,k),work(1:last_block_size),lwork,info)
    elseif (near_field_precon) then
    elseif (group_precon) then
    end if
    return
  end subroutine invert_preconditioner

  subroutine initialize_r0
    implicit none
    integer::jran,i
    real(kind=dp)::r0_dummy(1:2*nglunk)
    allocate(r0_initial(1:nglunk))
!    jran=1211
    jran=-3
    do i=1,2*nglunk
       r0_dummy(i)=2.0_dp*real(ran1(jran),dp)-1.0_dp
    end do
    do i=1,nglunk
       r0_initial(i)=&
         cmplx(r0_dummy(i),r0_dummy(nglunk+i),dp)
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





