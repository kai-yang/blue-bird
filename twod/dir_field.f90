subroutine dir_field_mom_only
! This subroutine fills znow and zpast
! also puts znow,nznow into compressed column storage format                 
  use global_com!,only:dp,dtd,eps0d,max_INF,CFIE,alpha,ndtmax,lin_tim,&
!       nprcd_blocks,prcd_blocksize,last_block_size,EFIE_only,MFIE_only,&
!       use_dc_loops,use_mot_prolates,ndtmax_near,integral_form,&
!       myid,master_id,CART_COMM,ierr,int_mem,real_mem,precon_distx
  use global_dim,only:pmatrix,zpast,prcdin,IPIV
  use global_geom,only:nglunk,nsuunk,nsuinf
 !  use loop_tree,only:receive_rwg_loops,form_loops!,send_rwg_loops
!
  implicit none
  real(kind=dp)::wghts,wghts_phi
  integer::me,ne,count
  integer::dummy
  real(kind=dp)::mem_est(3)
  real(kind=dp)::tim_mom,tim_dummy,tim_prcd
  integer::Itim_mom,Itim_dummy,Itim_prcd
  integer::zinteract_counter
  integer::prcd_block,prcd_so,prcd_ob    
!************************************************
!**************** PRE-CONDITIONER ***************
!************************************************
  call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
  if (is_iter) then
     if (diag_precon) then
        prcdin(:,:,:)=cmplx(0.0_dp,0.0_dp,dp)
     elseif(block_diag_precon) then
        prcdin(:,:,:)=cmplx(0.0_dp,0.0_dp,dp)
     elseif(near_field_precon) then
        ! TO-BE-CODED
     elseif(group_precon) then
        ! TO-BE-CODED
     end if
  end if
call system_clock(Itim_prcd); tim_prcd=real(Itim_prcd)/real(Itim_rate)-tim_dummy
!************************************************
  mem_est(1:3)=0.0_dp
  zinteract_counter=(nsuinf(1)+nsuinf(2))**2
  if (is_iter) then
     allocate(zpast(0:zinteract_counter-1))     
  else
     allocate(pmatrix(1:nsuinf(1)+nsuinf(2),1:nsuinf(1)+nsuinf(2)))     
  end if
  tim_mom=0.0_dp
  zinteract_counter=0
  if (is_iter) then
     source_real_iter:do dummy=1,nsuinf(1)+nsuinf(2)
        if (modulo(dummy-1,30)==0) then
           print*,dummy,'of',nsuinf(1)+nsuinf(2)
        end if
        ne=dummy
        call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
        observer_mom_iter:do count=1,nsuinf(1)+nsuinf(2)
           !for each observer
           me=count
           !************************************************
           !***************** GET MOM MATRIX ***************
           !************************************************
           if (count<=nsuinf(1)) then ! conductor
              call field(me,ne,wghts_phi)
           else ! conformal dielectric
!              call field_dmg(me,ne,wghts_phi)
           end if
           wghts=wghts_phi
!           print*,dummy,count!wghts_phi
           zpast(zinteract_counter+count-1)=wghts
        end do observer_mom_iter
call system_clock(Itim_mom); tim_mom=tim_mom+real(Itim_mom)/real(Itim_rate)-tim_dummy
        zinteract_counter=zinteract_counter+nsuinf(1)+nsuinf(2)
     end do source_real_iter
!************************************************
!**************** PRE-CONDITIONER ***************
!************************************************
     tim_prcd=0.d0
call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
     call do_preconditioner
call system_clock(Itim_prcd); tim_prcd=real(Itim_prcd)/real(Itim_rate)-tim_dummy+tim_prcd
  else
     source_real:do dummy=1,nsuinf(1)+nsuinf(2)
        if (modulo(dummy-1,30)==0) then
           print*,dummy,'of',nsuinf(1)+nsuinf(2)
        end if
        ne=dummy
        call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
        observer_mom:do count=1,nsuinf(1)+nsuinf(2) ! conductor testing
           !for each observer
           me=count
           !************************************************
           !***************** GET MOM MATRIX ***************
           !************************************************
           if (count<=nsuinf(1)) then
              call field(me,ne,wghts_phi)
           else
              call field_dmg(me,ne,wghts_phi)
           end if
           wghts=wghts_phi
!           print*,ne,count,wghts_phi
           pmatrix(me,ne)=wghts
        end do observer_mom
!        print*,ne,sum(pmatrix(:,ne))
call system_clock(Itim_mom); tim_mom=tim_mom+real(Itim_mom)/real(Itim_rate)-tim_dummy
     end do source_real
  end if
!************************************************
  print*,'DIR-FIELD TIMES(mom/precon):',tim_mom,tim_prcd
  return
contains
  subroutine do_preconditioner
    integer::prcd_info
    integer::preconi,preconj,precondummy,preconsearch
    integer::z,y,x,total,maxblocksize
    if (diag_precon) then
       do dummy=1,nsuinf(1)+nsuinf(2)
          ne=dummy
           if (ne<=nsuinf(1)) then ! conductor
              call field(ne,ne,wghts_phi)
           else ! conformal dielectric
              call field_dmg(ne,ne,wghts_phi)
           end if
          wghts=wghts_phi
          prcdin(dummy,1,1)=wghts
       end do
    elseif(block_diag_precon) then
       source_precon:do dummy=1,nprecon_ids
          if (modulo(dummy-1,100)==0) then
             print*,'PRECOND',dummy,'of',nprecon_ids
          end if
          
          prcd_block=floor(dble(dummy-1)/dble(prcd_blocksize))+1
          prcd_so=dummy-(prcd_block-1)*prcd_blocksize
          
          ne=dummy
          ! we will be updating
          ! prcdin(prcd_so,:,prcd_block)<-this is the transpose of the matrix!
          obs_precon: do count=(prcd_block-1)*prcd_blocksize+1,&
               min(prcd_block*prcd_blocksize,nprecon_ids)
             ! go over the observers in the block of this source
             prcd_ob=count-(prcd_block-1)*prcd_blocksize
             me=count
             call field(me,ne,wghts_phi)
             wghts=wghts_phi
             prcdin(prcd_so,prcd_ob,prcd_block)=wghts
          end do obs_precon
       end do source_precon
    elseif(near_field_precon) then
     ! TO-BE-CODED
    elseif(group_precon) then
     ! TO-BE-CODED
    end if
    return
  end subroutine do_preconditioner
end subroutine dir_field_mom_only

