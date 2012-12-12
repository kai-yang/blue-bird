subroutine dir_field_mom_only
! This subroutine fills znow and zpast
! also puts znow,nznow into compressed column storage format                 
  use global_com
  use global_dim,only:pmatrix,zpast,prcdin
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
!************************************************
!**************** PRE-CONDITIONER ***************
!************************************************
  call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
  if (is_iter) then
     prcdin(:)=0.d0
  end if
call system_clock(Itim_prcd); tim_prcd=real(Itim_prcd)/real(Itim_rate)-tim_dummy
!************************************************
  mem_est(1:3)=0.0_dp
  zinteract_counter=nglunk**2
  if (is_iter) then
     allocate(zpast(0:zinteract_counter-1))     
  else
     allocate(pmatrix(1:nglunk,1:nglunk))
  end if
  tim_mom=0.0_dp
  zinteract_counter=0
  if (is_iter) then
     source_real_iter:do dummy=1,nglunk
        if (modulo(dummy-1,100)==0) then
#ifdef ENABLED_FOR_MAIN_PROGRAM       
           print*,dummy,'of',nglunk
#endif
        end if
        ne=dummy
        call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
        observer_mom_iter:do count=1,nglunk
           !for each observer
           me=count
           !************************************************
           !***************** GET MOM MATRIX ***************
           !************************************************
           if (me<=nsuinf(1)) then ! observer conductor
              call field(me,ne,wghts_phi)
           else ! observer dielectric
              call field_diel(me,ne,wghts_phi)
           end if
           wghts=wghts_phi
!           print*,dummy,count!wghts_phi
           zpast(zinteract_counter+count-1)=wghts
        end do observer_mom_iter
call system_clock(Itim_mom); tim_mom=tim_mom+real(Itim_mom)/real(Itim_rate)-tim_dummy
        zinteract_counter=zinteract_counter+nglunk
     end do source_real_iter
!************************************************
!**************** PRE-CONDITIONER ***************
!************************************************
     tim_prcd=0.d0
call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
     call do_preconditioner
call system_clock(Itim_prcd); tim_prcd=real(Itim_prcd)/real(Itim_rate)-tim_dummy+tim_prcd
  else
     source_real:do dummy=1,nglunk
        if (modulo(dummy-1,100)==0) then
#ifdef ENABLED_FOR_MAIN_PROGRAM
           print*,dummy,'of',nglunk
#endif
        end if
        ne=dummy
        call system_clock(Itim_dummy); tim_dummy=real(Itim_dummy)/real(Itim_rate)
        observer_mom:do count=1,nglunk ! conductor testing
           !for each observer
           me=count
           !************************************************
           !***************** GET MOM MATRIX ***************
           !************************************************
           if (me<=nsuinf(1)) then ! conductor      
!              print*,'cond'
              call field(me,ne,wghts_phi)
           else ! dielectric
!              print*,'diel'
              call field_diel(me,ne,wghts_phi)
           end if
           wghts=wghts_phi
           if (isnan(wghts)) then
              print*,'NaN',me,ne
              stop
           end if
           pmatrix(me,ne)=wghts
        end do observer_mom
call system_clock(Itim_mom); tim_mom=tim_mom+real(Itim_mom)/real(Itim_rate)-tim_dummy
     end do source_real
  end if

#ifdef ENABLED_FOR_MAIN_PROGRAM
  print*,'DIR-FIELD TIMES(mom/precon):',tim_mom,tim_prcd
#endif
  return
contains
  subroutine do_preconditioner
    do dummy=1,nglunk
       ne=dummy
       if (ne<=nsuinf(1)) then ! conductor
          call field(ne,ne,wghts_phi)
       else ! dielectric
          call field_diel(ne,ne,wghts_phi)
       end if
       wghts=wghts_phi
       prcdin(dummy)=wghts
    end do
    return
  end subroutine do_preconditioner
end subroutine dir_field_mom_only

