!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
subroutine QP_apply_keep_occupations(en,k)
 !
 ! KfnQP_KeepOcc: the BSE uses the QP energies with the occupations of the reference
 ! levels, the ones of the screening and of the GW (setup.F, OCCUPATIONS_Fermi "X").
 ! QP_apply then skips its Fermi-Dirac refill of the QP levels ("K+QP").
 !
 ! Here that refill is done on a copy, only to report what it would have changed:
 ! en is not touched.
 !
 use pars,          ONLY:SP,lchlen
 use units,         ONLY:HA2EV
 use drivers,       ONLY:Finite_Tel
 use com,           ONLY:msg
 use electrons,     ONLY:levels,n_sp_pol,spin_occ,filled_tresh,E_reset,E_duplicate
 use D_lattice,     ONLY:Tel,T_elecs,T_holes,Bose_Temp
 use R_lattice,     ONLY:bz_samp
 use interfaces,    ONLY:OCCUPATIONS_Fermi
 !
 implicit none
 !
 type(levels)  ::en
 type(bz_samp) ::k
 !
 ! Work Space
 !
 integer, parameter :: n_rows_max=200
 real(SP),parameter :: df_report=1.E-3_SP
 type(levels)       :: en_QP
 integer            :: ib,ik,is,ib_range(2),i_max(3),n_rows,n_changed,n_blocked(2)
 real(SP)           :: Tel_SAVE(4),N_moved,df,df_max,Eo_eV
 logical            :: Finite_Tel_SAVE
 character(lchlen)  :: row
 character(1)       :: blocked(2)
 !
 ! The refill of the QP levels, on a copy. OCCUPATIONS_Fermi also updates the temperatures.
 !
 Tel_SAVE=(/Tel,T_elecs,T_holes,Bose_Temp/)
 Finite_Tel_SAVE=Finite_Tel
 !
 call E_reset(en_QP)
 call E_duplicate(en,en_QP)
 call OCCUPATIONS_Fermi(en_QP,k,"K+QP",mode="OCCUPATIONS")
 !
 Tel      =Tel_SAVE(1)
 T_elecs  =Tel_SAVE(2)
 T_holes  =Tel_SAVE(3)
 Bose_Temp=Tel_SAVE(4)
 Finite_Tel=Finite_Tel_SAVE
 !
 ! Differences
 !
 N_moved=0._SP
 df_max =0._SP
 i_max  =0
 n_changed=0
 do is=1,n_sp_pol
   do ik=1,en%nk
     do ib=1,en%nb
       df=en_QP%f(ib,ik,is)-en%f(ib,ik,is)
       N_moved=N_moved+abs(df)*k%weights(ik)/2._SP
       if (abs(df)/spin_occ>df_report) n_changed=n_changed+1
       if (abs(df)>abs(df_max)) then
         df_max=df
         i_max =(/ib,ik,is/)
       endif
     enddo
   enddo
 enddo
 !
 call msg('r','[K+QP] KfnQP_KeepOcc: occupations of the reference levels kept, no refill of the QP levels')
 call msg('r','[K+QP] Fermi level of a QP refill, from the reference one [eV]',en_QP%E_Fermi*HA2EV)
 call msg('r','[K+QP] Electrons per cell a QP refill would move',N_moved)
 call msg('r','[K+QP] Largest occupation change of a QP refill (f/spin_occ)',df_max/spin_occ)
 if (i_max(1)>0) call msg('r','[K+QP] ... at band, IBZ k-point',i_max(1:2))
 call msg('r','[K+QP] States a QP refill changes by more than 0.001',n_changed)
 !
 ! The states that differ: reference and QP energy, kept and refilled occupation.
 ! "full": no room for an electron, (spin_occ-f)/spin_occ <= filled_tresh. K_Transition_check
 ! drops every transition into such a state (Pauli blocking).
 !
 if (n_changed>0) then
   !
   ib_range=(/min(minval(en%nbf(:n_sp_pol)),minval(en_QP%nbf(:n_sp_pol)))+1,&
&             max(maxval(en%nbm(:n_sp_pol)),maxval(en_QP%nbm(:n_sp_pol)))/)
   ib_range(2)=min(ib_range(2),en%nb)
   !
   call msg('r','[K+QP] band, IBZ k, spin | E ref, E QP [eV] | f kept, f QP refill (f/spin_occ) | weight | full: kept, refill')
   n_rows=0
   n_blocked=0
   do is=1,n_sp_pol
     do ik=1,en%nk
       do ib=ib_range(1),ib_range(2)
         if (abs(en_QP%f(ib,ik,is)-en%f(ib,ik,is))/spin_occ<=df_report) cycle
         n_rows=n_rows+1
         blocked='n'
         if ((spin_occ-en%f(ib,ik,is)   )/spin_occ<=filled_tresh) blocked(1)='y'
         if ((spin_occ-en_QP%f(ib,ik,is))/spin_occ<=filled_tresh) blocked(2)='y'
         if (blocked(1)=='y') n_blocked(1)=n_blocked(1)+1
         if (blocked(2)=='y') n_blocked(2)=n_blocked(2)+1
         if (n_rows>n_rows_max) cycle
         Eo_eV=en%E(ib,ik,is)*HA2EV
         if (allocated(en%Eo)) Eo_eV=en%Eo(ib,ik,is)*HA2EV
         write (row,'(a,3i6,2f11.5,2f10.6,f10.6,2a4)') '[K+QP] ',ib,ik,is,Eo_eV,en%E(ib,ik,is)*HA2EV,&
&                en%f(ib,ik,is)/spin_occ,en_QP%f(ib,ik,is)/spin_occ,k%weights(ik),blocked
         call msg('r',trim(row))
       enddo
     enddo
   enddo
   if (n_rows>n_rows_max) call msg('r','[K+QP] ... rows not written',n_rows-n_rows_max)
   call msg('r','[K+QP] Of these states, full with the kept occupations, with a QP refill',n_blocked)
   !
   call warning('[K+QP] KfnQP_KeepOcc: the kept occupations are not a Fermi-Dirac filling of the QP levels')
   !
 endif
 !
 call msg('r','')
 !
 call E_reset(en_QP)
 !
end subroutine
