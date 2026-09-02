!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details): 
!
subroutine W_potential(q,X,Xw)
 !
 ! Build the screened Coulomb potential W(q,G,G') in the plane-wave basis from the
 ! inverse dielectric matrix stored in ndb.em1s (static) or ndb.pp (static + plasmon-pole)
 ! and write it to ndb.W (see mod_Wpot for the conventions).
 !
 ! The construction is the same used in QP_ppa_cohsex, K_screened_interaction and LSEX_potential:
 !
 !   W(q,G,G') = [ delta_GG' + X_mat(q,G,G';w=0) ] * gamp(q,G,G') * 4 pi / q_weight
 !
 ! with gamp from scatter_Gamp (mode 'c'), i.e. including the RIM regularization of the q->0
 ! limit when RIM is active, and q_weight = d3q_factor/(2 pi)^3, so that away from q->0
 !
 !   gamp * 4 pi / q_weight = 4 pi / ( |q+G| |q+G'| )
 !
 ! Plasmon-pole energies (Godby-Needs), as in QP_ppa_cohsex:
 !
 !   Omega(q,G,G') = ppaE / sqrt( X_mat(w=0)/X_mat(w=i ppaE) - 1 )
 !
 ! Note that the -4 pi/spin_occ prefactor applied by yambo when W is contracted in the
 ! self-energy is NOT included: W is stored in Hartree atomic units.
 !
 use pars,          ONLY:SP,pi,cZERO,cONE
 use units,         ONLY:HA2EV
 use drivers,       ONLY:l_rim_w
 use com,           ONLY:msg
 use LIVE_t,        ONLY:live_timing
 use parallel_m,    ONLY:master_cpu
 use parallel_int,  ONLY:PP_wait
 use frequency,     ONLY:w_samp
 use R_lattice,     ONLY:bz_samp,d3q_factor,RIM_qpg
 use X_m,           ONLY:X_t,X_mat
 use ALLOC,         ONLY:X_ALLOC_elemental
 use collision_el,  ONLY:elemental_collision,elemental_collision_free,elemental_collision_alloc
 use Wpot_m,        ONLY:W_is_PPA,W_ng,W_n_freqs,W_screening,W_alloc,W_free,&
&                        W_freqs,W_vcoul,W_epsm1,W_static,W_ppa_E,&
&                        W_PPcond_rate,W_TO_rate,W_PP_err
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:OP_RD_CL,OP_RD,OP_WR_CL,OP_APP_CL,RD_CL_IF_END,REP,NONE,VERIFY,&
&                        manage_action,deliver_IO_error_message,io_RESPONSE
 use timing_m,      ONLY:timing
 !
#include<y_memory.h>
 !
 type(bz_samp), intent(in) :: q
 type(X_t)                 :: X
 type(w_samp)              :: Xw
 !
 ! Work Space
 !
 type(elemental_collision) :: isc
 integer                   :: iq,ig1,ig2,ID_X,ID_W,IO_ACT,io_err,iq_first,iq_last,PPcond_rate,TO_rate
 integer,      external    :: io_X,io_W
 real(SP)                  :: q_weight,gamp_fac,X_ppaE,PP_err,X_norm
 complex(SP)               :: X1,X2,ratio,sq_ratio
 !
 call section('*','Screened Coulomb potential W(q,G,G'')')
 !==========================================================
 !
 if (X%whoami/=2.and.X%whoami/=4) call error('[W] W(q,G,Gp) needs a static (em1s) or plasmon-pole (ppa) screening')
 if (.not.io_RESPONSE) call error('[W] The screening database I/O is switched off (DBsIOoff)')
 !
 call timing('W_potential',OPR='start')
 !
 W_is_PPA   =(X%whoami==4)
 W_screening="em1s"
 if (W_is_PPA) W_screening="pp"
 W_ng       =X%ng
 W_n_freqs  =Xw%n_freqs
 X_ppaE     =X%ppaE
 !
 ! gamp = q_weight/(|q+G||q+G'|), see scatter_Gamp
 !
 q_weight=d3q_factor/(2._SP*pi)**3
 gamp_fac=4._SP*pi/q_weight
 !
 iq_first=X%iq(1)
 iq_last =min(X%iq(2),q%nibz)
 !
 call msg('r','[W] Screening database           ',trim(W_screening))
 call msg('r','[W] Matrix size                  ',W_ng)
 call msg('r','[W] Screening frequencies        ',W_n_freqs)
 call msg('r','[W] Q-points range               ',(/iq_first,iq_last/))
 if (W_is_PPA) call msg('r','[W] PPA imaginary energy         ',X_ppaE*HA2EV,"[eV]")
 call msg('r','[W] RIM regularized q->0 head    ',allocated(RIM_qpg))
 if (l_rim_w) call warning('[W] The RIM-W corrections (ndb.RIM_W) are not included in W(q,G,Gp)')
 !
 ! Screening DB check (VERIFY mode)
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=ID_X)
 io_err=io_X(X,Xw,ID_X)
 if (io_err<0) call error('[W] Incomplete and/or broken PPA/Static dielectric function database')
 !
 ! The construction is cheap and the I/O is serial: the master does everything
 !
 if (master_cpu) then
   !
   ! Allocations
   !
   call X_ALLOC_elemental('X')
   call X_ALLOC_elemental('X',(/X%ng,X%ng,Xw%n_freqs/))
   call W_alloc(W_ng,W_n_freqs,W_is_PPA)
   !
   call elemental_collision_free(isc)
   isc%ngrho=X%ng
   isc%iqref=0
   call elemental_collision_alloc(isc,NG=X%ng,GAMP_NG=(/X%ng,X%ng/),TITLE="W")
   !
   ! Screening DB header (kept open during the q-loop)
   !
   call io_control(ACTION=OP_RD,COM=NONE,SEC=(/1/),ID=ID_X)
   io_err=io_X(X,Xw,ID_X)
   !
   ! ndb.W header
   !
   call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),ID=ID_W)
   io_err=io_W(X,Xw,ID_W)
   !
   call live_timing('[W] W(q,G,Gp)',iq_last-iq_first+1)
   !
   do iq=iq_first,iq_last
     !
     ! Coulomb factors
     !
     isc%qs(2)=iq
     call scatter_Gamp(isc,'c')
     !
     ! eps^-1_sym - 1 from disk
     !
     IO_ACT=manage_action(RD_CL_IF_END,iq,iq_first,iq_last)
     call io_control(ACTION=IO_ACT,COM=NONE,SEC=(/2*iq,2*iq+1/),ID=ID_X)
     io_err=io_X(X,Xw,ID_X)
     call deliver_IO_error_message(io_err,'PP/Em1s',STOP_it=.TRUE.)
     !
     W_freqs(:W_n_freqs)=Xw%p(:W_n_freqs)
     !
     ! Static part
     !
     W_epsm1=X_mat(:,:,1)
     forall(ig1=1:X%ng) W_epsm1(ig1,ig1)=W_epsm1(ig1,ig1)+cONE
     !
     !$omp parallel do default(shared), private(ig1,ig2), collapse(2)
     do ig2=1,X%ng
       do ig1=1,X%ng
         W_static(ig1,ig2)=W_epsm1(ig1,ig2)*isc%gamp(ig1,ig2)*gamp_fac
       enddo
     enddo
     !$omp end parallel do
     !
     forall(ig1=1:X%ng) W_vcoul(ig1)=isc%gamp(ig1,ig1)*gamp_fac
     !
     ! Plasmon-pole energies
     !
     PPcond_rate=0
     TO_rate    =0
     PP_err     =0._SP
     !
     if (W_is_PPA) then
       !
       !$omp parallel do default(shared), private(ig1,ig2,X1,X2,ratio,sq_ratio,X_norm), &
       !$omp &   reduction(+:PPcond_rate,TO_rate,PP_err), collapse(2)
       do ig2=1,X%ng
         do ig1=1,X%ng
           !
           X1=X_mat(ig1,ig2,1)
           X2=X_mat(ig1,ig2,2)
           !
           ratio=cZERO
           if (abs(X2)>0._SP) ratio=X1/X2
           !
           if (real(ratio)<=1._SP) then
             !
             ! Plasmon-pole condition failure: the pole is set to 1 Ha, as in QP_ppa_cohsex
             !
             W_ppa_E(ig1,ig2)=cONE
             X_norm=max(abs(X1),abs(X2),epsilon(1._SP))
             PP_err=PP_err+abs(X1/(X_ppaE**2+1._SP)-X2)/X_norm
             PPcond_rate=PPcond_rate+1
             !
           else
             !
             sq_ratio=sqrt(ratio-cONE)
             W_ppa_E(ig1,ig2)=X_ppaE/sq_ratio
             if (real(1._SP/sq_ratio,SP)*aimag(1._SP/X1)>0._SP) TO_rate=TO_rate+1
             !
           endif
           !
         enddo
       enddo
       !$omp end parallel do
       !
     endif
     !
     W_PPcond_rate=PPcond_rate
     W_TO_rate    =TO_rate
     W_PP_err     =PP_err
     !
     ! Report
     !
     call msg('r',' Current Q-pt index             ',iq)
     call msg('r',' :: W(q,G=0,Gp=0)   [Ha]        ',(/real(W_static(1,1)),aimag(W_static(1,1))/))
     if (W_is_PPA) then
       call msg('r',' :: PP condition fails/total    ',real(PPcond_rate,SP)/real(X%ng,SP)**2)
       call msg('r',' :: Time ordering fails/rest    ',real(TO_rate,SP)/max(real(X%ng,SP)**2-real(PPcond_rate,SP),1._SP))
       call msg('r',' :: Mean rel dev of PP cond     ',PP_err/real(X%ng,SP)**2)
     endif
     !
     ! I/O
     !
     call io_control(ACTION=OP_APP_CL,COM=REP,SEC=(/2*iq,2*iq+1/),ID=ID_W)
     io_err=io_W(X,Xw,ID_W)
     !
     call live_timing(steps=1)
     !
   enddo
   !
   call live_timing()
   !
   ! Clean
   !
   call elemental_collision_free(isc)
   call X_ALLOC_elemental('X')
   call W_free()
   !
 endif
 !
 call PP_wait()
 !
 call timing('W_potential',OPR='stop')
 !
end subroutine W_potential
