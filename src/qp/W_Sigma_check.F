!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
subroutine W_Sigma_check(X,Xk,E,k,q,qp,Xw)
 !
 ! Recompute the diagonal correlation self-energy Sc(n,k;w=E0) of the QP states from the
 ! screened Coulomb potential W(q,G,G') stored in ndb.W and compare it with QP_Sc, computed
 ! by QP_ppa_cohsex from ndb.pp/ndb.em1s. This validates ndb.W, and the formulas below,
 ! against yambo's internal GW self-energy. Nothing computed here is written in ndb.QP.
 !
 ! Definitions (Hartree atomic units, s=spin_occ, f_m the occupations in [0,s]):
 !
 !   Wc(q,G,G')   = W(q,G,G';w=0) - v(q,G) delta_GG'       (correlation part of W)
 !   rho_m(G)     = < n k | e^{i(q+G).r} | m k-q >         (scatter_Bamp, includes the q-rotation)
 !   q_w          = d3q_factor/(2 pi)^3                    (BZ measure of one q-point)
 !   Om(q,G,G')   = Re[ Omega(q,G,G') ]                    (plasmon-pole energies)
 !
 ! Plasmon-pole (as QP_ppa_cohsex, l_ppa):
 !
 !   Sc(w) = -(q_w/s) Sum_{q in BZ} Sum_m Sum_GG' rho_m(G) rho_m^*(G') Wc(G,G') Om(G,G')/2 *
 !           [ (s-f_m+n_B)/(w-e_m-Om) + (f_m+n_B)/(w^*-e_m+Om) ],  w = E0 + i GDamping
 !
 ! COHSEX (as QP_ppa_cohsex, l_cohsex):
 !
 !   SEX   = -(q_w/s) Sum_q Sum_m f_m Sum_GG' rho_m(G) Wc(G,G') rho_m^*(G')
 !   COH   =  (q_w/2) Sum_q Sum_GG' rho_nn(G-G') Wc(G,G')      (closure relation, no empty states)
 !           [ with empty states (COHSEX_use_empties): (q_w/2) Sum_q Sum_m rho_m Wc rho_m^* ]
 !
 ! q runs over the whole BZ: W is read at the IBZ point and the rotation is carried by rho_m,
 ! except when q is obtained by time-reversal without spatial inversion: then
 ! W(-q,G,G') = W(q,-G',-G) and W is transposed (the -1 is included in the collisions).
 !
 ! The Bruneval-Gonze/BRS terminator (l_GW_terminator) is included as in QP_ppa_cohsex. Note
 ! that there yambo uses the plain X_mat*gamp residual (no RIM-W) and ppaE/Re[ppaE/Omega] as
 ! pole: both are reproduced here from W_EPSm1 and W_PPA_E.
 !
 use pars,          ONLY:SP,DP,pi,schlen,cZERO,cONE,rZERO
 use units,         ONLY:HA2EV
 use stderr,        ONLY:intc
 use frequency,     ONLY:w_samp
 use electrons,     ONLY:levels,spin_occ,spin,n_sp_pol
 use LIVE_t,        ONLY:live_timing
 use com,           ONLY:msg
 use drivers,       ONLY:l_ppa,l_cohsex,Finite_Tel,l_rt_carriers_in_use
 use D_lattice,     ONLY:nsym,i_time_rev,i_space_inv,mag_syms,Bose_Temp
 use functions,     ONLY:bose_E_cut
 use parallel_int,  ONLY:PP_wait,PP_redux_wait,PARALLEL_WF_index,PARALLEL_WF_distribute
 use parallel_m,    ONLY:PAR_IND_Xk_ibz,PAR_IND_G_b,PAR_IND_QP,PAR_IND_Q_ibz,PAR_IND_Q_ibz_ID,PAR_IND_QP_ID,&
&                        n_WF_bands_to_load,PAR_COM_G_b_INDEX,master_cpu
 use collision_el,  ONLY:elemental_collision,elemental_collision_free,elemental_collision_alloc
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:manage_action,OP_RD_CL,REP,VERIFY,NONE,OP_RD,RD_CL_IF_END,deliver_IO_error_message
 use QP_m,          ONLY:QP_t,QP_n_G_bands,QP_table,QP_Sc,QP_n_states,QP_G_damp,COHSEX_use_empties,&
&                        l_GW_terminator,l_QP_symmetrize,QP_states_simmetrize
 use X_m,           ONLY:X_t
 use wave_func,     ONLY:WF
 use R_lattice,     ONLY:qindx_S,bz_samp,G_m_G,G_m_G_maxval,d3q_factor
 use Wpot_m,        ONLY:W_is_PPA,W_ng,W_n_freqs,W_screening,W_free,W_vcoul,W_static,W_ppa_E,W_epsm1,W_RIM_W
 use interfaces,    ONLY:WF_load,WF_free
 use timing_m,      ONLY:timing
 !
#include<y_memory.h>
 !
 type(levels) ::E
 type(bz_samp)::Xk,k,q
 type(X_t)    ::X
 type(QP_t)   ::qp
 type(w_samp) ::Xw
 !
 ! Work Space
 !
 type(elemental_collision), target :: isc
 integer          :: i_qp,iqbz,iqibz,iqs,ib,ig1,ig2,ig_mG,ng,io_err,ID,IO_ACT,is,iq_loaded,timing_steps
 integer          :: os_save(4),qs_save(3)
 logical          :: l_TR,W_is_T,l_COH_closure,l_term
 real(SP)         :: q_weight,gamp_fac,E_kmq,f_kmq,PPA_E,bose_PPA_E,eet_cutoff0,eet_cutoff1,&
&                    eet_cutoff0_sp(n_sp_pol),eet_cutoff1_sp(n_sp_pol),max_dev,dev,X_ppaE
 complex(SP)      :: W_i,QP_ppa_EET,pre_factor
 complex(DP)      :: dp_dummy,ctmp
 character(schlen):: ch
 !
 real(SP),    allocatable :: Om(:,:),Om_t(:,:)
 complex(SP), allocatable :: Wc(:,:),Wc_t(:,:),Sc_W(:),eet_factor(:,:),local_rhotw(:)
 !
 integer,     external :: io_W
 complex(SP), external :: QP_ppa_EET_terminator
 !
 if (.not.(l_ppa.or.l_cohsex)) then
   call warning('[W-check] Sc from ndb.W is available for plasmon-pole and COHSEX runs only')
   return
 endif
 !
 call section('=','Correlation self-energy from ndb.W (check)')
 !===============================================================
 !
 call timing('W_Sigma_check',OPR='start')
 !
 ng      =X%ng
 X_ppaE  =X%ppaE
 q_weight=d3q_factor/(2._SP*pi)**3
 gamp_fac=4._SP*pi/q_weight
 !
 l_COH_closure=l_cohsex.and..not.COHSEX_use_empties
 l_term       =l_ppa.and.l_GW_terminator
 !
 ! ndb.W check (VERIFY mode): the DB must contain at least X%ng RL vectors
 !
 W_is_PPA   =l_ppa
 W_n_freqs  =Xw%n_freqs
 W_ng       =ng
 W_screening="em1s"
 if (l_ppa) W_screening="pp"
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=ID)
 io_err=io_W(X,Xw,ID)
 if (io_err<0) call error('[W-check] ndb.W not found or not compatible with the present run')
 !
 call msg('r','[W-check] Screening database     ',trim(W_screening))
 call msg('r','[W-check] Matrix size            ',ng)
 call msg('r','[W-check] Bands range            ',QP_n_G_bands)
 call msg('r','[W-check] RIM-W corrections in W ',W_RIM_W)
 if (l_ppa) call msg('r','[W-check] GW terminator          ',l_GW_terminator)
 if (l_cohsex) call msg('r','[W-check] COH from empty states  ',COHSEX_use_empties)
 !
 ! WF distribution and load, as in QP_ppa_cohsex (the WFs were freed at its end)
 !
 isc%ngrho=ng
 if (l_COH_closure.or.l_term) then
   if (.not.allocated(G_m_G)) call error('[W-check] G-G table not available')
   isc%ngrho=G_m_G_maxval
 endif
 !
 call PARALLEL_WF_distribute(K_index=PAR_IND_Xk_ibz,B_index=PAR_IND_G_b,CLEAN_UP=.TRUE.)
 call PARALLEL_WF_distribute(QP_index=PAR_IND_QP)
 call PARALLEL_WF_index( )
 !
 call WF_load(WF,isc%ngrho,maxval(qindx_S(:,:,2)),(/1,n_WF_bands_to_load/),(/1,k%nibz/),title='-W-check')
 !
 ! Allocations
 !
 call elemental_collision_free(isc)
 isc%iqref=0
 call elemental_collision_alloc(isc,NG=isc%ngrho,TITLE="W-check")
 !
 YAMBO_ALLOC(Wc,(ng,ng))
 YAMBO_ALLOC(Sc_W,(QP_n_states))
 YAMBO_ALLOC(local_rhotw,(ng))
 Sc_W=cZERO
 if (l_ppa) then
   YAMBO_ALLOC(Om,(ng,ng))
 endif
 if (l_term) then
   YAMBO_ALLOC(Wc_t,(ng,ng))
   YAMBO_ALLOC(Om_t,(ng,ng))
   YAMBO_ALLOC(eet_factor,(ng,ng))
   do is=1,n_sp_pol
     eet_cutoff0_sp(is)=minval(E%E(E%nbf(is)+1,:,is))
     eet_cutoff1_sp(is)=minval(E%E(QP_n_G_bands(2),:,is))
   enddo
   eet_cutoff0=minval(eet_cutoff0_sp(:))
   eet_cutoff1=minval(eet_cutoff1_sp(:))
 endif
 !
 ! ndb.W header (kept open during the q-loop)
 !
 call io_control(ACTION=OP_RD,COM=NONE,SEC=(/1/),ID=ID)
 io_err=io_W(X,Xw,ID)
 !
 timing_steps=PAR_IND_QP%n_of_elements(PAR_IND_QP_ID+1)*&
&             PAR_IND_Q_ibz%n_of_elements(PAR_IND_Q_ibz_ID+1)*&
&             count( PAR_IND_G_b%element_1D(QP_n_G_bands(1):QP_n_G_bands(2)) )
 call live_timing('[W-check] Sc from ndb.W',timing_steps)
 !
 iq_loaded=0
 W_is_T   =.false.
 !
 Q_loop: do iqbz=1,q%nbz
   !
   iqibz=q%sstar(iqbz,1)
   iqs  =q%sstar(iqbz,2)
   !
   if (.not.PAR_IND_Q_ibz%element_1D(iqibz)) cycle
   !
   isc%qs(2:)=(/iqibz,iqs/)
   !
   if (iqibz/=iq_loaded) then
     !
     ! I/O
     !
     IO_ACT=manage_action(RD_CL_IF_END,iqibz,1,q%nibz,PAR_INDEX=PAR_IND_Q_ibz)
     call io_control(ACTION=IO_ACT,COM=NONE,SEC=(/2*iqibz,2*iqibz+1/),ID=ID)
     io_err=io_W(X,Xw,ID)
     call deliver_IO_error_message(io_err,'W',STOP_it=.TRUE.)
     !
     if (size(W_static,1)<ng) call error('[W-check] ndb.W matrix smaller than the response block size')
     !
     ! Wc = W - v  and  Re[Omega]
     !
     Wc(:,:)=W_static(:ng,:ng)
     forall(ig1=1:ng) Wc(ig1,ig1)=Wc(ig1,ig1)-W_vcoul(ig1)
     !
     if (l_ppa) Om(:,:)=real(W_ppa_E(:ng,:ng),SP)
     !
     if (l_term) then
       !
       ! Terminator: plain X_mat*gamp residual, gamp rebuilt from scatter_Gamp, and
       ! pole ppaE/Re[ppaE/Omega], as in QP_ppa_cohsex
       !
       call elemental_collision_alloc(isc,GAMP_NG=(/ng,ng/))
       call scatter_Gamp(isc,'c')
       !$omp parallel do default(shared), private(ig1,ig2), collapse(2)
       do ig2=1,ng
         do ig1=1,ng
           Wc_t(ig1,ig2)=(W_epsm1(ig1,ig2)-cONE*merge(1,0,ig1==ig2))*isc%gamp(ig1,ig2)*gamp_fac
           Om_t(ig1,ig2)=X_ppaE/real(X_ppaE/W_ppa_E(ig1,ig2),SP)
         enddo
       enddo
       !$omp end parallel do
       YAMBO_FREE(isc%gamp)
       !
     endif
     !
     iq_loaded=iqibz
     W_is_T   =.false.
     !
   endif
   !
   ! Time-reversal without spatial inversion: W(-q,G,G') = W(q,-G',-G)
   !
   l_TR=(iqs>nsym/(i_time_rev+1)).and.(i_space_inv==0.or.mag_syms)
   if (l_TR.neqv.W_is_T) then
     Wc=transpose(Wc)
     if (l_ppa)  Om=transpose(Om)
     if (l_term) then
       Wc_t=transpose(Wc_t)
       Om_t=transpose(Om_t)
     endif
     W_is_T=l_TR
   endif
   !
   QP_loop: do i_qp=1,QP_n_states
     !
     if (.not.PAR_IND_QP%element_1D(i_qp)) cycle
     !
     isc%is=(/QP_table(i_qp,1),QP_table(i_qp,3),1,spin(QP_table(i_qp,:))/)
     isc%os(2:)=(/k%sstar(qindx_S(isc%is(2),iqbz,1),:),spin(QP_table(i_qp,:))/)
     isc%qs(1)=qindx_S(QP_table(i_qp,3),iqbz,2)
     !
     W_i=cmplx(qp%E_bare(i_qp),QP_G_damp,kind=SP)
     !
     ! COH (closure relation) and terminator initialization: rho_nn'(G-G') at q=0
     !
     if (l_COH_closure.or.l_term) then
       !
       os_save=isc%os
       qs_save=isc%qs
       isc%os=(/QP_table(i_qp,2),QP_table(i_qp,3),1,spin(QP_table(i_qp,:))/)
       isc%qs=(/1,1,1/)
       !
       call scatter_Bamp(isc)
       !
       if (l_term) eet_factor=cZERO
       !
       if (PAR_COM_G_b_INDEX%CPU_id==0) then
         dp_dummy=0._DP
         do ig2=1,ng
           do ig1=1,ng
             if (W_is_T) then
               ig_mG=G_m_G(ig2,ig1)
             else
               ig_mG=G_m_G(ig1,ig2)
             endif
             if (l_term)        eet_factor(ig1,ig2)=isc%rhotw(ig_mG)
             if (l_COH_closure) dp_dummy=dp_dummy+cmplx(isc%rhotw(ig_mG)*Wc(ig1,ig2),kind=DP)
           enddo
         enddo
         if (l_COH_closure) Sc_W(i_qp)=Sc_W(i_qp)+cmplx(q_weight/2._SP*dp_dummy,kind=SP)
       endif
       !
       isc%os=os_save
       isc%qs=qs_save
       !
     endif
     !
     do ib=QP_n_G_bands(1),QP_n_G_bands(2)
       !
       if (.not.PAR_IND_G_b%element_1D(ib)) cycle
       !
       if (iqs==1) call live_timing(steps=1)
       !
       isc%os(1)=ib
       call scatter_Bamp(isc)
       !
       E_kmq=E%E(isc%os(1),isc%os(2),isc%os(4))
       f_kmq=E%f(isc%os(1),isc%os(2),isc%os(4))
       !
       if (l_ppa) then
         !
         if (l_term) then
           !$omp parallel do default(shared), private(ig1,ig2), collapse(2)
           do ig2=1,ng
             do ig1=1,ng
               eet_factor(ig1,ig2)=eet_factor(ig1,ig2)-isc%rhotw(ig1)*conjg(isc%rhotw(ig2))
             enddo
           enddo
           !$omp end parallel do
         endif
         !
         dp_dummy=0._DP
         !
         !$omp parallel do default(shared), private(ig1,ig2,PPA_E,bose_PPA_E,ctmp), &
         !$omp &           reduction(+:dp_dummy), collapse(2)
         do ig2=1,ng
           do ig1=1,ng
             !
             PPA_E=Om(ig1,ig2)
             !
             bose_PPA_E=0._SP
             if (PPA_E<0._SP) bose_PPA_E=-spin_occ
             if (Finite_Tel) then
               if (abs(PPA_E)>epsilon(1._SP)) then
                 if (abs(PPA_E)<=bose_E_cut*Bose_Temp) bose_PPA_E=spin_occ*Bose_Temp/PPA_E
                 if (abs(PPA_E)> bose_E_cut*Bose_Temp) bose_PPA_E=spin_occ/(exp(PPA_E/Bose_Temp)-1._SP)
               else
                 bose_PPA_E=spin_occ*Bose_Temp/epsilon(1._SP)
               endif
             endif
             !
             ctmp=cmplx(-q_weight/spin_occ*isc%rhotw(ig1)*conjg(isc%rhotw(ig2))*Wc(ig1,ig2)*PPA_E/2._SP*&
&                       ( (spin_occ-f_kmq+bose_PPA_E)/(W_i-E_kmq-PPA_E)+&
&                         (f_kmq+bose_PPA_E)/(conjg(W_i)-E_kmq+PPA_E) ),kind=DP)
             !
             dp_dummy=dp_dummy+ctmp
             !
           enddo
         enddo
         !$omp end parallel do
         !
         Sc_W(i_qp)=Sc_W(i_qp)+cmplx(dp_dummy,kind=SP)
         !
       else
         !
         ! COHSEX: pre_factor = rho . Wc . rho^*  (no conjugation on the first rho, as V_dot_V)
         !
         !$omp parallel do default(shared), private(ig1,ig2)
         do ig1=1,ng
           local_rhotw(ig1)=cZERO
           do ig2=1,ng
             local_rhotw(ig1)=local_rhotw(ig1)+Wc(ig1,ig2)*conjg(isc%rhotw(ig2))
           enddo
         enddo
         !$omp end parallel do
         !
         pre_factor=cZERO
         do ig1=1,ng
           pre_factor=pre_factor+isc%rhotw(ig1)*local_rhotw(ig1)
         enddo
         !
         ! SEX
         !
         Sc_W(i_qp)=Sc_W(i_qp)-q_weight/spin_occ*pre_factor*f_kmq
         !
         ! COH (with empty states)
         !
         if (COHSEX_use_empties) Sc_W(i_qp)=Sc_W(i_qp)+q_weight/2._SP*pre_factor
         !
       endif
       !
     enddo
     !
     ! Terminator
     !
     if (l_term) then
       !
       dp_dummy=0._DP
       !
       do ig2=1,ng
         do ig1=1,ng
           !
           PPA_E=Om_t(ig1,ig2)
           QP_ppa_EET=QP_ppa_EET_terminator(W_i,E,isc%is,PPA_E,ig1,ig2,iqibz,eet_cutoff0,eet_cutoff1)
           !
           ! 4 pi/s * PPA_R * eet_factor * EET, PPA_R = -X gamp Om/2 = -(q_w/4 pi) Wc_t Om/2
           !
           ctmp=cmplx(-q_weight/spin_occ*Wc_t(ig1,ig2)*PPA_E/2._SP*eet_factor(ig1,ig2)*QP_ppa_EET,kind=DP)
           dp_dummy=dp_dummy+ctmp
           !
         enddo
       enddo
       !
       Sc_W(i_qp)=Sc_W(i_qp)+cmplx(dp_dummy,kind=SP)
       !
     endif
     !
   enddo QP_loop
   !
 enddo Q_loop
 !
 call live_timing()
 !
 ! Clean
 !
 call WF_free(WF)
 call elemental_collision_free(isc)
 call W_free()
 YAMBO_FREE(Wc)
 YAMBO_FREE(Wc_t)
 YAMBO_FREE(Om)
 YAMBO_FREE(Om_t)
 YAMBO_FREE(eet_factor)
 YAMBO_FREE(local_rhotw)
 !
 ! ALL 2 ALL
 !
 call PP_redux_wait(Sc_W)
 !
 ! Same symmetrization of degenerate states applied by QP_ppa_cohsex
 !
 if (.not.l_rt_carriers_in_use.and.l_QP_symmetrize) then
   call QP_states_simmetrize(E,what="Sc from ndb.W",V_complex=Sc_W,warn_me=.FALSE.)
 endif
 !
 ! Report
 !
 max_dev=rZERO
 do i_qp=1,QP_n_states
   max_dev=max(max_dev,abs(Sc_W(i_qp)-QP_Sc(i_qp,1)))
 enddo
 !
 call msg('r','')
 call msg('r','[W-check] Sc(n,k;Eo) [eV]: QP_ppa_cohsex (ndb.pp/em1s) vs ndb.W')
 call msg('r','[W-check] Max |Sc(ndb.W)-Sc(yambo)| [eV]',max_dev*HA2EV)
 call msg('r','')
 !
 if (master_cpu) then
   !
   ch='Sc_W'
   call msg('o '//trim(ch),'# Correlation self-energy Sc(n,k;Eo) from ndb.W vs QP_ppa_cohsex')
   call msg('o '//trim(ch),'#',(/' K-point','    Band','  Eo[eV]',' ReSc[Y]',' ImSc[Y]',' ReSc[W]',' ImSc[W]','   |dSc|'/),&
&           INDENT=0,USE_TABS=.true.)
   !
   do i_qp=1,QP_n_states
     dev=abs(Sc_W(i_qp)-QP_Sc(i_qp,1))*HA2EV
     call msg('o '//trim(ch),'#',(/real(QP_table(i_qp,3),SP),real(QP_table(i_qp,1),SP),qp%E_bare(i_qp)*HA2EV,&
&             real(QP_Sc(i_qp,1))*HA2EV,aimag(QP_Sc(i_qp,1))*HA2EV,real(Sc_W(i_qp))*HA2EV,aimag(Sc_W(i_qp))*HA2EV,dev/),&
&             INDENT=0,USE_TABS=.true.)
     call msg('r',' K B Sc[Y] Sc[W] |dSc| ',(/real(QP_table(i_qp,3),SP),real(QP_table(i_qp,1),SP),&
&             real(QP_Sc(i_qp,1))*HA2EV,real(Sc_W(i_qp))*HA2EV,dev/))
   enddo
   !
 endif
 !
 YAMBO_FREE(Sc_W)
 !
 call PP_wait()
 !
 call timing('W_Sigma_check',OPR='stop')
 !
end subroutine W_Sigma_check
