!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details): 
!
integer function io_W(X,Xw,ID)
 !
 ! Screened Coulomb potential W(q,G,G') database: ndb.W (header) + ndb.W_fragment_<iq> (per q-point).
 ! See mod_Wpot for the conventions.
 !
 ! Sections:  1        -> header, global descriptors, G-vectors and |q+G|
 !            2*iq     -> per-q descriptors (plasmon-pole statistics)
 !            2*iq+1   -> per-q matrices
 !
 use pars,          ONLY:SP,schlen
 use R_lattice,     ONLY:g_vec,bare_qpg,nqibz
 use X_m,           ONLY:X_t
 use frequency,     ONLY:w_samp
 use Wpot_m,        ONLY:W_desc,W_is_PPA,W_ng,W_n_freqs,W_alloc,W_freqs,W_vcoul,W_epsm1,W_static,W_ppa_E,&
&                        W_PPcond_rate,W_TO_rate,W_PP_err
 use IO_int,        ONLY:io_connect,io_disconnect,io_elemental,io_header,io_fragment,&
&                        def_variable_bulk,io_variable_bulk,def_variable_elemental,io_variable_elemental
 use IO_m,          ONLY:io_sec,io_com,io_status,write_is_on,read_is_on,NONE,frag_SCREEN,&
&                        IO_NO_DATABASE,io_mode,DUMP
 use stderr,        ONLY:intc
 use timing_m,      ONLY:timing
 !
#include<y_memory.h>
 !
 type(X_t)   :: X
 type(w_samp):: Xw
 integer     :: ID
 !
 ! Work Space
 !
 integer                 :: i1,iq,io_com_save,io_err,ID_frag,ng_qpg
 character(schlen)       :: ch
 real(SP),   allocatable :: RL_vecs_disk(:,:)
 complex(SP),allocatable :: qpg_disk(:,:)
 integer,    external    :: variables_W
 !
 call timing('io_W',OPR='start')
 !
 ID_frag=0
 io_err =0
 !
 iq=0
 do i1=1,size(io_sec(ID,:))
   if (io_sec(ID,i1)==0) cycle
   iq=io_sec(ID,i1)/2
 enddo
 !
 io_W=io_connect(desc="W",type=2,ID=ID,FORCE_READ_MODE=(.not.any((/io_sec(ID,:)==1/))).and.frag_SCREEN)
 if (io_W/=0) goto 1
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   ! Header
   !
   io_W=io_header(ID,QPTS=.true.,R_LATT=.true.,WF=.false.,CUTOFF=.true.,FRAG=frag_SCREEN)
   if (io_W/=0) goto 1
   !
   ! Global Descriptions
   !=====================
   !
   io_W=variables_W(ID,W_desc,X)
   if (io_W/=0) goto 1
   !
   ! G-vectors, in the ordering used for the W matrices
   !
   YAMBO_ALLOC(RL_vecs_disk,(W_ng,3))
   if (write_is_on(ID)) RL_vecs_disk(:,:)=g_vec(:W_ng,:)
   call def_variable_bulk(ID,'W_RL_vecs',2,(/W_ng,3/),SP)
   call  io_variable_bulk(ID,2,R2=RL_vecs_disk)
   YAMBO_FREE(RL_vecs_disk)
   !
   ! |q+G| (bare, before any regularization) for all the q-points in the IBZ
   !
   YAMBO_ALLOC(qpg_disk,(nqibz,W_ng))
   qpg_disk=(0._SP,0._SP)
   if (write_is_on(ID)) then
     ng_qpg=min(W_ng,size(bare_qpg,2))
     qpg_disk(:,:ng_qpg)=bare_qpg(:nqibz,:ng_qpg)
   endif
   call def_variable_bulk(ID,'W_QPG',3,(/2,nqibz,W_ng/),SP)
   call  io_variable_bulk(ID,3,C2=qpg_disk)
   YAMBO_FREE(qpg_disk)
   !
   ! Close & Menu lower line
   !
   if (.not.any((/io_sec(ID,:)==2/))) then
     call io_elemental(ID,VAR="",VAR_SZ=0,MENU=1)
     io_W=io_status(ID)
     if (io_W/=0) goto 1
   else
     call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   endif
   !
 endif
 !
 if (iq==0) goto 1
 !
 call io_fragment(ID,ID_frag,i_fragment=iq,ierr=io_err)
 io_W=io_err
 if (io_err/=0) goto 1
 !
 if (any((/io_sec(ID,:)==2*iq/))) then
   !
   ch="W_PARS_sec_iq"//trim(intc(iq))
   !
   io_com_save=io_com(ID_frag)
   if (iq>2) io_com(ID_frag)=NONE
   !
   ! 5 slots: iq, n_freqs, PP condition fails, time-ordering fails, PP mean relative deviation
   !
   call def_variable_elemental(ID_frag,trim(ch),5,SP,1)
   call io_variable_elemental(ID_frag,VAR="Current Q-pt index      ",I0=iq,CHECK=.true.,OP=(/"=="/))
   call io_variable_elemental(ID_frag,VAR="Number of frequencies   ",I0=W_n_freqs,CHECK=.true.,OP=(/"=="/))
   call io_variable_elemental(ID_frag,VAR="PP condition fails      ",I0=W_PPcond_rate)
   call io_variable_elemental(ID_frag,VAR="Time ordering fails     ",I0=W_TO_rate)
   call io_variable_elemental(ID_frag,VAR="PP mean rel deviation   ",R0=W_PP_err)
   call def_variable_elemental(ID_frag,"",0,0,1)
   io_W=io_status(ID_frag)
   !
   io_com(ID_frag)=io_com_save
   !
   if (io_W/=0) goto 1
   !
 endif
 !
 if (any((/io_sec(ID,:)==2*iq+1/))) then
   !
   if (read_is_on(ID)) call W_alloc(W_ng,W_n_freqs,W_is_PPA)
   !
   ! Frequencies of the screening DB
   !
   ch="W_FREQ_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),1,(/2,W_n_freqs/),SP)
   call  io_variable_bulk(ID_frag,1,C1=W_freqs)
   !
   ! Coulomb potential v(q,G) (diagonal), as used by yambo (RIM regularized when RIM is on)
   !
   ch="W_VCOUL_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),2,(/2,W_ng/),SP)
   call  io_variable_bulk(ID_frag,2,C1=W_vcoul)
   !
   ! Symmetrized inverse dielectric matrix eps^-1_sym = 1 + X_mat (dimensionless)
   !
   ch="W_EPSm1_STATIC_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),3,(/2,W_ng,W_ng/),SP)
   call  io_variable_bulk(ID_frag,3,C2=W_epsm1)
   !
   ! W(q,G,G';w=0) = v + v chi v [Ha]
   !
   ch="W_STATIC_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),4,(/2,W_ng,W_ng/),SP)
   call  io_variable_bulk(ID_frag,4,C2=W_static)
   !
   ! Plasmon-pole energies Omega(q,G,G') [Ha]
   !
   if (W_is_PPA) then
     ch="W_PPA_E_Q_"//trim(intc(iq))
     call def_variable_bulk(ID_frag,trim(ch),5,(/2,W_ng,W_ng/),SP)
     call  io_variable_bulk(ID_frag,5,C2=W_ppa_E)
   endif
   !
 endif
 !
 ! Sections 1 & 2 are used in VERIFY mode.
 ! Now that the menu is closed I can return with io_W/=0 if there is a restart point
 !
 if (io_mode(ID)/=DUMP) then
   io_W=IO_NO_DATABASE
   if (io_err==0) io_W=iq
 endif
 !
1 call io_fragment_disconnect(ID,ID_frag)
 call io_disconnect(ID)
 !
 call timing('io_W',OPR='stop')
 !
end function
!
integer function variables_W(ID,local_desc,X)
 !
 use pars,          ONLY:SP
 use units,         ONLY:HA2EV
 use descriptors,   ONLY:IO_desc,IO_desc_reset
 use X_m,           ONLY:X_t
 use R_lattice,     ONLY:RIM_qpg
 use Wpot_m,        ONLY:W_is_PPA,W_ng,W_n_freqs,W_screening,W_convention,W_RIM_W,W_RIM_W_ng
 use IO_int,        ONLY:io_elemental
 use IO_m,          ONLY:io_status
 !
 implicit none
 !
 integer            ::ID
 type(IO_desc)      ::local_desc
 type(X_t)          ::X
 !
 ! Work Space
 !
 integer            ::tmp_ib(2)
 logical            ::l_RIM_head
 !
 call IO_desc_reset(local_desc)
 !
 l_RIM_head=allocated(RIM_qpg)
 tmp_ib=X%ib
 !
 ! W_PARS_1: ng(1) + n_freqs(1) + bands(2) + ppaE(1) + PPA(1) + RIM(1) + RIM-W(1) + RIM-W ng(1) = 9 slots
 !
 call io_elemental(ID,VAR="W_PARS_1",VAR_SZ=9,MENU=0)
 call io_elemental(ID,I0=W_ng,&
&     VAR="W matrix size                   ",CHECK=.true.,OP=(/"<="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,I0=W_n_freqs,&
&     VAR="Screening frequencies           ",CHECK=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,I1=tmp_ib,&
&     VAR="Screening band range            ",CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,UNIT=HA2EV,R0=X%ppaE,&
&     VAR="PPA imaginary energy            ",TERMINATOR="eV",CHECK=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,L0=W_is_PPA,&
&     VAR="Plasmon-pole energies included  ",CHECK=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,L0=l_RIM_head,&
&     VAR="RIM regularized q->0 head       ",WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,L0=W_RIM_W,&
&     VAR="RIM-W corrections included      ",WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,I0=W_RIM_W_ng,&
&     VAR="RIM-W RL components             ",WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
 !
 ! W_SCREENING: the screening DB W has been built from
 !
 call io_elemental(ID,VAR="W_SCREENING",CH0='',VAR_SZ=1,MENU=0)
 call io_elemental(ID,CH0=W_screening,&
&     VAR="Screening database              ",CHECK=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
 !
 ! W_CONVENTION
 !
 call io_elemental(ID,VAR="W_CONVENTION",CH0='',VAR_SZ=1,MENU=0)
 call io_elemental(ID,CH0=W_convention,&
&     VAR="Convention                      ",DESCRIPTOR=local_desc)
 call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
 !
 variables_W=io_status(ID)
 !
end function
