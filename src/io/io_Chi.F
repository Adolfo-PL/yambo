!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
integer function io_Chi(X,Xw,ID)
 !
 ! ndb.Chi: chi0 (Kohn-Sham), P (quasiparticle) and fxc=chi0^-1-P^-1 in the plane-wave
 ! basis. Header in ndb.Chi, per-q data in ndb.Chi_fragment_<iq> when the response
 ! databases are fragmented. See mod_Chi for the conventions.
 !
 ! Sections:  1        -> header, descriptors, G-vectors and |q+G|
 !            2*iq     -> per-q descriptors
 !            2*iq+1   -> per-q matrices
 !
 use pars,          ONLY:SP,schlen
 use units,         ONLY:HA2EV
 use R_lattice,     ONLY:g_vec,bare_qpg,nqibz
 use X_m,           ONLY:X_t
 use frequency,     ONLY:w_samp
 use Chi_m,         ONLY:Chi_desc,Chi_ng,Chi_n_freqs,Chi_alloc,Chi_freqs,Chi_0,Chi_P,Chi_fxc,&
&                        Chi_QP_applied,Chi_QP_action,Chi_convention
 use IO_int,        ONLY:io_connect,io_disconnect,io_elemental,io_header,io_fragment,&
&                        def_variable_bulk,io_variable_bulk,def_variable_elemental,io_variable_elemental
 use IO_m,          ONLY:io_sec,io_com,io_status,write_is_on,read_is_on,NONE,frag_RESPONSE,&
&                        IO_NO_DATABASE,io_mode,DUMP
 use descriptors,   ONLY:IO_desc_reset
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
 integer                 :: i1,iq,io_com_save,io_err,ID_frag,ng_qpg,tmp_ib(2)
 character(schlen)       :: ch
 real(SP),   allocatable :: RL_vecs_disk(:,:)
 complex(SP),allocatable :: qpg_disk(:,:)
 !
 call timing('io_Chi',OPR='start')
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
 io_Chi=io_connect(desc="Chi",type=2,ID=ID,FORCE_READ_MODE=(.not.any((/io_sec(ID,:)==1/))).and.frag_RESPONSE)
 if (io_Chi/=0) goto 1
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   io_Chi=io_header(ID,QPTS=.true.,R_LATT=.true.,WF=.false.,CUTOFF=.true.,FRAG=frag_RESPONSE)
   if (io_Chi/=0) goto 1
   !
   call IO_desc_reset(Chi_desc)
   tmp_ib=X%ib
   !
   ! CHI_PARS_1: ng(1) + n_freqs(1) + bands(2) + QP applied(1) = 5 slots
   !
   call io_elemental(ID,VAR="CHI_PARS_1",VAR_SZ=5,MENU=0)
   call io_elemental(ID,I0=Chi_ng,&
&       VAR="Matrix size                     ",CHECK=.true.,OP=(/"<="/),DESCRIPTOR=Chi_desc)
   call io_elemental(ID,I0=Chi_n_freqs,&
&       VAR="Number of frequencies           ",CHECK=.true.,OP=(/"=="/),DESCRIPTOR=Chi_desc)
   call io_elemental(ID,I1=tmp_ib,&
&       VAR="Bands range                     ",CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=Chi_desc)
   call io_elemental(ID,L0=Chi_QP_applied,&
&       VAR="QP corrections applied to P     ",WARN=.true.,OP=(/"=="/),DESCRIPTOR=Chi_desc)
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   !
   call io_elemental(ID,VAR="CHI_QP_ACTION",CH0='',VAR_SZ=1,MENU=0)
   call io_elemental(ID,CH0=Chi_QP_action,&
&       VAR="QP action (XfnQPdb)             ",DESCRIPTOR=Chi_desc)
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   !
   call io_elemental(ID,VAR="CHI_CONVENTION",CH0='',VAR_SZ=1,MENU=0)
   call io_elemental(ID,CH0=Chi_convention,&
&       VAR="Convention                      ",DESCRIPTOR=Chi_desc)
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   !
   io_Chi=io_status(ID)
   if (io_Chi/=0) goto 1
   !
   ! G-vectors, in the ordering used for the matrices
   !
   YAMBO_ALLOC(RL_vecs_disk,(Chi_ng,3))
   if (write_is_on(ID)) RL_vecs_disk(:,:)=g_vec(:Chi_ng,:)
   call def_variable_bulk(ID,'CHI_RL_vecs',2,(/Chi_ng,3/),SP)
   call  io_variable_bulk(ID,2,R2=RL_vecs_disk)
   YAMBO_FREE(RL_vecs_disk)
   !
   ! |q+G| as used in the Coulomb potential (cutoff included when active)
   !
   YAMBO_ALLOC(qpg_disk,(nqibz,Chi_ng))
   qpg_disk=(0._SP,0._SP)
   if (write_is_on(ID)) then
     ng_qpg=min(Chi_ng,size(bare_qpg,2))
     qpg_disk(:,:ng_qpg)=bare_qpg(:nqibz,:ng_qpg)
   endif
   call def_variable_bulk(ID,'CHI_QPG',3,(/2,nqibz,Chi_ng/),SP)
   call  io_variable_bulk(ID,3,C2=qpg_disk)
   YAMBO_FREE(qpg_disk)
   !
   if (.not.any((/io_sec(ID,:)==2/))) then
     call io_elemental(ID,VAR="",VAR_SZ=0,MENU=1)
     io_Chi=io_status(ID)
     if (io_Chi/=0) goto 1
   else
     call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   endif
   !
 endif
 !
 if (iq==0) goto 1
 !
 call io_fragment(ID,ID_frag,i_fragment=iq,ierr=io_err)
 io_Chi=io_err
 if (io_err/=0) goto 1
 !
 if (any((/io_sec(ID,:)==2*iq/))) then
   !
   ch="CHI_PARS_sec_iq"//trim(intc(iq))
   !
   io_com_save=io_com(ID_frag)
   if (iq>2) io_com(ID_frag)=NONE
   !
   ! 6 slots: iq, n_freqs, energy range(2), damping range(2)
   !
   call def_variable_elemental(ID_frag,trim(ch),6,SP,1)
   call io_variable_elemental(ID_frag,VAR="Current Q-pt index      ",I0=iq,CHECK=.true.,OP=(/"=="/))
   call io_variable_elemental(ID_frag,VAR="Number of frequencies   ",I0=Chi_n_freqs,CHECK=.true.,OP=(/"=="/))
   call io_variable_elemental(ID_frag,UNIT=HA2EV,VAR="X energy range ",TERMINATOR="eV",R1=Xw%er,WARN=.true.,OP=(/"==","=="/))
   call io_variable_elemental(ID_frag,UNIT=HA2EV,VAR="X damping range",TERMINATOR="eV",R1=Xw%dr,WARN=.true.,OP=(/"==","=="/))
   call def_variable_elemental(ID_frag,"",0,0,1)
   io_Chi=io_status(ID_frag)
   !
   io_com(ID_frag)=io_com_save
   !
   if (io_Chi/=0) goto 1
   !
 endif
 !
 if (any((/io_sec(ID,:)==2*iq+1/))) then
   !
   if (read_is_on(ID)) call Chi_alloc(Chi_ng,Chi_n_freqs)
   !
   ch="CHI_FREQ_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),1,(/2,Chi_n_freqs/),SP)
   call  io_variable_bulk(ID_frag,1,C1=Chi_freqs)
   !
   ch="CHI0_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),2,(/2,Chi_ng,Chi_ng,Chi_n_freqs/),SP)
   call  io_variable_bulk(ID_frag,2,C3=Chi_0)
   !
   ch="P_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),3,(/2,Chi_ng,Chi_ng,Chi_n_freqs/),SP)
   call  io_variable_bulk(ID_frag,3,C3=Chi_P)
   !
   ch="FXC_Q_"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(ch),4,(/2,Chi_ng,Chi_ng,Chi_n_freqs/),SP)
   call  io_variable_bulk(ID_frag,4,C3=Chi_fxc)
   !
 endif
 !
 if (io_mode(ID)/=DUMP) then
   io_Chi=IO_NO_DATABASE
   if (io_err==0) io_Chi=iq
 endif
 !
1 call io_fragment_disconnect(ID,ID_frag)
 call io_disconnect(ID)
 !
 call timing('io_Chi',OPR='stop')
 !
end function
