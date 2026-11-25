!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Stub modules for compiling src/qp/QP_ppa_cohsex.F outside Yambo (test harness only)
module pars
 integer,parameter::SP=kind(1.0),DP=kind(1.0d0),schlen=100
 real(SP),parameter::pi=3.14159265358979323846_SP
 complex(SP),parameter::cZERO=(0._SP,0._SP),cI=(0._SP,1._SP),cONE=(1._SP,0._SP)
end module
module units
 use pars; real(SP),parameter::HA2EV=27.211386245988_SP
end module
module stderr
 contains
 character(10) function intc(i); integer,intent(in)::i; write(intc,'(i0)') i; end function
end module
module frequency
 use pars
 type w_samp
   integer::n_freqs=0
   complex(SP),allocatable::p(:)
   real(SP)::er(2)=0,dr(2)=0
   character(2)::grid_type='ra'
 end type
 contains
 subroutine W_reset(W); type(w_samp)::W; if(allocated(W%p)) deallocate(W%p); W%n_freqs=0; end subroutine
end module
module electrons
 use pars
 type levels
   integer::nb=0,nbf(1)=13
   real(SP),allocatable::E(:,:,:),f(:,:,:)
 end type
 real(SP)::spin_occ=2._SP
 integer::n_met_bands(1)=0,n_sp_pol=1
 contains
 integer function spin(v); integer,intent(in)::v(:); spin=1; end function
end module
module LIVE_t
 contains
 subroutine live_timing(message,steps); character(*),optional::message; integer,optional::steps; end subroutine
end module
module com
 use pars
 interface msg
   module procedure msg_s,msg_i,msg_iv,msg_r
 end interface
 contains
 subroutine msg_s(m,s); character(*)::m,s; end subroutine
 subroutine msg_i(m,s,i); character(*)::m,s; integer::i; end subroutine
 subroutine msg_iv(m,s,i); character(*)::m,s; integer::i(:); end subroutine
 subroutine msg_r(m,s,r,u); character(*)::m,s; real(SP)::r; character(*),optional::u; end subroutine
end module
module drivers
 logical::l_ppa=.true.,l_cohsex=.false.,l_sc_srpa=.false.,l_sc_sex=.false.,l_sc_coh=.false.,l_sc_run=.false.
 logical::l_rt_carriers_in_use=.false.,l_rim_w=.true.,Finite_Tel=.false.
end module
module parallel_m
 type PP_index
   logical,allocatable::element_1D(:)
   integer,allocatable::n_of_elements(:)
 end type
 type PP_comm
   integer::CPU_id=0,n_CPU=1
 end type
 type(PP_index)::PAR_IND_Xk_ibz,PAR_IND_G_b,PAR_IND_QP,PAR_IND_Q_ibz
 integer::PAR_IND_Q_ibz_ID=0,PAR_IND_QP_ID=0,n_WF_bands_to_load=0
 integer,allocatable::PAR_Q_ibz_index(:)
 logical::HEAD_QP_cpu=.true.
 type(PP_comm)::PAR_COM_RL_A2A,PAR_COM_QP_A2A,PAR_COM_G_b_INDEX
end module
module R_lattice
 use pars
 type bz_samp
   integer::nbz=0,nibz=0
   integer,allocatable::sstar(:,:)
 end type
 integer,allocatable::qindx_S(:,:,:),G_m_G(:,:)
 integer::G_m_G_maxval=0,nqibz=0,RIM_W_ng=0
 logical::RIM_W_is_diagonal=.false.
 real(SP),allocatable::RIM_W_E(:,:),RIM_W(:,:,:,:)
end module
module collision_el
 use pars
 type elemental_collision
   integer::ngrho=0,iqref=0,is(4)=0,os(4)=0,qs(3)=0
   complex(SP),allocatable::rhotw(:),gamp(:,:)
 end type
 contains
 subroutine elemental_collision_free(c); type(elemental_collision)::c
   if(allocated(c%rhotw)) deallocate(c%rhotw); if(allocated(c%gamp)) deallocate(c%gamp)
 end subroutine
 subroutine elemental_collision_alloc(c,NG,TITLE,GAMP_ROWS,GAMP_COLS)
   type(elemental_collision)::c; integer,optional::NG,GAMP_ROWS(2),GAMP_COLS(2); character(*),optional::TITLE
   if(present(NG)) then; if(allocated(c%rhotw)) deallocate(c%rhotw); allocate(c%rhotw(NG)); endif
   if(present(GAMP_ROWS)) then; if(.not.allocated(c%gamp)) allocate(c%gamp(GAMP_ROWS(1):GAMP_ROWS(2),GAMP_COLS(1):GAMP_COLS(2))); endif
 end subroutine
end module
module functions
 use pars; real(SP)::bose_E_cut=0.1_SP
 contains
 real(SP) function bose_f(e); real(SP)::e; bose_f=0; end function
end module
module IO_int
 contains
 subroutine io_control(ACTION,COM,SEC,MODE,ID); integer,optional::ACTION,COM,SEC(:),MODE,ID; end subroutine
end module
module IO_m
 integer,parameter::OP_RD_CL=1,REP=2,VERIFY=3,NONE=0,OP_RD=4,RD_CL_IF_END=5
 logical::io_RESPONSE=.false.
 contains
 integer function manage_action(a,i,i0,n,PAR_INDEX); use parallel_m; integer::a,i,i0,n; type(PP_index),optional::PAR_INDEX; manage_action=a; end function
 subroutine deliver_IO_error_message(e,w,STOP_it); integer::e; character(*)::w; logical,optional::STOP_it; end subroutine
end module
module QP_m
 use pars
 use electrons
 type QP_t
   integer::n_states=0
   real(SP),allocatable::E_bare(:)
 end type
 character(1)::QP_solver='g'
 integer::QP_n_G_bands(2)=0,QP_dSc_steps=0,QP_Sc_steps=0,QP_n_states=0
 integer,allocatable::QP_table(:,:)
 logical::l_QP_symmetrize=.false.,COHSEX_use_empties=.false.,l_GW_terminator=.true.,QP_retarded_G=.true.
 complex(SP),allocatable::QP_Sc(:,:)
 real(SP)::QP_G_damp=0,QP_dSc_delta=0.1_SP/27.2114_SP,GW_terminator_E=1.5_SP
 character(10)::GW_terminator_Kind='BG'
 contains
 subroutine QP_states_simmetrize(E,what,V_complex,warn_me)
   type(levels)::E; character(*),optional::what; complex(SP),optional::V_complex(:); logical,optional::warn_me
 end subroutine
 subroutine QP_prepare_G_grid(p,grid_type)
   complex(SP),intent(inout)::p(:); character(*),intent(in)::grid_type
   QP_G_damp=0._SP
 end subroutine
end module
module matrix
 use pars
 type PAR_matrix
   complex(SP),allocatable::blc(:,:,:)
   integer::rows(2)=0,cols(2)=0
 end type
 contains
 subroutine MATRIX_reset(M); type(PAR_matrix)::M; end subroutine
end module
module ALLOC
 contains
 subroutine X_ALLOC_elemental(w); character(*)::w; end subroutine
end module
module X_m
 use matrix
 type X_t
   integer::ng=0
   real(SP)::ppaE=1._SP
 end type
 type(PAR_matrix),allocatable,target::X_par(:)
end module
module wave_func
 type WAVEs; integer::i=0; end type
 type(WAVEs)::WF
end module
module D_lattice
 use pars
 integer::nsym=48,i_time_rev=0,i_space_inv=1
 logical::mag_syms=.false.
 real(SP)::Bose_Temp=0
end module
module wrapper
 use pars
 contains
 complex(SP) function V_dot_V_gpu(n,a,b); integer::n; complex(SP)::a(:),b(:); V_dot_V_gpu=sum(a*b); end function
 subroutine devxlib_xGEMV_gpu(t,m,n,al,A,lda,x,ix,be,y,iy)
   character::t; integer::m,n,lda,ix,iy; complex(SP)::al,be,A(:,:),x(:),y(:); y=al*matmul(A,x)+be*y
 end subroutine
end module
module interfaces
 contains
 subroutine QP_state_print(); end subroutine
 subroutine WF_load(W,ng,nq,b,k,title); use wave_func; type(WAVEs)::W; integer::ng,nq,b(2),k(2); character(*),optional::title
 end subroutine
 subroutine WF_free(W); use wave_func; type(WAVEs)::W; end subroutine
 subroutine MATRIX_transpose(op,M_out); use matrix; character(*)::op; type(PAR_matrix)::M_out
   stop 'MATRIX_transpose not expected'
 end subroutine
end module
module timing_m
 contains
 subroutine timing(n,OPR); character(*)::n; character(*),optional::OPR; end subroutine
end module
module devxlib
 use pars
 interface devxlib_memcpy_d2d; module procedure cp1; end interface
 interface devxlib_memcpy_d2h; module procedure cp2; end interface
 interface devxlib_memcpy_h2d; module procedure cp3; end interface
 interface devxlib_memset_d;   module procedure set2; end interface
 interface devxlib_conjg_d;    module procedure cj1; end interface
 contains
 subroutine cp1(d,s); complex(SP)::d(:),s(:); d=s(1:size(d)); end subroutine
 subroutine cp2(d,s); complex(SP)::d(:,:),s(:,:); d=s; end subroutine
 subroutine cp3(d,s); complex(SP)::d(:,:,:),s(:,:,:); d=s; end subroutine
 subroutine set2(d,v); complex(SP)::d(:,:),v; d=v; end subroutine
 subroutine cj1(d); complex(SP)::d(:); d=conjg(d); end subroutine
end module
module gpu_m
 logical::have_gpu=.false.
end module
module parallel_int
 contains
 subroutine PP_wait(); end subroutine
 subroutine PP_redux_wait(a); complex(kind(1.0))::a(:,:); end subroutine
 subroutine PARALLEL_global_indexes(E,k,q,what)
   use electrons; use R_lattice; type(levels)::E; type(bz_samp)::k,q; character(*)::what
 end subroutine
 subroutine PARALLEL_WF_index(); end subroutine
 subroutine PARALLEL_WF_distribute(K_index,B_index,CLEAN_UP,QP_index)
   use parallel_m; type(PP_index),optional::K_index,B_index,QP_index; logical,optional::CLEAN_UP
 end subroutine
end module
