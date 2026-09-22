! Small host-model fixtures. Production bubble, inversion and frequency
! switching are compiled unchanged; only the material/IO infrastructure is mocked.
module pars
 use iso_fortran_env,only:real64,real32
 implicit none
 integer,parameter::DP=real64,schlen=256
#ifdef _TEST_SINGLE
 integer,parameter::SP=real32
#else
 integer,parameter::SP=real64
#endif
 real(SP),parameter::pi=acos(-1._SP)
 complex(SP),parameter::cZERO=(0._SP,0._SP),cONE=(1._SP,0._SP)
 complex(SP),parameter::cI=(0._SP,1._SP)
 real(SP),parameter::rONE=1._SP,rZERO=0._SP
 real(SP),parameter::zero_dfl=1.E-5_SP
end module
module units
 use pars
 real(SP),parameter::HA2EV=27.211386_SP
 real(SP),parameter::HBAR_eVfs=.6582119569_SP
end module
module descriptors
 type IO_desc
   integer::n=0
 end type
contains
 subroutine IO_desc_reset(d)
   type(IO_desc)::d
   d%n=0
 end subroutine
 subroutine IO_desc_duplicate(a,b)
   type(IO_desc)::a,b
   b=a
 end subroutine
end module
module electrons
 use pars
 integer::n_sp_pol=1
 real(SP)::spin_occ=2._SP
 type levels
   integer::nb=2,nk=2,nbf(2)=1,nbm(2)=1
   real(SP)::E_Fermi=0._SP
   real(SP),allocatable::E(:,:,:),f(:,:,:),Eo(:,:,:)
 end type
contains
 integer function spin(table)
   integer::table(:)
   spin=1
 end function
end module
module R_lattice
 use pars
 type bz_samp
   integer::nibz=2,nbz=2
   real(SP),allocatable::pt(:,:)
   integer,allocatable::sstar(:,:)
 end type
 integer,allocatable::qindx_X(:,:,:)
 integer,allocatable::qindx_S(:,:,:)
 complex(SP),allocatable::bare_qpg(:,:)
 real(SP)::q0_def_norm=1.E-5_SP
 integer::g_rot(2,1)=reshape([1,2],[2,1])
end module
module D_lattice
 use pars
 real(SP)::DL_vol=1._SP
 integer::nsym=1,i_time_rev=0,sop_inv(1)=[1]
end module
module vec_operate
 use pars
contains
 real(SP) function v_norm(v)
   real(SP),intent(in)::v(:)
   v_norm=sqrt(sum(v*v))
 end function
end module
module frequency
 use pars
 type w_samp
   integer::n_freqs=2
   real(SP)::er(2)=0._SP,dr(2)=0._SP
   character(16)::grid_type='ra'
   complex(SP),allocatable::p(:)
 end type
contains
 subroutine W_reset(w)
   type(w_samp)::w
   if(allocated(w%p))deallocate(w%p)
 end subroutine
 subroutine W_duplicate(a,b)
   type(w_samp)::a,b
   b=a
 end subroutine
end module
module matrix
 use pars
 type PAR_matrix
   integer::rows(2)=[1,2],cols(2)=[1,2]
   complex(SP),allocatable::blc(:,:,:)
 end type
end module
module X_m
 use pars
 use matrix
 type X_t
   integer::ng=2,ib(2)=[1,2],whoami=2
   real(SP)::q0(3)=[1._SP,0._SP,0._SP]
 end type
 type(PAR_matrix),allocatable::X_par(:)
 type(PAR_matrix)::X_par_lower_triangle
 logical::X_lower_triangle_matrix_in_use=.FALSE.
 integer::current_iq=0
end module
module DIPOLES
 use pars
 complex(SP),allocatable::DIP_iR(:,:,:,:,:)
 type DIPOLE_t
   real(SP)::q0(3)=[1._SP,0._SP,0._SP]
 end type
end module
module ALLOC
contains
 subroutine DIPOLE_ALLOC_global()
 end subroutine
end module
#ifndef _TEST_REAL_QP_MODULE
module QP_m
 use pars
 type QP_t
   integer::nk=2,nb=2,n_states=0,GreenF_n_steps=0
   logical::GreenF_retarded=.FALSE.
   real(SP)::GreenF_mu=0._SP
   integer,allocatable::table(:,:)
   real(SP),allocatable::k(:,:),E_bare(:)
   complex(SP),allocatable::GreenF(:,:),GreenF_W(:,:),S_total(:,:)
 end type
 real(SP)::QP_G_er(2)=0._SP,QP_G_dr(2)=0._SP,QP_G_Zoom_treshold=0._SP
 logical::GF_is_causal=.FALSE.
contains
 subroutine QP_reset(qp)
   type(QP_t)::qp
   if(allocated(qp%table))deallocate(qp%table)
   if(allocated(qp%k))deallocate(qp%k)
   if(allocated(qp%E_bare))deallocate(qp%E_bare)
   if(allocated(qp%GreenF))deallocate(qp%GreenF)
   if(allocated(qp%GreenF_W))deallocate(qp%GreenF_W)
   if(allocated(qp%S_total))deallocate(qp%S_total)
   qp%GreenF_retarded=.FALSE.
 end subroutine
end module
#endif
module drivers
 logical::Finite_Tel=.FALSE.
end module
module functions
 use pars
 use, intrinsic::ieee_arithmetic
contains
 logical function NAN(x)
   real(SP)::x
   NAN=.not.ieee_is_finite(x)
 end function
 real(SP) function bose_f(x)
   real(SP)::x
   bose_f=0._SP
 end function
end module
module QP_CTL_m
 use pars
 type ctl
   character(schlen)::action='none'
 end type
 type(ctl)::QP_ctl_DB_user(1)
end module
module global_XC
 integer,parameter::QP_SE_COHSEX=7
 integer::QP_DB_kind=0
end module
module IO_m
 integer,parameter::OP_RD_CL=1,OP_WR_CL=2,OP_APP_CL=3,DUMP=1,REP=1
 logical::io_DIP=.TRUE.
end module
module IO_int
contains
 subroutine io_control(ACTION,COM,SEC,MODE,ID)
   integer::ACTION,COM,SEC(:),ID
   integer,optional::MODE
   ID=1
 end subroutine
end module
module com
 use pars
 interface msg
   module procedure msg_scalar,msg_vector,msg_integer
 end interface
contains
 subroutine msg_scalar(a,b,c)
   character(*)::a,b
   real(SP)::c
 end subroutine
 subroutine msg_vector(a,b,c,d)
   character(*)::a,b
   character(*),optional::d
   real(SP)::c(:)
 end subroutine
 subroutine msg_integer(a,b,c)
   character(*)::a,b
   integer::c
 end subroutine
end module
module stderr
contains
 function intc(i)result(s)
   integer::i
   character(24)::s
   write(s,'(i0)')i
 end function
end module
module LIVE_t
contains
 subroutine live_timing()
 end subroutine
end module
module timing_m
contains
 subroutine timing(s,OPR)
   character(*)::s,OPR
 end subroutine
end module
module parser_m
contains
 subroutine parser(s,l)
   character(*)::s
   logical::l
   l=.FALSE.
 end subroutine
end module
module parallel_m
 integer::ncpu=1,myid=0
 logical::master_cpu=.TRUE.
 type comm_t
   integer::COMM=0,CPU_id=0,n_CPU=1
 end type
 type index_t
   logical,allocatable::element_1D(:)
 end type
 type(comm_t)::PAR_COM_X_WORLD,PAR_COM_X_WORLD_RL_resolved
 type(comm_t)::PAR_COM_Q_INDEX,PAR_COM_RL_INDEX,PAR_COM_CON_INDEX_X(5),PAR_COM_VAL_INDEX_X(5)
 type(index_t)::PAR_IND_Xk_bz
end module
module parallel_int
 use pars
#ifdef _TEST_MPI
 use mpi
#endif
 interface PP_redux_wait
   module procedure reduce_sp
#ifdef _TEST_SINGLE
   module procedure reduce_dp
#endif
 end interface
contains
 subroutine PP_wait(COMM)
   integer::COMM
#ifdef _TEST_MPI
   integer::ierr
   call MPI_BARRIER(COMM,ierr)
   if(ierr/=MPI_SUCCESS)call MPI_ABORT(COMM,1,ierr)
#endif
 end subroutine
 subroutine reduce_sp(a,IN_PLACE,COMM)
   complex(SP)::a(:,:)
   logical,optional::IN_PLACE
   integer::COMM
#ifdef _TEST_MPI
   integer::ierr,datatype
   datatype=MPI_DOUBLE_COMPLEX
#ifdef _TEST_SINGLE
   datatype=MPI_COMPLEX
#endif
   call MPI_ALLREDUCE(MPI_IN_PLACE,a,size(a),datatype,MPI_SUM,COMM,ierr)
   if(ierr/=MPI_SUCCESS)call MPI_ABORT(COMM,1,ierr)
#endif
 end subroutine
#ifdef _TEST_SINGLE
 subroutine reduce_dp(a,COMM)
   complex(DP)::a(:,:)
   integer::COMM
#ifdef _TEST_MPI
   integer::ierr
   call MPI_ALLREDUCE(MPI_IN_PLACE,a,size(a),MPI_DOUBLE_COMPLEX,MPI_SUM,COMM,ierr)
   if(ierr/=MPI_SUCCESS)call MPI_ABORT(COMM,1,ierr)
#endif
 end subroutine
#endif
end module
module collision_el
 use pars
 type elemental_collision
   integer::is(4),os(4),qs(3)
   complex(SP),allocatable::rhotw(:)
 end type
contains
 subroutine elemental_collision_alloc(s,NG,CPU_ONLY,TITLE)
   type(elemental_collision)::s
   integer::NG
   logical::CPU_ONLY
   character(*)::TITLE
   allocate(s%rhotw(NG))
 end subroutine
 subroutine elemental_collision_free(s)
   type(elemental_collision)::s
   if(allocated(s%rhotw))deallocate(s%rhotw)
 end subroutine
end module
module pipeline_output
 use pars
 complex(SP),allocatable::saved_c0(:,:,:),saved_p(:,:,:),saved_fxc(:,:,:),saved_freq(:)
 real(SP),allocatable::saved_rcond(:,:)
end module
subroutine scatter_Bamp(s)
 use pars
 use collision_el
 type(elemental_collision)::s
 s%rhotw=0._SP
 if(s%is(1)==s%os(1))return
 s%rhotw=[cmplx(1._SP,0._SP,SP),cmplx(.2_SP,0._SP,SP)]
 if(s%is(2)==2)s%rhotw=s%rhotw(2:1:-1)
end subroutine
subroutine DIPOLE_rotate(a,b,ikbz,i_spin,what,k,dipole)
 use pars
 use R_lattice,only:bz_samp
 integer::a,b,ikbz,i_spin
 character(*)::what
 type(bz_samp)::k
 complex(SP)::dipole(3)
 if(a<=b)error stop 'DIP_iR mock stores only conduction/valence dipoles'
 dipole=0._SP
 dipole(1)=cmplx(0._SP,1._SP,SP)
end subroutine
subroutine DIPOLE_dimensions(e,d,bands,q0)
 use pars
 use electrons,only:levels
 use DIPOLES,only:DIPOLE_t
 type(levels)::e
 type(DIPOLE_t)::d
 integer::bands(2)
 real(SP)::q0(3)
 d%q0=q0
end subroutine
subroutine DIPOLE_IO(k,e,d,action,io_err,scheme)
 use R_lattice,only:bz_samp
 use electrons,only:levels
 use DIPOLES,only:DIPOLE_t
 type(bz_samp)::k
 type(levels)::e
 type(DIPOLE_t)::d
 character(*)::action,scheme
 integer::io_err
 io_err=0
end subroutine
subroutine FREQUENCIES_reset(w,s)
 use frequency
 type(w_samp)::w
 character(*)::s
end subroutine
subroutine X_ALLOC_parallel(m,ng,nw,s)
 use matrix
 type(PAR_matrix)::m
 integer::ng,nw
 character(*)::s
 if(allocated(m%blc))deallocate(m%blc)
 allocate(m%blc(ng,ng,nw))
 m%blc=0._SP
end subroutine
subroutine X_irredux(iq,s,m,e,k,w,x,d)
 use pars
 use matrix
 use electrons
 use R_lattice
 use frequency
 use X_m,only:X_t
 use DIPOLES
 integer::iq,iw
 character(*)::s
 type(PAR_matrix)::m
 type(levels)::e
 type(bz_samp)::k
 type(w_samp)::w
 type(X_t)::x
 type(DIPOLE_t)::d
 complex(SP)::factor
 do iw=1,w%n_freqs
   factor=1._SP/(w%p(iw)-2._SP)-1._SP/(w%p(iw)+2._SP)
   m%blc(:,:,iw)=factor*reshape([1.04_SP,.4_SP,.4_SP,1.04_SP],[2,2])
 enddo
end subroutine
subroutine error(s)
#ifdef _TEST_MPI
 use mpi
 integer::ierr
#endif
 character(*)::s
 print *,trim(s)
#ifdef _TEST_MPI
 call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
 error stop 1
end subroutine
