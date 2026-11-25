!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Stub external routines for the QP_ppa_cohsex harness
module harness_data
 use pars
 real(SP),allocatable::Eo(:)
end module
subroutine error(m); character(*)::m; print *,'ERROR: ',m; error stop 1; end subroutine
subroutine section(m,n); character(*)::m,n; end subroutine
subroutine k_expand(k); use R_lattice; type(bz_samp)::k; end subroutine
subroutine WF_spatial_inversion(E,k); use electrons; use R_lattice; type(levels)::E; type(bz_samp)::k; end subroutine
subroutine X_ALLOC_parallel(M,ng,nw,w); use matrix; type(PAR_matrix)::M; integer::ng,nw; character(*)::w
 stop 'X_ALLOC_parallel not expected'
end subroutine
integer function io_X(X,Xw,ID); use X_m; use frequency; type(X_t)::X; type(w_samp)::Xw; integer::ID; io_X=0; end function
subroutine FREQUENCIES_Green_Function(i_qp,W,E,l)
 use pars; use frequency; use harness_data; use QP_m, only:QP_Sc_steps
 integer::i_qp,i; type(w_samp)::W; real(SP)::E(:,:,:); logical::l; real(SP)::ha
 ha=27.211386245988_SP
 W%n_freqs=QP_Sc_steps; allocate(W%p(W%n_freqs))
 do i=1,W%n_freqs
   W%p(i)=cmplx(Eo(i_qp)+(-40._SP+80._SP*real(i-1,SP)/real(W%n_freqs-1,SP))/ha,0.136_SP/ha,SP)
 enddo
end subroutine
subroutine scatter_Gamp_gpu(c,mode)
 use collision_el; type(elemental_collision)::c; character(*)::mode; integer::i,j
 do j=lbound(c%gamp,2),ubound(c%gamp,2); do i=lbound(c%gamp,1),ubound(c%gamp,1)
   c%gamp(i,j)=cmplx((1._SP+0.1_SP*c%qs(2))/(1._SP+0.1_SP*real(i+j,SP)),0._SP,SP)
 enddo; enddo
end subroutine
subroutine scatter_Bamp_gpu(c)
 use collision_el; type(elemental_collision)::c; integer::g; real(SP)::a
 do g=1,size(c%rhotw)
   a=1.3_SP*g+0.7_SP*c%os(1)+0.31_SP*c%is(1)+0.17_SP*c%qs(1)+0.05_SP*c%os(2)+0.011_SP*c%qs(2)
   c%rhotw(g)=cmplx(sin(a),cos(1.7_SP*a),SP)/(1._SP+0.05_SP*g)
 enddo
end subroutine
complex(kind(1.0)) function QP_ppa_EET_terminator(W,E,vs,Ep,ig1,ig2,iq,cutoff0,cutoff1)
 ! BG branch of src/qp/QP_ppa_EET_terminator.F
 use pars; use electrons; use drivers, only:Finite_Tel; use QP_m, only:GW_terminator_E
 complex(SP)::W,lW; type(levels)::E; integer::vs(4),ig1,ig2,iq; real(SP)::Ep,cutoff0,cutoff1,eet_energy
 eet_energy=cutoff1+GW_terminator_E
 lW=W; if (Finite_Tel) lW=conjg(W)
 QP_ppa_EET_terminator=real(spin_occ)/(lW-eet_energy-Ep)
end function
