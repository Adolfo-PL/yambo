!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
program harness
 ! QP_ppa_cohsex (binned retarded poles) against QP_ppa_cohsex_direct (frequency loop)
 use pars; use electrons; use R_lattice; use X_m; use frequency; use QP_m; use parallel_m
 use drivers, only:l_rim_w; use harness_data
 implicit none
 type(levels)::E; type(bz_samp)::k,q; type(X_t)::X; type(QP_t)::qp; type(w_samp)::Xw
 complex(SP),allocatable::blc0(:,:,:,:),S_bin(:,:),S_dir(:,:)
 real(SP)::ha,u(4); integer::ib,i,j,iq,ic,seed(64); real(DP)::err
 character(40)::case_name(3)=['RIM-W + BG terminator','RIM-W, no terminator ','no RIM-W, BG         ']
 ha=27.211386245988_SP
 seed=7; call random_seed(put=seed)
 E%nb=20; allocate(E%E(20,1,1),E%f(20,1,1))
 do ib=1,20
   E%E(ib,1,1)=(-12._SP+real(ib,SP))/ha
   E%f(ib,1,1)=0._SP; if (ib<=13) E%f(ib,1,1)=spin_occ
 enddo
 k%nbz=1; k%nibz=1; allocate(k%sstar(1,2)); k%sstar=1
 q%nbz=4; q%nibz=2; allocate(q%sstar(4,2)); q%sstar(:,1)=[1,1,2,2]; q%sstar(:,2)=[1,2,1,2]
 allocate(qindx_S(1,4,2)); qindx_S=1
 QP_n_states=3; qp%n_states=3; allocate(QP_table(3,3),qp%E_bare(3),Eo(3))
 do i=1,3; QP_table(i,:)=[12+i,12+i,1]; Eo(i)=E%E(12+i,1,1); qp%E_bare(i)=Eo(i); enddo
 QP_n_G_bands=[1,20]
 allocate(PAR_IND_G_b%element_1D(20)); PAR_IND_G_b%element_1D=.true.
 allocate(PAR_IND_QP%element_1D(3),PAR_IND_QP%n_of_elements(1)); PAR_IND_QP%element_1D=.true.; PAR_IND_QP%n_of_elements=3
 allocate(PAR_IND_Q_ibz%element_1D(2),PAR_IND_Q_ibz%n_of_elements(1)); PAR_IND_Q_ibz%element_1D=.true.; PAR_IND_Q_ibz%n_of_elements=2
 allocate(PAR_Q_ibz_index(2)); PAR_Q_ibz_index=[1,2]
 X%ng=30; X%ppaE=1._SP; Xw%n_freqs=2; nqibz=2
 allocate(X_par(2),blc0(30,30,2,2))
 do iq=1,2
   allocate(X_par(iq)%blc(30,30,2)); X_par(iq)%rows=[1,30]; X_par(iq)%cols=[1,30]
   do j=1,30; do i=1,30
     call random_number(u)
     blc0(i,j,1,iq)=cmplx(-(0.1_SP+u(1)),0.05_SP*(u(2)-0.5_SP),SP)
     blc0(i,j,2,iq)=blc0(i,j,1,iq)*(0.02_SP+0.98_SP*u(3))        ! PPA fit: E = ppaE/sqrt(X0/X1-1)
     if (u(4)>0.95_SP) blc0(i,j,2,iq)=blc0(i,j,1,iq)*1.2_SP       ! some failing fits
   enddo; enddo
 enddo
 allocate(G_m_G(30,30)); G_m_G_maxval=60
 do j=1,30; do i=1,30; call random_number(u); G_m_G(i,j)=1+int(u(1)*59.999_SP); enddo; enddo
 RIM_W_ng=5; allocate(RIM_W(1,2,5,5),RIM_W_E(5,5))
 do j=1,5; do i=1,5; call random_number(u); RIM_W(1,:,i,j)=-(0.2_SP+u(1)); RIM_W_E(i,j)=0.3_SP+2._SP*u(2); enddo; enddo
 QP_Sc_steps=589; QP_dSc_steps=589
 allocate(QP_Sc(3,589),S_bin(3,589),S_dir(3,589))
 do ic=1,3
   l_rim_w=(ic/=3); l_GW_terminator=(ic/=2)
   call set_X(); call QP_ppa_cohsex(X,k,E,k,q,qp,Xw,0); S_bin=QP_Sc
   call set_X(); call QP_ppa_cohsex_direct(X,k,E,k,q,qp,Xw,0); S_dir=QP_Sc
   err=maxval(abs(S_bin-S_dir))/maxval(abs(S_dir))
   write(*,'(a,a,es9.2,a,3f10.5,a,f9.5)') case_name(ic),': max|bin-direct|/max|Sc| ',err,&
&    '  Re Sc(Eo) [eV]',real(S_dir(:,295))*ha,'  max|Sc| [eV]',maxval(abs(S_dir))*ha
   if (.not.(err<1.E-4_DP)) then; print *,'FAIL'; error stop 1; endif
 enddo
 print *,'PASS: binned retarded PPA self-energy = frequency loop'
contains
 subroutine set_X()
   ! QP_ppa_cohsex fits the PPA in place and frees X_par at the end
   integer :: jq
   if (allocated(X_par)) deallocate(X_par)
   allocate(X_par(2))
   do jq=1,2
     allocate(X_par(jq)%blc(30,30,2)); X_par(jq)%blc=blc0(:,:,:,jq)
     X_par(jq)%rows=[1,30]; X_par(jq)%cols=[1,30]
   enddo
 end subroutine
end program
