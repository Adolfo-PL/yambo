!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Reading copy: the numerics only (BSEChiDyn, coupling, as run). No checks, text dump or timing.
!
! LEGEND (names as in the full code)
!   n, nh    N resonant transitions, 2N rows of BS_mat
!   H0       BS_mat of this CHI run: diag(E) + exchange, fxc left out by K_kernel
!   Vf, sf   (Chi_ng,nh), (nh)  vertices rho_A(G) and s_A = sqrt(f_A) of every row, as in K_Chi_export
!   df       (nh)    optical dipoles d_A (BSS_dipoles_opt)
!   Chi_fxc  (Chi_ng,Chi_ng,nfreq)  kernel read from ndb.Chi; i0 = 1 skips its z = 0 sample
!   W%p      BSE frequencies w + i eta (must be the ndb.Chi grid)
!   F        fxc(z) at one frequency;   resp_w  head response at every frequency
!
subroutine K_Chi_dynamic(iq,Xk,W)
 !
 ! Resp(w) = sum_AB s_A conj(d_A) [ z - H0 - Kf(z) ]^-1_AB s_B d_B,   Resp_ii = -Co_factor Resp
 !
 use pars,           ONLY:SP,DP
 use frequency,      ONLY:w_samp
 use R_lattice,      ONLY:bz_samp
 use BS,             ONLY:BS_K_dim,BS_H_dim
 use BS_solvers,     ONLY:BS_mat,BSS_dipoles_opt,BSS_eh_f,Co_factor,para_term_w0
 use X_m,            ONLY:Resp_ii,X_OUTPUT_driver
 use Chi_m,          ONLY:Chi_freqs,Chi_fxc,Chi_n_freqs,Chi_ng,Chi_free,Chi_dyn_V,Chi_dyn_c,Chi_dyn_free
 use TDDFT_Chi_dynamic_m, ONLY:TDDFT_Chi_full_dynamic_response
 !
 implicit none
 !
 integer       :: iq
 type(bz_samp) :: Xk
 type(w_samp)  :: W
 !
 integer                  :: n,nh,nw,iw,it,i0
 complex(DP)              :: resp
 complex(DP), allocatable :: H0(:,:),F(:,:),resp_w(:),Vf(:,:),sf(:),df(:)
 !
 call X_OUTPUT_driver("DEF INIT",solver="ChiDyn")
 !
 nw=W%n_freqs
 i0=0
 if (Chi_n_freqs==nw+1) then
   if (abs(Chi_freqs(1))==0._SP) i0=1
 endif
 !
 call K_components_folded_in_serial_arrays(iq)
 call K_stored_in_a_big_matrix(1,iq,"DIAGO")
 !
 n =BS_K_dim(1)
 nh=BS_H_dim
 allocate(H0(nh,nh),resp_w(nw))
 H0=cmplx(BS_mat,kind=DP)
 !
 allocate(Vf(Chi_ng,nh),sf(nh),df(nh))
 Vf(:,:n)=Chi_dyn_V(:,:,1)
 Vf(:,n+1:)=Chi_dyn_V(:,:,2)
 do it=1,nh
   sf(it)=sqrt(cmplx(BSS_eh_f(it),kind=DP))
   df(it)=cmplx(BSS_dipoles_opt(1,it),kind=DP)
 enddo
 !
 !$omp parallel default(shared) private(iw,F,resp)
 allocate(F(Chi_ng,Chi_ng))
 !$omp do schedule(dynamic)
 do iw=1,nw
   F=cmplx(Chi_fxc(:,:,iw+i0),kind=DP)
   call TDDFT_Chi_full_dynamic_response(cmplx(W%p(iw),kind=DP),H0,F,Vf,sf,Chi_dyn_c,df,resp)
   resp_w(iw)=resp
 enddo
 !$omp end do
 deallocate(F)
 !$omp end parallel
 !
 Resp_ii(:,2)=cmplx(-real(Co_factor,DP)*resp_w,kind=SP)
 para_term_w0(2)=para_term_w0(1)
 call K_OUTPUT(iq,W,"open dump close ChiDyn","ALL")
 !
 deallocate(H0,resp_w,Vf,sf,df)
 call Chi_dyn_free()
 call Chi_free()
 !
end subroutine K_Chi_dynamic
