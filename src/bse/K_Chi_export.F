!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Reading copy: the numerics only (BSEChiOut = "exc" or "full", any q). No guards, reports or timing.
!
! LEGEND (names as in the full code)
!   n, nh    N resonant transitions, 2N rows of BS_mat (resonant + anti-resonant)
!   ng, ns   plane waves of the kernel; of those in the kernel set S (all G, or G_z = 0 with BSEGinplane)
!   H0       BS_mat of this BSE run (energies, exchange, W), double precision
!   Vf       (ng,nh) vertices rho_A(G): Chi_dyn_V resonant | de-excitation, collected in K_kernel
!   sf       (nh)    s_A = sqrt(f_A); i sqrt|f_A| for f_A < 0
!   E_ref    (nh)    energies of the reference bubble: BSE energies ("exc") or KS ones ("full")
!   d        (ng)    D_G = sqrt(4 pi)/|q+G|, so that v_G = D_G^2
!   zz       (nw)    z = 0, then the BSE frequencies w + i eta
!   chib     BSE response chib(G,G') (exchange included)
!   chi0     reference bubble
!   chs,c0s  chib, chi0 restricted to S;   fs,Ps  kernel and polarization on S;   f,Pm  the same on all G
!   l_full_L Lkind = "full": G = 0 exchange kept in vbar
!
subroutine K_Chi_export(iq,Ken,Xk,W)
 !
 ! fxc(z) = chi0(z)^-1 - ( chib(z)^-1 + vbar ),   chib = c U (z - H0)^-1 U^H,   chi0 = c U (z - E_ref)^-1 U^H
 !
 use pars,           ONLY:SP,DP,pi
 use frequency,      ONLY:w_samp,W_reset,W_duplicate
 use electrons,      ONLY:levels
 use R_lattice,      ONLY:bz_samp,bare_qpg,qindx_X,g_vec
 use stderr,         ONLY:STRING_match
 use X_m,            ONLY:X_t,X_reset
 use BS,             ONLY:BS_K_dim,BS_H_dim,BS_n_g_fxc,BS_bands,BS_K_Chi_out,BSE_L_kind,l_BS_ares_from_res,BS_G_inplane
 use BS_solvers,     ONLY:BS_mat,BSS_eh_E,BSS_eh_table,BSS_eh_f
 use Chi_m,          ONLY:Chi_ng,Chi_n_freqs,Chi_alloc,Chi_free,Chi_freqs,Chi_0,Chi_P,Chi_fxc,Chi_dyn_V,Chi_dyn_c,Chi_dyn_free
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:OP_WR_CL,OP_APP_CL,REP
 use TDDFT_Chi_dynamic_m, ONLY:TDDFT_Chi_full_response,TDDFT_Chi_full_bubble,TDDFT_Chi_kernel_from_response
 !
 implicit none
 !
 integer       :: iq
 type(levels)  :: Ken
 type(bz_samp) :: Xk
 type(w_samp)  :: W
 !
 logical                  :: l_exc,l_qp,l_full_L
 logical, save            :: header_done=.FALSE.
 integer                  :: n,nh,ng,ns,nw,iw,it,ig,ik,ikv,ID,io_err
 integer,     external    :: io_Chi
 integer,     allocatable :: S(:)
 real(DP),    allocatable :: E_ref(:)
 complex(DP), allocatable :: d(:),H0(:,:),zz(:),Vf(:,:),sf(:)
 complex(DP), allocatable :: chib(:,:),chi0(:,:),f(:,:),Pm(:,:),chs(:,:),c0s(:,:),fs(:,:),Ps(:,:)
 type(X_t)                :: X_io
 type(w_samp)             :: W_io
 !
 l_exc=trim(BS_K_Chi_out)=='exc'
 l_full_L=STRING_match(BSE_L_kind,"full")
 !
 call K_components_folded_in_serial_arrays(iq)
 call K_stored_in_a_big_matrix(1,iq,"DIAGO")
 !
 n =BS_K_dim(1)
 nh=BS_H_dim
 ng=BS_n_g_fxc
 !
 allocate(E_ref(nh),d(ng),Vf(ng,nh),sf(nh))
 Vf(:,:n)=Chi_dyn_V(:,:,1)
 Vf(:,n+1:)=Chi_dyn_V(:,:,2)
 do it=1,nh
   sf(it)=sqrt(cmplx(BSS_eh_f(it),kind=DP))
 enddo
 l_qp=.FALSE.
 if (allocated(Ken%Eo)) l_qp=maxval(abs(Ken%Eo-Ken%E))>1.E-6_SP
 do it=1,nh
   E_ref(it)=real(BSS_eh_E(it),DP)
   if (.not.l_exc.and.l_qp) then
     if (l_BS_ares_from_res.and.it>n) then
       E_ref(it)=-E_ref(it-n)
       cycle
     endif
     ik =Xk%sstar(BSS_eh_table(it,1),1)
     ikv=Xk%sstar(qindx_X(iq,BSS_eh_table(it,1),1),1)
     E_ref(it)=real(Ken%Eo(BSS_eh_table(it,3),ik ,BSS_eh_table(it,4))-&
&                   Ken%Eo(BSS_eh_table(it,2),ikv,BSS_eh_table(it,5)),DP)
   endif
 enddo
 do ig=1,ng
   d(ig)=sqrt(4._DP*pi)/cmplx(bare_qpg(iq,ig),kind=DP)
 enddo
 !
 allocate(S(ng))
 ns=0
 do ig=1,ng
   if (BS_G_inplane.and.abs(g_vec(ig,3))>1.E-5_SP) cycle
   ns=ns+1
   S(ns)=ig
 enddo
 !
 nw=W%n_freqs+1
 allocate(zz(nw))
 zz(1)=(0._DP,0._DP)
 zz(2:)=cmplx(W%p(:),kind=DP)
 !
 call Chi_alloc(ng,nw)
 Chi_freqs=cmplx(zz,kind=SP)
 allocate(H0(nh,nh))
 H0=cmplx(BS_mat,kind=DP)
 !
 !$omp parallel default(shared) private(iw,chib,chi0,f,Pm,chs,c0s,fs,Ps)
 allocate(chib(ng,ng),chi0(ng,ng),f(ng,ng),Pm(ng,ng))
 allocate(chs(ns,ns),c0s(ns,ns),fs(ns,ns),Ps(ns,ns))
 !$omp do schedule(dynamic)
 do iw=1,nw
   call TDDFT_Chi_full_response(zz(iw),H0,Vf,sf,Chi_dyn_c,chib)
   call TDDFT_Chi_full_bubble(zz(iw),E_ref,Vf,sf,Chi_dyn_c,chi0)
   chs=chib(S(:ns),S(:ns))
   c0s=chi0(S(:ns),S(:ns))
   call TDDFT_Chi_kernel_from_response(c0s,chs,d(S(:ns)),l_full_L,fs,Ps)
   f=(0._DP,0._DP)
   Pm=(0._DP,0._DP)
   f(S(:ns),S(:ns))=fs
   Pm(S(:ns),S(:ns))=Ps
   Chi_0(:,:,iw)=cmplx(chi0,kind=SP)
   Chi_P(:,:,iw)=cmplx(Pm,kind=SP)
   Chi_fxc(:,:,iw)=cmplx(f,kind=SP)
 enddo
 !$omp end do
 deallocate(chib,chi0,f,Pm,chs,c0s,fs,Ps)
 !$omp end parallel
 !
 Chi_ng=ng
 Chi_n_freqs=nw
 call X_reset(X_io)
 X_io%ib=BS_bands
 X_io%ng=ng
 call W_reset(W_io)
 call W_duplicate(W,W_io)
 if (.not.header_done) then
   call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),ID=ID)
   io_err=io_Chi(X_io,W_io,ID)
   header_done=.TRUE.
 endif
 call io_control(ACTION=OP_APP_CL,COM=REP,SEC=(/2*iq,2*iq+1/),ID=ID)
 io_err=io_Chi(X_io,W_io,ID)
 !
 call W_reset(W_io)
 call Chi_free()
 call Chi_dyn_free()
 deallocate(E_ref,d,Vf,sf,zz,H0,S)
 !
end subroutine K_Chi_export
