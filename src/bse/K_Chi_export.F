!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
subroutine K_Chi_export(iq,Ken,Xk,W)
 !
 ! BSEChiOut: the xc kernel that reproduces this BSE, written as an ndb.Chi.
 !
 ! The coupled BSE matrix M (its transition energies, the exchange vbar and, with
 ! BSKmod= "SEX", the screened e-h attraction W) gives the G-resolved response
 !
 !   chib(G,G';z) = c B_G^T [z - M]^-1 A_G'              (TDDFT_Chi_pair_response)
 !
 ! on the Chi vertices that K_kernel collected. With P^-1 = chib^-1 + vbar the kernel
 !
 !   "exc" : fxc = Pqp^-1 - P^-1,  Pqp = bubble of the BSE transition energies
 !   "full": fxc = chi0^-1 - P^-1, chi0 = bubble of the KS transition energies
 !
 ! makes the Casida response with the same reference energies and vbar equal to chib,
 ! and so reproduces this BSE, exciton included (tests/chi_dynamic). ndb.Chi stores
 ! chi0/Pqp as "chi0", P and fxc. Its first frequency is z = 0 exactly (static Casida),
 ! followed by the BSE frequencies (BSEChiDyn with the same BSE grid).
 !
 use pars,           ONLY:SP,DP,pi,schlen,cZERO
 use units,          ONLY:HA2EV
 use frequency,      ONLY:w_samp,W_reset,W_duplicate
 use electrons,      ONLY:levels
 use R_lattice,      ONLY:bz_samp,bare_qpg
 use com,            ONLY:msg,isec
 use stderr,         ONLY:intc,real2ch
 use parallel_m,     ONLY:ncpu
 use X_m,            ONLY:X_t,X_reset,global_gauge
 use BS,             ONLY:BS_K_coupling,BS_K_dim,BS_H_dim,BS_res_ares_n_mat,BS_n_g_fxc,BS_bands,BS_K_Chi_out
 use BS_solvers,     ONLY:BS_mat,BSS_eh_E,BSS_eh_W,BSS_eh_Z,BSS_eh_table,BSS_q0
 use Chi_m,          ONLY:Chi_ng,Chi_n_freqs,Chi_alloc,Chi_free,Chi_freqs,Chi_0,Chi_P,Chi_fxc,Chi_rcond,&
&                         Chi_QP_applied,Chi_QP_action,Chi_convention,Chi_G_mode,Chi_G_axis,Chi_G_db,&
&                         Chi_optical_q0,Chi_optical_q0_known,Chi_G0_error,Chi_rcond_floor,&
&                         Chi_dyn_V,Chi_dyn_s,Chi_dyn_filled,Chi_dyn_c,Chi_dyn_free
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:OP_WR_CL,OP_APP_CL,REP
 use TDDFT_Chi_dynamic_m, ONLY:TDDFT_Chi_pair_response,TDDFT_Chi_pair_bubble,TDDFT_Chi_kernel_from_response
 use timing_m,       ONLY:timing
 !
 implicit none
 !
 integer       :: iq
 type(levels)  :: Ken
 type(bz_samp) :: Xk
 type(w_samp)  :: W
 !
 ! Work Space
 !
 character(1)             :: sec
 logical                  :: l_exc,l_qp
 integer                  :: n,nh,ng,nw,iw,it,ig,jg,ik,ierr,n_bad,ID,io_err
 integer,     external    :: io_Chi
 integer,     allocatable :: err_w(:)
 real(DP)                 :: herm,scale_max,rc(2)
 real(DP),    allocatable :: E_ref(:)
 complex(DP), allocatable :: d(:),H0(:,:),zz(:),chib(:,:),chi0(:,:),f(:,:),Pm(:,:)
 type(X_t)                :: X_io
 type(w_samp)             :: W_io
 !
 call timing('ChiExport',OPR='start')
 !
 if (isec(2)/=0) sec="="
 if (isec(2)==0) sec="+"
 call section(sec,'Chi kernel from the BSE @q'//trim(intc(iq)))
 !
 l_exc=trim(BS_K_Chi_out)=='exc'
 !
 if (iq/=1) call error('[Chi/BSK] BSEChiOut: optical q only (BSEQptR 1)')
 if (ncpu/=1) call error('[Chi/BSK] BSEChiOut is serial (MPI); use OpenMP threads')
 if (.not.BS_K_coupling) call error('[Chi/BSK] BSEChiOut needs BSEmod= "coupling"')
 if (BS_res_ares_n_mat/=1) call error('[Chi/BSK] BSEChiOut: one BSE matrix only')
 if (trim(global_gauge)/='length') call error('[Chi/BSK] BSEChiOut: length gauge only')
 if (.not.allocated(Chi_dyn_V)) call error('[Chi/BSK] BSEChiOut: no Chi vertices (the kernel must be built in this run)')
 if (.not.all(Chi_dyn_filled)) call error('[Chi/BSK] BSEChiOut: some transitions have no Chi vertex')
 !
 call K_components_folded_in_serial_arrays(iq)
 if (allocated(BSS_eh_W).or.allocated(BSS_eh_Z)) &
&  call error('[Chi/BSK] BSEChiOut: transition widths/Z factors are not supported')
 call K_stored_in_a_big_matrix(1,iq,"DIAGO")
 !
 n =BS_K_dim(1)
 nh=BS_H_dim
 ng=BS_n_g_fxc
 if (size(BS_mat,1)/=nh.or.nh/=2*n.or.size(Chi_dyn_V,1)/=ng) call error('[Chi/BSK] BSEChiOut: inconsistent dimensions')
 !
 ! Reference transition energies of the bubble the kernel is measured from
 !
 allocate(E_ref(n),d(ng))
 l_qp=.FALSE.
 if (allocated(Ken%Eo)) l_qp=maxval(abs(Ken%Eo-Ken%E))>1.E-6_SP
 do it=1,n
   E_ref(it)=real(BSS_eh_E(it),DP)
   if (.not.l_exc.and.l_qp) then
     ik=Xk%sstar(BSS_eh_table(it,1),1)
     E_ref(it)=real(Ken%Eo(BSS_eh_table(it,3),ik,BSS_eh_table(it,4))-&
&                   Ken%Eo(BSS_eh_table(it,2),ik,BSS_eh_table(it,5)),DP)
   endif
 enddo
 do ig=1,ng
   ! as the BSE exchange (1/bare_qpg**2): complex and signed with a Coulomb cutoff
   d(ig)=sqrt(4._DP*pi)/cmplx(bare_qpg(iq,ig),kind=DP)
 enddo
 !
 call msg('r','[Chi/BSK] BSEChiOut kernel reference',trim(BS_K_Chi_out))
 call msg('r','[Chi/BSK] QP-corrected BSE transitions',l_qp)
 call msg('r','[Chi/BSK] lowest BSE transition energy (gap of the window) [eV]',real(minval(BSS_eh_E(1:n))*HA2EV,SP))
 if (l_exc.and..not.l_qp) call warning('[Chi/BSK] BSEChiOut= "exc" on KS transitions: the reference bubble is the KS one')
 !
 ! Frequencies: z=0 first, then the BSE grid
 !
 nw=W%n_freqs+1
 allocate(zz(nw),err_w(nw))
 zz(1)=(0._DP,0._DP)
 zz(2:)=cmplx(W%p(:),kind=DP)
 !
 call Chi_alloc(ng,nw)
 Chi_freqs=cmplx(zz,kind=SP)
 allocate(H0(nh,nh))
 H0=cmplx(BS_mat,kind=DP)
 err_w=0
 !
 !$omp parallel default(shared) private(iw,chib,chi0,f,Pm,rc,ierr)
 allocate(chib(ng,ng),chi0(ng,ng),f(ng,ng),Pm(ng,ng))
 !$omp do schedule(dynamic)
 do iw=1,nw
   call TDDFT_Chi_pair_response(zz(iw),H0,Chi_dyn_V(:,:,1),Chi_dyn_V(:,:,2),Chi_dyn_s,Chi_dyn_c,chib,ierr)
   if (ierr==0) then
     call TDDFT_Chi_pair_bubble(zz(iw),E_ref,Chi_dyn_V(:,:,1),Chi_dyn_V(:,:,2),Chi_dyn_s,Chi_dyn_c,chi0)
     call TDDFT_Chi_kernel_from_response(chi0,chib,d,.TRUE.,f,Pm,rc,ierr)
   endif
   err_w(iw)=ierr
   if (ierr==0) then
     Chi_0(:,:,iw)=cmplx(chi0,kind=SP)
     Chi_P(:,:,iw)=cmplx(Pm,kind=SP)
     Chi_fxc(:,:,iw)=cmplx(f,kind=SP)
     Chi_rcond(:,iw)=real(rc,SP)
   endif
 enddo
 !$omp end do
 deallocate(chib,chi0,f,Pm)
 !$omp end parallel
 !
 n_bad=count(err_w/=0)
 if (n_bad>0) then
   do iw=1,nw
     if (err_w(iw)/=0) exit
   enddo
   call error('[Chi/BSK] BSEChiOut failed at '//trim(intc(n_bad))//' frequencies; first at '//&
&    trim(real2ch(real(real(zz(iw))*HA2EV,SP)))//' eV with code '//trim(intc(err_w(iw)))//' (31/32: reduce BSENGfxc)')
 endif
 call msg('r','[Chi/BSK] minimum rcond (reference bubble, BSE response)',(/minval(Chi_rcond(1,:)),minval(Chi_rcond(2,:))/))
 if (minval(Chi_rcond)<Chi_rcond_floor) &
&  call error('[Chi/BSK] BSEChiOut: ill-conditioned response; reduce BSENGexx/BSENGfxc for this band window')
 !
 ! The static kernel must be Hermitian (on the scale |q+G| fxc |q+G'| where all blocks are O(1))
 !
 herm=0._DP
 scale_max=0._DP
 do jg=1,ng
   do ig=1,ng
     herm=max(herm,abs(Chi_fxc(ig,jg,1)-conjg(Chi_fxc(jg,ig,1)))/abs(d(ig)*d(jg)))
     scale_max=max(scale_max,abs(Chi_fxc(ig,jg,1))/abs(d(ig)*d(jg)))
   enddo
 enddo
 call msg('r','[Chi/BSK] static fxc head [Ha] and relative anti-Hermitian part',&
&         (/real(Chi_fxc(1,1,1)),real(herm/max(scale_max,tiny(1._DP)),SP)/))
 !
 ! ndb.Chi in this job's directory
 !
 Chi_ng=ng
 Chi_n_freqs=nw
 Chi_QP_applied=l_qp
 Chi_G_db=''
 Chi_G_axis='REAL'
 Chi_G0_error=-1._SP
 Chi_optical_q0=BSS_q0
 Chi_optical_q0_known=.TRUE.
 if (l_exc) then
   Chi_G_mode='BSE-EXC'
   Chi_QP_action='P from the coupled BSE; chi0 = bubble of its (QP) transition energies'
 else
   Chi_G_mode='BSE-KS'
   Chi_QP_action='P from the coupled BSE; chi0 = KS bubble'
 endif
 Chi_convention='Bare vertices; P^-1 = chib^-1 + vbar from the BSE; fxc=chi0^-1-P^-1, unsymmetrized [Ha]; w1 = 0'
 call X_reset(X_io)
 X_io%ib=BS_bands
 X_io%ng=ng
 call W_reset(W_io)
 call W_duplicate(W,W_io)
 call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),ID=ID)
 io_err=io_Chi(X_io,W_io,ID)
 if (io_err<0) call error('[Chi/BSK] cannot write the ndb.Chi header')
 call io_control(ACTION=OP_APP_CL,COM=REP,SEC=(/2*iq,2*iq+1/),ID=ID)
 io_err=io_Chi(X_io,W_io,ID)
 if (io_err<0) call error('[Chi/BSK] cannot write the ndb.Chi q fragment')
 call msg('r','[Chi/BSK] ndb.Chi written: G vectors, frequencies (z=0 first)',(/ng,nw/))
 !
 call W_reset(W_io)
 call Chi_free()
 call Chi_dyn_free()
 deallocate(E_ref,d,zz,err_w,H0)
 !
 call timing('ChiExport',OPR='stop')
 !
end subroutine K_Chi_export
