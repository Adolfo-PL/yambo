!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
subroutine K_Chi_dynamic(iq,W)
 !
 ! Casida response with the frequency-dependent kernel fxc(G,G';w) of a real-axis
 ! ndb.Chi (BSKmod="CHI" + BSEChiDyn). The static part of the BSE matrix (transition
 ! energies and exchange) comes from BS_mat as for the diagonalization; K_kernel
 ! left fxc out of it and kept the transition vertices (Chi_dyn_*). At every BSE
 ! frequency z = w + i eta the fxc(z) of the same ndb.Chi sample is projected and
 !
 !   Resp_ii(w) = -Co_factor * B^T [z - M(z)]^-1 A
 !
 ! is solved directly (TDDFT_Chi_dynamic_m). With coupling, the same bands and the
 ! same G basis as the ndb.Chi this reproduces its QP RPA+LF response at every w.
 ! For a z-independent Hermitian fxc it equals the diagonalization result.
 !
 use pars,           ONLY:SP,DP,cZERO,schlen
 use units,          ONLY:HA2EV
 use frequency,      ONLY:w_samp
 use com,            ONLY:msg,isec
 use stderr,         ONLY:intc,real2ch
 use parallel_m,     ONLY:ncpu
 use BS,             ONLY:BS_K_coupling,BS_K_dim,BS_H_dim,BS_res_ares_n_mat,l_BS_abs,l_BS_kerr,&
&                         l_BS_magnons,l_BS_dichroism,l_BS_photolum
 use BS_solvers,     ONLY:BS_mat,BSS_dipoles_opt,BSS_eh_f,BSS_eh_W,BSS_eh_Z,Co_factor,para_term_w0,BSS_desc
 use X_m,            ONLY:Resp_ii,X_OUTPUT_driver,global_gauge
 use Chi_m,          ONLY:Chi_freqs,Chi_fxc,Chi_n_freqs,Chi_ng,Chi_free,&
&                         Chi_dyn_V,Chi_dyn_s,Chi_dyn_filled,Chi_dyn_c,Chi_dyn_free
 use descriptors,    ONLY:IO_desc_add
 use TDDFT_Chi_dynamic_m, ONLY:TDDFT_Chi_dynamic_response
 use timing_m,       ONLY:timing
 !
 implicit none
 !
 integer      :: iq
 type(w_samp) :: W
 !
 ! Work Space
 !
 character(1)             :: sec
 integer                  :: n,nh,nw,iw,it,ierr,iw_bad,n_bad
 integer,     allocatable :: err_w(:)
 real(SP)                 :: tol
 complex(DP)              :: resp
 complex(DP), allocatable :: H0(:,:),F(:,:),a(:),b(:),resp_w(:)
 complex(SP)              :: sqrt_f
 real(DP)                 :: ares_err
 !
 call timing('ChiDyn Solver',OPR='start')
 !
 if (isec(2)/=0) sec="="
 if (isec(2)==0) sec="+"
 call section(sec,'Dynamic Chi solver @q'//trim(intc(iq)))
 !
 call X_OUTPUT_driver("DEF INIT",solver="ChiDyn")
 !
 ! Checks
 !========
 if (ncpu/=1) call error('[Chi/BSK] the dynamic solver is serial (MPI); use OpenMP threads')
 if (.not.allocated(Chi_dyn_V).or..not.allocated(Chi_fxc)) &
&  call error('[Chi/BSK] no transition vertices or fxc(w): the dynamic kernel must be built in this run')
 if (.not.all(Chi_dyn_filled)) call error('[Chi/BSK] some transitions have no fxc vertex')
 if (BS_res_ares_n_mat/=1) call error('[Chi/BSK] dynamic solver: one BSE matrix only')
 if (trim(global_gauge)/='length') call error('[Chi/BSK] dynamic solver: length gauge only')
 if (allocated(BSS_eh_W).or.allocated(BSS_eh_Z)) &
&  call error('[Chi/BSK] dynamic solver: transition widths/Z factors are not supported')
 if (.not.l_BS_abs.or.l_BS_kerr.or.l_BS_magnons.or.l_BS_dichroism.or.l_BS_photolum) &
&  call error('[Chi/BSK] dynamic solver: set BSEprop= "abs" only')
 !
 ! The BSE frequencies must be the ndb.Chi ones: fxc is not interpolated in w
 !
 nw=W%n_freqs
 if (any(aimag(Chi_freqs)<=0._SP)) &
&  call error('[Chi/BSK] ndb.Chi is not on the damped real axis: rebuild it with ChiGAxis="REAL"')
 iw_bad=0
 if (nw/=Chi_n_freqs) then
   iw_bad=-1
 else
   do iw=1,nw
     tol=max(1.E-6_SP,1.E-5_SP*abs(Chi_freqs(iw)))
     if (abs(W%p(iw)-Chi_freqs(iw))>tol) then
       iw_bad=iw
       exit
     endif
   enddo
 endif
 if (iw_bad/=0) then
   call msg('r','[Chi/BSK] ndb.Chi frequencies (n, first/last Re [eV], first/last damping [eV])')
   call msg('r','[Chi/BSK]   n=',Chi_n_freqs)
   call msg('r','[Chi/BSK]   Re w  ',(/real(Chi_freqs(1))*HA2EV,real(Chi_freqs(Chi_n_freqs))*HA2EV/))
   call msg('r','[Chi/BSK]   eta   ',(/aimag(Chi_freqs(1))*HA2EV,aimag(Chi_freqs(Chi_n_freqs))*HA2EV/))
   call error('[Chi/BSK] BSE and ndb.Chi frequencies differ: set BEnRange= '//&
&    trim(real2ch(real(Chi_freqs(1))*HA2EV))//' | '//trim(real2ch(real(Chi_freqs(Chi_n_freqs))*HA2EV))//&
&    ' eV, BEnSteps= '//trim(intc(Chi_n_freqs))//', BDmRange= '//trim(real2ch(aimag(Chi_freqs(1))*HA2EV))//&
&    ' | '//trim(real2ch(aimag(Chi_freqs(Chi_n_freqs))*HA2EV))//' eV')
 endif
 !
 ! Static matrix: energies and exchange (fxc was left out by K_kernel)
 !====================================================================
 call K_components_folded_in_serial_arrays(iq)
 call K_stored_in_a_big_matrix(1,iq,"DIAGO")
 !
 n=BS_K_dim(1)
 nh=n
 if (BS_K_coupling) nh=BS_H_dim
 if (size(BS_mat,1)/=nh.or.size(Chi_dyn_s)/=n) call error('[Chi/BSK] dynamic solver: inconsistent BSE dimensions')
 !
 allocate(H0(nh,nh),a(n),b(n),resp_w(nw),err_w(nw))
 H0=cmplx(BS_mat,kind=DP)
 !
 ! Residual vectors as in K_diago_left/right_residuals: a = d sqrt(f), b = conj(d) sqrt(f)
 !
 do it=1,n
   sqrt_f=sqrt(cmplx(BSS_eh_f(it),kind=SP))
   a(it)=cmplx(BSS_dipoles_opt(1,it)*sqrt_f,kind=DP)
   b(it)=cmplx(conjg(BSS_dipoles_opt(1,it))*sqrt_f,kind=DP)
 enddo
 if (BS_K_coupling) then
   ! The solver builds the anti-resonant vectors as (i b, i a); Yambo's own must agree
   ares_err=0._DP
   do it=1,n
     sqrt_f=sqrt(cmplx(BSS_eh_f(n+it),kind=SP))
     ares_err=max(ares_err,abs(cmplx(BSS_dipoles_opt(1,n+it)*sqrt_f,kind=DP)-(0._DP,1._DP)*b(it)),&
&                          abs(cmplx(conjg(BSS_dipoles_opt(1,n+it))*sqrt_f,kind=DP)-(0._DP,1._DP)*a(it)))
   enddo
   if (ares_err>1.E-4_DP*max(maxval(abs(a)),tiny(1._DP))) &
&    call error('[Chi/BSK] anti-resonant dipoles differ from the (i b, i a) convention of the dynamic solver')
 endif
 !
 call msg('r','[Chi/BSK] dynamic Casida: transitions',n)
 ! fxc(w) of ndb.Chi inverts the full (resonant + anti-resonant) bubbles. Without the
 ! coupling blocks the KS poles are not cancelled: weights and positions come out wrong.
 if (.not.BS_K_coupling) call warning('[Chi/BSK] TDA with fxc(w) is inconsistent with ndb.Chi; use BSEmod= "coupling"')
 call msg('r','[Chi/BSK] dynamic Casida: frequencies (fxc taken on the ndb.Chi grid)',nw)
 call msg('r','[Chi/BSK] fxc head at first/last frequency [Ha]',&
&         (/real(Chi_fxc(1,1,1)),aimag(Chi_fxc(1,1,1)),real(Chi_fxc(1,1,nw)),aimag(Chi_fxc(1,1,nw))/))
 !
 ! One linear solve per frequency (threads over frequencies)
 !===========================================================
 err_w=0
 resp_w=(0._DP,0._DP)
 !$omp parallel default(shared) private(iw,F,resp,ierr)
 allocate(F(Chi_ng,Chi_ng))
 !$omp do schedule(dynamic)
 do iw=1,nw
   F=cmplx(Chi_fxc(:,:,iw),kind=DP)
   call TDDFT_Chi_dynamic_response(cmplx(W%p(iw),kind=DP),H0,F,Chi_dyn_V(:,:,1),Chi_dyn_V(:,:,2),&
&                                  Chi_dyn_s,Chi_dyn_c,a,b,BS_K_coupling,resp,ierr)
   resp_w(iw)=resp
   err_w(iw)=ierr
 enddo
 !$omp end do
 deallocate(F)
 !$omp end parallel
 !
 n_bad=count(err_w/=0)
 if (n_bad>0) then
   do iw=1,nw
     if (err_w(iw)/=0) exit
   enddo
   call error('[Chi/BSK] dynamic solve failed at '//trim(intc(n_bad))//' frequencies; first at '//&
&    trim(real2ch(real(W%p(iw))*HA2EV))//' eV with code '//trim(intc(err_w(iw))))
 endif
 !
 Resp_ii(:,2)=cmplx(-real(Co_factor,DP)*resp_w,kind=SP)
 !
 ! The f-sum weight A of the header has no meaning for a w-dependent kernel:
 ! carry the IP value so that the Drude/diamagnetic bookkeeping stays neutral.
 para_term_w0(2)=para_term_w0(1)
 call msg('r','[Chi/BSK] diamagnetic "A" of the output carries the IP value (undefined for fxc(w))')
 !
 call IO_desc_add(BSS_desc,kind="s",str="[Chi/BSK] Kernel",S="fxc(w) from ndb.Chi, per-frequency solve")
 call K_OUTPUT(iq,W,"open dump close ChiDyn","ALL")
 !
 deallocate(H0,a,b,resp_w,err_w)
 call Chi_dyn_free()
 call Chi_free()
 !
 call timing('ChiDyn Solver',OPR='stop')
 !
end subroutine K_Chi_dynamic
