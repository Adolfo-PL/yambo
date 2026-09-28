!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
subroutine K_Chi_Kdump(iq,Xk,nh,H0)
 !
 ! BSEChiKoutW/BSEChiKoutB: the transition-space Casida matrix of the dynamic solver,
 ! written as text at the ndb.Chi frequencies z with Re z inside BSEChiKoutW.
 !
 !   o-<job>.Ktt_q<iq>_transitions : label, k point, bands, energy E_t and residual a_t
 !   o-<job>.Ktt_q<iq>_w<iw>        : t t' X_rr F_rr X_rc F_rc X_cr F_cr X_cc F_cc   [eV]
 !
 ! With coupling the matrix of K_Chi_dynamic is, block by block,
 !
 !   M(z) = [[ diag(E) + X_rr + F_rr ,           X_rc + F_rc ],
 !           [           X_cr + F_cr , -diag(E) + X_cc + F_cc ]]
 !
 ! X = the static part of the BSE matrix without the energies (exchange, and anything
 ! else the kernel holds), F = the fxc(z) part: F_rr = Kf_rr, F_rc = i Kf_rc,
 ! F_cr = i Kf_cr, F_cc = -Kf_cc with Kf_xy = Co/(4 pi) s V_x^H fxc(z) V_y s. The
 ! response is Resp_ii(z) = -Co_factor B^T [z - M(z)]^-1 A with A = (a, i b),
 ! B = (b, i a), b = conj(a). In TDA only the rr block is written.
 !
 ! BSEChiKoutB restricts the file to one valence/conduction pair: the rows and columns
 ! of those transitions (a sub-matrix, not a separate Casida problem).
 !
 use pars,           ONLY:SP,DP,lchlen,schlen
 use units,          ONLY:HA2EV
 use R_lattice,      ONLY:bz_samp
 use vec_operate,    ONLY:c2a
 use com,            ONLY:msg
 use stderr,         ONLY:intc
 use IO_int,         ONLY:get_name
 use BS,             ONLY:BS_K_coupling,BS_K_dim,BS_K_Chi_dump_w,BS_K_Chi_dump_b
 use BS_solvers,     ONLY:BSS_eh_E,BSS_eh_table,BSS_eh_f,BSS_dipoles_opt
 use Chi_m,          ONLY:Chi_freqs,Chi_fxc,Chi_n_freqs,Chi_ng,Chi_dyn_V,Chi_dyn_s,Chi_dyn_c
 use TDDFT_Chi_dynamic_m, ONLY:TDDFT_Chi_dynamic_kernel
 !
 implicit none
 !
 integer       :: iq,nh
 type(bz_samp) :: Xk
 complex(DP)   :: H0(nh,nh)
 !
 ! Work Space
 !
 integer                  :: n,m,it,jt,iw,n_w,ierr,iu,ib,nb
 integer,     allocatable :: sel(:)
 real(SP)                 :: k_rlu(3)
 real(DP)                 :: E_t
 real(DP),    allocatable :: s(:)
 complex(DP), allocatable :: Vr(:,:),Vc(:,:),F(:,:),X(:,:,:),Kf(:,:,:)
 complex(SP)              :: a_t
 character(lchlen)        :: fname
 complex(DP), parameter   :: ci=(0._DP,1._DP)
 !
 if (BS_K_Chi_dump_w(2)<0._SP) return
 !
 n=BS_K_dim(1)
 nb=1
 if (BS_K_coupling) nb=4
 if ((BS_K_coupling.and.nh/=2*n).or.(.not.BS_K_coupling.and.nh/=n)) &
&  call error('[Chi/BSK] BSEChiKoutW: inconsistent BSE matrix size')
 !
 ! Transitions kept: one valence/conduction pair or all
 !
 m=0
 do it=1,n
   if (keep(it)) m=m+1
 enddo
 if (m==0) then
   call warning('[Chi/BSK] BSEChiKoutB: no transition with valence band '//trim(intc(BS_K_Chi_dump_b(1)))//&
&               ' and conduction band '//trim(intc(BS_K_Chi_dump_b(2))))
   return
 endif
 allocate(sel(m))
 m=0
 do it=1,n
   if (.not.keep(it)) cycle
   m=m+1
   sel(m)=it
 enddo
 !
 n_w=0
 do iw=1,Chi_n_freqs
   if (in_window(iw)) n_w=n_w+1
 enddo
 if (n_w==0) then
   call warning('[Chi/BSK] BSEChiKoutW: no ndb.Chi frequency in the window, no Casida matrix written')
   return
 endif
 call msg('r','[Chi/BSK] Casida matrix written: transitions, frequencies',(/m,n_w/))
 if (real(m,DP)**2*n_w>1.E7_DP) call warning('[Chi/BSK] BSEChiKoutW: large text output ('//&
&  trim(intc(m))//'^2 lines per frequency); BSEChiKoutB selects one band pair')
 !
 ! Transitions
 !
 fname=get_name(desc='Ktt_q'//trim(intc(iq))//'_transitions',type=3,CORE_IO=.false.,MORE_IO=.false.,COM_IO=.true.)
 open(newunit=iu,file=trim(fname))
 write(iu,'(a)') '# Transitions of the Casida matrix (BSEChiKoutW). E_t [eV]: reference transition energy;'
 write(iu,'(a)') '# a_t = d_t sqrt(f_t): resonant residual (b_t = conj(a_t)); k in reduced (rlu) units'
 write(iu,'(a)') '#    t    T  ik_bz ik_ibz      k1        k2        k3     v    c        E_t [eV]'//&
&                '            f_t          Re a_t          Im a_t'
 do it=1,m
   jt=sel(it)
   k_rlu=0._SP
   if (allocated(Xk%ptbz)) call c2a(v_in=Xk%ptbz(BSS_eh_table(jt,1),:),v_out=k_rlu,mode='ki2a')
   a_t=BSS_dipoles_opt(1,jt)*sqrt(cmplx(BSS_eh_f(jt),kind=SP))
   write(iu,'(4i6,3f10.5,2i5,4es16.8)') it,jt,BSS_eh_table(jt,1),Xk%sstar(BSS_eh_table(jt,1),1),k_rlu,&
&        BSS_eh_table(jt,2),BSS_eh_table(jt,3),real(BSS_eh_E(jt),DP)*HA2EV,real(BSS_eh_f(jt),DP),&
&        real(a_t,DP),aimag(a_t)
 enddo
 close(iu)
 !
 ! Static part without the energies
 !
 allocate(X(m,m,nb),Kf(m,m,nb),F(Chi_ng,Chi_ng),Vr(Chi_ng,m),Vc(Chi_ng,m),s(m))
 do it=1,m
   Vr(:,it)=Chi_dyn_V(:,sel(it),1)
   Vc(:,it)=Chi_dyn_V(:,sel(it),2)
   s(it)=Chi_dyn_s(sel(it))
   do jt=1,m
     X(it,jt,1)=H0(sel(it),sel(jt))
     if (BS_K_coupling) then
       X(it,jt,2)=H0(sel(it),n+sel(jt))
       X(it,jt,3)=H0(n+sel(it),sel(jt))
       X(it,jt,4)=H0(n+sel(it),n+sel(jt))
     endif
   enddo
   E_t=real(BSS_eh_E(sel(it)),DP)
   X(it,it,1)=X(it,it,1)-E_t
   if (BS_K_coupling) X(it,it,4)=X(it,it,4)+E_t
 enddo
 !
 ! fxc(z) part at every frequency of the window
 !
 do iw=1,Chi_n_freqs
   if (.not.in_window(iw)) cycle
   F=cmplx(Chi_fxc(:,:,iw),kind=DP)
   call TDDFT_Chi_dynamic_kernel(F,Vr,Vr,s,Chi_dyn_c,Kf(:,:,1),ierr)
   if (ierr==0.and.BS_K_coupling) then
     call TDDFT_Chi_dynamic_kernel(F,Vr,Vc,s,Chi_dyn_c,Kf(:,:,2),ierr)
     Kf(:,:,2)=ci*Kf(:,:,2)
     if (ierr==0) call TDDFT_Chi_dynamic_kernel(F,Vc,Vr,s,Chi_dyn_c,Kf(:,:,3),ierr)
     Kf(:,:,3)=ci*Kf(:,:,3)
     if (ierr==0) call TDDFT_Chi_dynamic_kernel(F,Vc,Vc,s,Chi_dyn_c,Kf(:,:,4),ierr)
     Kf(:,:,4)=-Kf(:,:,4)
   endif
   if (ierr/=0) call error('[Chi/BSK] BSEChiKoutW: kernel projection failed with code '//trim(intc(ierr)))
   !
   fname=get_name(desc='Ktt_q'//trim(intc(iq))//'_w'//trim(intc(iw)),type=3,CORE_IO=.false.,MORE_IO=.false.,COM_IO=.true.)
   open(newunit=iu,file=trim(fname))
   write(iu,'(a,i6,a,2f14.8,a)') '# Casida matrix at ndb.Chi frequency',iw,', z =',&
&        real(Chi_freqs(iw),DP)*HA2EV,aimag(Chi_freqs(iw))*real(HA2EV,DP),' eV'
   write(iu,'(a)') '# M(z) = [[diag(E)+X_rr+F_rr, X_rc+F_rc], [X_cr+F_cr, -diag(E)+X_cc+F_cc]]  [eV]'
   write(iu,'(a)') '# X: static part without the energies (exchange); F: fxc(z) part, as it enters M'
   write(iu,'(a)') '# t, t'' label rows/columns of o-<job>.Ktt_q'//trim(intc(iq))//'_transitions'
   if (BS_K_coupling) then
     write(iu,'(a)') '#    t     t''    Re X_rr    Im X_rr    Re F_rr    Im F_rr    Re X_rc    Im X_rc'//&
&                    '    Re F_rc    Im F_rc    Re X_cr    Im X_cr    Re F_cr    Im F_cr    Re X_cc    Im X_cc'//&
&                    '    Re F_cc    Im F_cc'
   else
     write(iu,'(a)') '# TDA: resonant block only'
     write(iu,'(a)') '#    t     t''    Re X_rr    Im X_rr    Re F_rr    Im F_rr'
   endif
   do it=1,m
     do jt=1,m
       write(iu,'(2i6,16es16.8)') it,jt,(real(X(it,jt,ib))*HA2EV,aimag(X(it,jt,ib))*HA2EV,&
&                                       real(Kf(it,jt,ib))*HA2EV,aimag(Kf(it,jt,ib))*HA2EV,ib=1,nb)
     enddo
   enddo
   close(iu)
 enddo
 !
 deallocate(sel,X,Kf,F,Vr,Vc,s)
 !
 contains
   !
   logical function keep(i_t)
     integer, intent(in) :: i_t
     keep=(BS_K_Chi_dump_b(1)==0.or.BSS_eh_table(i_t,2)==BS_K_Chi_dump_b(1)).and.&
&         (BS_K_Chi_dump_b(2)==0.or.BSS_eh_table(i_t,3)==BS_K_Chi_dump_b(2))
   end function keep
   !
   logical function in_window(i_w)
     integer, intent(in) :: i_w
     real(SP) :: w_re,tol
     w_re=real(Chi_freqs(i_w))
     tol=1.E-6_SP*max(1._SP,abs(BS_K_Chi_dump_w(2)))
     in_window=w_re>=BS_K_Chi_dump_w(1)-tol.and.w_re<=BS_K_Chi_dump_w(2)+tol
   end function in_window
   !
end subroutine K_Chi_Kdump
