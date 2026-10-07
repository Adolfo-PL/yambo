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
 !   o-<job>.Ktt_q<iq>_transitions   : label, k point, bands, energy E_t and residual a_t
 !   o-<job>.Ktt_q<iq>_antiresonant  : the same for row N+t (explicit anti-resonant
 !                                     transitions, finite q, coupling only), with A and B
 !   o-<job>.Ktt_q<iq>_w<iw>          : t t' X_rr F_rr X_rc F_rc X_cr F_cr X_cc F_cc   [eV]
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
 ! At finite q the anti-resonant rows are explicit transitions (K_driver_init): row N+t
 ! is anti-resonant transition t (electron in v at k, hole in c at k-q), with its own
 ! energy E_N+t < 0, occupation f < 0, vertex and oscillator. Then, in every block,
 ! F_xy = Co/(4 pi) s_x s_y V_x^H fxc(z) V_y with s = sqrt(f) (i sqrt|f| for f < 0), the
 ! cc block carries diag(E_N+t) instead of -diag(E), and A = s d, B = s conj(d) row by
 ! row (both listed in the _antiresonant file). With derived rows (q=0) these are the
 ! formulas above: s_N+t = i s_t, d_N+t = conj(d_t), E_N+t = -E_t.
 !
 ! BSEChiKoutB restricts the file to one valence/conduction pair: the rows and columns
 ! of those transitions (a sub-matrix, not a separate Casida problem).
 !
 use pars,           ONLY:SP,DP,lchlen,schlen
 use units,          ONLY:HA2EV
 use R_lattice,      ONLY:bz_samp,rl_sop
 use vec_operate,    ONLY:c2a
 use com,            ONLY:msg
 use stderr,         ONLY:intc
 use IO_int,         ONLY:get_name
 use BS,             ONLY:BS_K_coupling,BS_K_dim,BS_K_Chi_dump_w,BS_K_Chi_dump_b,l_BS_ares_from_res
 use BS_solvers,     ONLY:BSS_eh_E,BSS_eh_table,BSS_eh_f,BSS_dipoles_opt
 use Chi_m,          ONLY:Chi_freqs,Chi_fxc,Chi_n_freqs,Chi_ng,Chi_dyn_V,Chi_dyn_c
 !
 implicit none
 !
 integer       :: iq,nh
 type(bz_samp) :: Xk
 complex(DP)   :: H0(nh,nh)
 !
 ! Work Space
 !
 integer                  :: n,m,ma,it,jt,iw,n_w,iu,ib,nb,mm
 integer,     allocatable :: sel(:),sela(:),rows(:)
 real(DP),    allocatable :: E_row(:)
 complex(DP), allocatable :: V(:,:),sv(:),F(:,:),X(:,:,:),Kf(:,:,:),FV(:,:),K2(:,:)
 complex(SP)              :: a_t
 logical                  :: l_explicit
 character(lchlen)        :: fname
 !
 if (BS_K_Chi_dump_w(2)<0._SP) return
 !
 l_explicit=BS_K_coupling.and..not.l_BS_ares_from_res
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
 allocate(sel(m),sela(m))
 m=0
 do it=1,n
   if (.not.keep(it)) cycle
   m=m+1
   sel(m)=it
 enddo
 ! Anti-resonant rows: the mirrors of sel (q=0), or the explicit anti-resonant transitions
 ! of the same band pair (finite q; their table holds c in the v slot and v in the c slot)
 sela=sel
 if (l_explicit) then
   ma=0
   do it=1,n
     if (.not.keep_ares(n+it)) cycle
     ma=ma+1
     if (ma<=m) sela(ma)=it
   enddo
   if (ma/=m) call error('[Chi/BSK] BSEChiKoutB: resonant and anti-resonant transitions of the pair differ in number')
 endif
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
 if (iq>1) write(iu,'(a)') '# finite q: transition t is v at k-q -> c at k; k is the conduction-band k'
 do it=1,m
   jt=sel(it)
   a_t=BSS_dipoles_opt(1,jt)*sqrt(cmplx(BSS_eh_f(jt),kind=SP))
   call write_row(it,jt,BSS_eh_table(jt,2),BSS_eh_table(jt,3),a_t,.FALSE.)
 enddo
 close(iu)
 if (l_explicit) then
   fname=get_name(desc='Ktt_q'//trim(intc(iq))//'_antiresonant',type=3,CORE_IO=.false.,MORE_IO=.false.,COM_IO=.true.)
   open(newunit=iu,file=trim(fname))
   write(iu,'(a)') '# Anti-resonant rows N+t of the Casida matrix (explicit transitions, finite q).'
   write(iu,'(a)') '# Transition t: electron in v at k, hole in c at k-q (k as listed); E_t < 0, f_t < 0.'
   write(iu,'(a)') '# A_t = d_t sqrt(f_t), B_t = conj(d_t) sqrt(f_t), sqrt(f) = i sqrt|f|: Resp = -Co B^T [z-M]^-1 A'
   write(iu,'(a)') '#    t    T  ik_bz ik_ibz      k1        k2        k3     v    c        E_t [eV]'//&
&                  '            f_t          Re A_t          Im A_t          Re B_t          Im B_t'
   do it=1,m
     jt=n+sela(it)
     a_t=BSS_dipoles_opt(1,jt)*sqrt(cmplx(BSS_eh_f(jt),kind=SP))
     call write_row(it,jt,BSS_eh_table(jt,3),BSS_eh_table(jt,2),a_t,.TRUE.)
   enddo
   close(iu)
 endif
 !
 ! Static part without the energies
 !
 ! Rows of BS_mat kept: sel (resonant) and, with coupling, N+sela (anti-resonant). Every
 ! row has its vertex V, s = sqrt(f) (complex) and energy; F = Co/(4 pi) s V^H fxc V s.
 mm=m
 if (BS_K_coupling) mm=2*m
 allocate(X(m,m,nb),Kf(m,m,nb),F(Chi_ng,Chi_ng),V(Chi_ng,mm),sv(mm),rows(mm),E_row(mm),FV(Chi_ng,mm),K2(mm,mm))
 rows(:m)=sel
 if (BS_K_coupling) rows(m+1:)=n+sela
 do it=1,mm
   if (rows(it)<=n) V(:,it)=Chi_dyn_V(:,rows(it),1)
   if (rows(it)> n) V(:,it)=Chi_dyn_V(:,rows(it)-n,2)
   sv(it)=sqrt(cmplx(BSS_eh_f(rows(it)),kind=DP))
   E_row(it)=real(BSS_eh_E(rows(it)),DP)
 enddo
 do ib=1,nb
   do jt=1,m
     do it=1,m
       X(it,jt,ib)=H0(rows(it+m*((ib-1)/2)),rows(jt+m*mod(ib-1,2)))
     enddo
   enddo
 enddo
 do it=1,m
   X(it,it,1)=X(it,it,1)-E_row(it)
   if (BS_K_coupling) X(it,it,4)=X(it,it,4)-E_row(m+it)
 enddo
 !
 ! fxc(z) part at every frequency of the window
 !
 do iw=1,Chi_n_freqs
   if (.not.in_window(iw)) cycle
   F=cmplx(Chi_fxc(:,:,iw),kind=DP)
   FV=matmul(F,V)
   K2=Chi_dyn_c*matmul(conjg(transpose(V)),FV)
   do it=1,mm
     K2(:,it)=K2(:,it)*sv(it)
     K2(it,:)=K2(it,:)*sv(it)
   enddo
   do ib=1,nb
     Kf(:,:,ib)=K2(1+m*((ib-1)/2):m*(1+(ib-1)/2),1+m*mod(ib-1,2):m*(1+mod(ib-1,2)))
   enddo
   !
   fname=get_name(desc='Ktt_q'//trim(intc(iq))//'_w'//trim(intc(iw)),type=3,CORE_IO=.false.,MORE_IO=.false.,COM_IO=.true.)
   open(newunit=iu,file=trim(fname))
   write(iu,'(a,i6,a,2f14.8,a)') '# Casida matrix at ndb.Chi frequency',iw,', z =',&
&        real(Chi_freqs(iw),DP)*HA2EV,aimag(Chi_freqs(iw))*real(HA2EV,DP),' eV'
   if (.not.l_explicit) write(iu,'(a)') '# M(z) = [[diag(E)+X_rr+F_rr, X_rc+F_rc], [X_cr+F_cr, -diag(E)+X_cc+F_cc]]  [eV]'
   if (l_explicit) write(iu,'(a)') '# M(z) = [[diag(E)+X_rr+F_rr, X_rc+F_rc], [X_cr+F_cr, diag(E_ares)+X_cc+F_cc]]  [eV]'//&
&    ' (E_ares, t'' of the c blocks: o-<job>.Ktt_q'//trim(intc(iq))//'_antiresonant)'
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
 deallocate(sel,sela,X,Kf,F,V,sv,rows,E_row,FV,K2)
 !
 contains
   !
   logical function keep(i_t)
     integer, intent(in) :: i_t
     keep=(BS_K_Chi_dump_b(1)==0.or.BSS_eh_table(i_t,2)==BS_K_Chi_dump_b(1)).and.&
&         (BS_K_Chi_dump_b(2)==0.or.BSS_eh_table(i_t,3)==BS_K_Chi_dump_b(2))
   end function keep
   !
   logical function keep_ares(i_t)
     ! explicit anti-resonant row: the table holds c in the v slot, v in the c slot
     integer, intent(in) :: i_t
     keep_ares=(BS_K_Chi_dump_b(1)==0.or.BSS_eh_table(i_t,3)==BS_K_Chi_dump_b(1)).and.&
&              (BS_K_Chi_dump_b(2)==0.or.BSS_eh_table(i_t,2)==BS_K_Chi_dump_b(2))
   end function keep_ares
   !
   subroutine write_row(i_t,j_t,i_v,i_c,a_r,l_ab)
     ! BZ point from its IBZ image and star symmetry, as k_ibz2bz (Xk%ptbz is not set here)
     integer,     intent(in) :: i_t,j_t,i_v,i_c
     complex(SP), intent(in) :: a_r
     logical,     intent(in) :: l_ab
     real(SP)    :: k_rlu(3),k_iku(3)
     complex(SP) :: b_r
     k_iku=matmul(rl_sop(:,:,Xk%sstar(BSS_eh_table(j_t,1),2)),Xk%pt(Xk%sstar(BSS_eh_table(j_t,1),1),:))
     call c2a(v_in=k_iku,v_out=k_rlu,mode='ki2a')
     if (.not.l_ab) then
       write(iu,'(4i6,3f10.5,2i5,4es16.8)') i_t,j_t,BSS_eh_table(j_t,1),Xk%sstar(BSS_eh_table(j_t,1),1),k_rlu,&
&            i_v,i_c,real(BSS_eh_E(j_t),DP)*HA2EV,real(BSS_eh_f(j_t),DP),real(a_r,DP),aimag(a_r)
     else
       b_r=conjg(BSS_dipoles_opt(1,j_t))*sqrt(cmplx(BSS_eh_f(j_t),kind=SP))
       write(iu,'(4i6,3f10.5,2i5,6es16.8)') i_t,j_t,BSS_eh_table(j_t,1),Xk%sstar(BSS_eh_table(j_t,1),1),k_rlu,&
&            i_v,i_c,real(BSS_eh_E(j_t),DP)*HA2EV,real(BSS_eh_f(j_t),DP),real(a_r,DP),aimag(a_r),&
&            real(b_r,DP),aimag(b_r)
     endif
   end subroutine write_row
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
