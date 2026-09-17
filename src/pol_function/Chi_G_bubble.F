!
! License-Identifier: GPL
!
subroutine Chi_G_bubble(iq,X,Xk,Xw,is_G0,response)
 ! Full band-pair density bubble. qindx_X maps k to k-q. scatter_Bamp
 ! supplies the same rho*(G) rho(G') convention as X_irredux_residuals.
 ! This first implementation deliberately requires serial CPU, finite q,
 ! an insulating KS reference, and no transition-energy filtering/terminator.
 use pars, ONLY:SP,DP
 use electrons, ONLY:n_sp_pol,spin_occ
 use Chi_m, ONLY:Chi_KS_levels,Chi_G_mode,Chi_G_axis,Chi_G_energy,Chi_G_weight,Chi_G_occupation
 use Chi_Green_m, ONLY:Chi_G0_pair,Chi_spectral_pair
 use X_m, ONLY:X_t
 use frequency, ONLY:w_samp
 use R_lattice, ONLY:bz_samp,qindx_X
 use D_lattice, ONLY:DL_vol
 use collision_el, ONLY:elemental_collision,elemental_collision_alloc,elemental_collision_free
 use, intrinsic :: ieee_arithmetic, ONLY:ieee_is_finite
 implicit none
 integer, intent(in) :: iq
 type(X_t), intent(in) :: X
 type(bz_samp), intent(in) :: Xk
 type(w_samp), intent(in) :: Xw
 logical, intent(in) :: is_G0
 complex(SP), intent(out) :: response(X%ng,X%ng,Xw%n_freqs)
 type(elemental_collision) :: scattering
 complex(DP) :: pair(Xw%n_freqs),z(Xw%n_freqs),vertex(X%ng,X%ng)
 complex(DP), allocatable :: accumulated(:,:,:)
 real(DP) :: normalization,fa,fb
 integer :: ikbz,ikpbz,ik,ikp,is,isp,i_spin,a,b,iw,ig,igp,info
 ! Validate before using an incomplete symmetry/distributed mapping.
 if (iq<=1) call error('[Chi/G] finite q required; set QpntsRXd(1)>1')
 z=cmplx(Xw%p,kind=DP)
 if (.not.all(ieee_is_finite(real(z,DP))).or..not.all(ieee_is_finite(aimag(z))))&
&  call error('[Chi/G] response frequencies must be finite')
 if (any(aimag(z)<0._DP)) call error('[Chi/G] response frequencies must be in the upper half plane')
 if (trim(Chi_G_axis)=='REAL'.and.any(aimag(z)<=0._DP))&
&  call error('[Chi/G] use positive response damping with ChiGAxis="REAL"')
 normalization=1._DP/(real(Xk%nbz,DP)*real(DL_vol,DP))
 allocate(accumulated(X%ng,X%ng,Xw%n_freqs))
 accumulated=(0._DP,0._DP)
 call elemental_collision_alloc(scattering,NG=X%ng,CPU_ONLY=.TRUE.,TITLE='Chi/G')
 do i_spin=1,n_sp_pol
 do ikbz=1,Xk%nbz
   ikpbz=qindx_X(iq,ikbz,1)
   ik =Xk%sstar(ikbz,1); is =Xk%sstar(ikbz,2)
   ikp=Xk%sstar(ikpbz,1); isp=Xk%sstar(ikpbz,2)
   do b=X%ib(1),X%ib(2)
   do a=X%ib(1),X%ib(2)
     if (is_G0.or.trim(Chi_G_mode)=='G0') then
       fa=real(Chi_KS_levels%f(a,ik,i_spin),DP)
       fb=real(Chi_KS_levels%f(b,ikp,i_spin),DP)
       if (fa==fb) cycle
       do iw=1,Xw%n_freqs
         pair(iw)=Chi_G0_pair(z(iw),real(Chi_KS_levels%E(a,ik,i_spin),DP),&
&            real(Chi_KS_levels%E(b,ikp,i_spin),DP),fa,fb)
       enddo
     else
       call Chi_spectral_pair(z,Chi_G_energy(a,ik,i_spin,:),Chi_G_energy(b,ikp,i_spin,:),&
&        Chi_G_weight(a,ik,i_spin,:),Chi_G_weight(b,ikp,i_spin,:),&
&        Chi_G_occupation(a,ik,i_spin,:),Chi_G_occupation(b,ikp,i_spin,:),pair,info)
       if (info/=0) call error('[Chi/G] invalid spectral convolution or undamped pole')
       pair=pair*real(spin_occ,DP)
       if (maxval(abs(pair))==0._DP) cycle
     endif
     scattering%is=(/a,ik,is,i_spin/)
     scattering%os=(/b,ikp,isp,i_spin/)
     scattering%qs=(/qindx_X(iq,ikbz,2),iq,1/)
     call scatter_Bamp(scattering)
     do igp=1,X%ng
     do ig=1,X%ng
       vertex(ig,igp)=conjg(cmplx(scattering%rhotw(ig),kind=DP))*&
&                            cmplx(scattering%rhotw(igp),kind=DP)*normalization
     enddo
     enddo
     do iw=1,Xw%n_freqs
       accumulated(:,:,iw)=accumulated(:,:,iw)+pair(iw)*vertex
     enddo
   enddo
   enddo
 enddo
 enddo
 response=cmplx(accumulated,kind=SP)
 deallocate(accumulated)
 call elemental_collision_free(scattering)
end subroutine
!
subroutine Chi_G_load(X,Xk)
 ! Load conventional retarded Dyson spectra without QP fitting/remapping.
 use pars, ONLY:SP,DP
 use Chi_m, ONLY:Chi_KS_levels,Chi_G_db,Chi_G_energy,Chi_G_weight,Chi_G_occupation,&
&                 Chi_G_norm,Chi_G_norm_tol,Chi_G_mu,Chi_G_kbt
 use Chi_Green_m, ONLY:Chi_spectral_weights,Chi_fermi,Chi_Dyson
 use X_m, ONLY:X_t
 use R_lattice, ONLY:bz_samp
 use electrons, ONLY:spin,n_sp_pol
 use QP_m, ONLY:QP_t,QP_reset,QP_G_er,QP_G_dr,GF_is_causal,QP_G_Zoom_treshold
 use global_XC, ONLY:QP_DB_kind
 use IO_int, ONLY:io_control
 use IO_m, ONLY:OP_RD_CL,DUMP,REP
 use com, ONLY:msg
 use stderr, ONLY:intc
 implicit none
 type(X_t), intent(in) :: X
 type(bz_samp), intent(in) :: Xk
 type(QP_t) :: qp
 real(SP) :: saved_er(2),saved_dr(2),saved_zoom
 integer :: saved_kind,ID,io_err,i,a,ik,s,iw,n,info
 logical :: saved_causal
 logical, allocatable :: covered(:,:,:)
 real(DP) :: norm,reference_tolerance
 complex(DP), allocatable :: reconstructed(:)
 integer, external :: io_QP_and_GF
 if (allocated(Chi_G_energy)) return
 if (len_trim(Chi_G_db)==0) call error('[Chi/G] ChiGDb must name a GRetarded ndb.G database')
 ! The generic reader modifies GW globals; preserve them in this response run.
 saved_er=QP_G_er; saved_dr=QP_G_dr; saved_zoom=QP_G_Zoom_treshold
 saved_causal=GF_is_causal; saved_kind=QP_DB_kind
 call QP_reset(qp)
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2,3/),MODE=DUMP,ID=ID)
 io_err=io_QP_and_GF(trim(Chi_G_db),qp,ID)
 QP_G_er=saved_er; QP_G_dr=saved_dr; QP_G_Zoom_treshold=saved_zoom
 GF_is_causal=saved_causal; QP_DB_kind=saved_kind
 if (io_err/=0) call error('[Chi/G] cannot read '//trim(Chi_G_db))
 if (.not.qp%GreenF_retarded) call error('[Chi/G] rebuild ndb.G with GRetarded; legacy spectra are unsupported')
 if (qp%nk/=Xk%nibz) call error('[Chi/G] ndb.G and response must have identical k meshes')
 if (qp%nb<X%ib(2)) call error('[Chi/G] ndb.G does not cover the response band range')
 if (maxval(abs(qp%k-Xk%pt))>1.E-5_SP) call error('[Chi/G] ndb.G k-point ordering differs from X')
 if (.not.allocated(qp%E_bare)) call error('[Chi/G] ndb.G lacks KS reference energies')
 Chi_G_mu=qp%GreenF_mu
 if (abs(Chi_G_mu-Chi_KS_levels%E_Fermi)>1.E-5_SP) call error('[Chi/G] KS/G energy references differ')
 n=qp%GreenF_n_steps
 if (n<2) call error('[Chi/G] at least two spectral energy samples required')
 allocate(reconstructed(n))
 allocate(covered(X%ib(1):X%ib(2),Xk%nibz,n_sp_pol))
 covered=.FALSE.
 allocate(Chi_G_energy(X%ib(2),Xk%nibz,n_sp_pol,n))
 allocate(Chi_G_weight(X%ib(2),Xk%nibz,n_sp_pol,n))
 allocate(Chi_G_occupation(X%ib(2),Xk%nibz,n_sp_pol,n))
 allocate(Chi_G_norm(X%ib(2),Xk%nibz,n_sp_pol))
 Chi_G_energy=0._DP; Chi_G_weight=0._DP; Chi_G_occupation=0._DP; Chi_G_norm=0._SP
 reference_tolerance=1.E-5_DP
 do i=1,qp%n_states
   a=qp%table(i,1); ik=qp%table(i,3); s=spin(qp%table(i,:))
   if (qp%table(i,2)/=a) call error('[Chi/G] off-diagonal self-energy is unsupported')
   if (a<X%ib(1).or.a>X%ib(2)) cycle
   if (ik<1.or.ik>Xk%nibz.or.s<1.or.s>n_sp_pol) call error('[Chi/G] invalid state index')
   if (covered(a,ik,s)) call error('[Chi/G] duplicate spectral state in ndb.G')
   if (abs(real(qp%E_bare(i),DP)-real(Chi_KS_levels%E(a,ik,s),DP))>reference_tolerance)&
&      call error('[Chi/G] ndb.G starting energies differ from the saved KS levels')
   Chi_G_energy(a,ik,s,:)=real(qp%GreenF_W(i,:),DP)
   if (any(aimag(qp%GreenF_W(i,:))<=0._SP)) call error('[Chi/G] spectral G must have positive damping')
   do iw=1,n
     reconstructed(iw)=Chi_Dyson(cmplx(qp%GreenF_W(i,iw),kind=DP),&
&                                real(qp%E_bare(i),DP),cmplx(qp%S_total(i,iw),kind=DP))
   enddo
   if (maxval(abs(reconstructed-cmplx(qp%GreenF(i,:),kind=DP)))>&
&       100._DP*real(epsilon(1._SP),DP)*max(maxval(abs(reconstructed)),tiny(1._DP)))&
&       call error('[Chi/G] exported G is inconsistent with its self-energy and KS reference')
   call Chi_spectral_weights(Chi_G_energy(a,ik,s,:),cmplx(qp%GreenF(i,:),kind=DP),&
&                           Chi_G_weight(a,ik,s,:),norm,info)
   if (info/=0) call error('[Chi/G] invalid/noncausal spectrum at state '//trim(intc(i)))
   Chi_G_norm(a,ik,s)=real(norm,SP)
   if (abs(norm-1._DP)>real(Chi_G_norm_tol,DP)) then
     call msg('r','[Chi/G] spectral integral at rejected state',real(norm,SP))
     call error('[Chi/G] expand/refine GEnRnge/GEnSteps; spectral weight is not normalized')
   endif
   do iw=1,n
     Chi_G_occupation(a,ik,s,iw)=Chi_fermi(Chi_G_energy(a,ik,s,iw),real(Chi_G_mu,DP),real(Chi_G_kbt,DP))
   enddo
   covered(a,ik,s)=.TRUE.
 enddo
 if (.not.all(covered)) call error('[Chi/G] QPkrange must cover every response band, k point and spin')
 call msg('r','[Chi/G] spectral integral min/max',&
& (/minval(Chi_G_norm(X%ib(1):,:,:)),maxval(Chi_G_norm(X%ib(1):,:,:))/))
 deallocate(covered,reconstructed)
 call QP_reset(qp)
end subroutine
