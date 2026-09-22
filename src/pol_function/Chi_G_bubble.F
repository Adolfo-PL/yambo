!
! License-Identifier: GPL
!
subroutine Chi_G_bubble(iq,X,Xk,Xw,Dip,is_G0,response)
 ! Full band-pair density bubble. qindx_X maps k to k-q. scatter_Bamp
 ! supplies the same rho*(G) rho(G') convention as X_irredux_residuals.
 ! This implementation requires CPU, a k-only MPI layout,
 ! an insulating KS reference, and no transition-energy filtering/terminator.
 use pars, ONLY:SP,DP
 use electrons, ONLY:n_sp_pol,spin_occ
 use Chi_m, ONLY:Chi_KS_levels,Chi_G_mode,Chi_G_axis,Chi_G_energy,Chi_G_weight,Chi_G_occupation
 use Chi_m, ONLY:Chi_G_static_energy,Chi_G_static_occupation
 use Chi_Green_m, ONLY:Chi_G0_pair,Chi_spectral_pair
 use X_m, ONLY:X_t
 use frequency, ONLY:w_samp
 use R_lattice, ONLY:bz_samp,qindx_X,q0_def_norm,g_rot
 use D_lattice, ONLY:DL_vol,nsym,i_time_rev,sop_inv
 use DIPOLES, ONLY:DIPOLE_t,DIP_iR
 use vec_operate, ONLY:v_norm
 use collision_el, ONLY:elemental_collision,elemental_collision_alloc,elemental_collision_free
 use parallel_m, ONLY:PAR_IND_Xk_bz,PAR_COM_X_WORLD
 use parallel_int, ONLY:PP_redux_wait
 use, intrinsic :: ieee_arithmetic, ONLY:ieee_is_finite
 implicit none
 integer, intent(in) :: iq
 type(X_t), intent(in) :: X
 type(bz_samp), intent(in) :: Xk
 type(w_samp), intent(in) :: Xw
 type(DIPOLE_t), intent(in) :: Dip
 logical, intent(in) :: is_G0
 complex(SP), intent(out) :: response(X%ng,X%ng,Xw%n_freqs)
 type(elemental_collision) :: scattering
 complex(DP) :: pair(Xw%n_freqs),z(Xw%n_freqs),vertex(X%ng,X%ng)
 complex(SP) :: density(X%ng),dipole(3),field_dir(3)
 complex(DP), allocatable :: accumulated(:,:,:)
 real(DP) :: normalization,fa,fb
 integer :: ikbz,ikpbz,ik,ikp,is,isp,i_spin,a,b,iw,ig,igp,info,sop_inv_l
 ! Validate before using an incomplete symmetry/distributed mapping.
 if (iq<1) call error('[Chi/G] invalid q index')
 if (iq==1) then
   if (trim(Chi_G_mode)=='DYSON')&
&    call error('[Chi/G] optical Dyson bubble needs an intraband spectral/Drude treatment')
   if (v_norm(Dip%q0)<=0._SP) call error('[Chi/G] optical direction is undefined')
   field_dir=cmplx(Dip%q0/v_norm(Dip%q0)*q0_def_norm,kind=SP)
 endif
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
   if (.not.PAR_IND_Xk_bz%element_1D(ikbz)) cycle
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
     else if (trim(Chi_G_mode)=='COHSEX') then
       fa=Chi_G_static_occupation(a,ik,i_spin)*real(spin_occ,DP)
       fb=Chi_G_static_occupation(b,ikp,i_spin)*real(spin_occ,DP)
       if (fa==fb) cycle
       do iw=1,Xw%n_freqs
         pair(iw)=Chi_G0_pair(z(iw),Chi_G_static_energy(a,ik,i_spin),&
&                                  Chi_G_static_energy(b,ikp,i_spin),fa,fb)
       enddo
     else
       call Chi_spectral_pair(z,Chi_G_energy(a,ik,i_spin,:),Chi_G_energy(b,ikp,i_spin,:),&
&        Chi_G_weight(a,ik,i_spin,:),Chi_G_weight(b,ikp,i_spin,:),&
&        Chi_G_occupation(a,ik,i_spin,:),Chi_G_occupation(b,ikp,i_spin,:),pair,info)
       if (info/=0) call error('[Chi/G] invalid spectral convolution or undamped pole')
       pair=pair*real(spin_occ,DP)
       if (maxval(abs(pair))==0._DP) cycle
     endif
     if (iq==1) then
       ! Native X_irredux_residuals treats the head as the dipole limit of
       ! the density vertex.  The ordinary q=0 overlap vanishes for a/=b.
       ! G>0 uses the q=0 scattering vertex rotated from the IBZ k star.
       if (X%ng>1) then
         scattering%is=(/a,ik,1,i_spin/)
         scattering%os=(/b,ik,1,i_spin/)
         scattering%qs=(/1,1,1/)
         call scatter_Bamp(scattering)
         sop_inv_l=sop_inv(is)
         do ig=2,X%ng
           density(ig)=scattering%rhotw(g_rot(ig,sop_inv_l))
         enddo
         if (is>nsym/(i_time_rev+1)) density(2:)=conjg(density(2:))
       endif
       ! DIP_iR is stored only for conduction/valence band order.  A
       ! reverse pair uses the same saved dipole and the opposite conjugate
       ! head, because DIP_iR(v,c)=-conjg(DIP_iR(c,v)).
       if (.not.allocated(DIP_iR)) call error('[Chi/G] optical DIP_iR is unavailable')
       if (fa<fb) then
         call DIPOLE_rotate(a,b,ikbz,i_spin,'DIP_iR',Xk,dipole)
         density(1)=-conjg(dot_product(field_dir,dipole))
       else
         call DIPOLE_rotate(b,a,ikbz,i_spin,'DIP_iR',Xk,dipole)
         density(1)=dot_product(field_dir,dipole)
       endif
     else
       scattering%is=(/a,ik,is,i_spin/)
       scattering%os=(/b,ikp,isp,i_spin/)
       scattering%qs=(/qindx_X(iq,ikbz,2),iq,1/)
       call scatter_Bamp(scattering)
       density=scattering%rhotw(1:X%ng)
     endif
     do igp=1,X%ng
     do ig=1,X%ng
       vertex(ig,igp)=conjg(cmplx(density(ig),kind=DP))*&
&                            cmplx(density(igp),kind=DP)*normalization
     enddo
     enddo
     do iw=1,Xw%n_freqs
       accumulated(:,:,iw)=accumulated(:,:,iw)+pair(iw)*vertex
     enddo
   enddo
   enddo
 enddo
 enddo
 do iw=1,Xw%n_freqs
   call PP_redux_wait(accumulated(:,:,iw),COMM=PAR_COM_X_WORLD%COMM)
 enddo
 response=cmplx(accumulated,kind=SP)
 deallocate(accumulated)
 call elemental_collision_free(scattering)
end subroutine
!
subroutine Chi_G_parallel_check(X)
 ! The native k partition includes k-q partners in PAR_IND_Xk_ibz. With
 ! c/v/g/q roles equal to one, every local k has both complete band ranges.
 ! Other layouts would duplicate pairs or require unavailable wavefunctions.
 use X_m, ONLY:X_t
 use parallel_m, ONLY:ncpu,PAR_COM_Q_INDEX,PAR_COM_RL_INDEX,PAR_COM_X_WORLD,&
&                     PAR_COM_CON_INDEX_X,PAR_COM_VAL_INDEX_X,PAR_IND_Xk_bz
 implicit none
 type(X_t), intent(in) :: X
 if (PAR_COM_Q_INDEX%n_CPU/=1.or.PAR_COM_RL_INDEX%n_CPU/=1.or.&
&     PAR_COM_CON_INDEX_X(X%whoami)%n_CPU/=1.or.PAR_COM_VAL_INDEX_X(X%whoami)%n_CPU/=1.or.&
&     PAR_COM_X_WORLD%n_CPU/=ncpu)&
&   call error('[Chi/G] use k-only MPI: X_and_IO_CPU="1 1 N 1 1", X_and_IO_ROLEs="q g k c v"')
 if (.not.allocated(PAR_IND_Xk_bz%element_1D)) call error('[Chi/G] native BZ k ownership is undefined')
end subroutine
!
subroutine Chi_G_load_static(qp,X,Xk)
 ! A native COHSEX QP energy already includes SEX+COH-Vxc once. Its static
 ! Dyson G has a single unit-weight pole; no broadened spectral grid is needed.
 use pars, ONLY:SP,DP
 use QP_m, ONLY:QP_t
 use Chi_m, ONLY:Chi_KS_levels,Chi_G_mu,Chi_G_static_energy,Chi_G_static_occupation,Chi_G_norm
 use Chi_Green_m, ONLY:Chi_fermi
 use X_m, ONLY:X_t
 use R_lattice, ONLY:bz_samp
 use electrons, ONLY:n_sp_pol,spin,spin_occ
 use, intrinsic :: ieee_arithmetic, ONLY:ieee_is_finite
 implicit none
 type(QP_t), intent(in) :: qp
 type(X_t), intent(in) :: X
 type(bz_samp), intent(in) :: Xk
 logical, allocatable :: covered(:,:,:)
 real(DP) :: pole,occupation,tolerance
 integer :: i,a,ik,s
 if (.not.allocated(qp%E).or..not.allocated(qp%Z)) call error('[Chi/G] missing COHSEX QP energies/residues')
 tolerance=max(1.E-6_DP,100._DP*real(epsilon(1._SP),DP))
 allocate(covered(X%ib(1):X%ib(2),Xk%nibz,n_sp_pol))
 allocate(Chi_G_static_energy(X%ib(2),Xk%nibz,n_sp_pol))
 allocate(Chi_G_static_occupation(X%ib(2),Xk%nibz,n_sp_pol))
 allocate(Chi_G_norm(X%ib(2),Xk%nibz,n_sp_pol))
 covered=.FALSE.
 Chi_G_static_energy=0._DP; Chi_G_static_occupation=0._DP; Chi_G_norm=0._SP
 Chi_G_mu=Chi_KS_levels%E_Fermi
 do i=1,qp%n_states
   a=qp%table(i,1); ik=qp%table(i,3); s=spin(qp%table(i,:))
   if (qp%table(i,2)/=a) call error('[Chi/G] off-diagonal COHSEX states are unsupported')
   if (a<X%ib(1).or.a>X%ib(2)) cycle
   if (ik<1.or.ik>Xk%nibz.or.s<1.or.s>n_sp_pol) call error('[Chi/G] invalid COHSEX state index')
   if (covered(a,ik,s)) call error('[Chi/G] duplicate COHSEX state')
   if (.not.ieee_is_finite(real(qp%E_bare(i),DP)).or.&
&      abs(real(qp%E_bare(i),DP)-real(Chi_KS_levels%E(a,ik,s),DP))>1.E-5_DP)&
&     call error('[Chi/G] COHSEX starting energies differ from the saved KS levels')
   pole=real(qp%E(i),DP)
   if (.not.ieee_is_finite(pole).or..not.ieee_is_finite(real(aimag(qp%E(i)),DP)))&
&     call error('[Chi/G] nonfinite COHSEX energy')
   if (abs(real(aimag(qp%E(i)),DP))>tolerance) call error('[Chi/G] static COHSEX cannot have a QP linewidth')
   if (.not.ieee_is_finite(real(qp%Z(i),DP)).or..not.ieee_is_finite(real(aimag(qp%Z(i)),DP)))&
&     call error('[Chi/G] nonfinite COHSEX residue')
   if (abs(cmplx(qp%Z(i),kind=DP)-1._DP)>tolerance)&
&     call error('[Chi/G] static COHSEX requires unit residues')
   occupation=Chi_fermi(pole,real(Chi_G_mu,DP),0._DP)
   if (abs(occupation-real(Chi_KS_levels%f(a,ik,s),DP)/real(spin_occ,DP))>tolerance)&
&     call error('[Chi/G] COHSEX poles cross the fixed KS chemical potential; occupations require a new treatment')
   Chi_G_static_energy(a,ik,s)=pole
   Chi_G_static_occupation(a,ik,s)=occupation
   Chi_G_norm(a,ik,s)=1._SP
   covered(a,ik,s)=.TRUE.
 enddo
 if (.not.all(covered)) call error('[Chi/G] COHSEX QPkrange must cover every response band, k point and spin')
 deallocate(covered)
end subroutine
!
subroutine Chi_G_load(X,Xk)
 ! Load static COHSEX poles or retarded Dyson spectra without QP fitting/remapping.
 use pars, ONLY:SP,DP
 use Chi_m, ONLY:Chi_KS_levels,Chi_G_mode,Chi_G_db,Chi_G_energy,Chi_G_weight,Chi_G_occupation,&
&                 Chi_G_norm,Chi_G_norm_tol,Chi_G_mu,Chi_G_kbt
 use Chi_m, ONLY:Chi_G_static_energy
 use Chi_Green_m, ONLY:Chi_spectral_weights,Chi_fermi,Chi_Dyson
 use X_m, ONLY:X_t
 use R_lattice, ONLY:bz_samp
 use electrons, ONLY:spin,n_sp_pol
 use QP_m, ONLY:QP_t,QP_reset,QP_G_er,QP_G_dr,GF_is_causal,QP_G_Zoom_treshold
 use global_XC, ONLY:QP_DB_kind,QP_SE_COHSEX
 use IO_int, ONLY:io_control
 use IO_m, ONLY:OP_RD_CL,DUMP,REP
 use com, ONLY:msg
 use stderr, ONLY:intc
 implicit none
 type(X_t), intent(in) :: X
 type(bz_samp), intent(in) :: Xk
 type(QP_t) :: qp
 real(SP) :: saved_er(2),saved_dr(2),saved_zoom
 integer :: saved_kind,disk_kind,ID,io_err,i,a,ik,s,iw,n,info
 logical :: saved_causal
 logical, allocatable :: covered(:,:,:)
 real(DP) :: norm,reference_tolerance
 complex(DP), allocatable :: reconstructed(:)
 integer, external :: io_QP_and_GF
 if (allocated(Chi_G_energy).or.allocated(Chi_G_static_energy)) return
 if (len_trim(Chi_G_db)==0) call error('[Chi/G] ChiGDb must name a COHSEX ndb.QP or a retarded ndb.G')
 ! The generic reader modifies GW globals; preserve them in this response run.
 saved_er=QP_G_er; saved_dr=QP_G_dr; saved_zoom=QP_G_Zoom_treshold
 saved_causal=GF_is_causal; saved_kind=QP_DB_kind
 call QP_reset(qp)
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2,3/),MODE=DUMP,ID=ID)
 io_err=io_QP_and_GF(trim(Chi_G_db),qp,ID)
 disk_kind=QP_DB_kind
 QP_G_er=saved_er; QP_G_dr=saved_dr; QP_G_Zoom_treshold=saved_zoom
 GF_is_causal=saved_causal; QP_DB_kind=saved_kind
 if (io_err/=0) call error('[Chi/G] cannot read '//trim(Chi_G_db))
 if (qp%nk/=Xk%nibz) call error('[Chi/G] ChiGDb and response must have identical k meshes')
 if (qp%nb<X%ib(2)) call error('[Chi/G] ChiGDb does not cover the response band range')
 if (maxval(abs(qp%k-Xk%pt))>1.E-5_SP) call error('[Chi/G] ChiGDb k-point ordering differs from X')
 if (.not.allocated(qp%E_bare)) call error('[Chi/G] ChiGDb lacks KS reference energies')
 if (trim(Chi_G_mode)=='COHSEX') then
   if (disk_kind/=QP_SE_COHSEX) call error('[Chi/G] COHSEX mode requires an actual COHSEX ndb.QP')
   call Chi_G_load_static(qp,X,Xk)
   call QP_reset(qp)
   return
 endif
 if (.not.qp%GreenF_retarded) call error('[Chi/G] rebuild ndb.G with GRetarded; legacy spectra are unsupported')
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
