!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine X_pre_setup(Xen,Xk,X)
 !
 use pars,          ONLY:SP
 use X_m,           ONLY:X_t,use_X_DbGd,X_DbGd_nkpts,l_X_terminator,X_use_lin_sys,&
&                        X_use_gpu,Chi_linalg_mode,X_ng,l_write_disk_Xo,l_recalculate_X
 use Chi_m,         ONLY:l_Chi_fxc,Chi_G_mode,Chi_KS_levels,Chi_G_norm_tol,Chi_G_check_tol,&
&                        Chi_rcond_floor,Chi_G_kbt,Chi_G_mu
 use Chi_m,         ONLY:Chi_G_axis,Chi_G_nu_steps,Chi_G_nu_range
 use electrons,     ONLY:levels,n_sp_pol
 use drivers,       ONLY:l_bs_fxc,l_col_cut,Finite_Tel,l_cohsex,l_rt_carriers_in_use,l_rpa_IP
 use R_lattice,     ONLY:nqibz,bz_samp
 use QP_m,          ONLY:COHSEX_use_empties
 use TDDFT,         ONLY:FXC_n_g_corr
 use interfaces,    ONLY:eval_G_minus_G
 use parallel_m,    ONLY:PAR_COM_WORLD
 use stderr,        ONLY:intc
 use parser_m,      ONLY:parser
 use gpu_m,         ONLY:have_gpu
 !
 implicit none
 !
 type(X_t)    :: X 
 type(bz_samp):: Xk
 type(levels) :: Xen
 ! 
 ! Work Space
 !
 real(SP), external :: G2E
 !
 if (X%ib(1)<0) X%ib(1)=1
 if (X%ib(2)<0) X%ib(2)=1
 if (X%iq(2)<0.or.X%iq(2)>nqibz) X%iq(2)=nqibz
 if (X%iq(1)<0) X%iq(1)=1
 if (X%iq(1)>X%iq(2)) X%iq(1)=X%iq(2)
 !
 ! In these two lines I ensure the at least all the metallic
 ! bands are included in the X summation
 !
 if (X%ib(1)>minval(Xen%nbf)+1) X%ib(1)=minval(Xen%nbf)+1
 if (X%ib(2)<maxval(Xen%nbm)  ) X%ib(2)=maxval(Xen%nbm)
 !
 ! Other checks
 !
 if (X%ib(1)<0.or.X%ib(1)>maxval(Xen%nbm).or.X%ib(2)<X%ib(1))  X%ib(1)=1
 if (X%ib(2)>Xen%nb.or.X%ib(2)<minval((/Xen%nbf(:n_sp_pol)/))) X%ib(2)=Xen%nb
 !
 ! Here I set up the Double Grid flag to (for example) skip the coarse 
 ! grid generator in the Xo calculation
 !
 use_X_DbGd=Xen%FG%nb/=0.and.X%ib(1)<Xen%FG%nb.and.(.not.l_bs_fxc)
 !
 if (use_X_DbGd) X_DbGd_nkpts=Xk%nbz+Xk%FGbz%N
 !
 ! COH check: when the X terminator is used or the empty bands are NOT used the
 ! response function is accessed in G-G' mode. Here I ensure that all G-G' are in the
 ! range provided externally
 !
 if ((l_cohsex.and..not.COHSEX_use_empties).or.l_X_terminator) then
     X%ng=eval_G_minus_G(X%ng,0,COMM=PAR_COM_WORLD)
     call warning(' Response block size reduced to '//intc(X%ng) &
&     //'RL ('//intc(int(1000.*G2E(X%ng)))//'mHa)')
 endif
 !
 ! Fxc Checks
 !
 if (FXC_n_g_corr>X%ng) FXC_n_g_corr=X%ng
 !
 ! In case the columbian is cutoffed a resonant/antiresonant 
 ! only calculation is meaningless because of the dependence 
 ! of the alpha from X
 !
 if (l_col_cut) then
   if ( trim(X%ordering) == 'a' .or. trim(X%ordering) == 'r' )  then
      X%ordering='T'
      call warning(" Anti/Resonant ordering not compatible with Coulomb cutoff, set to 'T' ")
   endif
 endif
 !
 ! Finite temperature X MUST be causal. This is because in X_eh_setup the
 ! "resonant" or "antiresonant" par of X (in the many body language) include
 ! negative poles as well. Only the sum of the two is meaningful, then.
 !
 if (Finite_Tel .or. l_rt_carriers_in_use) X%ordering="R"
 !
 ! Transfer the GF's ordering 
 !
 if ( allocated(Xen%GreenF).and.     Xen%GreenF_is_causal ) X%ordering='R'
 if ( allocated(Xen%GreenF).and..not.Xen%GreenF_is_causal ) then
   if ( .not.(Finite_Tel.or.l_rt_carriers_in_use) ) X%ordering='T'
   if (      (Finite_Tel.or.l_rt_carriers_in_use) ) then
     call error('Fractional occupations but GF is not retarded (causal in yambo language)')
   endif 
 endif
 !
 ! Dyson solver
 !
 X_use_lin_sys=.true.
 X_use_gpu=have_gpu
 if (index(Chi_linalg_mode,"INVERSION")>0) X_use_lin_sys=.false.
 if (index(Chi_linalg_mode,"LIN_SYS")>0  ) X_use_lin_sys=.true.
 if (index(Chi_linalg_mode,"CPU")>0      ) X_use_gpu=.false.
 if (index(Chi_linalg_mode,"GPU")>0      ) X_use_gpu=.true..and.have_gpu
 !
 ! Copy of X dimension for Parallel default environment solver
 !
 X_ng=X%ng
 !
 ! Check if Xo must be written on-the-fly
 !
 call parser('WriteXo',l_write_disk_Xo)
 !
 ! ndb.Chi needs BOTH bubbles, and they are built inside the q-loop of
 ! X_dielectric_matrix. A complete response database on disk would make that loop
 ! return early and no ndb.Chi would be written, so force the recalculation.
 !
 call parser('ChiFxc',l_Chi_fxc)
 if (l_Chi_fxc) then
   if (trim(Chi_G_mode)/='LEVELS'.and.trim(Chi_G_mode)/='G0'.and.&
&      trim(Chi_G_mode)/='DYSON'.and.trim(Chi_G_mode)/='COHSEX')&
&     call error('[Chi] ChiGMode must be LEVELS, G0, COHSEX or DYSON')
   if (Chi_G_norm_tol<=0._SP.or.Chi_G_norm_tol>=1._SP.or.Chi_G_check_tol<=0._SP.or.&
&      Chi_rcond_floor<0._SP.or.Chi_rcond_floor>=1._SP) call error('[Chi] invalid validation tolerances')
   if (trim(Chi_G_mode)/='LEVELS') then
     if (trim(Chi_G_axis)/='REAL'.and.trim(Chi_G_axis)/='IMAG') call error('[Chi/G] ChiGAxis must be REAL or IMAG')
     if (Chi_G_nu_steps<1.or.Chi_G_nu_range(1)<0._SP.or.Chi_G_nu_range(2)<Chi_G_nu_range(1))&
&       call error('[Chi/G] invalid imaginary frequency range/steps')
     if (have_gpu) call error('[Chi/G] CPU execution required; GPU bubble kernels are not implemented')
     if (X%iq(1)<=1) call error('[Chi/G] select finite q only: QpntsRXs or QpntsRXd must start after q=1')
     if (l_X_terminator.or.use_X_DbGd) call error('[Chi/G] disable XTermKind and double-grid interpolation')
     if (X%ehe(1)>0._SP.or.X%ehe(2)>0._SP)&
&       call error('[Chi/G] transition-energy filtering is unsupported; set EhEngyXs or EhEngyXd to zero')
     if (any(Chi_KS_levels%nbf/=Chi_KS_levels%nbm))&
&       call error('[Chi/G] first implementation requires an insulating KS reference')
     Chi_G_kbt=0._SP
     Chi_G_mu=Chi_KS_levels%E_Fermi
     if (Finite_Tel.or.l_rt_carriers_in_use) call error('[Chi/G] finite temperature/nonequilibrium is unsupported')
     if (allocated(Xen%W).or.allocated(Xen%Z).or.allocated(Xen%GreenF))&
&       call error('[Chi/G] use ChiGDb alone; disable external X widths, residues and Green functions')
     if (maxval(abs(Xen%E-Chi_KS_levels%E))>1.E-5_SP)&
&       call error('[Chi/G] disable external X energy corrections; G0 must use the KS reference')
     X%ordering='R'
     X%cg_percentual=0._SP
     if (trim(Chi_G_mode)=='DYSON'.or.trim(Chi_G_mode)=='COHSEX') call Chi_G_load(X,Xk)
   endif
 endif
 if (l_Chi_fxc) then
   l_recalculate_X=.TRUE.
   call warning('[Chi] the response matrix is recomputed: ndb.Chi needs chi0 and P together')
 endif
 !
 if (l_rpa_IP) l_write_disk_Xo=.FALSE.
 !
end subroutine
