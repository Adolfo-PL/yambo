program test_pipeline
 use pars
 use Chi_m
 use X_m
 use R_lattice
 use frequency
 use DIPOLES
 use pipeline_output
 use QP_m,only:QP_t,QP_Sc_steps,QP_G_amplitude_integral,QP_reset,QP_table,QP_n_states,&
&              QP_Sc,QP_Vnl_xc,QP_Vxc,use_GreenF_Zoom,use_GreenF_to_eval_QP
 implicit none
 type(X_t)::x
 type(bz_samp)::k
 type(w_samp)::w
 type(DIPOLE_t)::d
 type(QP_t)::qp
 complex(SP)::original(2,2,2),factor0,factor1,expected(2,2),gram_inverse(2,2)
 integer::iw
 character(24)::argument
 call get_command_argument(1,argument)
 allocate(k%sstar(2,2),k%pt(2,3),qindx_X(2,2,2),bare_qpg(2,2))
 k%sstar(:,1)=[1,2];k%sstar(:,2)=1;k%pt=0._SP
 qindx_X=1;qindx_X(2,:,1)=[2,1];bare_qpg=1._SP
 allocate(Chi_KS_levels%E(2,2,1),Chi_KS_levels%f(2,2,1))
 Chi_KS_levels%E(1,:,:)= -1._SP;Chi_KS_levels%E(2,:,:)=1._SP
 Chi_KS_levels%f(1,:,:)=2._SP;Chi_KS_levels%f(2,:,:)=0._SP
 allocate(w%p(2),X_par(1))
 allocate(X_par(1)%blc(2,2,2))
 w%p=[cmplx(.1_SP,.2_SP,SP),cmplx(.8_SP,.2_SP,SP)]
 original=cmplx(.3_SP,.4_SP,SP)
 X_par(1)%blc=original
 Chi_G_mode='G0';Chi_G_axis='IMAG';Chi_G_nu_steps=3;Chi_G_nu_range=[0._SP,2._SP]
 if(trim(argument)=='conditioning') Chi_rcond_floor=0.99_SP
 call Chi_fxc_eval(2,1,x,Chi_KS_levels,k,w,d)
 call require(size(saved_freq)==3,'separate imaginary export grid')
 call require(maxval(abs(real(saved_freq)))<1.E-12_SP,'pure imaginary frequencies')
 call require(maxval(abs(aimag(saved_freq)-[0._SP,1._SP,2._SP]))<1.E-12_SP,'frequency endpoints')
 call require(maxval(abs(saved_fxc))<1.E-12_SP,'G0 kernel vanishes')
 call require(maxval(abs(X_par(1)%blc-original))<1.E-12_SP,'screening matrix restored')
 call require(size(w%p)==2.and.abs(w%p(2)-cmplx(.8_SP,.2_SP,SP))<1.E-12_SP,'native grid preserved')
 call require(minval(saved_rcond)>0._SP,'LAPACK conditioning diagnostics')
 Chi_G_mode='DYSON'
 allocate(Chi_G_energy(2,2,1,1),Chi_G_weight(2,2,1,1),Chi_G_occupation(2,2,1,1))
 Chi_G_energy(1,:,:,:)=-1.2_DP;Chi_G_energy(2,:,:,:)=1.4_DP
 Chi_G_weight=1._DP
 Chi_G_occupation(1,:,:,:)=1._DP;Chi_G_occupation(2,:,:,:)=0._DP
 X_par(1)%blc=original
 call Chi_fxc_eval(2,1,x,Chi_KS_levels,k,w,d)
 gram_inverse=reshape([1.04_SP,-.4_SP,-.4_SP,1.04_SP],[2,2])/(1.04_SP**2-.4_SP**2)
 do iw=1,3
   factor0=1._SP/(saved_freq(iw)-2._SP)-1._SP/(saved_freq(iw)+2._SP)
   factor1=1._SP/(saved_freq(iw)-2.6_SP)-1._SP/(saved_freq(iw)+2.6_SP)
   expected=(1._SP/factor0-1._SP/factor1)*gram_inverse
   call require(maxval(abs(saved_fxc(:,:,iw)-expected))<max(1.E-11_SP,100._SP*epsilon(1._SP)),&
&               'Dyson bubble inverse difference')
 enddo
 call require(maxval(abs(X_par(1)%blc-original))<1.E-12_SP,'Dyson export preserves screening')
 call Chi_G_free()
 Chi_G_db='fixture.ndb.G'
 if(trim(argument)=='missing_state')Chi_G_db='missing.ndb.G'
 if(trim(argument)=='legacy')Chi_G_db='legacy.ndb.G'
 if(trim(argument)=='reference')Chi_G_db='reference.ndb.G'
 call Chi_G_load(x,k)
 call require(minval(Chi_G_norm)>.98_SP,'loader retains unnormalized spectral weight')
 call require(maxval(Chi_G_norm)<1._SP,'loader does not hide finite spectral tails')
 call Chi_G_free()
 ! Exercise production QP_Green_Function, including the XC subtraction.
 QP_Sc_steps=3;QP_n_states=1;use_GreenF_Zoom=.FALSE.;use_GreenF_to_eval_QP=.FALSE.
 allocate(QP_G_amplitude_integral(1),QP_table(1,3),QP_Sc(1,3),QP_Vnl_xc(1),QP_Vxc(1))
 QP_table(1,:)=[2,2,1];QP_Vnl_xc=.2_SP;QP_Vxc=.1_SP
 QP_Sc=cmplx(.3_SP,-.1_SP,SP)
 qp%GreenF_retarded=.TRUE.
 allocate(qp%GreenF(1,3),qp%GreenF_W(1,3),qp%S_total(1,3))
 call QP_Green_Function(qp,Chi_KS_levels,-1)
 do iw=1,3
   factor1=1._SP/(qp%GreenF_W(1,iw)-1._SP-cmplx(.4_SP,-.1_SP,SP))
   call require(abs(qp%GreenF(1,iw)-factor1)<100._SP*epsilon(1._SP),'retarded Dyson and XC subtraction')
 enddo
 call require(maxval(aimag(qp%GreenF))<0._SP,'retarded Green spectral sign')
 call QP_reset(qp)
 deallocate(QP_G_amplitude_integral,QP_table,QP_Sc,QP_Vnl_xc,QP_Vxc)
 print *,'PASS: native bubble/inversion pipeline with imaginary grid and frozen screening'
contains
 subroutine require(ok,label)
   logical::ok
   character(*)::label
   if(.not.ok)then
     print *,'FAIL: ',label
     error stop 1
   endif
 end subroutine
end program
