integer function io_Chi(x,w,id)
 use Chi_m
 use X_m,only:X_t
 use frequency
 use pipeline_output
 type(X_t)::x
 type(w_samp)::w
 integer::id
 saved_c0=Chi_0; saved_p=Chi_P; saved_fxc=Chi_fxc
 saved_freq=Chi_freqs; saved_rcond=Chi_rcond
 io_Chi=2
end function
subroutine FREQUENCIES_Green_Function(i,w,e,spectral)
 use pars
 use frequency
 integer::i
 type(w_samp)::w
 real(SP)::e(*)
 logical::spectral
 allocate(w%p(3))
 w%p=[cmplx(-2._SP,.2_SP,SP),cmplx(0._SP,.2_SP,SP),cmplx(2._SP,.2_SP,SP)]
end subroutine
subroutine QPartilize()
 error stop 'GreenF2QP must be disabled for the retarded test'
end subroutine
function RIntegrate(a,e,n) result(value)
 use pars
 integer::n,i
 real(SP)::a(n),e(n),value
 value=0._SP
 do i=2,n
   value=value+(a(i)+a(i-1))*(e(i)-e(i-1))/2._SP
 enddo
end function
integer function io_QP_and_GF(s,qp,id)
 use pars
 use QP_m
 character(*)::s
 type(QP_t)::qp
 integer::id,i,iw,band,kpt,n
 real(SP)::energy,delta
 qp%nk=2;qp%nb=2;qp%n_states=4;qp%GreenF_n_steps=2001;qp%GreenF_retarded=.TRUE.
 qp%GreenF_mu=0._SP
 if(index(s,'missing')>0)qp%n_states=3
 if(index(s,'legacy')>0)qp%GreenF_retarded=.FALSE.
 n=qp%GreenF_n_steps
 allocate(qp%table(qp%n_states,3),qp%k(2,3),qp%E_bare(qp%n_states))
 allocate(qp%GreenF(qp%n_states,n),qp%GreenF_W(qp%n_states,n),qp%S_total(qp%n_states,n))
 qp%k=0._SP
 do i=1,qp%n_states
   band=mod(i-1,2)+1;kpt=(i-1)/2+1
   qp%table(i,:)=[band,band,kpt]
   qp%E_bare(i)=real(2*band-3,SP)
   delta=-.2_SP
   if(band==2)delta=.4_SP
   do iw=1,n
     energy=-20._SP+.02_SP*real(iw-1,SP)
     qp%GreenF_W(i,iw)=cmplx(energy,.2_SP,SP)
     qp%S_total(i,iw)=delta
     qp%GreenF(i,iw)=1._SP/(qp%GreenF_W(i,iw)-qp%E_bare(i)-delta)
   enddo
 enddo
 if(index(s,'reference')>0)qp%E_bare(1)=qp%E_bare(1)+1._SP
 io_QP_and_GF=0
end function
