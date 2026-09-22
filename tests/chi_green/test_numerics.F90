program test_chi_green
 use pars, only:DP
 use Chi_Green_m
 implicit none
 complex(DP) :: z(3),p(3),expected(3),g(20001),sum_m,gp,gm
 real(DP) :: e(20001),weight(20001),norm,omega,beta,fa,fb
 real(DP) :: ea(2),eb(2),wa(2),wb(2),f_a(2),f_b(2)
 integer :: i,m,info
 z=[cmplx(0._DP,0._DP,DP),cmplx(0._DP,0.6_DP,DP),cmplx(1.2_DP,0.2_DP,DP)]
 ea=[1._DP,3._DP]; eb=[-1._DP,-4._DP]
 wa=[0.7_DP,0.3_DP]; wb=[0.8_DP,0.2_DP]
 f_a=0._DP; f_b=1._DP
 call Chi_spectral_pair(z,ea,eb,wa,wb,f_a,f_b,p,info)
 call require(info==0,'two-pole convolution status')
 expected=0.56_DP/(z-2._DP)+0.14_DP/(z-5._DP)+&
          0.24_DP/(z-4._DP)+0.06_DP/(z-7._DP)
 call close(p,expected,1.E-13_DP,'satellites and independent band indices')
 call Chi_spectral_pair(z,eb,ea,wb,wa,f_b,f_a,p,info)
 expected=-0.56_DP/(z+2._DP)-0.14_DP/(z+5._DP)-&
          0.24_DP/(z+4._DP)-0.06_DP/(z+7._DP)
 call close(p,expected,1.E-13_DP,'antiresonant poles')
 wa=[1._DP,0._DP]; wb=[1._DP,0._DP]
 call Chi_spectral_pair(z,ea,eb,wa,wb,f_a,f_b,p,info)
 do i=1,3
   expected(i)=Chi_G0_pair(z(i),1._DP,-1._DP,0._DP,1._DP)
 enddo
 call close(p,expected,1.E-13_DP,'G0 delta-spectrum limit')
 call Chi_spectral_pair(z,ea,ea,wa,wa,f_a,f_a,p,info)
 call require(info==0.and.maxval(abs(p))==0._DP,'equal-energy static 0/0 skipped')
 ! Independent explicit Matsubara sum fixes the susceptibility sign.
 beta=8._DP
 fa=Chi_fermi(1._DP,0._DP,1._DP/beta)
 fb=Chi_fermi(-1._DP,0._DP,1._DP/beta)
 sum_m=0._DP
 do m=-100000,99999
   omega=(2*m+1)*acos(-1._DP)/beta
   gp=1._DP/cmplx(-1._DP,omega+2*acos(-1._DP)/beta,DP)
   gm=1._DP/cmplx(1._DP,omega,DP)
   sum_m=sum_m+gp*gm/beta
 enddo
 expected(1)=Chi_G0_pair(cmplx(0._DP,2*acos(-1._DP)/beta,DP),1._DP,-1._DP,fa,fb)
 call require(abs(sum_m-expected(1))<5.E-6_DP,'explicit fermionic Matsubara sum')
 call require(Chi_fermi(1000._DP,0._DP,0.01_DP)==0._DP,'stable Fermi positive tail')
 call require(Chi_fermi(-1000._DP,0._DP,0.01_DP)==1._DP,'stable Fermi negative tail')
 ! Retarded Dyson spectral integral and high-frequency asymptote.
 do i=1,size(e)
   e(i)=-100._DP+200._DP*real(i-1,DP)/real(size(e)-1,DP)
   g(i)=Chi_Dyson(cmplx(e(i),0.1_DP,DP),0.3_DP,cmplx(0.2_DP,-0.1_DP,DP))
 enddo
 call Chi_spectral_weights(e,g,weight,norm,info)
 call require(info==0.and.abs(norm-1._DP)<0.002_DP,'causal Dyson spectral normalization')
 gp=Chi_Dyson(cmplx(0._DP,1.E6_DP,DP),0.3_DP,cmplx(0.2_DP,-0.1_DP,DP))
 call require(abs(cmplx(0._DP,1.E6_DP,DP)*gp-1._DP)<1.E-6_DP,'G high-frequency normalization')
 call Chi_spectral_weights(e,conjg(g),weight,norm,info)
 call require(info==2,'advanced/noncausal spectrum rejected')
 e(2)=e(1)
 call Chi_spectral_weights(e,g,weight,norm,info)
 call require(info==1,'nonmonotonic spectral grid rejected')
 print *, 'PASS: numerical G0, Dyson, spectral, Matsubara and causality checks'
contains
 subroutine require(ok,label)
   logical,intent(in)::ok
   character(*),intent(in)::label
   if (.not.ok) then
     print *, 'FAIL: ',label
     error stop 1
   endif
 end subroutine
 subroutine close(actual,reference,tolerance,label)
   complex(DP),intent(in)::actual(:),reference(:)
   real(DP),intent(in)::tolerance
   character(*),intent(in)::label
   call require(maxval(abs(actual-reference))<tolerance,label)
 end subroutine
end program
