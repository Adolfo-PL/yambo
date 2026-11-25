!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
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
 ! small negative dips (PPA Im Sigma > 0) are clipped and reported, deeper ones rejected
 call Chi_spectral_weights(e,g,weight,norm,info)
 g(7000)=cmplx(real(g(7000)),5.E-4_DP*maxval(-aimag(g)),DP)
 call Chi_spectral_weights(e,g,weight,norm,info,omega)
 call require(info==0.and.omega>0._DP.and.weight(7000)==0._DP,'small negative dip clipped')
 g(7000)=cmplx(real(g(7000)),2.E-2_DP*maxval(-aimag(g)),DP)
 call Chi_spectral_weights(e,g,weight,norm,info)
 call require(info==2,'deep negative dip rejected')
 g(6900:7100)=cmplx(real(g(6900:7100)),9.E-4_DP*maxval(-aimag(g)),DP)
 call Chi_spectral_weights(e,g,weight,norm,info)
 call require(info==2,'wide shallow negative region rejected (integrated weight)')
 e(2)=e(1)
 call Chi_spectral_weights(e,g,weight,norm,info)
 call require(info==1,'nonmonotonic spectral grid rejected')
 print *, 'PASS: numerical G0, Dyson, spectral, Matsubara and causality checks'
 call test_pole_bins()
contains
 subroutine test_pole_bins()
   ! Binned pole sum against the direct sum on the production grid shape:
   ! Eo +- 40 eV, 589 points, eta = 0.136 eV (Hartree units).
   integer, parameter :: nw=589,np=20000
   real(DP), parameter :: ha=27.211386245988_DP
   type(Chi_pole_bins) :: L
   complex(DP), allocatable :: b(:)
   complex(DP) :: w(nw),sd(nw),sb(nw),r(np)
   real(DP) :: x(np),u(4),eo,eta,scale(nw),err,err_rel
   integer :: i,ip,seed(64)
   seed=12345
   call random_seed(put=seed(1:size(seed)))
   eo=0.1_DP; eta=0.136_DP/ha
   do i=1,nw
     w(i)=cmplx(eo+(-40._DP+80._DP*real(i-1,DP)/real(nw-1,DP))/ha,eta,DP)
   enddo
   do ip=1,np
     call random_number(u)
     if (ip<=12000) then                      ! across the grid and its margin
       x(ip)=eo+(-60._DP+120._DP*u(1))/ha
     else if (ip<=16000) then                 ! clusters within eta of grid points
       x(ip)=real(w(1+int(u(1)*(nw-1))),DP)+(u(2)-0.5_DP)*eta
     else if (ip<=19000) then                 ! coarse mesh
       x(ip)=eo+sign(1._DP,u(1)-0.5_DP)*(50._DP+(1000._DP*u(2)))/ha
     else                                     ! far: series
       x(ip)=eo+sign(1._DP,u(1)-0.5_DP)*(1.E3_DP+1.E6_DP*u(2))/ha
     endif
     r(ip)=cmplx(u(3)-0.5_DP,u(4)-0.5_DP,DP)
   enddo
   ! strong cancellation: pairs of opposite residues 1 meV apart
   do ip=1,2000,2
     x(ip+1)=x(ip)+0.001_DP/ha; r(ip+1)=-r(ip)
   enddo
   call Chi_bins_setup(L,real(w(1),DP),real(w(nw),DP),eta)
   allocate(b(Chi_bins_size(L)))
   b=(0._DP,0._DP)
   do ip=1,np
     call Chi_bins_add(L,b,x(ip),r(ip))
   enddo
   call Chi_bins_eval(L,b,w,sb)
   sd=(0._DP,0._DP); scale=0._DP
   do ip=1,np
     sd=sd+r(ip)/(w-x(ip))
     scale=scale+abs(r(ip)/(w-x(ip)))
   enddo
   err=maxval(abs(sb-sd)/scale)
   err_rel=maxval(abs(sb-sd))/maxval(abs(sd))
   write(*,'(a,es10.2,a,es10.2,a,i0)') ' pole bins: max |error|/sum|terms| ',err,&
&    ', max |error|/max|S| ',err_rel,', bins ',size(b)
   call require(err<2.E-4_DP,'binned pole sum (term scale)')
   call require(err_rel<1.E-3_DP,'binned pole sum (sum scale)')
   ! infinite and NaN-free far poles add nothing
   b=(0._DP,0._DP)
   call Chi_bins_add(L,b,huge(1._DP),(1._DP,0._DP))
   call Chi_bins_eval(L,b,w,sb)
   call require(maxval(abs(sb))<1.E-30_DP,'pole at infinity')
   deallocate(b)
 end subroutine
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
