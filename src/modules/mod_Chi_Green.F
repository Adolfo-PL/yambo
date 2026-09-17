!
! License-Identifier: GPL
!
module Chi_Green_m
 ! Independent numerical kernels for a diagonal, bare-vertex dressed bubble.
 ! Energies and frequencies are in Hartree. Spectral weights include dE.
 use pars, ONLY:DP
 use, intrinsic :: ieee_arithmetic, ONLY:ieee_is_finite
 implicit none
 private
 public :: Chi_G0_pair,Chi_spectral_weights,Chi_spectral_pair,Chi_fermi,Chi_Dyson
 contains
 !
 pure function Chi_Dyson(z,e,delta_sigma) result(g)
   complex(DP), intent(in) :: z,delta_sigma
   real(DP), intent(in) :: e
   complex(DP) :: g
   g=1._DP/(z-e-delta_sigma)
 end function
 !
 pure function Chi_fermi(e,mu,kbt) result(f)
   real(DP), intent(in) :: e,mu,kbt
   real(DP) :: f,x
   if (kbt>0._DP) then
     x=(e-mu)/kbt
     ! Stable also for large positive/negative arguments.
     f=exp(-max(x,0._DP))/(1._DP+exp(-abs(x)))
   else
     f=0.5_DP
     if (e<mu) f=1._DP
     if (e>mu) f=0._DP
   endif
 end function
 !
 pure function Chi_G0_pair(z,ea,eb,fa,fb) result(p)
   complex(DP), intent(in) :: z
   real(DP), intent(in) :: ea,eb,fa,fb
   complex(DP) :: p
   p=(0._DP,0._DP)
   if (fa==fb) return
   p=(fb-fa)/(z-ea+eb)
 end function
 !
 subroutine Chi_spectral_weights(e,g,weight,norm,info)
   real(DP), intent(in) :: e(:)
   complex(DP), intent(in) :: g(:)
   real(DP), intent(out) :: weight(size(e)),norm
   integer, intent(out) :: info
   real(DP) :: a(size(e)),de(size(e)),negative_tolerance
   integer :: n
   info=1
   weight=0._DP
   norm=0._DP
   n=size(e)
   if (n<2.or.size(g)/=n) return
   if (.not.all(ieee_is_finite(e))) return
   if (.not.all(ieee_is_finite(real(g,DP)))) return
   if (.not.all(ieee_is_finite(aimag(g)))) return
   if (any(e(2:)<=e(:n-1))) return
   a=-aimag(g)/acos(-1._DP)
   ! Reject a noncausal spectrum. Only roundoff-sized negatives may be clipped.
   negative_tolerance=100._DP*epsilon(1._DP)*max(1._DP,maxval(abs(a)))
   info=2
   if (minval(a)<-negative_tolerance) return
   a=max(a,0._DP)
   de(1)=(e(2)-e(1))/2._DP
   de(n)=(e(n)-e(n-1))/2._DP
   if (n>2) de(2:n-1)=(e(3:)-e(:n-2))/2._DP
   weight=a*de
   norm=sum(weight)
   info=3
   if (.not.ieee_is_finite(norm).or.norm<=0._DP) return
   ! Never normalize away missing satellites or inadequate energy coverage.
   info=0
 end subroutine
 !
 subroutine Chi_spectral_pair(z,ea,eb,wa,wb,fa,fb,p,info)
   complex(DP), intent(in) :: z(:)
   real(DP), intent(in) :: ea(:),eb(:),wa(:),wb(:),fa(:),fb(:)
   complex(DP), intent(out) :: p(size(z))
   integer, intent(out) :: info
   complex(DP) :: denominator(size(z))
   real(DP) :: coefficient
   integer :: ia,ib
   p=(0._DP,0._DP)
   info=1
   if (size(wa)/=size(ea).or.size(fa)/=size(ea)) return
   if (size(wb)/=size(eb).or.size(fb)/=size(eb)) return
   if (any(aimag(z)<0._DP)) return
   if (any(wa<0._DP).or.any(wb<0._DP)) return
   if (any(fa<0._DP).or.any(fa>1._DP).or.any(fb<0._DP).or.any(fb>1._DP)) return
   do ib=1,size(eb)
     if (wb(ib)==0._DP) cycle
     do ia=1,size(ea)
       coefficient=wa(ia)*wb(ib)*(fb(ib)-fa(ia))
       ! Includes equal-energy/equal-occupation terms; avoid 0/0 at z=0.
       if (coefficient==0._DP) cycle
       denominator=z-ea(ia)+eb(ib)
       info=2
       if (any(abs(denominator)<=tiny(1._DP))) return
       p=p+coefficient/denominator
     enddo
   enddo
   info=0
 end subroutine
end module Chi_Green_m
