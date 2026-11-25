!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
module Chi_Green_m
 ! Independent numerical kernels for a diagonal, bare-vertex dressed bubble.
 ! Energies and frequencies are in Hartree. Spectral weights include dE.
 use pars, ONLY:DP
 use, intrinsic :: ieee_arithmetic, ONLY:ieee_is_finite
 implicit none
 private
 public :: Chi_G0_pair,Chi_spectral_weights,Chi_spectral_pair,Chi_fermi,Chi_Dyson
 public :: Chi_pole_bins,Chi_bins_setup,Chi_bins_size,Chi_bins_add,Chi_bins_eval
 !
 ! Binned sum of simple poles, S(w) = sum_p r_p/(w - x_p), with real pole energies x_p,
 ! evaluated at complex frequencies w with Im w >= eta > 0 (the retarded PPA self-energy
 ! on a Green-function grid). Each residue is split linearly between the two nearest
 ! nodes of an energy mesh, then the mesh is evaluated once at every w:
 !   fine nodes, step eta/bins_per_eta, over the real part of the grid plus a margin;
 !   coarse nodes, step margin/bins_per_margin, further out;
 !   a series in (w-c)/(x-c) about the grid centre c beyond the coarse mesh.
 ! Error of a pole relative to its own term: (step/|w-x|)^2/4 on the meshes,
 ! i.e. <= (1/bins_per_eta)^2/4 near the grid and (1/bins_per_margin)^2/4 beyond the
 ! margin; ratio^bins_n_mom (ratio <= 1/5) in the series.
 ! Packed storage b(:): fine nodes 0:n_f, coarse nodes 0:n_c, then bins_n_mom moments.
 !
 integer,  parameter :: bins_per_eta=80,bins_per_margin=50,bins_n_mom=8
 ! Negative spectral weight tolerated (and set to zero) in Chi_spectral_weights: dips down to
 ! Chi_clip_peak of the peak, adding up to at most Chi_clip_rel of the weight. The PPA self-energy
 ! (one pole per G,G' with weights not positive definite) can have Im Sigma > 0 at isolated
 ! frequencies (WS2 9x9: +17 eV at 61 eV from Eo, -1.1e-3 of the peak, 6e-4 of the weight);
 ! deeper or larger negative weight is still rejected as a noncausal spectrum.
 real(DP), parameter, public :: Chi_clip_rel=1.E-3_DP,Chi_clip_peak=1.E-2_DP
 real(DP), parameter :: bins_margin=10._DP/27.211386245988_DP   ! 10 eV, in Ha
 type Chi_pole_bins
   real(DP) :: x0_f=0._DP,d_f=1._DP,x0_c=0._DP,d_c=1._DP,c=0._DP
   integer  :: n_f=0,n_c=0
 end type
 contains
 !
 pure subroutine Chi_bins_setup(L,wr_min,wr_max,eta)
   type(Chi_pole_bins), intent(out) :: L
   real(DP), intent(in) :: wr_min,wr_max,eta
   real(DP) :: margin,extent
   margin=max(bins_margin,50._DP*eta)
   extent=max(40._DP*margin,4._DP*(wr_max-wr_min+2._DP*eta))
   L%c=(wr_min+wr_max)/2._DP
   L%d_f=eta/real(bins_per_eta,DP)
   L%x0_f=wr_min-margin
   L%n_f=ceiling((wr_max-wr_min+2._DP*margin)/L%d_f)
   L%d_c=margin/real(bins_per_margin,DP)
   L%x0_c=wr_min-extent
   L%n_c=ceiling((wr_max-wr_min+2._DP*extent)/L%d_c)
 end subroutine
 !
 pure function Chi_bins_size(L) result(n)
   type(Chi_pole_bins), intent(in) :: L
   integer :: n
   n=L%n_f+1+L%n_c+1+bins_n_mom
 end function
 !
 pure subroutine Chi_bins_add(L,b,x,r)
   type(Chi_pole_bins), intent(in) :: L
   complex(DP), intent(inout) :: b(:)
   real(DP), intent(in) :: x
   complex(DP), intent(in) :: r
   real(DP) :: s,t,y
   complex(DP) :: p
   integer :: j,k,o
   s=(x-L%x0_f)/L%d_f
   if (s>=0._DP.and.s<real(L%n_f,DP)) then
     j=int(s); t=s-real(j,DP)
     b(1+j)=b(1+j)+r*(1._DP-t)
     b(2+j)=b(2+j)+r*t
     return
   endif
   s=(x-L%x0_c)/L%d_c
   if (s>=0._DP.and.s<real(L%n_c,DP)) then
     o=L%n_f+1
     j=int(s); t=s-real(j,DP)
     b(o+1+j)=b(o+1+j)+r*(1._DP-t)
     b(o+2+j)=b(o+2+j)+r*t
     return
   endif
   ! Beyond the coarse mesh (also infinite or NaN x): moments sum r/(x-c)^(k+1)
   o=L%n_f+1+L%n_c+1
   y=1._DP/(x-L%c)
   p=r*y
   do k=1,bins_n_mom
     b(o+k)=b(o+k)+p
     p=p*y
   enddo
 end subroutine
 !
 pure subroutine Chi_bins_eval(L,b,w,s)
   type(Chi_pole_bins), intent(in) :: L
   complex(DP), intent(in) :: b(:),w(:)
   complex(DP), intent(out) :: s(:)
   complex(DP) :: u,acc
   real(DP) :: e
   integer :: iw,j,k,o
   s=(0._DP,0._DP)
   do j=0,L%n_f
     if (b(1+j)==(0._DP,0._DP)) cycle
     e=L%x0_f+real(j,DP)*L%d_f
     s=s+b(1+j)/(w-e)
   enddo
   o=L%n_f+1
   do j=0,L%n_c
     if (b(o+1+j)==(0._DP,0._DP)) cycle
     e=L%x0_c+real(j,DP)*L%d_c
     s=s+b(o+1+j)/(w-e)
   enddo
   ! 1/(w-x) = -sum_k (w-c)^k/(x-c)^(k+1)
   o=L%n_f+1+L%n_c+1
   do iw=1,size(w)
     u=w(iw)-L%c
     acc=(0._DP,0._DP)
     do k=bins_n_mom,1,-1
       acc=acc*u+b(o+k)
     enddo
     s(iw)=s(iw)-acc
   enddo
 end subroutine
 !
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
 subroutine Chi_spectral_weights(e,g,weight,norm,info,clipped)
   real(DP), intent(in) :: e(:)
   complex(DP), intent(in) :: g(:)
   real(DP), intent(out) :: weight(size(e)),norm
   integer, intent(out) :: info
   real(DP), intent(out), optional :: clipped
   real(DP) :: a(size(e)),de(size(e)),negative_tolerance,neg_weight
   integer :: n
   info=1
   weight=0._DP
   norm=0._DP
   if (present(clipped)) clipped=0._DP
   n=size(e)
   if (n<2.or.size(g)/=n) return
   if (.not.all(ieee_is_finite(e))) return
   if (.not.all(ieee_is_finite(real(g,DP)))) return
   if (.not.all(ieee_is_finite(aimag(g)))) return
   if (any(e(2:)<=e(:n-1))) return
   a=-aimag(g)/acos(-1._DP)
   de(1)=(e(2)-e(1))/2._DP
   de(n)=(e(n)-e(n-1))/2._DP
   if (n>2) de(2:n-1)=(e(3:)-e(:n-2))/2._DP
   ! Reject a noncausal spectrum: negatives beyond round-off and beyond Chi_clip_peak
   ! of the peak, or adding up to more than Chi_clip_rel. Tolerated ones are set to zero.
   negative_tolerance=max(100._DP*epsilon(1._DP)*max(1._DP,maxval(abs(a))),&
&                         Chi_clip_peak*max(maxval(a),0._DP))
   info=2
   if (minval(a)<-negative_tolerance) return
   neg_weight=-sum(min(a,0._DP)*de)
   if (neg_weight>Chi_clip_rel*max(sum(max(a,0._DP)*de),0._DP)) return
   if (present(clipped)) clipped=neg_weight
   a=max(a,0._DP)
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
