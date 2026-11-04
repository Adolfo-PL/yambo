!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Reading copy: the numerics only (no argument checks, error codes or condition numbers).
!
! LEGEND (names as in the full code)
!   z     complex frequency w + i eta
!   nh    number of rows of the coupled matrix, 2N: N resonant + N anti-resonant transitions
!   ng    number of plane waves G
!   v     (ng,nh)  vertices v(G,A) = rho_A(G), the bare transition densities
!   s     (nh)     s_A = sqrt(f_A); i sqrt|f_A| for the anti-resonant rows (f_A < 0)
!   c     spin_occ/(Omega Nk) = Co/4pi
!   E     (nh)     reference transition energies E_A (negative for the anti-resonant rows)
!   M,H0  (nh,nh)  BSE matrix in Yambo's storage; H0 = diag(E) + exchange, i.e. without the fxc part
!   d     (nh)     optical dipoles d_A         [full_dynamic_response]
!         (ng)     D_G = sqrt(4 pi)/|q+G|      [kernel_from_response]
!   chi0  bare bubble;  chib  BSE response (exchange included);  P  polarization (exchange removed)
!   fxc   exact kernel,  resp  head response
!
module TDDFT_Chi_dynamic_m
 use pars, only: DP
 implicit none
 private
 public :: TDDFT_Chi_full_response,TDDFT_Chi_full_bubble
 public :: TDDFT_Chi_kernel_from_response,TDDFT_Chi_full_dynamic_response
 complex(DP), parameter :: zero=(0._DP,0._DP),one=(1._DP,0._DP)
contains
 !
 recursive subroutine TDDFT_Chi_full_bubble(z,E,v,s,c,chi)
   ! chi0(G,G') = c sum_A s_A^2 v(G,A) conj(v(G',A)) / (z - E_A)
   complex(DP), intent(in)  :: z,v(:,:),s(:)
   real(DP),    intent(in)  :: E(:),c
   complex(DP), intent(out) :: chi(:,:)
   complex(DP), allocatable :: Y(:,:)
   integer :: nh,ng,i
   nh=size(s); ng=size(v,1)
   allocate(Y(ng,nh))
   do i=1,nh
     Y(:,i)=v(:,i)*s(i)**2/(z-E(i))
   enddo
   call zgemm('N','C',ng,ng,nh,cmplx(c,0._DP,kind=DP),Y,ng,v,ng,zero,chi,ng)
   deallocate(Y)
 end subroutine TDDFT_Chi_full_bubble
 !
 recursive subroutine TDDFT_Chi_full_response(z,M,v,s,c,chi)
   ! chib(G,G') = c sum_AB s_A v(G,A) [z - M]^-1_AB s_B conj(v(G',B))
   complex(DP), intent(in)  :: z,M(:,:),v(:,:),s(:)
   real(DP),    intent(in)  :: c
   complex(DP), intent(out) :: chi(:,:)
   complex(DP), allocatable :: W(:,:),X(:,:),B(:,:)
   integer,     allocatable :: ipiv(:)
   integer :: nh,ng,ig,i,info
   nh=size(s); ng=size(v,1)
   allocate(W(nh,nh),X(nh,ng),B(ng,nh),ipiv(nh))
   do ig=1,ng
     X(:,ig)=s*conjg(v(ig,:))
     B(ig,:)=s*v(ig,:)
   enddo
   W=-M
   do i=1,nh
     W(i,i)=W(i,i)+z
   enddo
   call zgesv(nh,ng,W,nh,ipiv,X,nh,info)
   call zgemm('N','N',ng,ng,nh,cmplx(c,0._DP,kind=DP),B,ng,X,nh,zero,chi,ng)
   deallocate(W,X,B,ipiv)
 end subroutine TDDFT_Chi_full_response
 !
 recursive subroutine TDDFT_Chi_kernel_from_response(chi0,chib,d,exch_head,fxc,P)
   ! fxc = chi0^-1 - P^-1,   P^-1 = chib^-1 + vbar
   ! done on D chi D (D = diag(d)), where vbar -> J = diag(0,1,..,1); J(1,1) = 1 too if exch_head
   complex(DP), intent(in)  :: chi0(:,:),chib(:,:),d(:)
   logical,     intent(in)  :: exch_head
   complex(DP), intent(out) :: fxc(:,:),P(:,:)
   complex(DP), allocatable :: A(:,:),B(:,:)
   integer :: ng,i,j
   ng=size(d)
   allocate(A(ng,ng),B(ng,ng))
   do j=1,ng
     do i=1,ng
       A(i,j)=d(i)*chi0(i,j)*d(j)
       B(i,j)=d(i)*chib(i,j)*d(j)
     enddo
   enddo
   call invert(A)
   call invert(B)
   do i=2,ng
     B(i,i)=B(i,i)+one
   enddo
   if (exch_head) B(1,1)=B(1,1)+one
   do j=1,ng
     do i=1,ng
       fxc(i,j)=d(i)*(A(i,j)-B(i,j))*d(j)
     enddo
   enddo
   A=B
   call invert(A)
   do j=1,ng
     do i=1,ng
       P(i,j)=A(i,j)/(d(i)*d(j))
     enddo
   enddo
   deallocate(A,B)
 end subroutine TDDFT_Chi_kernel_from_response
 !
 recursive subroutine TDDFT_Chi_full_dynamic_response(z,H0,fxc,v,s,c,d,resp)
   ! resp = sum_AB s_A conj(d_A) [z - H0 - Kf(z)]^-1_AB s_B d_B
   ! Kf(A,B) = c s_A s_B sum_GG' conj(v(G,A)) fxc(G,G') v(G',B)
   complex(DP), intent(in)  :: z,H0(:,:),fxc(:,:),v(:,:),s(:),d(:)
   real(DP),    intent(in)  :: c
   complex(DP), intent(out) :: resp
   complex(DP), allocatable :: M(:,:),fv(:,:),rhs(:)
   integer,     allocatable :: ipiv(:)
   integer :: nh,ng,i,info
   nh=size(s); ng=size(fxc,1)
   allocate(M(nh,nh),fv(ng,nh),rhs(nh),ipiv(nh))
   call zgemm('N','N',ng,nh,ng,one,fxc,ng,v,ng,zero,fv,ng)
   call zgemm('C','N',nh,nh,ng,cmplx(c,0._DP,kind=DP),v,ng,fv,ng,zero,M,nh)
   do i=1,nh
     M(:,i)=-M(:,i)*s(i)
     M(i,:)=M(i,:)*s(i)
   enddo
   M=M-H0
   do i=1,nh
     M(i,i)=M(i,i)+z
   enddo
   rhs=s*d
   call zgesv(nh,1,M,nh,ipiv,rhs,nh,info)
   resp=sum(s*conjg(d)*rhs)
   deallocate(M,fv,rhs,ipiv)
 end subroutine TDDFT_Chi_full_dynamic_response
 !
 subroutine invert(Mx)
   complex(DP), intent(inout) :: Mx(:,:)
   integer :: ipv(size(Mx,1)),info,m,lw
   complex(DP) :: wq(1)
   complex(DP), allocatable :: wk(:)
   m=size(Mx,1)
   call zgetrf(m,m,Mx,m,ipv,info)
   call zgetri(m,Mx,m,ipv,wq,-1,info)
   lw=max(1,int(real(wq(1))))
   allocate(wk(lw))
   call zgetri(m,Mx,m,ipv,wk,lw,info)
   deallocate(wk)
 end subroutine invert
end module TDDFT_Chi_dynamic_m
