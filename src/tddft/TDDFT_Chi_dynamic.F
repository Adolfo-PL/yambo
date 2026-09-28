!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
! Casida response with a frequency-dependent kernel fxc(G,G';z), one frequency at
! a time. Only BLAS/LAPACK is used, so tests/chi_dynamic builds this file alone.
! The routines are recursive: K_Chi_dynamic calls them from OpenMP threads.
!
! Conventions are those of K_kernel/K_stored_in_a_big_matrix/K_diago_*:
!
!  * H0 is Yambo's BS_mat without the fxc part: transition energies plus the
!    static kernel (exchange). N x N in TDA, 2N x 2N with coupling, where the
!    coupling blocks carry the factor i of K_kernel ("C" block times cI).
!  * The fxc part of a block is  Kf_xy = c diag(s) V_x^H fxc(z) V_y diag(s),
!    with V_res the excitation vertices (O_tddft_R in mode "R"), V_cpl the
!    de-excitation ones (mode "C"), s = sqrt(f) and c = Co/(4 pi) of K_kernel.
!    For a Hermitian static fxc this is exactly what K_kernel adds to BS_mat.
!  * a = d sqrt(f) and b = conj(d) sqrt(f) are the resonant residual vectors of
!    K_diago_left/right_residuals. The anti-resonant ones are (i b, i a).
!
! Returned is  resp = B^T [z - M(z)]^-1 A  with
!
!    coupling: M = H0 + [[ Kf_rr, i Kf_rc ], [ i Kf_cr, -Kf_cc ]],  A=(a,i b), B=(b,i a)
!    TDA     : the two diagonal blocks alone, the anti-resonant one being
!              -conj(H0) - Kf_cc (A = -R* of K_stored_in_a_big_matrix).
!
! Resp_ii = -Co_factor*resp is what K_diago_response_functions accumulates as
! -Co_factor*sum_l R_L R_R/(z-E_l) (TDA: plus the anti-resonant mirror). With a
! z-independent Hermitian fxc the two agree identically.
!
! Every block uses the same retarded fxc(z), z = w + i eta, as the pair-space
! Dyson equation chi = [chi0^-1 - K(z)]^-1 requires: the anti-resonant pairs see
! fxc at +w through their own vertices, not a mirror of the resonant block.
!
module TDDFT_Chi_dynamic_m
 use pars, only: DP
 implicit none
 private
 public :: TDDFT_Chi_dynamic_kernel,TDDFT_Chi_dynamic_response
 complex(DP), parameter :: zero=(0._DP,0._DP),one=(1._DP,0._DP),ci=(0._DP,1._DP)
contains
 !
 recursive subroutine TDDFT_Chi_dynamic_kernel(fxc,v_left,v_right,s,c,kernel,ierr)
   ! kernel(t,t') = c s_t [ sum_GG' conj(v_left(G,t)) fxc(G,G') v_right(G',t') ] s_t'
   complex(DP), intent(in)  :: fxc(:,:),v_left(:,:),v_right(:,:)
   real(DP),    intent(in)  :: s(:),c
   complex(DP), intent(out) :: kernel(:,:)
   integer,     intent(out) :: ierr
   complex(DP), allocatable :: fv(:,:)
   integer :: ng,nt,it
   ! ierr: 1 fxc not square, 2 vertex G size, 3 vertex/s/kernel transition size
   ierr=1
   ng=size(fxc,1)
   if (ng<1.or.size(fxc,2)/=ng) return
   ierr=2
   if (size(v_left,1)/=ng.or.size(v_right,1)/=ng) return
   ierr=3
   nt=size(v_right,2)
   if (nt<1.or.size(v_left,2)/=nt.or.size(s)/=nt) return
   if (size(kernel,1)/=nt.or.size(kernel,2)/=nt) return
   allocate(fv(ng,nt))
   call zgemm('N','N',ng,nt,ng,one,fxc,ng,v_right,ng,zero,fv,ng)
   call zgemm('C','N',nt,nt,ng,cmplx(c,0._DP,kind=DP),v_left,ng,fv,ng,zero,kernel,nt)
   deallocate(fv)
   do it=1,nt
     kernel(:,it)=kernel(:,it)*s(it)
     kernel(it,:)=kernel(it,:)*s(it)
   enddo
   ierr=0
 end subroutine TDDFT_Chi_dynamic_kernel
 !
 recursive subroutine TDDFT_Chi_dynamic_response(z,H0,fxc,v_res,v_cpl,s,c,a,b,coupling,resp,ierr)
   complex(DP), intent(in)  :: z,H0(:,:),fxc(:,:),v_res(:,:),v_cpl(:,:),a(:),b(:)
   real(DP),    intent(in)  :: s(:),c
   logical,     intent(in)  :: coupling
   complex(DP), intent(out) :: resp
   integer,     intent(out) :: ierr
   complex(DP), allocatable :: M(:,:),K(:,:),rhs(:)
   integer,     allocatable :: ipiv(:)
   integer :: n,nh,i,info
   ! ierr: 10+ kernel size error, 20 H0 size, 30+info singular z-M (LAPACK info>0), 40 LAPACK argument
   resp=zero
   n=size(a)
   nh=n
   if (coupling) nh=2*n
   ierr=20
   if (size(b)/=n.or.size(H0,1)/=nh.or.size(H0,2)/=nh) return
   allocate(K(n,n),M(nh,nh),rhs(nh),ipiv(nh))
   if (coupling) then
     M=-H0
     call TDDFT_Chi_dynamic_kernel(fxc,v_res,v_res,s,c,K,ierr)
     if (ierr/=0) goto 1
     M(:n,:n)=M(:n,:n)-K
     call TDDFT_Chi_dynamic_kernel(fxc,v_res,v_cpl,s,c,K,ierr)
     if (ierr/=0) goto 1
     M(:n,n+1:)=M(:n,n+1:)-ci*K
     call TDDFT_Chi_dynamic_kernel(fxc,v_cpl,v_res,s,c,K,ierr)
     if (ierr/=0) goto 1
     M(n+1:,:n)=M(n+1:,:n)-ci*K
     call TDDFT_Chi_dynamic_kernel(fxc,v_cpl,v_cpl,s,c,K,ierr)
     if (ierr/=0) goto 1
     M(n+1:,n+1:)=M(n+1:,n+1:)+K
     do i=1,nh
       M(i,i)=M(i,i)+z
     enddo
     rhs(:n)=a
     rhs(n+1:)=ci*b
     call solve(nh,M,rhs,ierr)
     if (ierr/=0) goto 1
     resp=sum(b*rhs(:n))+sum(ci*a*rhs(n+1:))
   else
     ! Resonant block: b^T (z - H0 - Kf_rr)^-1 a
     call TDDFT_Chi_dynamic_kernel(fxc,v_res,v_res,s,c,K,ierr)
     if (ierr/=0) goto 1
     M=-H0-K
     do i=1,n
       M(i,i)=M(i,i)+z
     enddo
     rhs=a
     call solve(n,M,rhs,ierr)
     if (ierr/=0) goto 1
     resp=sum(b*rhs)
     ! Anti-resonant block: (i a)^T (z + conj(H0) + Kf_cc)^-1 (i b)
     call TDDFT_Chi_dynamic_kernel(fxc,v_cpl,v_cpl,s,c,K,ierr)
     if (ierr/=0) goto 1
     M=conjg(H0)+K
     do i=1,n
       M(i,i)=M(i,i)+z
     enddo
     rhs=b
     call solve(n,M,rhs,ierr)
     if (ierr/=0) goto 1
     resp=resp-sum(a*rhs)
   endif
   ierr=0
1  if (ierr>0.and.ierr<10) ierr=10+ierr
   deallocate(K,M,rhs,ipiv)
 contains
   subroutine solve(m_dim,mat,x,err)
     integer,     intent(in)    :: m_dim
     complex(DP), intent(inout) :: mat(:,:),x(:)
     integer,     intent(out)   :: err
     call zgesv(m_dim,1,mat,m_dim,ipiv,x,m_dim,info)
     err=0
     if (info>0) err=30
     if (info<0) err=40
   end subroutine solve
 end subroutine TDDFT_Chi_dynamic_response
end module TDDFT_Chi_dynamic_m
