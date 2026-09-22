!
! License-Identifier: GPL
!
! Project an unsymmetrized finite-q density kernel onto transition densities.
! The stored ndb.Chi matrix has units of energy times cell volume. Vectors
! passed from the BSE kernel are the raw Fourier overlap coefficients in
! Yambo's O_tddft_L/R convention, without Coulomb or volume scaling.
!
module TDDFT_Chi_Fourier_m
 use pars, only: DP
 implicit none
 private
 public :: TDDFT_Chi_pair_G,TDDFT_Chi_project_G,TDDFT_Chi_project_R
contains
 !
 subroutine TDDFT_Chi_pair_G(fxc,d_left,d_right,value,ierr)
   ! The BSE caller applies its own spin/k/volume factor. In K_kernel that
   ! factor includes 4*pi, so divide this result by 4*pi before H_x*Co.
   complex(DP), intent(in)  :: fxc(:,:),d_left(:),d_right(:)
   complex(DP), intent(out) :: value
   integer,     intent(out) :: ierr
   complex(DP) :: work(size(d_left))
   integer :: ng

   ierr=1
   value=(0._DP,0._DP)
   ng=size(d_left)
   if (ng<1) return
   if (size(d_right)/=ng.or.size(fxc,1)/=ng.or.size(fxc,2)/=ng) return
   work=matmul(fxc,d_right)
   value=dot_product(d_left,work)
   ierr=0
 end subroutine TDDFT_Chi_pair_G
 !
 subroutine TDDFT_Chi_project_G(fxc,vertex,volume,kernel,ierr)
   complex(DP), intent(in)  :: fxc(:,:),vertex(:,:)
   real(DP),    intent(in)  :: volume
   complex(DP), intent(out) :: kernel(:,:)
   integer,     intent(out) :: ierr
   complex(DP) :: work(size(vertex,1))
   integer :: ng,nt,it,jt

   ierr=1
   kernel=(0._DP,0._DP)
   ng=size(vertex,1)
   nt=size(vertex,2)
   if (ng<1.or.nt<1.or.volume<=0._DP) return
   if (size(fxc,1)/=ng.or.size(fxc,2)/=ng) return
   if (size(kernel,1)/=nt.or.size(kernel,2)/=nt) return

   do jt=1,nt
     work=matmul(fxc,vertex(:,jt))
     do it=1,nt
       kernel(it,jt)=dot_product(vertex(:,it),work)/volume
     enddo
   enddo
   ierr=0
 end subroutine TDDFT_Chi_project_G
 !
 subroutine TDDFT_Chi_project_R(fxc,qg_cart,r_cart,volume,transition_r,kernel,ierr)
   !
   ! This is the inverse-transform action of fxc on each ket transition:
   ! d_j(G') = integral dr' exp[-i(q+G')r'] T_j(r')
   ! h_j(r) = (1/volume) sum_G exp[i(q+G)r] sum_G' fxc(G,G') d_j(G')
   ! K_ij = integral dr conjg(T_i(r)) h_j(r).
   ! The uniform cell grid uses dV=volume/nr. Only O(ng+nr) scratch is used;
   ! the two-position real-space kernel is never materialized.
   !
   complex(DP), intent(in)  :: fxc(:,:),transition_r(:,:)
   real(DP),    intent(in)  :: qg_cart(:,:),r_cart(:,:),volume
   complex(DP), intent(out) :: kernel(:,:)
   integer,     intent(out) :: ierr
   complex(DP) :: density_g(size(fxc,1)),field_g(size(fxc,1)),field_r
   real(DP) :: dvol,phase
   integer :: ng,nr,nt,ig,ir,it,jt

   ierr=1
   kernel=(0._DP,0._DP)
   ng=size(fxc,1)
   nr=size(r_cart,2)
   nt=size(transition_r,2)
   if (ng<1.or.nr<1.or.nt<1.or.volume<=0._DP) return
   if (size(fxc,2)/=ng.or.size(qg_cart,1)/=3.or.size(qg_cart,2)/=ng) return
   if (size(r_cart,1)/=3.or.size(transition_r,1)/=nr) return
   if (size(kernel,1)/=nt.or.size(kernel,2)/=nt) return
   dvol=volume/real(nr,DP)

   do jt=1,nt
     density_g=(0._DP,0._DP)
     do ig=1,ng
       do ir=1,nr
         phase=dot_product(qg_cart(:,ig),r_cart(:,ir))
         density_g(ig)=density_g(ig)+&
&          exp(cmplx(0._DP,-phase,kind=DP))*transition_r(ir,jt)*dvol
       enddo
     enddo
     field_g=matmul(fxc,density_g)
     do ir=1,nr
       field_r=(0._DP,0._DP)
       do ig=1,ng
         phase=dot_product(qg_cart(:,ig),r_cart(:,ir))
         field_r=field_r+exp(cmplx(0._DP,phase,kind=DP))*field_g(ig)
       enddo
       field_r=field_r/volume
       do it=1,nt
         kernel(it,jt)=kernel(it,jt)+conjg(transition_r(ir,it))*field_r*dvol
       enddo
     enddo
   enddo
   ierr=0
 end subroutine TDDFT_Chi_project_R
end module TDDFT_Chi_Fourier_m
