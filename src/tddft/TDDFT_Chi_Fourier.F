!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Reading copy: the numerics only.
!
! LEGEND
!   fxc      (ng,ng)  kernel fxc(G,G') of ndb.Chi at z = 0
!   d_left   (ng)     left vertex  of transition t  (O_tddft_L)
!   d_right  (ng)     right vertex of transition t' (O_tddft_R)
!   value             sum_GG' conj(d_left(G)) fxc(G,G') d_right(G')
!
module TDDFT_Chi_Fourier_m
 use pars, only: DP
 implicit none
 private
 public :: TDDFT_Chi_pair_G
contains
 !
 subroutine TDDFT_Chi_pair_G(fxc,d_left,d_right,value)
   ! value = d_left^H fxc d_right.  K_kernel divides it by 4 pi and multiplies by Co
   complex(DP), intent(in)  :: fxc(:,:),d_left(:),d_right(:)
   complex(DP), intent(out) :: value
   complex(DP) :: work(size(d_left))
   work=matmul(fxc,d_right)
   value=dot_product(d_left,work)
 end subroutine TDDFT_Chi_pair_G
end module TDDFT_Chi_Fourier_m
