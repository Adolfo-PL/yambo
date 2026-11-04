!
! License-Identifier: GPL
!
! Copyright (C) 2026 Adolfo Partida
!
! Authors: Adolfo Partida, with AI assistance (Claude Code)
!
! Reading copy: only the data the exact-kernel path uses (ndb.Chi arrays and the transition vertices).
!
! LEGEND
!   Chi_ng, Chi_n_freqs   number of plane waves G and of frequencies of ndb.Chi
!   Chi_freqs             (n_freqs)        z = 0 first, then the BSE frequencies w + i eta
!   Chi_0, Chi_P, Chi_fxc (ng,ng,n_freqs)  chi0, P and fxc, unsymmetrized
!   Chi_dyn_V             (ng,N,2)         vertices rho_t(G) of the N resonant transitions: [:,:,1] excitation, [:,:,2] de-excitation
!   Chi_dyn_s             (N)              s_t = sqrt(f_t)
!   Chi_dyn_c             spin_occ/(Omega Nk) = Co/4pi
!
module Chi_m
 !
 use pars,        ONLY:SP,DP
 !
#include<y_memory.h>
 !
 integer           :: Chi_ng          =0
 integer           :: Chi_n_freqs     =0
 complex(SP), allocatable :: Chi_freqs(:)
 complex(SP), allocatable :: Chi_0(:,:,:)
 complex(SP), allocatable :: Chi_P(:,:,:)
 complex(SP), allocatable :: Chi_fxc(:,:,:)
 complex(DP), allocatable :: Chi_dyn_V(:,:,:)
 real(DP),    allocatable :: Chi_dyn_s(:)
 logical,     allocatable :: Chi_dyn_filled(:,:)
 real(DP)                 :: Chi_dyn_c=0._DP
 !
 contains
   !
   subroutine Chi_alloc(ng,nw)
     integer, intent(in) :: ng,nw
     call Chi_free()
     YAMBO_ALLOC(Chi_freqs,(nw))
     YAMBO_ALLOC(Chi_0,(ng,ng,nw))
     YAMBO_ALLOC(Chi_P,(ng,ng,nw))
     YAMBO_ALLOC(Chi_fxc,(ng,ng,nw))
   end subroutine Chi_alloc
   !
   subroutine Chi_free()
     YAMBO_FREE(Chi_freqs)
     YAMBO_FREE(Chi_0)
     YAMBO_FREE(Chi_P)
     YAMBO_FREE(Chi_fxc)
   end subroutine Chi_free
   !
   subroutine Chi_dyn_free()
     if (allocated(Chi_dyn_V)) deallocate(Chi_dyn_V)
     if (allocated(Chi_dyn_s)) deallocate(Chi_dyn_s)
     if (allocated(Chi_dyn_filled)) deallocate(Chi_dyn_filled)
     Chi_dyn_c=0._DP
   end subroutine Chi_dyn_free
   !
end module Chi_m
