!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details): 
!
module Wpot_m
 !
 ! Screened Coulomb potential W(q,G,G') in the plane-wave basis, stored on disk (ndb.W).
 !
 ! Conventions (Hartree atomic units):
 !
 !   v_G(q)   = 4 pi / |q+G|^2
 !   eps^-1   = 1 + v chi                      (chi is the reducible response function)
 !   W        = eps^-1 v = v + v chi v          (Hermitian)
 !            = v^1/2 [ 1 + X_mat ] v^1/2
 !
 ! where X_mat = v^1/2 chi v^1/2 = eps^-1_sym - 1 is what is stored in ndb.em1s/ndb.pp.
 !
 ! The Coulomb factors are the ones used by yambo in the self-energy (scatter_Gamp, mode 'c'),
 ! i.e. they include the RIM regularization of the q->0 head/wings when RIM is active. The
 ! bare |q+G| are stored as well so that the analytic 4 pi/(|q+G||q+G'|) form can be rebuilt.
 !
 ! Plasmon-pole (Godby-Needs) representation, when W is built from ndb.pp:
 !
 !   W(q,G,G';w) = v(q,G) delta_GG' + [ W_static(q,G,G') - v(q,G) delta_GG' ] Omega^2/(Omega^2 - w^2)
 !
 ! with Omega(q,G,G') the pole energies (yambo uses Re[Omega] in the GW self-energy).
 !
 use pars,        ONLY:SP,schlen
 use descriptors, ONLY:IO_desc
 !
#include<y_memory.h>
 !
 ! Header
 !
 type(IO_desc)     :: W_desc
 logical           :: W_is_PPA   =.FALSE.
 integer           :: W_ng       =0
 integer           :: W_n_freqs  =0
 character(schlen) :: W_screening="none"       ! "em1s" or "pp"
 character(schlen) :: W_convention="W = eps^-1 v = v + v chi v, eps^-1 = 1 + v chi [Ha]"
 !
 ! RIM-W (ndb.RIM_W) corrections folded in W_static (real part) and W_ppa_E (q=0) for
 ! G,G' <= W_RIM_W_ng, exactly as QP_ppa_cohsex uses them in the plasmon-pole self-energy.
 ! They are applied to the plasmon-pole W only (yambo does not use RIM-W in COHSEX).
 !
 logical           :: W_RIM_W    =.FALSE.
 integer           :: W_RIM_W_ng =0
 !
 ! Per q-point (IBZ) quantities
 !
 complex(SP), allocatable :: W_freqs(:)        ! frequencies of the screening DB (n_freqs)
 complex(SP), allocatable :: W_vcoul(:)        ! v(q,G) (ng)
 complex(SP), allocatable :: W_epsm1(:,:)      ! eps^-1_sym(q,G,G';w=0) = 1 + X_mat (ng,ng)
 complex(SP), allocatable :: W_static(:,:)     ! W(q,G,G';w=0) (ng,ng)
 complex(SP), allocatable :: W_ppa_E(:,:)      ! Plasmon-pole energies Omega(q,G,G') (ng,ng)
 !
 ! Plasmon-pole statistics (per q)
 !
 integer           :: W_PPcond_rate=0
 integer           :: W_TO_rate    =0
 real(SP)          :: W_PP_err     =0._SP
 !
 contains
   !
   subroutine W_alloc(ng,n_freqs,l_PPA)
     integer, intent(in) :: ng,n_freqs
     logical, intent(in) :: l_PPA
     call W_free()
     YAMBO_ALLOC(W_freqs,(n_freqs))
     YAMBO_ALLOC(W_vcoul,(ng))
     YAMBO_ALLOC(W_epsm1,(ng,ng))
     YAMBO_ALLOC(W_static,(ng,ng))
     if (l_PPA) then
       YAMBO_ALLOC(W_ppa_E,(ng,ng))
     endif
   end subroutine W_alloc
   !
   subroutine W_free()
     YAMBO_FREE(W_freqs)
     YAMBO_FREE(W_vcoul)
     YAMBO_FREE(W_epsm1)
     YAMBO_FREE(W_static)
     YAMBO_FREE(W_ppa_E)
   end subroutine W_free
   !
end module Wpot_m
