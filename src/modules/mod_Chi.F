!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
module Chi_m
 !
 ! Irreducible polarizability in the plane-wave basis, with and without quasiparticle
 ! dressing, and the exchange-correlation kernel that maps one onto the other (ndb.Chi).
 !
 !   chi0(q,G,G';w)  irreducible bubble with the Kohn-Sham levels
 !   P   (q,G,G';w)  irreducible bubble with the quasiparticle levels applied to the
 !                   X energies through XfnQPdb (energies, and Z factors when the
 !                   database provides them: see X_irredux_residuals, Z_eh_occ)
 !   fxc (q,G,G';w)  = chi0^-1 - P^-1
 !
 ! Both bubbles are computed by X_irredux, so chi0 and P carry yambo's own convention:
 ! Adler-Wiser sum over transitions normalized by 1/(N_kbz V_cell), rho^*(G) rho(G') in
 ! the residual, GreenF from X_GreenF_analytical with the requested time ordering.
 ! The q=0 head and wings use the dipole q->0 limit at the small vector q0 (X%q0), as
 ! everywhere else in yambo.
 !
 ! The inversion is done in double precision on the symmetrized matrices
 !
 !   chi~ = D chi D,  D = diag( sqrt(4 pi)/|q+G| )     (|q+G| = bare_qpg, cutoff included)
 !
 ! whose head stays O(1) at q->0, and fxc = D (chi0~^-1 - P~^-1) D. The stored fxc is the
 ! bare (unsymmetrized) G-space kernel in Hartree atomic units. Its q->0 head is evaluated
 ! at q0 and diverges as 1/q^2, as the physical long-range kernel does.
 !
 ! Self-test: with no XfnQPdb the two bubbles are computed from identical levels through
 ! identical code and fxc must vanish exactly.
 !
 use pars,        ONLY:SP,schlen
 use descriptors, ONLY:IO_desc
 use electrons,   ONLY:levels
 !
#include<y_memory.h>
 !
 type(IO_desc)     :: Chi_desc
 !
 logical           :: l_Chi_fxc           =.FALSE.  ! ChiFxc input flag
 logical           :: Chi_header_done =.FALSE.
 logical           :: Chi_QP_applied  =.FALSE.
 integer           :: Chi_ng          =0
 integer           :: Chi_n_freqs     =0
 character(schlen) :: Chi_QP_action   ="none"
 character(schlen) :: Chi_convention  ="chi0,P: yambo Xo [1/(Nk V)]; fxc=chi0^-1-P^-1 [Ha], bare G-space"
 !
 ! Kohn-Sham levels, duplicated before the QP corrections are applied to the X levels
 !
 type(levels)      :: Chi_KS_levels
 !
 ! Per q-point (IBZ)
 !
 complex(SP), allocatable :: Chi_freqs(:)      ! (n_freqs)
 complex(SP), allocatable :: Chi_0(:,:,:)      ! chi0 (ng,ng,n_freqs)
 complex(SP), allocatable :: Chi_P(:,:,:)      ! P    (ng,ng,n_freqs)
 complex(SP), allocatable :: Chi_fxc(:,:,:)    ! fxc  (ng,ng,n_freqs)
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
end module Chi_m
