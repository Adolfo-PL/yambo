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
 ! In LEVELS mode both bubbles are computed by X_irredux. G0/DYSON mode
 ! instead uses explicit finite-q density vertices and spectral band pairs.
 ! Both routes use yambo's normalization and vertex convention:
 ! Adler-Wiser sum over transitions normalized by 1/(N_kbz V_cell), rho^*(G) rho(G') in
 ! the residual. LEVELS uses X_GreenF_analytical with the requested ordering.
 ! In LEVELS mode the q=0 head and wings use the dipole limit at q0 (X%q0), as
 ! everywhere else in yambo.
 !
 ! The inversion is done in double precision on the symmetrized matrices
 !
 !   chi~ = D chi D,  D = diag( sqrt(4 pi)/|q+G| )     (|q+G| = bare_qpg, cutoff included)
 !
 ! whose head stays O(1) at q->0, and fxc = D (chi0~^-1 - P~^-1) D. The stored fxc is the
 ! bare (unsymmetrized) G-space kernel in Hartree atomic units. Its q->0 head is evaluated
 ! at q0 in LEVELS mode. G0/DYSON currently excludes q=0 and a response vertex.
 !
 ! Self-test: with no XfnQPdb the two bubbles are computed from identical levels through
 ! identical code and fxc must vanish exactly.
 !
 use pars,        ONLY:SP,DP,schlen
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
 ! LEVELS preserves the existing QP-energy path. G0/DYSON are independent
 ! finite-q band-pair bubbles; DYSON reads conventional retarded ndb.G.
 character(schlen) :: Chi_G_mode="LEVELS",Chi_G_db=""
 character(schlen) :: Chi_G_axis="IMAG"
 integer :: Chi_G_nu_steps=16
 real(SP) :: Chi_G_nu_range(2)=(/0._SP,1._SP/)
 real(SP) :: Chi_G_norm_tol=0.05_SP,Chi_G_check_tol=1.E-4_SP
 real(SP) :: Chi_rcond_floor=1.E-12_SP,Chi_G_mu=0._SP,Chi_G_kbt=0._SP
 real(SP) :: Chi_G0_error=0._SP
 real(SP), allocatable :: Chi_rcond(:,:),Chi_G_norm(:,:,:)
 real(DP), allocatable :: Chi_G_energy(:,:,:,:),Chi_G_weight(:,:,:,:),Chi_G_occupation(:,:,:,:)
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
     YAMBO_ALLOC(Chi_rcond,(2,nw))
   end subroutine Chi_alloc
   !
   subroutine Chi_free()
     YAMBO_FREE(Chi_freqs)
     YAMBO_FREE(Chi_0)
     YAMBO_FREE(Chi_P)
     YAMBO_FREE(Chi_fxc)
     YAMBO_FREE(Chi_rcond)
   end subroutine Chi_free
   ! Spectra live across q points and are released by the driver.
   subroutine Chi_G_free()
     if (allocated(Chi_G_energy)) deallocate(Chi_G_energy)
     if (allocated(Chi_G_weight)) deallocate(Chi_G_weight)
     if (allocated(Chi_G_occupation)) deallocate(Chi_G_occupation)
     if (allocated(Chi_G_norm)) deallocate(Chi_G_norm)
   end subroutine
   !
end module Chi_m
