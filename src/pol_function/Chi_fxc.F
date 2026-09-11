!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details):
!
subroutine Chi_fxc_eval(iq,iq_mem,X,Xen,Xk,Xw,Dip)
 !
 ! Called by X_dielectric_matrix right after X_irredux has filled X_par(iq_mem) with the
 ! irreducible polarizability built from the X levels Xen. When XfnQPdb is active those
 ! levels carry the quasiparticle corrections, so X_par holds P. This routine
 !
 !   1. gathers P from the distributed X_par
 !   2. re-runs X_irredux on the same matrix with the Kohn-Sham levels saved before the
 !      corrections were applied (Chi_KS_levels), after resetting the pole tables so the
 !      transition energies are re-tabulated from those levels; this gives chi0
 !   3. puts P back into X_par (and its lower-triangle companion) so X_redux proceeds
 !      exactly as if nothing had happened
 !   4. on the master: fxc = chi0^-1 - P^-1 in double precision on the symmetrized
 !      matrices (mod_Chi), then writes ndb.Chi
 !
 use pars,          ONLY:SP,DP,pi,cZERO,cONE,schlen
 use units,         ONLY:HA2EV
 use com,           ONLY:msg
 use stderr,        ONLY:intc
 use LIVE_t,        ONLY:live_timing
 use frequency,     ONLY:w_samp
 use electrons,     ONLY:levels
 use R_lattice,     ONLY:bz_samp,bare_qpg
 use DIPOLES,       ONLY:DIPOLE_t
 use X_m,           ONLY:X_t,X_par,current_iq,X_par_lower_triangle,X_lower_triangle_matrix_in_use
 use matrix,        ONLY:PAR_matrix
 use Chi_m,         ONLY:Chi_KS_levels,Chi_ng,Chi_n_freqs,Chi_alloc,Chi_free,Chi_freqs,Chi_0,Chi_P,Chi_fxc,&
&                        Chi_header_done,Chi_QP_applied,Chi_QP_action
 use QP_CTL_m,      ONLY:QP_ctl_DB_user
 use parallel_m,    ONLY:master_cpu,PAR_COM_X_WORLD,PAR_COM_X_WORLD_RL_resolved
 use parallel_int,  ONLY:PP_redux_wait,PP_wait
 use parser_m,      ONLY:parser
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:OP_WR_CL,OP_APP_CL,REP
 use timing_m,      ONLY:timing
 !
#include<y_memory.h>
 !
 integer,       intent(in) :: iq,iq_mem
 type(X_t)                 :: X
 type(levels)              :: Xen
 type(bz_samp)             :: Xk
 type(w_samp)              :: Xw
 type(DIPOLE_t)            :: Dip
 !
 ! Work Space
 !
 integer                  :: ng,nw,iw,ig1,ig2,r1,r2,c1,c2,ID,io_err,info,lwork
 integer,     external    :: io_Chi
 real(SP)                 :: dev_P_chi0,fxc_body_max
 complex(SP), allocatable :: P_full(:,:,:),C0_full(:,:,:)
 complex(DP), allocatable :: A(:,:),B(:,:),work(:)
 complex(DP)              :: vs(X%ng),wq(1)
 integer,     allocatable :: ipiv(:)
 logical                  :: l_QP_db
 !
 call timing('Chi_fxc',OPR='start')
 !
 ng=X%ng
 nw=Xw%n_freqs
 !
 r1=X_par(iq_mem)%rows(1); r2=X_par(iq_mem)%rows(2)
 c1=X_par(iq_mem)%cols(1); c2=X_par(iq_mem)%cols(2)
 !
 YAMBO_ALLOC(P_full,(ng,ng,nw))
 YAMBO_ALLOC(C0_full,(ng,ng,nw))
 !
 ! 1. P from the distributed matrix (same gather as X_mat_filling)
 !
 call gather(X_par(iq_mem),P_full)
 !
 ! 2. chi0 with the Kohn-Sham levels, on the same matrix
 !
 call FREQUENCIES_reset(Xw,"coarse_grid")
 current_iq=-1
 call X_ALLOC_parallel(X_par(iq_mem),ng,nw,"XUP")
 call X_irredux(iq,"Xo[KS]",X_par(iq_mem),Chi_KS_levels,Xk,Xw,X,Dip)
 call gather(X_par(iq_mem),C0_full)
 !
 ! 3. P back in place for X_redux
 !
 X_par(iq_mem)%blc(r1:r2,c1:c2,:nw)=P_full(r1:r2,c1:c2,:)
 if (X_lower_triangle_matrix_in_use) X_par_lower_triangle%blc(c1:c2,r1:r2,:nw)=P_full(c1:c2,r1:r2,:)
 !
 ! 4. fxc and I/O, master only
 !
 if (master_cpu) then
   !
   Chi_ng     =ng
   Chi_n_freqs=nw
   call Chi_alloc(ng,nw)
   Chi_freqs(:)=Xw%p(:nw)
   Chi_0 (:,:,:)=C0_full
   Chi_P (:,:,:)=P_full
   !
   ! D = sqrt(4 pi)/|q+G|, the same Coulomb factor X_redux uses in the RPA kernel
   !
   do ig1=1,ng
     vs(ig1)=sqrt(4._DP*pi)/cmplx(bare_qpg(iq,ig1),kind=DP)
   enddo
   !
   YAMBO_ALLOC(A,(ng,ng))
   YAMBO_ALLOC(B,(ng,ng))
   YAMBO_ALLOC(ipiv,(ng))
   call zgetri(ng,A,ng,ipiv,wq,-1,info)
   lwork=max(int(real(wq(1))),ng)
   YAMBO_ALLOC(work,(lwork))
   !
   do iw=1,nw
     !
     do ig2=1,ng
       do ig1=1,ng
         A(ig1,ig2)=vs(ig1)*cmplx(C0_full(ig1,ig2,iw),kind=DP)*vs(ig2)
         B(ig1,ig2)=vs(ig1)*cmplx(P_full (ig1,ig2,iw),kind=DP)*vs(ig2)
       enddo
     enddo
     !
     call zgetrf(ng,ng,A,ng,ipiv,info)
     if (info/=0) call error('[Chi] chi0 is singular at q='//trim(intc(iq))//' w='//trim(intc(iw)))
     call zgetri(ng,A,ng,ipiv,work,lwork,info)
     call zgetrf(ng,ng,B,ng,ipiv,info)
     if (info/=0) call error('[Chi] P is singular at q='//trim(intc(iq))//' w='//trim(intc(iw)))
     call zgetri(ng,B,ng,ipiv,work,lwork,info)
     !
     do ig2=1,ng
       do ig1=1,ng
         Chi_fxc(ig1,ig2,iw)=cmplx(vs(ig1)*(A(ig1,ig2)-B(ig1,ig2))*vs(ig2),kind=SP)
       enddo
     enddo
     !
   enddo
   !
   YAMBO_FREE(work)
   YAMBO_FREE(ipiv)
   YAMBO_FREE(B)
   YAMBO_FREE(A)
   !
   ! Report
   !
   dev_P_chi0=maxval(abs(P_full-C0_full))/max(maxval(abs(C0_full)),tiny(1._SP))
   fxc_body_max=0._SP
   if (ng>1) fxc_body_max=maxval(abs(Chi_fxc(2:,2:,1)))
   !
   call msg('r',' [Chi] Current Q-pt index                ',iq)
   call msg('r',' [Chi] max|P-chi0|/max|chi0|             ',dev_P_chi0)
   call msg('r',' [Chi] fxc(q,0,0;w1)   [Ha]              ',(/real(Chi_fxc(1,1,1)),aimag(Chi_fxc(1,1,1))/))
   if (ng>1) call msg('r',' [Chi] max|fxc| body (G,Gp>0;w1) [Ha]    ',fxc_body_max)
   !
   ! I/O
   !
   if (.not.Chi_header_done) then
     call parser('XfnQPdb',l_QP_db)
     Chi_QP_applied=l_QP_db
     Chi_QP_action ="none"
     if (l_QP_db) Chi_QP_action=trim(QP_ctl_DB_user(1)%action)
     call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),ID=ID)
     io_err=io_Chi(X,Xw,ID)
     Chi_header_done=.TRUE.
   endif
   !
   call io_control(ACTION=OP_APP_CL,COM=REP,SEC=(/2*iq,2*iq+1/),ID=ID)
   io_err=io_Chi(X,Xw,ID)
   !
   call Chi_free()
   !
 endif
 !
 YAMBO_FREE(P_full)
 YAMBO_FREE(C0_full)
 !
 call PP_wait(COMM=PAR_COM_X_WORLD%COMM)
 !
 call timing('Chi_fxc',OPR='stop')
 !
 contains
   !
   subroutine gather(M,full)
     type(PAR_matrix)          :: M
     complex(SP)               :: full(ng,ng,nw)
     integer                   :: iw_
     full=cZERO
     if (PAR_COM_X_WORLD_RL_resolved%CPU_id==0) then
       full(M%rows(1):M%rows(2),M%cols(1):M%cols(2),:)=M%blc(:,:,:nw)
     endif
     do iw_=1,nw
       call PP_redux_wait(full(:,:,iw_),IN_PLACE=.true.,COMM=PAR_COM_X_WORLD%COMM)
     enddo
   end subroutine gather
   !
end subroutine Chi_fxc_eval
