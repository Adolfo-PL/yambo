      program casida_simple

!     Exciton energies from the Casida matrix written by fkq_simple.sh (or fkq_dump.sh).
!
!     The matrix at frequency w, in eV (N transitions, size 2N):
!
!        M(w) = | diag(E)  + X_rr + F_rr(w)      X_rc + F_rc(w)          |
!               | X_cr + F_cr(w)                diag(E_ar) + X_cc + F_cc(w) |
!
!     E   = transition energies (resonant rows, from the _transitions file)
!     E_ar = energies of the anti-resonant rows (= -E at q1, from the _antiresonant file at q>1)
!     X   = Hartree coupling, F = f_xc coupling (your F_kq), both from the _w file
!
!     Do NOT drop the rc, cr, cc blocks (Tamm-Dancoff): the kernel was made from the full
!     (coupled) BSE, and without them it does not give the exciton.
!
!     Usage:   ./casida_simple  transitions_file  w_file  [more w_files ...]
!
!     One w file : eigenvalues -> eigen.dat, absorption -> absorption.dat
!     Several    : (give them in frequency order: ls -v) lowest bright eigenvalue
!                  lambda(w) of each -> scan.dat, and the
!                  exciton = the w where lambda(w) = w (linear interpolation)
!
!     Build:   gfortran -O2 -ffpe-summary=none casida_simple.f90 -o casida_simple -llapack -lblas

      implicit none
      integer, parameter :: dp=kind(1.0d0)
      integer :: n, nn, i, j, l, it, jt, iw, n_w, info, lwork, ios, l_bright
      character(len=512) :: f_tr, f_ar, f_w, line
      logical :: have_ar
      real(dp) :: v(16), col(15), z_re, z_im, gap, s_max, omega, eps2, eta
      real(dp), allocatable :: E(:), E_ar(:), strength(:), w_list(:), lam_list(:)
      complex(dp) :: ii, a_l, b_l, n_l
      complex(dp), allocatable :: a(:), A_ar(:), B_ar(:), Avec(:), Bvec(:)
      complex(dp), allocatable :: M(:,:), VL(:,:), VR(:,:), lam(:), st(:), work(:)
      real(dp), allocatable :: rwork(:)
      integer, allocatable :: order(:)

      ii=dcmplx(0.0d0,1.0d0)
      eta=0.05d0                  ! Lorentzian width of absorption.dat [eV]

!     ---------------- transitions: E_t and the dipoles a_t ----------------

      call get_command_argument(1,f_tr)
      n_w=command_argument_count()-1

      n=0
      open(1,File=trim(f_tr),status='old')
      do
         read(1,'(a)',iostat=ios) line
         if (ios.ne.0) exit
         if (line(1:1).ne.'#' .and. len_trim(line).gt.0) n=n+1
      end do
      close(1)
      nn=2*n
      write(*,*) 'transitions N =',n,'   matrix size 2N =',nn

      allocate(E(n),E_ar(n),a(n),A_ar(n),B_ar(n),Avec(nn),Bvec(nn))
      allocate(M(nn,nn),VL(nn,nn),VR(nn,nn),lam(nn),st(nn),strength(nn),order(nn))
      allocate(w_list(n_w),lam_list(n_w))
      lwork=4*nn
      allocate(work(lwork),rwork(2*nn))

!     columns: t T ik_bz ik_ibz k1 k2 k3 v c E_t f_t Re(a_t) Im(a_t)
      open(1,File=trim(f_tr),status='old')
      i=0
      do
         read(1,'(a)',iostat=ios) line
         if (ios.ne.0) exit
         if (line(1:1).eq.'#' .or. len_trim(line).eq.0) cycle
         i=i+1
         read(line,*) col(1:13)
         E(i)=col(10)
         a(i)=dcmplx(col(12),col(13))
      end do
      close(1)
      gap=minval(E)

!     anti-resonant rows: own file at finite q, mirror of the resonant ones at q1
      f_ar=f_tr(1:index(f_tr,'_transitions')-1)//'_antiresonant'
      inquire(File=trim(f_ar),exist=have_ar)
      if (have_ar) then
!        columns: t T ik_bz ik_ibz k1 k2 k3 v c E_t f_t Re(A) Im(A) Re(B) Im(B)
         open(1,File=trim(f_ar),status='old')
         i=0
         do
            read(1,'(a)',iostat=ios) line
            if (ios.ne.0) exit
            if (line(1:1).eq.'#' .or. len_trim(line).eq.0) cycle
            i=i+1
            read(line,*) col(1:15)
            E_ar(i)=col(10)
            A_ar(i)=dcmplx(col(12),col(13))
            B_ar(i)=dcmplx(col(14),col(15))
         end do
         close(1)
         write(*,*) 'anti-resonant rows read from ',trim(f_ar)
      else
         E_ar=-E
         A_ar=ii*dconjg(a)
         B_ar=ii*a
      end if
!     light couples to the rows through A (in) and B (out)
      Avec(1:n)=a
      Avec(n+1:nn)=A_ar
      Bvec(1:n)=dconjg(a)
      Bvec(n+1:nn)=B_ar

      write(*,*) 'gap (lowest transition) =',gap,' eV'

!     ---------------- loop over the frequency files ----------------

      do iw=1,n_w
         call get_command_argument(iw+1,f_w)

!        header: '# Casida matrix at ndb.Chi frequency  iw, z =  Re  Im eV'
         M=dcmplx(0.0d0,0.0d0)
         open(1,File=trim(f_w),status='old')
         do
            read(1,'(a)') line
            if (line(1:1).ne.'#') then
               backspace(1)
               exit
            end if
            if (index(line,'z =').gt.0) read(line(index(line,'z =')+3:),*) z_re, z_im
         end do
!        each line: t t'  then Re/Im of X_rr F_rr X_rc F_rc X_cr F_cr X_cc F_cc
         do l=1,n*n
            read(1,*) it, jt, v
            M(it  ,jt  )=dcmplx(v(1) +v(3) ,v(2) +v(4))
            M(it  ,n+jt)=dcmplx(v(5) +v(7) ,v(6) +v(8))
            M(n+it,jt  )=dcmplx(v(9) +v(11),v(10)+v(12))
            M(n+it,n+jt)=dcmplx(v(13)+v(15),v(14)+v(16))
         end do
         close(1)
         do i=1,n
            M(i,i)=M(i,i)+E(i)
            M(n+i,n+i)=M(n+i,n+i)+E_ar(i)
         end do

!        diagonalize (non-Hermitian: left and right eigenvectors)
         call ZGEEV('V','V',nn,M,nn,lam,VL,nn,VR,nn,work,lwork,rwork,info)
         if (info.ne.0) write(*,*) 'ZGEEV failed, info =',info

!        oscillator strength of eigenvalue l:  (B.R_l) (L_l^H A) / (L_l^H R_l)
         do l=1,nn
            b_l=sum(Bvec*VR(:,l))
            a_l=sum(dconjg(VL(:,l))*Avec)
            n_l=sum(dconjg(VL(:,l))*VR(:,l))
            st(l)=b_l*a_l/n_l
            strength(l)=abs(st(l))
         end do

!        positive eigenvalues, lowest first
         j=0
         do l=1,nn
            if (dreal(lam(l)).gt.1.0d-6) then
               j=j+1
               order(j)=l
            end if
         end do
         call sort_by_energy(order,j,lam,nn)
         s_max=maxval(strength(order(1:j)))
         l_bright=order(1)
         do i=1,j
            if (strength(order(i)).gt.0.01d0*s_max) then
               l_bright=order(i)
               exit
            end if
         end do

         if (n_w.eq.1) then
            write(*,*) 'frequency z =',z_re,z_im,' eV'
            open(2,File='eigen.dat')
            write(2,*) '# Re E [eV]   Im E [eV]   strength/max   E - gap [eV]'
            do i=1,j
               l=order(i)
               write(2,'(4f14.6)') dreal(lam(l)), dimag(lam(l)), strength(l)/s_max, dreal(lam(l))-gap
               if (i.le.8) write(*,'(a,4f12.6)') '   E, Im E, strength, E-gap:', &
                  dreal(lam(l)), dimag(lam(l)), strength(l)/s_max, dreal(lam(l))-gap
            end do
            close(2)
            write(*,*) 'lowest bright eigenvalue =',dreal(lam(l_bright)),' eV'
            write(*,*) 'binding energy (gap - bright) =',gap-dreal(lam(l_bright)),' eV'

!           absorption, up to a constant: eps2(w) = -Im sum_l st_l / (w + i eta - E_l)
            open(2,File='absorption.dat')
            do i=0,800
               omega=0.005d0*i
               eps2=0.0d0
               do l=1,j
                  eps2=eps2-dimag(st(order(l))/(omega+ii*eta-lam(order(l))))
               end do
               write(2,*) omega, eps2
            end do
            close(2)
         else
            w_list(iw)=z_re
            lam_list(iw)=dreal(lam(l_bright))
            write(*,'(a,f10.5,a,f10.5,a,f10.5)') ' w =',z_re,'   lambda_bright =',lam_list(iw), &
                                                 '   lambda - w =',lam_list(iw)-z_re
         end if
      end do

!     ---------------- several frequencies: solve lambda(w) = w ----------------

      if (n_w.gt.1) then
         if (abs(z_im).gt.1.0d-6) write(*,*) 'warning: damped frequencies; use undamped ones (BDmRange 0|0)'
         open(2,File='scan.dat')
         write(2,*) '# w [eV]   lambda_bright(w) [eV]'
         do iw=1,n_w
            write(2,'(2f14.6)') w_list(iw), lam_list(iw)
         end do
         close(2)
!        files given in frequency order; a sign change of lambda - w is the exciton,
!        unless lambda jumps up there (a pole of f_xc(w), not a solution)
         do iw=1,n_w-1
            if ((lam_list(iw)-w_list(iw))*(lam_list(iw+1)-w_list(iw+1)).le.0.0d0) then
               if (lam_list(iw+1)-lam_list(iw).gt.10.0d0*(w_list(iw+1)-w_list(iw))) then
                  write(*,*) 'pole of f_xc between',w_list(iw),' and',w_list(iw+1),' (not a solution)'
               else
                  omega=w_list(iw)-(lam_list(iw)-w_list(iw))*(w_list(iw+1)-w_list(iw))   &
                        /((lam_list(iw+1)-w_list(iw+1))-(lam_list(iw)-w_list(iw)))
                  write(*,*) 'exciton: lambda(w) = w at',omega,' eV'
                  write(*,*) 'binding energy (gap - exciton) =',gap-omega,' eV'
               end if
            end if
         end do
      end if

      stop
      end


      subroutine sort_by_energy(order,m,lam,nn)
!     insertion sort of order(1:m) by Re lam
      implicit none
      integer :: m, nn, order(nn), i, j, k
      complex(kind(1.0d0)) :: lam(nn)
      do i=2,m
         k=order(i)
         j=i-1
         do while (j.ge.1)
            if (dreal(lam(order(j))).le.dreal(lam(k))) exit
            order(j+1)=order(j)
            j=j-1
         end do
         order(j+1)=k
      end do
      end subroutine sort_by_energy
