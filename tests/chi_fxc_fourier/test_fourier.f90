program test_chi_fourier
  use pars, only: DP
  use TDDFT_Chi_Fourier_m, only: TDDFT_Chi_pair_G,TDDFT_Chi_project_G,TDDFT_Chi_project_R
  implicit none
  integer, parameter :: ng=3,nt=2,nr=16,nw=2
  real(DP), parameter :: pi=acos(-1._DP),volume=8.5_DP
  real(DP) :: qg(3,ng),r(3,nr),q,phase,dvol,max_error
  complex(DP) :: vertex(ng,nt),transition_r(nr,nt),fxc(ng,ng,nw)
  complex(DP) :: from_g(nt,nt),from_r(nt,nt),from_double(nt,nt)
  complex(DP) :: real_kernel,pair,previous(nt,nt)
  integer :: ig,jg,ir,jr,it,jt,iw,ierr

  ! Nonzero q checks that its Bloch phase cancels only after both transforms.
  q=0.37_DP*2._DP*pi/volume
  qg=0._DP
  r=0._DP
  do ig=1,ng
    qg(1,ig)=q+2._DP*pi*real(ig-2,DP)/volume
  enddo
  do ir=1,nr
    r(1,ir)=volume*real(ir-1,DP)/real(nr,DP)
  enddo
  dvol=volume/real(nr,DP)

  do it=1,nt
    do ig=1,ng
      vertex(ig,it)=cmplx(0.17_DP*ig+0.23_DP*it,&
&       0.31_DP*ig-0.13_DP*it+0.07_DP*ig*it,kind=DP)
    enddo
  enddo
  do ir=1,nr
    do it=1,nt
      transition_r(ir,it)=(0._DP,0._DP)
      do ig=1,ng
        phase=dot_product(qg(:,ig),r(:,ir))
        transition_r(ir,it)=transition_r(ir,it)+&
&        vertex(ig,it)*exp(cmplx(0._DP,phase,kind=DP))/volume
      enddo
    enddo
  enddo

  do iw=1,nw
    do jg=1,ng
      do ig=1,ng
        ! A complex nonsymmetric matrix detects transposes and conjugation.
        fxc(ig,jg,iw)=cmplx(0.4_DP*ig-0.27_DP*jg+0.21_DP*iw,&
&         0.19_DP*ig*jg-0.13_DP*jg+0.17_DP*iw*ig,kind=DP)
      enddo
    enddo
  enddo

  previous=(0._DP,0._DP)
  do iw=1,nw
    call TDDFT_Chi_project_G(fxc(:,:,iw),vertex,volume,from_g,ierr)
    if (ierr/=0) error stop 'reciprocal projection rejected valid data'
    do jt=1,nt
      do it=1,nt
        call TDDFT_Chi_pair_G(fxc(:,:,iw),vertex(:,it),vertex(:,jt),pair,ierr)
        if (ierr/=0) error stop 'pair projection rejected valid data'
        if (abs(pair-volume*from_g(it,jt))>3.e-13_DP) &
&         error stop 'pair contraction has incorrect normalization'
      enddo
    enddo
    call TDDFT_Chi_project_R(fxc(:,:,iw),qg,r,volume,transition_r,from_r,ierr)
    if (ierr/=0) error stop 'streamed inverse transform rejected valid data'
    max_error=maxval(abs(from_r-from_g))
    if (max_error>3.e-13_DP) error stop 'real-space and reciprocal projections differ'

    ! Independently form the two-position kernel only in this tiny test.
    from_double=(0._DP,0._DP)
    do jr=1,nr
      do ir=1,nr
        real_kernel=(0._DP,0._DP)
        do jg=1,ng
          do ig=1,ng
            real_kernel=real_kernel+&
&            exp(cmplx(0._DP, dot_product(qg(:,ig),r(:,ir)),kind=DP))*&
&            fxc(ig,jg,iw)*&
&            exp(cmplx(0._DP,-dot_product(qg(:,jg),r(:,jr)),kind=DP))/volume
          enddo
        enddo
        do jt=1,nt
          do it=1,nt
            from_double(it,jt)=from_double(it,jt)+&
&            conjg(transition_r(ir,it))*real_kernel*transition_r(jr,jt)*dvol*dvol
          enddo
        enddo
      enddo
    enddo
    if (maxval(abs(from_double-from_g))>3.e-13_DP) &
&      error stop 'explicit inverse transform differs from reciprocal projection'
    if (iw>1.and.maxval(abs(from_g-previous))<1.e-2_DP) &
&      error stop 'frequency-dependent matrices were collapsed'
    previous=from_g
  enddo

  ! Invalid geometry must fail before any Fourier operation.
  call TDDFT_Chi_project_R(fxc(:,:,1),qg,r,0._DP,transition_r,from_r,ierr)
  if (ierr==0) error stop 'zero volume accepted'
  call TDDFT_Chi_project_G(fxc(:2,:2,1),vertex,volume,from_g,ierr)
  if (ierr==0) error stop 'mismatched density basis accepted'
  call TDDFT_Chi_pair_G(fxc(:2,:2,1),vertex(:,1),vertex(:,2),pair,ierr)
  if (ierr==0) error stop 'pair contraction accepted mismatched density basis'
  print *, 'PASS finite-q Fourier projection at two frequencies'
end program test_chi_fourier
