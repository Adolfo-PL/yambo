! Real MPI collectives around the production bubble; synthetic host/IO only.
program test_mpi
 use mpi
 use pars
 use Chi_m
 use X_m
 use R_lattice
 use frequency
 use DIPOLES
 use parallel_m
 implicit none
 type(X_t)::x
 type(bz_samp)::k
 type(w_samp)::w
 type(DIPOLE_t)::d
 complex(SP)::response(2,2,3),expected(2,2),factor
 real(SP)::gap,tolerance
 integer::ierr,ik,iw,mode
 call MPI_INIT(ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
 master_cpu=myid==0
 PAR_COM_X_WORLD%COMM=MPI_COMM_WORLD
 PAR_COM_X_WORLD%n_CPU=ncpu
 PAR_COM_X_WORLD%CPU_id=myid
 allocate(PAR_IND_Xk_bz%element_1D(2))
 do ik=1,2
   PAR_IND_Xk_bz%element_1D(ik)=mod(ik-1,ncpu)==myid
 enddo
 call Chi_G_parallel_check(x)
 allocate(k%sstar(2,2),k%pt(2,3),qindx_X(2,2,2))
 k%sstar(:,1)=[1,2];k%sstar(:,2)=1;k%pt=0._SP
 qindx_X=1;qindx_X(2,:,1)=[2,1]
 allocate(Chi_KS_levels%E(2,2,1),Chi_KS_levels%f(2,2,1))
 Chi_KS_levels%E(1,:,:)=-1._SP;Chi_KS_levels%E(2,:,:)=1._SP
 Chi_KS_levels%f(1,:,:)=2._SP;Chi_KS_levels%f(2,:,:)=0._SP
 allocate(w%p(3))
 w%n_freqs=3;w%p=[cmplx(0._SP,0._SP,SP),cmplx(0._SP,1._SP,SP),cmplx(0._SP,2._SP,SP)]
 Chi_G_axis='IMAG'
 tolerance=200._SP*epsilon(1._SP)
 do mode=1,3
   select case(mode)
   case(1)
     Chi_G_mode='G0';gap=2._SP
   case(2)
     Chi_G_mode='COHSEX';Chi_G_db='cohsex.ndb.QP';gap=2.6_SP
     call Chi_G_load(x,k)
   case(3)
     Chi_G_mode='DYSON';gap=2.6_SP
     allocate(Chi_G_energy(2,2,1,1),Chi_G_weight(2,2,1,1),Chi_G_occupation(2,2,1,1))
     Chi_G_energy(1,:,:,:)=-1.2_DP;Chi_G_energy(2,:,:,:)=1.4_DP
     Chi_G_weight=1._DP
     Chi_G_occupation(1,:,:,:)=1._DP;Chi_G_occupation(2,:,:,:)=0._DP
   end select
   call Chi_G_bubble(2,x,k,w,d,.FALSE.,response)
   do iw=1,w%n_freqs
     factor=1._SP/(w%p(iw)-gap)-1._SP/(w%p(iw)+gap)
     expected=factor*reshape([1.04_SP,.4_SP,.4_SP,1.04_SP],[2,2])
     if(maxval(abs(response(:,:,iw)-expected))>tolerance)then
       print *,'FAIL: MPI ',trim(Chi_G_mode),' rank ',myid
       call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     endif
   enddo
   call Chi_G_free()
 enddo
 if(master_cpu)print *,'PASS: G0/COHSEX/Dyson production bubble with real MPI, ranks=',ncpu
 call MPI_FINALIZE(ierr)
end program
