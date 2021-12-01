Module mod_spectra_heom
!! Reproducing figure 5 of Strumpfer Schulten, JCTC 8, 2808 (2012)
!! Chen, Zheng, Shi, Tan, JCP 131, 094502 (2009)
implicit none
real*8, parameter :: clight=137.d0,av=6.0221367D23,hbar=1.d0!1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.d0!1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: wave_to_J=1.98644568d-23
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2s=2.418884326505d-17
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi

!! HEOM variables
integer LL,KK
integer nn_tot
integer,allocatable :: nn(:,:,:)
integer,allocatable :: map_nplus(:,:,:),map_nneg(:,:,:),map_sum(:)
integer,allocatable :: zero(:)
complex*16,allocatable :: rho(:,:,:),rho_t0(:,:)
complex*16,allocatable :: c_matsubara(:),omg_matsubara(:)
real*8 tolerance
real*8 sum_c_over_nu

!! Output
complex*16,allocatable :: dip_mom_corr(:)

!! System specific
integer nquant
real*8 gamma,eta,temperature,lambda
complex*16,allocatable :: Hamil_e(:,:),dip_moment(:,:)
real*8,allocatable :: Q_op(:,:,:)

!! Evolution
integer nsteps,nsteps_filter
real*8 dt,tim_tot,curr_time
integer flag_scale,flag_truncate,flag_spectra

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i
  real*8,allocatable:: Ham_site(:,:),dip_real(:,:)

  pi=dacos(-1.d0)
  open(10,file="spectra_heom.inp")
  read(10,*) nquant
  read(10,*) LL
  read(10,*) KK
  read(10,*) tolerance
  read(10,*) dt
  read(10,*) tim_tot
  read(10,*) nsteps_filter
  read(10,*) flag_scale
  read(10,*) flag_truncate
  read(10,*) flag_spectra
  read(10,*) gamma
  read(10,*) lambda
  read(10,*) temperature
  allocate(Hamil_e(nquant,nquant),Ham_site(nquant,nquant))
  read(10,*)
  read(10,*)
  do i=1,nquant
    read(10,*) Ham_site(i,:)
  enddo
  allocate(Q_op(nquant,nquant,nquant))
  read(10,*)
  read(10,*)
  do i=1,nquant
    read(10,*) Q_op(:,:,i)
    read(10,*)
  enddo
  allocate(dip_moment(nquant,nquant),dip_real(nquant,nquant))
  if(flaG_spectra==1) then
    read(10,*)
    read(10,*)dip_real
    dip_moment=dip_real
  endif
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  
  if(flag_scale.ne.0.and.flag_scale.ne.1) then
    write(6,*) "flag_scale should be either 1 or 0. terminating ..."
    stop
  endif

  if(flag_truncate.ne.0.and.flag_truncate.ne.1) then
    write(6,*) "flag_truncate should be either 1 or 0. terminating ..."
    stop
  endif

  if(flag_spectra.ne.0.and.flag_spectra.ne.1) then
    write(6,*) "flag_spectra should be either 1 or 0. terminating ..."
    stop
  endif
  !---------------------------------------------------------- 

  Ham_site=Ham_site*wave_to_J/au2J
  dt=dt/au2s
  tim_tot=tim_tot/au2s
  gamma=1.d15/gamma*au2s
  eta=2*lambda/hbar*wave_to_J/au2J
!  eta=2*lambda/hbar
  temperature=temperature*1.38064852d-23/au2J
  Hamil_e=Ham_site

  !nn_tot=(factorial(LL+nquant*KK)/factorial(LL))/factorial(nquant*KK)
  call compute_nn_tot
  write(6,*) "total number of auxillary density matrices = ",nn_tot
  nsteps=nint(tim_tot/dt)

  allocate(dip_mom_corr(nsteps))
  allocate(nn(nquant,0:KK,nn_tot),map_nplus(nquant,0:KK,nn_tot),map_nneg(nquant,0:KK,nn_tot))
  allocate(map_sum(0:LL),zero(nn_tot))
  allocate(rho(nquant,nquant,nn_tot),rho_t0(nquant,nquant))
  allocate(c_matsubara(0:KK),omg_matsubara(0:KK))

  call compute_nn
  call compute_map

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer clck_counts_beg, clck_counts_end, clck_rate
  real*8 t1,t2

  call system_clock ( clck_counts_beg, clck_rate )
  call cpu_time(t1)

  open(100,file="rho.out")
  if(flag_spectra==1)open(101,file="abs_spectra.out")
  call setup_parameters
  call init_cond
  call evolve(nsteps)
  if(flag_spectra==1)call fft_dip_mom_corr

  call system_clock ( clck_counts_end, clck_rate )
  call cpu_time(t2)
  write(6,*) "total wall clock time=",(clck_counts_end - clck_counts_beg) / real(clck_rate)
  write(6,*) "total cpu time=",t2-t1

  close(100)
  if(flag_spectra==1)close(101)

end subroutine main
!---------------------------------------------------------- 

subroutine setup_parameters
  implicit none
  integer i,j,k
  real*8 tmp,omg,cc,eps

  !JJ=1.0!*wave_to_J
  omg_matsubara(0)=gamma 
  c_matsubara(0)=eta*gamma/2.d0* (1.d0/tan(hbar*gamma/(2.d0*kb*temperature))-iota)
  do k=1,KK
    omg_matsubara(k)=2*k*pi*kb*temperature/hbar
    c_matsubara(k)=2*eta*gamma*kb*temperature/hbar * omg_matsubara(k)/(omg_matsubara(k)**2-gamma**2)
  enddo

  sum_c_over_nu=0.d0
  do k=KK+1,KK+200
    omg=2*k*pi*kb*temperature/hbar
    cc=2*eta*gamma*kb*temperature/hbar * omg/(omg**2-gamma**2)
    sum_c_over_nu=sum_c_over_nu+cc/omg
  enddo
  !write(6,*) sum_c_over_nu,eta*kb*temperature/gamma-real(sum(c_matsubara/omg_matsubara))
  !stop

  eps=0.d0
  do i=1,nquant
    eps=eps+Hamil_e(i,i)
  enddo
  eps=eps/real(nquant)
  do i=1,nquant
    Hamil_e(i,i)=Hamil_e(i,i)-eps
  enddo

  zero=0

!  Hamil_e=Hamil_e*wave_to_J/au2J

write(6,*) "Parameters Set ..."

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i

  rho=0.d0
  rho(1,1,1)=1.d0

  if(flag_spectra==1) rho(:,:,1)=matmul(dip_moment,rho(:,:,1))

  rho_t0=rho(:,:,1)

  curr_time=0.d0
  dip_mom_corr=0.d0

write(6,*) "Intitial Conditions set ... "

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in)::nsteps
  integer i

  do i=1,nsteps
    if(flag_spectra==1)call compute_dip_mom_corr(i)
    call write_output(i)
    if(mod(i,nsteps_filter)==0)call filter
    call rk4(rho,dt)
    curr_time=curr_time+dt
  enddo

end subroutine evolve
!-----------------------------------------------------------------  

subroutine write_output(i)
  implicit none
  integer,intent(in):: i
  integer j

  if(flag_spectra==0)write(100,'(f15.7$)') curr_time*au2s*1.d15
  do j=1,nquant
    if(flag_spectra==0)write(100,'(f15.7$)')real(rho(j,j,1))
  enddo
  write(100,*)

end subroutine write_output
!-----------------------------------------------------------------  

subroutine compute_dip_mom_corr(j)
  implicit none
  integer,intent(in) :: j
  integer i
  complex*16 mat(nquant,nquant)

  mat=matmul(dip_moment,rho(:,:,1))

  do i=1,nquant
    dip_mom_corr(j)=dip_mom_corr(j)+mat(i,i)
  enddo

end subroutine compute_dip_mom_corr
!-----------------------------------------------------------------  

subroutine rk4(rho,dt)
  implicit none
  complex*16,intent(inout)::rho(nquant,nquant,nn_tot)
  real*8,intent(in)::dt
  complex*16,dimension(nquant,nquant,nn_tot) :: k1,k2,k3,k4

  call compute_deriv(rho,k1)
  call compute_deriv(rho+0.5*dt*k1,k2)
  call compute_deriv(rho+0.5*dt*k2,k3)
  call compute_deriv(rho+dt*k3,k4)

  rho=rho+dt/6.d0*(k1+2*k2+2*k3+k4)

end subroutine rk4
!-----------------------------------------------------------------  

subroutine compute_deriv(rho,drho_dt)
  implicit none
  complex*16,intent(in)::rho(nquant,nquant,nn_tot)
  complex*16,intent(out)::drho_dt(nquant,nquant,nn_tot)
  complex*16 tmp(nquant,nquant)
  integer n
  integer tid,OMP_GET_THREAD_NUM

  call omp_set_num_threads(4)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,tmp)
  !$OMP DO SCHEDULE(STATIC)
  do n=1,nn_tot
    call compute_deriv_n(n,rho,tmp)
    drho_dt(:,:,n)=tmp
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

end subroutine compute_deriv
!-----------------------------------------------------------------  

subroutine compute_deriv_n(n,rho,drho_dt)
  !! implicit none
  !! Eq. 15, JCP 131, 094502 (2009)
  !! Scaled HEOM approach from JCP 130, 084105 (2009)
  integer,intent(in) :: n
  complex*16,intent(in)::rho(nquant,nquant,nn_tot)
  complex*16,intent(out)::drho_dt(nquant,nquant)
  complex*16 tmp(nquant,nquant),mat_tmp(nquant,nquant)
  integer m,k,nplus,nminus
  integer nvec(nquant,0:KK)!,nvec_plus(nquant,KK),nvec_neg(nquant,KK)

    m=n
    call index(m,nvec,0)
    tmp=-iota*commute(Hamil_e,rho(:,:,n))/hbar

    if(zero(n)==0) then !! matrix at n is not filtered out

      do m=1,nquant
        do k=0,KK
          tmp=tmp-nvec(m,k)*omg_matsubara(k)*rho(:,:,n)
        enddo
      enddo

      !! Ihizaki-Tanimura scheme for truncation
      !! JPSJ 74 3131, 2005
      if(flag_truncate==1) then
        do m=1,nquant
          !mat_tmp=0.d0;mat_tmp(m,m)=1.d0 !! matrix |m><m|
          mat_tmp=Q_op(:,:,m)
          !tmp=tmp-(eta*kb*temperature/gamma-real(sum(c_matsubara/omg_matsubara)))*commute(mat_tmp,commute(mat_tmp,rho(:,:,n)))
          tmp=tmp-sum_c_over_nu*commute(mat_tmp,commute(mat_tmp,rho(:,:,n)))
        enddo
      endif
    endif

    do m=1,nquant
      !mat_tmp=0.d0;mat_tmp(m,m)=1.d0 !! matrix |m><m|
      mat_tmp=Q_op(:,:,m)
      do k=0,KK
        nplus=map_nplus(m,k,n)

        if(nplus>0.and.nplus<=nn_tot) then
          !! Scaled HEOM approach from JCP 130, 084105 (2009)
          if(zero(nplus)==0.and.flag_scale==1) tmp=tmp - iota*commute(mat_tmp,rho(:,:,nplus))*&
          sqrt((nvec(m,k)+1.d0)*abs(c_matsubara(k)))
          if(zero(nplus)==0.and.flag_scale==0) tmp=tmp - iota*commute(mat_tmp,rho(:,:,nplus))
        endif
      enddo
    enddo

    do m=1,nquant
      !mat_tmp=0.d0;mat_tmp(m,m)=1.d0 !! matrix |m><m|
      mat_tmp=Q_op(:,:,m)
      do k=0,KK
        nminus=map_nneg(m,k,n)
        if(nminus>0.and.nminus<=nn_tot) then
          if(zero(nminus)==0)then
            !! Scaled HEOM approach from JCP 130, 084105 (2009)
            if(flag_scale==1)tmp=tmp-iota*(c_matsubara(k)*matmul(mat_tmp,rho(:,:,nminus)) &
- dconjg(c_matsubara(k))*matmul(rho(:,:,nminus),mat_tmp) )*sqrt(nvec(m,k)/abs(c_matsubara(k)))
            if(flag_scale==0)tmp=tmp-iota*nvec(m,k)*(c_matsubara(k)*matmul(mat_tmp,rho(:,:,nminus)) &
- dconjg(c_matsubara(k))*matmul(rho(:,:,nminus),mat_tmp) )

          endif
        endif
      enddo
    enddo

    drho_dt=tmp
    !stop

end subroutine compute_deriv_n
!-----------------------------------------------------------------  

subroutine filter
  !! Shi, Chen, Nan, Xu, Yan, JCP 130, 084105 (2009)
  integer n

  zero=0
  do n=1,nn_tot
    if(maxval(abs(rho(:,:,n)))<tolerance) then
      rho(:,:,n)=0.d0
      zero(n)=1
    endif
  enddo

!write(6,*) curr_time,nn_tot-sum(zero),nn_tot

end subroutine filter
!-----------------------------------------------------------------  

subroutine index(n,nvec,iflag)
  implicit none
  integer,intent(inout)::n
  integer,intent(inout)::nvec(nquant,0:KK)
  integer,intent(in)::iflag
  integer m,k,ind,temp
  integer n_beg,n_end,l_sum
  integer clck_counts_beg,clck_counts_end, clck_rate
  real*8 t1,t2

  if(iflag==0) then  !n is input, nvec is output
    nvec=nn(:,:,n)
  endif

  if(iflag==1) then  !nvec is input, n is output
    n=0
    l_sum=sum(nvec)
    if(l_sum<=LL) then
      n_beg=map_sum(l_sum)
      if(l_sum==LL) n_end=nn_tot
      if(l_sum<LL) n_end=map_sum(l_sum+1)-1
      do m=1,nn_tot!n_beg,n_end
        if(all(nvec==nn(:,:,m))) then
          n=m
          exit
        endif
      enddo
    endif
  endif
        
end subroutine index
!-----------------------------------------------------------------  

subroutine binary(n,nvec,iflag,base)
  implicit none
  integer,intent(inout)::n
  integer,intent(inout)::nvec(nquant,0:KK)
  integer,intent(in)::iflag,base
  integer m,k,ind,temp

  if(iflag==0) then  !n is input, nvec is output
    temp = n-1
    nvec=0
    m=1;k=0
    DO WHILE(temp > 0)
      nvec(m,k)=mod(temp,base)
      k=k+1
      if(k>KK) then
        k=1;m=m+1
      endif
      temp = temp / base
    END DO
  endif

  if(iflag==1) then  !nvec is input, n is output
    n=0
    ind=0
    do m=1,nquant
      do k=0,KK
        if(nvec(m,k).ne.0)n=n+nvec(m,k)*base**ind
        ind=ind+1
      enddo
    enddo
    n=n+1
  endif

end subroutine binary
!-----------------------------------------------------------------  

subroutine compute_nn
  implicit none
  integer n,m,n_tot,tot
  integer nvec(nquant,0:KK)
  integer n_beg,n_end

!! This method works, but inefficient for large LL
!  nvec=0;nvec(nquant,KK)=LL
!  call binary(n_tot,nvec,1,LL+1)
!
!  tot=0
!  do n=1,n_tot
!    m=n
!    call binary(m,nvec,0,LL+1)
!    if(sum(nvec)<=LL) then
!      tot=tot+1
!      nn(:,:,tot)=nvec
!    endif
!  enddo

!! Alternate recursive method

  n_tot=0
  n_beg=0;n_end=0
  write(6,*) "hierarchy level, start index, end index"
  do n=0,LL
    call compute_nn_sum_L(n,n_beg,n_end)
    map_sum(n)=n_beg
    write(6,*) n,n_beg,n_end
  enddo
  !write(6,*) n_end,nn_tot

end subroutine compute_nn
!-----------------------------------------------------------------  

subroutine compute_nn_sum_L(L,n_beg,n_end)
  implicit none
  integer,intent(in)::L
  integer,intent(inout)::n_beg,n_end!! at input, state/end point of entries for sum L-1, at output, for entries for sum L
  integer i,j,m,k,tot_L
  integer flag_new,nvec(nquant,0:KK)

  if(L==0) then
    nn(:,:,1)=0
    n_beg=1;n_end=1
  else
    tot_L=n_end
    do i=n_beg,n_end
      do m=1,nquant
        do k=0,KK
          
          nvec=nn(:,:,i);nvec(m,k)=nvec(m,k)+1
          flag_new=0
          if(tot_L>n_end) then
            do j=n_end+1,tot_L
              if(all(nvec==nn(:,:,j))) then
                flag_new=1;exit
              endif
            enddo
          endif

          if(flag_new==0) then
            tot_L=tot_L+1
            nn(:,:,tot_L)=nn(:,:,i)
            nn(m,k,tot_L)=nn(m,k,i)+1
!            map_nneg(m,k,tot_L)=i
!            map_nplus(m,k,i)=tot_L
          endif
        enddo
      enddo
    enddo
    n_beg=n_end+1
    n_end=tot_L
  endif

end subroutine compute_nn_sum_L
!-----------------------------------------------------------------  

subroutine compute_map
  implicit none
  integer n,m,k
  integer nvec(nquant,0:KK),nvec_plus(nquant,0:KK),nvec_neg(nquant,0:KK)
  integer nplus,nneg

  map_nneg=0
  map_nplus=0
  do n=1,nn_tot
    m=n
    call index(m,nvec,0)
    do m=1,nquant
      do k=0,KK
        nvec_plus=nvec;nvec_plus(m,k)=nvec_plus(m,k)+1
        call index(nplus,nvec_plus,1)
        map_nplus(m,k,n)=nplus
       
        if(nvec(m,k)>0) then
          nvec_neg=nvec;nvec_neg(m,k)=nvec_neg(m,k)-1
          call index(nneg,nvec_neg,1)
          map_nneg(m,k,n)=nneg
        endif
      enddo
    enddo
  enddo

write(6,*) "map completed ..."

end subroutine compute_map
!-----------------------------------------------------------------  

subroutine compute_nn_tot
  !nn_tot=(factorial(LL+nquant*KK)/factorial(LL))/factorial(nquant*KK)
  implicit none
  integer i
  real*8 tmp

  tmp=1.d0
  do i=1,nquant*(KK+1)
    tmp=tmp*(LL+i)/real(i)
  enddo
  nn_tot=nint(tmp)

end subroutine compute_nn_tot
!-----------------------------------------------------------------  

function commute(mat1,mat2) result(mat3)
  complex*16,intent(in) :: mat1(:,:),mat2(:,:)
  complex*16,allocatable :: mat3(:,:)
  integer i1

  i1=size(mat1,1)
  allocate(mat3(i1,i1))

  mat3=matmul(mat1,mat2)-matmul(mat2,mat1)

end function commute
!-----------------------------------------------------------------  

integer function factorial(n)
  implicit none
  integer,intent(in) :: n
  integer i

  if(n==0) then
    factorial=1
  else
    factorial=1
    do i=1,n
      factorial=factorial*i
    enddo
  endif

end function factorial
!-----------------------------------------------------------------  

subroutine fft_dip_mom_corr
  implicit none
  integer i,j,k
  real*8 w1,tim(nsteps),dat_r(2*nsteps),dat(4*nsteps),tim_doub(2*nsteps)
  complex*16 dat_c(nsteps),tmp(3)
  real*8,allocatable :: freq(:),dat_fft(:)
  integer nst,m_nst,k1,k2
  real*8 cos_damp(nsteps),cc,ss,re,im
  real*8 ww(nsteps)
  complex*16 ft(nsteps)

  ft=0.d0
  do i=1,nsteps
    ww(i)=-5.d0+10.d0*i/real(nsteps)
    tim(i)=(i-1)*dt
    writE(100,*)tim(i),real(dip_mom_corr(i)),imag(dip_mom_corr(i))
  enddo

  ft=0.d0
  do i=1,nsteps
    do j=1,nsteps
      ft(i)=ft(i)+exp(iota*ww(i)*tim(j))*dip_mom_corr(j)
    enddo
    write(101,*)ww(i),real(ft(i))*dt
  enddo

!    k=0
!    do j=1,nsteps
!      k=k+1
!      !read(100,*)tim(k),dat_r(2*k-1),dat_r(2*k)  !! real part, complex part
!      tim(k)=dt*(j-1)
!      dat_r(2*k-1)=real(dip_mom_corr(k))
!      dat_r(2*k)=imag(dip_mom_corr(k))
!    enddo
!    
!    nst=2*nsteps
!    
!    do j=1,nsteps
!      k1=nsteps+j
!      k2=nsteps-j+1
!      dat(nst+2*j-1)=(dat_r(2*j-1))!*dcos(pi*tim(j)/(2*tim(nsteps)))
!      dat(nst+2*j)=(dat_r(2*j))!*dcos(pi*tim(j)/(2*tim(nsteps)))
!      tim_doub(k1)=tim(j)
!      tim_doub(k2)=-tim(j)
!    enddo
!    
!    do j=1,nst
!      k1=2*j-1;k2=2*(nst-j+1)-1
!      dat(k1)=dat(k2)
!      k1=2*j;k2=2*(nst-j+1)
!      dat(k1)=-dat(k2)
!    enddo
!    
!    do j=1,nst
!      write(102,*)tim_doub(j),dat(2*j-1),dat(2*j)
!    enddo
!    
!    !tim_doub=tim_doub*1.d-15
!    call fft(tim_doub,dat,2*nst,freq,dat_fft,m_nst)
!    do j=m_nst/4+1,m_nst/2
!      cc=dcos(freq(j)*tim(nsteps))
!      ss=-dsin(freq(j)*tim(nsteps))
!      re=dat_fft(2*j-1);im=dat_fft(2*j)
!      write(101,'(4es15.5)')freq(j),re,im!re*cc-im*ss,im*cc+re*ss
!    enddo
!    do j=1,m_nst/4
!      cc=dcos(freq(j)*tim(nsteps))
!      ss=-dsin(freq(j)*tim(nsteps))
!      re=dat_fft(2*j-1);im=dat_fft(2*j)
!      write(101,'(4es15.5)')freq(j),re,im!re*cc-im*ss,im*cc+re*ss
!    enddo
!    write(101,*)
!
end subroutine fft_dip_mom_corr
!-----------------------------------------------------------------  
!
!subroutine fft(t,dat,ndat,w,dat_f,m)
!  !! computes fft of dat[]
!  !! tim in fs; freq in cm-1
!  implicit none
!  integer,intent(in)::ndat
!  real*8,intent(in)::t(ndat/2),dat(ndat)
!  integer,intent(out)::m
!  real*8,intent(out),allocatable::w(:),dat_f(:)
!  real*8 delt
!  complex*16 ff,tmp
!  integer i,j
!
!  m=1
!  do while(m<ndat)
!   m=m*2
!   if(m>=ndat) exit
!  enddo
!
!write(6,*) m
!
!  allocate(dat_f(m),w(m/2))
!  dat_f=0.d0
!  do i=1,ndat
!    dat_f(i)=dat(i)
!  enddo
!
!  delt = t(2)-t(1)
!
!  do i=1,m/4
!   w(i) = (i-1)/(m/2*delt)
!  enddo
!  do i=m/4+1,m/2
!   w(i) = -(m/2-i+1)/(m/2*delt)
!  enddo
!  !w = w/clight
!  w=w*2*pi
!
!write(6,*)
!
!do i=1,m/2
!  write(103,*) dat_f(2*i-1),dat_f(2*i)
!enddo
!
!  !call fft_gateway(dat_f(:),m)
!  dat_f=0.d0
!  do i=1,m/2
!    tmp=0.d0
!    do j=1,ndat/2
!      ff=dat(2*j-1)+iota*dat(2*j)
!      tmp=tmp+exp(iota*w(i)*t(j))*ff
!    enddo
!    dat_f(2*i-1)=real(tmp)
!    dat_f(2*i)=imag(tmp)
!  enddo
!  dat_f = dat_f*delt
!
!end subroutine fft
!!-----------------------------------------------------------------  

End Module mod_spectra_heom
