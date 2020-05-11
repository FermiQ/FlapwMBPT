      subroutine output_xi_rax_mt(ch,f,wrx,nrax)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      character*2, intent(in) :: ch
      integer, intent(in) :: nrax
      real*8, intent(in) :: wrx(0:nrax)
      complex*16, intent(in) :: f(0:nrax,nsort,2)
      integer :: i_len,n,i,nk,isort
      real*8 :: omeg,ev,hh,t,hh1
      real*8, allocatable :: sm(:,:,:),sm1(:,:)
      complex*16, allocatable :: f1(:,:,:)
      ev=evolt/2.d0
      i_len=len_trim(allfile)
      nk=min(3,nsort)
      allocate(sm(0:nrax,nsort,2))
      allocate(sm1(0:nrax,nsort))
      allocate(f1(0:nrax,nsort,2))
      sm=0.d0
      sm1=0.d0
      f1=(0.d0,0.d0)
      do isort=1,nsort
        do i=nrax,1,-1
c --------------- Thermal factor ---------------------------------
          t=wrx(i)
          t=0.5d0*betta_t*t
          t=1.d0/tanh(t)
          f1(i,isort,:)=t*f(i,isort,:)
        enddo
        f1(0,isort,:)=f1(1,isort,:)
        do i=1,nrax
          hh=wrx(i)-wrx(i-1)
c --------------- Thermal factor ---------------------------------
          t=(wrx(i)+wrx(i-1))/2
          t=0.5d0*betta_t*t
          t=1.d0/tanh(t)
          hh1=hh*t
c ----------------------------------------------------------------
          if(ch(1:1)=='Z') then
            sm(i,isort,1)=sm(i-1,isort,1)
     &                   +hh*imag(f(i-1,isort,1)+f(i,isort,1))/2/pi
            sm1(i,isort)=sm1(i-1,isort)
     &                  +hh1*imag(f(i-1,isort,1)+f(i,isort,1))/2/pi
          else if(ch(1:1)=='+') then
            sm(i,isort,1)=sm(i-1,isort,1)
     &                   +hh*(imag(f(i-1,isort,1)+f(i,isort,1))
     &                      +imag(f(i-1,isort,2)+f(i,isort,2)))/4/pi
            sm(i,isort,2)=sm(i-1,isort,2)
     &                   +hh*(imag(f(i-1,isort,1)+f(i,isort,1))
     &                      -imag(f(i-1,isort,2)+f(i,isort,2)))/8/pi
            sm1(i,isort)=sm1(i-1,isort)
     &                  +hh1*(imag(f(i-1,isort,1)+f(i,isort,1))
     &                       +imag(f(i-1,isort,2)+f(i,isort,2)))/4/pi
          endif
        enddo
      enddo
      if(maswrk) then
        if(ch(1:1)=='Z') then
          open(3,file=allfile(1:i_len)//'_'//ch//'MT_dif_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(imag(f(i,n,1)/ev/pi),
     &                           imag(f1(i,n,1)/ev/pi),n=1,nk)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_'//ch//'MT_sum_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(sm(i,n,1),sm1(i,n),n=1,nk)
          enddo
          close(3)
        else if(ch=='+-') then
          open(3,file=allfile(1:i_len)//'_XX+YY'//'MT_dif_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(imag(f(i,n,1)+f(i,n,2))/2/ev/pi,
     &                           imag(f1(i,n,1)+f1(i,n,2))/2/ev/pi,
     &                           n=1,nk)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_XX+YY'//'MT_sum_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(sm(i,n,1),sm1(i,n),n=1,nk)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_Z'//'MT_dif_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(imag(f(i,n,1)-f(i,n,2))/4/ev/pi,n=1,nk)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_Z'//'MT_sum_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(sm(i,n,2),n=1,nk)
          enddo
          close(3)
        else if(ch=='+0') then
          open(3,file=allfile(1:i_len)//'_XX0YY'//'MT_dif_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(imag(f(i,n,1)+f(i,n,2))/2/ev/pi,
     &                           imag(f1(i,n,1)+f1(i,n,2))/2/ev/pi,
     &                           n=1,nk)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_XX0YY'//'MT_sum_dft.rax')
          do i=0,nrax
            omeg=wrx(i)
            write(3,100)omeg*ev,(sm(i,n,1),sm1(i,n),n=1,nk)
          enddo
          close(3)
        endif
      endif
      deallocate(sm,f1,sm1)
100   format(f11.5,1x,6(e11.4,1x))
      end
