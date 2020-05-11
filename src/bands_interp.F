      subroutine bands_interp(ebnd)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: ebnd(nbndf,npnt,nspin)
      integer :: ispin,k,ind_k,info,i,ind,m,lmax,lmin,nn,nb,n,ns,nstr
      real*8 :: pi2,phase,c1,c2,xx,ee
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: rhom(:)
      complex*16, allocatable :: smk(:,:),x(:,:),a(:,:),sm(:),b(:,:)
      pi2=pi+pi
      nstr=nstar_r
      allocate(smk(2:nstr,npnt))
      smk=(0.d0,0.d0)
      do k=1,npnt
        do m=2,nstr
          do i=indstar_r(m-1)+1,indstar_r(m)
            phase=pi2*dot_product(rbs(:,i),pnt(:,k))
            smk(m,k)=smk(m,k)+dcmplx(cos(phase),sin(phase))
          enddo
          smk(m,k)=smk(m,k)/(indstar_r(m)-indstar_r(m-1))
        enddo
      enddo
c ------------------------------------------------------------------
      allocate(rhom(2:nstr))
      c1=0.75d0
      c2=0.75d0
      do i=2,nstr
        xx=(rcostar_r(i)/rcostar_r(2))**2
        rhom(i)=(1.d0-c1*xx)**2+c2*xx**3
      enddo
c ------------------------------------------------------------------
      do k=1,npnt-1
        smk(:,k)=smk(:,k)-smk(:,npnt)
      enddo
      n=npnt-1
      allocate(x(2:nstr,n))
      do i=2,nstr
        x(i,:)=smk(i,1:n)/rhom(i)
      enddo
      allocate(a(n,n))
      call zgemm('c','n',n,n,nstr-1,(1.d0,0.d0),x(2,1),nstr-1,
     &           smk(2,1),nstr-1,(0.d0,0.d0),a,n)
      deallocate(x)
c ------------------------------------------------------------------
      lmin=10000
      lmax=0
      do ispin=1,nspin
        do k=1,npnt
          nb=n_bnd(k,ispin)
          do i=1,nb
            ee=ebnd(i,k,ispin)-chem_pot
            if(ee>emindos) lmin=min(lmin,i)
            if(ee<emaxdos) lmax=max(lmax,i)
          enddo
        enddo
      enddo
      nn=lmax-lmin+1
      allocate(b(n,nn*nspin))
      do ispin=1,nspin
        do i=1,nn
          ind=(ispin-1)*nn+i
          do k=1,n
            b(k,ind)=ebnd(lmin-1+i,k,ispin)-ebnd(lmin-1+i,npnt,ispin)
          enddo
        enddo
      enddo
      allocate(ipiv(n))
      call zgesv(n,nn*nspin,a,n,ipiv,b,n,info)
      deallocate(ipiv,a)
      ns=nstr-1
      allocate(a(nstr,nn*nspin))
      call zgemm('n','n',ns,nn*nspin,n,(1.d0,0.d0),smk(2,1),ns,b,n,
     &           (0.d0,0.d0),a(2,1),nstr)
      deallocate(b)
      do i=2,nstr
        a(i,:)=conjg(a(i,:))/rhom(i)
      enddo
      deallocate(rhom)
      do ispin=1,nspin
        do i=1,nn
          ind=(ispin-1)*nn+i
          a(1,ind)=ebnd(lmin-1+i,npnt,ispin)
          do m=2,nstr
            a(1,ind)=a(1,ind)-a(m,ind)*smk(m,npnt)
          enddo
        enddo
      enddo
c --------------------- Check --------------------------------------
c      do k=1,n
c        xx=0.d0
c        do m=2,nstr
c          xx=xx+a(m,1)*smk(m,k)
c        enddo
c        print *,k,ebnd(lmin,npnt,1)+xx,ebnd(lmin,k,1)
c      enddo
c ------------------------------------------------------------------
      deallocate(smk)
      allocate(e_plot(nn,num_k_plot,nspin))
      e_plot=0.d0
      allocate(sm(2:nstr))
      do ispin=1,nspin
        do ind_k=1,ndim_kplot(me_k+1)
          k=n_mpi_kplot(me_k+1)+ind_k
          sm=(0.d0,0.d0)
          do m=2,nstr
            do i=indstar_r(m-1)+1,indstar_r(m)
              phase=pi2*dot_product(rbs(:,i),k_plot(:,k))
              sm(m)=sm(m)+dcmplx(cos(phase),sin(phase))
            enddo
            sm(m)=sm(m)/(indstar_r(m)-indstar_r(m-1))
          enddo
          do i=1,nn
            ind=(ispin-1)*nn+i
            e_plot(i,k,ispin)=a(1,ind)
            do m=2,nstr
              e_plot(i,k,ispin)=e_plot(i,k,ispin)+a(m,ind)*sm(m)
            enddo
          enddo
        enddo
      enddo
      deallocate(a,sm)
      k=nn*num_k_plot*nspin
      if(nproc_k/=1) call DGOP(e_plot,k,'  +',comm_k)
      call bands_gnuplot(nn,nn,e_plot,1)
      deallocate(e_plot)
      end
