      subroutine str_con_matrix(s,nlma,nlmb,hl,khbr)
c     makes lattice structure constants. make hl
c     in rcnke0 but note that there dr is in alat units, q in tpiba.
c     exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
c     functions are real hankels and bessels without extra factors. E=0
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer, intent(in) :: nlma,nlmb,khbr
      integer :: lmax,l,mlm,lm,lget,klm,lk,ii,icg,llm,ll0
      real*8 :: p32,fac2
      complex*16 :: s(nlma,nlmb),hl(*),sm
      real*8, allocatable :: sig(:)
      p32=32.d0*pi*pi
      lmax=lget(nlma)+lget(nlmb)
      allocate(sig(0:lmax))
      sig(0)=1.d0
      do l=1,lmax
        sig(l)=-sig(l-1)
      enddo
c ---------------------------------------
      do mlm=1,nlma
        lm=lget(mlm)
        do klm=1,nlmb
          lk=lget(klm)
          sm=(0.d0,0.d0)
          ii=max0(mlm,klm)
          ii=(ii*(ii-1))/2+min0(mlm,klm)
          do icg=indxcg(ii),indxcg(ii+1)-1
            llm=jcg(icg)
            ll0=lget(llm)
            if(khbr==0.and.ll0.eq.lm+lk) sm=sm+cg(icg)*hl(llm)
            if(khbr==1.and.mod(lm+lk-ll0,2)==0) sm=sm+cg(icg)*hl(llm)
     &                                         *omega_hse**(lm+lk-ll0)
          enddo
          if(khbr==0) s(mlm,klm)=p32*sm*sig(lk)/fac2(lm+lm+1)
     &                                         /fac2(lk+lk+1)
          if(khbr==1) s(mlm,klm)=p32*sm*sig(lm)
        enddo
      enddo
      deallocate(sig)
      end
