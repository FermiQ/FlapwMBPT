      subroutine small_star_q(iq,nq0,indk0,wk0)
c     Finds irreducible k-points with respect to
c     the symmetry operations which left q-point intact
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq
      integer, intent(out) :: nq0,indk0(nqdiv)
      real*8, intent(out) :: wk0(nqdiv)
      integer :: k,ig,k0,ig0
      integer, allocatable :: nst(:),flag(:)
      allocate(nst(nqdiv))
      allocate(flag(nqdiv))
      nq0=0
      nst=0
      flag=0
      do k=1,nqdiv
        if(flag(k)==1) cycle
        nq0=nq0+1
        nst(nq0)=1
        indk0(nq0)=k
        flag(k)=1
        do ig0=1,num_sym_k(iq)
          ig=list_sym_k(ig0,iq)
          k0=ip_k(k,ig)
          if(flag(k0)==1) cycle
          flag(k0)=1
          nst(nq0)=nst(nq0)+1
        enddo
      enddo
      do k=1,nq0
        wk0(k)=dfloat(nst(k))/nqdiv
      enddo
      deallocate(nst,flag)
      end
