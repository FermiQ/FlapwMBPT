      subroutine hovb_lucb(eb,ho,lu)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      real*8, intent(in) :: eb(nbndf,npnt,nspin)
      real*8, intent(out) :: ho,lu
      integer :: ispin,k,n,nocc
      integer, allocatable :: iref(:)
      real*8, allocatable :: e(:)
      nocc=idnint(nelec)
      if(irel/=2.and.nspin==1) nocc=nocc/2
      ho=-1.d5
      lu=1.d5
      do k=1,npnt
        n=n_bnd(k,1)
        if(nspin==2) n=n+n_bnd(k,2)
        allocate(e(n))
        allocate(iref(n))
        ii=0
        do ispin=1,nspin
          do i=1,n_bnd(k,ispin)
            ii=ii+1
            e(ii)=eb(i,k,ispin)
          enddo
        enddo
        call ordering_all(n,e,iref)
        ho=max(ho,e(iref(nocc)))
        lu=min(lu,e(iref(nocc+1)))
        deallocate(e,iref)
      enddo
      end
