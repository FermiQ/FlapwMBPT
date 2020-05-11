      subroutine memory_solid(maxbv,maxbj,maxw2)
	  use atom_mod
	  use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
      integer, intent(in) :: maxbv,maxbj,maxw2
	  integer :: m21,nlmb
	  nlmb=nrel*maxb+1
      m21=max0(maxb,maxpb)
	  call memory_atom
      allocate(sym(limlpb,nsort))
	  sym=.false.
      allocate(nsym(nsort))
      allocate(lmsym(limlpb,nsort))
	  lmsym=0
      if(magn.eq.2) then
        allocate(symb(limlpb,3,nsort))
	    symb=.false.
        allocate(nsymb(nsort))
        allocate(lmsymb(limlpb,3,nsort))
	    lmsymb=0
      endif
      allocate(eny(maxntle,nlmb,nsort,nspin_0))
      allocate(augm_coef(3,nlmb,nsort,nspin_0))
      allocate(fi(maxntle,nlmb,nsort,nspin))
      allocate(dfi(maxntle,nlmb,nsort,nspin))
      allocate(dny(maxntle,nlmb,nsort,nspin))
      allocate(fidot(maxntle,nlmb,nsort,nspin))
      allocate(dfidot(maxntle,nlmb,nsort,nspin))
      allocate(fidot2(maxntle,nlmb,nsort,nspin))
      allocate(dfidot2(maxntle,nlmb,nsort,nspin))
      allocate(ffsmt(2,2,maxntle,maxntle,nlmb,nsort,nspin))
      allocate(ffhmt(2,2,maxntle,maxntle,nlmb,nsort,nspin))
      allocate(a0val(maxntle,nlmb,nsort,nspin_0))
      a0val=1.d0
      allocate(a0cor(ncormax,nsort,nspin_0))
      a0cor=1.d0
	  allocate(vxc0(nsort))
	  allocate(vxc1(nsort))
c	 *****************************************************
	  allocate(bett((maxylm+1)**2))
      call cnsbet(maxylm,bett)
      allocate(rcostar(maxplwro/2))
      allocate(indstar(0:maxplwro/2))
      allocate(igbs(3,maxplw))
      allocate(igbsf(maxplwro))
      allocate(gbs(3,maxplw))
      end
