      subroutine sym_sig_k_from_k0_mm(s,k,ig)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: ig,k
	  complex*16, intent(inout) :: s(nfun,nfun)
	  integer :: ia,ib,iatom,isort,i
	  real*8 :: pi2,tt(3),phase
      complex*16, allocatable :: tmp(:,:),tr(:)
	  pi2=pi+pi
      allocate(tmp(nfun,nfun))
      allocate(tr(natom))
	  do iatom=1,natom
	    tt=tshift(:,iatom,ig)
	    phase=pi2*dot_product(pnt(:,k),tt)
	    tr(iatom)=dcmplx(cos(phase),sin(phase))
	  enddo
	  tmp=(0.d0,0.d0)
	  do i=1,nfun
c ------------ Index from the right  ------------------------------
	    call sym_val_left(ig,s(i,:),tmp(:,i),2)
	    do iatom=1,natom
	      isort=is(iatom)
	      ia=io_lem(iatom)
	      ib=ia+lfunm(isort)-1
	      tmp(ia:ib,i)=tmp(ia:ib,i)*conjg(tr(iatom))
	    enddo
	  enddo
	  tmp=transpose(tmp)
	  do i=1,nfun
c ------------ Index from the left ------------------------------
        call sym_val_left(ig,tmp(:,i),s(:,i),1)
	    do iatom=1,natom
	      isort=is(iatom)
	      ia=io_lem(iatom)
	      ib=ia+lfunm(isort)-1
	      s(ia:ib,i)=s(ia:ib,i)*tr(iatom)
	    enddo
	  enddo
	  deallocate(tmp,tr)
      end
