      subroutine tos_lda(x,ttos,tdos)
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      use vertex_mod
      use units_mod
      implicit none
	  integer :: ispin,k,ibnd,n,ii
	  real*8 :: const,ttos,tdos,st,sd,ffd,fermi_dirac,x,ffd1,f,ee
	  const=2.d0/dfloat(nspin*nrel)
      ttos=0.d0
      tdos=0.d0
      do ispin=1,nspin
        do k=1,npnt
	      n=n_bnd(k,ispin)
          st=0.d0
          sd=0.d0
          do ibnd=1,n
            ffd=fermi_dirac(e_bnd(ibnd,k,ispin)-x)
            ffd1=betta_t*ffd*(1.d0-ffd)
            if(k_integral=='O1') then
              ffd=0.d0
              ffd1=0.d0
              do ii=1,nn_k_int
                ee=e_bnd(ibnd,k,ispin)
     &            +de_bnd(ibnd,k,ispin,1)*k_int_vec(1,ii)
     &            +de_bnd(ibnd,k,ispin,2)*k_int_vec(2,ii)
     &            +de_bnd(ibnd,k,ispin,3)*k_int_vec(3,ii)
                f=fermi_dirac(ee-x)
                ffd=ffd+f
                ffd1=ffd1+betta_t*f*(1.d0-f)
              enddo
              ffd=ffd/nn_k_int
              ffd1=ffd1/nn_k_int
            endif
            st=st+ffd
            sd=sd+ffd1
          enddo
          ttos=ttos+wgt(k)*st*const
          tdos=tdos+wgt(k)*sd*const
        enddo  !!! over k
      enddo  !!! over ispin
      if(maswrk) write(iun,'(3(1x,f13.8))')x,ttos,tdos
      end
