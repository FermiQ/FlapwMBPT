      subroutine tos_x_solid(x,ttos,tdos)
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: x
	  real*8, intent(inout) :: ttos,tdos
	  integer :: me1,ispin,k,i,n,ii
	  real*8 :: const,st,sd,de,ffd,dffd,fermi_dirac,ee,f
	  me1=me+1
	  const=2.d0/dfloat(nspin*nrel)
      do ispin=1,nspin
        do k=1,npnt
	      n=n_bnd(k,ispin)
          st=0.d0
          sd=0.d0
	      do i=1,n
            de=e_bnd(i,k,ispin)-x
            ffd=fermi_dirac(de)
            dffd=betta_t*ffd*(1.d0-ffd)
            if(k_integral=='O1') then
              ffd=0.d0
              dffd=0.d0
              do ii=1,nn_k_int
                ee=e_bnd(i,k,ispin)
     &            +de_bnd(i,k,ispin,1)*k_int_vec(1,ii)
     &            +de_bnd(i,k,ispin,2)*k_int_vec(2,ii)
     &            +de_bnd(i,k,ispin,3)*k_int_vec(3,ii)
                f=fermi_dirac(ee-x)
                ffd=ffd+f
                dffd=dffd+betta_t*f*(1.d0-f)
              enddo
              ffd=ffd/nn_k_int
              dffd=dffd/nn_k_int
            endif
            st=st+ffd
            sd=sd+dffd
          enddo
          ttos=ttos+wgt(k)*st*const
          tdos=tdos+wgt(k)*sd*const
        enddo  !!! over k
      enddo  !!! over ispin
      if(maswrk.and.iter<=iter_dft+iter_hf) 
     &   write(iun,'(3(1x,f13.8))')x,ttos,tdos
      end
