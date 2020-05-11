      subroutine optics_loc_fld(nrax,wrx,e,mode)
 	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
      character*1, intent(in) :: mode
      integer, intent(in) :: nrax
      complex*16, intent(in) :: wrx(0:nrax)
      complex*16, intent(out) :: e(3,3,0:nrax)
      integer :: n,npw,ne,i
      real*8, allocatable :: ve(:)
      complex*16, allocatable :: vie(:,:),b11(:,:,:),b12(:,:,:),
     &                           b21(:,:,:),b22(:,:,:),
     &                           thet_m(:,:,:),thet(:,:),tmp(:,:)
      npw=nplwgw(1)
      n=n_pbmt+npw
      allocate(ve(n))
      allocate(vie(n,n))
      call v_spectrum(1,n,ve,vie,ne)
c --------------------------------------------------------------------
      allocate(b11(0:nrax,3,3))
      allocate(b21(ne,0:nrax,3))
      allocate(b12(ne,0:nrax,3))
      allocate(b22(ne,ne,0:nrax))
      call opt_chi0_q0(b11,b12,b21,b22,ne,nrax,n,vie,wrx)
      allocate(thet(ne,ne))
      thet=(0.d0,0.d0)
      if(mode/='1') then  !! if it is not RPA
        call resp_lda_prepare('100',.false.)
        allocate(thet_m(n,n,n_ixc_0))
        call theta_q_ss(1,n,thet_m)
        allocate(tmp(n,ne))
        call zgemm('n','n',n,ne,n,(1.d0,0.d0),thet_m,n,vie,n,
     &             (0.d0,0.d0),tmp,n)
        deallocate(thet_m)
        call zgemm('c','n',ne,ne,n,(1.d0,0.d0),vie,n,tmp,n,(0.d0,0.d0),
     &             thet,ne)
        deallocate(tmp)
      endif
      do i=1,ne
        thet(i,i)=thet(i,i)+ve(i)
      enddo
      call opt_enh(b11,b12,b21,b22,ne,nrax,thet,e)
      deallocate(ve,vie,b11,b12,b21,b22)
      end
