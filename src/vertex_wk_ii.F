      subroutine vertex_wk_ii(wii,vx,nrr3,indrr3)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: nrr3,
     &                       indrr3(nr_full_red,nr_full_red,nqdiv_c)
	  real*8, intent(in) :: wii(nrr3)
      complex*16, intent(inout) :: vx(nrel*nr_full_red,nrel*nr_full_red,
     &                                nqdiv_c)
      integer :: i,j,ir,irl,jrl,ii,jj
      real*8 :: wi
      do ir=1,nqdiv_c
        do j=1,nr_full_red
          do i=1,nr_full_red
            ii=indrr3(i,j,ir)
            wi=wii(ii)
            if(irel/=2) then
              vx(i,j,ir)=wi*vx(i,j,ir)
            else
              do jrl=1,2
                jj=(jrl-1)*nr_full_red+j
                do irl=1,2
                  ii=(irl-1)*nr_full_red+i
                  vx(ii,jj,ir)=wi*vx(ii,jj,ir)
                enddo
              enddo
            endif
          enddo   !! over j
        enddo   !! over ir
      enddo   !! over i
      end

