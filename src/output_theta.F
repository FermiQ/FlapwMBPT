      subroutine output_theta(n,th,ch)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      character*3, intent(in) :: ch
      integer, intent(in) :: n
      complex*16, intent(in) :: th(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,i10,i20,i30,i40,i_len,i,j,k,l,i1,
     &           j1,k1,l1
	i_len=len_trim(allfile)
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_Theta_loc_'//ch)
        do isp4=1,2
          do isp3=1,2
            do isp2=1,2
              do isp1=1,2
                i40=n*(isp4-1)
                i30=n*(isp3-1)
                i20=n*(isp2-1)
                i10=n*(isp1-1)
                do l=1,n
                  l1=i40+l
                  do k=1,n
                    k1=i30+k
                    do j=1,n
                      j1=i20+j
                      do i=1,n
                        i1=i10+i
                        write(3,'(8i3,1x,e12.5,1x,e12.5)')
     &                    i,j,k,l,isp1,isp2,isp3,isp4,th(i1,j1,k1,l1)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
	  enddo
	  close(3)	
	endif   
      end