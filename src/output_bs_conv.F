      subroutine output_bs_conv(conv_table,iter_lad)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_lad
      real*8, intent(in) :: conv_table(0:iter_lad,0:7,nc_line)
      integer :: k,iq,it,i,i_nu
      if(maswrk) then
        do iq=1,nc_line
          k=kline_in_npnt_c(iq)
          write(iun,'(a34,3f12.5)')
     &	      ' Convergency of Ladder eqn. for Q= ',
     &	              (pnt_c(i,k),i=1,3)
          write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
          do it=0,iter_lad
            write(iun,'(i3,8e9.2)')it,
     & 	        (conv_table(it,i_nu,iq),i_nu=0,7)
          enddo
        enddo
      endif
      end
      
      
      
      
      subroutine output_bs_conv_1q(conv_table,iter_lad)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_lad
      real*8, intent(in) :: conv_table(0:iter_lad,0:7)
      integer :: it,i,i_nu
      if(maswrk) then
        write(iun,'(a34,3f12.5)')' Convergency of Ladder eqn. for Q= ',
     &	                       (pnt_c(i,q_suscept),i=1,3)
        write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
        do it=0,iter_lad
          write(iun,'(i3,8e9.2)')it,
     & 	        (conv_table(it,i_nu),i_nu=0,7)
        enddo
      endif
      end