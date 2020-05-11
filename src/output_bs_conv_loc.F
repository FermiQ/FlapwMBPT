      subroutine output_bs_conv_loc(conv_table,iter_lad)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_lad
      real*8, intent(in) :: conv_table(0:iter_lad,0:7,natom)
      integer :: it,i_nu,iatom
      if(maswrk) then
        do iatom=1,natom
          write(iun,'(a40,i5)')
     &        ' Convergency of Ladder eqn. for IATOM = ',iatom
          write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
          do it=0,iter_lad
            write(iun,'(i3,8e9.2)')it,
     &          (conv_table(it,i_nu,iatom),i_nu=0,7)
          enddo
        enddo
      endif
      end
