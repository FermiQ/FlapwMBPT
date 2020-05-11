      subroutine write_vrt_scale(key)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      integer :: i_len
      i_len=len_trim(tmpfile)
      open(8,file=tmpfile(1:i_len)//'_vrt_scale',form='unformatted')
      if(key==0) write(8)vrt_scale
      if(key==1) read(8)vrt_scale
      close(8)
      end
