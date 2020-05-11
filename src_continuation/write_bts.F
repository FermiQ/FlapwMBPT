      subroutine write_bts(key)
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
      character*5 :: cha
c ------------------------------------------------------------
      cha='     '
      if(me<10) then
        write(cha(1:1),'(i1)')me
      else if(me<100) then
        write(cha(1:2),'(i2)')me
      else if(me<1000) then
        write(cha(1:3),'(i3)')me
      else if(me<10000) then
        write(cha(1:4),'(i4)')me
      else if(me<100000) then
        write(cha(1:5),'(i5)')me
      endif
      i_len=len_trim(tmpfile)
      open(8,file=tmpfile(1:i_len)//'_bts'//cha,form='unformatted')
      if(key==0) write(8)fxc_bts
      if(key==1) read(8)fxc_bts
      close(8)
      end
