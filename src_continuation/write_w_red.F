      subroutine write_w_red(key,g)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(inout) :: g(n_pbtot_red,n_pbtot_red,ndim3_nu,
     & 	                       ndim3_k(me_k+1))
      integer :: ich,i_len
      character*5 :: cha
c ------------------------------------------------------------
      cha='     '
      ich=1000*me_t+me_k
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
      endif
      i_len=len_trim(tmpfile)
      open(8,file=tmpfile(1:i_len)//'_wq'//cha,form='unformatted')
      if(key==0) write(8)g
      if(key==1) read(8)g
      close(8)
      end
