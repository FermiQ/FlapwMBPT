      subroutine write_sig_x(key,key1,g)
c     key = 0 - Write
c     key = 1 - Read
c     key1 = 1 - V_H_BB_OLD
c     key1 = 2 - SIGX_SOLID_K_OLD
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
	  integer, intent(in) :: key,key1
      complex*16, intent(inout) :: g(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),
     &                               nspin)
	  integer :: ich,i_len
	  character*5 :: cha
      i_len=len_trim(tmpfile)
	  cha='     '
	  ich=me
	  if(ich<10) then
	    write(cha(1:1),'(i1)')ich
        if(key1==1) then
          open(8,file=tmpfile(1:i_len)//'_vh'//cha(1:1),
     &         form='unformatted')
        else if(key1==2) then
          open(8,file=tmpfile(1:i_len)//'_sx'//cha(1:1),
     &         form='unformatted')
        endif
	  else if(ich<100) then
	    write(cha(1:2),'(i2)')ich
        if(key1==1) then
          open(8,file=tmpfile(1:i_len)//'_vh'//cha(1:2),
     &         form='unformatted')
        else if(key1==2) then
          open(8,file=tmpfile(1:i_len)//'_sx'//cha(1:2),
     &         form='unformatted')
        endif
	  else if(ich<1000) then
	    write(cha(1:3),'(i3)')ich
        if(key1==1) then
          open(8,file=tmpfile(1:i_len)//'_vh'//cha(1:3),
     &         form='unformatted')
        else if(key1==2) then
          open(8,file=tmpfile(1:i_len)//'_sx'//cha(1:3),
     &         form='unformatted')
        endif
	  else if(ich<10000) then
	    write(cha(1:4),'(i4)')ich
        if(key1==1) then
          open(8,file=tmpfile(1:i_len)//'_vh'//cha(1:4),
     &         form='unformatted')
        else if(key1==2) then
          open(8,file=tmpfile(1:i_len)//'_sx'//cha(1:4),
     &         form='unformatted')
        endif
	  else if(ich<100000) then
	    write(cha(1:5),'(i5)')ich
        if(key1==1) then
          open(8,file=tmpfile(1:i_len)//'_vh'//cha(1:5),
     &         form='unformatted')
        else if(key1==2) then
          open(8,file=tmpfile(1:i_len)//'_sx'//cha(1:5),
     &         form='unformatted')
        endif
	  endif
      if(key==0) write(8)g
	  if(key==1) read(8)g
	  close(8)
      end
