      subroutine write_pw_q_tau_mm(key,gmm)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(inout) :: gmm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                             ndim3_k(me_k+1))
      integer :: ich,i_len
      character*8 :: cha
c ------------------------------------------------------------
      i_len=len_trim(tmpfile)
      cha='        '
      ich=1000000*me_b+1000*me_t+me_k
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:1),
     &       form='unformatted')
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:2),
     &       form='unformatted')
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:3),
     &       form='unformatted')
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:4),
     &       form='unformatted')
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:5),
     &       form='unformatted')
      else if(ich<1000000) then
        write(cha(1:6),'(i6)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:6),
     &       form='unformatted')
      else if(ich<10000000) then
        write(cha(1:7),'(i7)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:7),
     &       form='unformatted')
      else if(ich<100000000) then
        write(cha(1:8),'(i8)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mm'//cha(1:8),
     &       form='unformatted')
      endif
      if(key==0) write(8)gmm
      if(key==1) read(8)gmm
      close(8)
      end




      subroutine write_pw_q_tau_mi(key,gmi)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(inout) :: gmi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                             ndim3_tn,ndim3_k(me_k+1))
      integer :: ich,i_len
      character*8 :: cha
c ------------------------------------------------------------
      i_len=len_trim(tmpfile)
      cha='        '
      ich=1000000*me_b+1000*me_t+me_k
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:1),
     &       form='unformatted')
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:2),
     &       form='unformatted')
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:3),
     &       form='unformatted')
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:4),
     &       form='unformatted')
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:5),
     &       form='unformatted')
      else if(ich<1000000) then
        write(cha(1:6),'(i6)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:6),
     &       form='unformatted')
      else if(ich<10000000) then
        write(cha(1:7),'(i7)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:7),
     &       form='unformatted')
      else if(ich<100000000) then
        write(cha(1:8),'(i8)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_mi'//cha(1:8),
     &       form='unformatted')
      endif
      if(key==0) write(8)gmi
      if(key==1) read(8)gmi
      close(8)
      end




      subroutine write_pw_q_tau_ii(key,gii)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(inout) :: gii(ncmpl*nd_b_pbint(me_b+1),ndim3_tn,
     &                             ndim3_k(me_k+1))
      integer :: ich,i_len
      character*8 :: cha
c ------------------------------------------------------------
      i_len=len_trim(tmpfile)
      cha='        '
      ich=1000000*me_b+1000*me_t+me_k
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:1),
     &       form='unformatted')
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:2),
     &       form='unformatted')
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:3),
     &       form='unformatted')
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:4),
     &       form='unformatted')
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:5),
     &       form='unformatted')
      else if(ich<1000000) then
        write(cha(1:6),'(i6)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:6),
     &       form='unformatted')
      else if(ich<10000000) then
        write(cha(1:7),'(i7)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:7),
     &       form='unformatted')
      else if(ich<100000000) then
        write(cha(1:8),'(i8)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_q_tau_ii'//cha(1:8),
     &       form='unformatted')
      endif
      if(key==0) write(8)gii
      if(key==1) read(8)gii
      close(8)
      end




      subroutine write_ws_head(key)
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
      open(8,file=tmpfile(1:i_len)//'_ws_head',form='unformatted')
      if(ncmpl==1) then
        if(key==0) write(8)ws_head,ws_head_nu,ws01_met_r
        if(key==1) read(8)ws_head,ws_head_nu,ws01_met_r
      else
        if(key==0) write(8)ws_head,ws_head_nu,ws01_met
        if(key==1) read(8)ws_head,ws_head_nu,ws01_met
      endif
      close(8)
      end




      subroutine write_w11_red(key)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      integer :: i_len,ich
      character*5 :: cha
      i_len=len_trim(tmpfile)
      cha='        '
      ich=me_t
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
        open(8,file=tmpfile(1:i_len)//'_w11_red'//cha(1:1),
     &       form='unformatted')
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
        open(8,file=tmpfile(1:i_len)//'_w11_red'//cha(1:2),
     &       form='unformatted')
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
        open(8,file=tmpfile(1:i_len)//'_w11_red'//cha(1:3),
     &       form='unformatted')
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
        open(8,file=tmpfile(1:i_len)//'_w11_red'//cha(1:4),
     &       form='unformatted')
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
        open(8,file=tmpfile(1:i_len)//'_w11_red'//cha(1:5),
     &       form='unformatted')
      endif
      if(key==0) write(8)w11_red
      if(key==1) read(8)w11_red
      close(8)
      end




      subroutine write_ws_head_par(key)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      integer :: i_len,ich
      character*5 :: cha
      i_len=len_trim(tmpfile)
      cha='        '
      ich=me_t
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
        open(8,file=tmpfile(1:i_len)//'_ws_par'//cha(1:1),
     &       form='unformatted')
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
        open(8,file=tmpfile(1:i_len)//'_ws_par'//cha(1:2),
     &       form='unformatted')
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
        open(8,file=tmpfile(1:i_len)//'_ws_par'//cha(1:3),
     &       form='unformatted')
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
        open(8,file=tmpfile(1:i_len)//'_ws_par'//cha(1:4),
     &       form='unformatted')
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
        open(8,file=tmpfile(1:i_len)//'_ws_par'//cha(1:5),
     &       form='unformatted')
      endif
      if(ncmpl==1) then
        if(key==0) write(8)ws_head_tens_r,ws_wing_r
        if(key==1) read(8)ws_head_tens_r,ws_wing_r
      else
        if(key==0) write(8)ws_head_tens,ws_wing
        if(key==1) read(8)ws_head_tens,ws_wing
      endif
      close(8)
      end


      subroutine write_pw_sclp(key)
c     key = 0 - Write
c     key = 1 - Read
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      integer :: ich,i_len
      character*8 :: cha
c ------------------------------------------------------------
      i_len=len_trim(tmpfile)
      cha='        '
      ich=1000000*me_b+1000*me_t+me_k
      if(ich<10) then
        write(cha(1:1),'(i1)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:1),
     &       form='unformatted')
      else if(ich<100) then
        write(cha(1:2),'(i2)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:2),
     &       form='unformatted')
      else if(ich<1000) then
        write(cha(1:3),'(i3)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:3),
     &       form='unformatted')
      else if(ich<10000) then
        write(cha(1:4),'(i4)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:4),
     &       form='unformatted')
      else if(ich<100000) then
        write(cha(1:5),'(i5)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:5),
     &       form='unformatted')
      else if(ich<1000000) then
        write(cha(1:6),'(i6)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:6),
     &       form='unformatted')
      else if(ich<10000000) then
        write(cha(1:7),'(i7)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:7),
     &       form='unformatted')
      else if(ich<100000000) then
        write(cha(1:8),'(i8)')ich
        open(8,file=tmpfile(1:i_len)//'_pw_sclp'//cha(1:8),
     &       form='unformatted')
      endif
      if(key==0) write(8)pw_sclp
      if(key==1) read(8)pw_sclp
      close(8)
      end

