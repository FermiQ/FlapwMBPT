      subroutine read_boltzman(key,n,nel_b,ebnd,e_don,s_don)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: key
      integer, intent(inout) :: n
      real*8, intent(out) :: ebnd(n,npnt,nspin),nel_b,e_don,s_don
      character*3 :: b_meth
      integer :: k,i,i_len,ispin
      real*8 :: v(3)
      b_meth='_HF'
      if(maswrk) then
        i_len=len_trim(allfile)
        if(b_meth=='DFT') then
          open(3,file=allfile(1:i_len)//'_DFT.energy')
        else if(b_meth=='_HF') then
          open(3,file=allfile(1:i_len)//'_HF.energy')
        else if(b_meth=='_QP') then
          open(3,file=allfile(1:i_len)//'_scQSGW.energy')
        endif
        read(3,*)
        read(3,*)npnt,nspin,chem_pot,nel_b,e_don,s_don
        if(key==0) then
          read(3,'(3(1x,e17.10),i5)')(v(i),i=1,3),n
        else
          do k=1,npnt
            read(3,'(3(1x,e17.10),i5)')(v(i),i=1,3),n
            do i=1,n
              read(3,'(2(1x,e17.10))')(ebnd(i,k,ispin),ispin=1,nspin)
            enddo
          enddo
        endif
        close(3)
      endif
      if(goparr) then
#ifdef MPI
        call brdcst(msgint,nel_b,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,chem_pot,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,e_don,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,s_don,8,master,MPI_COMM_WORLD)
        if(key/=0) then
          call brdcst(msgdbl,ebnd,8*n*npnt*nspin,master,MPI_COMM_WORLD)
        endif
#endif
      endif
      end
