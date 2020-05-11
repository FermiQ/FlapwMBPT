      subroutine output_map_mpi
      use manager_mod
      use parallel_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: i_len,i,ii
      integer, allocatable :: map(:,:)
      allocate(map(0:nproc-1,6))
      map=0
      map(me,1)=me_t
      map(me,2)=me_k
      map(me,3)=me_b
      map(me,4)=me_t_k
      map(me,5)=me_t_b
      map(me,6)=me_k_b
#ifdef MPI
      if(nproc/=1) call igop(map,6*nproc,'  +',MPI_COMM_WORLD)
#endif
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_MPI_MAP')
        write(3,'(7a9)')'       ME','     ME_T','     ME_K','     ME_B',
     &                  '   me_t_b','   ME_T_B','   ME_K_B'
        do i=0,nproc-1
          write(3,'(7i9)')i,(map(i,ii),ii=1,6)
        enddo
      endif
      deallocate(map)
      end
