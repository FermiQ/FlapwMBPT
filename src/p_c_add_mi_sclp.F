      subroutine p_c_add_mi_sclp(pwsclp,s1,npw,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw,k
      complex*16, intent(in) :: s1(n_pbmt,npw)
      complex*16, intent(inout) :: pwsclp(nrow_pw_max,ncol_pw_max)
      integer :: j,i,ii,jj
      do j=1,ncol_pw(k)
        jj=refcol_pw(j,k)-n_pbmt
        if(jj>0) then
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)
            if(ii<=n_pbmt) then
              pwsclp(i,j)=pwsclp(i,j)+s1(ii,jj)
            endif
          enddo
        endif
        jj=refcol_pw(j,k)
        if(jj<=n_pbmt) then
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)-n_pbmt
            if(ii>0) then
              pwsclp(i,j)=pwsclp(i,j)+conjg(s1(jj,ii))
            endif
          enddo
        endif
      enddo
      end



      subroutine p_r_add_mi_sclp(pwsclp,s1,npw,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw,k
      real*8, intent(in) :: s1(n_pbmt,npw)
      real*8, intent(inout) :: pwsclp(nrow_pw_max,ncol_pw_max)
      integer :: j,i,ii,jj
      do j=1,ncol_pw(k)
        jj=refcol_pw(j,k)-n_pbmt
        if(jj>0) then
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)
            if(ii<=n_pbmt) then
              pwsclp(i,j)=pwsclp(i,j)+s1(ii,jj)
            endif
          enddo
        endif
        jj=refcol_pw(j,k)
        if(jj<=n_pbmt) then
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)-n_pbmt
            if(ii>0) then
              pwsclp(i,j)=pwsclp(i,j)+s1(jj,ii)
            endif
          enddo
        endif
      enddo
      end

