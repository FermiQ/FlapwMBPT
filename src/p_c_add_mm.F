      subroutine p_c_add_mm(pmm,tc,n1pbm,n2pbm,ind0,jnd0,n0,n1)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1pbm,n2pbm,ind0,jnd0,n0,n1
      complex*16, intent(in) :: tc(n1pbm,n2pbm)
      complex*16, intent(inout) :: pmm(ndim_pbmt(me_b+1))
      integer :: j,km1,m,i,km,ind,ij
      do j=1,n2pbm
        km1=jnd0+j
        m=(km1-1)*km1/2
        do i=1,n1pbm
          km=ind0+i
          if(km>km1) cycle
          ind=m+km
          if(ind<=n0) cycle
          if(ind>n1) cycle
          ij=ind-n0
          pmm(ij)=pmm(ij)+tc(i,j)
        enddo
      enddo
      end


      subroutine p_c_add_mm_sclp(pwsclp,tc,n1pbm,n2pbm,ind0,jnd0,k,key)
c     key = 0 for right upper triangle
c     key = 1 for left lower triangle
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1pbm,n2pbm,ind0,jnd0,k,key
      complex*16, intent(in) :: tc(n1pbm,n2pbm)
      complex*16, intent(inout) :: pwsclp(nrow_pw_max,ncol_pw_max)
      integer :: j,i,ii,jj
      if(key==0) then
        do j=1,ncol_pw(k)
          jj=refcol_pw(j,k)-jnd0
          if(jj<1.or.jj>n2pbm) cycle
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)-ind0
            if(ii<1.or.ii>n1pbm) cycle
            pwsclp(i,j)=pwsclp(i,j)+tc(ii,jj)
          enddo
        enddo
      endif
      if(key==1) then
        do j=1,ncol_pw(k)
          jj=refcol_pw(j,k)-ind0
          if(jj<1.or.jj>n1pbm) cycle
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)-jnd0
            if(ii<1.or.ii>n2pbm) cycle
            pwsclp(i,j)=pwsclp(i,j)+conjg(tc(jj,ii))
          enddo
        enddo
      endif
      end
