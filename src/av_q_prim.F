      subroutine av_q_prim
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,i,j,nn
      real*8 :: pi2a,pi8,pi2,vol0,qs,qb(3,3),qq0(3),qq(3),qqs,
     &          rad_inscribed_sphere
      pi2=pi+pi
      pi2a=pi2/par
      pi8=8.d0*pi
      vol0=pi2**3/amega/nqdiv
      qs=(3.d0*vol0/4.d0/pi)**(1.d0/3.d0)
      q0aver=4.d0*pi*qs**3/3.d0/vol0
      q1aver=pi2*qs*qs/vol0
      q2aver=4.d0*pi*qs/vol0
      q0aver_c=q0aver
      q1aver_c=q1aver
      q2aver_c=q2aver
      logq_aver=log(qs)-1.d0/3.d0
      if(mesh_k_c.and.nqdiv_c/=nqdiv) then
        vol0=pi2**3/amega/nqdiv_c
        qs=(3.d0*vol0/4.d0/pi)**(1.d0/3.d0)
        q0aver_c=4.d0*pi*qs**3/3.d0/vol0
        q1aver_c=pi2*qs*qs/vol0
        q2aver_c=4.d0*pi*qs/vol0
        logq_aver_c=log(qs)-1.d0/3.d0
      endif
      return
c ----- Average over parallelepiped --------------------------------
      nn=30
      qs=rad_inscribed_sphere(qb0)*pi2a
      q1aver=pi2*qs*qs
      q2aver=4.d0*pi*qs
      qb=qb0/nn
      qq0=0.5d0*(qb0(:,1)+qb0(:,2)+qb0(:,3))
      do k=1,nn
        do j=1,nn
          do i=1,nn
            qq=pi2a*(qb(:,1)*(i-0.5d0)+qb(:,2)*(j-0.5d0)
     &              +qb(:,3)*(k-0.5d0)-qq0)
            qqs=sqrt(dot_product(qq,qq))
            if(qqs>qs) then
              q1aver=q1aver+vol0/qqs/nn**3
              q2aver=q2aver+vol0/qqs/qqs/nn**3
            endif
          enddo
        enddo
      enddo
      q1aver=q1aver/vol0
      q2aver=q2aver/vol0
      if(mesh_k_c.and.nqdiv_c/=nqdiv) then
        vol0=pi2**3/amega/nqdiv_c
        qs=rad_inscribed_sphere(qb0_c)*pi2a
        q1aver_c=pi2*qs*qs
        q2aver_c=4.d0*pi*qs
        qb=qb0_c/nn
        qq0=0.5d0*(qb0_c(:,1)+qb0_c(:,2)+qb0_c(:,3))
        do k=1,nn
          do j=1,nn
            do i=1,nn
              qq=pi2a*(qb(:,1)*(i-0.5d0)+qb(:,2)*(j-0.5d0)
     &                +qb(:,3)*(k-0.5d0)-qq0)
              qqs=sqrt(dot_product(qq,qq))
              if(qqs>qs) then
                q1aver_c=q1aver_c+vol0/qqs/nn**3
                q2aver_c=q2aver_c+vol0/qqs/qqs/nn**3
              endif
            enddo
          enddo
        enddo
        q1aver_c=q1aver_c/vol0
        q2aver_c=q2aver_c/vol0
      endif
      end
