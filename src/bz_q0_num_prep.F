      subroutine bz_q0_num_prep
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i,n_k,ii,j,iq,id,k
      real*8 :: pi2a,v(3),d,d0
      pi2a=(pi+pi)/par
      if(istruc<=3) then
c ------------------------ Cubic -----------------------------------
        n_phead_dir=1
        allocate(khead(n_phead_pnt,n_phead_dir))
        d0=1.d5
        do iq=2,npnt
          d=dot_product(pnt(:,iq),pnt(:,iq))
          if(d<d0) then
            d0=d
            khead(1,1)=iq
          endif
        enddo
      else if((istruc>3.and.istruc<8).or.istruc==14) then
c  --------- Tetragonal, Hexagonal or Rhombohedral -----------------
        n_phead_dir=2
        allocate(khead(n_phead_pnt,n_phead_dir))
c ----- First point: shortest in XY plane --------------------------
        d0=1.d5
        do iq=2,npnt
          if(abs(pnt(3,iq))>1.d-5) cycle
          d=dot_product(pnt(:,iq),pnt(:,iq))
          if(d<d0) then
            d0=d
            khead(1,1)=iq
          endif
        enddo
c ----- First point: shortest in Z direction -----------------------
        d0=1.d5
        do iq=2,npnt
          d=pnt(1,iq)**2+pnt(2,iq)**2
          if(d>1.d-7) cycle
          d=dot_product(pnt(:,iq),pnt(:,iq))
          if(d<d0) then
            d0=d
            khead(1,2)=iq
          endif
        enddo
      else if(istruc==8) then
c ---------------- Simple Orthorhombic -----------------------------
        n_phead_dir=3
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,3
          if(i==1) v=qb0(:,1)
          if(i==2) v=qb0(:,2)
          if(i==3) v=qb0(:,3)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      else if(istruc==9) then
c ---------------- Body Centered Orthorhombic ----------------------
        n_phead_dir=3
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,3
          if(i==1) v=-qb0(:,1)+qb0(:,2)+qb0(:,3)
          if(i==2) v=qb0(:,1)-qb0(:,2)+qb0(:,3)
          if(i==3) v=qb0(:,1)+qb0(:,2)-qb0(:,3)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      else if(istruc==10) then
c ---------------- Face Centered Orthorhombic ----------------------
        n_phead_dir=3
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,3
          if(i==1) v=qb0(:,2)+qb0(:,3)
          if(i==2) v=qb0(:,1)+qb0(:,3)
          if(i==3) v=qb0(:,1)+qb0(:,2)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      else if(istruc==11) then
c ---------------- Base Centered Orthorhombic ----------------------
        n_phead_dir=3
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,3
          if(i==1) v=-qb0(:,1)+qb0(:,2)
          if(i==2) v=qb0(:,1)+qb0(:,2)
          if(i==3) v=qb0(:,3)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      else if(istruc==12) then
c --------------- Simple monoclinic --------------------------------
c -------     C axis is unique along Z; A axis is along X ----------
        n_phead_dir=4
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,4
          if(i==1) v=qb0(:,1)
          if(i==2) v=qb0(:,1)-qb0(:,2)
          if(i==3) v=qb0(:,2)
          if(i==4) v=qb0(:,3)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      else if(istruc==121) then
c ------------------ Simple monoclinic -----------------------------
c -------------    a,b<c;alpha<90;beta,gamma=90 --------------------
        n_phead_dir=4
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,4
          if(i==1) v=qb0(:,1)
          if(i==2) v=qb0(:,2)
          if(i==3) v=qb0(:,2)+qb0(:,3)
          if(i==4) v=qb0(:,3)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      else
        n_phead_dir=6
        allocate(khead(n_phead_pnt,n_phead_dir))
        do i=1,6
          if(i==1) v=qb0(:,1)
          if(i==2) v=qb0(:,2)
          if(i==3) v=qb0(:,3)
          if(i==4) v=qb0(:,1)+qb0(:,2)
          if(i==5) v=qb0(:,1)+qb0(:,3)
          if(i==6) v=qb0(:,2)+qb0(:,3)
          call zone1_number(v,rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          khead(1,i)=iq
        enddo
      endif
c ------- Add second point in each direction if necessary ----------
      if(n_phead_pnt==2) then
        do i=1,n_phead_dir
          v=2.d0*pnt(:,khead(1,i))
          call zone1_number(v,rb0,ndiv,iq)  !! in 1-zone
          iq=index_k1(iq)    !! in 0-zone
          khead(2,i)=i_kref(iq)
        enddo
      endif
c ------------------------------------------------------------------
      allocate(khead0(n_phead_pnt*n_phead_dir+1))
      n_k=1
      khead0(1)=1
      do id=1,n_phead_dir
        do i=1,n_phead_pnt
          ii=0
          do j=1,n_k
            if(khead0(j)==khead(i,id)) ii=1
          enddo
          if(ii==0) then
            n_k=n_k+1
            khead0(n_k)=khead(i,id)
          endif
        enddo
      enddo
      n_k_head_0=n_k
      do id=1,n_phead_dir
        do i=1,n_phead_pnt
          do j=1,n_k_head_0
            if(khead(i,id)==khead0(j)) khead(i,id)=j
          enddo
        enddo
      enddo
      allocate(w_head_dir(n_phead_dir))
      do id=1,n_phead_dir
        k=khead(1,id)
        iq=khead0(k)
        w_head_dir=dfloat(ngroup)/num_sym_k0(iq)
      enddo
      d=sum(w_head_dir)
      w_head_dir=w_head_dir/d
      end
