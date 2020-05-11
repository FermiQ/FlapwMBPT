      integer function zone1_c_int(iv1,iv2,iv3)
c     Returns kb-point given its integer (in terms of qb0_c) coordinates
c     in the extended first BZ
      use models_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iv1,iv2,iv3
      integer :: i1,i2,i3,ind
      i1=mod(iv1,ndiv_c(1))
      i2=mod(iv2,ndiv_c(2))
      i3=mod(iv3,ndiv_c(3))
      ind=ndiv_c(1)*ndiv_c(2)*i3+ndiv_c(1)*i2+i1+1
      zone1_c_int=index_k1_c(ind)
      end
