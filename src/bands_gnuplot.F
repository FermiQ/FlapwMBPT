      subroutine bands_gnuplot(nbn,nbn0,epl,mode)
c     mode = 0 regular bands (now only for DFT)
c     mode = 1 interpolated bands
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: nbn,nbn0,mode
      real*8, intent(in) :: epl(nbn0,num_k_plot,nspin)
      character*1 :: ap
	  integer :: i_len,ispin,k,i
	  character*72 :: f_bnd
	  character*3, allocatable :: glabel(:)
	  real*8, allocatable :: e(:)
	  allocate(e(2:n_k_div))
	  if(.not.maswrk) return
	  ap="'"
	  i_len=len_trim(allfile)
	  if(ubi=='dft') then
	    if(mode==0) f_bnd=allfile(1:i_len)//'_dft_band_LAPW.dat'
        if(mode==1) f_bnd=allfile(1:i_len)//'_dft_band_INTRP.dat'
	  else if(ubi==' hf') then
	    f_bnd=allfile(1:i_len)//'_hf_band.dat'
	  else if(ubi==' qp') then
	    f_bnd=allfile(1:i_len)//'_qp_band.dat'
	  else if(ubi==' gw') then
	    f_bnd=allfile(1:i_len)//'_gw_band.dat'
	  else if(ubi=='psi') then
	    f_bnd=allfile(1:i_len)//'_psi_band.dat'
	  else if(ubi=='bsp') then
	    f_bnd=allfile(1:i_len)//'_bsp_band.dat'
	  endif
	  open(3,file=f_bnd)
c ----------------------------------------------------------
      do i=1,nbn
        do k=1,num_k_plot
          write(3,'(3E16.8)')k_dist(k),
     &      ((epl(i,k,ispin)-chem_pot)*evolt/2,ispin=1,nspin)
        enddo
        write(3,*) ' '
      enddo
      close(3)
c --------------------------------------------------------------
	  if(ubi=='dft') then
        if(mode==0) open(3,file=allfile(1:i_len)//'_dft_band_LAPW.gnu')
        if(mode==1) open(3,file=allfile(1:i_len)//'_dft_band_INTRP.gnu')
	  else if(ubi==' hf') then
	    open(3,file=allfile(1:i_len)//'_hf_band.gnu')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_qp_band.gnu')
	  else if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_gw_band.gnu')
	  else if(ubi=='psi') then
	    open(3,file=allfile(1:i_len)//'psi_band.gnu')
	  else if(ubi=='bsp') then
	    open(3,file=allfile(1:i_len)//'bsp_band.gnu')
	  endif
c ----------------------------------------------------------------
	  allocate(glabel(num_k_path+1))
      glabel(1)=' '//k_label(1,1)//' '
      do i=2,num_k_path
        if(k_label(2,i-1)/=k_label(1,i)) then
          glabel(i)=k_label(2,i-1)//'/'//k_label(1,i)
        else
          glabel(i)=' '//k_label(2,i-1)//' '
        end if
      end do
      glabel(num_k_path+1)=' '//k_label(2,num_k_path)//' '
c ---------------------------------------------------------------------
      write(3,*)
     & 'set loadpath '//ap//
     & 'C:/P4Win/Gnuplot/gnuplot/bin/share/PostScript/'//ap
      write(3,*)
     & 'set terminal postscript eps enhanced color "fontsize" 22'
	  if(ubi=='dft') then
        if(mode==0) write(3,*)'set output '//ap//allfile(1:i_len)//
     &                       '_lda_bnd_LAPW.eps'//ap
        if(mode==1) write(3,*)'set output '//ap//allfile(1:i_len)//
     &                       '_lda_bnd_INTRP.eps'//ap
	  else if(ubi==' hf') then
        write(3,*)'set output '//ap//allfile(1:i_len)//'_hf_bnd.eps'//ap
	  else if(ubi==' qp') then
        write(3,*)'set output '//ap//allfile(1:i_len)//'_qp_bnd.eps'//ap
	  else if(ubi==' gw') then
        write(3,*)'set output '//ap//allfile(1:i_len)//'_gw_bnd.eps'//ap
	  else if(ubi=='psi') then
        write(3,*)'set output '//ap//allfile(1:i_len)//
     &                        '_psi_bnd.eps'//ap
	  else if(ubi=='bsp') then
        write(3,*)'set output '//ap//allfile(1:i_len)//
     &                        '_bsp_bnd.eps'//ap
      endif
      write(3,701)k_dist(num_k_plot),emindos*evolt/2,emaxdos*evolt/2
      do i = 1, num_k_path-1
        write(3,705)sum(kpath_len(1:i)),emindos*evolt/2,
     &	            sum(kpath_len(1:i)),emaxdos*evolt/2
      enddo
      write(3,705)0.d0,0.d0,k_dist(num_k_plot),0.d0
      write(3,702, advance="no")glabel(1),0.d0,
     &	            (glabel(i+1),sum(kpath_len(1:i)),i=1,num_k_path-1)
      write(3,703)glabel(num_k_path+1),sum(kpath_len(:))
      write(3,*)'plot \'
      if(nspin==1) then
        write(3,*)ap//trim(f_bnd)//ap//' using 1:2 with lines ls 1'
      else if(nspin==2) then
        write(3,*)ap//trim(f_bnd)//ap//' using 1:2 with lines ls 1, \'
        write(3,*)ap//trim(f_bnd)//ap//' using 1:3 with lines ls 2'
      endif
      write(3,*)'unset ylabel'
      write(3,*)'unset output'
      write(3,*)'unset terminal'
      write(3,*)'unset loadpath'   
      close(3)
      deallocate(glabel,e)
701   format('unset colorbox',/,
     &       'set style line 1 lt 1 lc rgb "black" lw 3',/,
     &       'set style line 2 lt 1 lc rgb "red" lw 3',/,
     &       'set nokey',/,
     &       'set ylabel "Energy (eV)"',/,
     &       'set xrange [0:',F8.5,']',/,
     &       'set yrange [',F9.5,' :',F9.5,']')
702   format('set xtics (',:20('"',A3,'" ',F8.5,','))
703   format(A3,'" ',F8.5,')')
705   format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5,
     &       ' nohead')
      end 
