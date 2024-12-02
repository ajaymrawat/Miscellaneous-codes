! @Ajay M. Rawat - Dec 2, 2024. (Tested with He + LiH+ -----> H + LiHe+ reaction)
! This module contains the subroutines which calculates the following observables:
!     1. Total integral cross section: Ecol vs ICS (Subroutine: total_ics) 2-D
!     2. Product vibrational level resolved ICS: Ecol vs ICS(v') (Subroutine: vibres_ics) 2-D
!     3. Product rotational level resolved ICS: Ecol vs ICS(v',j') (Subroutine: rot_res_ics) 2-D
!     4. Product rotational and vibrational level resolved ICS: Ecol vs j' vs ICS(Ecol, v', j') 3-D
module ics_mods

        integer, parameter:: njmx=16,nvabmx=6,nkmx=16,nemax=1001, &
                             jbmax=51,jbstart=0,jbend=50,isym=1
        real*8, parameter :: autev=27.211648d0,b2asq=0.280028557d0,&
                             evib=0.0153471010632d0/autev
        integer :: njab,nvab,ne,ip,jbig,nk,jstart,kstart,i,j,jb,ie, &
                   k,iv,gk,jqm(njmx),inj(nkmx)
        real*8  :: eray(nemax),prob(nvabmx,njmx,nkmx),const, &
                 pa(nemax,nvabmx,njmx),tpa(nemax,nvabmx,njmx),ics,evj0

        contains

subroutine calc_prob     ! Calculates the S-matrix and reaction probability 
      implicit none
      real*8  ::      emur,pi,p,le,gs
      complex*16   :: ss(nvabmx,njmx,nkmx)
      character*7  :: ipc
      character*4  :: jbc
      character*28 :: outputfile1,outputfile2
      character*19 :: outputfile3

      pi=4.d0*datan(1.d0)

      open(80,file='pedat_lihhep_gr_v0j0_J000-050.p',form='unformatted')
      !open(80,file='big_pedat_mod.p',form='unformatted')
! all the memory locations are declared zero at the starting 
      ss(:,:,:)=0.d0
      prob(:,:,:)=0.d0
      pa(:,:,:)=0.d0
      tpa(:,:,:)=0.d0

! starting the j (small j) loop
      do i=jbstart,jbend
      read(80) ipc,ip,jbc,jb         !unfomatted files are read without format specification
      read(80) ne,nvab,njab,nk,jstart,kstart,evj0,emur
      read(80) (jqm(j), j=1,njab)
      read(80) (inj(k), k=1,nk)
      print*,'i',i,jb
! t/(eray(le)-evj0))*tpa(le,iv,j)*b2asqtarting energy loop
      do ie=1,ne
      read(80) eray(ie)
      read(80)(((ss(iv,j,k),iv=1,nvab),j=inj(k),njab),k=1,nk)
      prob(:,:,:)=dble(isym)*(ss(:,:,:)*conjg(ss(:,:,:)))
! starting vibration loop       
      do iv=1,nvab
      do j=1,njab
      if(mod((j-1),2).eq.0) then
      gs=(1.d0/2.d0)
      else
      gs=(3.d0/2.d0)
      endif
      gs=1.d0
      p=0.0d0
!      p=sum(prob(iv,j,:))+p   !commented when k loop is there. 
      do k=1,min(i+1,j)
      p=p+prob(iv,j,k)
      enddo  !---k loop ends here
!      pa(ie,iv,j)=gs*dble((2*i+1))*p      ! with gs
      pa(ie,iv,j)=dble((2*i+1))*p        ! without gs
      enddo !---j loop ends here
      enddo !---v (iv) loop ends here
      enddo !---energy (ie) loop ends here 
      tpa(:,:,:)=pa(:,:,:)+tpa(:,:,:)  !  adding the prob for different J (jbig) values
      enddo !---J (i) loop ends here

      if(kstart==0) then
        gk=1
      else
        gk=2
      endif

      const=(pi*dble(gk))/(2.0d0*emur*dble(2*jstart+1)) ! constant in the formula 

      end subroutine

subroutine total_ics   ! Ecol vs ICS(Ecol) : ICS is summed over v' and j' 
        do ie=1,ne
        ics=0.d0
        do j=1,njab
        do iv=1,nvab
        ics=ics+(const/(eray(ie)-evib))*tpa(ie,iv,j)*b2asq
        enddo !j loop
        enddo ! iv loop
        write(140,*) ((eray(ie)-evib)*autev),ics
        enddo ! ne loop
end subroutine

subroutine vibres_ics ! Ecol vs ICS(Ecol,v') : ICS is summed over j'
      character*14 :: outputfile2
      real*8       :: tpvib

      outputfile2="stsics_v00.dat"

      do iv=1,nvab   !writing the results in files
      write(unit=outputfile2(9:9),fmt='(i1)') ((iv-1)/10)
      write(unit=outputfile2(10:10),fmt='(i1)') mod((iv-1),10)
      open(11,file=trim(outputfile2))

      do ie=1,ne
      tpvib=0.d0
      do j=1,njab
      tpvib=tpvib+tpa(ie,iv,j)
      enddo  ! ending ij loop
      if (ie.ne.1) then
      write(11,*) (eray(ie)-evib)*autev, (const/(eray(ie)-evib))*tpvib*b2asq
      endif
      enddo  ! ending ie loop
      close(11)
      enddo  ! ending iv loop

end subroutine


subroutine vib_rot_res_ics ! Ecol vs vs j' vs ICS(Ecol,v',j') : Each v' have different files

      character*28 :: outputfile3
      outputfile3='s2sics_v00_vs_j.dat'
      do iv=1,nvab   !writing the results in files
      write(unit=outputfile3(9:9),fmt='(i1)') ((iv-1)/10)
      write(unit=outputfile3(10:10),fmt='(i1)') mod((iv-1),10)
      open(15,file=trim(outputfile3))
      do j=1,njab
      do ie=1,ne
!--- for specific vib level all the j (3-d plot)
      if (ie.ne.1) then
      write(15,*) (eray(ie)-evib)*autev, j-1, (const/(eray(ie)-evib))*tpa(ie,iv,j)*b2asq
      endif
!--- for specific vib and rot level (2-d plots)
      enddo  !--- end of ie loop
      write(15,*)
      enddo  !-- end of j loop
      close(15)
      enddo  !---end of iv loop
end subroutine

subroutine rot_res_ics             ! Ecol vs ICS(v',j') : Each combination of v' and j' has separate file.

      character*28::outputfile1
      outputfile1="s2sics_aj_v00_j00_gs.dat"

      do iv=1,nvab   !writing the results in files
      do j=1,njab
      write(unit=outputfile1(12:12),fmt='(i1)') ((iv-1)/10)
      write(unit=outputfile1(13:13),fmt='(i1)') mod((iv-1),10)
      write(unit=outputfile1(16:16),fmt='(i1)') ((j-1)/10)
      write(unit=outputfile1(17:17),fmt='(i1)') mod((j-1),10)
      open(10,file=trim(outputfile1))
      do ie=1,ne
!--- for specific vib level all the j (3-d plot)
      if (ie.ne.1) then
      write(10,*) (eray(ie)-evib)*autev, (const/(eray(ie)-evib))*tpa(ie,iv,j)*b2asq
      endif
      enddo  !--- end of ie loop
      enddo  !-- end of j loop
      enddo  !---end of iv loop

end subroutine

subroutine en_res_ics(ren)
      real*8 :: ren
      character*12 :: enoutfile

      do i=1,ne
      if(int((eray(i)*autev)-ren).eq.0) then
      ie=i

      print*,ie,eray(ie)
      endif
      enddo

      enoutfile= 'ics_e00.dat'
      write(enoutfile(6:6),'(i1)') ie/10
      write(enoutfile(7:7),'(i1)') mod(ie,10)
      open(14,file=trim(enoutfile))
      do iv=1,nvab
      do j=1,njab
      write(14,*) (iv-1),(j-1), (const/(eray(le)-evj0))*tpa(le,iv,j)*b2asq
      enddo  !--- end of j loop
      write(14,*)
      enddo  !-- end of iv loop

end subroutine

end module ics_mods

