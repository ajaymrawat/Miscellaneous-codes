! 


program ics_jne0
        implicit none
        integer, parameter :: njmax=61,nvmax=5,nkmax=61,nemax=601,jbstart=0,&
                              jbmax=81,jinit=1,jbend=80,isym=2
        integer      :: ip(4),kstart(4),i,ij,iv,j,k,ne,nvab,njab,nk,kkstart,jjstart,&
                        iip,ibig,jqm(njmax),inj(nkmax),ie,iq
        real*8       :: evj0,emur,pi,ics,gs,gk,gkp,tpa(nemax,nvmax,njmax), &
                        prob(nvmax,njmax,nkmax),ttpa(nemax,nvmax,njmax),pa(nemax,nvmax,njmax),pk,&
                        eray(nemax),const
        complex*16   :: ss(nvmax,njmax,nkmax)
        real*8, parameter  :: autev=27.211648d0,b2asq=0.280028557d0,evjstart=0.097350d0
        character*7  :: ciparity
        character*4  :: cjbig
        character*30 :: pedat1,pedat2,pedat3,pedat4
        real*8       :: t1,t2

        call cpu_time(t1)

        pi=4.d0*datan(1.d0)

        pedat1='bigpedat_o3_gr_j1_k0p0.p'
        pedat2='bigpedat_o3_gr_j1_k0p1.p'
        pedat3='bigpedat_o3_gr_j1_k1p0.p'
        pedat4='bigpedat_o3_gr_j1_k1p1.p'

        ip(1)=0;ip(2)=1;ip(3)=0;ip(4)=1
        kstart(1)=0;kstart(2)=0;kstart(3)=1;kstart(4)=1

        open(81,file=pedat1,form='unformatted')
        open(82,file=pedat2,form='unformatted')
        open(83,file=pedat3,form='unformatted')
        open(84,file=pedat4,form='unformatted')
! initialization of arrays
        ttpa(:,:,:)=0.d0
        tpa(:,:,:)=0.d0
        ss(:,:,:)=0.d0
        pa(:,:,:)=0.d0

        do i=1,2*(jinit+1)
           tpa(:,:,:)=0.d0
           print*,''
           print*,'pedat    :',kstart(i),ip(i)
           do ij=0,jbend
              ss(:,:,:)=0.d0
!              if (ij.eq.22) stop
              if (((kstart(i)+ip(i)).ne.0).and.(ij.eq.0)) goto 1
              print*,'J :',ij
              read(80+i)ciparity,iip,cjbig,ibig
              read(80+i)ne,nvab,njab,nk,jjstart,kkstart,evj0,emur
              read(80+i)(jqm(iq),iq=1,njab)
              read(80+i)(inj(k),k=1,nk)
1             continue

              do ie=1,ne
                 
                 if (((kstart(i)+ip(i)).ne.0).and.(ij.eq.0)) then
                         ss(:,:,:)=0.d0
                         if (ie.eq.1) print*,'ij inside      :',ij
                         pa(:,:,:)=0.d0
                 else
                         read(80+i)eray(ie)
                         read(80+i)(((ss(iv,j,k),iv=1,nvab),j=inj(k),njab),k=1,nk)
                         if ((ij.eq.20).and.(ie.eq.5)) print*,ss(1,2,1)                         
                 prob(:,:,:)=dble(isym)*(ss(:,:,:)*conjg(ss(:,:,:)))
                 do iv=1,nvab
                    do j=1,njab
                       ! gs=post-antisymmetry factor ( only when ...
                       ! ... the nuclei in product diatom are same)
                       if(mod((j-1),2).eq.0) then
                               gs=(1.d0/2.d0)
                       else
                               gs=(3.d0/2.d0)
                       endif
                       pk=0.d0
                       do k=1,min(ij+1,j)
                          if ((k-1).eq.0) then
                                  gkp=1.d0
                          else
                                  gkp=1.d0/2.d0
                          endif
                          pk=pk+(gkp*prob(iv,j,k))

!                       if (ie.eq.4) then
!                       if ((prob(iv,j,k).ne.0)) write(13,*)ij,ibig,prob(iv,j,k),iv,j,k

!                       endif
                       enddo     ! k loop end (omega')
                       pa(ie,iv,j)=dble((2*ij+1))*pk
                    enddo        ! j loop end (j')
                 enddo           ! iv loop end (v')
                 endif
              enddo              ! ie loop end (E)
              tpa(:,:,:)=pa(:,:,:)+tpa(:,:,:)
           enddo                 ! ij loop end (J)
           
           if (kstart(i).eq.0) then
                   gk=1.d0/(2*jinit+1)
           else
                   gk=2.d0/(2*jinit+1)
           endif
           ttpa(:,:,:)=gk*tpa(:,:,:)+ttpa(:,:,:)
        enddo                    ! loop over omega and p

! total state-selected ics
        do ie=1,ne
           ics=0.d0
           const=(pi)/(2*emur*(eray(ie)-(evjstart/autev)))
           do iv=1,nvab
              do j=1,njab
                 ics=ics+(const*ttpa(ie,iv,j))!*b2asq
              enddo
           enddo
           if (ie.ne.1) write(10,*) ((eray(ie)*autev)-evjstart),ics
        enddo
        call cpu_time(t2)
        stop
        end
