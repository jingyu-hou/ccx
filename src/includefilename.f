!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine includefilename(text,includefn,lincludefn)
!
!     determines the name of an include file
!
      implicit none
!
      character*132 includefn
      character*1320 text
!
      integer nstart,nend,ii,jj,kk,lincludefn
!
      nstart=0
      nend=0
!
      loop: do ii=1,lincludefn
        if(text(ii:ii).eq.'=') then
           jj=ii+1
           if(text(jj:jj).eq.'"') then
              nstart=jj+1
              do kk=jj+1,lincludefn
                 if(text(kk:kk).eq.'"') then
                    nend=kk-1
                    exit loop
                 endif
              enddo
              write(*,*)'*ERROR in includefilename: ',
     &             'finishing quotes are lacking'
              write(*,*) '*ERROR in the input deck. Card image:'
              write(*,'(132a1)') 
     &             (text(kk:kk),kk=1,lincludefn)
              call exit(201)
           else
              nstart=jj
              nend=lincludefn
              exit
           endif
        endif
      enddo loop
      if(ii.eq.lincludefn+1) then
         write(*,*) '*ERROR in includefilename: syntax error'
         write(*,*) '*ERROR in the input deck. Card image:'
         write(*,'(132a1)')
     &             (text(kk:kk),kk=1,lincludefn)
         call exit(201)
      endif
!
      if(nend.ge.nstart) then
         if(nend-nstart+1.le.132) then
            includefn(1:nend-nstart+1)=text(nstart:nend)
            lincludefn=nend-nstart+1
         else
            write(*,*) '*ERROR in includefilename: file name too long'
            write(*,*) '*ERROR in the input deck. Card image:'
            write(*,'(132a1)')
     &             (text(kk:kk),kk=1,lincludefn)
            call exit(201)
         endif
      else
         write(*,*) '*ERROR in includefilename: file name is lacking'
         write(*,*) '*ERROR in the input deck. Card image:'
         write(*,'(132a1)')
     &             (text(kk:kk),kk=1,lincludefn)
         call exit(201)
      endif
!
      return
      end




