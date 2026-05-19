!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine qsorti(n,list,key)
!     
      implicit none
!     
      integer list(*),key(*),n,ll,lr,lm,nl,nr,ltemp,stktop,maxstk,guess
!     
      parameter(maxstk=32)
!     
      integer lstack(maxstk),rstack(maxstk)
!     
      ll= 1
      lr=n
      stktop=0
 10   if(ll.lt.lr) then
         nl=ll
         nr=lr
         lm=(ll+lr)/2
         guess=key(list(lm))
 20      if (key(list(nl)).lt.guess) then
            nl=nl+1
            goto 20
         end if
 30      if (guess.lt.key(list(nr))) then
            nr=nr-1
            goto 30
         end if
         if(nl.lt.(nr-1)) then
            ltemp=list(nl)
            list(nl)=list(nr)
            list(nr)=ltemp
            nl=nl+1
            nr=nr-1
            goto 20
         end if
         if(nl.le.nr) then
            if(nl.lt.nr) then
               ltemp=list(nl)
               list(nl)=list(nr)
               list(nr)=ltemp
            end if
            nl=nl+1
            nr=nr-1
         end if
         stktop=stktop+1
         if(nr.lt.lm) then
            lstack(stktop)=nl
            rstack(stktop)=lr
            lr=nr
         else
            lstack(stktop)=ll
            rstack(stktop)=nr
            ll=nl
         end if
         goto 10
      end if
      if (stktop.ne.0) then
         ll=lstack(stktop)
         lr=rstack(stktop)
         stktop=stktop-1
         goto 10
      end if
      end
      
