
C************************************************************************
C                            潜热的计�?
C************************************************************************
      SUBROUTINE phaselatentheat(shcon,nshcon,rhcon,
     &  nrhcon,v,vold,vini,nk,ne,lakon,xstateini,xstate,
     &  nstate,mi,ipkon,kon,ntmat,cphase,phase_inf,nphase,iout)
        IMPLICIT NONE
 
        character*8 lakon(*),lakonl
        integer phase_inf(4),nshcon(*),nrhcon(*),mi(*),
     &     ipkon(*),kon(*),nstate,nk,ne,ntmat,iout

       real*8   shcon(0:3,ntmat,*),rhcon(0:1,ntmat,*),
     &   v(0:mi(2),*),vold(0:mi(2),*),vini(0:mi(2),*),
     &   xstate(nstate,mi(1),*),xstateini(nstate,mi(1),*),
     &    cphase(13+phase_inf(4),phase_inf(1),*),rho,sph

      integer nonei20(3,12),nonei15(3,9),nonei10(3,6),
     &inum(nk),konl(20),nphase(11,*)

      integer I,J,l,k,indexe,node_index,nope,mint3d,two,
     & four,id,imat,nfield

      real*8  yn(nstate,nk),yn0(nstate,nk),a27(20,27),a9(6,9),
     & a8(8,8),a6(6,6),a4(4,4),a2(6,2),field(999,20*mi(3)),
     &field0(999,20*mi(3)),dxstate,dtemp

C      open(unit=10000,file='10000.dat')

      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
      data a2 /1.1455d0,-0.1455d0,1.1455d0,-0.1455d0,1.1455d0,-0.1455d0,
     &        -0.1455d0,1.1455d0,-0.1455d0,1.1455d0,-0.1455d0,1.1455d0/
      data a4 /  1.92705d0, -0.30902d0, -0.30902d0, -0.30902d0,
     &          -0.30902d0,  1.92705d0, -0.30902d0, -0.30902d0,
     &          -0.30902d0, -0.30902d0,  1.92705d0, -0.30902d0,
     &          -0.30902d0, -0.30902d0, -0.30902d0,  1.92705d0/
!
!     extrapolation from a 6 integration point scheme in a wedge to
!     the vertex nodes
!
      data a6 / 2.04904d0, 0.00000d0, 0.00000d0,-0.54904d0, 0.00000d0,
     &          0.00000d0,
     &         -0.34151d0, 1.70753d0,-0.34151d0, 0.09151d0,-0.45753d0,
     &          0.09151d0,
     &         -0.34151d0,-0.34151d0, 1.70753d0, 0.09151d0, 0.09151d0,
     &         -0.45753d0,
     &         -0.54904d0, 0.00000d0, 0.00000d0, 2.04904d0, 0.00000d0,
     &          0.00000d0,
     &          0.09151d0,-0.45753d0, 0.09151d0,-0.34151d0, 1.70753d0,
     &         -0.34151d0,
     &          0.09151d0, 0.09151d0,-0.45753d0,-0.34151d0,-0.34151d0,
     &          1.70753d0/
!
      data a9 / 1.63138d0,-0.32628d0,-0.32628d0,-0.52027d0, 0.10405d0,
     &          0.10405d0,
     &         -0.32628d0, 1.63138d0,-0.32628d0, 0.10405d0,-0.52027d0,
     &          0.10405d0,
     &         -0.32628d0,-0.32628d0, 1.63138d0, 0.10405d0, 0.10405d0,
     &         -0.52027d0,
     &          0.55556d0,-0.11111d0,-0.11111d0, 0.55556d0,-0.11111d0,
     &         -0.11111d0,
     &         -0.11111d0, 0.55556d0,-0.11111d0,-0.11111d0,0.55556d0,
     &         -0.11111d0,
     &         -0.11111d0,-0.11111d0, 0.55556d0,-0.11111d0,-0.11111d0,
     &          0.55556d0,
     &         -0.52027d0, 0.10405d0, 0.10405d0, 1.63138d0,-0.32628d0,
     &         -0.32628d0,
     &          0.10405d0,-0.52027d0, 0.10405d0,-0.32628d0, 1.63138d0,
     &         -0.32628d0,
     &          0.10405d0, 0.10405d0,-0.52027d0,-0.32628d0,-0.32628d0,
     &          1.63138d0/
!
!     extrapolation from a 2x2x2=8 integration point scheme in a hex to
!     the vertex nodes
!    
      data a8 /2.549d0,-.683d0,.183d0,-.683d0,-.683d0,.183d0,
     &        -.04904d0,.183d0,-.683d0,2.549d0,-.683d0,.183d0,
     &        .183d0,-.683d0,.183d0,-.04904d0,-.683d0,.183d0,
     &        -.683d0,2.549d0,.183d0,-.04904d0,.183d0,-.683d0,
     &        .183d0,-.683d0,2.549d0,-.683d0,-.04904d0,.183d0,
     &        -.683d0,.183d0,-.683d0,.183d0,-.04904d0,.183d0,
     &        2.549d0,-.683d0,.183d0,-.683d0,.183d0,-.683d0,
     &        .183d0,-.04904d0,-.683d0,2.549d0,-.683d0,.183d0,
     &        .183d0,-.04904d0,.183d0,-.683d0,-.683d0,.183d0,
     &        -.683d0,2.549d0,-.04904d0,.183d0,-.683d0,.183d0,
     &        .183d0,-.683d0,2.549d0,-.683d0/  
!
!     extrapolation from a 3x3x3=27 integration point scheme in a hex to
!     the all nodes in a 20-node element
!    
      data a27 /
     &  2.37499d0,-0.12559d0,-0.16145d0,-0.12559d0,-0.12559d0,
     & -0.16145d0, 0.11575d0,
     & -0.16145d0, 0.32628d0, 0.11111d0, 0.11111d0, 0.32628d0,
     &  0.11111d0,-0.10405d0,
     & -0.10405d0, 0.11111d0, 0.32628d0, 0.11111d0,-0.10405d0,
     &  0.11111d0,-0.31246d0,
     & -0.31246d0, 0.31481d0, 0.31481d0, 0.31481d0, 0.31481d0,
     & -0.16902d0,-0.16902d0,
     &  1.28439d0,-0.27072d0,-0.19444d0,-0.27072d0,-0.19444d0,
     &  0.15961d0,-0.00661d0,
     &  0.15961d0,-0.27072d0,-0.27072d0, 0.15961d0, 0.15961d0,
     & -0.12559d0, 2.37499d0,
     & -0.12559d0,-0.16145d0,-0.16145d0,-0.12559d0,-0.16145d0,
     &  0.11575d0, 0.32628d0,
     &  0.32628d0, 0.11111d0, 0.11111d0, 0.11111d0, 0.11111d0,
     & -0.10405d0,-0.10405d0,
     &  0.11111d0, 0.32628d0, 0.11111d0,-0.10405d0,-0.31246d0,
     &  0.31481d0, 0.31481d0,
     & -0.31246d0, 0.31481d0,-0.16902d0,-0.16902d0, 0.31481d0,
     & -0.27072d0,-0.19444d0,
     & -0.27072d0, 1.28439d0, 0.15961d0,-0.00661d0, 0.15961d0,
     & -0.19444d0,-0.27072d0,
     &  0.15961d0, 0.15961d0,-0.27072d0,-0.48824d0,-0.48824d0,
     & -0.48824d0,-0.48824d0,
     &  0.22898d0, 0.22898d0, 0.22898d0, 0.22898d0, 0.05556d0,
     &  0.05556d0, 0.05556d0,
     &  0.05556d0, 0.05556d0, 0.05556d0, 0.05556d0, 0.05556d0,
     & -0.22222d0,-0.22222d0,
     & -0.22222d0,-0.22222d0, 0.31481d0,-0.31246d0,-0.31246d0,
     &  0.31481d0,-0.16902d0,
     &  0.31481d0, 0.31481d0,-0.16902d0,-0.27072d0, 1.28439d0,
     & -0.27072d0,-0.19444d0,
     &  0.15961d0,-0.19444d0, 0.15961d0,-0.00661d0, 0.15961d0,
     & -0.27072d0,-0.27072d0,
     &  0.15961d0,-0.12559d0,-0.16145d0,-0.12559d0, 2.37499d0,
     & -0.16145d0, 0.11575d0,
     & -0.16145d0,-0.12559d0, 0.11111d0, 0.11111d0, 0.32628d0,
     &  0.32628d0,-0.10405d0,
     & -0.10405d0, 0.11111d0, 0.11111d0, 0.11111d0,-0.10405d0,
     &  0.11111d0, 0.32628d0,
     &  0.31481d0, 0.31481d0,-0.31246d0,-0.31246d0,-0.16902d0,
     & -0.16902d0, 0.31481d0,
     &  0.31481d0,-0.19444d0,-0.27072d0, 1.28439d0,-0.27072d0,
     & -0.00661d0, 0.15961d0,
     & -0.19444d0, 0.15961d0, 0.15961d0, 0.15961d0,-0.27072d0,
     & -0.27072d0,-0.16145d0,
     & -0.12559d0, 2.37499d0,-0.12559d0, 0.11575d0,-0.16145d0,
     & -0.12559d0,-0.16145d0,
     &  0.11111d0, 0.32628d0, 0.32628d0, 0.11111d0,-0.10405d0,
     &  0.11111d0, 0.11111d0,
     & -0.10405d0,-0.10405d0, 0.11111d0, 0.32628d0, 0.11111d0,
     & -0.31246d0, 0.31481d0,
     & -0.16902d0, 0.31481d0,-0.31246d0, 0.31481d0,-0.16902d0,
     &  0.31481d0,-0.27072d0,
     &  0.15961d0, 0.15961d0,-0.27072d0,-0.27072d0, 0.15961d0,
     &  0.15961d0,-0.27072d0,
     &  1.28439d0,-0.19444d0,-0.00661d0,-0.19444d0,-0.48824d0,
     & -0.48824d0, 0.22898d0,
     &  0.22898d0,-0.48824d0,-0.48824d0, 0.22898d0, 0.22898d0,
     &  0.05556d0,-0.22222d0,
     &  0.05556d0,-0.22222d0, 0.05556d0,-0.22222d0, 0.05556d0,
     & -0.22222d0, 0.05556d0,
     &  0.05556d0, 0.05556d0, 0.05556d0, 0.31481d0,-0.31246d0,
     &  0.31481d0,-0.16902d0,
     &  0.31481d0,-0.31246d0, 0.31481d0,-0.16902d0,-0.27072d0,
     & -0.27072d0, 0.15961d0,
     &  0.15961d0,-0.27072d0,-0.27072d0, 0.15961d0, 0.15961d0,
     & -0.19444d0, 1.28439d0,
     & -0.19444d0,-0.00661d0,-0.48824d0, 0.22898d0, 0.22898d0,
     & -0.48824d0,-0.48824d0,
     &  0.22898d0, 0.22898d0,-0.48824d0,-0.22222d0, 0.05556d0,
     & -0.22222d0, 0.05556d0,
     & -0.22222d0, 0.05556d0,-0.22222d0, 0.05556d0, 0.05556d0,
     &  0.05556d0, 0.05556d0,
     &  0.05556d0,-0.29630d0,-0.29630d0,-0.29630d0,-0.29630d0,
     & -0.29630d0,-0.29630d0,
     & -0.29630d0,-0.29630d0,-0.11111d0,-0.11111d0,-0.11111d0,
     & -0.11111d0,-0.11111d0,
     & -0.11111d0,-0.11111d0,-0.11111d0,-0.11111d0,-0.11111d0,
     & -0.11111d0,-0.11111d0,
     &  0.22898d0,-0.48824d0,-0.48824d0, 0.22898d0, 0.22898d0,
     & -0.48824d0,-0.48824d0,
     &  0.22898d0,-0.22222d0, 0.05556d0,-0.22222d0, 0.05556d0,
     & -0.22222d0, 0.05556d0,
     & -0.22222d0, 0.05556d0, 0.05556d0, 0.05556d0, 0.05556d0,
     &  0.05556d0, 0.31481d0,
     & -0.16902d0, 0.31481d0,-0.31246d0, 0.31481d0,-0.16902d0,
     &  0.31481d0,-0.31246d0,
     &  0.15961d0, 0.15961d0,-0.27072d0,-0.27072d0, 0.15961d0,
     &  0.15961d0,-0.27072d0,
     & -0.27072d0,-0.19444d0,-0.00661d0,-0.19444d0, 1.28439d0,
     &  0.22898d0, 0.22898d0,
     & -0.48824d0,-0.48824d0, 0.22898d0, 0.22898d0,-0.48824d0,
     & -0.48824d0, 0.05556d0,
     & -0.22222d0, 0.05556d0,-0.22222d0, 0.05556d0,-0.22222d0,
     &  0.05556d0,-0.22222d0,
     &  0.05556d0, 0.05556d0, 0.05556d0, 0.05556d0,-0.16902d0,
     &  0.31481d0,-0.31246d0,
     &  0.31481d0,-0.16902d0, 0.31481d0,-0.31246d0, 0.31481d0,
     &  0.15961d0,-0.27072d0,
     & -0.27072d0, 0.15961d0, 0.15961d0,-0.27072d0,-0.27072d0
     & , 0.15961d0,-0.00661d0,
     & -0.19444d0, 1.28439d0,-0.19444d0,-0.12559d0,-0.16145d0,
     &  0.11575d0,-0.16145d0,
     &  2.37499d0,-0.12559d0,-0.16145d0,-0.12559d0, 0.11111d0,
     & -0.10405d0,-0.10405d0,
     &  0.11111d0, 0.32628d0, 0.11111d0, 0.11111d0, 0.32628d0,
     &  0.32628d0, 0.11111d0,
     & -0.10405d0, 0.11111d0, 0.31481d0, 0.31481d0,-0.16902d0,
     & -0.16902d0,-0.31246d0,
     & -0.31246d0, 0.31481d0, 0.31481d0,-0.19444d0, 0.15961d0,
     & -0.00661d0, 0.15961d0,
     &  1.28439d0,-0.27072d0,-0.19444d0,-0.27072d0,-0.27072d0,
     & -0.27072d0, 0.15961d0,
     &  0.15961d0,-0.16145d0,-0.12559d0,-0.16145d0, 0.11575d0,
     & -0.12559d0, 2.37499d0,
     & -0.12559d0,-0.16145d0, 0.11111d0, 0.11111d0,-0.10405d0,
     & -0.10405d0, 0.32628d0,
     &  0.32628d0, 0.11111d0, 0.11111d0, 0.11111d0, 0.32628d0,
     &  0.11111d0,-0.10405d0,
     &  0.31481d0,-0.16902d0,-0.16902d0, 0.31481d0,-0.31246d0,
     &  0.31481d0, 0.31481d0,
     & -0.31246d0, 0.15961d0,-0.00661d0, 0.15961d0,-0.19444d0,
     & -0.27072d0,-0.19444d0,
     & -0.27072d0, 1.28439d0,-0.27072d0, 0.15961d0, 0.15961d0,
     & -0.27072d0, 0.22898d0,
     &  0.22898d0, 0.22898d0, 0.22898d0,-0.48824d0,-0.48824d0,
     & -0.48824d0,-0.48824d0,
     &  0.05556d0, 0.05556d0, 0.05556d0, 0.05556d0, 0.05556d0,
     &  0.05556d0, 0.05556d0,
     &  0.05556d0,-0.22222d0,-0.22222d0,-0.22222d0,-0.22222d0,
     & -0.16902d0, 0.31481d0,
     &  0.31481d0,-0.16902d0, 0.31481d0,-0.31246d0,-0.31246d0,
     &  0.31481d0, 0.15961d0,
     & -0.19444d0, 0.15961d0,-0.00661d0,-0.27072d0, 1.28439d0,
     & -0.27072d0,-0.19444d0,
     &  0.15961d0,-0.27072d0,-0.27072d0, 0.15961d0,-0.16145d0,
     &  0.11575d0,-0.16145d0,
     & -0.12559d0,-0.12559d0,-0.16145d0,-0.12559d0, 2.37499d0,
     & -0.10405d0,-0.10405d0,
     &  0.11111d0, 0.11111d0, 0.11111d0, 0.11111d0, 0.32628d0,
     &  0.32628d0, 0.11111d0,
     & -0.10405d0, 0.11111d0, 0.32628d0,-0.16902d0,-0.16902d0,
     &  0.31481d0, 0.31481d0,
     &  0.31481d0, 0.31481d0,-0.31246d0,-0.31246d0,-0.00661d0,
     &  0.15961d0,-0.19444d0,
     &  0.15961d0,-0.19444d0,-0.27072d0, 1.28439d0,-0.27072d0,
     &  0.15961d0, 0.15961d0,
     & -0.27072d0,-0.27072d0, 0.11575d0,-0.16145d0,-0.12559d0,
     & -0.16145d0,-0.16145d0,
     & -0.12559d0, 2.37499d0,-0.12559d0,-0.10405d0, 0.11111d0,
     &  0.11111d0,-0.10405d0,
     &  0.11111d0, 0.32628d0, 0.32628d0, 0.11111d0,-0.10405d0,
     &  0.11111d0, 0.32628d0,
     &  0.11111d0/

      if(nphase(5,1).ne.1)return
      imat=1;two=2;four=4
      nfield=nstate

      do i=1,nk
         inum(i)=0
      enddo

      do i=1,nk
         do j=1,nstate
            yn(j,i)=0.d0
            yn0(j,i)=0.d0
         enddo
      enddo
      
      do i=1,ne
          if(ipkon(i).lt.0) cycle

          if(lakon(i)(4:4).eq.'2') then
             nope=20
             if(lakon(i)(6:7).eq.'RA')then
               mint3d=4
             else
               mint3d=8
             endif
          elseif(lakon(i)(4:4).eq.'8') then
             nope=8
          elseif(lakon(i)(4:5).eq.'10') then
             nope=10
             mint3d=4
          elseif(lakon(i)(4:4).eq.'4') then
             nope=4
             mint3d=1
          elseif(lakon(i)(4:5).eq.'15') then
             nope=15
             mint3d=9
          elseif(lakon(i)(4:4).eq.'6') then
             nope=6
             mint3d=2
          else
             nope=0
             mint3d=0
          endif

          indexe=ipkon(i)!第i个单元第一个节点的数组编号
          do j=1,nope
              konl(j)=kon(indexe+j)
              inum(konl(j))=inum(konl(j))+1
          enddo

          lakonl=lakon(i)
C************************************************************************
C               积分点结果外推至节点 
C************************************************************************

!        determining the field values in the vertex nodes
!        for C3D20R and C3D8: trilinear extrapolation (= use of the
!                             C3D8 linear interpolation functions)
!        for C3D8R: constant field value in each element
!        for C3D10: use of the C3D4 linear interpolation functions
!        for C3D4: constant field value in each element
!        for C3D15: use of the C3D6 linear interpolation functions
!        for C3D6: use of a linear interpolation function
            if((lakonl(4:6).eq.'20R').or.
     &         (lakonl(4:5).eq.'8 ').or.(lakonl(4:5).eq.'8I')) then
               if(lakonl(7:8).ne.'LC') then
                  do j=1,8
                     do k=1,nfield
                        field(k,j)=0.d0
                        field0(k,j)=0.d0
                        do l=1,8
                           field(k,j)=field(k,j)+a8(j,l)*xstate(k,l,i)
                           field0(k,j)=field0(k,j)+a8(j,l)*
     &                      xstateini(k,l,i)
                        enddo
                     enddo
                  enddo
               endif
            elseif(lakonl(4:4).eq.'8') then
               do j=1,8
                  do k=1,nfield
                     field(k,j)=xstate(k,1,i)
                     field0(k,j)=xstateini(k,1,i)
                  enddo
               enddo
            elseif(lakonl(4:5).eq.'10') then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=0.d0
                     field0(k,j)=0.d0
                     do l=1,4
                        field(k,j)=field(k,j)+a4(j,l)*xstate(k,l,i)
                        field0(k,j)=field0(k,j)+a4(j,l)*
     &                      xstateini(k,l,i)
                     enddo
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'2') then
               do j=1,20
                  do k=1,nfield
                     field(k,j)=0.d0
                     field0(k,j)=0.d0
                     do l=1,27
                        field(k,j)=field(k,j)+a27(j,l)*xstate(k,l,i)
                        field0(k,j)=field0(k,j)+a27(j,l)*
     &                      xstateini(k,l,i)
                     enddo
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'4') then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=xstate(k,1,i)
                     field0(k,j)=xstateini(k,1,i)
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'1') then
               if(lakonl(7:8).ne.'LC') then
                  do j=1,6
                     do k=1,nfield
                        field(k,j)=0.d0
                        field0(k,j)=0.d0
                        do l=1,9
                           field(k,j)=field(k,j)+a9(j,l)*xstate(k,l,i)
                           field0(k,j)=field0(k,j)+a9(j,l)*
     &                          xstateini(k,l,i)
                        enddo
                     enddo
                  enddo
               endif
            else
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     field0(k,j)=0.d0
                     do l=1,2
                        field(k,j)=field(k,j)+a2(j,l)*xstate(k,l,i)
                        field0(k,j)=field0(k,j)+a2(j,l)*
     &                       xstateini(k,l,i)
                     enddo
                  enddo
               enddo
            endif
!
!        determining the field values in the midside nodes
!
         if(lakonl(4:6).eq.'20R') then
            if(lakonl(7:8).ne.'LC') then
               do j=9,20
                  do k=1,nfield
                     field(k,j)=(field(k,nonei20(2,j-8))+
     &                    field(k,nonei20(3,j-8)))/2.d0
                     field0(k,j)=(field0(k,nonei20(2,j-8))+
     &                    field0(k,nonei20(3,j-8)))/2.d0
                  enddo
               enddo
            endif
         elseif(lakonl(4:5).eq.'10') then
            do j=5,10
               do k=1,nfield
                  field(k,j)=(field(k,nonei10(2,j-4))+
     &                 field(k,nonei10(3,j-4)))/2.d0
                  field0(k,j)=(field0(k,nonei10(2,j-4))+
     &                 field0(k,nonei10(3,j-4)))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
               do j=7,15
                  do k=1,nfield
                     field(k,j)=(field(k,nonei15(2,j-6))+
     &                    field(k,nonei15(3,j-6)))/2.d0
                     field0(k,j)=(field0(k,nonei15(2,j-6))+
     &                    field0(k,nonei15(3,j-6)))/2.d0
                  enddo
               enddo
            endif
         endif
!
!        transferring the field values into yn
!
         if(lakonl(7:8).ne.'LC') then
            do j=1,nope
               do k=1,nfield
                  yn(k,kon(indexe+j))=yn(k,kon(indexe+j))+
     &                 field(k,j)
                  yn0(k,kon(indexe+j))=yn0(k,kon(indexe+j))+
     &                 field0(k,j)
               enddo
            enddo
         endif
c     Bernhardi start
c        incompatible modes elements
         if(lakonl(1:5).eq.'C3D8I') then
            do j=1,3
               do k=1,nfield
                  yn(k,kon(indexe+nope+j))=0.0d0
                  yn0(k,kon(indexe+nope+j))=0.0d0
               enddo
            enddo
         endif
      enddo

      do i=1,nk
        if(inum(i).gt.0) then
C************************************************************************
C                     确定密度 
C************************************************************************
          call ident2(rhcon(0,1,imat),v(0,i),nrhcon(imat),two,id)
          if(nrhcon(imat).eq.0) then
            continue
          elseif(nrhcon(imat).eq.1) then
           rho=rhcon(1,1,imat)
          elseif(id.eq.0) then
           rho=rhcon(1,1,imat)
          elseif(id.eq.nrhcon(imat)) then
           rho=rhcon(1,id,imat)
          else
           rho=rhcon(1,id,imat)+
     &     (rhcon(1,id+1,imat)-rhcon(1,id,imat))*
     &            (v(0,i)-rhcon(0,id,imat))/
     &     (rhcon(0,id+1,imat)-rhcon(0,id,imat))
          endif
C************************************************************************
C                     确定比热 
C************************************************************************
          call ident2(shcon(0,1,imat),v(0,i),nshcon(imat),four,id)
          if(nshcon(imat).eq.0) then
            continue
          elseif(nshcon(imat).eq.1) then
            sph=shcon(1,1,imat)
          elseif(id.eq.0) then
            sph=shcon(1,1,imat)
          elseif(id.eq.nshcon(imat)) then
            sph=shcon(1,id,imat)
          else
            sph=shcon(1,id,imat)+
     &      (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &      (v(0,i)-shcon(0,id,imat))/
     &      (shcon(0,id+1,imat)-shcon(0,id,imat))
          endif
C************************************************************************
C                      计算节点相变结果 
C************************************************************************
           do j=2,phase_inf(1)
             yn(j,i)=yn(j,i)/inum(i)
             yn0(j,i)=yn0(j,i)/inum(i)
             dxstate=dabs(yn(j,i)-yn0(j,i))
             if(dxstate.gt.0.d0)then
C************************************************************************
C                计算节点温度增量以及更新温度 
C************************************************************************
               dtemp=dxstate*cphase(2,j,imat)/sph/rho
                v(0,i)=v(0,i)+dtemp
                vold(0,i)=v(0,i)
C                 write(10000,*)yn(j,i),yn0(j,i)
C                 write(10000,*)"打印�?0",xstateini(j,1,1)
C                 write(10000,*)"打印�?,xstate(j,1,1)
C                 write(10000,*)dxstate
C               endif
             endif
           enddo
        endif
!         write(10000,*)v(0,i)
      enddo

      RETURN
      END
