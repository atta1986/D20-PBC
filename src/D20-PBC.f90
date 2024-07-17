   program D20
   implicit real*8 (a-h,o-z)
!   implicit none
   Integer :: nsa,nsb,typ,nsmax,vers,maxat,max_elem,maxc,nn
   Real*8 :: aswitch,bswitch,conv,a0,zero,autoang,autokcal,c6conv,autoev
   parameter (nsmax=1000,aswitch=1.0,bswitch=8.8,conv=0.529177209,vers=3)

   parameter (maxat=20000)
   parameter (max_elem=94)
! maximum coordination number references per element
   parameter (maxc=5)
! coversion factors
   parameter (autoang=0.52917726d0)
   parameter (autokcal=627.509541d0)
   parameter (autoev= 27.21138505)

   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr,disper,rad2d,pi
!   Integer :: nsa,nsb,typ
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   Real*8,dimension(:,:) :: sita(5,nsmax),sitb(5,nsmax),sitaa(5,nsmax),sitbb(5,nsmax),&
   & ama(nsmax),amb(nsmax),tran(3,3)
   Real, dimension(:) :: aol(3), bol(3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(5,nsmax) :: latom(5,nsmax),atom(5,nsmax)
   Integer,dimension(:):: chara(nsmax),charb(nsmax)
   dimension:: bweight(nsmax,nsmax),sit_bnd(3)
   data a0 / 0.529177249d0/, zero /0.d0/
   DATA doit/'DOIT'/,done/'DONE'/,doitr6/'DOR6'/
   data printed /.false./
   Logical :: switch_on,onesystem
!   pi = dacos(-1.d0)
!   rad2d = 180.d0/pi
!      implicit none             
!   integer :: maxat,max_elem,maxc                      
! conversion factors
!   real*8 :: autoang,autokcal,c6conv,autoev
!   parameter (maxat=20000)
!   parameter (max_elem=94)
! maximum coordination number references per element
!   parameter (maxc=5)
! coversion factors
!   parameter (autoang=0.52917726d0)
!   parameter (autokcal=627.509541d0)
!   parameter (autoev= 27.21138505)

! DFT-D version
   integer version
! number of atoms
!   integer n
! coordinates in au
   real*8,dimension(:,:), allocatable :: xyz,abc
! fixed atoms in geometry opt
   logical fix(maxat)
! lattice in au
   real*8 lat(3,3)
! gradient
   real*8,dimension(:,:), allocatable :: g      
   real*8 g_lat(3,3)
! cardinal numbers of elements
   integer,dimension(:), allocatable :: iz  
! cut-off radii for all element pairs
   real*8 r0ab(max_elem,max_elem)
! C6 for all element pairs 
   real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
! how many different C6 for one element
   integer mxc(max_elem)
! C6810 
   real*8 c6,c8,c10
! coordination numbers of the atoms
   real*8,dimension(:), allocatable :: cn  
! covalent radii
   real*8 rcov(max_elem)
! atomic <r^2>/<r^4> values
   real*8 r2r4(max_elem)
! energies
   real*8 e6, e8, e10, e12, disp, e6abc        
! THE PARAMETERS OF THE METHOD (not all a "free")
   real*8 rs6, rs8, rs10, s6, s18, alp6, alp8, alp10, s42, rs18, alp
! printout option
   logical echo
! grad ?
   logical grad
! analyse results ?
   logical anal
! third-order term?
   logical noabc
! gradient calctype
   logical numgrad
! special parameters
   logical tz
! periodic boundary conditions
   logical pbc
! repetitions of the unitcell to match the rthr and c_thr
   integer rep_vdw(3),rep_cn(3)
! R^2 distance neglect threshold (important for speed in case of large systems)
   real*8 rthr,rthr2
! R^2 distance to cutoff for CN_calculation
   real*8 cn_thr
! Integer for assigning only max/min cn C6 (0=normal, 1=min, 2=max)
! local and dummy variables
   character*80 atmp,btmp,ctmp,dtmp,etmp,ftmp,func
   character*2  esym 
   integer i,j,k,z,iat,jat,i1,i2
   integer ida(max_elem),ipot
   real*8  x,y,dispr,displ,gdsp,dum,xx(10),dum6(86)
   real*8  dum1,dum2,dum3(3)
   logical ex,pot,fdum
   logical minc6list(max_elem),maxc6list(max_elem),minc6,maxc6
! these new data are scaled with k2=4./3.  and converted a_0 via
! autoang=0.52917726d0
   data rcov/&
   & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,&
   & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,&
   & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,&
   & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,&
   & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,&
   & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,&
   & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,&
   & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,&
   & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,&
   & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,&
   & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,&
   & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,&
   & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,&
   & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,&
   & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,&
   & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,&
   & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,&
   & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,&
   & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

! k1-k3
!   include 'param'
   echo=.true. 
   grad=.false.
   pot =.false.
   anal=.false.
   noabc=.true. 
   numgrad=.false.
   tz=.false.
   func=' none (read from parameter file)'
   version=0
!   pbc=.false.
   pbc=.true.
   fix=.false.
   minc6=.false.
   maxc6=.false.
   minc6list=.false.
   maxc6list=.false.
   fdum=.false.
! Cutoff r^2 thresholds for the gradient in bohr^2.
! rthr influences N^2 part of the gradient.
! rthr2 influences the N^3 part of the gradient. When using
! dftd3 in combination with semi-empirical methods or FFs, and large
! (>1000 atoms) systems, rthr2 is crucial for speed:
! Recommended values are 20^2 to 25^2 bohr.

   rthr=18000.0d0   ! UR, SE
   rthr2=1600.0d0
   cn_thr=2000.0d0

! J/mol nm^6 - > au
   c6conv=1.d-3/2625.4999d0/((0.052917726d0)**6)
   pi = dacos(-1.d0)
   rad2d = 180.d0/pi

! get coord filename
   call getarg(1,etmp)
   inquire(file=etmp,exist=ex)
!   if(.not.ex) call printoptions       
   ex=.false.
   ipot=0


   switch_on=.true.
   onesystem=.true.
!   onesystem=.false.
!   Read in the sites of monomer A
!   WRITE(fileA,'(a)') "fileA.dat"
!   if(onesystem) then
!     open(unit=1,file='file-nbc10ext-fortran.dat',form='formatted',status='old')
!     read(1,*)nsa
!     do i=1,nsa
!         read(1,*)(sita(j,i),j=1,4)
!     end do
!     nsb=nsa
!     sitaa=sita
!     sitbb=sita
!   else
!!     print*,'else=',nsa
!     open(unit=1,file='file.dat',form='formatted',status='old')
!     read(1,*)nsa
!   read(1,*)
!   open(99,file='file.dat')
!!   write(99,*)nsa

!     do i=1,nsa
!         read(1,*)(sita(j,i),j=1,4)
!!         print*,(sita(j,i),j=1,4)
!!         print*,'****************'
!!         print*,sita(:,1)
!     end do

!     read(1,*)nsb
   
!     do i=1,nsb
!         read(1,*)(sitb(j,i),j=1,4)
!!         print*,(sitb(j,i),j=1,4)
!!         print*,ltempb(i)
!     end do  
!     close(1)
!     sitaa=sita
!     sitbb=sitb
!   end if
!   call analyze_hydrogens(sita,nsa,sitaa)
!   sita=sitaa
!   do i=1,nsa
!      print*,'sita=',sita(:,i)
!   end do
!   print*,'*****************'
!   call analyze_carbons(sita,nsa,sitaa)
!   sita=sitaa
!   do i=1,nsa
!      print*,'sita=',sita(:,i)
!   end do
!   print*,'*****************'
!   do i=1,nsa      
!      print*,sita(5,i)
!   end do
!   call analyze_hydrogens(sitb,nsb,sitbb)
!   sitb=sitbb
!   do i=1,nsb
!      print*,sitb(:,i)
!   end do
!   print*,'*****************'
!   call analyze_carbons(sitb,nsb,sitbb)
!   sitb=sitbb
!   do i=1,nsb
!      print*,sitb(:,i)
!   end do
!   print*,sita(:,1)
!   Call fitted_dispersion(sita,sitb,nsa,nsb,switch_on,aswitch,bswitch,disper)
!   print*,disper
!
! For Periodic systems (from Grimme D3 dispersion function)
!
   if (pbc) then
     call pbcrdatomnumber(etmp,nn)
!   else
!     call rdatomnumber(etmp,n)
   end if
   !print*,nn
   allocate(xyz(3,nn))
   allocate(g(3,nn))
   allocate(iz(nn))
   allocate(cn(nn))
! reading coordinates and cell in VASP.5.2-format
! determing repetitions of unitcell
   if (pbc) then
           call pbcrdcoord(etmp,lat,nn,xyz,iz,autoang,sita)
           call set_criteria(rthr,lat,dum3)
           rep_vdw=int(dum3)+1
           !print*,'rep_vdw(1-3)',rep_vdw(1),rep_vdw(2),rep_vdw(3)
           call set_criteria(cn_thr,lat,dum3)
           rep_cn=int(dum3)+1
           !print*,'rep_cn(1-3)',rep_cn(1),rep_cn(2),rep_cn(3)
   end if
   if(nn.lt.1)     call stoprun( 'no atoms' )
   if(nn.gt.maxat) call stoprun( 'too many atoms' )
!
!   print*,'sita',sita(1,1),sita(2,1),sita(3,1)
   call analyze_hydrogens(sita,nn,lat,rep_cn,cn_thr,sitaa)
   sita=sita+sitaa
!   Do i=1,nn
!     print*,sita(4,i),sitaa(5,i),sita(5,i)
     !sita(5,i)=sitaa(5,i)
!   end do
!
   call analyze_carbons(sita,nn,lat,rep_cn,cn_thr,sitbb)
   sita=sita+sitbb
   !print*,'Number of atoms in unit cell: ',nn
   !Do i=1,nn
   !   print*,sita(4,i),sitbb(5,i),sita(5,i)
!     sita(5,i)=sitbb(5,i)
   !end do
!
!   print*,'sita',sita(1,1),sita(2,1),sita(3,1),sita(4,1)
!  Call PBC subroutine:
   call  pbcedisp(max_elem,maxc,nn,xyz,sita,iz,c6ab,mxc,r2r4,r0ab,rcov, &
   &     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
   &     e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn,switch_on,aswitch,bswitch,disper)

!   e6   = e6   *s6

   ! e6abc= e6abc*s6 ! old and wrong !

!   e8   = e8   *s18

!   disp =-e6-e8-e6abc

   print*,'Dispersion (kcal/mol)in PBC= ',disper


   end program D20

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C compute energy
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
   subroutine pbcedisp(max_elem,maxc,n,xyz,laa,iz,c6ab,mxc,r2r4,r0ab,&
  &    rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
  &    e6,e8,e10,e12,e63,lat,rthr,rep_vdw,cn_thr,rep_cn,switch_on,aswitch,bswitch,disper)
   implicit none  
   integer max_elem,maxc
   real*8 r2r4(max_elem),rcov(max_elem)
   real*8 rs6,rs8,rs10,alp6,alp8,alp10
   real*8 rthr,cn_thr,crit_cn
   integer rep_vdw(3),rep_cn(3)
   integer n,iz(*),version,mxc(max_elem)
!   integer rep_v(3)=rep_vdw!,rep_cn(3)
   real*8 xyz(3,*),r0ab(max_elem,max_elem),lat(3,3)!,r2r4(*)
!   real*8 rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
   real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
   real*8 e6, e8, e10, e12, e63!,crit_vdw,crit_cn
   logical noabc,switch_on
!
   Integer, Parameter:: vers=3,nsmax=1000
   Real*8 :: dispfit,x,y,disper,laa,lab, &
   & disp_total,aswitch,bswitch
   Real :: dict1,dict81,beta1,dict2,dict82,beta2
   Integer :: at1,at11,at2,at22,nsa,nsb,i,j
   Dimension :: x(5,nsmax),y(5,nsmax),laa(5,nsmax),lab(5,nsmax)
   logical :: onesystem
!
   integer iat,jat,kat
   real*8 r,r2,r6,r8,tmp,dx,dy,dz,c6,c8,c10,ang,rav,R0
   real*8 damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14
   real*8 cn(n),rxyz(3),dxyz(3)
   real*8 r2ab(n*n),cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,tau(3)
   integer ij,ik,jk !lin
   integer taux,tauy,tauz,counter
   real*8 a1,a2  !BJ-parameter
   real*8 bj_dmp6,bj_dmp8
   real tmp1,tmp2,atta,length,a0
   data a0 / 0.529177249d0/
   !Real :: dict1,dict81,beta1,dict2,dict82,beta2

   x=laa
   e6 =0
   e8 =0
   e10=0
   e12=0
   e63=0
   tau=(/0.0,0.0,0.0/)
   counter=0
   dispfit=0.0
   crit_cn=cn_thr
   atta=0.0
   !rthr=100.0d0
!
!   print*,'laa',laa(4,1)
   !print*,'lat(1,1-3)',lat(1,1),lat(1,2),lat(1,3)
   !print*,'lat(2,1-3)',lat(2,1),lat(2,2),lat(2,3)
   !print*,'lat(3,1-3)',lat(3,1),lat(3,2),lat(3,3)
   Do iat=1,n
      at1=x(4,iat)
      !print*,'x(4,iat),x(:iat)',x(4,iat),x(1,iat)*a0,x(2,iat)*a0,x(3,iat)*a0
      Call Coefficients(x,dict1,dict81,beta1,iat,vers)
      !print*,'iat,dict1,dict81,beta1',iat,dict1,dict81,beta1 
      do jat=iat+1,n
         !length=sqrt((x(1,iat)-x(1,jat))**2+(x(2,iat)-x(2,jat))**2+(x(3,iat)-x(3,jat))**2)
         !if(at1.eq.1.0.and.length.le.3.0) then
           !print*,'length=',length
         !end if
         at2=x(4,jat)
         Call Coefficients(x,dict2,dict82,beta2,jat,vers)
! get C6
!         call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
!        &                              cn(iat),cn(jat),c6)

         rxyz=xyz(:,iat)-xyz(:,jat)
!         r42=r2r4(iz(iat))*r2r4(iz(jat))
!         bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
!         bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

!         if(.not.noabc)then
!           ij=lin(jat,iat)
! store C6 for C9, calc as sqrt
!           cc6ab(ij)=sqrt(c6)
!         endif
!         print*,'rep_vdw',rep_vdw(1),rep_vdw(2),rep_vdw(3)
         do taux=-rep_vdw(1),rep_vdw(1)
         do tauy=-rep_vdw(2),rep_vdw(2)
         do tauz=-rep_vdw(3),rep_vdw(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            
            dxyz=rxyz+tau

            r2=sum(dxyz*dxyz)
! cutoff
            !print*,'r2.gt.rthr=',rthr
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            !if(r.le.3.0d0) cycle
!            print*,'r',r
            dispfit=dispfit + disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
            !if(r.le.3.0d0) then
              !print*,'r,dispfit',r,disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
              !atta=atta+disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
            !end if
!            rr=r0ab(iz(jat),iz(iat))/r

!            r6=r2**3      

!            e6 =e6+c6/(r6+bj_dmp6)

! stored in main as sqrt
!            c8 =3.0d0*c6*r42
!            r8 =r6*r2

!            e8 =e8+c8/(r8+bj_dmp8)

            counter=counter+1

         enddo !tauz
         enddo !tauy
         enddo !taux
         !print*,'dispfit',dispfit
      enddo !jat

! Now the self interaction
      !print*,'at1-at2=',at1,at2
      !print*,'r',r
      !print*,'dispfit-before',dispfit
      jat=iat
      at2=x(4,jat)
      Call Coefficients(x,dict2,dict82,beta2,jat,vers)
! get C6
!      call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
!     &                                  cn(iat),cn(jat),c6)
!      r42=r2r4(iz(iat))*r2r4(iz(iat))
!      bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
!      bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8
           
!      if(.not.noabc)then
!         ij=lin(jat,iat)
! store C6 for C9, calc as sqrt
!         cc6ab(ij)=dsqrt(c6)
!      endif

      do taux=-rep_vdw(1),rep_vdw(1)
      do tauy=-rep_vdw(2),rep_vdw(2)
      do tauz=-rep_vdw(3),rep_vdw(3)
         if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
         tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

         r2=sum(tau*tau)
            !print*,'rthr',rthr
! cutoff
         if(r2.gt.rthr) cycle
         !print*,'r2',r2
         r =sqrt(r2)
         dispfit=dispfit + 0.50d0*disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
         rr=r0ab(iz(jat),iz(iat))/r


!            r6=r2**3      

!            e6 =e6+c6/(r6+bj_dmp6)*0.50d0

! stored in main as sqrt
!            c8 =3.0d0*c6*r42
!            r8 =r6*r2

!            e8 =e8+c8/(r8+bj_dmp8)*0.50d0
         counter=counter+1

      enddo !tauz
      enddo !tauy
      enddo !taux
      !print*,'dispfit-self',dispfit
   enddo !iat
   !print*,'atta',atta
   !disper=dispfit*0.50d0
   !print*,'counter=',counter
   disper=dispfit
   end subroutine pbcedisp
!
   subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

      character(len=*),intent(inout) :: str,before
      character(len=*),intent(in) :: delims
      character,optional :: sep
      logical :: pres
      character :: ch,cha

      pres=present(sep)
      str=adjustl(str)
      call compact(str)
      lenstr=len_trim(str)
      if(lenstr == 0) return        ! string str is empty
      k=0
      ibsl=0                        ! backslash initially inactive
      before=' '
      do i=1,lenstr
         ch=str(i:i)
         if(ibsl == 1) then          ! backslash active
            k=k+1
            before(k:k)=ch
            ibsl=0
            cycle
         end if
         if(ch == '\') then          ! backslash with backslash inactive
            k=k+1
            before(k:k)=ch
            ibsl=1
            cycle
         end if
         ipos=index(delims,ch)         
         if(ipos == 0) then          ! character is not a delimiter
            k=k+1
            before(k:k)=ch
            cycle
         end if
         if(ch /= ' ') then          ! character is a delimiter that is not aspace
            str=str(i+1:)
            if(pres) sep=ch
            exit
         end if
         cha=str(i+1:i+1)            ! character is a space delimiter
         iposa=index(delims,cha)
         if(iposa > 0) then          ! next character is a delimiter
            str=str(i+2:)
            if(pres) sep=cha
            exit
         else
            str=str(i+1:)
            if(pres) sep=ch
            exit
         end if
      end do
      if(i >= lenstr) str=''
      str=adjustl(str)              ! remove initial spaces
      return

   end subroutine split
!
   subroutine compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control
! characters;
! removes initial spaces.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str)):: outstr
      
      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      isp=0
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)
  
        select case(ich)
  
          case(9,32)     ! space or tab character
            if(isp==0) then
              k=k+1
              outstr(k:k)=' '
            end if
            isp=1
            
          case(33:)      ! not a space, quote, or control character
            k=k+1
            outstr(k:k)=ch
            isp=0
      
        end select
        
      end do

      str=adjustl(outstr)

   end subroutine compact

!**********************************************************************

   subroutine removesp(str)

      ! Removes spaces, tabs, and control characters in string str

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str))::outstr

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)
        select case(ich)    
          case(0:32)  ! space, tab, or control character
               cycle       
          case(33:)  
            k=k+1
            outstr(k:k)=ch
        end select
      end do
      
      str=adjustl(outstr)
      
   end subroutine removesp
!
   subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str))::outstr

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      k=0
      ibsl=0                        ! backslash initially inactive
      
      do i=1,lenstr
        ch=str(i:i)
        if(ibsl == 1) then          ! backslash active
         k=k+1
         outstr(k:k)=ch
         ibsl=0
         cycle
        end if
        if(ch == '\') then          ! backslash with backslash inactive
         ibsl=1
         cycle
        end if
        k=k+1
        outstr(k:k)=ch              ! non-backslash with backslash inactive
      end do
      
      str=adjustl(outstr)
      
   end subroutine removebksl



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       
!c            string pars procedures
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.
   interface 
       subroutine split(str,delims,before,sep)
         character(len=*),intent(inout) :: str,before
         character(len=*),intent(in) :: delims
         character,optional,intent(inout) :: sep
       end subroutine split
   end interface

      character(len=*),intent(inout) :: str
      character(len=*),intent(in) :: delims
      character(len=len_trim(str)) :: strsav
      character(len=*),dimension(:),intent(inout) :: args
      integer, intent(out) :: nargs
      
      strsav=str
      call compact(str)
      na=size(args)
      do i=1,na
        args(i)=' '
      end do  
      nargs=0
      lenstr=len_trim(str)
      if(lenstr==0) return
      k=0

      do
         if(len_trim(str) == 0) exit
         nargs=nargs+1
         call split(str,delims,args(nargs))
         call removebksl(args(nargs))
      end do   
      str=strsav

   end subroutine parse




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c read coordinates in Angst and converts them to au 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   subroutine pbcrdcoord(fname,lattice,n,xyz,iat,autoang,sita)
   implicit none             
   interface
     subroutine parse(str,delims,args,nargs)
     character(len=*),intent(inout) :: str
     character(len=*),intent(in)  :: delims
     character(len=*),dimension(:),intent(inout) :: args
     integer, intent(out) :: nargs
     end subroutine parse
   end interface
   
   Real*8, Dimension(:,:):: sita(5,1000)
   real*8                :: xyz(3,*),atmnum
   real*8, INTENT(OUT)   ::lattice(3,3)
   integer, INTENT(out)               :: iat(*) 
   integer, INTENT(in)               :: n 
   character*(*), INTENT(IN)          :: fname
   logical              :: selective=.FALSE. ! Selective dynamics
   logical              :: cartesian=.TRUE.  ! Cartesian or direct
   real*8, INTENT(IN)   ::autoang

   real*8 xx(10),scalar,a0
   character*200 line,bb
   character*80 args(90),args2(90)
   data a0 / 0.529177249d0/
   integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck


   lattice=0
   
   ich=142
   open(unit=ich,file=fname)
   rewind(ich)
   ncheck=0
   ntype=0
   read(ich,'(a)',end=200)line !first line must contain Element Info
   call parse(line,' ',args,ntype)
   read(ich,'(a)',end=200)line !second line contains global scaling factor
   call readl(line,xx,nn)
   scalar=xx(1)/autoang        !the Ang->au conversion is included in the scaling factor
!  write(*,'(F8.6)')scalar
   DO i=1,3            ! reading the lattice constants
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn < 3) call stoprun( 'Error reading unit cell vectors' )
     lattice(1,i)=xx(1)*scalar
     lattice(2,i)=xx(2)*scalar
     lattice(3,i)=xx(3)*scalar
   !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
   ENDDO
   read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
   line=adjustl(line)
   call readl(line,xx,nn)
   IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
     call parse(line,' ',args,ntype)
     read(ich,'(a)',end=200)line
     line=adjustl(line)
     call readl(line,xx,nn)
   ENDIF
!       call elem(args(1),i_dummy2)
!       IF (i_dummy2<1 .OR. i_dummy2>94) THEN
!          args=args2
!       ENDIF
   IF (nn.NE.ntype ) THEN
     call stoprun( 'Error reading number of atomtypes')
   ENDIF
   ncheck=0
   DO i=1,nn
     i_dummy1=INT(xx(i))
     call elem(args(i),i_dummy2, atmnum)
!     print*,atmnum,i_dummy2
     IF (i_dummy2<1 .OR. i_dummy2>94) then 
        call stoprun( 'Error: unknown element.')
     end IF
     DO j=1,i_dummy1
       ncheck=ncheck+1
       sita(4,ncheck)=atmnum
       iat(ncheck)=i_dummy2
     ENDDO
   ENDDO
   if (n.ne.ncheck) call stoprun('Error reading Number of Atoms')

   read(ich,'(a)',end=200)line
   line=adjustl(line)
   IF (line(:1).EQ.'s' .OR. line(:1).EQ.'S') THEN
     selective=.TRUE.
     read(ich,'(a)',end=200)line
     line=adjustl(line)
   ENDIF
!   print*,'n',n
!   write(*,*)line(:1)
   cartesian=(line(:1).EQ.'c' .OR. line(:1).EQ.'C' .OR. line(:1).EQ.'k' .OR. line(:1).EQ.'K')
   !Open(unit=21,file='dimer1.xyz')
   DO i=1,n
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn.NE.3) call stoprun( 'Error reading coordinates.')

     IF (cartesian) THEN
       xyz(1,i)=xx(1)*scalar
       xyz(2,i)=xx(2)*scalar
       xyz(3,i)=xx(3)*scalar
     ELSE
       xyz(1,i)=lattice(1,1)*xx(1)+lattice(1,2)*xx(2)+lattice(1,3)*xx(3)
       xyz(2,i)=lattice(2,1)*xx(1)+lattice(2,2)*xx(2)+lattice(2,3)*xx(3)
       xyz(3,i)=lattice(3,1)*xx(1)+lattice(3,2)*xx(2)+lattice(3,3)*xx(3)
     ENDIF
     !read(21,*) bb,(sita(j,i),j=1,3)
     !sita(1,i)=sita(1,i)/a0
     !sita(2,i)=sita(2,i)/a0
     !sita(3,i)=sita(3,i)/a0
     sita(1,i)=xyz(1,i)
     sita(2,i)=xyz(2,i)
     sita(3,i)=xyz(3,i)
     !print*,sita(1,i)*a0,sita(2,i)*a0,sita(3,i)*a0 !,sita(4,i)
!      write(*,'(3F20.10,1X,I3)')xyz(:,i),iat(i)   !debug printout
      
   ENDDO
   !close(21)
!   print*,sita(1,38),sita(2,38),sita(3,38)
      
      
 200  continue

   close(ich)
   end subroutine pbcrdcoord


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Returns the number of a given element string (h-pu, 1-94)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   SUBROUTINE ELEM(KEY1, NAT, atmnum)
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   CHARACTER*(*) KEY1
   CHARACTER*2 ELEMNT(94),E
   Real*8 :: atmnum

   DATA ELEMNT/'h ','he',&
  & 'li','be','b ','c ','n ','o ','f ','ne',&
  & 'na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr',&
  & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
  & 'cd','in','sn','sb','te','i ','xe',&
  & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
  & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
  & 'au','hg','tl','pb','bi','po','at','rn',&
  & 'fr','ra','ac','th','pa','u ','np','pu'/
!   print*,key1
   if (key1.eq.'H' .or. key1.eq.'h') then
     atmnum=1.00
   elseif (key1.eq.'He' .or. key1.eq.'HE' .or. key1.eq.'he') then
     atmnum=2.00
   elseif (key1.eq.'Li' .or. key1.eq.'LI' .or. key1.eq.'li') then
     atmnum=3.00
   elseif (key1.eq.'B' .or. key1.eq.'b') then
     atmnum=5.00
   elseif (key1.eq.'C' .or. key1.eq.'c') then
     atmnum=6.00
   elseif (key1.eq.'N' .or. key1.eq.'n') then
     atmnum=7.00
   elseif (key1.eq.'O' .or. key1.eq.'o') then
     atmnum=8.00
   elseif (key1.eq.'F' .or. key1.eq.'f') then
     atmnum=9.00
   elseif (key1.eq.'Ne' .or. key1.eq.'NE' .or. key1.eq.'ne') then
     atmnum=10.00
   elseif (key1.eq.'Al' .or. key1.eq.'AL' .or. key1.eq.'al') then
     atmnum=13.00
   elseif (key1.eq.'Si' .or. key1.eq.'SI' .or. key1.eq.'si') then
     atmnum=14.00
   elseif (key1.eq.'P' .or. key1.eq.'p') then
     atmnum=15.00
   elseif (key1.eq.'S' .or. key1.eq.'s') then
     atmnum=16.00
   elseif (key1.eq.'Cl' .or. key1.eq.'CL' .or. key1.eq.'cl') then
     atmnum=17.00
   elseif (key1.eq.'Ar' .or. key1.eq.'AR' .or. key1.eq.'ar') then
     atmnum=18.00
   elseif (key1.eq.'Br' .or. key1.eq.'BR' .or. key1.eq.'br') then
     atmnum=35.00
   elseif (key1.eq.'I' .or. key1.eq.'i') then
     atmnum=53.00
   end if
!
   nat=0
   e='  '
   k=1
   DO J=1,len(key1)
         if (k.gt.2)exit       
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
            k=k+1
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
   enddo

   DO I=1,94
      if(e.eq.elemnt(i))then
         NAT=I
         RETURN
      ENDIF
   ENDDO
!   print*,key1
!   if (key1.eq.'H' .or. key1.eq.'h') then
!     atmnum=1.00
!   elseif (key1.eq.'He' .or. key1.eq.'HE' .or. key1.eq.'he') then
!     atmnum=2.00
!   elseif (key1.eq.'Li' .or. key1.eq.'LI' .or. key1.eq.'li') then
!     atmnum=3.00
!   elseif (key1.eq.'B' .or. key1.eq.'b') then
!     atmnum=5.00
!   elseif (key1.eq.'C' .or. key1.eq.'c') then
!     atmnum=6.00
!   elseif (key1.eq.'N' .or. key1.eq.'n') then
!     atmnum=7.00
!   elseif (trim(adjustl(key1)).eq.'O' .or. key1.eq.'o') then
!     atmnum=8.00
!   elseif (key1.eq.'F' .or. key1.eq.'f') then
!     atmnum=9.00
!   elseif (key1.eq.'Ne' .or. key1.eq.'NE' .or. key1.eq.'ne') then
!     atmnum=10.00
!   elseif (key1.eq.'Al' .or. key1.eq.'AL' .or. key1.eq.'al') then
!     atmnum=13.00
!   elseif (key1.eq.'Si' .or. key1.eq.'SI' .or. key1.eq.'si') then
!     atmnum=14.00
!   elseif (key1.eq.'P' .or. key1.eq.'p') then
!     atmnum=15.00
!   elseif (key1.eq.'S' .or. key1.eq.'s') then
!     atmnum=16.00
!   elseif (key1.eq.'Cl' .or. key1.eq.'CL' .or. key1.eq.'cl') then
!     atmnum=17.00
!   elseif (key1.eq.'Ar' .or. key1.eq.'AR' .or. key1.eq.'ar') then
!     atmnum=18.00
!   elseif (key1.eq.'Br' .or. key1.eq.'BR' .or. key1.eq.'br') then
!     atmnum=35.00
!   elseif (key1.eq.'I' .or. key1.eq.'i') then
!     atmnum=53.00
!   end if
    
   end SUBROUTINE ELEM

!C     *****************************************************************         

   FUNCTION ESYM(I)
   CHARACTER*2 ESYM
   CHARACTER*2 ELEMNT(94)
   DATA ELEMNT/'h ','he',&
  & 'li','be','b ','c ','n ','o ','f ','ne',&
  & 'na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr',&
  & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
  & 'cd','in','sn','sb','te','i ','xe',&
  & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
  & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
  & 'au','hg','tl','pb','bi','po','at','rn',&
  & 'fr','ra','ac','th','pa','u ','np','pu'/
   ESYM=ELEMNT(I)
   RETURN
   END FUNCTION ESYM



   SUBROUTINE SET_CRITERIA(rthr,lat,tau_max)

     REAL*8 :: r_cutoff,rthr
     REAL*8 :: lat(3,3)
     REAL*8 :: tau_max(3)
     REAL*8 :: norm1(3),norm2(3),norm3(3)
     REAL*8 :: cos10,cos21,cos32
     real*8,external :: vectorsize

     r_cutoff=sqrt(rthr)
!       write(*,*) 'lat',lat
       !c find normal to the plane...
       ! kreuzprodukt means cross product
     call kreuzprodukt(lat(:,2),lat(:,3),norm1)
     call kreuzprodukt(lat(:,3),lat(:,1),norm2)
     call kreuzprodukt(lat(:,1),lat(:,2),norm3)
!     write(*,*) 'norm2',norm2
     !c ...normalize it...
     norm1=norm1/VECTORSIZE(norm1)
     norm2=norm2/VECTORSIZE(norm2)
     norm3=norm3/VECTORSIZE(norm3)
!     write(*,*) 'norm2_',norm2
       !c cos angles between normals and lattice vectors
     cos10=SUM(norm1*lat(:,1))
     cos21=SUM(norm2*lat(:,2))
     cos32=SUM(norm3*lat(:,3))
       !write(*,*) 'cos32',cos32
       !tau_max(1)=abs(2*r_cutoff/cos10)
       !tau_max(2)=abs(2*r_cutoff/cos21)
       !tau_max(3)=abs(2*r_cutoff/cos32)
       !write(*,*) 'r_cutoff',r_cutoff
     tau_max(1)=abs(r_cutoff/cos10)
     tau_max(2)=abs(r_cutoff/cos21)
     tau_max(3)=abs(r_cutoff/cos32)
!     write(*,'(3f8.4)')tau_max(1),tau_max(2),tau_max(3)
   END SUBROUTINE SET_CRITERIA


   SUBROUTINE kreuzprodukt(A,B,C)
     IMPLICIT NONE
  
     REAL*8 :: A(3),B(3)
     REAL*8 :: X,Y,Z
     REAL*8 :: C(3)
     
     X=A(2)*B(3)-B(2)*A(3)
     Y=A(3)*B(1)-B(3)*A(1)
     Z=A(1)*B(2)-B(1)*A(2)
     C=(/X,Y,Z/)
   END SUBROUTINE kreuzprodukt

   FUNCTION VECTORSIZE(VECT)

      REAL*8 :: VECT(3)
      REAL*8 :: SVECT(3)
      REAL*8 :: VECTORSIZE

      SVECT=VECT*VECT
      VECTORSIZE=SUM(SVECT)
      VECTORSIZE=VECTORSIZE**(0.5)
   END FUNCTION VECTORSIZE



  
!  These subroutines are copied and modified from Grimme D3 dispersion function
!  This subroutine reads POSCAR (of VASP package) file:
   subroutine pbcrdatomnumber(fname,n)
   implicit none
   interface
     subroutine parse(str,delims,args,nargs)
     character(len=*),intent(inout) :: str
     character(len=*),intent(in)  :: delims
     character(len=*),dimension(:),intent(inout) :: args
     integer, intent(out) :: nargs
     end subroutine parse
   end interface

   integer, INTENT(out)               :: n
   character*(*), INTENT(IN)          :: fname
   logical              :: selective=.FALSE. ! Selective dynamics
   logical              :: cartesian=.TRUE.  ! Cartesian or direct

   real*8 xx(10),scalar,fdum
   character*80 line,args(90),args2(90)

   integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2

   ich=142
   open(unit=ich,file=fname)
   n=0
   ntype=0
   read(ich,'(a)',end=200)line !first line must contain Element Info
   call parse(line,' ',args,ntype)
   read(ich,'(a)',end=200)line !second line contains global scaling factor
   call readl(line,xx,nn)
!   write(*,'(F8.6)')scalar
   DO i=1,3            ! reading the lattice constants
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn < 3) call stoprun( 'Error reading unit cell vectors' )
   !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
   ENDDO
   read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
   line=adjustl(line)
   call readl(line,xx,nn)
   IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
     call parse(line,' ',args,ntype)
     read(ich,'(a)',end=200)line
     line=adjustl(line)
     call readl(line,xx,nn)
   ENDIF
!    call elem(args(1),i_dummy2)
!    IF (i_dummy2<1 .OR. i_dummy2>94) THEN
!       args=args2
!    ENDIF
   IF (nn.NE.ntype ) THEN
!      IF(nn.NE.ntype2) THEN
     call stoprun( 'Error reading number of atomtypes')
!      ELSE
!        ntype=ntype2
!      ENDIF
   ENDIF
   n=0
   DO i=1,nn
     i_dummy1=INT(xx(i))
       n=n+i_dummy1
   ENDDO

 200  continue

   close(ich)
   end subroutine pbcrdatomnumber
!
!     *****************************************************************
!     Reads a given line
!     *****************************************************************

   subroutine readl(A1,X,N)
   IMPLICIT REAL*8 (A-H,O-Z)
   CHARACTER*(*) A1
   DIMENSION X(*)
   I=0
   IS=1
10 I=I+1
   X(I)=READAA(A1,IS,IB,IE)
   IF(IB.GT.0 .AND. IE.GT.0) THEN
           IS=IE
           GOTO 10
   ENDIF
   N=I-1
   RETURN
   end subroutine readl

   FUNCTION READAA(A,ISTART,IEND,IEND2)                                   
   IMPLICIT REAL*8 (A-H,O-Z)                                                 
   REAL*8 READAA                                                             
   CHARACTER*(*) A                                                      
   NINE=ICHAR('9')                                                           
   IZERO=ICHAR('0')                                                          
   MINUS=ICHAR('-')                                                          
   IDOT=ICHAR('.')                                                           
   ND=ICHAR('D')                                                             
   NE=ICHAR('E')                                                             
   IBL=ICHAR(' ')                                                            
   IEND=0                                                                    
   IEND2=0                                                                   
   IDIG=0                                                                    
   C1=0                                                                      
   C2=0                                                                      
   ONE=1.D0                                                                  
   X = 1.D0                                                                  
   NL=LEN(A) 
   DO 10 J=ISTART,NL-1                                                       
      N=ICHAR(A(J:J))                                                          
      M=ICHAR(A(J+1:J+1)) 
      IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
      IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20                                                 
10 CONTINUE                                                                  
   READAA=0.D0                                                               
   RETURN                                                                    
20 CONTINUE                                                                  
   IEND=J                                                                    
   DO 30 I=J,NL                                                              
         N=ICHAR(A(I:I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
30 CONTINUE                                                                  
40 CONTINUE                                                                  
   IDIG=0                                                                    
   DO 50 II=I+1,NL                                                           
         N=ICHAR(A(II:II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
50 CONTINUE                                                                  
!                                                                               
! PUT THE PIECES TOGETHER                                                       
!                                                                               
60 CONTINUE                                                                  
   READAA= ONE * ( C1 + C2 * X)                                              
   DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J:J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
   RETURN                                                                    
                                                                                
57 C1=0.0D0                                                                  
   ONE=1.0D0                                                                 
   DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I:I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
31 CONTINUE                                                                  
61 CONTINUE                                                                  
70 READAA=READAA*10**(ONE*C1)                                                
   RETURN                                                                    
   END FUNCTION READAA




   subroutine stoprun(s)
   character*(*) s
   write(*,*)'program stopped due to: ',s
   call system('touch dscf_problem')
   stop 'must stop!'
   end subroutine stoprun
!
!

!
!
   Subroutine Coefficients(sita,dict,dict8,beta,i,vers)

   parameter (nsmax=100,conv=0.529177209)
   Real*8 :: a,b,sita
   Real, Intent(out) :: dict,dict8,beta
   Integer :: i,vers
   dimension :: sita(5,nsmax),sitb(5,nsmax),sitout(5,nsmax),sitbb(5,nsmax),&
   & ama(nsmax),amb(nsmax),dicta(100),dicta8(100),betaa(100)
   Data dicta/ 0.199730d0,0.625058d0,0.216319d0,1.5290d0,0.117956d0,0.155250d0,0.028773d0,&
         &  0.984668d0,0.957855d0,0.338996d0,0.405991d0,0.09588d0,0.185489d0,&
         &  0.185411d0,0.071774d0,0.003104d0,1.297139d0,0.914964d0,0.461576d0,&
         &  1.427069d0,0.72359d0,0.900101d0,0.261125d0,0.415583d0,0.009843d0,&
         &  0.890484d0,4.793738d0,6.254113d0,3.33151d0,9.995193d0,21.167203d0,69*0/
   Data dicta8/ 0.010225d0,0.004286d0,0.002841d0,0.1649d0,0.012842d0,0.003195d0,0.005962d0,&
         &  0.060477d0,0.011257d0,0.29648d0,0.011234d0,0.008433d0,0.020926d0,&
         &  0.045112d0,0.004783d0,0.264031d0,0.171713d0,0.268786d0,0.235496d0,&
         &  0.070811d0,0.052314d0,0.015629d0,0.028508d0,0.900131d0,0.894086d0,&
         &  2.415552d0,0.825826d0,0.691238d0,0.509268d0,1.169732d0,2.745814d0,69*0/
   Data betaa/ 1.898977d0,2.022739d0,1.737714d0,1.9509d0,1.456906d0,1.623847d0,1.618463d0,&
         &    1.828202d0,1.62796d0,0.589819d0,1.416878d0,1.66046d0,1.386093d0,&
         &    1.310652d0,2.225149d0,1.695239d0,1.951721d0,1.65896d0,1.780635d0,&
         &    2.396253d0,2.368908d0,2.458106d0,2.13224d0,1.133178d0,2.726511d0,&
         &    1.634346d0,2.279597d0,1.643112d0,1.843843d0,1.572187d0,1.463387d0,69*0/
   a=sita(4,i)
   b=sita(5,i)
   !print*,'vers',vers
!   print*,'vers',vers
!  version 3 corresponds to Das_2020
   if (vers.eq.3) then
      if (a.eq.1.0.and.b.eq.1.0) then
         dict=dicta(1)
         dict8=dicta8(1)
         beta=betaa(1)
      elseif (a.eq.1.0.and.b.eq.5.0) then
         dict=dicta(2)
         dict8=dicta8(2)
         beta=betaa(2)
      elseif (a.eq.1.0.and.b.eq.6.0) then
         dict=dicta(3)
         dict8=dicta8(3)
         beta=betaa(3)
      elseif (a.eq.6.0.and.b.eq.1.0) then
         dict=dicta(4)
         dict8=dicta8(4)
         beta=betaa(4)
      elseif (a.eq.1.0.and.b.eq.7.0) then
         dict=dicta(5)
         dict8=dicta8(5)
         beta=betaa(5)
      elseif (a.eq.1.0.and.b.eq.8.0) then
         dict=dicta(6)
         dict8=dicta8(6)
         beta=betaa(6)
      elseif (a.eq.1.0.and.b.eq.9.0) then
         dict=dicta(7)
         dict8=dicta8(7)
         beta=betaa(7)
      elseif (a.eq.1.0.and.b.eq.13.0) then
         dict=dicta(8)
         dict8=dicta8(8)
         beta=betaa(8)
      elseif (a.eq.1.0.and.b.eq.14.0) then
         dict=dicta(9)
         dict8=dicta8(9)
         beta=betaa(9)
      elseif (a.eq.1.0.and.b.eq.15.0) then
         dict=dicta(10)
         dict8=dicta8(10)
         beta=betaa(10)
      elseif (a.eq.1.0.and.b.eq.16.0) then
         dict=dicta(11)
         dict8=dicta8(11)
         beta=betaa(11)
      elseif (a.eq.1.0.and.b.eq.17.0) then
         dict=dicta(12)
         dict8=dicta8(12)
         beta=betaa(12)
      elseif (a.eq.1.0.and.b.eq.35.0) then
         dict=dicta(13)
         dict8=dicta8(13)
         beta=betaa(13)
      elseif (a.eq.1.0.and.b.eq.53.0) then
         dict=dicta(14)
         dict8=dicta8(14)
         beta=betaa(14)
      elseif (a.eq.2.0.and.b.eq.0.0) then
         dict=dicta(15)
         dict8=dicta8(15)
         beta=betaa(15)
      elseif (a.eq.5.0.and.b.eq.0.0) then
         dict=dicta(16)
         dict8=dicta8(16)
         beta=betaa(16)
      elseif (a.eq.6.0.and.b.eq.2.0) then
         dict=dicta(17)
         dict8=dicta8(17)
         beta=betaa(17)
      elseif (a.eq.6.0.and.b.eq.3.0) then
         dict=dicta(18)
         dict8=dicta8(18)
         beta=betaa(18)
      elseif (a.eq.6.0.and.b.eq.4.0) then
         dict=dicta(19)
         dict8=dicta8(19)
         beta=betaa(19)
      elseif (a.eq.7.0.and.b.eq.0.0) then
         dict=dicta(20)
         dict8=dicta8(20)
         beta=betaa(20)
      elseif (a.eq.8.0.and.b.eq.0.0) then
         dict=dicta(21)
         dict8=dicta8(21)
         beta=betaa(21)
      elseif (a.eq.9.0.and.b.eq.0.0) then
         dict=dicta(22)
         dict8=dicta8(22)
         beta=betaa(22)
      elseif (a.eq.10.0.and.b.eq.0.0) then
         dict=dicta(23)
         dict8=dicta8(23)
         beta=betaa(23)
      elseif (a.eq.13.0.and.b.eq.0.0) then
         dict=dicta(24)
         dict8=dicta8(24)
         beta=betaa(24)
      elseif (a.eq.14.0.and.b.eq.0.0) then
         dict=dicta(25)
         dict8=dicta8(25)
         beta=betaa(25)
      elseif (a.eq.15.0.and.b.eq.0.0) then
         dict=dicta(26)
         dict8=dicta8(26)
         beta=betaa(26)
      elseif (a.eq.16.0.and.b.eq.0.0) then
         dict=dicta(27)
         dict8=dicta8(27)
         beta=betaa(27)
      elseif (a.eq.17.0.and.b.eq.0.0) then
         dict=dicta(28)
         dict8=dicta8(28)
         beta=betaa(28)
      elseif (a.eq.18.0.and.b.eq.0.0) then
         dict=dicta(29)
         dict8=dicta8(29)
         beta=betaa(29)
      elseif (a.eq.35.0.and.b.eq.0.0) then
         dict=dicta(30)
         dict8=dicta8(30)
         beta=betaa(30)
      elseif (a.eq.53.0.and.b.eq.0.0) then
         dict=dicta(31)
         dict8=dicta8(31)
         beta=betaa(31)
      end if      
   end if
   !print*,'dict-dict8-beta',dict,dict8,beta
   end Subroutine

   Function catom(dict)
!  """Returns C6 coefficient for atom at, converted to kcal*bohr^6/mol."""
   Real :: dict,catom

   catom=dict*10884.3258684610

   End Function

   Function catomatom(dicta,dictb)

   Real*8 :: c6a,c6b,catomatom
   Real :: dicta,dictb

   c6a = catom(dicta)
   c6b = catom(dictb)
   catomatom = sqrt(c6a*c6b)

   End Function catomatom 


   Function catom8(dict8)
!  """Returns C8 coefficient for atom at, converted to kcal*bohr^8/mol."""
   Real*8 :: catom8
   Real :: dict8

   catom8=dict8*3886863.35457252

   End Function

   Function catomatom8(dict8a,dict8b)
!  """Returns C8 coefficient between atoms with charges at1 and at2."""
   Real*8 :: c8a,c8b,catomatom8
   Real :: dict8a,dict8b
   c8a=catom8(dict8a)
   c8b=catom8(dict8b)

   catomatom8=sqrt(c8a*c8b)

   End Function

   Function damp_TT(r,beta1,beta2)
!   """Returns the value of Tang-Toennies damping function at separation r."""
!  it is assumed below that the damping factor is the geometric mean of atomic
!  damping factors:
   Real*8 :: alpha,r,br,damp_TT,sm,term
   Real :: beta1,beta2

   alpha=sqrt(abs(beta1*beta2))
   br=alpha*r
!   print*,alpha
   sm=1.0
   term=1.0
   Do i=1,6
     term=term*br/i
     sm=sm+term
   end do
   damp_TT=1.0 - Exp(-br)*sm
   
   End Function


   Function damp_TT8(r,beta1,beta2)
!  """Returns the value of Tang-Toennies damping function f_8 at separation r."""
! it is assumed below that the damping factor is the geometric mean of atomic
! damping factors:
   Real*8 :: alpha,sm,br,damp_TT8,r,term
   Real :: beta1,beta2
!   Real :: beta1,beta2
!   print*,beta1,beta2
   alpha=sqrt(abs(beta1*beta2))
   br=alpha*r
!   print*,beta1,beta2
!   print*,alpha
   sm=1.0
   term=1.0
   Do i=1,8
     term=term*br/i
     sm=sm+term
   end do
!   print*,sm
   damp_TT8=1.0 - Exp(-br)*sm

   End Function
 

   Function covalent_radius(att)
!  """Returns the covalent (Bragg) radius of the atom at, converted to bohr.
!  The radii are taken from the Cambridge Structural Database."""
!  print('at[0]=',at)
   Real*8 :: r,covalent_radius
   Integer :: att,at
    at=att
!    print*,'covalent_radius=',at
!    at=att
    if (at==2) then
      r=1.50d0
    elseif (at==3) then
      r=1.28d0
    elseif (at==4) then
      r=0.96d0
    elseif (at==5) then
      r=0.83d0
    elseif (at==6) then
      r=0.68d0
    elseif (at==7) then
      r=0.68d0
    elseif (at==8) then
      r=0.68d0
    elseif (at==9) then
      r=0.64d0
    elseif (at==10) then
      r=1.50d0
    elseif (at==11) then
      r=1.66d0
    elseif (at==12) then
      r=1.41d0
    elseif (at==13) then
      r=1.21d0
    elseif (at==14) then
      r=1.20d0
    elseif (at==15) then
      r=1.05d0
    elseif (at==16) then
      r=1.02d0
    elseif (at==17) then
      r=0.99d0
    elseif (at==18) then
      r=1.51d0
    elseif (at==35) then
      r=1.21d0
    elseif (at==53) then
      r=1.40d0
    elseif (at==1) then
      r=0.23d0
!    else
!      print*,'Wrong atom type in covalent_radius'
    end if
    covalent_radius=r/0.529177209
   End Function

   Function switching_function(r,at1,at2,aswitch,bswitch)
!  """Returns the value of the dispersion switching (Fermi) function for
!  atoms at1,at2 at separation r."""
   Real*8 :: r,aswitch,bswitch,r0,switching_function
   Integer*8 :: at1,at2
!   dimension(5,nsmax) :: at1,at2
!   print*,covalent_radius(at1),covalent_radius(at2)
   r0=aswitch*(covalent_radius(at1)+covalent_radius(at2))
   switching_function=1.0/(1.0 + Exp(-bswitch*(r/r0-1.0)))

   End Function


   Function disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)

!  """Returns the overall contribution to dispersion
!  energy (in kcal/mol) for atoms at1,at2 at separation r."""
!  if not(at1 in dict.keys()):
!    print ("Dispersion coefficients not defined for atom",at1)
!    raise "Error!"
!  if not(at2 in dict.keys()):
!    print ("Dispersion coefficients not defined for atom",at2)
!    raise "Error!"

   Real*8 :: disp_total,edisp,r,aswitch,bswitch, &
   & damp_TT
   Real :: dict1,dict81,beta1,dict2,dict82,beta2,e6,e8,disp
   Logical :: switch_on
   Integer :: at1,at2
!   dimension(5,nsmax) :: at1,at2
   !print*,dict1,dict2,beta1,dict81,dict82,beta2
   edisp=-catomatom(dict1,dict2)/(r*r*r*r*r*r)*damp_TT(r,beta1,beta2) - &
   &      catomatom8(dict81,dict82)/(r*r*r*r*r*r*r*r)*damp_TT8(r,beta1,beta2)
   !e6=-catomatom(dict1,dict2)/(r*r*r*r*r*r)*damp_TT(r,beta1,beta2)
   !e8=-catomatom8(dict81,dict82)/(r*r*r*r*r*r*r*r)*damp_TT8(r,beta1,beta2)
   !if (abs(e8).ge.abs(e6)) then
   !   print*,'r',r
   !   print*,'e6-e8',e6,e8
   !   print*,'disp=',(e6+e8)*switching_function(r,at1,at2,aswitch,bswitch)
   !end if
   !print*,'e6=',-catomatom(dict1,dict2)/(r*r*r*r*r*r)*damp_TT(r,beta1,beta2)
   !print*,'e8=',-catomatom8(dict81,dict82)/(r*r*r*r*r*r*r*r)*damp_TT8(r,beta1,beta2)
!   print*,'edisp=',edisp
!   print*,damp_TT(r,beta1,beta2)
!   if (switch_on) then

   disp_total=edisp*switching_function(r,at1,at2,aswitch,bswitch)
   !print*,'disp_total',disp_total
!   end if

   End Function



   Subroutine fitted_dispersion(laa,lab,nsa,nsb,switch_on,aswitch,bswitch,disper)
!  """Calculates the fitted dispersion energy, either within a molecule
!     or between two molecules."""
!  print('laa-fitted_dispersion',laa)
!  print('lab-fitted_dispersion',lab)
   Parameter (vers=3,nsmax=100)
   Real*8 :: dispfit,x,y,disper,laa,lab,r, &
   & disp_total
   Real :: dict1,dict81,beta1,dict2,dict82,beta2
   Integer :: at1,at11,at2,at22,nsa,nsb,i,j
   Dimension :: x(5,nsmax),y(5,nsmax),laa(5,nsmax),lab(5,nsmax)
   logical :: onesystem
   onesystem=.true.
!   onesystem=.false.
!   intent(in) laa,lab

!   print*,lab(:,1)
   x=laa
   y=lab
   dispfit=0.0
   if(onesystem) then
      Do i = 1,nsa
         at1=x(4,i)
         Call Coefficients(x,dict1,dict81,beta1,i,vers)
         Do j = 1,nsa
            at2=x(4,j)
            Call Coefficients(x,dict2,dict82,beta2,j,vers)
            r=sqrt((x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2 + (x(3,i)-x(3,j))**2)
            if(r.gt.0.0) then
               !print*,r
               !print*,dict1,dict81
               dispfit=dispfit + disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
               !print*,disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
            else
               !print*,r 
               !dispfit=dispfit+0.0
               cycle
            end if
         end Do
      end Do
!   if onesystem
!     for x in laa
!       for y in laa
!          r=math.sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]))
!          if r>0.0
!            dispfit+=disp_total(r,x[3],y[3],dict,dict8,dicd,switch_on,aswitch,bswitch)
!     dispfit*=0.5
!   elseif
!  if(twosystem)
   else
    Do i = 1,nsa
      at1=x(4,i)
!      print*,x(:,i)
      Call Coefficients(x,dict1,dict81,beta1,i,vers)
      Do j = 1,nsb
         at2=y(4,j)
!         
         Call Coefficients(y,dict2,dict82,beta2,j,vers)
         r=sqrt((x(1,i)-y(1,j))**2 + (x(2,i)-y(2,j))**2 + (x(3,i)-y(3,j))**2)
!         print*,'r =',r
         dispfit=dispfit + disp_total(r,at1,at2,dict1,dict81,beta1,dict2,dict82,beta2,switch_on,aswitch,bswitch)
      end Do
!   print*,'dispfit ',dispfit
    end Do
   end if
   disper=dispfit/2.0
   
   End subroutine 



   
   Subroutine analyze_carbons(latom,nat,lat,rep_cn,cn_thr,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=1000,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr,cn_thr,tau(3),lat(3,3),dx,dy,dz
   Integer :: nat,typ,rep_cn(3),taux,tauy,tauz
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   Real*8,dimension(:) :: sita(5,nsmax),sitb(5,nsmax),sitout(5,nsmax),sitbb(5,nsmax),&
   & ama(nsmax),amb(nsmax)
!   Integer, Intent(out) :: a
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(5,nsmax) :: latom(5,nsmax),atom(5,nsmax)

   !print*,'cn_thr: ',cn_thr
   !print*,rep_cn(1),rep_cn(2),rep_cn(3)
   Do i=1,nat
      atom=latom
      if (atom(4,i).ne.6.0) then
         cycle
      end if
      nneighbors=0
      Do j=1,nat
        do taux=-rep_cn(1),rep_cn(1)
        do tauy=-rep_cn(2),rep_cn(2)
        do tauz=-rep_cn(3),rep_cn(3)
          if(j.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
          dx=atom(1,j)-atom(1,i)+tau(1)
          dy=atom(2,j)-atom(2,i)+tau(2)
          dz=atom(3,j)-atom(3,i)+tau(3)
          r=(dx*dx+dy*dy+dz*dz)
          if (r.gt.cn_thr) cycle
          r=sqrt(r)
!         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
!         & (atom(3,i)-latom(3,j))**2)
          if (r.eq.0.0) then
             cycle
          elseif (r.lt.3.0) then
             nneighbors=nneighbors+1
          end if
        enddo !tauz
        enddo !tauy
        enddo !taux
      end Do ! j
      sitout(5,i)=nneighbors
   end Do ! i

   end subroutine



   subroutine analyze_hydrogens(latom,nat,lat,rep_cn,crit_cn,sitout)

   implicit real*8 (a-h,o-z)
   parameter (nsmax=1000,conv=0.529177209)
   character(128) buffer,f
   Real*8 :: mass,cha,t,r,closestr,crit_cn,dx,dy,dz
   Integer :: nat,typ,rep_cn(3),taux,tauy,tauz
   character*4 doit,done,status,doitr6,dummy,N,C,O,H
   character*3 itypa(nsmax),itypb(nsmax),itypc(nsmax)
   character(200) :: filename,file,fileB
   LOGICAL printed,eqv,ok
   Real*8,dimension(:) :: sita(5,nsmax),sitb(5,nsmax),sitout(5,nsmax),sitbb(5,nsmax),&
   & ama(nsmax),amb(nsmax),lat(3,3),tau(3)
   dimension aol(3), bol(3),tran(3,3),com(3)
   Real, dimension(nsmax,nsmax) :: ltempa,ltempb
   Real*8, dimension(5,nsmax) :: latom(5,nsmax),atom(5,nsmax)

!   print*,'latom',latom(1,1),latom(2,1),latom(3,1)
!   sita=latom
!   sitb=latom
   !print*,'nat=',nat
   do i=1,nat
      atom=latom
      if (atom(4,i).ne.1.0) then
!         print*,'atom(1,i)=',atom(1,i)
         cycle
      end if
!      else
      closestr=200.0
!      print*,atom(:,i)
      Do j=1,nat
!         print*,'latom=',latom(1,j)
        do taux=-rep_cn(1),rep_cn(1)
        do tauy=-rep_cn(2),rep_cn(2)
        do tauz=-rep_cn(3),rep_cn(3)
          if(j.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
          dx=atom(1,j)-atom(1,i)+tau(1)
          dy=atom(2,j)-atom(2,i)+tau(2)
          dz=atom(3,j)-atom(3,i)+tau(3)
          r=(dx*dx+dy*dy+dz*dz)
          if (r.gt.crit_cn) cycle
          r=sqrt(r)
!         r=sqrt((atom(1,i)-latom(1,j))**2 + (atom(2,i)-latom(2,j))**2 + &
!         & (atom(3,i)-latom(3,j))**2)
!         print*,atom(1,i)
!         print*,'r=',r
!         print*,'closestr=',closestr
          if (r.eq.0.0) then
!            print*,atom(1,i)
            cycle
!         print*,'closestr=',closestr
!         print*,'r=',r
          elseif (r.lt.closestr) then
!            print*,'closestr=',closestr
            !print*,'r=',r
            closestr=r
            typ=latom(4,j)
            !print*,'type***************=',typ
          end if
!         print*,'closestr=',closestr
!         print*,'r=',r
        enddo !tauz
        enddo !tauy
        enddo !taux
      end do ! j
      !print*,'closestr=',closestr
      !print*,'type=',typ
      if (typ.eq.1) then
          !print*,"H connected to H"
          sitout(5,i)=typ
!          print*,'atom(5,i)*********************',sitaa(5,i)
      elseif (typ.eq.3) then
!     print "H connected to Li"
          sitout(5,i)=typ
      elseif (typ.eq.4) then
!     print "H connected to Be"
          sitout(5,i)=typ
      elseif (typ.eq.5) then
!     print "H connected to B"
          sitout(5,i)=typ
      elseif (typ.eq.6) then
!     print "H connected to C"
          sitout(5,i)=typ
!          sitaa(5,i)=int(sitaa(5,i))
!          print*,'sitaa(4,i)****',sitaa(4,i)
!          print*,'sitaa(5,i)**********************',sitaa(5,i)
      elseif (typ.eq.7) then
!     print "H connected to N"
          sitout(5,i)=typ
      elseif (typ.eq.8) then
          !print*,"H connected to O"
          sitout(5,i)=typ
!          print*,'sitaa(5,i)**********************',sitaa(5,i)
      elseif (typ.eq.9) then
!     print "H connected to F"
          sitout(5,i)=typ
      elseif (typ.eq.11) then
!     print "H connected to Na"
          sitout(5,i)=typ
      elseif (typ.eq.12) then
!     print "H connected to Mg"
          sitout(5,i)=typ
      elseif (typ.eq.13) then
!     print "H connected to Al"
          sitout(5,i)=typ
      elseif (typ.eq.14) then
!     print "H connected to Si"
          sitout(5,i)=typ
      elseif (typ.eq.15) then
!     print "H connected to P"
          sitout(5,i)=typ
      elseif (typ.eq.16) then
!     print "H connected to S"
          sitout(5,i)=typ
      elseif (typ.eq.17) then
!     print "H connected to Cl"
          sitout(5,i)=typ
!      elseif (typ.eq.) 1,3,4,5,6,7,8,9,11,12,13,14,15,16,17 then
!          print*,'Wrong neighbour! ',typ
      end if
   !print*,'analyze-hydrogen',sitout(5,i)
   end do ! i
   end subroutine
