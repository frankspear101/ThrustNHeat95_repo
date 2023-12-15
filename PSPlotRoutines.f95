      subroutine psopcl(itype)
      implicit none
!      routine to perform various postscript open/close operations
!      such as
!      1)toggle postscript output on and off
!      2)Open a postscript file and write header
!      3)Close a postscript file and write trailer
!      itype is an option mode variable
!      itype = 1 user interaction mode
!      itype = 2 shut down postscript file only
!      itype = 3 Automatic open of postscript file
!      itype = 4 Open Postscript scratch file
!      itype = 5 Save PostScript scratch file
      	character*132 readit
      	integer itype,status
! ------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
! ------------------------------------------------
!      	JPS is the postscript on/off toggle
!             0 = off
!             1 = on
!      
!      	psfile is the postscript file on/off flag
!            0 = file is not open
!            1 = file is open illustrator 88
! 	    2 = file is open illustrator 3.0 and 5.0
! 
!      	toggle PS output off if file is closed, to avoid any possible file errors
      	if(psfile.eq.0)jps=0 


! ***************************************
!        automatic open mode
        if(itype.eq.3)then
!        get new PS file name        
        psfile=0
      	OPEN(psunit,FILE='',STATUS='UNKNOWN',iostat=status)
      	if(status.ne.0)then
      		psfile=0
      		jps=0
      		return
      		endif
        INQUIRE (psunit, NAME=psfilename)
!        write header for Adobe illustrator
	call WritepsHeader
!          write(psunit,*)'%!PS-Adobe-3.0'
!           write(psunit,*)'%%Creator: Adobe Illustrator(TM) 3.0'
!           write(psunit,*)'%%For: (FSS) (RPI)'
!           write(psunit,*)'%%Title: (Adobe test 1)'
!           write(psunit,*)'%%CreationDate: (1/1/90) (12:00 PM)'
!           write(psunit,*)'%%BoundingBox: 310 415 359 470'
!           write(psunit,*)'%AI3_TileBox: 30 31 582 761'
!           write(psunit,*)'%%EndComments'
!           write(psunit,*)'%%EndProlog'
!           write(psunit,*)'/_Helvetica 14 12 0 1 z'
!          write(psunit,*)'/_Times-Roman 14 12 0 1 z'
!           write(psunit,*)'0 O'
!           write(psunit,*)'0 g'
!           write(psunit,*)'0 i 0 J 1 w 4 M []0 d'
          psfile=2      			! 2 means opened for illustrator 3.0 (works in Illustrator 5.0)
          jps=1
       return
        endif           			! end  automatic open mode


! ***************************************
!        Open Postscript Scratch file
        if(itype.eq.4)then
!        call volset(volrefnum)      		! open file on startup volume
!         OPEN(psunit,FILE='PostScript.tmp',STATUS='NEW')
        OPEN(psunit,FILE='PostScriptTemp.ai',STATUS='UNKNOWN')
!         OPEN(psunit,STATUS='SCRATCH')
!        write header for Adobe illustrator
	call writepsheader
          psfile=2            			! 2 means opened for illustrator 3.0 (works in Illustrator 5.0)
          jps=1

        endif           			! end  automatic open mode

! ***************************************
!       Save Postscript Scratch file
        if(itype.eq.5)then
!        get new PS file name        
!        iokL=.false.
!        iokL = stdfil(2,vref,'PostScript filename',psfilename,1,'TEXT')	! 2 is putfile and set vref
!        if (.not.iokL)then     						! no file name was specified
!            return
!            endif
!      	OPEN(21,FILE=psfilename,STATUS='UNKNOWN')
!      	OPEN(21,FILE='',STATUS='UNKNOWN',ACTION='WRITE',iostat=status)
      	OPEN(21,FILE='',STATUS='NEW',ACTION='WRITE',iostat=status)
      	if(status.ne.0)return
        INQUIRE (21, NAME=psfilename)
      	rewind(psunit)
      	do
      	read(psunit,'(A)',end=923)readit
      	write(21,'(A)')ADJUSTL(trim(readit))
      	repeat
923   	continue
      	backspace(psunit)
!       write(21,*)'%%Trailer'
      	write(21,924)
924	format('%%Trailer')
      	close(21)
      	endif           ! end  automatic save
        
        return
      end
      
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine writepsheader
      implicit none
! ------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
! ------------------------------------------------
	write(psunit,100)
100	format('%!PS-Adobe-3.0')
        write(psunit,101)
101	format('%%Creator: Adobe Illustrator(TM) 3.0')
        write(psunit,102)
102	format('%%For: (FSS) (RPI)')
        write(psunit,103)
103	format('%%Title: (Adobe test 1)')
        write(psunit,104)
104	format('%%CreationDate: (1/1/90) (12:00 PM)')
        write(psunit,105)
105	format('%%BoundingBox: 310 415 359 470')
        write(psunit,106)
106	format('%AI3_TileBox: 30 31 582 761')
        write(psunit,107)
107	format('%%EndComments')
        write(psunit,108)
108	format('%%EndProlog')
        write(psunit,109)
109	format('/_Helvetica 14 12 0 1 z')
!          write(psunit,*)'/_Times-Roman 14 12 0 1 z')
        write(psunit,110)
110	format('0 O')
        write(psunit,111)
111	format('0 g')
        write(psunit,112)
112	format('0 i 0 J 1 w 4 M []0 d')
	return
	end
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE puser(iauto)
	implicit none
!C	Routine to scale the user units for PS output
	integer iauto
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!	Include "FSS_Ppplot.inc"
!c------------------------------------------------
	real*4 fact
	
	SAVE

!c	do autoscaling if desired
	if(iauto.eq.1)then
	  if(psfile.eq.1) then     !Illustrator '88 scaling
	    pylen=15
	    pxlen=pylen*1.33333 
 	    pxor=70-0
 	    pyor=70+420
	  endif
	  if(psfile.eq.2) then    !Illustrator v 3.0 scaling  (also v 5)
	    pylen=15
	    pxlen=pylen*1.33333 
 	    pxor=100
 	    pyor=-200
	  endif
	else
	  if(psfile.eq.1) then          !Illustrator '88 scaling
	    pxlen=xlen*2
	    pylen=ylen*2
  	    pxor=70-0
 	    pyor = 70+420
	  endif
	  if(psfile.eq.2) then          !Illustrator v 3.0 scaling (Also v 5)
!c       I changed numbers here for the Illustrator v 5.0
	    pxlen=xlen
	    pylen=ylen
  	    pxor=100
 	    pyor = -200
	  endif
	endif
	fact=.035278
	pxconv=pxlen/(fact*(xmax-xmin))
	pyconv=pylen/(fact*(ymax-ymin))
	RETURN
	END
      
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine paxis()
	implicit none
!	integer nxstep,nystep,nxdec,nydec,i
	integer i
	real*4 xunits,yunits,xstep,ystep,x,y,xnum,ynum
	character*32 npout1
!	character*80 xlab,ylab,pltitle
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!	Include "FSS_Ppplot.inc"
!c------------------------------------------------

	SAVE
	xunits=xmax-xmin
	yunits=ymax-ymin
	write(psunit,'(A)')'1 w'	    !Weight of 1; can be changed if desired
!c
!c	-----------------------DRAW BOX---------------------
!c
	call pplot (xmin,ymin,0)
	call pplot (xmax,ymin,1)
	call pplot (xmax,ymax,1)
	call pplot (xmin,ymax,1)
	call pplot (xmin,ymin,1)
	write(psunit,*)'s'
!c
!c	--------------------draw tic marks-------------------
!c
	write(psunit,'(A)')'u'
	xstep=xunits/nxstep
	ystep=yunits/nystep
	do 10 i=1,nxstep-1
	  x=xmin+i*xstep
	  y=ymin
	  call pplot(x,y,0)
	  call pplot(x,y+.02*(ymax-ymin),1)
	  write(psunit,*)'S'
10	continue
	do 20 i=1,nystep-1
	  x=xmax
	  y=ymin+i*ystep
	  call pplot(x,y,0)
	  call pplot(x-.02*(xmax-xmin),y,1)
	  write(psunit,*)'S'
20	continue
	do 30 i=1,nxstep-1
	  x=xmin+i*xstep
	  y=ymax
	  call pplot(x,y,0)
	  call pplot(x,y-.02*(ymax-ymin),1)
	  write(psunit,*)'S'
30	continue
	do 40 i=1,nystep-1
	  x=xmin
	  y=ymin+i*ystep
	  call pplot(x,y,0)
	  call pplot(x+.02*(xmax-xmin),y,1)
	  write(psunit,*)'S'
40	continue
	write(psunit,'(A)')'U'
	write(psunit,*)'1 w'	 !Weight of 1; can be changed if desired
!c
!c	-----------------Put numbers on X and Y axes----------------
!c
2000	FORMAT(I13)
2001	FORMAT(F13.1)
2002	FORMAT(F13.2)
2003	FORMAT(F13.3)
2004	FORMAT(F13.4)
2005	FORMAT(F13.5)
2006	FORMAT(F13.6)
2110	format('[1 0 0 1 ',2f10.1,']e')
2111	format('[0 1 -1 0 ',2f10.1,']e')
2120	format(1i5,'(',a,')t')
!c
!c	------------------Put Numbers on X axis------------------
!c
	write(psunit,'(A)')'u'
	write(psunit,*)'/_Helvetica 18 18 0 1 z'
	do 50 i=1,nxstep+1,2
	  y=ymin
	  x=xmin+(i-1)*xunits/nxstep
	  xnum=x
	  if(nxdec.eq.0)write(npout1,2000)int(xnum)
	  if(nxdec.eq.1)write(npout1,2001)xnum
	  if(nxdec.eq.2)write(npout1,2002)xnum
	  if(nxdec.eq.3)write(npout1,2003)xnum
	  if(nxdec.eq.4)write(npout1,2004)xnum
	  if(nxdec.eq.5)write(npout1,2005)xnum
	  if(nxdec.ge.6)write(npout1,2006)xnum
!c
!c	trim leading and trailing blanks
!c
	  npout1=trim(adjustl(npout1))
!c
!c	5 pixels width per character in monaco 9 pt
!c	7 pixels height
!c	8 pixels width per character in 14 pt. Times
!c	about 12 pixels height
!c
!c	  write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-16
	  write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-20
!c	  write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-25
!c	  write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-32
!c	  write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-48
	  write(psunit,2120)len(trim(npout1)),trim(npout1)
	  write(psunit,*)'T'
50	CONTINUE
	write(psunit,'(A)')'U'
!c
!c	------------------PUT NUMBERS ON Y AXIS-----------------
!c
	write(psunit,'(A)')'u'
	write(psunit,*)'/_Helvetica 18 18 0 2 z'
	do 60 i=1,nystep+1,2
	  x=xmin
	  y=ymin+(i-1)*yunits/nystep
	  ynum=y
	  if(nydec.eq.0)write(npout1,2000)int(ynum)
	  if(nydec.eq.1)write(npout1,2001)ynum
	  if(nydec.eq.2)write(npout1,2002)ynum
	  if(nydec.eq.3)write(npout1,2003)ynum
	  if(nydec.eq.4)write(npout1,2004)ynum
	  if(nydec.eq.5)write(npout1,2005)ynum
	  if(nydec.ge.6)write(npout1,2006)ynum
!c
!c	Trim leading and trailing blanks
!c
	npout1=TRIM(ADJUSTL(npout1))
	  write(psunit,2110)pxor-8,-(pyor-(y-ymin)*pyconv)-8
!c	  write(psunit,2110)pxor-4,-(pyor-(y-ymin)*pyconv)-10
!c	  write(psunit,2110)pxor-4,-(pyor-(y-ymin)*pyconv)-12
!c	  write(psunit,2110)pxor-4,-(pyor-(y-ymin)*pyconv)-16
!c	  write(psunit,2110)pxor-4,-(pyor-(y-ymin)*pyconv)-32
	  write(psunit,2120)len(trim(npout1)),trim(npout1)
	  write(psunit,*)'T'
60	CONTINUE
	write(psunit,'(A)')'U'
	write(psunit,*)'/_Helvetica 26 26 0 1 z'
!c
!c	---------------------WRITE LABELS----------------------
!c
!c	---------------------WRITE X LABEL---------------------
!c
	x=xmin+.5*xunits
!c	write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-32
	write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-64
!c	write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor-100
       
	write(psunit,2120)len(trim(xlab)),trim(xlab)
	write(psunit,'(A1)')'T'
!c	---------------------WRITE Y LABEL---------------------
!c
	y=ymin+.5*yunits
!c	write(psunit,2111)pxor-32,-pyor+(y-ymin)*pyconv
	write(psunit,2111)pxor-64,-pyor+(y-ymin)*pyconv
!c	write(psunit,2111)pxor-100,-pyor+(y-ymin)*pyconv
	write(psunit,2120)len(trim(ylab)),trim(ylab)
	write(psunit,'(A1)')'T'
!c
!c	-------------------WRITE PLOT LABEL--------------------
	x=xmin+.5*xunits
	write(psunit,2110)pxor+(x-xmin)*pxconv,-pyor+(ymax-ymin)*pyconv+24
	write(psunit,2120)len(trim(pltitle)),trim(pltitle)
	write(psunit,'(A1)')'T'
	return
	end

!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE pplot (x,y,ipen)
	implicit none
!c	routine to move pen to x,y in user units
!c	ipen = 0 move with pen up
!c	ipen = 1 move with pen down
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!	Include "FSS_Ppplot.inc"
!c------------------------------------------------
	integer*4 ipen,i
	real*4 x,y
	save

	IF (ipen.EQ.0)then
		write(psunit,5)(CMYK(CurrentColor,i),i=1,4)
5		format(4f10.5,'    K')			! Capital "K" sets the line color (lower case "k" sets fill color)
		write(psunit,10)pxor+(x-xmin)*pxconv,-pyor+(y-ymin)*pyconv
10		format(2f10.1,' m')	
		else
		write(psunit,20)pxor+(x-xmin)*pxconv,-pyor+(y-ymin)*pyconv
20		format(2f10.1,' L')
		endif
	return
	END

!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	subroutine pslab(overstring,xplot,yplot,isize,ipos)
!c	subroutine to write out a text string to a Postscript file	
!c	overstring is the character string with text to write
!c	x and y are the position of string in user units
!c	isize is the point size of the text
!c	ipos	0 = left-justified
!c		1 = centered
!c		2 = right-justified


	implicit none
	character*(*) overstring
	real*4 x,y,xplot,yplot
	integer isize,ipos

!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!	Include "FSS_Ppplot.inc"
!c------------------------------------------------
	
10	format('[1 0 0 1',2f10.1,']e')
20	format(1i5,'(',A,')t')
30	format('/_Helvetica ',4i3,' z')
	if(jps.eq.1) then
	  x=pxor+(xplot-xmin)*pxconv
	  y=-pyor+(yplot-ymin)*pyconv
	  write(psunit,30)isize,isize,0,ipos
	  write(psunit,10)x,y
!c	set color of text
!c        1      2       3         4     5        6      7         8  
!c      black  yellow  magenta   red   cyan   green     blue  white
!c	select case (icolor)
!c	case(1)		!black
!c	write(psunit,*)'     0  0  0  1   K'
!c	case(2)		!Yellow
!c	write(psunit,*)'     0  0  1  0   K'
!c	case(3)		!magenta
!c	write(psunit,*)'     0  1  0  0   K'
!c	case(4)		!red
!c	write(psunit,*)'     0  1  1  0   K'
!c	case(5)		!cyan
!c	write(psunit,*)'     1  0  0  0   K'
!c	case(6)		!green
!c	write(psunit,*)'     1  0  1  0   K'
!c	case(7)		!blue
!c	write(psunit,*)'     1  .5  0  0   K'
!c	case(8)		!white
!c	write(psunit,*)'     0  0  0  0   K'
!c	case default
!c	end select

	  write(psunit,20)len(trim(overstring)),trim(overstring)
	  write(psunit,'(A1)')'T'
	endif
	return
	end
	
	
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
	Subroutine ppenup
	implicit none
!c	routine to write an end of path command to postscript file
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!c------------------------------------------------
	if(jps.eq.1)write(psunit,*)'S'
	return
	end
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
 	SUBROUTINE psCircle(xp,yp,radp,SFB)
!	Routine to draw a circle centered at centerx,centery with radius=radius
!		Units are user units
!		ColorS is color stroke CMYK
!		ColorF is color fill   CMYK
!		SFB 1=stroke = s; 2=fill = f; 3=both = b
!		radius is the radius in the X user direction
!			e.g. if cx = T = 500 and rad = 3 then the radius is 3 degrees
	implicit none
	INCLUDE "PlotStuff.inc"	
	integer*4 SFB,i
	real*4 xp,yp,radp,cx,cy,rad,x,y,t1x,t1y,t2x,t2y,tangk
	tangk = 0.55222966

!	write(psunit,*)'     0  0  0  1   K'		!black stroke
!	write(psunit,*)'     .611765 0 1 0 k'		! green fill
!	write(psunit,5)(colorS(i),i=1,4)
!	write(psunit,6)(colorF(i),i=1,4)

	write(psunit,5)(CMYK(CurrentColor,i),i=1,4)
	write(psunit,6)(CMYK(CurrentColor,i),i=1,4)
5	format(4f10.5,'    K')			! Capital "K" sets the line color (lower case "k" sets fill color)
6	format(4f10.5,'    k')

!												=centerx+radius	=centery	m
!=centerx+radius	=centery+radius*tangk	=centerx+radius*tangk	=centery+radius		=centerx	=centery+radius	c
!=centerx-radius*tangk	=centery+radius		=centerx-radius		=centery+radius*tangk	=centerx-radius	=centery	c
!=centerx-radius	=centery-radius*tangk	=centerx-radius*tangk	=centery-radius		=centerx	=centery-radius	c
!=centerx+radius*tangk	=centery-radius		=centerx+radius		=centery-radius*tangk	=centerx+radius	=centery	c

10	format(2f15.4,' m')	
20	format(6f15.4,' c')
! For things to be scaled correctly, cx, cy and radius need to first be converted to psuser units
!	write(*,*)pxconv,pyconv
!	write(*,*)xp,yp,radp
	cx =  pxor+(xp-xmin)*pxconv
	cy = -pyor+(yp-ymin)*pyconv
	rad = radp*pxconv
!	write(*,*)cx,cy,rad
! =centerx+radius	=centery	m
	x = cx + rad
	y = cy
!	x =  pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,10)x,y
!=centerx+radius	=centery+radius*tangk	=centerx+radius*tangk	=centery+radius		=centerx	=centery+radius	c
	t1x = cx + rad
	t1y = cy + rad*tangk
	t2x = cx + rad*tangk
	t2y = cy + rad
	x = cx
	y = cy + rad
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x =  pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y
!=centerx-radius*tangk	=centery+radius		=centerx-radius		=centery+radius*tangk	=centerx-radius	=centery	c
	t1x = cx - rad*tangk
	t1y = cy + rad
	t2x = cx - rad
	t2y = cy + rad*tangk
	x = cx - rad
	y = cy
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x = pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y
!=centerx-radius	=centery-radius*tangk	=centerx-radius*tangk	=centery-radius		=centerx	=centery-radius	c
	t1x = cx - rad
	t1y = cy - rad*tangk
	t2x = cx - rad*tangk
	t2y = cy - rad
	x = cx
	y = cy - rad
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x = pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y
!=centerx+radius*tangk	=centery-radius		=centerx+radius		=centery-radius*tangk	=centerx+radius	=centery	c
	t1x = cx + rad*tangk
	t1y = cy - rad
	t2x = cx + rad
	t2y = cy - rad*tangk
	x = cx + rad
	y = cy
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x = pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y

	select case(SFB)
	case(1)		! stroke only
	write(psunit,*)'s'
	case(2)		! fill only
	write(psunit,*)'f'
	case(3)		! stroke and fill
	write(psunit,*)'b'
	case default
	end select

	return
	end	
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
 	SUBROUTINE psymb(xus,yus,isymb,size)
	implicit none
!c	Routine to plot a centered symbol at x,y
!c	isymb = symbol type
!c	Size = half size of symbol in pixels (72 dpi)
!c	symbols are:
!c		     0 = blank
!c		     1 = vertical line
!c		     2 = horizontal line
!c		     3 = +
!c		     4 = cross x
!c		     5 = * (star)
!c		     6 = triangle
!c		     7 = upside down triangle
!c		     8 = square
!c		     9 = diamond
	integer*4 isymb
	integer*4 size
	real*4 x,y,xus,yus
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!	Include "FSS_PPPlot.inc"
!c------------------------------------------------

1	format('u')
10	format(2f10.1,' m')	
20	format(2f10.1,' L')	
30	format('S')
99	format('U')
	x = pxor+(xus-xmin)*pxconv
	y = -pyor+(yus-ymin)*pyconv
	if(isymb.eq.0)return	!no symbol drawn
	if(isymb.eq.1)then	!vertical line
	write(psunit,1)
	write(psunit,10)x,y+size
	write(psunit,20)x,y-size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.2)then	!horizontal line
	write(psunit,1)
	write(psunit,10)x+size,y
	write(psunit,20)x-size,y
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.3)then	!plus +
	write(psunit,1)
	write(psunit,10)x+size,y
	write(psunit,20)x-size,y
	write(psunit,30)
	write(psunit,10)x,y+size
	write(psunit,20)x,y-size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.4)then	!cross (X)
	write(psunit,1)
	write(psunit,10)x+size,y+size
	write(psunit,20)x-size,y-size
	write(psunit,30)
	write(psunit,10)x+size,y-size
	write(psunit,20)x-size,y+size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.5)then	!star (*)
	write(psunit,1)
	write(psunit,10)x+size,y
	write(psunit,20)x-size,y
	write(psunit,30)
	write(psunit,10)x,y+size
	write(psunit,20)x,y-size
	write(psunit,30)
	write(psunit,10)x+size,y+size
	write(psunit,20)x-size,y-size
	write(psunit,30)
	write(psunit,10)x+size,y-size
	write(psunit,20)x-size,y+size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.6)then	!triangle
	write(psunit,1)
	write(psunit,10)x,y
	write(psunit,20)x,y-size
	write(psunit,20)x-size,y-size
	write(psunit,20)x,y+size
	write(psunit,20)x+size,y-size
	write(psunit,20)x,y-size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.7)then	!upside down triangle
	write(psunit,1)
	write(psunit,10)x,y
	write(psunit,20)x,y+size
	write(psunit,20)x-size,y+size
	write(psunit,20)x,y-size
	write(psunit,20)x+size,y+size
	write(psunit,20)x,y+size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.8)then	!square
	write(psunit,1)
	write(psunit,10)x,y
	write(psunit,20)x,y+size
	write(psunit,20)x-size,y+size
	write(psunit,20)x-size,y-size
	write(psunit,20)x+size,y-size
	write(psunit,20)x+size,y+size
	write(psunit,20)x,y+size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif
	if(isymb.eq.9)then	!diamond
	write(psunit,1)
	write(psunit,10)x,y
	write(psunit,20)x,y+size
	write(psunit,20)x-size,y
	write(psunit,20)x,y-size
	write(psunit,20)x+size,y
	write(psunit,20)x,y+size
	write(psunit,30)
	write(psunit,99)
	go to 999
	endif

999	continue
	return
	END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotBoxInPSFile(T,dT,P,dP)
	real*4 T,dT,P,dP
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"

		write(psunit,*)'u'	! this should "group" all of the same assemblage for later coloring
		call pplot(T-dT,(P+dP)/1000,0)
		call pplot(T+dT,(P+dP)/1000,1)
		call pplot(T+dT,(P-dP)/1000,1)
		call pplot(T-dT,(P-dP)/1000,1)
		call pplot(T-dT,(P+dP)/1000,1)
		call ppenup				! writes out an "S"
		write(psunit,*)'U'		! Turn off "group"

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine FillBoxInPSFile(T,dT,P,dP)
	real*4 T,dT,P,dP
	integer*4 i
	INCLUDE "PlotStuff.inc"	

!	write(psunit,5)(color(i),i=1,4)
	write(psunit,5)(CMYK(CurrentColor,i),i=1,4)
5	format(4f10.5,'    k')				! lower case "k" sets fill color

	write(psunit,*)'u'	! this should "group" all of the same assemblage for later coloring
	call pplot(T-dT,(P+dP)/1000,0)		! note that these calls will set a stroke color, but we're not stroking
	call pplot(T+dT,(P+dP)/1000,1)
	call pplot(T+dT,(P-dP)/1000,1)
	call pplot(T-dT,(P-dP)/1000,1)
	call pplot(T-dT,(P+dP)/1000,1)
	write(psunit,*)'f'		! fill box
!		call ppenup				! writes out an "S"
	write(psunit,*)'U'		! Turn off "group"

	return
	end
	

