! 	F77 -c -f Graphics.f
!	SUBROUTINE axis(xtic,nxdec,xlab,ytic,nydec,ylab,pltitle)
!	SUBROUTINE user(xxor,xxmin,xxmax,xxlen,yyor,yymin,yymax,yylen)
!      	SUBROUTINE al2sio5(itype,isw)
!      	SUBROUTINE abqtz (isw)
!      	subroutine setplt(Xor,Xmin,Xmax,Xlen,Xtic,Xdec,Xlab,Yor,Ymin,Ymax,Ylen,Ytic,Ydec,Ylab,Pltitle,iok)
!      	subroutine settriangle(Xor,Yor,Llab,Rlab,Tlab,scale,itic,Pltitle,iok)
!      	subroutine triaxis(Llab,Rlab,Tlab,scale,itic,Pltitle)
!      	subroutine tripaxis(Llab,Rlab,Tlab,scale,itic,Pltitle)
!	subroutine triplot(left,right,top,scale,ipen)
!	subroutine tripplot(left,right,top,scale,ipen)
!	subroutine trisymb(left,right,top,scale,isymb,isize)
!	subroutine tripsymb(left,right,top,scale,isymb,isize)
!	subroutine listsy
! 	SUBROUTINE symb(xus,yus,isymb,size)
!	SUBROUTINE symbpl (x,y,itype,xsize,ysize,scaleit,wnumT)
!	subroutine psopcl(itype)
!      	subroutine WritePSHeader
!	SUBROUTINE puser(iauto)
!	subroutine paxis(xtic,nxdec,xlab,ytic,nydec,ylab,pltitle)
!	SUBROUTINE pplot (x,y,ipen)
!	subroutine pslab(overstring,xplot,yplot,isize,ipos)
!	Subroutine ppenup
! 	SUBROUTINE psymb(xus,yus,isymb,size)
!---------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    	Subroutine Setplt(xor,xmin,xmax,xlen,nxstep,nxdec,xlab,
!     &			yor,ymin,ymax,ylen,nystep,nydec,ylab,
!     &			pltitle,iok)
    	Subroutine Setplot(iok)
    	use AWE_Interfaces
   	implicit none
	INCLUDE "PlotStuff.inc"	
!	INCLUDE "FSS_PPPlot.inc"	
    	integer :: result, iok
!    	character*2 num
    	character *16 editChars
!    	character *16 mytext
!
! dialog and item types
!    
    	type(AWE_FormDialog)     :: Plotdialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! stastep text
!    	type(AWE_FormComboBox)   :: comboBox        ! drop down menu
	type(AWE_FormLineEdit)   :: lineEditPlTitle       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditXmin       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditXmax       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditXstep       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditXdec       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditXlab       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYmin       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYmax       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYstep       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYdec       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYlab       ! text box with default text

! set the dialog title and create it
    Plotdialog%title = "Plot values"
    call AWE_createDialog(Plotdialog)

! add stastuc text

!    dialoglabel%text = "Plotter Variables"
!    call AWE_addToDialog(dialoglabel, Plotdialog)

	lineEditPlTitle%title = "Plot Title"
	lineEditPlTitle%text = Pltitle
	call AWE_addToDialog(lineEditPlTitle, Plotdialog)

	lineEditXmin%title = "Xmin"
	write(editChars,*)Xmin
	lineEditXmin%text = editChars
	call AWE_addToDialog(lineEditXmin, Plotdialog)
	lineEditXmax%title = "Xmax"
	write(editChars,*)Xmax
	lineEditXmax%text = editChars
	call AWE_addToDialog(lineEditXmax, Plotdialog)
	lineEditXstep%title = "X step size"
	write(editChars,*)nxstep
	lineEditXstep%text = editChars
	call AWE_addToDialog(lineEditXstep, Plotdialog)
	lineEditXdec%title = "X decimals"
	write(editChars,*)nXdec
	lineEditXdec%text = editChars
	call AWE_addToDialog(lineEditXdec, Plotdialog)
	lineEditXlab%title = "X label"
	lineEditXlab%text = Xlab
	call AWE_addToDialog(lineEditXlab, Plotdialog)

	lineEditYmin%title = "Ymin"
	write(editChars,*)Ymin
	lineEditYmin%text = editChars
	call AWE_addToDialog(lineEditYmin, Plotdialog)
	lineEditYmax%title = "Ymax"
	write(editChars,*)Ymax
	lineEditYmax%text = editChars
	call AWE_addToDialog(lineEditYmax, Plotdialog)
	lineEditYstep%title = "Y step size"
	write(editChars,*)nystep
	lineEditYstep%text = editChars
	call AWE_addToDialog(lineEditYstep, Plotdialog)
	lineEditYdec%title = "Y decimals"
	write(editChars,*)nYdec
	lineEditYdec%text = editChars
	call AWE_addToDialog(lineEditYdec, Plotdialog)
	lineEditYlab%title = "Y label"
	lineEditYlab%text = Ylab
	call AWE_addToDialog(lineEditYlab, Plotdialog)





! show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(Plotdialog)




    	if (result == 1) then
		iok = 0
 		read(lineEditPltitle%text,*)Pltitle
      		print *,"Plot title =  ", Pltitle
 		read(lineEditXmin%text,*)Xmin
      		print *,"Xmin =      ", Xmin
 		read(lineEditXmax%text,*)Xmax
      		print *,"Xmax =      ", Xmax
 		read(lineEditXstep%text,*)nxstep
      		print *,"Xstep =      ", nxstep
 		read(lineEditXdec%text,*)nxdec
      		print *,"Xdec =      ", nXdec
 		read(lineEditXlab%text,*)Xlab
      		print *,"X label =  ", Xlab

 		read(lineEditYmin%text,*)Ymin
      		print *,"Ymin =      ", Ymin
 		read(lineEditYmax%text,*)Ymax
      		print *,"Ymax =      ", Ymax
 		read(lineEditYstep%text,*)nystep
      		print *,"Ystep =      ", nystep
 		read(lineEditYdec%text,*)nYdec
      		print *,"Ydec =      ", nYdec
 		read(lineEditYlab%text,*)Ylab
      		print *,"Y label =  ", Ylab
		else
		iok = 1
   		endif
	
	end
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotAxes(Canvas)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: Canvas

!	include "MMAD_Plotter_2.inc"
	INCLUDE "PlotStuff.inc"	
!	include "Plotted.inc"
!	Include "FSS_PostScript.inc"


      	close(PSUNIT)         			! close old PostScript scratch file
      	call psopcl(4)    			! open a new scratch file

      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot

!	xlen and ylen are in cm. 
!	the canvas width and height are in pixels
!	to make sure the canvas is large enough we need to convert
!	The scaling in subroutine USER assumes 72 dpi so
!	cm = pixels * 2.54/72 or
!	pixels = cm * 72/2.54 = cm * 28.34646
!	add 25% just to be safe
	canvas%width = 1.25*(xlen * 28.34646)			! *2 is needed to make room for the assemblage information
	canvas%height = 1.5*(ylen * 28.34646)
!	canvas%backgroundColor = AWE_teal
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = canvas%height - 100.
	CALL AWE_createCanvas(canvas)	

      	call axis(canvas)			! draws the axis in AWE
      	call paxis()			! draws the axis in postscript file

	return
	end
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE axis(canvas)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	TYPE(AWE_CanvasPen) :: pen
	TYPE(AWE_Rect) :: rect
!	CHARACTER(LEN=*) text
!	INTEGER, optional :: flags
	INTEGER(KIND=4) :: flags
	TYPE(AWE_Font) :: font
	INTEGER :: textColor = AWE_black
!	INTEGER, OPTIONAL:: textColor = AWE_black
	integer*4 ticsize,i,nxmove,nychar
	real*4 xstep,ystep,xunits,yunits,x,y,xnum,ynum
	CHARACTER*32 npout1
	character*1 char1
!	character*80 xlab,ylab,pltitle
	real xpix,ypix
!c------------------------------------------------
!	Include "FSS_Ppplot.inc"
	INCLUDE "PlotStuff.inc"	
!c------------------------------------------------
	SAVE
	xunits=xmax-xmin
	yunits=ymax-ymin
	call plot(canvas,xmin,ymin,0)
	call plot(canvas,xmax,ymin,1)
	call plot(canvas,xmax,ymax,1)
	call plot(canvas,xmin,ymax,1)
	call plot(canvas,xmin,ymin,1)
!c	 draw tic marks
	 ticsize=4		!this should be some % of window size
	 xstep=xunits/nxstep
	 ystep=yunits/nystep
	 Do 10 i=1,nxstep-1
	 x=xmin+i*xstep
	 y=ymin
	 call plot(canvas,x,y,0)
!c	CALL LINE (0,-ticsize)
	  call plot(canvas,x,y+.02*(ymax-ymin),1)
10	continue
	do 20 i=1,nystep-1
	 x=xmax
	 y=ymin+i*ystep
	 call plot(canvas,x,y,0)
!c	CALL LINE (-ticsize,0)
	  call plot(canvas,x-.02*(xmax-xmin),y,1)
20	continue
	Do 30 i=1,nxstep-1
	 x=xmin+i*xstep
	 y=ymax
	 call plot(canvas,x,y,0)
!c	 CALL LINE (0,ticsize)
	 call plot(canvas,x,y-.02*(ymax-ymin),1)
30	continue
	do 40 i=1,nystep-1
	  x=xmin
	  y=ymin+i*ystep
	  call plot(canvas,x,y,0)
!c	CALL LINE (ticsize,0)
	  call plot(canvas,x+.02*(xmax-xmin),y,1)
40	continue
!c	 Put numbers on X and Y axes
	font%pointSize = 14
!	flags = AWE_TextFlag_AlignHCenter		! this makes text disappear (I suspect it's outside of the rect)
	do 50 i=1,nxstep+1,2
	 y=ymin
	 x=xmin+(i-1)*xunits/nxstep
	 xnum=x
!	call plot(canvas,x,y,0)
	 npout1=' '
	IF(nxdec.EQ.0)WRITE(npout1,2000)INT(xnum)
2000	FORMAT(I13)
	IF(nxdec.EQ.1)WRITE(npout1,2001)xnum
2001	FORMAT(F13.1)
	IF(nxdec.EQ.2)WRITE(npout1,2002)xnum
2002	FORMAT(F13.2)
	IF(nxdec.EQ.3)WRITE(npout1,2003)xnum
2003	FORMAT(F13.3)
	IF(nxdec.EQ.4)WRITE(npout1,2004)xnum
2004	FORMAT(F13.4)
	IF(nxdec.EQ.5)WRITE(npout1,2005)xnum
2005	FORMAT(F13.5)
	IF(nxdec.GE.6)WRITE(npout1,2006)xnum
2006	FORMAT(F13.6)
!c	Trim leading and trailing blanks
	npout1=ADJUSTL(npout1)
	nxmove=5*(LEN(TRIM(npout1))/2.)
!c	5 pixels width per character in monaco 9 pt
!c	7 pixels height
!c	8 pixels width per character in 14 pt. Times
!c	about 12 pixels height
	rect%origin%x = xpix(x) - 10
	rect%origin%y = ypix(y) + 5		! rectangle is 10 pixels below point
!	rect%size%width = nxmove
	rect%size%width = 40
	rect%size%height = 15
	CALL AWE_canvasDrawText(canvas, rect, npout1, flags, font, textColor)
!	call fss_move(-nxmove,10)
!	call fss_drawstring(npout1)
!c	WRITE(*,'(A)')TRIM(npout1)
50	CONTINUE

!c	 PUT NUMBERS ON Y AXIS
	 Do 60 i=1,nystep+1,2
	 x=xmin
	 y=ymin+(i-1)*yunits/nystep
	 ynum=y
!	 call plot(canvas,x,y,0)
	 npout1=' '
	IF(nydec.EQ.0)WRITE(npout1,2000)int(ynum)
	IF(nydec.EQ.1)WRITE(npout1,2001)ynum
	IF(nydec.EQ.2)WRITE(npout1,2002)ynum
	IF(nydec.EQ.3)WRITE(npout1,2003)ynum
	IF(nydec.EQ.4)WRITE(npout1,2004)ynum
	IF(nydec.EQ.5)WRITE(npout1,2005)ynum
	IF(nydec.GE.6)WRITE(npout1,2006)ynum
!c	Trim leading blanks
	npout1=ADJUSTL(npout1)
	nxmove=5*(LEN(TRIM(npout1)))
!c	5 pixels width per character in monaco 9 pt
!c	7 pixels height
	rect%origin%x = xpix(xmin) - nxmove - 8
	rect%origin%y = ypix(y) - 8		! move the rect up 8 pixels from the value
	rect%size%width = nxmove + 5
	rect%size%height = 15
!	flags = AWE_TextFlag_AlignRight
	CALL AWE_canvasDrawText(canvas, rect, npout1, flags, font, textColor)
!	call fss_move(-nxmove,3)
!	call fss_drawstring(npout1)
!c	WRITE(*,'(A)')TRIM(npout1)
60	CONTINUE

!c	WRITE LABELS
!c	WRITE X LABEL
!c	FIRST LOCATE THE MIDPOINT OF X AXIS
	 x=xmin+.5*xunits
!	call plot(canvas,x,ymin,0)
!c	shift left by 1/2 of label length
!c	and down by 3 character heights (21 pixels)
	nxmove=12*(LEN(TRIM(xlab))/2.)  !12 pixels per character
	rect%origin%x = xpix(x) - nxmove		! shift the origin to the left
	rect%origin%y = ypix(ymin) + 25 
	rect%size%width = nxmove*2 + 20
	rect%size%height = 28
	font%pointSize = 28
	CALL AWE_canvasDrawText(canvas, rect, xlab, flags, font, textColor)
!	CALL fss_move (-nxmove,21)
!	call fss_drawstring(xlab)
!c	write(*,'(A)')TRIM(xlab)


!c	print ylabel
!c	here, npout1 still contains the largest ynum value
	nxmove=5*(LEN(TRIM(npout1))) + 25
	nychar=len(trim(ylab))
	 x=xmin
	 y=ymin+yunits/2.
	do 110 i=1,nychar
!	 call plot(canvas,x,y,0)
	rect%origin%x = xpix(xmin) - 50
	rect%origin%y = ypix(y) + i*28 - ypix(y)/3
	rect%size%width = 25
	rect%size%height = 28
!	flags = AWE_TextFlag_TextWordWrap
	char1 = ylab(i:i)
	CALL AWE_canvasDrawText(canvas, rect, char1, flags, font, textColor)
!	CALL AWE_canvasDrawText(canvas, rect, text, flags, font, textColor)

!	 call fss_move(-nxmove,-nychar*5+i*10)	!move 10 pixels per character
!	call fss_drawstring(ylab(i:i))
!c	write(*,'(A)')ylab(i:i)
110	continue
!c
!c	print Plot label
!c	first locate midpoint of Plot	
	 x=xmin + .5*xunits
!	call plot(canvas,x,ymax,0)
	nxmove=15*(LEN(TRIM(pltitle))/2.)   !15 pixels per character
	rect%origin%x = xpix(x) -  nxmove
	rect%origin%y = ypix(ymax) - 30
	rect%size%width = nxmove*2 + 20
	rect%size%height = 28
	CALL AWE_canvasDrawText(canvas, rect, pltitle, flags, font, textColor)
!	 CALL fss_move(-nxmove,-15)
!	call fss_drawstring(pltitle)
!c	write(*,'(A)')pltitle
	return
	end
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE user()
!	SUBROUTINE user(xxor,xxmin,xxmax,xxlen,yyor,yymin,yymax,yylen)
!C	ROUTINE TO SCALE THE USER UNITS
	implicit none
!	real*4 xxor,xxmin,xxmax,xxlen,yyor,yymin,yymax,yylen,xunits,yunits,fact
	real*4 xunits,yunits,fact
	integer iauto
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!	Include "FSS_PostScript.inc"
!	Include "FSS_Ppplot.inc"
!c------------------------------------------------

	SAVE

	fact=.035278     !cm/pixel (at 72 dpi)

!	xor=xxor
!	yor = yyor + (yylen / fact)
!	ylen=yylen
!	xlen=xxlen

!	xmin=xxmin
!	xmax=xxmax
!	ymin=yymin
!	ymax=yymax
	xunits=xmax-xmin
	yunits=ymax-ymin
	xconv=xlen/(fact*xunits)
	yconv=ylen/(fact*yunits)
	write(*,*)xconv,yconv
	iauto = 0
 	call puser(iauto)
	RETURN
	END
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c*******************************************************
!c*******************************************************
	Real function XPIX(x)
	INCLUDE "PlotStuff.inc"	
!	Include "fss_ppplot.inc"
	real*4 x
!	integer*4 xpix
!	xpix = int(xor+(x-xmin)*xconv + .5)
	xpix = xor + (x-xmin)*xconv
	return
	end
!c*******************************************************
!c*******************************************************
	real function YPIX(y)
	INCLUDE "PlotStuff.inc"	
!	Include "fss_ppplot.inc"
	real*4 y
!	integer*4 ypix
!	ypix = int(yor-(y-ymin)*yconv + .5)
	ypix = yor - (y-ymin)*yconv
	return
	end
!c*******************************************************
!c*******************************************************
	SUBROUTINE plot (canvas,x,y,ipen)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	TYPE(AWE_Line) :: line(1)
	TYPE(AWE_CanvasPen) :: pen

!c	routine to move pen to x,y in user units
!c	ipen = 0 move with pen up
!c	ipen = 1 move with pen down
!	integer*2 ix2,iy2,ixold2,iyold2
	integer*4 ix,iy,ipen,ixold,iyold
	real*4 x,y
	INCLUDE "PlotStuff.inc"	
!	Include "fss_ppplot.inc"
!c***********************************************************
!c      real*4 xor,xconv,xmin,xmax,xlen,yor,yconv,ymin,ymax,ylen
!c      real*4 pxor,pxconv,pxlen,pyor,pyconv,pylen
!c	COMMON /PPPlot/
!c     +	xor,xconv,xmin,xmax,xlen,
!c     +	yor,yconv,ymin,ymax,ylen,
!c     +  pxor,pyor,pxconv,pyconv,pxlen,pylen
!c***********************************************************
	save
	ix=int(xor+(x-xmin)*xconv + .5)
	iy=int(yor-(y-ymin)*yconv + .5)
!	ix2 = ix		! we don't need to use integer*2
!	iy2 = iy
!	pen%penColor = 16711680		! red line for testing (it works)
	pen%penColor = AWEcolorNumber(CurrentColor)		! Color for AWE line - in PlotStuff.inc common block

	IF (ipen.EQ.0) then
!		CALL MoveTo (val2(ix2),val2(iy2))
!		If this is a move, then the start and end of the line are the same
		line%start%x = ix
		line%start%y = iy
		line%end%x = ix
		line%end%y = iy
		CALL AWE_canvasDrawLines(canvas, line, pen)		
		else
!		endif
!	IF (ipen.EQ.1) THEN 
!		CALL MoveTo (val2(ixold2),val2(iyold2))
!		CALL LineTo (val2(ix2),val2(iy2))
!		If this is a line, then the start is the last point called
		line%start%x = ixold
		line%start%y = iyold
		line%end%x = ix
		line%end%y = iy
		CALL AWE_canvasDrawLines(canvas, line, pen)		
		ENDIF
	ixold=ix
	iyold=iy
!c	xold=xor+(x-xmin)*xconv
!c	yold=-yor+(y-ymin)*yconv
!	ixold2 = ixold
!	iyold2 = iyold
	return
	END
!c*******************************************************
!c*******************************************************

	SUBROUTINE TextOnPlot(canvas,xplot,yplot,text,textSize)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	TYPE(AWE_CanvasPen) :: pen
	TYPE(AWE_Rect) :: rect
!	INTEGER, optional :: flags
	INTEGER :: flags
	TYPE(AWE_Font) :: font
	INTEGER :: textColor = AWE_black
	integer textSize
	real xplot,yplot,width
	real xpix,ypix                ! functions
	CHARACTER(LEN=*) text
	
	font%pointSize = textSize
	rect%origin%x = xpix(xplot)			! functions to return pixel values of x and y given plot dimensions
	rect%origin%y = ypix(yplot) 		! rectangle is 10 pixels below point
	text = adjustL(text)
	width = 1.4*float(textSize*(LEN(TRIM(text))))       !12 pixels per character

	rect%size%width = width
	rect%size%height = textSize
	flags = AWE_TextFlag_AlignLeft
	CALL AWE_canvasDrawText(canvas, rect, text, flags, font, textColor)
	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE abqtz (canvas,isw)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
!
!     subroutine to plot the alpha-beta quartz transition following Berman's database
!     isw is a switch
!	 isw = 0  y axis is in pressure
!	 isw = 1  y axis is in Km
!
!	T is in C
!	P is in kbar
!
!------------------------------------------------
	INCLUDE "PlotStuff.inc"	
!------------------------------------------------
	real*4 Pinc,Tstart,Pstart,Pb,Pkbar,T
	integer*4 isw,iup,i


!	line for a-b quartz transition from Gibbs subroutine ASV
!	TBETA	= 0.024 * (Pb-5000.d0) + 959.d0    (T in K)

	Pinc=100.
!	Tstart=565.85
	Tstart=573
	Pstart=0
	T=Tstart
	Pb=0
	Pkbar=0
	iup=0
	do 10 i=1,200
    	Pb=Pb+Pinc
!	T = 0.024 * (Pb-5000.d0) + 685.85
	T = 0.024 * Pb + 573						! The transition should be at 573 at 1 bar.
      	Pkbar=Pb/1000.
	if(T.LT.xmax.and.T.Gt.xmin.and.Pkbar.gt.ymin.and.Pkbar.lt.ymax)then      
      		IF(isw.EQ.1)Pkbar=Pkbar/.275
		call plot(canvas,T,Pkbar,iup)
      		if(jps.eq.1)call pplot(T,Pkbar,iup)
      		iup=1
		endif
10    CONTINUE
	if(jps.eq.1)call ppenup
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotBoxOnScreen(canvas,T,dT,P,dP)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	real*4 T,dT,P,dP

!	Include "FSS_PostScript.inc"
	!	write(*,*)T,P
		call plot(canvas,T-dT,(P+dP)/1000.,0)
		call plot(canvas,T+dT,(P+dP)/1000.,1)
		call plot(canvas,T+dT,(P-dP)/1000.,1)
		call plot(canvas,T-dT,(P-dP)/1000.,1)
		call plot(canvas,T-dT,(P+dP)/1000.,1)
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotRectOnScreen(canvas,T,dT,P,dP,brush,pen)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	TYPE(AWE_Rect) :: rect

	real xpix,ypix
	real*4 T,dT,P,dP
	!	write(*,*)T,dT,P,dP
	rect%origin%x = xpix(T)
	rect%origin%y = ypix(P)	
!	rect%size%width = xpix(dT)
!	rect%size%height = ypix(dP)
	rect%size%width = xpix(T+dT) - xpix(T)
	rect%size%height = ypix(P+dP) - ypix(P)
	!	write(*,*)rect%origin%x,rect%size%width,rect%origin%y,rect%size%height
	call AWE_canvasDrawRect(canvas,rect,pen,brush)	
	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotCenteredRectOnScreen(canvas,T,dT,P,dP,brush,pen)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	TYPE(AWE_Rect) :: rect

	real xpix,ypix
	real*4 T,dT,P,dP,halfwidth,halfheight,centerx,centery
!	write(*,*)T,dT,P,dP
!	brush%brushColor = AWE_gold

	centerx = xpix(T)
	centery = ypix(P)
	halfwidth = xpix(T + dT) - xpix(T)
	halfheight = ypix(P - dP) - ypix(P)
!	write(*,*)centerx,halfwidth,centery,halfheight
	rect%origin%x = centerx - halfwidth
	rect%origin%y = centery - halfheight	
	rect%size%width = 2.*halfwidth
	rect%size%height = 2.*halfheight
!	write(*,*)rect%origin%x,rect%size%width,rect%origin%y,rect%size%height
	call AWE_canvasDrawRect(canvas,rect,pen,brush)	

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotCenteredEllipseOnScreen(canvas,T,dT,P,dP,brush,pen)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	TYPE(AWE_Rect) :: rect

	real xpix,ypix
	real*4 T,dT,P,dP,halfwidth,halfheight,centerx,centery
!	write(*,*)T,dT,P,dP
!	brush%brushColor = AWE_gold

	centerx = xpix(T)
	centery = ypix(P)
	halfwidth = xpix(T + dT) - xpix(T)
	halfheight = ypix(P - dP) - ypix(P)
!	write(*,*)centerx,halfwidth,centery,halfheight
	rect%origin%x = centerx - halfwidth
	rect%origin%y = centery - halfheight	
	rect%size%width = 2.*halfwidth
	rect%size%height = 2.*halfheight
!	write(*,*)rect%origin%x,rect%size%width,rect%origin%y,rect%size%height
	call AWE_canvasDrawEllipse(canvas,rect,pen,brush)	

	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine ReadColorFile()
	implicit none
	integer i,j,R,G,B
	character*8 dummy
	INCLUDE "PlotStuff.inc"	
	write(*,*)ColorFileName
	open(16,file=ColorFileName,status = 'OLD')
	read(16,*)AWEnumColors
	read(16,*)dummy
	do 10 i = 1,AWEnumColors
!	read(16,101)AWEColorName(i),AWEColorNumber(i),RGBDecimal,R !,G,B !,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
	read(16,*)AWEColorName(i),AWEColorNumber(i),R,G,B,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
!	Write(*,*)AWEColorName(i),AWEColorNumber(i),RGBDecimal,R,G,B,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
101	format(A24,Z,4I,4F)
10	continue
	close (16)
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PickAWEColor(myColor)
    	use AWE_Interfaces
    	implicit none
    	integer :: result, myColor,i
	INCLUDE "PlotStuff.inc"	
!
! 	dialog and item types
    	type(AWE_FormDialog)     :: ColorDialog	        ! the dialog
    	type(AWE_FormLabel)      :: label           ! static text
    	type(AWE_FormComboBox)   :: comboBox        ! drop down menu

! set the dialog title and create it
    	ColorDialog%title = "Color Picker"
    	call AWE_createDialog(ColorDialog)

! add static text
    	label%text = "Pick your color"
    	call AWE_addToDialog(label, ColorDialog)

! add comboBox drop down menu
!    comboBox%title = "Select your poison"
    	allocate(comboBox%items(AWEnumColors))
	comboBox%selected = myColor			! this should set the selection at the last one chosen
	do 10 i = 1,AWEnumColors
	comboBox%items(i) = trim(AWEColorName(i))
10	continue	
    	call AWE_addToDialog(comboBox, ColorDialog)

! show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(ColorDialog)
    	if (result == 1) then
      		print *,"You picked...            ", comboBox%items(comboBox%selected)
		myColor = comboBox%selected
		endif
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!call symb(XXX,YYY,3,3)
	Subroutine symb(Canvas,x,y,dX,dY)
!	plots a simple symbol on the graphic at position x,y
!	size of symbol is Xsize = abs(xmax-xmin)*dX/100.		(dX is the % of the limits
!	size of symbol is Ysize = abs(ymax-ymin)*dy/100.		(dX is the % of the limits
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	INCLUDE "PlotStuff.inc"	
	real*4 x,y,dX,dY,xSize,ySize
	
	xSize = (abs(xmax-xmin)*dX/100.)/2.		! tic length should be dX or dY in % of axes
	ySize = (abs(ymax-ymin)*dy/100.)/2.		
	call plot(Canvas,x+xSize,y,0)
	call plot(Canvas,x-xSize,y,1)
	call plot(Canvas,x,y+ySize,0)
	call plot(Canvas,x,y-ySize,1)
	return
	end
