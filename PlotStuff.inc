!***********/PLOTTED/***********************************************
!	X-Y plotting variables
      INTEGER*4 nxstep,nystep,nxdec,nydec,NXPLT,NYPLT,idraw
      REAL*4 xor,yor,xmin,ymin,xmax,ymax,xlen,ylen,xconv,yconv
      REAL*4 pxor,pyor,pxlen,pylen,pxconv,pyconv
      CHARACTER*80 XLAB,YLAB,Pltitle
!      integer*4 AWEcolor
      common /PLOTTED/nxstep,nystep,nxdec,nydec,nxplt,nyplt,idraw,  &
          xor,yor,xmin,ymin,xmax,ymax,xlen,ylen,xconv,yconv,  &
		pxor,pyor,pxconv,pyconv,pxlen,pylen,         &				! Postscript (AI) values
		xlab,ylab,pltitle
!*********/AWE and PS Colors/************************
	Integer AWEnumColors,AWEColorNumber(200),RGB(200,3)
	Integer CurrentColor		! this is the color number to use for plot and fill routines
	real*4 CMYK(200,4)
	character*24 AWEColorName(200)
	Character*128 ColorFileName
	common /ALLcolors/ColorFileName,CurrentColor,AWEnumColors,AWEColorNumber,AWEColorName,RGB,CMYK
!*********/PS/************************
      character*32 psFileName
      integer jps,psFile,psPref,icolor
!           icolor is the color of the line to be plotted in the graphics window
      integer*4 psUnit	!postscript temporary file unit number
!      real*4 cmyk(100,4)
      parameter (psUnit=51)
      COMMON /ps/jps,psFile,psFileName,psPref,icolor
!*************************************
