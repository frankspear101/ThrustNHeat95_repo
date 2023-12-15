! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine ReadPositionFile(modelsteps,isostasy_flag,itopt,gofast)

!     Routine to write original positions of each grid point to make it easier to locate rocks
      implicit none
	Include "ThrustnHeat.inc"

! **********************************************************

      integer*4 i,i1,j,k,xpix,ypix,x,y,modelsteps,modelstart,posFilLen,id,xpix2
      integer*4 isostasy_flag,itopt,gofast,numberOfTimeSteps
      real*4 xdist,ydepth,saveTimeMa
      character*100 posFil,fileBase
      character*32 chnum
      integer*2 imageout(1000,1000)
      common /imagine/imageout

      REAL*4 doldFile(xGridMax,yGridMax),xgrid_oldFile(xGridMax,yGridMax)
      COMMON /filearrays/doldFile,xgrid_oldFile

! ************************************

!      positionFile=trim(modelinFull)//".PosFil.txt"    
!      open(23,file=tempFile,FORM='UNFORMATTED',status='unknown')
!      open(28,file=positionFile,FORM='formatted',status='old')
      open(28,file='',FORM='formatted',status='old')
	inquire(28,name = posFil)
	posFilLen = len(trim(posFil))
	do 5 i = posFilLen,1,-1
		if(posFil(i:i+6).eq.'PosFile')then
		fileBase = posFil(1:i-1)
!		modelIn = filBase(1:i-2)
		go to 6
		endif
5	continue
6	continue
	write(*,*)'posFil ',posFil
	write(*,*)'filBase',fileBase
	write(*,*)'modelIn',fileBase
	read(28,*)gridX,gridZ,nvert,gridbottom,exaggeration,saveTimeMa,modelSteps,numberOfTimeSteps
      do 10 i = 1,gridBottom       ! do a row
!         do 20 k=1,nvert      ! loop through each column
		read(28,*)(xGrid_OldFile(k,i),dOldFile(k,i),k=1,nvert)
!20       continue
10    continue
	do 12 i = 1,modelSteps
	read(28,*)i1,singleRockIndex(i)
12	continue

      close(28)
!      write(*,*) 'Input time in Ma between binary image files (saveTimeMa)'
!      read(*,*)saveTimeMa

!---------------------------------------------------------------
15	continue
	write(*,*)' Input X and Y coordinates for the point'
	write(*,*)' 0,0 to exit and return'
	read(*,*)x,y
	if(x.eq.0.and.y.eq.0)return
	x = int(float(x)/float(exaggeration))
	xpix = xGrid_OldFile(x,y)
	modelstart = 1
!	if(xpix.le.0)then
!		modelstart = modelSteps + xpix
!		write(*,*)'Xposition is off grid to left.'
!		write(*,*)' You must start rock at model step',modelstart
!		endif
	xdist = xpix*gridX
	ydepth = dOldFile(x,y)
	ypix = ydepth/gridZ
!	ypix = y
	write(*,*)'Dist   ',xdist,ydepth
	write(*,*)'Pixels ',xpix,ypix
	write(*,*)'Msteps,NumTimesteps',modelSteps,numberOfTimeSteps

!	Set up to track this "rock" through the model
!	First add the rock to the appropriate arrays
	nrocks = 0
	singleRockX(1) = xpix
	singleRockZ(1) = ypix + sealevel

	gofast = 1
	isostasy_flag = 0
	itopt = 4
	
	surface(0)=seaLevel
      do 701 i=1,nvert
	surface(i)=seaLevel
	baseOfCrust(i) = referenceBaseOfCrust_gridZ + sealevel
	crustThickness(i) = baseOfCrust(i) - surface(i)
	elevation(i) = 0
701   continue
	! set initial thrust depths to values from model file
	do 1611 k = 1,nthrust
	do 1610 j = 1,nvert
	dThrust(k,j) = dThrustModelFile(k,j)
1610	continue
1611	continue	
	nThrust = 0
	mct = 1
	call Auto2D               ! plot grid to graphics window
	pause 'Start -- Take a look at grid then hit return'
!**********************************************
	do 200 i = 1,modelSteps
!**********************************************

	mct = i
	nThrust = nThrustModelStep(i)
	if(i.gt.1)then
		singleRockZ(i) = singleRockZ(i-1)
		singleRockX(i) = singleRockX(i-1)
		endif

	call MoveThrusts(i,isostasy_flag,itopt)

!      write(*,*)'Plotting grid to graphics window'
!	if(goFast.eq.1)call Auto2D               ! plot grid to graphics window

200	continue

	call Auto2D               ! plot grid to graphics window

!	do 4321 i = 1,numberOfTimeSteps
	do 4321 i = 1,modelSteps
	write(13,*)i,singleRockX(i),singleRockZ(i)
4321	continue
	pause 'End -- Take a look at grid then hit return'


!     now open and read the binary image file
      pixelNumber = nvert*exaggeration
!	do 250 i = modelstart,numberOfTimeSteps
	do 250 i = 1,modelSteps
	if(singleRockIndex(i).eq.0)then
		go to 250
		endif
	tempNum = singleRockIndex(i)	
	if(singleRockX(i).le.0)then
		tRock(1,i) = 0
		dRock(1,i) = 0
		go to 251		
		endif
      write(chNum,'(i4)') tempNum
      tempFile=trim(modelinFull)//"."//trim(chNum)    
!	tempFile = trim(fileBase)//trim(chNum)
!	write(*,*)'Image file name ',tempFile
      open(23,file=tempFile,FORM='UNFORMATTED',status='old',ACCESS="Transparent")
      do 221 i1 = 1,gridBottom
      read(23)(imageout(j,i1),j=1,pixelNumber)
221   continue
      close(23)

	xpix = singleRockX(i)		!vertical section for this rock
	ypix = singleRockZ(i)
!	extract the temperature
	xpix2 = xpix*exaggeration
	tRock(1,i) = imageout(xpix2,ypix)
!	calculate the depth
!	to do this we must first figure out where the surface is (TC = 0 at surface)
	id = 0
	do 260 i1 = ypix,1,-1
	if(imageout(xpix2,i1).le.0)go to 251
	id = id + 1
260	continue
	write(*,*)'Did not find the surface. i = ',i
	pause 'hit return to continue'	
251	continue
	dRock(1,i) = float(id)*gridZ
	write(13,*)i,xpix,ypix,tRock(1,i),dRock(1,i)
!	pause 'Look at output window and hit return to read the next file'


250	continue


	go to 15
      end

