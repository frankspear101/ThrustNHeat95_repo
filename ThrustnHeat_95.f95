      PROGRAM ThrustnHeat

!     NUMERICAL SOLUTION OF CRUSTAL THICKENING
!     To sumulate the emplacement of thrusts with lateral variation
!     in thermal conditions
!     Model is a serial 1-dimensional model

!     Coded by F. S. Spear

	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: XYPlot
	TYPE(AWE_CanvasPen) :: pen
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_Line) :: line(1)

! *****************************************
	INCLUDE "PlotStuff.inc"				! in folder "AWE_Subs"
! *****************************************
	Include 'ThrustnHeat.inc'

! ---------------------------------------
!      local variables     
	Character*16 programName
      Character*25 WORDS
      Character*32 CONFIG
      integer*4 k,i,ii,j,jj,itherm,iter,iopt,singleStep,it,stepsLikeThis,stepLoop
      real*8 timeintMa,TimeLeft,dTimeMaMax,z,eros,rtot,ddd,showTime
	real*8 saveTimeMa,timeSinceSaveMa
      integer*4 iok,itopt,isyplot,status,isostasy_flag
      integer*4 numY,iheat,ialsitp,icontinue
      real*8 g,x,y,delH,dmdp,rocp,tChange,pChange,dmdt
      integer*4 niter,relaxOnOff
      real*8 tRate,ZNO,Rvalue

      integer*4 modelSteps,istart,iend,goFast
      real*8 modelTime,totalmodeltime

!      LOGICAL*4 stdfil,iokL
!      integer*2 vref

      data words/' Rock paths position =   '/

! ------Start program here-------------------

! -------Block to initialize windows and menus-------------

      PSFILE = 0

! -----Give the program name----------------
!      call dlog(2450)


! open configuration (preferences) file
!      CONFIG='ThrustnHeat.fig'
!      open(3,file=CONFIG,status='old')

!     read default window sizes
!      read(3,*) titleLine      
!      read(3,*)W1Ylow,W1Xlow,W1Yhigh,W1Xhigh
!      read(3,*)W2Ylow,W2Xlow,W2Yhigh,W2Xhigh
!      read(3,*)W3Ylow,W3Xlow,W3Yhigh,W3Xhigh

! ----- open windows------------------------
!     Now open windows
!      Call fss_init('Command',W1Xlow,W1Ylow,W1Xhigh,W1Yhigh)
!      Call fss_openwindow(12,1,'Graphics',W2Xlow,W2Ylow,W2Xhigh,W2Yhigh)
!      Call fss_openwindow(13,0,'Output',W3Xlow,W3Ylow,W3Xhigh,W3Yhigh)
!      Call fss_openwindow(14,1,'TD plot',W3Xlow+5,W3Ylow+5,W3Xhigh,W3Yhigh)
!      CALL fss_makeplotmenu
!      Call fss_selectwindow(1)
!      Call fss_setport(1)

	OPEN(13,FILE = 'OutPut',ACCESS = 'window, 400, 600')


! -----Store startup volume number--------
!      call volget(volRefNum)        ! get reference to starting directory

! ----program starts here


	ColorFileName = 'ThrustnHeatColors.txt'
	call ReadColorFile		! set up AWE and PS colors
	CurrentColor = 1		! default = Black


! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!     MAIN LOOP

!     iCount COUNTS TOTAL ITERATIONS
      iCount=0
!     eros COUNTS TOTAL EROSION
      eros=0
!     rtot counts total R value
      rtot=0

	xlen=30
	ylen = 20


3000  continue
! 	iter COUNTS ITERATIONS FOR OUTPUT TO DISK or plotter
      iter=0

! 	MAIN LOOP MENU
      WRITE(*,*)' '
      WRITE(*,*)' Main loop menu'
      WRITE(*,*)' 0 = Return'
      write(*,*)' 1 = start a model'
      write(*,*)' 2 = Auto run (1+2+3 together)'

      write(*,*)' 4 = Save 2D to binary file for IMAGE program'
      write(*,*)' 5 = Set scale for graphics window plot '

      WRITE(*,*)'10 = Plot current geotherm(s)'
      WRITE(*,*)'11 = Plot current rock T-t paths'
      WRITE(*,*)'12 = Plot 2-D section'
      WRITE(*,*)'13 = Output status of rocks'
      WRITE(*,*)'14 = Write a DEBUG file entry'
      WRITE(*,*)'15 = Plot STEADY STATE geotherm(s)'
      WRITE(*,*)'16 = Save geotherms to disk'
      WRITE(*,*)'17 = write rock P-T paths to disk'
      WRITE(*,*)'18 = Plot axes'
      write(*,*)'19 = Delta H reaction'
      write(*,*)'20 = List model parameters'
      write(*,*)'21 = Read position file to locate rock positions'
      write(*,*)'22 = List singleRock array'

      read(*,*)iopt

      if(iopt.eq.0)go to 3000
!---------------------------------------------------------

	if(iopt.eq.22)then
		do 4321 i = 1,mct
		write(13,*)i,singleRockX(i),singleRockZ(i)
4321		continue
		pause 'hit return to continue'
		go to 3000
		endif

!---------------------------------------------------------

	if(iopt.eq.21)then
	goFast = 1
	isostasy_flag = 0
	itopt = 4
	call readPositionFile(modelsteps,isostasy_flag,itopt,gofast)
	go to 3000
	endif

!---------------------------------------------------------

	if(iopt.eq.5)then
	write(*,*) 'xlen,ylen = ',xlen,ylen
	write(*,*)' input new values'
	read(*,*)xlen,ylen
	
	go to 3000
	endif
!---------------------------------------------------------

	if(iopt.eq.1)then
	

1     continue
!     surface temperature
      tSurface=20.
	tSky =    0.
	iopt = 0
	
      tempNum=0
      do 23 k=1,xGridMax			! vsectionMax
        do 23 i=1,yGridMax
          tC(k,i)=tSky
          rho(k,i)=0.0
          cp(k,i)=0.0
          kond(k,i)=0.0
          agen(k,i)=0.0
          kappa(k,i)=0.0
          told(k,i)=0.0
          dold(k,i)=0.0
          xGrid_old(k,i)=0.0
23    continue

!     default values
	nThrust = 0
      Rvalue = .1
      niter=1     ! number of iterations for thermal relaxation
      tRate=1     ! thrusting rate in cm/year
!     convert to meters per Ma
      tRate = tRate * 1.e4

!     erosion constants
      uconst=100.       ! meters/million years
      ZNO=0.            ! depth to zero strain
      uz=0.             ! strain rate below ZNO

!     Average density for P calculations
      density=2750.
      
!     open INPUT MODEL file:
!        This file contains names of i/O files
!        Plotter information
!        Model input data
!	call FSS_StdFile(2,modelinFull,directory)

!       open(50,file='',status='old', ACTION='READ', IOSTAT=status)
!       open(50,file=modelinFull,status='old', ACTION='READ', IOSTAT=status)
       open(50,file='',status='old', ACTION='READ', IOSTAT=status)
	if(status.ne.0)stop
        INQUIRE (50, NAME=modelin)
        modelinFull = modelin
	write(*,*)'Modelin',modelin
	write(*,*)'ModelinFull',modelinfull
	write(*,*)'directory',directory


!     read name of output data file
      read(50,'(A)')programName
      if(programName.ne."ThrustNHeat7")then
      	write(*,*)' This input file is an old one (for ThrustNHeat5 or 6)'
      	write(*,*)' Program will abort - please update file for new format'
		pause 'hit return to exit'
		stop
		endif
      read(50,'(A)')outfile
!      open(10,file=outfile,status='unknown')
      WRITE(*,'(A)')modelin

!     read plotter information
	write(*,*)' Reading plotting info'
      READ(50,'(A)')
      READ(50,'(A)')pltitleTD
      READ(50,'(A)')XLABTD
      read(50,*)xorTD,xminTD,xmaxTD,xlenTD,nxticTD,nxdecTD
      READ(50,'(A)')YLABTD
      read(50,*)yorTD,yminTD,ymaxTD,ylenTD,nyticTD,nydecTD

!     Read model input file
	write(*,*)' Reading model info'

      read(50,*)		!read a row of stars
      read(50,*)
      read(50,*)seaLevel	! starting grid point for surface of earth
!     Read Grid spacing
      READ(50,*)gridZ		! vertical grid spacing

!     Read the number of vertical sections nvert
      READ(50,*)nvert		! number of horizontal grid points (number of vertical sections)

! 	first thing is to set the surface of the domain at the grid point of
! 	nbegin. So, in fact, surface is the data of off-set of the surface

      do 701 i=1,nvert
        surface(i)=seaLevel
701   continue

!     Read the spacing between vertical sections
!     Note that each thrusting increment must correspond to the distance
!     between vertical sections.
!     If you need small thrusting increments, then you will need a small
!     distance between vertical sections and a large number of vertical
!     sections.
      READ(50,*)gridX		! Horizontal grid spacing in meters

      timeintMa = gridX/tRate      ! default time interval

!	read reference base of crust (for isostasy calculations)
      READ(50,*)referenceBaseOfCrust_meters
      referenceBaseOfCrust_gridZ = referenceBaseOfCrust_meters/gridZ
      READ(50,*)isostasyFactor		!crust density/mantle density
      READ(50,*)kappaFofT			!switch for kappa as a function of T (Whittington et al 2009) 0 = no; 1=yes

! Read Base of model (in meters)
      READ(50,*)baseOfModel
      gridBottom = int(baseOfModel/gridZ + .01) + seaLevel
!	the grid is set up to start at seaLevel for the surface conditions 
! 		and go to gridBottom
! 	Note that point surfact(k) contains the surface boundary condition,
! 		which is usually constant tSurface=0
	nRocks = 0		! Initialize nrocks
      read(50,*)		! read a row of stars or dashes
      read(50,*)		! read Header for crustal layers

	write(*,*)' Reading layer info'
	do 5 j = 1,nvert

! read: #  llower qstar noLay zLLay(1) zuLay(1) rhoLay(1) cpLay(1) kLay(1) AgenLay(1)  zLLay(2) zuLay(2) rhoLay(2) cpLay(2) kLay(2) AgenLay(2)
	read(50,*)k,LLower(j),qStar(j),noLay(j),(zLLay(j,i),zuLay(j,i),rhoLay(j,i),cpLay(j,i),kLay(j,i),agenLay(j,i),i=1,noLay(j))
	k = j

! 	convert basal heat flow to watts/meter
      qStar(k)=qStar(k)/1000.
!     compute thermal diffusivity (kappa) for each layer (m**2 Ma-1)
      do 13 i=1,noLay(k)
13    kappaLay(k,i)=(kLay(k,i)/(rhoLay(k,i)*cpLay(k,i)))*secMa
!	Set base of crust for each vertical section - units are grid points
	baseOfCrust(k) = referenceBaseOfCrust_gridZ + sealevel
	crustThickness(k) = baseOfCrust(k) - surface(k)
	elevation(k) = 0
5	continue


	write(*,*)' Reading rock info'
      Read(50,*)        ! read a line of stars
	do 22 i = 1,nModelMax
	rocksToAdd(i) = 0			! zero out array
22	continue
!	write(13,*)'     Read rocks to add'
      READ(50,*)		! dummy line
2034	continue
      read(50,*)i,rocksToAdd(i)		! i is the model step for adding these rocks
	if(i.eq.0)go to 2033
      do 17 k=1,rocksToAdd(i)
	read(50,*) rocksToAdd_rockGridX(i,k),rocksToAdd_dRock(i,k),rocksToAdd_rockColor(i,k)
17    continue    
      do 18 k=1,rocksToAdd(i)
      rocksToAdd_rockGridZ(i,k)=int(rocksToAdd_dRock(i,k)/gridZ + .01) + seaLevel
18    continue    
	go to 2034

2033	continue

	write(*,*)' Reading model steps'
	
	modelSteps = 0
	totalmodeltime = 0
4024	continue
!     Read the model step parameters
      read(50,*,end=4025)        ! read a line of Stars
      read(50,*,end=4025)stepsLikeThis        ! read the number of steps of this type - 
	write(*,*)' Steps Like This = ',stepsLikeThis
	
	do 4026 stepLoop = 1,stepsLikeThis
	modelSteps = modelSteps + 1
	i = modelSteps
	write(13,*)'Reading modelStep = ',i
	write(13,*)'     Read relax time'
	if(stepLoop.eq.1)then
		read(50,*)relaxTime(i)
		else
		relaxTime(i) = relaxTime(i-1)
		endif
      totalmodeltime = totalmodeltime + relaxTime(i)

	if(stepLoop.eq.1)then
	      read(50,*,end=4025)        ! read a line of dashes
		endif
	write(13,*)'     Read add a fault'
	if(stepLoop.eq.1)then
	      read(50,*)addAFault(i)
		else
		addAFault(i) = addAFault(i-1)
		endif

      if(addAFault(i).eq.1)then
	      nThrust = nThrust + 1
 		read(50,*)FaultReference(nThrust)		! Read whether fault is referenced to sealevel  or surface
   ! 		Read depth and thickness and shear stress (in bars)of each thrust
		if(stepLoop.eq.1)then
      			read(50,*)        ! This first line just lists the vertical sections
			endif
		write(13,*)'          Read depth of thrust'
		if(stepLoop.eq.1)then
      		read(50,*)(dThrustModelFile(nThrust,j),j=1,nvert)
			else
			do 16 j=1,nvert
			dThrustModelFile(nThrust,j) = dThrustModelFile(nThrust-1,j)
16			continue
			endif
		write(13,*)'          Read shear stress'
		if(stepLoop.eq.1)then
      		read(50,*)shear(nThrust)
			else
			shear(nThrust) = shear(nThrust-1)
			endif
! 		convert shear stress to Pa (j/m**3)
      	shear(nThrust)=shear(nThrust)*1.e5
		endif
	
	nThrustModelStep(i) = nThrust
   
      
	write(13,*)'          Read move thrust'
	if(stepLoop.eq.1)then
	      read(50,*)words,((thrustTopBot(i,j),thrustDirection(i,j)),j=1,nThrust)
		else
		do 19 j=1,nThrust
		thrustTopBot(i,j) = thrustTopBot(i-1,j)
		thrustDirection(i,j) = thrustDirection(i-1,j)
19		continue
		endif
		if(stepLoop.eq.1)then
		      read(50,*,end=4025)        ! read a line of dashes
			endif

	write(13,*)'     Read pluton switch'
	if(stepLoop.eq.1)then
	      	read(50,*)plutonSwitch(i)
		else
		plutonSwitch(i) = plutonSwitch(i-1)
		endif
      if(plutonSwitch(i).eq.1)then
		if(stepLoop.eq.1)then
	            	read(50,*)plutonT(i),plutonCp(i)
!	 		Read pluton position in meters (from left, top) and convert to pixel position
			read(50,*)plutonUpperXm(i),plutonUpperZm(i)
 			read(50,*)plutonLowerXm(i),plutonLowerZm(i)
			else
	            	plutonT(i) = plutonT(i-1)
	            	plutonCp(i) = plutonCp(i-1)
!	 		Read pluton position in meters (from left, top) and convert to pixel position
			plutonUpperXm(i) = plutonUpperXm(i-1)
			plutonUpperZm(i) = plutonUpperZm(i-1)
 			plutonLowerXm(i) = plutonLowerXm(i-1)
 			plutonLowerZm(i) = plutonLowerZm(i-1)
			endif
		! convert pluton position from km to pixels
		plutonUpperX(i) = INT(plutonUpperXm(i)/gridX + .01)
		plutonLowerX(i) = INT(plutonLowerXm(i)/gridX + .01)
		plutonUpperZ(i) = INT(plutonUpperZm(i)/gridZ + .01) + surface(plutonUpperX(i))
		plutonLowerZ(i) = INT(plutonLowerZm(i)/gridZ + .01) + surface(plutonLowerX(i))
            	endif

	write(13,*)'     Read erosion'
	if(stepLoop.eq.1)then
	      	read(50,*,end=4025)        ! read a line of dashes
	      	read(50,*)erosionSwitch(i),ErosionFactor(i)        ! read erode header
		else
		erosionSwitch(i) = erosionSwitch(i-1)
		erosionFactor(i) = erosionFactor(i-1)
		endif
!	if(erosionSwitch(i).eq.1)then
!		if(stepLoop.eq.1)then
!		      read(50,*)        ! read vertical section list
!		      read(50,*)(erode(i,j),j=1,nvert)
!			else
!			do 14 j = 1,nvert
!			erode(i,j) = erode(i-1,j)
!14			continue
!			endif
!		endif

4026 	continue		! end of stepsLikeThis loop

	go to 4024		! go back and read the next step

4025  continue		! done reading model steps


! 	Set Total increment counter

      mct=1

	singleRockZ(mct) = 200 
	singleRockX(mct) = 1
	singleRockTime(mct) = 0
      close(50)


	! set initial thrust depths to values from model file
	do 1611 k = 1,nthrust
	do 1610 j = 1,nvert
	dThrust(k,j) = dThrustModelFile(k,j)
1610	continue
1611	continue	
	nThrust = 0    ! there are no thrusts active until the first model step is run

	do 1612 i = 1,modelSteps
	singleRockIndex(i) = 0
1612	continue
				! this may screw things up
! 	Assign grid points to respective layers for each vertical section
      do 21 k=1,nvert
!      write(13,*)k,' = working on vertical section'
      ii=surface(k)
      rho(k,ii)=rhoLay(k,1)
      cp(k,ii)=cpLay(k,1)
      kond(k,ii)=kLay(k,1)
      kappa(k,ii)=kappaLay(k,1)
      agen(k,ii)=agenLay(k,1)

! 	here i is the number of space

      do 20 ii = surface(k)+1,gridBottom
      if(ii.gt.yGridMax) then
        	write(*,*) 'grid point is too large'
        	write(*,*) 'increase the array'
        	pause
      	endif
      z=gridZ*float(ii - surface(k))
      do 25 j=1,noLay(k)
      if(z.gt.zLLay(k,j).and.z.LE.zuLay(k,j))then
! 		if here, then the grid point is in layer j
        	rho(k,ii)=rhoLay(k,j)
        	cp(k,ii)=cpLay(k,j)
        	kond(k,ii)=kLay(k,j)
        	kappa(k,ii)=kappaLay(k,j)
        	agen(k,ii)=agenLay(k,j)
        	go to 20
      		ENDIF
25    	CONTINUE
      	write(13,*)k,j,noLay(k),' = vert section,layer,noLay(k) (no layer is found)'
      	pause
20    	continue
21	continue


56    	continue
! 	COMPUTE INITIAL TEMPERATURES
      	itherm=1
      	if (itherm.eq.1) then
        	do 65 k=1,nvert
        	if(kappaFofT.eq.0)then
	            	CALL GTSTEADY(k)
			else
			call GTSTEADY_K_of_T(k)
			endif
            	do 60 j=surface(k),gridBottom
60          	tC(k,j)=tSteady(j)
65      	continue
        	go to 90
      		endif

      if(itherm.ne.0.or.itherm.ne.1.or.itherm.ne.2) goto 56

90    continue
          

!     set initial depth and temperature in arrays told and dold 
!     (these store old depth and T for each V section for use in
!     H reaction calculations
      do 95 k=1,nvert
      do 95 j = surface(k),gridBottom
      told(k,j)=tC(k,j)
    	dold(k,j)=FLOAT(j - surface(k))*gridZ
      xGrid_old(k,j)=k
95	continue

      call SetTemp


      rockTime(1)=0.       ! time counter for rock paths

!     Write model parameters to text window
      	write(*,*)' Writing model parameters to text window'
      	modelTime = 0           !  here model time is a counter for the output
!      	call txend
!	go to 4011
      	do 4010 i=1,modelSteps
      	write(13,*)
      	WRITE(13,*)' Model step = ',i

      	modelTime = modelTime + relaxTime(i)
      	if(relaxTime(i).le.0)then
            tRate=0
            else
            tRate = (gridX/relaxTime(i))/10000.
            endif
      	WRITE(13,*)' Relax Time (Ma)   tRate(cm/y)  modelTime '
      	write(13,*)relaxTime(i),tRate,modelTime
	write(13,*)rocksToAdd(i),'    rocks to add'
	if(rocksToAdd(i).gt.0)then
		do 4006 j = 1,rocksToAdd(i)
		write(13,*) rocksToAdd_rockGridX(i,j),rocksToAdd_dRock(i,j)
4006		continue
		endif
	nThrust = nThrustModelStep(i)
      	write(13,*)' Thrust number   ',(j,j=1,nThrust)
      	WRITE(13,*)' Thrust direction = ',(thrustTopBot(i,j),thrustDirection(i,j),j=1,nThrust)

      	if(plutonSwitch(i).eq.0)then
		 write(13,*)' No pluton this model step'
		 else
		 write(13,*)' Pluton this step'
		 write(13,*)'Pluton T    Pluton cp'
		 write(13,*)plutonT(i),plutonCp(i)
		 write(13,*)'Pluton upper left  x,z ',plutonUpperX(i),plutonUpperZ(i)
		 write(13,*)'Pluton lower right x,z ',plutonLowerX(i),plutonLowerZ(i)
		 endif
!	if (erosionSwitch(i).ne.1)then
!		write(13,*) ' No erosion this step'
!		else	
      	write(13,*)' Erosion switch, ErosionFactor = ',erosionSwitch(i),ErosionFactor(i)
!      	write(13,4009)(j,j=1,nvert)
!      	write(13,4009)(erode(i,j),j=1,nvert)
!4009  	format(50I4)
!		endif
		
4010  	continue
4011	continue
      	write(13,*)'   '
      	exaggeration = int((gridX/gridZ)+.01)
      	pixelNumber = nvert*exaggeration
      	write(13,*)'Number of pixels written = ',pixelNumber
      	write(13,*)'Number of lines written  = ',gridBottom

! 	Write a file that contains the image dimensions
!	outfile = trim(modelin)//'.gridsize.txt'
	outfile = trim(modelinFull)//'.gridsize.txt'
	open(24,file=outfile,status='UNKNOWN')
      	write(24,*)'Number of pixels written = ',pixelNumber
      	write(24,*)'Number of lines written  = ',gridBottom
	close(24)
!     	write out the initial grid configuration to the binary file

	nthrust = 0
      	call Save2D()      ! save initial configuration to binary file

      	call Auto2D(XYPlot,0)                ! plot grid to graphics window - (0) means open a new canvas
      	call Auto2D(XYPlot,1)                ! plot grid to graphics window - (1) means clear canvas and write new data
	
	
	
	endif
	

!---------------------------------------------------------
! 	-------IOPT=4: save 2D grid ----------------------------
         if(iopt.eq.4)then

! 	Routine to save 2D grid to binary file for input to program IMAGE
! 	set the volume reference number to the folder where we opened the input file
!         call VolSet(volRefNum)
         call Save2D()
         go to 3000
         endif



! ***********#################**********************************########################
! ***********#################**********************************########################
! ***********#################**********************************########################

! --------------auto thrusting and relaxing------------------------
      if(iopt.eq.2)then
!     Auto thrust and relax 
	iend = 1
5555	continue
      write(*,*)' Run the model'
      write(*,*)
      write(*,*)' Number of modelSteps = ',modelSteps
      write(*,*)' TotalModelTime       = ',totalModelTime
      write(*,*)'  0 = return'
      write(*,*)'  1 = excute range of steps'
      write(*,*)'  2 = run all steps'
      write(*,*)'  3 = Do ISOSTASY'
      write(*,*)'  4 = Go to specified step without saving (testing thrust model)'
      write(*,*)'  5 = Single step through thrust model (testing thrust model)'
      write(*,*)'  6 = Write Fault File'
      read(*,*)itopt
      if(itopt.eq.0)go to 3000

      if(itopt.eq.6)then
	      write(*,*)' Writing Fault File'
		call WriteFaultFile(iend)
		go to 5555
            endif

      if(itopt.eq.3)then
	      write(*,*)' Doing isostasy'
		call isostasy(itopt)
		call Auto2D(XYPlot,1)			!(1) means clear the old canvas and redraw
         	call Save2D()
		go to 5555
            endif
	

      if(itopt.eq.1)then
		goFast = 1
	      	write(*,*)' Input starting and ending steps'
		read(*,*)istart,iend
		if(istart.eq.iend)then
      			singleStep = 1
      			else
      			singleStep = 0
      			endif
            	endif

	if (itopt.eq.2)then
		singleStep = 0
		istart = 1
            	iend = modelSteps
		write(*,*)'0 = no relax; 1 = relax'
		read(*,*)relaxOnOff
		write(*,*)' Input saveTimeMa (time between saving binary images)'
		read(*,*)saveTimeMa
		WRITe(*,*)' Do you want ISOSTASY? 0 = no, 1  = yes'
		read(*,*)isostasy_Flag
		write(*,*)'Do you want to 0 = go fast; 1 = watch steps'
		read(*,*)goFast
		endif

	if(itopt.eq.4)then
	      write(*,*)' Input end step (0 to abort)'
		read(*,*)iend
		if(iend.eq.0)go to 5555
		goFast = 0
		relaxOnOff = 0
		saveTimeMa = 1000
		isostasy_flag = 0
		singleStep = 0
		istart = 1
		endif
      if(itopt.eq.5)then		! we are doing a single step from whereever we start
		goFast = 1
		istart = iend		! from previous run=4
		iend = istart
      		singleStep = 1
		relaxOnOff = 0
		saveTimeMa = 1000
		isostasy_flag = 0
		singleStep = 1
      		endif


5000  continue

      write(*,*)
      write(*,*)
      modelTime = 0           ! modelTime stores actual model time when model is running
	timeSinceSaveMa = 0.0	! counts time between saves to binary image
5556	continue		! loop here for single steps
   !**************************************
      do 5010 i=istart,iend		! Main model step loop
   !**************************************
      modelTime = modelTime + relaxTime(i)
      if(relaxTime(i).le.0)then
            tRate=0
            else
            tRate = (gridX/relaxTime(i))/10000.
            endif

	nThrust = nThrustModelStep(i)
	if(itopt.ne.4)then
	      WRITE(*,*)' Model step = ',i
		endif

5001  continue

!	iok = 0
!      call trapit(iok)
!      if(iok.eq.1)then
!            write(*,*)' User hit escape'
!            pause
!            go to 3000
!            endif


	if(rocksToAdd(i).gt.0)then
		do 5011 j=1,rocksToAdd(i)
		jj=j+nRocks
		rockColor(jj) = rocksToAdd_rockColor(i,j) 
		rockGridZ(jj,mct) = rocksToAdd_rockGridZ(i,j)
		rockGridX(jj,mct) = rocksToAdd_rockGridX(i,j)
		tRock(jj,mct) = tC(rockGridX(jj,mct),rockGridZ(jj,mct))
		do 5012 k = 1,mct-1
		rockGridZ(jj,k) = 0
		rockGridX(jj,k) = 0
		tRock(jj,k) = 0
5012		continue
5011		continue
		nRocks = nRocks + rocksToAdd(i)
		endif	


      mct=mct+1
!     set the rock grid points for this new time step
      do 5018 j = 1,nRocks
      rockGridZ(j,mct) = rockGridZ(j,mct-1)
      rockGridX(j,mct) = rockGridX(j,mct-1)
      tRock(j,mct) = tRock(j,mct-1)
5018  continue
	singleRockZ(mct) = singleRockZ(mct-1)
	singleRockX(mct) = singleRockX(mct-1)


	call MoveThrusts(i,isostasy_flag,itopt)


!      write(*,*)'Plotting grid to graphics window'
	if(goFast.eq.1)then
		call Auto2D(XYPlot,1)               ! plot grid to graphics window
		endif
	if(itopt.eq.5)then
		! loop back and do another step if desired
		Write(*,*)'  0 = continue with next step'
		Write(*,*)'  1 = continue to end'
		write(*,*)' -1 = abort single step'
		read(*,*)icontinue
		select case (icontinue)
		case(-1)
			go to 5555
		case(0)
			istart = istart+1
			iend = iend + 1
			if(iend.eq.modelsteps)then
				write(*,*)' All done with model steps'
				go to 5555
				endif
			go to 5556
		case(1)
			istart = istart+1
			iend = modelsteps
			gofast = 0		! don't write out until the model is done
			itopt = 6		! this will allow loop to keep executing until we hit the end
			go to 5556
		case default
			write(*,*)'Bad choice'
			pause 'hit return to continue'
			go to 5555
		end select
		endif			


	if(relaxOnOff.eq.0)then
	      if(itopt.ne.4)then
		      call Save2D()      ! save this result to binary file
			endif
		go to 5201
		endif



	if(singleStep.eq.1) go to 5555
!     relax all vertical sections
!     compute the time increment for the given R value
      dTimeMa = (1./(gridZ*gridZ) + 1./(gridX*gridX))
      dTimeMaMax = Rvalue/(1.e-6*secMa*dTimeMa)

!     The problem must be solved for a total time timeintMa
!     with time steps no bigger than dTimeMa
      dTimeMa  = dTimeMaMax
      timeintMa = relaxTime(i)
      timeLeft = timeintMa
	showTime = 0.
5100  continue
      if(timeLeft.le.0)go to 5200

      if(dTimeMa.le.timeLeft)then
            dTimeMa = dTimeMaMax    ! use the maximum on this increment
            TimeLeft = TimeLeft - dTimeMa
            else        
            dTimeMa = TimeLeft      ! we must use a smaller increment
            TimeLeft = 0.
            endif

      iCount = iCount + 1

	showTime = showTime + dTimeMa
	if(showTime.gt.0.1)then
		showTime = 0.
		endif
	
      call RxExp1

      Call HeatGen

	timeSinceSaveMa = timeSinceSaveMa + dTimeMa
	if(timeSinceSaveMa.gt.saveTimeMa)then
	      call Save2D()      ! save this result to binary file
		timeSinceSaveMa = 0.0
		endif

      go to 5100
5200  continue

5201	continue		! skip to here if relaxOnOff = 0


!     Set new temperatures for rocks in array tRock
      CALL SetTemp

!     set rockTime
      rockTime(mct)=modelTime

5010  continue

!	write out original grid positions into two files
	call WritePositionFile(ModelSteps,saveTimeMa)

!     write out rock paths, if desired

	if(itopt.eq.1)go to 5555
	if(itopt.eq.4)then
		call auto2D(XYPlot,1)
		call save2D
		go to 5555
		endif

	if(goFast.eq.0)call auto2D(XYPlot,1)
	
	Write(*,*)' Specify name of file for storing rock PT paths'
	call fss_alert('Alert','Specify name of file for storing rock PT paths')
      open(45,file='',status='New',  IOSTAT=status)
	if(status.ne.0)go to 3000
      INQUIRE (45, NAME=outfile)
      call RockOut(45,modelTime)
      close(45)


        go to 3000

      endif      ! end auto thrust+relax option



! ------IOPT Plot geotherms-------------------------------------------
      if(iopt.eq.10) then
 !       CALL user(xorTD,xminTD,xmaxTD,xlenTD,yorTD,yminTD,ymaxTD,ylenTD)

        call GeoTpl(0,modelTime)

        go to 3000
        ENDIF

! ------IOPT: Plot rock paths ----------------------------
      if(iopt.eq.11)then
      write(*,*)' Do you want to plot symbols? 0 = no, 1 = yes'
      read(*,*)(isyplot)
!      CALL user(xorTD,xminTD,xmaxTD,xlenTD,yorTD,yminTD,ymaxTD,ylenTD)
!	call fss_openpicture(4)
      call ROCKPLOT(nThrust,isyplot)
!	call fss_closepicture(4)
      go to 3000
      endif
! ------IOPT: Make 2-D plot of grid on graphics window ----------------------------
      if(iopt.eq.12)then
      call Plot2D()
      go to 3000
      endif
! ------IOPT: Output current status of rocks ----------------------------
      if(iopt.eq.13)then

      write(13,*)' '
      WRITE(13,*)' Current status of ROCKS'
      do 3505 j=1,nRocks
      write(13,*)' '
      WRITE(13,'(A,I6)')WORDS
      write(13,*)'Rock #  Time Depth(m)  P(kb)     T'
      do 3505 i = 1,mct
      write(13,3506)j,rockTime(i),dRock(j,i),dRock(j,i)*density*9.8/1.e8,tRock(j,i)
3506  format(i5,f12.5,f12.2,f12.3,f12.1)
3505  continue

      go to 3000
      endif
! ------IOPT: save a DEBUG file entry ----------------------------
      if(iopt.eq.14)then
3612   continue
      WRITE(*,*)' Debugging output routine'
      WRITE(*,*)' Current vertical sections are:'
      WRITE(*,*)'    Twigland      Rootzone'
      WRITE(*,*)1,nvert
      WRITE(*,*)' 0 = Exit'
      WRITE(*,*)' k = values of vertical section'
      read(*,*)(k) 
      if(k.eq.0) go to 3000
      WRITE(*,*)'  Point Depth   ToC     rho     cp      k       kappa           agen'
      do 3610 i=surface(k),gridBottom
      ii=i-surface(k)
      WRITE(*,1056)i,gridZ*(ii),tC(k,i),rho(k,i),cp(k,i),kond(k,i),kappa(k,i),agen(k,i)
3610  continue
1056  format(I5,f8.0,f12.1,2f8.0,f8.2,2E10.3,f8.3)
      WRITE(*,*)' LLower     qStar  '
      WRITE(*,*)LLower(k),qStar(k)
      write (*,*)' Hit return when ready'
      read(*,*)(k)
      go to 3612
      endif
! ------IOPT=15: Plot a STEADY STATE geotherm ----------------------------
      if(iopt.eq.15)then
3700  continue
          WRITE(*,*)' Input number of geotherm to plot (0 to exit)'
          WRITE(*,*)' Current vertical sections are:'
          write(*,*)' Twigland.........Rootzone'
          WRITE(*,'(20I3)')(i,i=1,nvert)
          read(*,*)(k)
          if(k.eq.0)go to 3000
!          CALL user(xorTD,xminTD,xmaxTD,xlenTD,yorTD,yminTD,ymaxTD,ylenTD)
	    write(*,*)' Constant thermal diffusivity = 0; Kappa_function_of_T = 1'
	    read(*,*)it
	    if(it.eq.0)then
	        call GTSTEADY(k)
		  else
		  call GTSTEADY_K_of_T(k)
		  endif
!	    call fss_openpicture(4)
!	    call scolor(icolor)
	currentcolor = 3
          call GEOPLT(tSurface,gridZ,xmaxTD,ymaxTD,tSteady,modelTime,surface(k),gridBottom)
	currentcolor = 1
!	    call scolor(1)
!	    call fss_closepicture(4)

          go to 3700
          endif

! ------IOPT=16: pause ----------------------------
         if(iopt.eq.16)then
         call SaveGT()
         go to 3000
         endif

! ------IOPT: write rock P-T paths to disk ----------------------------
      if(iopt.eq.17)then

       	open(45,file='',status='unknown',  IOSTAT=status)
	if(status.ne.0)go to 3000
        INQUIRE (45, NAME=outfile)


!       outfile=trim(modelinFull)//'.PTt'

!       iokL = stdfil(2,vref,'Input save file name',outfile,1,'TEXT') ! 2 is setfile and set vref
!       if(.not.iokL)go to 3000
!       open(45,file=outfile,status='unknown')

!        write(*,*)' Input name of file for saving P-T path information'
!        outfile=trim(modelinFull)//'.PTt'
!        call saveas(45,iok,outfile,'Input save file name')
!        if(iok.eq.1)go to 3000
        call RockOut(45,modelTime)
        close(45)
        go to 3000
        endif
        
!  ---------IOPT Plot axes ---------------------------------------
      if(iopt.eq.18)then

!      call setplt(XORTD,XMINTD,XMAXTD,XLENTD,nXticTD,nXdecTD,XLABTD,			&
!                 YORTD,YMINTD,YMAXTD,YLENTD,nYticTD,nYdecTD,YLABTD,PltitleTD,iok)
!      call Setplt(26257,xor,xmin,xmax,xlen,nxtic,nxdec,XLAB,
!     * yor,ymin,ymax,ylen,nytic,nydec,YLAB,pltitle,iok)

!      WRITE(*,*)' Do you want to plot Al2SiO5 triple point?'
!      WRITE(*,*)' 0=NO, 1=YES'
!      READ(*,*)ialsitp

	CurrentColor = 1		! black is the default
	call Setplot(iOK)
	if(iok.eq.1)go to 3000  ! cancel button was selected

	! I need to make the plot with these parameters -- not the ones in the common block
      	!	CALL axis(nxticTD,nxdecTD,XLABTD,nyticTD,nydecTD,YLABTD,pltitleTD)
	call PlotAxes(XYPlot)

      

!	call clgs(4)
!	call fss_openpicture(4)
!     SET UP PLOTTER AXES
!      CALL user(xorTD,xminTD,xmaxTD,xlenTD,yorTD,yminTD,ymaxTD,ylenTD)
!      CALL axis(nxticTD,nxdecTD,XLABTD,nyticTD,nydecTD,YLABTD,pltitleTD)

!     draw pressure axis at x=xmax
!      pressure is calculated as

!            P = rho * g * h

! 		P(kbar) = rho * g * h(meters)/10^8

! 		h(meters)= P(kbar)*10^8/(rho*g)

! 	solve for h (depth) where P equals an integer value in kbar
! 	the variable i iterates for pressure in kbar
! 	draw tic marks first
      i=0
      g=9.8
30    continue
      i=i+1
 	xmaxTD = xmax
 	ymaxTD = ymax
      x=xmaxTD
      y=yminTD+float(i)*1.e8/(density*g)	! depth for pressure i in meters
	y=y/1000.					! depth in km
      if(y.gt.ymaxTD)go to 35
      call plot(XYPlot,x,y,0)
      call plot(XYPlot,x+10,y,1)		! I think this is supposed to draw a line at the depth in km. The original was 10 pixels. I'm not sure how long the line should be
      

!      call fss_line(10,0)

      go to 30
35    continue
!     PUT pressure NUMBERS on y AXIS
      i=0
40    continue
      i=i+5       				! label in 5 kbar increments
      x=xmaxTD
      y=   float(i)*1.e8/(density*g)      ! calculate depth to 5 kbar in meters
	y=y/1000.                           ! convert to km
      if(y.gt.ymaxTD)go to 49
      numY=i
!      CALL PLOT(x,y,0)
!      call fss_move(5,5)
      write(ppp,'(I4)')numY
!	call fss_drawstring(ppp)
	call TextOnPlot(XYPlot,x+5,y+5,ppp,12)

      go to 40

49    continue
!      if(ialsitp.eq.1)call al2sio5(1,1)

      go to 3000
      endif             ! end plot axes
! ------IOPT: HEAT of REACTION-------------------------------------
      if(iopt.eq.19)then
4100  write(*,*)' Options'
      write(*,*)' 0 = return'
      write(*,*)' 1 = List Changes in rock position'
      write(*,*)' 2 = Specify heat of reaction'
      write(*,*)' 3 = Compute T from H of reaction'
      read(*,*)iheat
      if (iheat.eq.0)go to 3000
      if(iheat.eq.1)then
      write(*,*)'  NGRID  dold-dnew told-tnew'
      do 4105 i=surface(1),gridBottom      ! this is a bug
      write(*,4110)i,(float(i-1)*gridZ-dold(k,i),tC(k,i)-told(k,i),k=1,nvert)
4110  format(1x,I5,20f10.1)
4105  continue
      go to 4100
      endif
      if(iheat.eq.2)then
      write(*,*)' Current value of delH is:',delH
      write(*,*)' Input heat of reaction (per mole of garnet produced)(e.g. 75E3 is low, 200E3 is high)'
      write(*,*)' units are kj/mole of garnet produced'
      read(*,*)delH
!     specify T and P derivatives of m garnet
!     temperature derivative (moles garnet / degC - m**3)
      dmdt = 1.0 * delH
!     Pressure derivative (moles garnet/ bar - m**3)
      dmdp = 3.E-2 * delH
!     dmdt and dmdp are now in units of j m-3 k-1 and j m-3 k-1
      go to 4100
      endif
      if(iheat.eq.3)then
      write(*,*)' k   i   tChange   pChange'
      do 4130 k=1,nvert
      do 4130 i=1+surface(k),gridBottom
      ddd=float(i-1)*gridZ
      rocp = rho(k,i)*cp(k,i)
      tChange= (dmdt*(tC(k,i)-told(k,i)))/rocp
      pChange=(dmdp*(ddd-dold(k,i))*(density*9.8/1.e5))/rocp
!     *(dmdt*(tC(k,i)-told(k,i)) + dmdp*(float(i-1)*grid-dold(k,i)))/
!     * (rho(k,i)*cp(k,i))    
      if(tChange.LT.0.)tChange=0.
      if(pChange.LT.0.)pChange=0.
      tC(k,i) = tC(k,i) - tChange - pChange
!     write(*,4135)k,i,tChange,pChange
!4135 format(2i5,2f10.2)
      dold(k,i)=ddd
      told(k,i)=tC(k,i)
4130  continue    
      endif
      go to 4100
      endif             ! end heat of reaction

      go to 3000
      END

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE ISOSTASY(itopt)
      implicit none
!	This routine computes the height of the mountains based on isntaneous isostatic compensation
!	The reference crust thickness is read in from the input file and it is assumed to be in equilibrium
!	This is the initial depth of the crust from sealevel
!	The act of thrusting thickens the crust and will depress the root
!	Here is the method:
!		Each vertical section is checked for the thickness of the crust.
!			Thickness of crust = base of crust - surface
!		The height that this vertical section should float above the mantle is calculated
!			If the crust is thickening then this usually means that the mountains are too high
!		Once the new elevation is calculated, the entire vertical section plus rocks and thrusts are moved up or down as needed.


	Include "ThrustnHeat.inc"
! ---local variables------------
      integer i,ii,j,k,move,itopt
	real*8 factor
	
!	isostasy factor
!	factor = 1.0 - 2700./3300.
	factor = 1.0-isostasyFactor
	! determine thickness of crust in each vertical section

	if(itopt.eq.3)then
	write(13,*)
	write(13,*)'---------------------------'
	write(13,*)' Initial conditions for ISOSTASY calculations'
	write(13,*)'Reference base of crust = ',referenceBaseOfCrust_gridZ
	write(13,*)'Surface grid'
	write(13,1022)(surface(k),k=1,nvert)
	write(13,*)'Elevation'
	write(13,1022)(elevation(k),k=1,nvert)
	write(13,*)'Base of crust'
	write(13,1022)(baseOfCrust(k),k=1,nvert)
	write(13,*)'CrustThickness'
	write(13,1022)((baseOfCrust(k) - surface(k)),k=1,nvert)
	endif

	do 10 k = 1,nvert
	! surface(K) and baseOfCrust are in gridpoints
	! crustThickness is in meters
	crustThickness(k) = (baseOfCrust(k) - surface(k))

	! do isostacy calculations
	elevation(k) = ((crustThickness(k)-referenceBaseOfCrust_gridZ)*factor)
	!elevation is in grid points

	!move section up or down
	move = -(elevation(k) - (sealevel - surface(k)))
	!move is in grid points

!	write(13,*)' vertSection K = ',k
!	write(13,*)'     elevation = ',elevation
!	write(13,*)'          move = ',move
	
	if(move.lt.0)then

		! moving the vertical section up - start from the top and work down
		do 20 i = surface(k),gridbottom
	        ii = i + move
	        tC(k,ii)    = tC(k,i)
	        rho(k,ii)   = rho(k,i)
	        cp(k,ii)    = cp(k,i)
	        kond(k,ii)  = kond(k,i)
	        agen(k,ii)  = agen(k,i)
	        kappa(k,ii) = kappa(k,i)
	        dold(k,ii)  = dold(k,i)
	        told(k,ii)  = told(k,i)
	        xGrid_old(k,ii) =  xGrid_old(k,i)
20		continue	
		! fill in the bottom of the vertical section with values
	! there is nothing under the grid to move up.
	! Not doing anything just leaves the last "move" grid points alone
	!      This probably is too far away from anything to matter - so give it a try.
!		do 21 ii = gridbottom-move,gridbottom
!		i = gridbottom-move-1
!      	tC(k,ii)    = tC(k,i)
!     	 	rho(k,ii)   = rho(k,i)
!        	cp(k,ii)    = cp(k,i)
!        	kond(k,ii)  = kond(k,i)
!        	agen(k,ii)  = agen(k,i)
!        	kappa(k,ii) = kappa(k,i)
!        	dold(k,ii)  = dold(k,i)
!        	told(k,ii)  = told(k,i)
!21		continue	
		elseif(move.gt.0)then
		! moving the section down. Start from bottom and work up
		do 23 i = gridbottom,surface(k)-move,-1
      	  ii = i + move
      	  if(i.ge.yGridMax.or.i.lt.1)go to 23
      	  if(ii.ge.yGridMax.or.ii.lt.1)go to 23
      	  tC(k,ii)    = tC(k,i)
      	  rho(k,ii)   = rho(k,i)
      	  cp(k,ii)    = cp(k,i)
      	  kond(k,ii)  = kond(k,i)
      	  agen(k,ii)  = agen(k,i)
      	  kappa(k,ii) = kappa(k,i)
      	  dold(k,ii)  = dold(k,i)
      	  told(k,ii)  = told(k,i)
	        xGrid_old(k,ii) = xGrid_old(k,i)
23		continue
		!fill in above the surface with zeros
!		do 24 ii = surface(k),surface(k)-move,-1
!      	  tC(k,ii)    = 0
!      	  rho(k,ii)   = 0
!      	  cp(k,ii)    = 0
!      	  kond(k,ii)  = 0
!      	  agen(k,ii)  = 0
!      	  kappa(k,ii) = 0
!      	  dold(k,ii)  = 0
!      	  told(k,ii)  = 0
!24		continue	
		else
		! move = 0 - don't do anything
		endif
		
!	move the rocks
      do 150 j=1,nRocks
      if(rockGridX(j,mct).eq.k)then        ! this one must be checked (it is on v section kk)
      	rockGridZ(j,mct) = rockGridZ(j,mct) + move
      	endif
150   continue

	! move the single rock position
	if(singleRockX(mct).eq.k)then
		singleRockZ(mct) = singleRockZ(mct) + move
		endif

! -----------------------------------------
!     Move thrusts
      do 160 j=1,nThrust
      	  thrustGridZ(j,k) = thrustGridZ(j,k) + move
160   continue

	baseOfCrust(k) = baseOfCrust(k) + move
	surface(k) = surface(k) + move
10	continue

	if(itopt.eq.3)then
	write(13,*)
	write(13,*)'Final conditions for ISOSTASY calculations'
	write(13,*)'Surface grid'
	write(13,1022)(surface(k),k=1,nvert)
1022	format(40I6)
	write(13,*)'Elevation'
	write(13,1022)(elevation(k),k=1,nvert)
	write(13,*)'Base of crust'
	write(13,1022)(baseOfCrust(k),k=1,nvert)
	write(13,*)'CrustThickness'
	write(13,1022)(crustThickness(k),k=1,nvert)
	pause 'Hit return to continue'
	endif
	
      return
      END

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE GTSTEADY(k)
      implicit none
!     THIS ROUTINE COMPUTES THE STEADY STATE GEOTHERM FOR THE GIVEN
!     SET OF BOUNDARY CONDITIONS, INTERNAL HEAT PRODUCTION
!     AND THERMAL PARAMETERS
!     THE CRUSTAL MODEL IS A NGRID LAYER MODEL

	Include "ThrustnHeat.inc"
! ---local variables------------
      integer jj,j,k
      real*8 qSum

!     DEPTH OF UPPER LAYER IS DSCALE
!     agen IS HEAT GENERATION
!     qStar IS MANTLE HEAT FLUX
!     k IS THERMAL CONDUCTIVITY
!     NOTE --- THIS ROUTINE IS DONE IN DIMENSIONALIZED PARAMETERS

      tSteady(surface(k))=tSurface       ! chen
      do 10 j = surface(k)+1, gridBottom
      qSum=qStar(k)
      do 20 jj = j,gridBottom
20    qSum=qSum+agen(k,jj)*gridZ
      tSteady(j)=(agen(k,j)*gridZ*gridZ)/(2.0*kond(k,j)) + qSum*gridZ/kond(k,j) + tSteady(j-1)
10    continue

      return
      END

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE GTSTEADY_K_of_T(k)
      implicit none
!     THIS ROUTINE COMPUTES THE STEADY STATE GEOTHERM FOR THE GIVEN
!     SET OF BOUNDARY CONDITIONS, INTERNAL HEAT PRODUCTION
!     AND THERMAL PARAMETERS
!     THE CRUSTAL MODEL IS A NGRID LAYER MODEL
!	
!	in this code, the model of Whittington is used where kappa is a function of T

	Include "ThrustnHeat.inc"
! ---local variables------------
      integer jj,j,k,newiter
      real*8 qSum,kond_of_T,kappa_of_T,tSteady_new(yGridMax),TK

!     DEPTH OF UPPER LAYER IS DSCALE
!     agen IS HEAT GENERATION
!     qStar IS MANTLE HEAT FLUX
!     k IS THERMAL CONDUCTIVITY
!     NOTE --- THIS ROUTINE IS DONE IN DIMENSIONALIZED PARAMETERS

	!first calculate the geotherm based on the input thermal conductivity
      tSteady(surface(k))=tSurface       ! chen
      do 10 j = surface(k)+1, gridBottom
      qSum=qStar(k)
      do 20 jj = j,gridBottom
20    qSum=qSum+agen(k,jj)*gridZ
      tSteady(j)=(agen(k,j)*gridZ*gridZ)/(2.0*kond(k,j)) + qSum*gridZ/kond(k,j) + tSteady(j-1)
10    continue

!	Now recalculate the thermal conductivity as a function of temperature and recompute the geotherm
!	Iterate until it converges (not sure it really does converge)


	do 50 newiter = 1,10		! do 10 iterations - more than enough to converge

	do 30 j = surface(k)+1,gridBottom
	TK = tSteady(j) + 273.15
		if(TK.le.846)then	! we are in the alpha quartz field
		kappa_of_T = 567.3/TK - 0.062
		else
		kappa_of_T = 0.732 - 0.000135*TK
		endif
!	units need to be 

	kond_of_T = kappa_of_T*rho(k,j)*cp(k,j)/1.0e6		! convert from mm^2/sec to m^2/sec

      qSum=qStar(k)
      do 35 jj = j,gridBottom
35    qSum=qSum+agen(k,jj)*gridZ
      tSteady_new(j)=(agen(k,j)*gridZ*gridZ)/(2.0*kond_of_T) + qSum*gridZ/kond_of_T + tSteady_new(j-1)
30	continue

!	write(*,*)'New Tsteady, Old Tsteady, New-Old '
	do 40 j = surface(k)+1,gridBottom
!	write(*,*)tSteady_new(j),tSteady(j),tSteady_new(j)-tSteady(j)
	tSteady(j) = tSteady_new(j)
40	continue

50	continue

	

      return
      END

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE HeatGen
!!     SUBROUTINE TO COMPUTE THE EFFECT OF HEAT GENERATION
!!     FOR EACH FINITE DIFFERENCE ITERATION
!     Note that each grid point has values of agen, k and rho

      implicit none
	Include "ThrustnHeat.inc"

      integer*4 j,k
      
      dTimeSec=secMa*dTimeMa
      do 5 k = 1, nvert
      DO 10 j = surface(k)+1,gridBottom-1
      tC(k,j)=tC(k,j)+(agen(k,j)/(rho(k,j)*cp(k,j)))*dTimeSec
10    CONTINUE
5     continue  
      RETURN
      END

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE SetTemp
      implicit none
	Include "ThrustnHeat.inc"
      integer k,j,i

!     set the temperature of the rock
      do 20 j=1,nRocks
        if(dRock(j,mct).le.0.) then
          tRock(j,mct)=0.0
          else
          k=rockGridX(j,mct)
          i=rockGridZ(j,mct)
          tRock(j,mct)=tC(k,i)
          endif
20    CONTINUE

      return
      end
! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE ROCKPLOT(XYPlot,isyplot)

	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: XYPlot
! *****************************************
	INCLUDE "PlotStuff.inc"				! in folder "AWE_Subs"
! *****************************************
	Include "ThrustnHeat.inc"

      real*4 x,y
      integer size,isymb,isyplot

      INTEGER*4 iup,j,m

      isymb = 3
      size=4
      do 20 j=1,nRocks
       ! set rock color rockColor
        iup=0
        do 20 m=1,mct
          x=tRock(j,m)
          y=dRock(j,m)/1000.
          call plot(XYPlot,x,y,iup)
          call symb(XYPlot,x,y,isymb,size)
          iup=1
20    CONTINUE

! plot symbols on rock P-T points, if desired

      if(isyplot.eq.1)then
        do 30 j=1,nRocks
        do 30 m=1,mct
        x=tRock(j,m)
        y=dRock(j,m)/1000.
        call symb(XYPlot,x,y,j+4,4)
30      CONTINUE
      endif

      return
      end
! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE RockOut(iout,modelTime)

! print data to terminal and disk

      implicit none
	Include "ThrustnHeat.inc"
! **********************************************
      Character*25 words
      real*8 modelTime
      integer j,m,iout
! *********************************************************
      data words/' Rock paths position =   '/

!     OUTPUT data to unit iout
1     continue
      write(iout,1044)modelinFull
1044  format(a40,'       , = Input file name') 
!      write(iout,*)'   '
      WRITE(iout,*)modelTime, '      = Total time of model (Ma)'
      WRITE(iout,*)nRocks   , '      = Number of rocks '
      WRITE(iout,*)mct,       '      = Number of iterations in this model (mct)'

!     Output Rock P-T paths to disk
      do 11 j=1,nRocks
!      write(iout,*)
      write(iout,*)'**************************************************'
      write(iout,*)'  Rock #  ',j
      write(iout,*)rockColor(j),'      Rock color'
      write(iout,*)'  Point #    Time   Depth(m)        P(kb)      T'
      do 10 m=1,mct
      write(iout,3506)m,rockTime(m),dRock(j,m),dRock(j,m)*density*9.8/(1.e8), tRock(j,m)
3506  format(i5,f12.5,f12.2,f12.3,f12.1)
10    continue
11    continue

      return
      end


! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine WritePositionFile(modelSteps,saveTimeMa)

!     Routine to write original positions of each grid point to make it easier to locate rocks
      implicit none
	Include "ThrustnHeat.inc"
	real*8 saveTimeMa
	integer modelSteps

! **********************************************************

      integer*4 i,k
! ************************************

      positionFile=trim(modelinFull)//".PosFile.txt"    
!      open(23,file=tempFile,FORM='UNFORMATTED',status='unknown')
      open(28,file=positionFile,FORM='formatted',status='unknown')

	write(28,*)GridX,GridZ,nvert,gridbottom,exaggeration,saveTimeMa,modelSteps,tempNum
      do 10 i = 1,gridBottom       ! do a row
!         do 20 k=1,nvert      ! loop through each column
		write(28,28)(xGrid_Old(k,i),dOld(k,i),k=1,nvert)
28		format(500(2F10.1))
!20       continue
10    continue
	do 12 i = 1,modelSteps
	write(28,*)i,singleRockIndex(i)
12	continue


      close(28)

      return
      end
! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine WriteFaultFile(modelStep)

!     Routine to write out location of faults at modelstep
      implicit none
	Include "ThrustnHeat.inc"
	real*8 saveTimeMa
	integer modelStep
	character*255 faultfile
! **********************************************************

      integer*4 i,k
! ************************************
	write(ppp,100)modelstep
100	format(I5)
      	FaultFile=trim(modelinFull)//".FaultFile_step_"//trim(ppp)//".txt"    
!      	open(23,file=tempFile,FORM='UNFORMATTED',status='unknown')
      	open(28,file=FaultFile,FORM='formatted',status='unknown')
	write(28,101)(k,k=1,nvert)
101	format(500I6)
	write(28,101)(surface(k),k=1,nvert)
	write(28,101)(surface(k)-sealevel,k=1,nvert)
	write(28,*)' '
	do 10 i = 1,nThrust
	write(28,*)' '
	write(28,*)'Fault = ',i
	write(28,101)(thrustGridZ(i,k),k=1,nvert)
	write(28,101)(thrustGridZ(i,k)-sealevel,k=1,nvert)
	write(28,101)(int((thrustGridZ(i,k)-sealevel)*gridZ),k=1,nvert)

10	continue


      close(28)

      return
      end

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine Save2D()

!     Routine to save 2 dimensional grid to a format that can be imported
!     into NIH Image

      implicit none
     
	Include "ThrustnHeat.inc"

! **********************************************************

!     Imageout is the array to make to write out.
!     The first dimension must be larger than the tC array because we write
!     the grid out with no vertical exaggeration.
!     tC is dimensioned tC(50,400).
!     To accomodate a 10x vertical/horizontal ratio, imageout needs to be
!     dimensioned 500,400 (e.g. 400k for this one array!)
      integer*2 imageout(1000,1000)
      common /imagine/imageout

      integer*4 i,j,jj,jjj,k,itemp,rkcolor
      character*32 chnum
! ************************************

!     option to save 2-D grid to binary file for opening in IMAGE program

!     open the file.  Give it an extension .number, where number is the number of file in this problem
      tempNum=tempNum+1
      write(chNum,'(i4)') tempNum
      tempFile=trim(modelinFull)//"."//trim(chNum)    
!      open(23,file=tempFile,FORM='UNFORMATTED',status='unknown')
      open(23,file=tempFile,FORM='UNFORMATTED',status='unknown',ACCESS="Transparent")

!     First make the integer*2 array we will write out
!     the "surface" is at grid point 50 

!     We must fill inbetween the vertical sections because the vertical grid spacing is 
!           not the same as the same as the horizontal grid spacint
!     exaggeration is the ratio of horizontal to vertical grid spacing
!      exaggeration = int((gridX/gridZ)+.01)
!      pixelNumber = nvert*exaggeration

      do 10 i = 1,gridBottom       ! do a row
         j = 0    ! pixel counter
         do 20 k=1,nvert      ! loop through each column
           itemp = int(tC(k,i))
           do 30 jj = 1,exaggeration     
           j = j + 1
           imageout(j,i) = itemp
30         continue
20       continue
10    continue


! ----------------------------------------------------------------------
!     Plot thrust positions
!     thrustGridZ(4,50) is the z position of each thrust at each vertical section

      do 510 jj = 1,nThrust       ! do a thrust
         j = 0    ! pixel counter
         do 520 k=1,nvert      ! loop through each column
           i = thrustGridZ(jj,k)
           do 530 jjj = 1,exaggeration     
           j = j + 1
!           imageout(j,i) = 0          ! this will force black in NIH Image
	    if(i.gt.1000)then
	    	write(*,*)' i = ',i
	    	pause ' this will bomb code'
	    	endif
 	    if(j.gt.1000)then
	    	write(*,*)' i = ',i
	    	pause ' this will bomb code'
	    	endif
          imageout(j,i) = 1000          ! this will force special color #1 in NIH Image
530         continue
520       continue
510    continue

! ----------------------------------------------------------------------
!     Plot rock positions
	rkcolor = 900
      do 400 k=1,nRocks

!     rockGridZ is the z grid point of this rock
      i = rockGridZ(k,mct)

!     get the horizontal grid point of this rock
      j = rockGridX(k,mct)
      j = (j-1)*exaggeration+1
!      imageout(j,i)=10000         		! This will force a "special color" in NIH image
      if(j.le.0) goto 400             	! do not plot the rock
      imageout(j,i)= rkcolor        		! This will force special color #2 in NIH image
	
	! Make the rock position 3x3
	imageout(j-1,i-1) = rkcolor
	imageout(j,i-1) = rkcolor
	imageout(j+1,i-1) = rkcolor
	imageout(j-1,i) = rkcolor
	imageout(j+1,i) = rkcolor
	imageout(j-1,i+1) = rkcolor
	imageout(j,i+1) = rkcolor
	imageout(j+1,i+1) = rkcolor


400   CONTINUE

! ----------------------------------------------------------------------
!     now write the temperature to a file
      do 120 i = 1,gridBottom
      write(23)(imageout(j,i),j=1,pixelNumber)
120   continue

      close(23)

      singleRockIndex(mct) = tempNum
      
      return
      end

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine RxExp1


      implicit none
	Include "ThrustnHeat.inc"


      real*8 curvex,curvez,TheCurve,DeltaT,TheKappa,OldT,theKond,TK
      real*8 gridz2,gridx2
      integer k,j

!     routine to do thermal relaxation over entire grid

!     This code is for a fully explicit 2-D grid 

!     Calculate the curvature at each grid point and store in array Curve(i,j)
!     The curvature is calculated in both x and y directions and summed:

!     curve(i,j) = 2T/x2 + 2T/y2

!     for most grid points the curvature is calculated as

!      2T/X2 = (T(i+1) - 2T(i) + T(i-1))/x^2

!      2T/Y2 = (T(i+1) - 2T(i) + T(i-1))/y^2

!     for the near-surface grid points we need a slightly different
!     formulation, but that is not yet implemented



!     Loop for each v. section

!      write(*,*)' Calculating curvature'
      gridz2 = gridZ*gridZ
      gridx2 = gridX*gridX
!     the edges are fixed temperature
      do 10 k = 2, nvert-1
      do 11 j = surface(k)+1,gridBottom-1
      curvez = (tC(k,j-1) - 2.0*tC(k,j) + tC(k,j+1))/(gridz2)
      curvex = (tC(k-1,j) - 2.0*tC(k,j) + tC(k+1,j))/(gridx2)
      curve(k,j) = curvex+curvez
11    continue
10    continue

!     do the left most vertical section
!     Note that here we assume that the "off grid" temperature is the
!     same as the first column
      k = 1
      do 15 j = surface(k)+1,gridBottom - 1
      curvez = (tC(k,j-1) - 2.0*tC(k,j) + tC(k,j+1))/(gridz2)
      curvex = (tC(k,j) - 2.0*tC(k,j) + tC(k+1,j))/(gridx2)
      curve(k,j) = curvex+curvez
15    continue

!     do the right most vertical sectin
!     Note that here we assume that the "off grid" temperature is the
!     same as the last (right-most) column
      k = nvert
      do 17 j = surface(k)+1,gridBottom - 1
      curvez = (tC(k,j-1) - 2.0*tC(k,j) + tC(k,j+1))/(gridz2)
      curvex = (tC(k-1,j) - 2.0*tC(k,j) + tC(k,j))/(gridx2)
      curve(k,j) = curvex+curvez
17    continue



!      Now solve for the new time step
      do 20 k = 1, nvert
      do 25 j = surface(k)+1,gridBottom-1
      TheCurve = curve(k,j)
      if(kappaFofT.eq.0)then
	      TheKappa = kappa(k,j)
		else
		TK = tC(k,j)+273.15
		TheKappa =  secMa*(567.3/TK - 0.062)/1.0e6
		endif
      DeltaT = TheCurve*TheKappa*dTimeMa
      OldT = tC(k,j)
      tC(k,j) = OldT + DeltaT
25    continue
      if(kappaFofT.eq.0)then
	      theKond = kond(k,gridbottom)
		else
		TK = tC(k,gridbottom-1)+273.15
		theKappa =  (567.3/TK - 0.062)/1.0e6			!convert from mm^2 to m^2 sec-1
		theKond = theKappa*rho(k,gridbottom)*cp(k,gridbottom)
		endif
      tC(k,gridBottom)=tC(k,gridBottom-1) + (qStar(k)/theKond)*gridZ

20    continue

      return
      end


! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine MoveThrusts (i,isostasy_flag,itopt)

	implicit none
	! i = the model step number
	Include 'ThrustnHeat.inc'

! ---------------------------------------
!      local variables     
      integer*4 k,i,j,it
      real*8 roundoff
      integer*4 itopt,isostasy_flag



!     First move the thrusts
!	Are we adding a thrust?
	if(addAFault(i).eq.1)then
!     set the grid points for the thrust to grid points
!     loop over all thrusts
!       Loop through nvert - 1 vertical sections
        do 41 k = 1,nvert
!          find the nearest grid point to the depth of the thrust at this vertical section
!          Assign this grid pointto thrustGridZ, which is regarded as the appoximation to the thrust line
	     if(dThrust(nThrust,k).ge.0.)then
			roundoff = 0.01
			else
			roundoff = -0.01
			endif
!          compute the depth of grid point thrustGridZ(i,k)
	if(FaultReference(nThrust).eq.'sealeve')then		! only 7 characters are read and checked
 		thrustGridZ(nThrust,k)=INT(dThrust(nThrust,k)/gridZ + roundoff) + sealevel
		else	! we are referencing to the 'surface'
		thrustGridZ(nThrust,k)=INT(dThrust(nThrust,k)/gridZ + roundoff) + surface(k)
		endif
 		! I removed the "+surface(k)" term so that thrusts are now entered relative to sealevel, rather than the surface
 		!  Jan 17, 2012
!          compute the depth of grid point thrustGridZ(i,k)
!           dThrust(nThrust,k)=(thrustGridZ(nThrust,k) - surface(k))*gridZ
	   ! note that dThrust is never used again in the program
           dThrust(nThrust,k) =(thrustGridZ(nThrust,k) - surface(k))*gridZ
41      continue

	endif		! done adding a fault
	
	!now plot the fault on the grid


      	do 5020 j=1,nThrust
	if(thrustdirection(i,j).eq.0)go to 5020		! we aren't moving this thrust
      	if(thrustTopBot(i,j).eq."top")then		!moving the top of the section
	      	if(thrustDirection(i,j).ge.1)then	! moving to the left
	      		do 5040 it = 1,abs(thrustdirection(i,j))
      	      		call ThrustLTop(j)
5040			continue  
      	      		else				! moving to the right
	      		do 5041 it = 1,abs(thrustdirection(i,j))
      	      		call ThrustRTop(j)
5041			continue  
      	      		endif
		else					! moving the bottom of the section "bot"
	      	if(thrustDirection(i,j).ge.1)then	! moving to the left
	      		do 5042 it = 1,abs(thrustdirection(i,j))
      	      		call ThrustLBot(j)
5042			continue  
      	      		else				! moving to the right
	      		do 5043 it = 1,abs(thrustdirection(i,j))
      	      		call ThrustRBot(j)
5043			continue  
      	      		endif
		endif	
5020  continue


!     now insert a pluton, if any
      if(plutonSwitch(i).eq.1)then
!      write(*,*)' Inserting pluton'
         do 6000 k = plutonUpperX(i),plutonLowerX(i)
         do 6000 j = plutonUpperZ(i),plutonLowerZ(i)
            tC(k,j) = plutonT(i)
            cp(k,j) = plutonCp(i)
6000     continue
         endif

!     Erode the surface
	Call Erosion(i)

!	if(erosionSwitch(i).eq.1)then
!		do 5030 k = 1,nvert
!		if(erode(i,k).eq.0)go to 5030		!no erosion
!		if(erode(i,k).gt.0)then	!erosion is positive - remove material
!			surface(k) = surface(k) + erode(i,k)
!			tC(k,surface(k)) = tSurface
!			do 5031 j = 1,surface(k)-1           ! set temperature of eroded cells to sky T
!			tC(k,j) = tSky
!5031	  		continue
!			else				! erosion is negative - add material to top
!			!note surface is still old surface for this loop
!			do 5033 j=surface(k)+erode(i,k),surface(k)
!			tC(k,j) = tSurface
!			rho(k,j)=rhoLay(k,1)
!			cp(k,j)=cpLay(k,1)
!			kond(k,j)=kLay(k,1)
!			kappa(k,j)=kappaLay(k,1)
!			agen(k,j)=agenLay(k,1)
!5033			continue
!			surface(k) = surface(k) + erode(i,k)
!			endif
!5030		continue	! end erosion loop
!		endif


!     now that all thrusts are moved, and surface is eroded
!      compute depth to each rock
      	do 5025 j=1,nRocks
        k = rockGridX(j,mct)        ! vertical section for this rock
        dRock(j,mct) = (float(rockGridZ(j,mct) - surface(k)))*gridZ
5025  	continue

	! if we are invoking isostasy then call the routine here

	if (isostasy_flag.eq.1)then
		call isostasy(itopt)
		endif
		
	return
	end



! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE ThrustRTop(ith)

!  variable   type    i/O      description
! -----------------------------------------------------------
!  mct        int      i       thrust "event" counter=number of times we have moved
!  ith        int      i       the ith thrust
! -----------------------------------------------------------
!     SUBROUTINE TO move thrust number ith
!     This routine moves upper plate to the RIGHT

      implicit none
	Include "ThrustnHeat.inc"
! ************************************************************
      integer i,ii,k,kk,j,deltaGrid,ith
      real*8 dtogo,DTC
! ************************************************************

!     Move the kth vertical section to the k-1 position
!     Only rocks above the fault move
!     Note that we are assuming that rocks that plot on the fault are part of
!     the lower plate and don't move

!     The RIGHT most part of the grid above the fault is removed from the grid
!     the LEFT most part is extrapolated from the LEFT side of the grid

!     To fix the LEFTmost (1) vertical section, we
!     first store the information for the 1 vertical section in the nvert+1 column.
!     We must do this because we can't be sure whether the shift will be up or down
      j = nvert+1
      do 30 i = surface(1),gridBottom
        tC(j,i)   = tC(1,i)
        rho(j,i)  =rho(1,i)
        cp(j,i)   =cp(1,i)
        kond(j,i) =kond(1,i)
        agen(j,i) =agen(1,i)
        kappa(j,i)=kappa(1,i)
        dold(j,i) =dold(1,i)
        told(j,i) =told(1,i)
        xGrid_old(j,i) = xGrid_old(1,i)-1		! this should generate negative numbers in xGrid_old if it starts off the grid
30    continue
      surface(j) = surface(1)
      deltaGrid = thrustGridZ(ith,2) - thrustGridZ(ith,1)
!      thrustGridZ(ith,j) = thrustGridZ(ith,1)+deltaGrid
      do 15 j = 1,nThrust
      thrustGridZ(j,nvert+1) = thrustGridZ(j,1)+deltaGrid
15    continue


      do 10 k = nvert-1,0,-1    ! work from right to left
      if(k.eq.0)then
            kk=nvert+1          ! this picks up the temporary storage arrays from the right of the grid
            else
            kk = k    ! otherwise the kth column shifts from the left
            endif

!     compute the number of vertical grid points to shift
      deltaGrid = thrustGridZ(ith,kk) - thrustGridZ(ith,k+1)
!     (Note that for a thrust that dips down to the east, deltaGrid is negative.
!     Thus when it says "up" below, it means "up" by a negative amount i.e. "down"
!     Also note we are moving "from" v section k "to" v section k+1
!           (kk = k except when we are at left hand side (k=0))
!     All temperatures and rock properties above the thrust must shift up by this number of grid points

!     Surface must shift up by this number of grid points
      	surface(k+1) = surface(kk) - deltaGrid
	ElevationMeters(k+1) = float(sealevel-surface(k+1))*gridZ
      do 20 i = surface(k+1),thrustGridZ(ith,kk)-1     ! only do loop for rocks above thrust
        ii = i - deltaGrid
        tC(k+1,ii)    = tC(kk,i)
        rho(k+1,ii)   = rho(kk,i)
        cp(k+1,ii)    = cp(kk,i)
        kond(k+1,ii)  = kond(kk,i)
        agen(k+1,ii)  = agen(kk,i)
        kappa(k+1,ii) = kappa(kk,i)
        dold(k+1,ii)  = dold(kk,i)
        told(k+1,ii)  = told(kk,i)
        xGrid_old(k+1,ii) = xGrid_old(kk,i)
20    continue

! 	set the temperature of the surface....just in case a new surface is exposed
	tC(k+1,surface(k+1)) = tSurface

      if(surface(k+1).eq.0) then
      	write(13,*)k,k+1,surface(k+1)
        write(*,*) 'redefine the parameter seaLevel to a bigger number'
        write(*,*) 'the parameter is defined in the input file'
        write(*,*) 'RETURN to stop'
        pause
        stop
      endif


! ----------------------------------------
!     Process the rocks in this vertical section

!     Rocks above the thrust must be moved the same amount as any grid point
!     Check to see if there are any rocks in this vertical section that must
!     be moved
!     note that mct is incremented before calling this subroutine
!           therefore, mct is where we are RIGHT NOW
      do 150 j=1,nRocks
      if(rockGridX(j,mct).eq.kk)then        ! this one must be checked (it is on v section kk)
!     check to see if we are above or below fault
      if(rockGridZ(j,mct).lt.thrustGridZ(ith,kk)) then
        ! if here then this rock is above the thrust and must move accordingly
        rockGridZ(j,mct) = rockGridZ(j,mct) - deltaGrid
        rockGridX(j,mct) = rockGridX(j,mct) + 1
        else
        !  if here, then rock is below thrust and doesn't move
        endif
      endif

150   continue


	! move the single rock position
	if(singleRockX(mct).eq.kk)then	! our rock is on this vertical section
		if(singleRockZ(mct).lt.thrustGridZ(ith,kk))then		!our rock is above the fault and must be moved 
			singleRockZ(mct) = singleRockZ(mct) - deltaGrid
			singleRockX(mct) = singleRockX(mct) +1
			endif
		endif

! -----------------------------------------
!     Move thrusts that lie above this thrust

!     Faults above the thrust must be moved the same amount as any grid point
!     Check to see if there are any faults in this vertical section that must be moved
!     note that mct is incremented before calling this subroutine
!        therefore, mct is where we are RIGHT NOW
      do 160 j=1,nThrust
!     check to see if we are above or below moving fault
      if(thrustGridZ(j,kk).lt.thrustGridZ(ith,kk)) then
        ! if here then this fault is above the thrust and must move accordingly
        thrustGridZ(j,k+1) = thrustGridZ(j,kk) - deltaGrid
!        if(thrustGridZ(j,k+1).lt.surface(k+1))then
!            thrustGridZ(j,k+1)=surface(k+1)
!            endif
        else
        ! if here, then fault is below thrust and doesn't move
        endif
160   continue


10    continue
! 	end loop on every vertical section
! ---------------------------------------
	if(singleRockX(mct).le.0)then		!our rock is off the grid to the left
		if(singleRockZ(mct).lt.thrustGridZ(ith,1))then
			singleRockX(mct) = singleRockX(mct) +1
			endif
		endif		


!     Zero temperatures above the surface
	do 50 k = 1,nvert
      do 50 i = 1,surface(k)-1
      tC(k,i) = tSky
50    continue




! -----------------------------------------
! Calculate the temperature increase from shear heating
! Note that the thrust is ON the thrustGridZ(k) now 
! in the Kth vertical section
! Note that qShear is the heat produced for one time increment 
! and has units of j m-2

	call ShearHeat(ith)
!      if(Shear(ith).eq.0)go to 111
!      do 110 k = 2,nvert
!      deltaGrid = thrustGridZ(ith,k-1) - thrustGridZ(ith,k)
!     dtogo is the distance traveled in 1 thrust step
!      dtogo = sqrt((deltaGrid)**2 + gridX**2)
!     qShear is the heat produced in 1 thrust step    
!      qShear(ith)= dtogo * shear(ith)

!      DTC = qShear(ith)/(rho(k,thrustGridZ(ith,k))*cp(k,thrustGridZ(ith,k))*2.*gridZ)
!      tC(k,thrustGridZ(ith,k)) = tC(k,thrustGridZ(ith,k)) + DTC
!110   continue
!111   continue
!	write(*,*)k,dtogo,qShear(ith),DTC



      RETURN
      END




! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE ThrustLTop(ith)

!  variable   type    i/O      description
! -----------------------------------------------------------
!  mct        int      i       thrust "event" counter=number of times we have moved
!  ith        int      i       the ith thrust
!-----------------------------------------------------------
!     SUBROUTINE TO move thrust number ith
!	move the top to the left

      implicit none
	Include "ThrustnHeat.inc"
! ************************************************************
      integer i,ii,k,j,deltaGrid,ith
      real*8 dtogo,DTC
! ************************************************************


!     Move the kth vertical section to the k-1 position
!     Only rocks above the fault move
!     Note that we are assuming that rocks that plot on the fault are part of
!     the lower plate and don't move

!     The left most part of the grid above the fault is removed from the grid
!     the right most part is extrapolated from the right side of the grid

!     To fix the rightmost (nvert) vertical section, we
!     first store the information for the nvert vertical section in the nvert+1 column.
!     We must do this because we can't be sure whether the shift will be up or down
      j = nvert+1
      do 30 i = surface(nvert),gridBottom
        tC(j,i)   = tC(nvert,i)
        rho(j,i)  =rho(nvert,i)
        cp(j,i)   =cp(nvert,i)
        kond(j,i) =kond(nvert,i)
        agen(j,i) =agen(nvert,i)
        kappa(j,i)=kappa(nvert,i)
        dold(j,i) =dold(nvert,i)
        told(j,i) =told(nvert,i)
        xGrid_old(j,i) = xGrid_old(nvert,i) + 1		!this should generate off grid numbers for the top on the right
30    continue
      surface(j) = surface(nvert)
      deltaGrid = thrustGridZ(ith,nvert) - thrustGridZ(ith,nvert-1)

      do 15 j = 1,nThrust
      thrustGridZ(j,nvert+1) = thrustGridZ(j,nvert)+deltaGrid
15    continue

      do 10 k = 2,nvert+1

	! check to see if the fault is below the surface
	if(thrustGridZ(ith,k).gt.surface(k)) then
!     	compute the number of grid points to shift
      	deltaGrid = thrustGridZ(ith,k) - thrustGridZ(ith,k-1)
!     	All temperatures and rock properties above the thrust must shift up by this number of grid points
!     	Surface must shift up by this number of grid points
      	surface(k-1) = surface(k) - deltaGrid
	ElevationMeters(k-1) = float(sealevel-surface(k-1))*gridZ
      	do 20 i = surface(k-1),thrustGridZ(ith,k)-1     ! only do loop for rocks above thrust
      	  ii = i - deltaGrid
      	  tC(k-1,ii)    = tC(k,i)
      	  rho(k-1,ii)   = rho(k,i)
      	  cp(k-1,ii)    = cp(k,i)
      	  kond(k-1,ii)  = kond(k,i)
      	  agen(k-1,ii)  = agen(k,i)
      	  kappa(k-1,ii) = kappa(k,i)
      	  dold(k-1,ii)  = dold(k,i)
      	  told(k-1,ii)  = told(k,i)
	        xGrid_old(k-1,ii) = xGrid_old(k,i)
20    	continue
		endif


! 	set the temperature of the surface....just in case a new surface is exposed
	tC(k-1,surface(k-1)) = tSurface

      if(surface(k-1).eq.0) then
        write(*,*) 'redefine the parameter NBEGIN to a biger number'
        write(*,*) 'the parameter is defined in thrust.inc'
        write(*,*) 'RETURN to stop'
        pause
        stop
      endif


! ----------------------------------------
!     Process the rocks in this vertical section

!     Rocks above the thrust must be moved the same amount as any grid point
!     Check to see if there are any rocks in this vertical section that must
!     be moved
!     note that mct is incremented before calling this subroutine
!        therefore, mct is where we are RIGHT NOW
      do 150 j=1,nRocks
      if(rockGridX(j,mct).eq.k)then        ! this one must be moved
!     check to see if we are above or below fault
      if(rockGridZ(j,mct).lt.thrustGridZ(ith,k)) then
        ! if here then this rock is above the thrust and must move accordingly
        rockGridZ(j,mct) = rockGridZ(j,mct) - deltaGrid
        rockGridX(j,mct) = rockGridX(j,mct) - 1
        else
        ! if here, then rock is below thrust and doesn't move
        endif
      ! in either case, the temperatures must stay the same as last step   
      endif
150   continue


	! move the single rock position
	if(singleRockX(mct).eq.k)then	! our rock is on this vertical section
		if(singleRockZ(mct).lt.thrustGridZ(ith,k))then		!our rock is above the fault and must be moved 
			singleRockZ(mct) = singleRockZ(mct) - deltaGrid
			singleRockX(mct) = singleRockX(mct) - 1
			endif
		endif

! -----------------------------------------
!     Move thrusts that lie above this thrust

!     Faults above the thrust must be moved the same amount as any grid point
!     Check to see if there are any faults in this vertical section that must be moved
!     note that mct is incremented before calling this subroutine
!        therefore, mct is where we are RIGHT NOW
      do 160 j=1,nThrust
!     check to see if we are above or below moving fault
      if(thrustGridZ(j,k).lt.thrustGridZ(ith,k)) then
        ! if here then this fault is above the thrust and must move accordingly
        thrustGridZ(j,k-1) = thrustGridZ(j,k) - deltaGrid
        if(thrustGridZ(j,k-1).lt.surface(k-1))then
            thrustGridZ(j,k-1)=surface(k-1)
            endif
        else
        ! if here, then fault is below thrust and doesn't move
        endif
160   continue


10    continue

	if(singleRockX(mct).gt.nvert)then		!our rock is off the grid to the right
		if(singleRockZ(mct).lt.thrustGridZ(ith,nvert))then
			singleRockX(mct) = singleRockX(mct) - 1   ! move it to the left
			endif
		endif		

!     Zero temperatures above the surface
	do 50 k = 1,nvert
      do 50 i = 1,surface(k)-1
      tC(k,i) = tSky
50    continue


! ----------------------------------------
! Calculate the temperature increase from shear heating
! Note that the thrust is ON the thrustGridZ(k) now 
! in the Kth vertical section
! Note that qShear is the heat produced for one time increment 
! and has units of j m-2

	call ShearHeat(ith)

!      if(Shear(ith).eq.0)go to 111
!      do 110 k = 2,nvert
!      deltaGrid = thrustGridZ(ith,k) - thrustGridZ(ith,k-1)
!     dtogo is the distance traveled in 1 thrust step
!      dtogo = sqrt((deltaGrid)**2 + gridX**2)

!     qShear is the heat produced in 1 thrust step    
!      qShear(ith)= dtogo * shear(ith)

!      DTC = qShear(ith)/(rho(k,thrustGridZ(ith,k))*cp(k,thrustGridZ(ith,k))*2.*gridZ)
!      tC(k,thrustGridZ(ith,k)) = tC(k,thrustGridZ(ith,k)) + DTC
!110   continue
!111	continue
!	write(*,*)k,dtogo,qShear(ith),DTC


      RETURN
      END


! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE ThrustRBot(ith)

!  variable   type    i/O      description
! -----------------------------------------------------------
!  mct        int      i       thrust "event" counter=number of times we have moved
!  ith        int      i       the ith thrust
! -----------------------------------------------------------
!     SUBROUTINE TO move thrust number ith
!     This routine moves lower plate to the RIGHT (ie the Bottom plate)

      implicit none
	Include "ThrustnHeat.inc"
! ************************************************************
      integer i,ii,k,kk,j,deltaGrid,ith
      real*8 dtogo,DTC
! ************************************************************

!     Move the kth vertical section to the k+1 position
!     Only rocks below the fault move
!     Note that we are assuming that rocks that plot on the fault are part of
!     the lower plate and therefore move (true for rocks but not for other faults???)

!     The RIGHT most part of the grid below the fault is removed from the grid
!     the LEFT most part is extrapolated from the LEFT side of the grid

!     To fix the LEFTmost (1) vertical section, we
!     first store the information for the 1 vertical section in the nvert+1 column.
!     We must do this because we can't be sure whether the shift will be up or down
      j = nvert+1
      do 30 i = surface(1),gridBottom
        tC(j,i)   = tC(1,i)
        rho(j,i)  =rho(1,i)
        cp(j,i)   =cp(1,i)
        kond(j,i) =kond(1,i)
        agen(j,i) =agen(1,i)
        kappa(j,i)=kappa(1,i)
        dold(j,i) =dold(1,i)
        told(j,i) =told(1,i)
        xGrid_old(j,i) = xGrid_old(1,i) - 1		!this should generate negative numbers for the bottom
30    continue
      surface(j) = surface(1)
      deltaGrid = thrustGridZ(ith,2) - thrustGridZ(ith,1)
!      thrustGridZ(ith,j) = thrustGridZ(ith,1)+deltaGrid
      do 15 j = 1,nThrust
      thrustGridZ(j,nvert+1) = thrustGridZ(j,1)+deltaGrid
15    continue


      do 10 k = nvert-1,0,-1    ! work from right to left
      if(k.eq.0)then
            kk=nvert+1          ! this picks up the temporary storage arrays from the right of the grid
            else
            kk = k    ! otherwise the kth column shifts from the left
            endif

!     compute the number of vertical grid points to shift
      deltaGrid = thrustGridZ(ith,kk) - thrustGridZ(ith,k+1)
!     (Note that for a thrust that dips down to the east, deltaGrid is negative.
!     Thus when it says "up" below, it means "up" by a negative amount i.e. "down"
!     Also note we are moving "from" v section k "to" v section k+1
!           (kk = k except when we are at left hand side (k=0))
!     All temperatures and rock properties below the thrust must shift up by this number of grid points

!     Surface must shift up by this number of grid points
!      surface(k+1) = surface(kk) - deltaGrid
!	ElevationMeters(k+1) = float((sealevel-surface(k+1))*gridZ
!     Surface shouldn't shift if we are moving the bottom down.

      do 20 i = thrustGridZ(ith,kk),gridBottom     ! only do loop for rocks from the thrust to the bottom of the grid
        ii = i - deltaGrid
        tC(k+1,ii)    = tC(kk,i)
        rho(k+1,ii)   = rho(kk,i)
        cp(k+1,ii)    = cp(kk,i)
        kond(k+1,ii)  = kond(kk,i)
        agen(k+1,ii)  = agen(kk,i)
        kappa(k+1,ii) = kappa(kk,i)
        dold(k+1,ii)  = dold(kk,i)
        told(k+1,ii)  = told(kk,i)
        xGrid_old(k+1,ii) = xGrid_old(kk,i)
20    continue

! 	set the temperature of the surface....just in case a new surface is exposed
!	tC(k+1,surface(k+1)) = tSurface

!      if(surface(k+1).eq.0) then
!        write(*,*) 'redefine the parameter seaLevel to a bigger number'
!        write(*,*) 'the parameter is defined in the input file'
!        write(*,*) 'RETURN to stop'
!        pause
!        stop
!      endif


! ----------------------------------------
!     Process the rocks in this vertical section

!     Rocks below the thrust must be moved the same amount as any grid point
!     Check to see if there are any rocks in this vertical section that must
!     be moved
!     note that mct is incremented before calling this subroutine
!           therefore, mct is where we are RIGHT NOW
      do 150 j=1,nRocks
      if(rockGridX(j,mct).eq.kk)then        ! this one must be checked (it is on v section kk)
!     check to see if we are above or below fault
      if(rockGridZ(j,mct).ge.thrustGridZ(ith,kk)) then
        ! if here then this rock is below the thrust and must move accordingly
        rockGridZ(j,mct) = rockGridZ(j,mct) - deltaGrid
        rockGridX(j,mct) = rockGridX(j,mct) + 1
        else
        !  if here, then rock is above thrust and doesn't move
        endif
      endif
150   continue

	! move the single rock position
	if(singleRockX(mct).eq.kk)then	! our rock is on this vertical section
		if(singleRockZ(mct).ge.thrustGridZ(ith,kk))then		!our rock is above the fault and must be moved 
			singleRockZ(mct) = singleRockZ(mct) - deltaGrid
			singleRockX(mct) = singleRockX(mct) + 1
			endif
		endif



! -----------------------------------------
!     Move thrusts that lie below this thrust

!     Faults below the thrust must be moved the same amount as any grid point
!     Check to see if there are any faults in this vertical section that must be moved
!     note that mct is incremented before calling this subroutine
!        therefore, mct is where we are RIGHT NOW
      do 160 j=1,nThrust
      if(j.eq.ith)go to 160	! don't try to move the one we're working on
!     check to see if we are above or below moving fault
!      if(thrustGridZ(j,kk).lt.thrustGridZ(ith,kk)) then
!      if(thrustGridZ(j,kk).ge.thrustGridZ(ith,kk)) then	! changed to gt so we won't move a thrust that is in the identical depth
      if(thrustGridZ(j,kk).gt.thrustGridZ(ith,kk)) then
        ! if here then this fault is below the thrust and must move accordingly
! write(*,*)k,kk,thrustGridZ(j,kk),thrustGridZ(ith,kk)
! pause 'hit return'
        thrustGridZ(j,k+1) = thrustGridZ(j,kk) - deltaGrid
!        thrustGridZ(j,kk) = thrustGridZ(j,kk) - deltaGrid
!        if(thrustGridZ(j,k+1).lt.surface(k+1))then
!            thrustGridZ(j,k+1)=surface(k+1)
!            endif
        else
        ! if here, then fault is above thrust and doesn't move
        endif
160   continue


10    continue
! 	end loop on every vertical section
! ---------------------------------------

	if(singleRockX(mct).le.0)then		!our rock is off the grid to the left
		if(singleRockZ(mct).ge.thrustGridZ(ith,1))then
			singleRockX(mct) = singleRockX(mct) + 1
			endif
		endif		


!     Zero temperatures above the surface
!	do 50 k = 1,nvert
!      do 50 i = 1,surface(k)-1
!      tC(k,i) = tSky
!50    continue




! -----------------------------------------
! Calculate the temperature increase from shear heating
! Note that the thrust is ON the thrustGridZ(k) now 
! in the Kth vertical section
! Note that qShear is the heat produced for one time increment 
! and has units of j m-2
	call ShearHeat(ith)

!      if(Shear(ith).eq.0)go to 111
!      do 110 k = 2,nvert
!      deltaGrid = thrustGridZ(ith,k-1) - thrustGridZ(ith,k)
!     dtogo is the distance traveled in 1 thrust step
!      dtogo = sqrt((deltaGrid)**2 + gridX**2)
!     qShear is the heat produced in 1 thrust step    
!      qShear(ith)= dtogo * shear(ith)

!      DTC = qShear(ith)/(rho(k,thrustGridZ(ith,k))*cp(k,thrustGridZ(ith,k))*2.*gridZ)
!      tC(k,thrustGridZ(ith,k)) = tC(k,thrustGridZ(ith,k)) + DTC
!110   continue
!111	continue
!	write(*,*)k,dtogo,qShear(ith),DTC


      RETURN
      END




! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE ThrustLBot(ith)

!  variable   type    i/O      description
! -----------------------------------------------------------
!  mct        int      i       thrust "event" counter=number of times we have moved
!  ith        int      i       the ith thrust
!-----------------------------------------------------------
!     SUBROUTINE TO move thrust number ith
!	move the lower plate to the left

      implicit none
	Include "ThrustnHeat.inc"
! ************************************************************
      integer i,ii,k,j,deltaGrid,ith
      real*8 dtogo,DTC
! ************************************************************


!     Move the kth vertical section to the k-1 position
!     Only rocks BELOW the fault move
!     Note that we are assuming that rocks that plot on the fault are part of
!     the lower plate and therefore do move

!     The left most part of the grid BELOW the fault is removed from the grid
!     the right most part is extrapolated from the right side of the grid

!     To fix the rightmost (nvert) vertical section, we
!     first store the information for the nvert vertical section in the nvert+1 column.
!     We must do this because we can't be sure whether the shift will be up or down
      j = nvert+1
      do 30 i = surface(nvert),gridBottom
        tC(j,i)   = tC(nvert,i)
        rho(j,i)  =rho(nvert,i)
        cp(j,i)   =cp(nvert,i)
        kond(j,i) =kond(nvert,i)
        agen(j,i) =agen(nvert,i)
        kappa(j,i)=kappa(nvert,i)
        dold(j,i) =dold(nvert,i)
        told(j,i) =told(nvert,i)
        xGrid_old(j,i) = xGrid_old(nvert,i) + 1		!should generate off grid numbers for the bottom right
30    continue
      surface(j) = surface(nvert)
      deltaGrid = thrustGridZ(ith,nvert) - thrustGridZ(ith,nvert-1)

      do 15 j = 1,nThrust
      thrustGridZ(j,nvert+1) = thrustGridZ(j,nvert)+deltaGrid
15    continue

      do 10 k = 2,nvert+1
!     compute the number of grid points to shift
      deltaGrid = thrustGridZ(ith,k) - thrustGridZ(ith,k-1)
!     All temperatures and rock properties BELOW the thrust must shift down by this number of grid points
!     Surface must shift up by this number of grid points
!      surface(k-1) = surface(k) - deltaGrid
!	ElevationMeters(k-1) = float((sealevel-surface(k-1))*gridZ
!	Surface DOES NOT shift
      do 20 i = thrustGridZ(ith,k),gridBottom     ! only do loop for rocks BELOW thrust
        ii = i - deltaGrid
        tC(k-1,ii)    = tC(k,i)
        rho(k-1,ii)   = rho(k,i)
        cp(k-1,ii)    = cp(k,i)
        kond(k-1,ii)  = kond(k,i)
        agen(k-1,ii)  = agen(k,i)
        kappa(k-1,ii) = kappa(k,i)
        dold(k-1,ii)  = dold(k,i)
        told(k-1,ii)  = told(k,i)
        xGrid_old(k-1,ii) = xGrid_old(k,i)
20    continue

! 	set the temperature of the surface....just in case a new surface is exposed
!	tC(k-1,surface(k-1)) = tSurface

!      if(surface(k-1).eq.0) then
!        write(*,*) 'redefine the parameter NBEGIN to a biger number'
!        write(*,*) 'the parameter is defined in thrust.inc'
!        write(*,*) 'RETURN to stop'
!        pause
!        stop
!      endif


! ----------------------------------------
!     Process the rocks in this vertical section

!     Rocks BELOW the thrust must be moved the same amount as any grid point
!     Check to see if there are any rocks in this vertical section that must
!     be moved
!     note that mct is incremented before calling this subroutine
!        therefore, mct is where we are RIGHT NOW
      do 150 j=1,nRocks
      if(rockGridX(j,mct).eq.k)then        ! this one must be moved
!     check to see if we are above or below fault
      if(rockGridZ(j,mct).ge.thrustGridZ(ith,k)) then
        ! if here then this rock is BELOW the thrust and must move accordingly
        rockGridZ(j,mct) = rockGridZ(j,mct) - deltaGrid
        rockGridX(j,mct) = rockGridX(j,mct) - 1
        else
        ! if here, then rock is above thrust and doesn't move
        endif
      ! in either case, the temperatures must stay the same as last step   
      endif
150   continue


	! move the single rock position
	if(singleRockX(mct).eq.k)then	! our rock is on this vertical section
		if(singleRockZ(mct).ge.thrustGridZ(ith,k))then		!our rock is above the fault and must be moved 
			singleRockZ(mct) = singleRockZ(mct) - deltaGrid
			singleRockX(mct) = singleRockX(mct) - 1
			endif
		endif

! -----------------------------------------
!     Move thrusts that lie above this thrust

!     Faults BELOW the thrust must be moved the same amount as any grid point
!     Check to see if there are any faults in this vertical section that must be moved
!     note that mct is incremented before calling this subroutine
!        therefore, mct is where we are RIGHT NOW
      do 160 j=1,nThrust
!     check to see if we are above or below moving fault
!      if(thrustGridZ(j,k).ge.thrustGridZ(ith,k)) then		! changed to gt so we don't move a thrust at the same depth
      if(thrustGridZ(j,k).gt.thrustGridZ(ith,k)) then		
        ! if here then this fault is BELOW the thrust and must move accordingly
        thrustGridZ(j,k-1) = thrustGridZ(j,k) - deltaGrid
        if(thrustGridZ(j,k-1).lt.surface(k-1))then
            thrustGridZ(j,k-1)=surface(k-1)
            endif
        else
        ! if here, then fault is above thrust and doesn't move
        endif
160   continue


10    continue

	if(singleRockX(mct).gt.nvert)then		!our rock is off the grid to the right
		if(singleRockZ(mct).ge.thrustGridZ(ith,nvert))then
			singleRockX(mct) = singleRockX(mct) - 1   ! move it to the left
			endif
		endif		

!     Zero temperatures above the surface
!	do 50 k = 1,nvert
!      do 50 i = 1,surface(k)-1
!      tC(k,i) = tSky
!50    continue


! ----------------------------------------
! Calculate the temperature increase from shear heating
! Note that the thrust is ON the thrustGridZ(k) now 
! in the Kth vertical section
! Note that qShear is the heat produced for one time increment 
! and has units of j m-2

	call ShearHeat(ith)
!      if(Shear(ith).eq.0)go to 111
!      do 110 k = 2,nvert
!      deltaGrid = thrustGridZ(ith,k) - thrustGridZ(ith,k-1)
!     dtogo is the distance traveled in 1 thrust step
!      dtogo = sqrt((deltaGrid)**2 + gridX**2)

!     qShear is the heat produced in 1 thrust step    
!      qShear(ith)= dtogo * shear(ith)

!      DTC = qShear(ith)/(rho(k,thrustGridZ(ith,k))*cp(k,thrustGridZ(ith,k))*2.*gridZ)
!      tC(k,thrustGridZ(ith,k)) = tC(k,thrustGridZ(ith,k)) + DTC
!110   continue
!111	continue
!	write(*,*)k,dtogo,qShear(ith),DTC


      RETURN
      END


! **********************************************************************
! **********************************************************************
! Calculate the temperature increase from shear heating
! Note that the thrust is ON the thrustGridZ(k) now 
! in the Kth vertical section
! Note that qShear is the heat produced for one time increment 
! 	and has units of j m-2
	Subroutine ShearHeat(ith)
	implicit none
	Include "ThrustnHeat.inc"
	real*8 dtogo,DTC
	integer*4 deltaGrid,ith,k

	if(shear(ith).eq.0.)return

	do 10 k = 2,nvert
	deltaGrid = thrustGridZ(ith,k) - thrustGridZ(ith,k-1)
!     dtogo is the distance traveled in 1 thrust step
	dtogo = sqrt((deltaGrid)**2 + gridX**2)

!     qShear is the heat produced in 1 thrust step   units of J m-2
	qShear(ith)= dtogo * shear(ith)

	DTC = qShear(ith)/(rho(k,thrustGridZ(ith,k))*cp(k,thrustGridZ(ith,k))*2.*gridZ)
	tC(k,thrustGridZ(ith,k)) = tC(k,thrustGridZ(ith,k)) + DTC
10   	continue
	write(*,*)k,dtogo,qShear(ith),DTC
	return
	end


! **********************************************************************
! **********************************************************************
      Subroutine Auto2D(XYPlot,inew)

!     Routine to plot 2-D thermal structure
!	inew = 0 the first time we call this routine (sets up a new canvas)
!	inew = 1 - clear canvas and rewrite with new information
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: XYPlot
	TYPE(AWE_CanvasPen) :: pen
	Type(AWE_CanvasBrush) :: brush

	Include "ThrustnHeat.inc"
	INCLUDE "PlotStuff.inc"				! in folder "AWE_Subs"
! *************************************************************

      integer i,j,k,num,ith,iup,inew
      real*4 x,y,dx,dy

! *********************************************************
!      real xor,xmin,xmax,xlen,yor,ymin,ymax,ylen,overlap
!      integer nxtic,nytic,nxdec,nydec,kk
!      character*40 xlab,ylab,pltitle

! Plot sections
! First set up scaling
! Total horizontal distance is nvert*gridX
! Vertical distance is set at 40 km
	select case (inew)
	Case (0)		! open a new canvas window

      xor=50
      xmin=0
!     note: Xdist is in meters
      xmax=gridX*(nvert-1)
!      xlen=18            ! 0 for autoscaling
      yor=200
!     Note:YDIST is in km
      ymin=-40
      ymax=0
!      ylen=14
      nxstep=0
      nystep=0
      nxdec=0
      nydec=0
      pltitle=' '
      xlab='  '
      ylab='  '

!      CALL user(xor,xmin,xmax,xlen,yor,ymin,ymax,ylen)

! plot the grid-------------------------
!	call clgs(2)
!	call fss_openpicture(2)
	xlen = 30.
	ylen= 20.
      	CALL USER()			! sets the user coordinates that were defined above
!	xlen and ylen are in cm. 
!	the canvas width and height are in pixels
!	to make sure the canvas is large enough we need to convert
!	The scaling in subroutine USER assumes 72 dpi so
!	cm = pixels * 2.54/72 or
!	pixels = cm * 72/2.54 = cm * 28.34646
!	add 25% just to be safe
	XYPlot%width = 1.25*(xlen * 28.34646)			! *2 is needed to make room for the assemblage information
	XYPlot%height = 1.5*(ylen * 28.34646)
!	canvas%backgroundColor = AWE_teal
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = XYPlot%height - 100.
	CALL AWE_createCanvas(XYPlot)	


	case(1)			!Clear canvas and redraw
	Call AWE_clearCanvas(XYPlot)
!	XYPlot%width = 1.25*(xlen * 28.34646)			! *2 is needed to make room for the assemblage information
!	XYPlot%height = 1.5*(ylen * 28.34646)

! plot surface, begin from the first point
      iup = 0
      do 101 k=1,nvert
        x = (k-1) * gridX
        y = float(seaLevel-surface(k))*gridZ/1000.
        call plot(XYPlot,x,y,iup)
        iup = 1
101    continue

! plot the horizontal line which is the original surface of the region

!      call plot(0.,0.,0)
!      x= gridX*(nvert-1)
!      call plot(x,0.,1)

! draw each vertical line. Each line may be of different height.

	go to 206		! skip this next part for now
      Do 205 k=1,nvert
      x = (k-1) * gridX
      call plot(XYPlot,x,0.,0)
      do 210 i = surface(k),gridBottom
!     find the depth of the point 
      y = float(i - seaLevel) * gridZ/1000.    ! plot in km, not meters
      y = -y
!     if the point is higher than the original surface, do not plot it
!      if(y.gt.0) then
!        goto 210
!        endif
!      call plot(x,y,1)
!	call symbplTh (x,y,2,1,1)
	pen%penStyle = CanvasPenStyle_SolidLine

!	skip plotting dots for now
!	call PlotCenteredEllipseOnScreen(XYPlot,x,dX,y,dY,brush,pen)
	
!      call symb(x,y,1,2)

!     draw tics on last line
!      if(k.eq.nvert)then
!        Draw tic mark
! 	   call plot(x,y,0)
!          call fss_line(8,0)
!       endif

210   continue
205   continue
206	continue

!     label y axis
      do 213 i=0,1000,5
      y = -gridZ*float(i)/1000.
      if(y.Lt.ymin)goto 217
      x=0
      WRITE(ppp,215)-y
215   format(F4.1)
	call TextOnPlot(XYPlot,x-400.,y+.7,ppp,11)
!     Draw tic mark
      call plot(XYPlot,x,y,0)
      call plot(XYPlot,x-50.,y,1)
!      call fss_line(5,0)    ! draw horizontal tic
!      call fss_move(-30,0)  ! move to left
!      num=-y
!      WRITE(ppp,215)num
!215   format(i4)
213   continue
217   continue


!     Draw and label first line
      x=0
      call plot(XYPlot,x,0.,0)
      do 212 i = surface(1),gridBottom
      y = float(i - seaLevel)*gridZ/1000.
      y = -y
      if(y.Lt.ymin)goto 214
      call plot(XYPlot,x,y,1)
!     Draw tic mark
!      call plot(XYPlot,x,y+5.,1)
!      call fss_line(0,5)
212   continue
214   continue

!     plot thrust fault

      do 312 ith = 1,nThrust
      x=0.0
      y = float(thrustGridZ(ith,1)-seaLevel)*gridZ/1000.
      y = -y
      call plot(XYPlot,x,y,0)
      do 303 k=1,nvert
        x=(k-1)*gridX
        y = float(thrustGridZ(ith,k)-seaLevel)*gridZ/1000.
        y = -y
!	write(*,*)x,y,xmin,xmax,ymin,ymax
        call plot(XYPlot,x,y,1)
303   continue
312   continue


! -----------------------------------------------------------------------
!     Plot rock positions

	!dx = 50.
	!dy = 0.2
	dx = (xmax-xmin)/400.
	dy = (ymax-ymin)/200.
	if(nRocks.eq.0)then
		currentcolor = 3
         	x = float(singleRockX(mct)-1)*gridX
         	if(x.lt.0.0) then
			X = 0		! plot the rock on the left
!            		goto 511    	! do not plot the rock
            		endif
         	if(x.gt.xmax) then
			X = xmax		! plot the rock on the right
!            		goto 511    ! do not plot the rock
            		endif
         	y = -float(singleRockZ(mct)-seaLevel)*gridZ/1000.
!         	y = -float(singleRockZ(mct))*gridZ/1000.
         	write(*,*)'PixelPlot',x,y
         	!call symb(XYPlot,x,y,6,4)
		Call PlotCenteredRectOnScreen(XYPlot,X,dx,Y,dy,brush,pen)
511   		continue
		else		
		currentcolor = 2
		do 420 j=1,nRocks
      		x = float(rockGridX(j,mct)-1)*gridX
         	if(x.lt.0.0) then
            		goto 510    ! do not plot the rock
            		endif
         	y = -float(rockGridZ(j,mct)-seaLevel)*gridZ/1000.
         	!call symb(XYPlot,x,y,6,4)
		Call PlotCenteredRectOnScreen(XYPlot,X,dx,Y,dy,brush,pen)
510   		continue
420     	CONTINUE
		endif
		currentcolor = 1

	case default
	write(*,*)' Illegal call to subroutine Save2D.... iNEW = ',iNew
	pause 'hit return to continue'
!	call fss_closepicture(2)
	end select
      return
      end

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      Subroutine Plot2D()

!     Routine to save/print/plot 2-D thermal structure
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: XYPlot
	TYPE(AWE_CanvasPen) :: pen
	Type(AWE_CanvasBrush) :: brush

	Include "ThrustnHeat.inc"
	INCLUDE "PlotStuff.inc"				! in folder "AWE_Subs"
! *************************************************************

      integer i,j,k,num,iwind,iplot,iok,ith,iup
      real*4 x,y
	integer*4 Tint,shift
      Character cnum*128
! *********************************************************
! ************************************
!      real xor,xmin,xmax,xlen,yor,ymin,ymax,ylen,overlap
!      integer nxtic,nytic,nxdec,nydec
!      character*40 xlab,ylab,pltitle

      data shift /0/

      iplot = 1
      go to 5

200   continue
      WRITE(*,*)' '
      WRITE(*,*)' Plotting options'
      WRITE(*,*)' 0 = Return'
      write(*,*)' 1 = Plot the grid'
      WRITE(*,*)' 2 = Draw thrust'
      write(*,*)' 3 = Plot rock positions'
	write(*,*)' 4 = Plot T at each grid point'
!      write(*,*)' 6 = Plot Temprature grades'
      read(*,*)(iplot)
      
      if(iplot.eq.0)return

! ----------------------------------------------------------------------
5     continue
      if(iplot.eq.1)then

!  
!  Plot sections
!  First set up scaling
!  Total horizontal distance is nvert*gridX
!  This is set in the main program so that all plots are drawn to the same
!  horizontal scale.
!  Vertical distance is set at 50 km
!  
      xor=25
      xmin=0
!     note: Xdist is in meters
      xmax=gridX*(nvert-1)
      xlen=18            ! 0 for autoscaling
      yor=70
!     Note:YDIST is in km
      ymin=-40
      ymax=0
      ylen=14
      nxstep=0
      nystep=0
      nxdec=0
      nydec=0
      pltitle=' '
      xlab='  '
      ylab='  '

!      call setplt(xor,xmin,xmax,xlen,nxtic,nxdec,XLAB,yor,ymin,ymax,ylen,nytic,nydec,YLAB,pltitle,iok)
      call setplot(iok)
      if(iok.eq.1)go to 200

!      CALL user(xor,xmin,xmax,xlen,yor,ymin,ymax,ylen)

      write(*,*)' '
      write(*,*)'  '
!      write(*,*)' Select window for plotting'
!      do 206 i=2,maxwindo
!      write(*,*)'Window number=',i
!206   continue
!      write(*,*)' New window  =',maxwindo+1
!      read(*,*)(iwind)
!      if(iwind.gt.maxwindo)then
!        if(maxwindo.gt.5)then
!          write(*,*)' Sorry, only 5 windows possible'
!          write(*,*)' Hit return when ready'
!          pause
!          go to 200
!	  endif
!        maxwindo=maxwindo+1

!      call opwndo(2,'Thrust''n Heat 1.0-Graphic',www,W2Ylow,W2Xlow,
!     &W2Yhigh,W2Xhigh,W2Ymax,W2Xmax)
       shift = shift + 10
!      Call fss_openwindow(maxwindo+11,1,'ThrustnHeat Graphic',W2Ylow+shift,W2Xlow+shift,W2Yhigh+shift,W2Xhigh+shift)

!        call opwndo(maxwindo,'Thrust Heat-Graphic',www,
!     &W2Ylow+shift,W2Xlow+shift,W2Yhigh+shift,W2Xhigh+shift,
!     &W2Ymax,W2Xmax)
!      endif
!      
! plot the grid-------------------------

!      call nuwind(iwind)
!      call nuport(0,iwind)
!      call cls
!	call clgs(iwind)
!	call fss_openpicture(iwind)

! plot surface, begin from the first point
      iup = 0
      do 101 k=1,nvert
        x = (k-1) * gridX
        y = float(seaLevel-surface(k))*gridZ/1000.
!        call plot(x,y,iup)
        call plot(XYPlot,x,y,iup)
        iup = 1
101    continue

! plot the horizontal line which is the original surface of the reggion

!      call plot(0.,0.,0)
      call plot(XYPlot,0.,0.,0)
      x= gridX*(nvert-1)
!      call plot(x,0.,1)
      call plot(XYPlot,x,0.,1)

! draw each vertical line. Each line may be of different height.

      Do 205 k=1,nvert

! Draw each vertical section

      x = (k-1) * gridX
      call plot(XYPlot,x,0.,0)
      do 210 i = surface(k),gridBottom

! find the depth of the point 
      y = float(i - seaLevel) * gridZ/1000.    ! plot in km, not meters
      y = -y

! if the point is higher than the original surface, do not plot it

      if(y.gt.0) then
        goto 210
      endif
!      call plot(x,y,1)
      call plot(XYPlot,x,y,1)

! draw tics on last line

      if(k.eq.nvert)then
! Draw tic mark
!         call fss_line(-5,0)
      endif

210   continue
205   continue

!     Draw and label first line
      x=0
!      call plot(x,0.,0)
      call plot(XYPlot,x,0.,0)
!      do 212 i=1,NGRID(1)   !chen
      do 212 i = surface(1),gridBottom
!      y = -gridZ*i/1000.
      y = float(i - seaLevel)*gridZ/1000.
      y = -y
      if(y.Lt.ymin)goto 214
!      call plot(x,y,1)
      call plot(XYPlot,x,y,1)
! Draw tic mark
!         call fss_line(0,5)
212   continue
214   continue
! label y axis
      do 213 i=0,1000,5
      y = -gridZ*float(i)/1000.
      if(y.Lt.ymin)goto 217
      x=0
!      call plot(x,y,0)
      call plot(XYPlot,x,y,0)
!      call fss_line(5,0)    ! draw horizontal tic
!      call fss_move(-30,0)  ! move to left
      num=-y
      WRITE(*,215)num
215   format(i4)
213   continue
217   continue

!     plot the surface of the earth
! plot surface, begin from the first point

!      call nuport(1,iwind)
!      call updscr(iwind)
	!call fss_closepicture(iwind)      
      go to 200      

      endif       ! end plotting grid

!      
! -----------------------------------------------------------------------
      
      if(iplot.eq.2)then

      do 312 ith = 1,nThrust

!      WRITE(*,*)' Input number of thrust fault to plot'
!      read(*,*)(ith)

!     plot thrust fault

!	call fss_openpicture(iwind)
      x=0.0
      y = float(thrustGridZ(ith,1)-seaLevel)*gridZ/1000.
      y = -y
      call plot(XYPlot,x,y,0)
      do 303 k=1,nvert
        x=(k-1)*gridX
        y = float(thrustGridZ(ith,k)-seaLevel)*gridZ/1000.
        y = -y
!        call plot(x,y,1)
        call plot(XYPlot,x,y,1)
303   continue
!      call nuport(1,iwind)
!      call updscr(iwind)
!	call fss_closepicture(iwind)
312   continue

      go to 200
!     end plotting thrust
      endif

! -----------------------------------------------------------------------
! 	 Plot all temperatures
      if(iplot.eq.4)then
301     continue
!	  call fss_openpicture(iwind)
!	  call fss_textsize(8)
      goto 333
333   continue

        Do 320 k=1,nvert
          x = (k-1) * gridX
          do 330 i=surface(k),gridBottom
            y = float(i - seaLevel)*gridZ/1000.
            y = -y
            CALL Plot(XYPlot,x,y,0)
! 		Label
	      Tint=tC(k,i)		! convert to integer
	      write(cnum,'(I4)')Tint
	      cnum=trim(cnum)
!      call plot(XYPlot,x,y,0)
	call TextOnPlot(XYPlot,x,y,cnum,12)
!	    call fss_drawstring(cnum)
330       continue

320     continue
!        call fss_textsize(12)
!	  call fss_closepicture(iwind)
      
        go to 200 

! end plotting temperatures

      endif

! -----------------------------------------------------------------------
!     Plot rock positions
      if(iplot.eq.3)then
!	  call fss_openpicture(iwind)
        do 420 j=1,nRocks
           x = float(rockGridX(j,mct)-1)*gridX
          if(x.lt.0.0) then
            goto 510    ! do not plot the rock
          endif
          y = -float(rockGridZ(j,mct)-seaLevel)*gridZ/1000.
         call symb(XYPlot,x,y,6,4)
510   continue
420     CONTINUE
!	  call fss_closepicture(iwind)
        go to 200 
      endif       ! end plotting rock positions



      goto 200


      end


! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
!      Subroutine SaveGTherm()
      Subroutine SaveGT()

!     Routine to save/geotherms

      implicit none
	Include "ThrustnHeat.inc"
! ************************************************************

!      real xor,xmin,xmax,xlen,yor,ymin,ymax,ylen
!      integer nxtic,nytic,nxdec,nydec
!      character*40 xlab,ylab,pltitle
      integer isave,iout,ivert,i,iend,istart,k,idone,status
!      LOGICAL*4 stdfil,iokL
!      integer*2 vref

1     continue
      WRITE(*,*)' Menu options:'
      WRITE(*,*)' 0 = EXIT'
      WRITE(*,*)' 1 = Save 1 geotherm to disk for future use'
      WRITE(*,*)' 2 = Save all geotherms to disk'
      read(*,*)(isave)
      if(isave.eq.0)return

! --------------------------------------------------------------
      if(isave.eq.1)then
      iout=20
      write(*,*)' Input name of file to save geotherm'
 
       	open(20,file='',status='unknown', IOSTAT=status)
	if(status.ne.0)go to 1
        INQUIRE (20, NAME=outfile)

!      outfile='Gtherm-'//modelin
!       iokL = stdfil(2,vref,'Input save file name',outfile,1,'TEXT') ! 2 is setfile and set vref
!       if(.not.iokL)go to 1
!       open(20,file=outfile,status='unknown')

      write(*,*)' Input number of vertical column to output (0 to exit)'
      read(*,*)ivert
      if(ivert.eq.0)then
         close(20)
         go to 1
         endif
          
      WRITE(20,*)' Geotherms'
      WRITE(20,*)' Vertical section number =  ',ivert
      do 110 i=surface(ivert),gridBottom
        write(20,112)(i-surface(k))*gridZ,tC(ivert,i)
112     format(' ',18F10.1)
110   continue

! Fill in the geotherm to gridpoint 200 assuming 
! a constant gradient below the lower boundary
! compute delz below base of grid
! (note base of grid is at depth of L(1)

!      ngridb=ngrid(ivert)   !chen
!      do 113 i = surface(ivert)+ngrid(ivert)+2,200
!      delz      = (i-1-surface(ivert))*grid - L(ivert)
!      write(20,112) (i-1-surface(k))*grid,
!     &tC(ivert,ngridb)+(qStar(ivert)/kond(ivert,NGRIDB))*DELZ
! 113   continue

!     fss modify 2/17/95
!      do 113 i = gridBottom + 1,200
!      delz  = (i-1-surface(ivert))*grid - L(ivert)    !delz is distance below base of model
!      write(20,112) (i-1-surface(k))*grid,
!     &tC(ivert,gridBottom)+(qStar(ivert)/kond(ivert,gridBottom))*DELZ
! 113   continue

      close(20)
      endif                   ! end isave=1


! -----------------------------------------------------------------------
      if(isave.eq.2)then

      iout=20
      write(*,*)' Input name of file to save geotherm'
       	open(20,file='',status='unknown', IOSTAT=status)
	if(status.ne.0)go to 1
        INQUIRE (20, NAME=outfile)

!      outfile='Gtherm-'//modelin
!      iokL = stdfil(2,vref,'Input save file name',outfile,1,'TEXT') ! 2 is setfile and set vref
!      if(.not.iokL)go to 1
!       open(20,file=outfile,status='unknown')
!      call saveas(20,iok,outfile,'Input file name')
!      if(iok.eq.1)go to 1
          
      iend=0
      istart=1
220   continue
      if(iend.ge.nvert)then
            close(20)
            go to 1
            endif

      istart=idone+1
      iend=istart+20
      if(iend.gt.nvert)iend=nvert
      WRITE(20,*)' Geotherms'
      WRITE(20,*)' Vertical sections = ',istart,iend
      
!      do 210 i=1+surface(1),ngrid(1)+surface(1)  !this is a bug
!      write(20,212) ((i-1)*GRID,tC(k,i),k=istart,iend)
! 212    format(' ',21F7.1)
! 210   continue

!     fss modify 2/17/95
!      do 210 i=1+surface(1),ngrid(1)+surface(1)  !this is a bug
!      write(20,212) ((i-1)*GRID,tC(k,i),k=istart,iend)
! 212    format(' ',21F7.1)
! 210   continue

! Fill in the geotherm to gridpoint 200 assuming 
! a constant gradient below the lower boundary
! compute delz below base of grid
! (note base of grid is at depth of L(1)

!      ngridb=ngrid(ivert)   !chen
!      do 213 i=surface(ivert)+ngrid(ivert)+1,200
!      delz      = (i-1-surface(ivert))*grid - L(ivert)
!      write(20,212)(i-1-surface(k))*grid,
!     &(tC(k,ngridb)+(qStar(1)/kond(1,NGRIDB))*DELZ,k=istart,iend)
! 213   continue

      go to 220

      endif                   ! end isave=2
              

      go to 1

      end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
!      Subroutine GeoTpl(isect,time,xmax,ymax)
      Subroutine GeoTpl(isect,time)
      
!     subroutine to plot a current geotherm for any vertical section
!     isect is the section to plot
!     if isect=0 then the user is prompted to plot any vertical section

      implicit none
	Include "ThrustnHeat.inc"
! *****************************************
	INCLUDE "PlotStuff.inc"				! in folder "AWE_Subs"
! *****************************************

! **********************************************************
      integer i,k,isect
	real*8 time
! **********************************************************

      if(isect.eq.0)then      ! user selected plot
10      continue
        WRITE(*,*)' Current vertical sections are:'
        write(*,*)' Twigland.........Rootzone'
        WRITE(*,'(20i3)')(i,i=1,nvert)
        WRITE(*,*)' 0 = Return'
        WRITE(*,*)' k = Specify geotherm to plot'
        read(*,*)(k)
        if(k.le.0)return

        if(k.gt.nvert)go to 10

	  do 20 i = surface(k),gridBottom
20      tSteady(i)=tC(k,i)

!	  call fss_openpicture(4)
	currentcolor = 3
!		call scolor(icolor)
!! need pass surface position and gridpoint here
        call GEOPLT(tSurface,gridZ,xmaxTD,ymaxTD,tSteady,time,surface(k),gridBottom)
	currentcolor = 1
!	call scolor(1)
!	  call fss_closepicture(4)
!	call scolor(1)
        go to 10
        endif

        if(isect.ne.0)then          ! autoplot
        if(isect.lt.1.or.isect.gt.nvert)return

        k=isect   
        do 30 i=surface(k),gridBottom
30      tSteady(i)=tC(k,i)
!	  call fss_openpicture(4)
        call GEOPLT(tSurface,gridZ,xmaxTD,ymaxTD,tSteady,time,surface(k),gridBottom)
!	  call fss_closepicture(4)
        return
      endif


      return
      end




! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************
      SUBROUTINE GEOPLT(TS,GRID,xmax,ymax,TPLOT,time,surface,gridBottom)

!      Subroutine to plot geotherm contained in array TPLOT
!      ATTENTION: THE NGRID IS MODIFIED WHEN THE ROUTINE IS CALLED
!            THE NGRID HERE IN THIS ROUTINE IS THE NUMBER OF 
!            SPACE (IT IS NOT THE NUMBER OF GRID POINT)

	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: XYPlot

      real*8 ts,grid,xmax,ymax,tplot(500),time
      integer surface,gridBottom,iup

      integer i
      real*4 x,y
      character*6 cnum
      
!      x=TS
!      y=0.
!      CALL PLOT(x,y,0)
!      DO 10 i=1,NGRID
	iup = 0
      do 10 i = surface,gridBottom
!      j=i+1
      x=Tplot(i)
      y=(i-surface)*GRID/1000.
      if(x.le.xmax.and.y.le.ymax)then      ! inside diagram, so plot
         call plot(XYPlot,x,y,iup)
 	   iup = 1
      else           ! outside of diagram -- draw label then return
         go to 20
      endif
         
10    CONTINUE

20    continue
!     done -- now label the curve
         x = Tplot(i-5)             ! pick a spot to draw label
         y = (i-surface-5)*grid/1000
!         call plot(XYPlot,x,y,0)
	   write(cnum,'(F5.0)')time
   	   cnum=adjustl(cnum)
	   cnum=trim(cnum)
!         call fss_overstring(cnum)
	call TextOnPlot(XYPlot,x,y,cnum,12)
         return         
      end

! **********************************************************************
! **********************************************************************
! **********************************************************************
! **********************************************************************

!	SUBROUTINE symbplTh (x,y,itype,xsize,ysize)
!  	routine to plot oval inside rect at point x,y (here, rect = size of symbol)
!     symbols are
!       0 = blank
!     1 = open oval
!     2 = filled oval
!     3 = open rectangle
!     4 = filled rectangle

!	implicit none
! *****************************************
!	INCLUDE "PlotStuff.inc"				! in folder "AWE_Subs"
! *****************************************

!	integer*2 rect(4)
!      integer*4 itype  
!      integer*4 xsize,ysize,ix,iy
!	real*8 x,y
!	save

!     calculate position of symbol in pixels
!	ix=int(xor+(x-xmin)*xconv)
!	iy=int(yor-(y-ymin)*yconv)

!     calculate new rect from ofset ix,iy
!      rect(1) = iy - ysize
!      rect(2) = ix - xsize
!      rect(3) = iy + ysize
!      rect(4) = ix + xsize


!      go to (1,2,3,4)itype

!1     call draovl(rect)
!      return

!2     call paiovl(rect)
!      return

!3     call drarec(rect)
!      return

!4     call pairec(rect)
!      return

!      end



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine FSS_alert(title, text)
    	use AWE_Interfaces
!	interface Alert
!	subroutine AWE_alertBox(title, text)
	character(len=*) :: title, text
	call AWE_alertBox(title,text)
	end
!	end subroutine AWE_alertBox
!	end interface
!	title is used as the title of the alert box text is the text that will be displayed in it.


