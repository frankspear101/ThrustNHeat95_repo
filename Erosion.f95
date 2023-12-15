	Subroutine Erosion(ith)
	implicit none
	Include "ThrustnHeat.inc"
! ************************************************************
	integer ith,k,surfaceOld,j
	real*8 dtime,ElevationMetersOld,sumTime,Erate
! ************************************************************
	
!	Erosion model is hyperbolic with elevation
!	y(erosion rate) = ((10/(10-X))-1)*factor
!	Asymptote is 10 km -- mountains can't get higher than that.

!	Tectonics can move the mountains above 10 km, so there is no solution to the equation.
!	Set a maximum erosion rate for any mountain over .9*10 = 9 km high

!	we also need to iterate so that all erosion isn't done in 1 time step with a huge erosion rate.
!	The time step for each modelstep is the relaxtime variable.
!	So after thrusting we need to (for each vertical section)
!	(1) calculate the new surface elevation = (surface - sealevel)*gridZ
!	(2) calculate the erosion rate for this elevation (in km or meters) (Erate)
!	(3) calculate a time step that doesn't result in too much erosion. 
!		To do this, I'll pick a nominal amount of erosion (say 1 meter) and calculate ∆t as 1meter/Erate
!		Keep iterating till the end of the modelstep time

!	Another issue is that erosion is in meters but the surface is in grid points. 
!	I think I'll keep track of where the surface is in real space (meters) and then add or subtract grid points if needed
!		but only after an integral value of the grid point is overstepped (e.g. when erosion>gridZ)

!     Erode the surface
	if(ErosionSwitch(ith).gt.0.)then
		do 10 k = 1,nvert
		ElevationMetersOld = ElevationMeters(k)
		sumtime = 0.0d0
20		continue
!		Erate is scaled for mm/year = km/Ma
		Erate = (10.0d0/(10.0d0-ElevationMeters(k)/1000.0d0) - 1.0d0)*ErosionFactor(ith)
		if(Erate.gt.1.0D-2)then		! only erode if the rate > 0.01 mm/year (km/Ma)
			Erate = Erate/1000.0d0		! erosion rate in m/year
			dtime = 1.0d0/Erate		! scale dtime to erode only 1 meter at a time - in years
			sumtime = sumtime + dtime	! sumtime is in years
			if (sumTime.gt.relaxtime(ith)*1.D6)then		! relaxtime is in Ma
				sumTime = sumTime - dtime
				dtime = relaxtime(ith)*1.0d6 - sumtime
				ElevationMeters(k) = ElevationMeters(k) - Erate*dtime
				go to 30
				endif
			ElevationMeters(k) = ElevationMeters(k) - 1.0d0
			go to 20
30			continue		
!			set the grid point for the surface
			surfaceOld = surface(k)
			surface(k) = sealevel - ElevationMeters(k)/GridZ
			tC(k,surface(k)) = tSurface
!			Since this is erosion, surfaceOld should always be a smaller number than surface(k)
!			Only reset the sky if the grid point has moved
			if(surface(k)-surfaceOld.gt.0)then
				do 40 j = surfaceOld,surface(k)-1           ! set temperature of eroded cells to sky T
				tC(k,j) = tSky
40	  			continue
				endif
			!write(*,*)ith,k,SurfaceOld,Surface(k),ElevationMetersOld,ElevationMeters(k)
			endif
10		continue
		endif

	return
	end
		