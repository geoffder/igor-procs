#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <FilterDialog> menus=0
#include <All IP Procedures>
#include <Image Saver>
Function miniSpotter (prefix, preFilter, postFilter)
	String prefix
	Variable preFilter, postFilter
	//This program creates a dF/F matrix for a given calcium imaging stack (name = prefix + trialNum).
	//Input arguments:
	//prefix = wave prefix in quotes ("dir")
	//trialNum = trial number of desired wave with the given prefix
	//preFilter = turn on/off spatial filtering of the raw image stacks before dF/F calculations on individual pixels (0 or 1)
	//postFilter = turn on/off spatial filtering on the 2D dF/F matrix at end of function (0 or 1)
	
	Variable filterSize = 3 //n x n size of filter (pre and post, median by default) (n = 3 by default)
	Variable threshold = 50 //to throwout non-dendrite regions
	Variable smthSize = 7 //smoothing function bin size
	Variable start = 6 //frame to start calculations from (skip laser response)
	
	Variable xWidth, yHeight, frames
	Variable xOffset, yOffset, xDelta, yDelta
	Variable xPos, yPos, zPos, z, darkROIavg
	Variable count, i, darkSubAvg
	String winTitle, pathToMatrix
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 20; i += 1)
		pathToMatrix = prefix + num2str(i)
		if(exists(pathToMatrix))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	WAVE Matrix = $(prefix + num2str(trials[0]))
	xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
	yHeight = DimSize(Matrix, 1)
	frames = DimSize(Matrix, 2)
	
	//replaced xOffset and yOffset with 0 to avoid inconsistencies (comment out to return)
	xOffset = 0
	yOffset = 0
	xDelta = .00000025 //3 pixels per micron
	yDelta = .00000025
	
	for	(i = 0; i < numpnts(trials); i += 1)
		WAVE Matrix = $(prefix + num2str(trials[i]))
		//**work around** set scaling to a default to make ROIs work
		SetScale/P x xOffset, xDelta, "m", Matrix 
		SetScale/P y yOffset, yDelta, "m", Matrix
	endfor
	//SetScale/P x 0,2.5e-07,"", Matrix //3.333e-07
	//SetScale/P y 0,2.5e-07,"", Matrix
							
	//Get scaling from original stacks so all output can be the same
	//xOffset = DimOffset(Matrix, 0)
	//xDelta = DimDelta(Matrix, 0)
	//yOffset = DimOffset(Matrix, 1)
	//yDelta = DimDelta(Matrix, 1)
	
	//###### Dark Fluorescence (using ROI) #########
	//target window for drawing (for both dark ROI and lines)
	GetWindow kwTopWin wtitle //get window title for active graph
	winTitle = S_value
	
	//Create continue button window and place it next to the active graph
	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
	DoWindow/C PauseForROI //Name the continue button window
	AutoPositionWindow/E/M=1/R=$("Graph" + num2str(str2num(winTitle[5,7])))
			
	DrawText 21,20,"Draw darkROI, then"	
	DrawText 21,40,"Click Continue."
	
	Button button0,pos={80,58},size={92,20},title="Continue" //make button
	Button button0,proc=PauseForROI_ContButtonProc //point button to procedure that kills the window

	GraphWaveDraw/F=3/O/L/T $("darkROI_y"), $("darkROI_x")
	PauseForUser PauseForROI, $("Graph" + num2str(str2num(winTitle[5,7]))) //needs to be just Graph# part of title
	GraphNormal //end draw mode and return graph to normal
	
	//Convert drawn "wave" into a mask
	ImageBoundaryToMask ywave=$("darkROI_y"),xwave=$("darkROI_x"),width=(dimSize(Matrix,0)),height=(dimSize(Matrix,1)),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
	WAVE ROIMask = $("M_ROIMask")	

	Make/O/N=(frames) ROIavg
	for(z = 0; z < frames; z += 1)
		ImageStats/M=1/P=(z)/R=ROImask Matrix
		ROIavg[z] = V_avg
	endfor
	darkROIavg = mean(ROIavg)
	KillWaves/Z ROIavg
	//################################
	
	for	(i = 0; i < numpnts(trials); i += 1)
		WAVE Matrix = $(prefix + num2str(trials[i]))
		Duplicate/O Matrix $("work_" + prefix + num2str(trials[i]))
		WAVE workMatrix = $("work_" + prefix + num2str(trials[i]))
		Make/O/N=(xWidth,yHeight) $("avg_"+ prefix + num2str(trials[i])), tempMatrix
		WAVE avgMatrix = $("avg_"+ prefix + num2str(trials[i]))
		Make/O/N=(xWidth,yHeight) $("darkSubAvg_"+ prefix + num2str(trials[i])), tempMatrix
		WAVE darkSubAvgMatrix = $("darkSubAvg_"+ prefix + num2str(trials[i]))
		
		if (preFilter == 1)
			for (zPos = 0; zPos < frames; zPos +=1)
				tempMatrix[][] = Matrix[p][q][zPos] //create 2D matrix on each frame of the stack
				MatrixFilter/N=(filterSize) avg tempMatrix //set prefilter type here (types: avg, median, etc)
				workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
			endfor
		endif
		
		//creation of the mask used for selecting which pixels have dF/F generated, and which have dF/F set to 0
		//also used as the baseline (F) for dendrite pixels that have their dF/F calculated
		for (yPos = 0; yPos < yHeight; yPos += 1)
			for(xPos = 0; xPos < xWidth; xPos += 1)
				MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
				avgMatrix[xPos][yPos] = mean(zWave,start,numpnts(zWave)-1) //Plug the mean intensity for this pixel in to avgMatrix
				darkSubAvgMatrix[xPos][yPos] = avgMatrix[xPos][yPos] - darkROIavg
				if(darkSubAvgMatrix[xPos][yPos] < 0)
					darkSubAvgMatrix[xPos][yPos] = 0
				endif
			endfor
		endfor
		
		//output 2D dF/F matrix
		Make/O/N=(xWidth,yHeight) $("dF_" + prefix + num2str(trials[i]))
		WAVE dFMatrix = $("dF_" + prefix + num2str(trials[i]))
		
		SetScale/P x xOffset, xDelta, "m", dFMatrix 
		SetScale/P y yOffset, yDelta, "m", dFMatrix
		
		//Holder wave for the single pixel/time zWave	
		Make/O/N=(frames) zWave  
			
		for (yPos = 0; yPos < yHeight; yPos += 1)
			for(xPos = 0; xPos < xWidth; xPos += 1)
				//if block to ensure these operations only occur on pixels that pass threshold
				if (darkSubAvgMatrix[xPos][yPos]  > threshold)
					//create beam of current pixel and clean it up
					MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
					zWave -= darkROIavg //dark baseline subtraction
					Smooth/S=2 smthSize, zWave
					
					//calculate largest deviation from baseline for whole beam 
					darkSubAvg = avgMatrix[xPos][yPos] - darkROIavg
					dFMatrix[xPos][yPos] = (wavemax(zWave,start,numpnts(zWave)-1)- darkSubAvg)/darkSubAvg
					//dFMatrix[xPos][yPos] = (wavemax(zWave,start,numpnts(zWave)-1)-avgMatrix[xPos][yPos])/avgMatrix[xPos][yPos]
				else
					dFMatrix[xPos][yPos] = 0
				endif
			endfor
		endfor
		
		//post dF/F calculation filtering
		if(postFilter == 1)	
			MatrixFilter/N=(filterSize) median dfMatrix //set filter type here (avg, median, etc)
		endif
	endfor
		
	KillWaves/Z zWave, workMatrix, tempMatrix
end