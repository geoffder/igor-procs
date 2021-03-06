#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <FilterDialog> menus=0
#include <All IP Procedures>
#include <Image Saver>
Function hotSpotter (prefix, waveNum, trialNum, preFilter, postFilter, multiThreshold, noiseThreshold, rawAVG, bslnStart, bslnEnd, peakStart, peakEnd)
	String prefix
	Variable waveNum, trialNum, preFilter, postFilter, multiThreshold, noiseThreshold, rawAVG
	Variable bslnStart, bslnEnd, peakStart, peakEnd
	//This program creates a dF/F matrix for a given calcium imaging stack (name = prefix + waveNum).
	//Input arguments:
	//prefix = wave prefix in quotes ("dir")
	//waveNum = number of desired wave (direction) with the given prefix
	//trialNum = trial number of desired wave (direction) with the given prefix e.g. prefixwaveNum_trialNum
	//preFilter = turn on/off spatial filtering of the raw image stacks before dF/F calculations on individual pixels (0 or 1)
	//postFilter = turn on/off spatial filtering on the 2D dF/F matrix at end of function (0 or 1)
	//multiThreshold = turn on/off threshold/mask that determines which pixels to calculate dF/F for based on the
	//intensity of the corresponding pixel in the multiplication matrix generated by hotSpotterDS(). If you intend
	//to use this option when calling hotSpotter() directly, make sure you have already generated this multiMatrix.
	//noiseThreshold = 1 or 0 (turn on/off the > bslnnoise requirment to register dF for each pixel)
	////rawAVG = 1 or 0 (disregard trial and use "_AVG" tag instead)
	
	Variable filterSize = 3 //n x n size of filter (pre and post, median by default) (n = 3 by default)
	
	NVAR threshold = $(GetDataFolder(1) + "threshold") //generated in hotSpotterDS
	
	Variable bslnThreshMode, bslnThreshold, rawBsln
	
	//experimental conditionals/parameters
	bslnThreshMode = 0
	bslnThreshold = 120
	
	
	if(!rawAVG)
		WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(waveNum) + "_" + num2str(trialNum))
	elseif(rawAVG)
		WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(waveNum) + "_AVG")
	endif
	
	SetScale/P x 0,2.5e-07,"", Matrix //3.333e-07
	SetScale/P y 0,2.5e-07,"", Matrix
							
	//point to the multiMatrix in this folder generated by hotSpotterDS()
	WAVE multiMatrix = $(GetDataFolder(1) + "multiMatrix")
	
	Duplicate/O Matrix $("work_" + prefix + num2str(waveNum) + "_" + num2str(trialNum)), rawWorkMatrix
	WAVE workMatrix = $(GetDataFolder(1) + "work_" + prefix + num2str(waveNum) + "_" + num2str(trialNum))
	
	Variable/G xWidth, yHeight, frames
	Variable xOffset, yOffset, xDelta, yDelta
	Variable xPos, yPos, zPos, z, darkROIavg
	
	xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
	yHeight = DimSize(Matrix, 1)
	frames = DimSize(Matrix, 2)
	
	//Get scaling from original stacks so all output can be the same
	//xOffset = DimOffset(Matrix, 0)
	//xDelta = DimDelta(Matrix, 0)
	//yOffset = DimOffset(Matrix, 1)
	//yDelta = DimDelta(Matrix, 1)
	
	//replaced xOffset and yOffset with 0 to avoid inconsistencies (comment out to return)
	xOffset = 0
	yOffset = 0
	xDelta = .00000025 //3 pixels per micron
	yDelta = .00000025
	
	//###### Dark Fluorescence #########
	//Convert drawn "wave" into a mask
	ImageBoundaryToMask ywave=$("darkROI_y"),xwave=$("darkROI_x"),width=(dimSize(Matrix,0)),height=(dimSize(Matrix,1)),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
	WAVE ROIMask = $(GetDataFolder(1) + "M_ROIMask")	

	Make/O/N=(frames) ROIavg
	for(z = 0; z < frames; z += 1)
		ImageStats/M=1/P=(z)/R=ROImask Matrix
		ROIavg[z] = V_avg
	endfor
	darkROIavg = mean(ROIavg)
	KillWaves/Z ROIavg
	//################################
	
	Make/O/N=(xWidth,yHeight) avgMatrix, tempMatrix

//	if (preFilter == 1)
//		for (zPos = 0; zPos < frames; zPos +=1)
//			for (yPos = 0; yPos < yHeight; yPos += 1)
//				for(xPos = 0; xPos < xWidth; xPos += 1) 
//					tempMatrix[xPos][yPos] = Matrix[xPos][yPos][zPos] //create 2D matrix on each frame of the stack
//				endfor
//			endfor
//			MatrixFilter/N=(filterSize) avg tempMatrix //set prefilter type here (types: avg, median, etc)
//			workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
//		endfor
//	endif
	
	if (preFilter == 1)
		for (zPos = 0; zPos < frames; zPos +=1)
			tempMatrix[][] = Matrix[p][q][zPos] //create 2D matrix on each frame of the stack
			MatrixFilter/N=(filterSize) avg tempMatrix //set prefilter type here (types: avg, median, etc)
			workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
		endfor
	endif
	
	Variable i
	
	//creation of the mask used for selecting which pixels have dF/F generated, and which have dF/F set to 0
	if(multiThreshold == 1)
		avgMatrix = multiMatrix  //if using multiThreshold, simply set avgMatrix equal to the provided multiMatrix
	else
		for (yPos = 0; yPos < yHeight; yPos += 1)
			for(xPos = 0; xPos < xWidth; xPos += 1)
				MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
				avgMatrix[xPos][yPos] = mean(zWave) //Plug the mean intensity for this pixel in to avgMatrix
			endfor
		endfor
	endif
	
	String deltaMatrix, baselineMatrix
	
	//names for output 2D dF/F matrix
	if(!rawAVG)
		deltaMatrix = "dF_" + prefix + num2str(waveNum) + "_"+ num2str(trialNum)  //giving the new output matrix a name
		baselineMatrix = "bsln_" + prefix + num2str(waveNum) + "_" + num2str(trialNum)
	elseif(rawAVG)
		deltaMatrix = "dF_" + prefix + num2str(waveNum) + "_AVG"
		baselineMatrix = "bsln_" + prefix + num2str(waveNum) + "_AVG" 
	endif
	//Killwaves/Z $deltaMatrix
	
	Make/O/N=(xWidth,yHeight) $baselineMatrix, peakMatrix, $deltaMatrix
	
	WAVE dFMatrix = $(GetDataFolder(1) + deltaMatrix)
	WAVE bslnMatrix = $(GetDataFolder(1) + baselineMatrix)
	SetScale/P x xOffset, xDelta, "m", dFMatrix 
	SetScale/P y yOffset, yDelta, "m", dFMatrix
		
	//KillWaves/Z zWave //kill and delete existing zWave in current folder
	Make/O/N=(frames) zWave  //Holder waves for the single pixel/time zWave and projected exponential baseline
	Make/O/N=(xWidth*yHeight) allDF
	Variable bsln, peak, bslnNoise
	
	
	//allDir avg usage goes in here
	//use trialNum to indicate which to use
	//TESTING: CALC BSLN OF ALLDIR AVG FOR THE CURRENT TRIAL FOR BSLNTHRESHOLD
	if(1 && bslnThreshMode)
		if(trialNum == 42) //if operating on rawAVG, use all scan AVG,trials don't exist for them
			WAVE rawMatrix = $(GetDataFolder(1) + "megaRawAvgMatrix") 
		else
			WAVE rawMatrix = $(GetDataFolder(1) + "allDir_" + prefix + "_" + num2str(trialNum))
		endif
	endif
	
	if(0 && bslnThreshMode)
		if(trialNum == 42)
			WAVE rawMatrix = $(GetDataFolder(1) + prefix + num2str(waveNum) + "_AVG")
		else
			WAVE rawMatrix = $(GetDataFolder(1) + "raw_" + prefix + num2str(waveNum) + "_AVG")
		endif
	endif
	
	if(1 && bslnThreshMode)
		for (zPos = 0; zPos < frames; zPos +=1)
			tempMatrix[][] = rawMatrix[p][q][zPos] //create 2D matrix on each frame of the stack
			MatrixFilter/N=(filterSize) avg tempMatrix //set prefilter type here (types: avg, median, etc)
			rawWorkMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
		endfor
	endif
		
	for (yPos = 0; yPos < yHeight; yPos += 1)
		for(xPos = 0; xPos < xWidth; xPos += 1)
			//if block to ensure these operations only occur on pixels that pass threshold
			if (avgMatrix[xPos][yPos]  > threshold)
				
				//begin caculation after first mask, need bsln for second mask
				MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
				//zWave -= darkROIavg //baseline subtraction
				bsln = mean(zWave,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
				
				if(1 && bslnThreshMode)
					MatrixOp/O zWave = beam(rawWorkMatrix, xPos, yPos)   //Create a zWave for current pixel over time
					zWave -= darkROIavg //baseline subtraction
					rawBsln = mean(zWave,bslnStart,bslnEnd) //baseline before stimulus //default 5,30	
				endif
				
				//If bslnThreshMode is engaged, current position must pass bslnThreshold
				//Otherwise OR statement is satisfied and additional threshold is not applied 
				if(!bslnThreshMode || rawBsln > bslnThreshold) //bsln or rawBsln	
					
					MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
					
					//THIS IS EXPERIMENTAL ALERT LOOK HERE###############################
					Smooth/S=2 9, zWave
					
					//rest of calculations
					bslnNoise = sqrt(variance(zWave,bslnStart,bslnEnd)) //standard deviation of pre-stimulus baseline
					bslnMatrix[xPos][yPos] = bsln
					peak = mean(zWave,peakStart,peakEnd) //make sure time window for peak integration is appropriate //default 36,42
					peakMatrix[xPos][yPos] = peak
					
					//ensure peak of signal is seperate from baseline noise (decrease coefficient to eliminate this filter)
					if(noiseThreshold)
						if(peak > 2.5*bslnNoise && bsln > 0) //2.5*
							dFMatrix[xPos][yPos] = (peak - bsln)/bsln //Plug the calculated dF/F into the output matrix
						else
							dFMatrix[xPos][yPos] = 0
						endif
					else
						dFMatrix[xPos][yPos] = (peak - bsln)/bsln //no thresholding, take dF for all pixels
					endif
					
					allDF[xPos+(yPos*xWidth)] = (peak - bsln)/bsln
				else
					dFMatrix[xPos][yPos] = 0
				endif
				
			else
				dFMatrix[xPos][yPos] = 0
			endif
		endfor
	endfor
	
	//post dF/F calculation filtering
	if(postFilter == 1)	
		MatrixFilter/N=(filterSize) median dfMatrix //set filter type here (avg, median, etc)
	endif
	
	KillWaves/Z zWave
	KillWaves/Z workMatrix, rawWorkMatrix
end