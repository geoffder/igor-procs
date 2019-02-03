#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//This function calculates the difference in theta at given distances to quantify the spatial acuity of DS
//Uses relThetaMatrix (user is prompted to specify suffix e.g. "_1" or "AVGraw" 
//Prompts user to draw a freehand ROI in which to perform the theta comparisons

Function spaceDiffer(importROIs)
	Variable importROIs
	
	String suffix, roiName, winTitle, ROIfolder
	Variable i, j, diffY, diffX, first, distMode, binSize, lastLevelCross
	Variable maxDist, numROIs, xWidth, yHeight, yPos, xPos, yPos1, xPos1
	
	Prompt suffix, "relThetaMatrix suffix: "
	Prompt maxDist, "Max distance to compare: "
	Prompt numROIs, "Number of ROIs: "
	DoPrompt "Enter parameters", suffix, maxDist, numROIs
	
	Print "Suffix: " + suffix + ", Maximum distance: " + num2str(maxDist)
	Print "Number of ROIs: " + num2str(numROIs)
	
	WAVE Matrix = $(GetDataFolder(1) +"relThetaMatrix" + suffix) //first trial of first direction as example
	xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
	yHeight = DimSize(Matrix, 1)
	
	if(!importROIs)
		for(j =1; j <= numROIs; j += 1) //start with differROI1
			roiName = "differROI" + num2str(j) 
			
			print "Draw differ ROI ..." + num2str(j) 
		
			GetWindow kwTopWin wtitle //get window title for active graph
			winTitle = S_value
			
			//Create continue button window and place it next to the active graph
			NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
			DoWindow/C PauseForROI //Name the continue button window
			AutoPositionWindow/E/M=1/R=$("Graph" + num2str(str2num(winTitle[5,7])))
			
			//Draw text and button on Continue window
			DrawText 21,20,"Draw ROI" + num2str(j) + ", then" 
			DrawText 21,40,"Click Continue."
			Button button0,pos={80,58},size={92,20},title="Continue" //make button
			Button button0,proc=PauseForROI_ContButtonProc //point button to procedure that kills the window

			GraphWaveDraw/F=3/O/L/T $(roiName + "_y"), $(roiName + "_x")
			PauseForUser PauseForROI, $("Graph" + num2str(str2num(winTitle[5,7]))) //needs to be just Graph# part of title
			GraphNormal //end draw mode and return graph to normal
		endfor
	elseif(importROIs)
		Prompt ROIfolder, "Folder with ROIs to import: "
		DoPrompt "Leave blank if in current folder", ROIfolder
		
		if(!stringmatch(ROIfolder, "")) //if not in current folder
			for(j = 1; j <= numROIs; j += 1)
				WAVE ROI = $("root:" + ROIfolder + ":differROI" + num2str(j) + "_y")
				Duplicate/O ROI $("differROI" + num2str(j) + "_y")
				WAVE ROI = $("root:" + ROIfolder + ":differROI" + num2str(j) + "_x")
				Duplicate/O ROI $("differROI" + num2str(j) + "_x")
			endfor
		endif
	endif
	
	//Loop through all dendrite (non-NaN) points within the ROI, skipping others, and compare each to
	//the surrounding points a certain distance away. 
	
	for(j = 1; j <= numROIs; j += 1)
		//Point to correct X and Y roi waves (define hand-drawn ROI)
		WAVE roiY = $(GetDataFolder(1) + "differROI" + num2str(j) + "_y")
		WAVE roiX = $(GetDataFolder(1) + "differROI" + num2str(j) + "_x")
		
		//Convert drawn "waves" into a mask
		ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(xWidth),height=(yHeight),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
		WAVE ROIMask = $(GetDataFolder(1) + "M_ROIMask")
		
		//Prepare output
		Make/O/N=(maxDist+1) $("spaceDiffAVG_ROI" + num2str(j))
		WAVE spaceDiffAVG = $("spaceDiffAVG_ROI" + num2str(j))
		Make/O/N=(maxDist+1) $("spaceDiffSEM_ROI" + num2str(j))
		WAVE spaceDiffSEM = $("spaceDiffSEM_ROI" + num2str(j))
		Make/O/N=(maxDist+1) $("spaceDiffSD_ROI" + num2str(j))
		WAVE spaceDiffSD = $("spaceDiffSD_ROI" + num2str(j))
		Make/O/N=(maxDist+1) $("spaceDiffMED_ROI" + num2str(j))
		WAVE spaceDiffMED = $("spaceDiffMED_ROI" + num2str(j))
			
		for(i = 1; i <= maxDist; i += 1)
			Make/O/N=1 $("spaceDiff" + suffix + "_ROI" + num2str(j) + "_d" + num2str(i))
			WAVE spaceDiffW = $("spaceDiff" + suffix + "_ROI" + num2str(j) + "_d" + num2str(i))
			first = 1
			
			for(yPos = 0; yPos < yHeight; yPos += 1)
				
				for(xPos = 0; xPos < xWidth; xPos += 1)
					//if pixel is NaN or not within ROI, skip it
					if(numtype(Matrix[xPos][yPos]) == 2 || ROIMask[xPos][yPos] == 1) 
						continue //advance to next pixel
					endif
					
					for(diffY = yPos-i; diffY <= yPos+i; diffY += 1)
						
						for(diffX = xPos-i; diffX<= xPos+i; diffX += 1)
							
							if(diffX < 0 || diffY < 0 || diffX >= xWidth || diffY >= yHeight)
								continue //skip if out of bounds
							elseif(abs(xPos - diffX) != i && abs(yPos - diffY) != i)
								continue //Skip all comparisons not on current 'concentric ring' (skip inside)
							elseif(numtype(Matrix[diffX][diffY]) == 2 || ROIMask[diffX][diffY] == 1)
							 	continue //as before, skip those not inside ROI (mask 1s) or on dendrite (NaNs)
							elseif(abs(xPos - diffX) == i &&  abs(yPos - diffY) == i)
								//continue //skip corner/diagonal points
							endif
							
							if (first == 1)
								spaceDiffW[0] = {abs(Matrix[xPos][yPos] - Matrix[diffX][diffY])}
								first = 0
							else
								spaceDiffW[numpnts(spaceDiffW)] = {abs(Matrix[xPos][yPos] - Matrix[diffX][diffY])}
							endif
						endfor // end of diffX loop (inner ring loop)
					endfor //end of diffY loop (outer ring loop)	
				endfor			
			endfor
			WaveStats/Q/W spaceDiffW 
			WAVE stats = M_WaveStats // W flag in WaveStats puts results in M_WaveStats
			spaceDiffAVG[i] = stats[3] //mean
			spaceDiffSEM[i] = stats[26] //standard error
			spaceDiffSD[i] = stats[4] //standard deviation
			
			Duplicate/O spaceDiffW spaceDiffWmed
			Sort spaceDiffWmed spaceDiffWmed
			spaceDiffMED[i] = spaceDiffWmed[(numpnts(spaceDiffWmed)-1)/2] //median
			
		endfor //end of dist loop (i)	
		
		binSize = 1
		distMode = 1 //originally optional, now just leaving on
		
		if(distMode == 1)
			Make/O/N=1 $("thetaDiffWave_ROI" + num2str(j))
			WAVE diffWave = $(GetDataFolder(1) + "thetaDiffWave_ROI" + num2str(j))
				
			Make/O/N=1 $("thetaDistWave_ROI" + num2str(j))
			WAVE distWave = $(GetDataFolder(1) +"thetaDistWave_ROI" + num2str(j))
			
			first = 1
			
			for (yPos = 0; yPos < yHeight; yPos += 1)
				for(xPos = 0; xPos < xWidth; xPos += 1)
					if(numtype(Matrix[xPos][yPos]) == 2 || ROIMask[xPos][yPos] == 1) 
						continue //advance to next pixel
					endif
					
					for(yPos1 = yPos; yPos1 < yHeight; yPos1 += 1)
						for(xPos1 = xPos+1; xPos1 < xWidth; xPos1 += 1)
							if(numtype(Matrix[xPos1][yPos1]) == 2 || ROIMask[xPos1][yPos1] == 1) 
								continue //advance to next pixel
							endif
							
							if (first == 1)
								diffWave[0] = {abs(Matrix[xPos][yPos] - Matrix[xPos1][yPos1])}
								//Calculate distance between the two points being compared
								distWave[0] = {sqrt((xPos1 - xPos)^2 + (yPos1 - yPos)^2)}
								first = 0
							else
								diffWave[numpnts(diffWave)] = {abs(Matrix[xPos][yPos] - Matrix[xPos1][yPos1])}
								distWave[numpnts(distWave)] = {sqrt((xPos1 - xPos)^2 + (yPos1 - yPos)^2)}
							endif
						endfor
					endfor
				endfor
			endfor
			
			Sort distWave, distWave, diffWave
			
			Make/O/N=1  $("diffBinsAVG_ROI" + num2str(j))
			WAVE diffBinsAVG =  $("diffBinsAVG_ROI" + num2str(j))
			Make/O/N=1  $("diffBinsMED_ROI" + num2str(j))
			WAVE diffBinsMED =  $("diffBinsMED_ROI" + num2str(j))
			Make/O/N=1  $("diffBinsSD_ROI" + num2str(j))
			WAVE diffBinsSD =  $("diffBinsSD_ROI" + num2str(j))
			Make/O/N=1  $("diffBinsSEM_ROI" + num2str(j))
			WAVE diffBinsSEM =  $("diffBinsSEM_ROI" + num2str(j))
			
			SetScale/P x binSize, binSize, "pixels", diffBinsAVG, diffBinsMED, diffBinsSD, diffBinsSEM
			
			Make/O/N=1  $("distBins_ROI" + num2str(j))
			WAVE distBins =  $("distBins_ROI" + num2str(j))
			
			first = 1
			lastLevelCross = 0
			for(i = binSize; i <= maxDist;i += binSize)
				
				FindLevel/Q/P/EDGE=1/R=(lastLevelCross) distWave, i
				
				if(first == 1)
					Duplicate/O/R=[lastLevelCross,V_LevelX] diffWave, diffBinSort 
					
					WaveStats/Q/W diffBinSort
					WAVE stats = M_WaveStats // W flag in WaveStats puts results in M_WaveStats
					diffBinsAVG[0] = stats[3] //mean
					diffBinsSEM[0] = stats[26] //standard error
					diffBinsSD[0] = stats[4] //standard deviation
					
					Sort diffBinSort,  diffBinSort
					diffBinsMED[0] = diffBinSort[(numpnts(diffBinSort)-1)/2] //median
			
					distBins[0] = i
					first = 0
				else
					Duplicate/O/R=[lastLevelCross+1,V_LevelX] diffWave, diffBinSort 
					
					WaveStats/Q/W diffBinSort
					WAVE stats = M_WaveStats // W flag in WaveStats puts results in M_WaveStats
					diffBinsAVG[numpnts(diffBinsAVG)] = {stats[3]} //mean
					diffBinsSEM[numpnts(diffBinsSEM)] = {stats[26]} //standard error
					diffBinsSD[numpnts(diffBinsSD)] = {stats[4]} //standard deviation
					
					Sort diffBinSort,  diffBinSort
					diffBinsMED[numpnts(diffBinsMED)] = {diffBinSort[(numpnts(diffBinSort)-1)/2]} //median
					
					distBins[numpnts(distBins)] = {i}
				endif
				
				lastLevelCross = V_LevelX
			endfor
			
		endif
		
	endfor //end of ROI loop
	
	if (numROIs > 1)
		Make/O/N=(numpnts(spaceDiffAVG)) $("spaceDiffAVG_AVGN" + num2str(numROIs))
		WAVE spaceDiffAVG_AVG = $("spaceDiffAVG_AVGN" + num2str(numROIs))
		Make/O/N=(numpnts(spaceDiffMED)) $("spaceDiffMED_AVGN" + num2str(numROIs))
		WAVE spaceDiffMED_AVG = $("spaceDiffMED_AVGN" + num2str(numROIs))
		Make/O/N=(numpnts(diffBinsAVG)) $("diffBinsAVG_AVGN" + num2str(numROIs))
		WAVE diffBinsAVG_AVG = $("diffBinsAVG_AVGN" + num2str(numROIs))
		Make/O/N=(numpnts(diffBinsMED))  $("diffBinsMED_AVGN" + num2str(numROIs))
		WAVE diffBinsMED_AVG = $("diffBinsMED_AVGN" + num2str(numROIs))
		
		SetScale/P x binSize, binSize, "pixels", diffBinsAVG_AVG, diffBinsMED_AVG
		
		for(j = 1; j <= numROIs; j +=1)
			WAVE spaceDiffAVG = $("spaceDiffAVG_ROI" + num2str(j))
			WAVE spaceDiffMED = $("spaceDiffMED_ROI" + num2str(j))
			WAVE diffBinsAVG = $("diffBinsAVG_ROI" + num2str(j))
			WAVE diffBinsMED = $("diffBinsMED_ROI" + num2str(j))
			
			spaceDiffAVG_AVG +=  spaceDiffAVG
			spaceDiffMED_AVG +=  spaceDiffMED
			diffBinsAVG_AVG += diffBinsAVG
			diffBinsMED_AVG += diffBinsMED
		endfor
		
		spaceDiffAVG_AVG /=  numROIs
		spaceDiffMED_AVG /=  numROIs
		diffBinsAVG_AVG /= numROIs
		diffBinsMED_AVG /= numROIs
	endif
	KillWaves spaceDiffWmed
End
