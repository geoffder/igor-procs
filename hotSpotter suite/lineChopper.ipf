#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function lineChopper(prefix, prefTheta, importLines, numLines, smthMode [,suffix])
	String prefix, suffix
	Variable prefTheta, importLines, numLines, smthMode

	Variable pt, h ,i, j, k, m, count1, count2, xDelta, yDelta, avgDelta, first
	Variable frames, bslnStart, bslnEnd, peakStart, peakEnd
	Variable  z, chunkStart, darkROIavg, bslnNoise, x1
	Variable bsln, peak, xsum, ysum, radius, theta, corrTheta, DSi
	Variable lineDivAVG, ptsPerPix, chunkSz
	String pathToMatrix, suffixStr, dirLabel
	String lineName, winTitle, graphname, lineFolder
	
	Variable maxLength = 80 //change to be more clear that it's 3 per micron?
	Variable lineWidth = 5
	
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45}
	Make/O/N=(numpnts(angle)) xpts, ypts, dirDFs
	
	Prompt bslnStart, "bsln start (default = 5): "
	Prompt bslnEnd, "bsln end (default = 28): "
	Prompt peakStart, "peak start (default = 36): "
	Prompt peakEnd, "peak end (default = 42): "
	DoPrompt "Leave blank to use defaults", bslnStart, bslnEnd, peakStart, peakEnd
	
	//when prompt inputs are left blank or user messes up, set integration parameters to defaults
	if(bslnStart >= bslnEnd)
		bslnStart = 5
		bslnEnd = 28
	endif
	if(peakStart >= peakEnd)
		peakStart = 36
		peakEnd = 42
	endif
	
	Print "bslnStart: " + num2str(bslnStart) + ", bslnEnd: "+ num2str(bslnEnd)
	Print "peakStart: " + num2str(peakStart) + ", peakEnd: " + num2str(peakEnd)
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count1 = 0
	for(i = 1; i < 12; i += 1)
		pathToMatrix = GetDataFolder(1) + prefix + "0_" + num2str(i)
		if(exists(pathToMatrix))
			trials[count1] = {i}
			count1 += 1
		endif
	endfor
	
	if(ParamIsDefault(suffix))
		suffixStr = num2str(trials[0])
	else
		suffixStr = suffix
	endif
			
	WAVE Matrix = $(GetDataFolder(1) + prefix + "0_" + suffixStr)
	xDelta = DimDelta(Matrix, 0)
	yDelta = DimDelta(Matrix, 1)
	avgDelta = (xDelta+yDelta)/2
	
	//target window for drawing (for both dark ROI and lines)
	GetWindow kwTopWin wtitle //get window title for active graph
	winTitle = S_value
	
	//################## Dark ROI ################
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
	//##############################################
	
	//############# Create/import lines for profiles ########################	
	if(!importLines)
		for(j = 1; j <= numLines; j += 1) 
			lineName = "line" + num2str(j) 
			
			//Create continue button window and place it next to the active graph
			NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
			DoWindow/C PauseForROI //Name the continue button window
			AutoPositionWindow/E/M=1/R=$("Graph" + num2str(str2num(winTitle[5,7])))
			
			DrawText 21,20,"Draw line" + num2str(j) + ", then" 
			DrawText 21,40,"Click Continue." //window dialogue
			Button button0,pos={80,58},size={92,20},title="Continue" //make button
			Button button0,proc=PauseForROI_ContButtonProc //point button to procedure that kills the window

			GraphWaveDraw/F=3/O/L/T $(lineName + "_y"), $(lineName + "_x") //removed /F=3 for smooth freehand
			PauseForUser PauseForROI, $("Graph" + num2str(str2num(winTitle[5,7]))) //needs to be just Graph# part of title
			GraphNormal //end draw mode and return graph to normal
		endfor
	elseif(importLines) //non-zero importLines specified in arguments
		Prompt lineFolder, "Folder with lines to import: "
		DoPrompt "Leave blank if in current folder", lineFolder
		
		if(!stringmatch(lineFolder, "")) //if not in current folder
			for(j = 0; j < numLines + 1; j += 1)
				WAVE ROI = $("root:" + lineFolder + ":line" + num2str(j) + "_y")
				Duplicate/O ROI $("line" + num2str(j) + "_y")
				WAVE ROI = $("root:" + lineFolder + ":line" + num2str(j) + "_x")
				Duplicate/O ROI $("line" + num2str(j) + "_x")
			endfor
		endif
	endif
	//###############################################################
	
	for(j = 1; j <= numLines; j+=1)
		WAVE lineY = $(GetDataFolder(1) + "line" + num2str(j) + "_y")
		WAVE lineX = $(GetDataFolder(1) + "line" + num2str(j) + "_x")
		
		lineDivAVG = 0
		for(pt = 0; pt <= 30; pt += 1)
			lineDivAVG += sqrt((lineX[pt+1] - lineX[pt])^2 + (lineY[pt+1] - lineY[pt])^2)
		endfor 
		lineDivAVG /= 30
		
		ptsPerPix = round(avgDelta/lineDivAVG)
		Print "------ Line " + num2str(j) + " ------"
		Print "Points per pixel = " + num2str(ptsPerPix)
		
		for(k = 0; k < numpnts(trials); k += 1)//will still run once if there are no trials (e.g. running on AVGs instead)
		
			//set suffix to trial number by default if no suffix (e.g. "AVG") specified
			if(ParamIsDefault(suffix))
				suffixStr = num2str(trials[k])
			else
				suffixStr = suffix
			endif	
			
			//############## Trial Specific 'pref/avg' theta ##################
			if(1) //add another optional argument to activate this process
				
				for(i = 0; i < 8; i += 1) //for each direction 
					
					WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + "_" + suffixStr)
					frames = DimSize(Matrix, 2)
					
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
					//################################
					
					//###### dF/F0 Calculations ##########
					//Make/O/N=(bslnEnd - bslnStart + 1) bsln
					//Make/O/N=(peakEnd - peakStart + 1) peak
					
					for(z = 0; z < frames; z += 1)
						ImageLineProfile/P=(z) xWave=lineX, yWave=lineY, srcwave=Matrix, width=lineWidth //width original =4
						WAVE profile = $(GetDataFolder(1) + "W_ImageLineProfile")
						ROIavg[z] = mean(profile)
					endfor
					
					ROIavg -= darkROIavg
					
					if(smthMode)
						Smooth/S=2 9, ROIavg
					endif 				
					
					bsln = mean(ROIavg,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
					bslnNoise = sqrt(variance(ROIavg,bslnStart,bslnEnd)) //standard deviation of pre-stimulus baseline
					peak = mean(ROIavg,peakStart,peakEnd) //response integration  //default 36,42
					
					//dF/F0 = (peakRegion - bslnRegion)/bslnRegion, if peak larger than noise
					if(peak > bslnNoise*1.5)
						dirDFs[i] = (peak - bsln)/bsln
					else
						dirDFs[i] = 0
					endif
					//################################
					
				endfor //dir loop (i)
				
				//############ Vector Calculations ###########
				xpts = dirDFs*cos(angle*pi/180) //convert each value to x and y  coordinates
				ypts = dirDFs*sin(angle*pi/180)
				xsum = sum(xpts) //sum across all 8 directions
				ysum = sum(ypts) 
				
				prefTheta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
				
				if (prefTheta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
					prefTheta += 360
				endif
				
				print "Trial " + num2str(trials[k]) + " prefTheta = " + num2str(prefTheta)
				//#########################################
			
			endif
			//###############################################################
			
			Make/O/N=(maxLength-1) $("thetaVar_VS_chunkSize_" + suffixStr + "_ln" + num2str(j))
			WAVE thetaVar_VS_chunkSize =  $("thetaVar_VS_chunkSize_" + suffixStr + "_ln" + num2str(j))
			Make/O/N=(maxLength-1) $("degOffAVG_VS_chunkSize_" + suffixStr + "_ln" + num2str(j))
			WAVE degOffAVG_VS_chunkSize = $("degOffAVG_VS_chunkSize_" + suffixStr + "_ln" + num2str(j))
			Make/O/N=(maxLength-1) $("degOffMED_VS_chunkSize_" + suffixStr + "_ln" + num2str(j))
			WAVE degOffMED_VS_chunkSize = $("degOffMED_VS_chunkSize_" + suffixStr + "_ln" + num2str(j))
			SetScale/P x .6666,.3333,"µm", thetaVar_VS_chunkSize, degOffAVG_VS_chunkSize, degOffMED_VS_chunkSize
			
			count1 = 0
			for(h = 2; h <= maxLength; h += 1)		
				
				Make/O/N=1 $("thetas_chnkSz" + num2str(h) + "_" + suffixStr + "_ln" + num2str(j))
				WAVE chunkThetas = $("thetas_chnkSz" + num2str(h) + "_" + suffixStr + "_ln" + num2str(j))
				Make/O/N=1 $("DSis_chnkSz" + num2str(h) + "_" + suffixStr + "_ln" + num2str(j))
				WAVE chunkDSis = $("DSis_chnkSz" + num2str(h) + "_" + suffixStr + "_ln" + num2str(j))
				first = 1 //for conditional used to guide filling of these growing result waves
				
				chunkSz = h * ptsPerPix
				
				//Chop up line from start to finish in chunks of size h and calculate dF/F0 
				for(chunkStart = 0; chunkStart+chunkSz < numpnts(lineX); chunkStart += chunkSz)
				
					Duplicate/O/R=[chunkStart,chunkStart+chunkSz-1] lineX, chopX
					Duplicate/O/R=[chunkStart,chunkStart+chunkSz-1] lineY, chopY
					
					for(i = 0; i < 8; i += 1) //for each direction 
					
						WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + "_" + suffixStr)
						frames = DimSize(Matrix, 2)
						
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
						//################################
						
						//###### dF/F0 Calculations ##########
						
						for(z = 0; z < frames; z += 1)
							ImageLineProfile/P=(z) xWave=chopX, yWave=chopY, srcwave=Matrix, width=4
							WAVE profile = $(GetDataFolder(1) + "W_ImageLineProfile")
							ROIavg[z] = mean(profile)
						endfor
						
						//ROIavg -= darkROIavg
						
						if(smthMode)
							Smooth/S=2 9, ROIavg
						endif 				
						
						bsln = mean(ROIavg,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
						bslnNoise = sqrt(variance(ROIavg,bslnStart,bslnEnd)) //standard deviation of pre-stimulus baseline
						peak = mean(ROIavg,peakStart,peakEnd) //response integration  //default 36,42
						
						//dF/F0 = (peakRegion - bslnRegion)/bslnRegion, if peak larger than noise
						if(peak > bslnNoise*1.5 && bsln > 0)
							dirDFs[i] = (peak - bsln)/bsln
						else
							dirDFs[i] = 0
						endif
						//################################
					
					endfor //dir loop (i)
				
					//############ Vector Calculations ###########
					xpts = dirDFs*cos(angle*pi/180) //convert each value to x and y  coordinates
					ypts = dirDFs*sin(angle*pi/180)
					xsum = sum(xpts) //sum across all 8 directions
					ysum = sum(ypts) 
					radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
					theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
					if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
						corrTheta = 360 + theta
					else
						corrTheta = theta
					endif
					
					if (first == 1)
						chunkThetas[0] = {corrTheta}
						chunkDSis[0] = {radius/sum(dirDFs)}
						first = 0
					else
						chunkThetas[numpnts(chunkThetas)] = {corrTheta}
						chunkDSis[numpnts(chunkDSis)] = {radius/sum(dirDFs)}
					endif
					//#########################################
					
				endfor //chopper loop (chunkStart)
				
				//############### OUTPUTS ###################
				
				//Convert to relative theta measurements for linear visualization/quantification
				for(m = 0; m < numpnts(chunkThetas); m += 1)
					chunkThetas[m] = prefTheta - chunkThetas[m]
					
					if(chunkThetas[m] > 180)
						chunkThetas[m] -= 360 //Make sure we're on -180 to +180 scale
					elseif(chunkThetas[m] < -180)
						chunkThetas[m] += 360
					endif
				endfor
				
				Duplicate/O chunkThetas chunkThetas_degsOff
				chunkThetas_degsOff = abs(chunkThetas)
				
				thetaVar_VS_chunkSize[count1] = variance(chunkThetas) //done on degsOff because it's linear (not circular)
				degOffAVG_VS_chunkSize[count1] = mean(chunkThetas_degsOff)
				
				//determine median by sorting and picking middle
				Sort chunkThetas_degsOff chunkThetas_degsOff
				degOffMED_VS_chunkSize[count1] = chunkThetas_degsOff[(numpnts(chunkThetas_degsOff)-1)/2]
				
				count1 += 1
				//#############################################
			endfor //length loop (h), go on to next chunk size
			
			//Add up results from all trials
			if(numpnts(trials) > 1)
				if(k == 0)
					Duplicate/O thetaVar_VS_chunkSize $("thetaVar_VS_chunkSize_trAV_ln" + num2str(j))
					Duplicate/O degOffAVG_VS_chunkSize $("degOffAVG_VS_chunkSize_trAV_ln" + num2str(j))
					Duplicate/O degOffMED_VS_chunkSize $("degOffMED_VS_chunkSize_trAV_ln" + num2str(j))
					WAVE thetaVar_VS_chunkSize_trialAV =  $("thetaVar_VS_chunkSize_trAV_ln" + num2str(j))
					WAVE degOffAVG_VS_chunkSize_trialAV = $("degOffAVG_VS_chunkSize_trAV_ln" + num2str(j))
					WAVE degOffMED_VS_chunkSize_trialAV = $("degOffMED_VS_chunkSize_trAV_ln" + num2str(j))
				else
					thetaVar_VS_chunkSize_trialAV +=  thetaVar_VS_chunkSize
					degOffAVG_VS_chunkSize_trialAV += degOffAVG_VS_chunkSize
					degOffMED_VS_chunkSize_trialAV += degOffMED_VS_chunkSize
				endif
			endif
			
		endfor //trial loop (k)
		
		//Divide results by number of trials (calculate average)
		if(numpnts(trials) > 1)
			thetaVar_VS_chunkSize_trialAV /= numpnts(trials)
			degOffAVG_VS_chunkSize_trialAV /= numpnts(trials)
			degOffMED_VS_chunkSize_trialAV /= numpnts(trials)
		endif
			
	endfor //numlines loop (j)
	
	KillWaves/Z bsln, peak, chunkThetas_degsOff, chopX, chopY
	KillWaves/Z profile, W_LineProfileX, W_LineProfileY, ROIavg
	KillWaves/Z xpts, ypts, dirDFs, M_ROImask, M_WaveStats
End