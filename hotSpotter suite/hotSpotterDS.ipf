#pragma rtGlobals=1		// Use modern global access method and strict wave access.
Function hotSpotterDS(prefix, prefTheta, preFilter, postFilter, multiThreshold, noiseThreshold)
	String prefix
	Variable prefTheta, preFilter, postFilter, multiThreshold, noiseThreshold
	
	//This program calls hotSpotter() to create dF/F matrices for each (8) direction then calculates
	//DSi and theta for each pixel (point in space) as you would to create polar plots. 
	//Input arguments:
	//prefix = wave prefix in quotes ("dir")
	//prefTheta = preferred direction of cell calculated from spikes
	//preFiler and postFilter = turn on/off pre and post filters in hotSpotter() (0 or 1)
	//multiThreshold = turn on/off multiplication of direction matrix averages to create a template
	//noiseThreshold = 1 or 0 (turn on/off the > bslnnoise requirment to register dF for each pixel)
	//representing the areas in which the dendrites overlap across all stimulations. Used for thresholding
	//in hotSpotter() [selecting pixels for which to calculate dF/F. Under threshold scores as 0]	
	
	//Added functionality 2016-02-16: looping over multiple trials. New naming convention is prefixwaveNum_trialNum
	//e.g. dir0_1, ..., dir7_1, dir0_2, ..., dir7_2, etc. Adhering to this naming scheme is important for compatibility with crossPlotter()
	
	String pathToMatrix, winTitle
	Variable i, j, xPos, yPos, xsum, ysum, radius, theta, DSi, DSi_PN
	Variable nullIndex, prefIndex, corrTheta, relTheta, relThetaAbs
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45} //order and direction of light stimulation
	Make/O/N=(numpnts(angle)) xpts, ypts, dFValue
	
	Variable xWidth, yHeight, frames
	Variable xOffset, yOffset, xDelta, yDelta
	Variable bslnStart, bslnEnd, peakStart, peakEnd
	
	Variable count //generic loop counter
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 1; i < 12; i += 1)
		pathToMatrix = GetDataFolder(1) + prefix + "0_" + num2str(i)
		if(exists(pathToMatrix))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Variable/G threshold //declared globally so that hotSpotter() can use the same value
	
	if(multiThreshold == 1)
		threshold = 1e+14 //1e+15
	else
		threshold = 20 //threshold to use on avg matrix thresholding (when multiThreshold is off)
	endif
	
	//To get appropriate matrix dimensions
	WAVE Matrix = $(GetDataFolder(1) + prefix + "0_" + num2str(trials[0])) //first trial of first direction as example
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
			
	//prompt to set integration times
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
	
	//create 'raw_ _AVG's of each direction here, so hotSpotter() may access them when called
	for(i = 0; i < 8; i += 1)
		Make/O/N=(xWidth, yHeight, frames) $("raw_" + prefix + num2str(i) + "_AVG")
		//now point to newly created husk
		WAVE avgRAWMatrix = $(GetDataFolder(1) + "raw_" + prefix + num2str(i) + "_AVG")
		
		SetScale/P x xOffset, xDelta, "m", avgRAWMatrix 
		SetScale/P y yOffset, yDelta, "m", avgRAWMatrix
		
		
		for(j = 0; j < numpnts(trials); j += 1)
			//point to current original 3D stack
			WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + "_" + num2str(trials[j]))

			//add current stack into avgRAWMatrix to be averaged (operated on later)
			avgRAWMatrix += Matrix 
		endfor
		
		//compute the average of Ca++ signal for the 8 directions (pre-dF/F0 calculation)
		avgRAWMatrix /= numpnts(trials)
	endfor
	
	Make/O/N=(xWidth, yHeight, frames) megaRawAvgMatrix
	SetScale/P x xOffset, xDelta, "m", megaRawAvgMatrix 
	SetScale/P y yOffset, yDelta, "m", megaRawAvgMatrix
	
	//create 'allDir_'s (8 directions collapsed for each trial) for bslnThresholdingin hotSpotter()
	for(j = 0; j < numpnts(trials); j += 1)
		Make/O/N=(xWidth, yHeight, frames) $("allDir_" + prefix + "_" + num2str(trials[j]))
		WAVE avgAllDirMatrix = $(GetDataFolder(1) + "allDir_" + prefix + "_" + num2str(trials[j]))
		
		SetScale/P x xOffset, xDelta, "m", avgAllDirMatrix 
		SetScale/P y yOffset, yDelta, "m", avgAllDirMatrix
		
		for(i = 0; i < 8; i += 1)
			//point to current original 3D stack
			WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + "_" + num2str(trials[j]))

			avgAllDirMatrix += Matrix 
		endfor
		megaRawAvgMatrix += avgAllDirMatrix
		avgAllDirMatrix /= 8	
		
	endfor
	
	megaRawAvgMatrix /= (8*numpnts(trials))
	
	for(j = 0; j < numpnts(trials); j += 1)		
		
		if(multiThreshold == 1)
			Make/O/N=(xWidth,yHeight) avgMatrix
			SetScale/P x xOffset, xDelta, "m", avgMatrix
			SetScale/P y yOffset, yDelta, "m", avgMatrix
			for(i = 0; i < numpnts(angle); i += 1)
				WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + "_" + num2str(trials[j]))
		
				for (yPos = 0; yPos < yHeight; yPos += 1)
					for(xPos = 0; xPos < xWidth; xPos += 1)
						MatrixOp/O zWave = beam(Matrix, xPos, yPos)   //Create a zWave for current pixel over time
						zWave -= 20 //baseline subtraction to keep dark areas with low numbers for thresholding
						avgMatrix[xPos][yPos] = mean(zWave) //Plug the mean intensity for this pixel in to avgMatrix
					endfor
				endfor
				
				if(i == 0)
					Duplicate/O avgMatrix multiMatrix //begin multiMatrix
				else
					MatrixOp/O multiMatrix = multiMatrix * avgMatrix //multiply current direction against growing multiMatrix
				endif
			endfor
			 
		endif
		
		//Generate dF/F0s
		for(i = 0; i < numpnts(angle); i += 1)
			//all arguments provided by arguments given to hotSpotterDS() [this function] except trials[j] which is current trial
			hotSpotter(prefix, i, trials[j], preFilter, postFilter, multiThreshold, noiseThreshold, 0, bslnStart, bslnEnd, peakStart, peakEnd) //rawAVG set to 0 (using trials)
		endfor
		
		Make/O/N=(xWidth, yHeight) thetaMatrix
		SetScale/P x xOffset, xDelta, "m", thetaMatrix
		SetScale/P y yOffset, yDelta, "m", thetaMatrix
		Duplicate/O thetaMatrix DSiPNMatrix, DSiMatrix, relThetaMatrix, relThetaMatrixAbs //same wave scaling
		
		//Make/O/N=(xWidth*yHeight) allDSi, allRelTheta, allRelThetaAbs, allTheta //1D arrays of output values for all pixels
		
		for(yPos = 0; yPos < yHeight; yPos += 1)
			for(xPos = 0; xPos < xWidth; xPos += 1)
				for(i = 0; i < numpnts(angle); i += 1)
					//for each pixel create a wave of all 8 dirs, then calculate DSi and theta 
				
					//point Matrix to current direction 
					WAVE Matrix = $(GetDataFolder(1) + "dF_" + prefix + num2str(i) + "_" + num2str(trials[j])) 
					
					if(Matrix[xPos][yPos] > 0) //discard negative dF/F values
						dFValue[i] = Matrix[xPos][yPos]
					else
						dFValue[i] = 0
					endif
				endfor
				
				xpts = dFValue*cos(angle*pi/180) //convert each value to x and y  coordinates
				ypts = dFValue*sin(angle*pi/180)
				xsum = sum(xpts) //sum across all 8 directions
				ysum = sum(ypts) 
				radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
				theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
				
				if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
					corrTheta = 360 + theta
				elseif(theta == 0)
					corrTheta = NaN //otherwise everything that does not have a value looks like pref direction
				else
					corrTheta = theta
				endif
				
				relTheta = prefTheta - corrTheta // = degrees off from spiking determined preferred direction
				if(relTheta > 180)
					relTheta -= 360 //Make sure we're on -180 to +180 scale
				elseif(relTheta < -180)
					relTheta += 360
				endif
				
				relThetaAbs = abs(relTheta)
				DSi = radius/sum(dFValue)
				
				Wavestats/Q dFValue
				nullIndex = V_minloc
				prefIndex = V_maxloc
				DSi_PN = (dFValue[prefIndex] - dFValue[nullIndex])/(dFValue[prefIndex] + dFValue[nullIndex])
				
				DSiPNMatrix[xPos][yPos] = DSi_PN
				DSiMatrix[xPos][yPos] = DSi
				thetaMatrix[xPos][yPos] = corrTheta //actual theta (non-centred)
				relThetaMatrix[xPos][yPos] = relTheta
				relThetaMatrixAbs[xPos][yPos] = relThetaAbs
				
				//the following are 1D waves of all DSi, relTheta and theta values (for all pixels) for scatter plots
				//allDSi[xPos+(yPos*xWidth)] = DSi
				//allRelTheta[xPos+(yPos*xWidth)] = relTheta
				//allRelThetaAbs[xPos+(yPos*xWidth)] = relTheta
				//allTheta[xPos+(yPos*xWidth)] = prefTheta - corrTheta
			endfor
		endfor
		
		//renaming all output matrices according to current trial (to avoid overwritting)
		Duplicate/O multiMatrix $("multiMatrix_" + num2str(trials[j]))
		Duplicate/O DSiPNMatrix $("DSiPNMatrix_" + num2str(trials[j]))
		Duplicate/O DSiMatrix $("DSiMatrix_" + num2str(trials[j]))
		Duplicate/O thetaMatrix $("thetaMatrix_" + num2str(trials[j]))
		Duplicate/O relThetaMatrix $("relThetaMatrix_" + num2str(trials[j]))
		Duplicate/O relThetaMatrixAbs $("relThetaMatrixAbs_" + num2str(trials[j]))
		
		//Duplicate/O allDSi $("allDSi_" + num2str(trials[j]))
		//Duplicate/O allRelTheta $("allRelTheta_" + num2str(trials[j]))
		//Duplicate/O allRelThetaAbs $("allRelThetaAbs_" + num2str(trials[j]))
		//Duplicate/O allTheta $("allTheta_" + num2str(trials[j]))
	endfor
	
	//############### variance of dF(across trials) for each pixel ################
				
	for (i = 0; i < numpnts(angle); i += 1)
		
		//3D matrix to stack up dFs in
		Make/O/N=(xWidth,yHeight,numpnts(trials)) $("stackDF_" + prefix + num2str(i))
		WAVE stackMatrix = $(GetDataFolder(1) + "stackDF_" + prefix + num2str(i))
		
		Make/O/N=(xWidth,yHeight) $("dFvar_" + prefix + num2str(i))
		WAVE varMatrix = $(GetDataFolder(1) + "dFvar_" + prefix + num2str(i))
		
		SetScale/P x xOffset, xDelta, "m", varMatrix
		SetScale/P y yOffset, yDelta, "m", varMatrix
		
		//assemble dFs of all trials in to a 3D matrix (stack)
		for(j = 0; j < numpnts(trials); j += 1) 
			WAVE Matrix = $(GetDataFolder(1) + "dF_" + prefix + num2str(i) + "_" + num2str(trials[j]))
			stackMatrix[][][j] = Matrix[p][q] //p and q specifies all rows and columns 
		endfor	
		
		for (yPos = 0; yPos < yHeight; yPos += 1)
			
			for (xPos = 0; xPos < xWidth; xPos += 1)
				if (multiMatrix[xPos][yPos] > threshold)
					MatrixOp/O zWave = beam(stackMatrix, xPos, yPos)
					Wavestats/Q zWave
					varMatrix[xPos][yPos] = V_sdev^2
				else
					varMatrix[xPos][yPos] = NaN
				endif
			endfor	
		endfor
		
	endfor
	
	//############### raw AVG theta/Dsi computations################
	
	//Generate dF/F0s from raw signal avgs
	for(i = 0; i < numpnts(angle); i += 1)
		//all arguments provided by arguments given to hotSpotterDS() [this function] except trials[j] which is current trial
		hotSpotter(("raw_"+prefix), i, 42, preFilter, postFilter, multiThreshold, noiseThreshold, 1, bslnStart, bslnEnd, peakStart, peakEnd)////rawAVG set to 1 (suffix = "_AVG", disregard trial)
	endfor
	
	//Do this again with the dF matrices calculated from raw avgs (one for each direction) to properly calculate thetas and DSi
	Make/O/N=(xWidth, yHeight) thetaMatrixAVGraw
	SetScale/P x xOffset, xDelta, "m", thetaMatrixAVGraw
	SetScale/P y yOffset, yDelta, "m", thetaMatrixAVGraw
	Duplicate/O thetaMatrixAVGraw DSiMatrixAVGraw, relThetaMatrixAVGraw, relThetaMatrixAbsAVGraw
	
	for(yPos = 0; yPos < yHeight; yPos += 1)
			for(xPos = 0; xPos < xWidth; xPos += 1)
				for(i = 0; i < numpnts(angle); i += 1)
					//for each pixel create a wave of all 8 dirs
					//then calculate DSi and theta 
					
					//point Matrix to current direction 
					WAVE Matrix = $(GetDataFolder(1) + "dF_" + "raw_" + prefix + num2str(i) + "_AVG") 
					if(Matrix[xPos][yPos] > 0) //discard negative dF/F values
						dFValue[i] = Matrix[xPos][yPos]
					else
						dFValue[i] = 0
					endif
				endfor
				
				xpts = dFValue*cos(angle*pi/180) //convert each value to x and y  coordinates
				ypts = dFValue*sin(angle*pi/180)
				xsum = sum(xpts) //sum across all 8 directions
				ysum = sum(ypts) 
				radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
				theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
				
				if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
					corrTheta = 360 + theta
				elseif(theta == 0)
					corrTheta = NaN //otherwise everything that does not have a value looks like pref direction
				else
					corrTheta = theta
				endif
				
				relTheta = prefTheta - corrTheta // = degrees off from spiking determined preferred direction
				if(relTheta > 180)
					relTheta -= 360 //Make sure we're on -180 to +180 scale
				elseif(relTheta < -180)
					relTheta += 360
				endif
				relThetaAbs = abs(relTheta)
				DSi = radius/sum(dFValue)
				
				Wavestats/Q dFValue
				
				DSiMatrixAVGraw[xPos][yPos] = DSi
				thetaMatrixAVGraw[xPos][yPos] = corrTheta //actual theta (non-centred)
				relThetaMatrixAVGraw[xPos][yPos] = relTheta
				relThetaMatrixAbsAVGraw[xPos][yPos] = relThetaAbs
			endfor
		endfor
		
		//clean up working waves
		Killwaves/Z angle, dFValue, ypts, xpts, avgMatrix, multiMatrix, tempMatrix, peakMatrix, bslnMatrix
		Killwaves/Z allDF, relThetaMatrixAbs, relThetaMatrix, DSiMatrix, thetaMatrix, DSiPNMatrix, allTheta
		Killwaves/Z allRelThetaAbs, allRelTheta, allDSi
end