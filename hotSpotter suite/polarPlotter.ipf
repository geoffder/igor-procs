#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <New Polar Graphs>

//Before running, ensure graph you would like to draw ROIs on is topmost/active
//Enter command (ex: polarPlotter("raw_dir","_AVG",0,3) for 3 ROIs).
//You'll be prompted to draw a dark region, then a number of signal ROIs you specified in the args.
//Program will spit out ROIavgs (z projections), dFs, thetas, DSis and whole polar plots.
//For now, run multiple times using suffix (e.g. "_1") to deal with multiple trials

Function polarPlotter(prefix, suffix, importROIs, numROIs, smthMode)
	String prefix, suffix
	Variable importROIs, numROIs, smthMode
	
	Variable i, j, k, frames, darkAvg, bsln, bslnNoise, peak, bsln_notSmth
	Variable xsum, ysum, radius, theta, corrTheta, DSi
	Variable bslnStart, bslnEnd, peakStart, peakEnd
	
	String roiName, winTitle, graphname, ROIfolder
	
	Variable radMax = 1.2
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45}
	Make/O/N=(numpnts(angle)) xpts, ypts, ROIdF
	
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
	
	//output preferred thetas and DSis for each ROI into waves
	//position in wave corresponds to number of ROI (0 is dark, don't use) 
	Make/O/N=(numROIs+1) $(prefix + suffix + "_ROIs" + "_thetas")
	WAVE thetas = $(GetDataFolder(1) + prefix + suffix + "_ROIs" + "_thetas")
	Make/O/N=(numROIs+1) $(prefix + suffix + "_ROIs" + "_rel0thetas")
	WAVE rel0thetas = $(GetDataFolder(1) + prefix + suffix + "_ROIs" + "_rel0thetas")
	Make/O/N=(numROIs+1) $(prefix + suffix + "_ROIs" + "_DSis")
	WAVE DSis = $(GetDataFolder(1) + prefix + suffix + "_ROIs" + "_DSis")
	
	if(!importROIs)
		for(j = 0; j < numROIs + 1; j += 1) //0 is dark ROI
			roiName = "ROI" + num2str(j) 
			if(j == 0)
				print "Draw dark fluorescence ROI..."
			else
				print "Draw signal ROI" + num2str(j) 
			endif
			GetWindow kwTopWin wtitle //get window title for active graph
			winTitle = S_value
			
			//Create continue button window and place it next to the active graph
			NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
			DoWindow/C PauseForROI //Name the continue button window
			AutoPositionWindow/E/M=1/R=$("Graph" + num2str(str2num(winTitle[5,7])))
			
			if(j == 0)
				DrawText 21,20,"Draw dark ROI, then"
			else
				DrawText 21,20,"Draw ROI" + num2str(j) + ", then" 
			endif
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
			for(j = 0; j < numROIs + 1; j += 1)
				WAVE ROI = $("root:" + ROIfolder + ":ROI" + num2str(j) + "_y")
				Duplicate/O ROI $("ROI" + num2str(j) + "_y")
				WAVE ROI = $("root:" + ROIfolder + ":ROI" + num2str(j) + "_x")
				Duplicate/O ROI $("ROI" + num2str(j) + "_x")
			endfor
		endif
	endif
	
	// make ROI Mask from ROI
	for(j = 0; j < numROIs+1; j += 1) //0 = dark ROI, 1+ are signal ROIs
		WAVE roiY = $(GetDataFolder(1) + "ROI" + num2str(j) + "_y")
		WAVE roiX = $(GetDataFolder(1) + "ROI" + num2str(j) + "_x")
		
		for(i = 0; i < 8; i += 1) //for each direction
			
			WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + suffix) //matrix mask is applied to
			
			//Convert drawn "wave" into a mask
			ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(dimSize(Matrix,0)),height=(dimSize(Matrix,1)),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
			WAVE ROIMask = $(GetDataFolder(1) + "M_ROIMask")
			
			//Calculate ROI averages (inside of drawn masks) using ImageStats
			frames = DimSize(Matrix, 2)
			if(frames > 0) //if target matrices are 3D stacks (e.g. Ca++ signal)
				Make/O/N=(frames) ROIavg
				for(k = 0; k < frames; k += 1)
					ImageStats/M=1/P=(k)/R=ROImask Matrix
					ROIavg[k] = V_avg
				endfor
			else //if target matrices are 2D (e.g. dF matrix)
				Make/O/N=1 ROIavg
				ImageStats/M=1/R=ROImask Matrix
				ROIavg[0] = V_avg
			endif
			
			if(j == 0)
				darkAvg = mean(ROIavg) //calculate background fluorescence
			else
				ROIavg -= darkAvg //subtrace background fluorescence
				
				if(smthMode) //take care of unsmoothed ROI trace and normalize
					Duplicate/O ROIavg ROIavg_notSmth
					Smooth/S=2 9, ROIavg //Savitzky-Golay
					bsln_notSmth = mean(ROIavg_notSmth,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
					ROIavg_notSmth = (ROIavg_notSmth - bsln_notSmth)/bsln_notSmth
				endif
				
				bsln = mean(ROIavg,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
				bslnNoise = sqrt(variance(ROIavg,bslnStart,bslnEnd)) //standard deviation of pre-stimulus baseline
				peak = mean(ROIavg,peakStart,peakEnd) //response integration  //default 36,42
				
				
				//dF/F0 = (peakRegion - bslnRegion)/bslnRegion, if peak larger than noise
				if(peak > bslnNoise*1.5)
					ROIdF[i] = (peak - bsln)/bsln
				else
					ROIdF[i] = 0
				endif
			endif
			
			//Now suppressing output of non-normalized waveforms, uncomment or add conditional
			//if needed in the future
			//Save with appropriate names
			//if(!smthMode)
			//	Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j))
			//else
			//	Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_smth")
			//endif
			
			if(j > 0)
				//dF/F0 normalization
				ROIavg = (ROIavg - bsln)/bsln
			endif
			
			if(smthMode)
				if(j == 0)
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j))
				else
					Duplicate/O ROIavg_notSmth $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrm")
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrmsm")
				endif
			else
				if(j == 0)
					 Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j))
				else
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrm")
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrmsm")
					Smooth/S=2 9,$(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrmsm")
				endif
			endif

			Duplicate/O ROIdF $(prefix  + suffix + "_ROI" + num2str(j) + "_dF")
			
		endfor
		
	endfor
	
	//Use dFs of each ROI (stored in waves of length = numpnts(angle)) to calculate theta and DSi
	for(j = 1; j < numROIs+1; j += 1)
	
		WAVE dirDFs = $(GetDataFolder(1) + prefix + suffix + "_ROI" + num2str(j) + "_dF")
		
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
			
		thetas[j] = corrTheta
		rel0thetas[j] = theta
		DSis[j] = radius/sum(dirDFs)
	endfor
	
	//Now that we have theta and DSi for each ROI, make polar plots
	for(j = 1; j < numROIs+1; j += 1)
		
		Make/O/N=2 $("ROI" + num2str(j) + "_DSiArrowY")
		WAVE DSiArrowY = $(GetDataFolder(1) + "ROI" + num2str(j) + "_DSiArrowY")
		Make/O/N=2 $("ROI" + num2str(j) + "_thetaArrowX")
		WAVE thetaArrowX = $(GetDataFolder(1) + "ROI" + num2str(j) + "_thetaArrowX")
		
		//set up vector for plot with values from current ROI
		DSiArrowY[0] = 0
		DSiArrowY[1] = DSis[j] * radMax//arrow pointing from origin on scale of 0 to 1 (visually)
		thetaArrowX[0] = thetas[j]
		thetaArrowX[1] = thetas[j] //direction that arrowy points in
		
		//point to directional DFs from current ROI
		WAVE dirDFs = $(GetDataFolder(1) + prefix + suffix + "_ROI" + num2str(j) + "_dF")
		
		//sort angle and corresponding dFs to be clockwise
		Duplicate/O angle angleSort 
		Sort angleSort, angleSort, dirDFs //sort angleSort ascending and sort dirDFs to angleSort
		
		angleSort[numpnts(angleSort)] = {angleSort[0]}
		dirDFs[numpnts(dirDFs)] = {dirDFs[0]}
		
		//Build polar plots
		graphname = "ROI" + num2str(j) + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, dirDFs, angleSort, 360)	
		WMPolarAppendTrace(graphname, DSiArrowY, thetaArrowX, 360)
		WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		WMPolarSetManualRadiusRange(graphname, -.01, radMax)
	endfor
	
	//trim dark roi point out of dsi and theta wave ouptuts
	DeletePoints 0,1, thetas, rel0thetas, DSis
	
	KillWaves/Z ROIavg, ROIavg_notSmth, xpts, ypts, ROIdF //$(GetDataFolder(1) + "M_ROIMask")
end	

Function calGrid(prefix, numROIs, importROIs)
	String prefix
	Variable numROIs, importROIs
	
	Variable i, j, k,f, count, perSpX, perX, perY, perSpY
	String botAxName, leftAxName, pathToMatrix
	Variable smthMode = 1
	
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
	Make/O trials = {1,2,3}
		
	perSpX = .02
	perX = (1-perSpx*7)/8
	perSpY = .03
	perY = (1-perSpY*(numROIs-1))/(numROIs)
	
	Variable frames, darkAvg, bsln, bslnNoise, peak, bsln_notSmth
	Variable xsum, ysum, radius, theta, corrTheta, DSi
	Variable bslnStart, bslnEnd, peakStart, peakEnd
	
	String roiName, winTitle, graphname, ROIfolder
	
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45}
	Make/O dispOrder = {6, 1, 5, 3, 7, 0, 4, 2} // display in this order so 0deg is centred
	Make/O/N=(numpnts(angle)) xpts, ypts, ROIdF
	
	Prompt bslnStart, "bsln start (default = 0): "
	Prompt bslnEnd, "bsln end (default = 10): "
	Prompt peakStart, "peak start (default = 20): "
	Prompt peakEnd, "peak end (default = 25): "
	DoPrompt "Leave blank to use defaults", bslnStart, bslnEnd, peakStart, peakEnd
	
	//when prompt inputs are left blank or user messes up, set integration parameters to defaults
	if(bslnStart >= bslnEnd)
		bslnStart = 0
		bslnEnd = 10
	endif
	if(peakStart >= peakEnd)
		peakStart = 20
		peakEnd = 25
	endif
	
	Print "bslnStart: " + num2str(bslnStart) + ", bslnEnd: "+ num2str(bslnEnd)
	Print "peakStart: " + num2str(peakStart) + ", peakEnd: " + num2str(peakEnd)
	
	if(!importROIs)
		for(j = 0; j < numROIs + 1; j += 1) //0 is dark ROI
			roiName = "ROI" + num2str(j) 
			if(j == 0)
				print "Draw dark fluorescence ROI..."
			else
				print "Draw signal ROI" + num2str(j) 
			endif
			GetWindow kwTopWin wtitle //get window title for active graph
			winTitle = S_value
			
			//Create continue button window and place it next to the active graph
			NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
			DoWindow/C PauseForROI //Name the continue button window
			AutoPositionWindow/E/M=1/R=$("Graph" + num2str(str2num(winTitle[5,7])))
			
			if(j == 0)
				DrawText 21,20,"Draw dark ROI, then"
			else
				DrawText 21,20,"Draw ROI" + num2str(j) + ", then" 
			endif
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
			for(j = 0; j < numROIs + 1; j += 1)
				WAVE ROI = $("root:" + ROIfolder + ":ROI" + num2str(j) + "_y")
				Duplicate/O ROI $("ROI" + num2str(j) + "_y")
				WAVE ROI = $("root:" + ROIfolder + ":ROI" + num2str(j) + "_x")
				Duplicate/O ROI $("ROI" + num2str(j) + "_x")
			endfor
		endif
	endif
	
	/// DARK ROI
	WAVE Matrix = $(prefix+num2str(0)+"_"+num2str(trials[0]))	
	WAVE roiY = $(GetDataFolder(1) + "ROI0_y")
	WAVE roiX = $(GetDataFolder(1) + "ROI0_x")
	//Convert drawn "wave" into a mask
	ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(dimSize(Matrix,0)),height=(dimSize(Matrix,1)),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
	WAVE ROIMask = $(GetDataFolder(1) + "M_ROIMask")
	//Calculate ROI averages (inside of drawn masks) using ImageStats
	frames = DimSize(Matrix, 2)
	Make/O/N=(frames) ROIavg
	for(f = 0; f < frames; f += 1)
		ImageStats/M=1/P=(f)/R=ROImask Matrix
		ROIavg[f] = V_avg
	endfor		
	darkAvg = mean(ROIavg) //calculate background fluorescence
					
	for(i = 1; i <= numROIs; i+=1)
		WAVE roiY = $(GetDataFolder(1) + "ROI" + num2str(i) + "_y")
		WAVE roiX = $(GetDataFolder(1) + "ROI" + num2str(i) + "_x")
		leftAxName = "r"+num2str(i)	
			
		for(j = 0; j < 8; j+=1)
			botAxName = "d"+num2str(j)
			for(k = 0; k < numpnts(trials); k+= 1)
				WAVE Matrix = $(prefix+num2str(dispOrder[j])+"_"+num2str(trials[k]))
			
				//Convert drawn "wave" into a mask
				ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(dimSize(Matrix,0)),height=(dimSize(Matrix,1)),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
				WAVE ROIMask = $(GetDataFolder(1) + "M_ROIMask")
			
				//Calculate ROI averages (inside of drawn masks) using ImageStats
				frames = DimSize(Matrix, 2)
				Make/O/N=(frames) ROIavg
				for(f = 0; f < frames; f += 1)
					ImageStats/M=1/P=(f)/R=ROImask Matrix
					ROIavg[f] = V_avg
				endfor
			
				ROIavg -= darkAvg //subtract background fluorescence
				
				if(smthMode) //take care of unsmoothed ROI trace and normalize
					Duplicate/O ROIavg ROIavg_notSmth
					Smooth/S=2 9, ROIavg //Savitzky-Golay
					bsln_notSmth = mean(ROIavg_notSmth,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
					ROIavg_notSmth = (ROIavg_notSmth - bsln_notSmth)/bsln_notSmth
				endif
				
				bsln = mean(ROIavg,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
				bslnNoise = sqrt(variance(ROIavg,bslnStart,bslnEnd)) //standard deviation of pre-stimulus baseline
				peak = mean(ROIavg,peakStart,peakEnd) //response integration  //default 36,42
			
				
				//dF/F0 normalization
				ROIavg = (ROIavg - bsln)/bsln
				Duplicate/O ROIavg $(prefix+num2str(j)+"_"+num2str(trials[k])+"_r"+num2str(i))
				WAVE dF = $(prefix+num2str(j)+"_"+num2str(trials[k])+"_r"+num2str(i))
				if(i==1 && !j && !k)
					Display/L=$leftAxName/B=$botAxName dF //vs timeAx
				else
					AppendToGraph/L=$leftAxName/B=$botAxName dF //vs timeAx
				endif
				ModifyGraph axisEnab($botAxName)={(j*(perX+perSpX)),(j*(perX+perSpX)+perX)}
				ModifyGraph freePos($botAxName)=0
				
			endfor
		endfor
		ModifyGraph axisEnab($leftAxName)={((i-1)*(perY+perSpY)),((i-1)*(perY+perSpY)+perY)}
		ModifyGraph freePos($leftAxName)=0
	endfor
	Killwaves/Z temp
End

Function newPolarPlotter(prefix, suffix, importROIs, numROIs, smthMode)
	String prefix, suffix
	Variable importROIs, numROIs, smthMode
	
	Variable i, j, k, count, frames, darkAvg, bsln, bslnNoise, peak, bsln_notSmth
	Variable xsum, ysum, radius, theta, corrTheta, DSi
	Variable bslnStart, bslnEnd, peakStart, peakEnd
	
	String roiName, winTitle, graphname, ROIfolder, pathToMatrix
	
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45}
	Make/O/N=(numpnts(angle)) xpts, ypts, ROIdF
	
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
	Make/O trials = {1,2,3}
	
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
	
	//output preferred thetas and DSis for each ROI into waves
	//position in wave corresponds to number of ROI (0 is dark, don't use) 
	Make/O/N=(numROIs+1) $(prefix + suffix + "_ROIs" + "_thetas")
	WAVE thetas = $(GetDataFolder(1) + prefix + suffix + "_ROIs" + "_thetas")
	
	Make/O/N=(numROIs+1) $(prefix + suffix + "_ROIs" + "_DSis")
	WAVE DSis = $(GetDataFolder(1) + prefix + suffix + "_ROIs" + "_DSis")
	
	if(!importROIs)
		for(j = 0; j < numROIs + 1; j += 1) //0 is dark ROI
			roiName = "ROI" + num2str(j) 
			if(j == 0)
				print "Draw dark fluorescence ROI..."
			else
				print "Draw signal ROI" + num2str(j) 
			endif
			GetWindow kwTopWin wtitle //get window title for active graph
			winTitle = S_value
			
			//Create continue button window and place it next to the active graph
			NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
			DoWindow/C PauseForROI //Name the continue button window
			AutoPositionWindow/E/M=1/R=$("Graph" + num2str(str2num(winTitle[5,7])))
			
			if(j == 0)
				DrawText 21,20,"Draw dark ROI, then"
			else
				DrawText 21,20,"Draw ROI" + num2str(j) + ", then" 
			endif
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
			for(j = 0; j < numROIs + 1; j += 1)
				WAVE ROI = $("root:" + ROIfolder + ":ROI" + num2str(j) + "_y")
				Duplicate/O ROI $("ROI" + num2str(j) + "_y")
				WAVE ROI = $("root:" + ROIfolder + ":ROI" + num2str(j) + "_x")
				Duplicate/O ROI $("ROI" + num2str(j) + "_x")
			endfor
		endif
	endif
	
	// make ROI Mask from ROI
	for(j = 0; j < numROIs+1; j += 1) //0 = dark ROI, 1+ are signal ROIs
		WAVE roiY = $(GetDataFolder(1) + "ROI" + num2str(j) + "_y")
		WAVE roiX = $(GetDataFolder(1) + "ROI" + num2str(j) + "_x")
		
		for(i = 0; i < 8; i += 1) //for each direction
			
			WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + suffix) //matrix mask is applied to
			
			//Convert drawn "wave" into a mask
			ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(dimSize(Matrix,0)),height=(dimSize(Matrix,1)),scalingwave=Matrix,seedx=(dimOffset(Matrix,0)+dimDelta(Matrix,0)),seedy=(dimOffset(Matrix,1)+dimDelta(Matrix,1))
			WAVE ROIMask = $(GetDataFolder(1) + "M_ROIMask")
			
			//Calculate ROI averages (inside of drawn masks) using ImageStats
			frames = DimSize(Matrix, 2)
			if(frames > 0) //if target matrices are 3D stacks (e.g. Ca++ signal)
				Make/O/N=(frames) ROIavg
				for(k = 0; k < frames; k += 1)
					ImageStats/M=1/P=(k)/R=ROImask Matrix
					ROIavg[k] = V_avg
				endfor
			else //if target matrices are 2D (e.g. dF matrix)
				Make/O/N=1 ROIavg
				ImageStats/M=1/R=ROImask Matrix
				ROIavg[0] = V_avg
			endif
			
			if(j == 0)
				darkAvg = mean(ROIavg) //calculate background fluorescence
			else
				ROIavg -= darkAvg //subtrace background fluorescence
				
				if(smthMode) //take care of unsmoothed ROI trace and normalize
					Duplicate/O ROIavg ROIavg_notSmth
					Smooth/S=2 9, ROIavg //Savitzky-Golay
					bsln_notSmth = mean(ROIavg_notSmth,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
					ROIavg_notSmth = (ROIavg_notSmth - bsln_notSmth)/bsln_notSmth
				endif
				
				bsln = mean(ROIavg,bslnStart,bslnEnd) //baseline before stimulus //default 5,30
				bslnNoise = sqrt(variance(ROIavg,bslnStart,bslnEnd)) //standard deviation of pre-stimulus baseline
				peak = mean(ROIavg,peakStart,peakEnd) //response integration  //default 36,42
				
				
				//dF/F0 = (peakRegion - bslnRegion)/bslnRegion, if peak larger than noise
				if(peak > bslnNoise*1.5)
					ROIdF[i] = (peak - bsln)/bsln
				else
					ROIdF[i] = 0
				endif
			endif
			
			//Now suppressing output of non-normalized waveforms, uncomment or add conditional
			//if needed in the future
			//Save with appropriate names
			//if(!smthMode)
			//	Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j))
			//else
			//	Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_smth")
			//endif
			
			if(j > 0)
				//dF/F0 normalization
				ROIavg = (ROIavg - bsln)/bsln
			endif
			
			if(smthMode)
				if(j == 0)
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j))
				else
					Duplicate/O ROIavg_notSmth $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrm")
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrmsm")
				endif
			else
				if(j == 0)
					 Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j))
				else
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrm")
					Duplicate/O ROIavg $(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrmsm")
					Smooth/S=2 9,$(prefix + num2str(i) + suffix + "_ROI" + num2str(j) + "_nrmsm")
				endif
			endif

			Duplicate/O ROIdF $(prefix  + suffix + "_ROI" + num2str(j) + "_dF")
			
		endfor
		
	endfor
	
	//Use dFs of each ROI (stored in waves of length = numpnts(angle)) to calculate theta and DSi
	for(j = 1; j < numROIs+1; j += 1)
	
		WAVE dirDFs = $(GetDataFolder(1) + prefix + suffix + "_ROI" + num2str(j) + "_dF")
		
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
			
		thetas[j] = corrTheta
		DSis[j] = radius/sum(dirDFs)
	endfor
	
	//Make/O/N=(numROIs) $(prefix+suffix+"_DSis") $(prefix+suffix+"_thetas")
	//Now that we have theta and DSi for each ROI, make polar plots
	for(j = 1; j < numROIs+1; j += 1)
		
		Make/O/N=2 $("ROI" + num2str(j) + "_DSiArrowY")
		WAVE DSiArrowY = $(GetDataFolder(1) + "ROI" + num2str(j) + "_DSiArrowY")
		Make/O/N=2 $("ROI" + num2str(j) + "_thetaArrowX")
		WAVE thetaArrowX = $(GetDataFolder(1) + "ROI" + num2str(j) + "_thetaArrowX")
		
		//set up vector for plot with values from current ROI
		DSiArrowY[0] = 0
		DSiArrowY[1] = DSis[j] //arrow pointing from origin with length DSis[i]
		thetaArrowX[0] = thetas[j]
		thetaArrowX[1] = thetas[j] //direction that arrowy points in
		
		//point to directional DFs from current ROI
		WAVE dirDFs = $(GetDataFolder(1) + prefix + suffix + "_ROI" + num2str(j) + "_dF")
		
		//sort angle and corresponding dFs to be clockwise
		Duplicate/O angle angleSort 
		Sort angleSort, angleSort, dirDFs //sort angleSort ascending and sort dirDFs to angleSort
		
		angleSort[numpnts(angleSort)] = {angleSort[0]}
		dirDFs[numpnts(dirDFs)] = {dirDFs[0]}
		
		//Build polar plots
		graphname = "ROI" + num2str(j) + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, dirDFs, angleSort, 360)	
		WMPolarAppendTrace(graphname, DSiArrowY, thetaArrowX, 360)
		WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		//WMPolarSetManualRadiusRange(graphname, -.01, WaveMax(dirDFs))
	endfor
	
	KillWaves/Z ROIavg, ROIavg_notSmth, xpts, ypts, ROIdF //$(GetDataFolder(1) + "M_ROIMask")
end	
