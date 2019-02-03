#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <New Polar Graphs>

//Before running, ensure graph you would like to draw ROIs on is topmost/active
//Enter command (ex: polarPlotter("raw_dir","_AVG",0,3) for 3 ROIs).
//You'll be prompted to draw a dark region, then a number of signal ROIs you specified in the args.
//Program will spit out ROIavgs (z projections), dFs, thetas, DSis and whole polar plots.
//For now, run multiple times using suffix (e.g. "_1") to deal with multiple trials

Function miniPlotter(prefix, importROIs, numROIs)
	String prefix
	Variable importROIs, numROIs
	
	Variable i, j, k, frames, darkAvg, bsln, count
	String roiName, winTitle, graphname, ROIfolder, pathToMatrix
	
	Variable smthSize = 7 //smoothing function bin size
	Variable start = 6 //frame to start calculations from (skip laser response)
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 12; i += 1)
		pathToMatrix = prefix + num2str(i)
		if(exists(pathToMatrix))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
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
		
		for(i = 0; i < numpnts(trials); i += 1) //for each trial
			
			WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i)) //matrix mask is applied to
			
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
			
			if(j == 0) //first ROI is dark spot
				darkAvg = mean(ROIavg) //calculate background fluorescence
				Duplicate/O ROIavg $(prefix + num2str(i) + "_ROI" + num2str(j))
			else
				ROIavg -= darkAvg //subtract background fluorescence
				//take care of unsmoothed ROI trace
				Duplicate/O ROIavg ROIavg_notSmth
				Smooth/S=2 smthSize, ROIavg //Savitzky-Golay
				bsln = mean(ROIavg,start,numpnts(ROIavg)-1) //baseline of entire time series	
				
				//Now suppressing output of non-normalized waveforms, uncomment or add conditional
				//if needed in the future
				//Save with appropriate names
				//Duplicate/O ROIavg $(prefix + num2str(i) + "_ROI" + num2str(j))
				//Duplicate/O ROIavg $(prefix + num2str(i) + "_ROI" + num2str(j) + "_smth")
			
				//dF/F0 normalization
				ROIavg = (ROIavg - bsln)/bsln
				ROIavg_notSmth = (ROIavg_notSmth - bsln)/bsln
				//saving output
				Duplicate/O ROIavg_notSmth $(prefix + num2str(i) + "_ROI" + num2str(j) + "_dF")
				Duplicate/O ROIavg $(prefix + num2str(i) + "_ROI" + num2str(j) + "_dFsm")
			endif
		
		endfor
		
	endfor
	
	KillWaves/Z ROIavg
end	

