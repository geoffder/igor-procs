#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//runs for all trials automatically
//expects naming conventionto be (prefix + trialNum) e.g. test0, test1 ...
//make sure the graph you want to draw on is the top window
//outputs space-time dF plots, where x is space, and y is time
Function linePlotterMini(prefix, importLines, numLines)
	String prefix
	Variable importLines, numLines
	
	Variable i, j, k, z, count, xDelta, yDelta, avgDelta, start
	Variable lineDivAVG, ptsPerPix, pt, lineWidth, xPos, yPos, bsln
	Variable darkROIavg, frames
	
	String pathToMatrix
	String lineName, winTitle, graphname, lineFolder
	
	lineWidth = 5 //width of pixels (centred on line) to avg over 
	start = 6 //first frame to calculate baseline avg from
	
	//find how many trials there are and create 
	//an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 12; i += 1)
		pathToMatrix = prefix + num2str(i)
		if(exists(pathToMatrix))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	WAVE Matrix = $(prefix + num2str(trials[0]))

	xDelta = DimDelta(Matrix, 0)
	yDelta = DimDelta(Matrix, 1)
	frames = DimSize(Matrix, 2)
	
	avgDelta = (xDelta+yDelta)/2
	
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
			
	for(j = 1; j <= numLines; j+=1)
		WAVE lineY = $(GetDataFolder(1) + "line" + num2str(j) + "_y")
		WAVE lineX = $(GetDataFolder(1) + "line" + num2str(j) + "_x")
		
		lineDivAVG = 0
		for(pt = 0; pt < numpnts(lineX); pt += 1)
			lineDivAVG += sqrt((lineX[pt+1] - lineX[pt])^2 + (lineY[pt+1] - lineY[pt])^2)
		endfor 
		lineDivAVG /= numpnts(lineX)
		
		ptsPerPix = round(avgDelta/lineDivAVG)
		Print "------ Line " + num2str(j) + " ------"
		Print "Points per pixel = " + num2str(ptsPerPix)
		
		for(k = 0; k < numpnts(trials); k += 1) //will still run once if there are no trials (e.g. running on AVGs instead)
			
			ImageLineProfile/P=-2 xWave=lineX, yWave=lineY, srcwave=Matrix, width=lineWidth
			WAVE profile = $("M_ImageLineProfile")
			SetScale/P x 0, avgDelta/ptsPerPix,"m", profile
			 
			for(xPos = 0; xPos < DimSize(profile, 0); xPos += 1)
				MatrixOP/O timeWave = row(profile, xPos)
				Redimension/N=(DimSize(profile, 1))/D timeWave
				timeWave -= darkROIavg
				bsln = mean(timeWave, start, numpnts(timeWave))
				timeWave = (timeWave - bsln)/bsln //convert to dF
				Redimension/N=(1,numpnts(timeWave))/D timeWave
				profile[xPos][] = timeWave[0][q]
			endfor
			Duplicate/O profile $("spTmDF" + num2str(j) + "_" + prefix + num2str(k))
		
		endfor //trial loop (k)
		
	endfor //line loop (j)
	
	KillWaves/z profile, timeWave
End