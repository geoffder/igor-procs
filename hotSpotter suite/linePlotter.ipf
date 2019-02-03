#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//runs for all trials automatically. Alternately,a suffix like "AVG" can be specified
//optional parameters such as [,suffix] are activated in args like so: (... , suffix="AVG") 
Function linePlotter(prefix, importLines, numLines [,suffix])
	String prefix, suffix
	Variable importLines, numLines
	
	Variable h,i, j, k, count, xDelta, yDelta, avgDelta,x1, first, maxDist, thetaMode, DSiMode
	Variable lineDivAVG, ptsPerPix, pt, lineWidth, is3D, xPos, yPos, bsln
	
	String pathToMatrix, suffixStr, dirLabel
	String lineName, winTitle, graphname, lineFolder
	
	lineWidth= 5
	maxDist = 30
	
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
	
	//engage theta mode if prefix contains the word 'theta'. Name output accordingly and do not deal in directions
	//only execute below i loop once if input is already multi-directional (one containing all directions)
	if(StringMatch(prefix,"*theta*"))
		thetaMode = 1
		DSiMode = 0
	elseif(StringMatch(prefix,"*DSi*"))
		thetaMode = 0
		DSiMode = 1
	else
		thetaMode = 0
		DSiMode = 0
	endif
	
	if(thetaMode||DSiMode)
		WAVE Matrix = $(GetDataFolder(1) + prefix + suffix)
	elseif(ParamIsDefault(suffix))
		WAVE Matrix = $(GetDataFolder(1) + prefix + "0_" + num2str(trials[0]))
	else
		WAVE Matrix = $(GetDataFolder(1) + prefix + "0_" + suffix)
	endif
	xDelta = DimDelta(Matrix, 0)
	yDelta = DimDelta(Matrix, 1)
	
	if(DimSize(Matrix, 2) > 1)
		is3D = 1
	else
		is3D = 0
	endif
	
	avgDelta = (xDelta+yDelta)/2
	
	if(!importLines)
		for(j = 1; j <= numLines; j += 1) 
			lineName = "line" + num2str(j) 
			
			GetWindow kwTopWin wtitle //get window title for active graph
			winTitle = S_value
			
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
		for(pt = 0; pt <= 30; pt += 1)
			lineDivAVG += sqrt((lineX[pt+1] - lineX[pt])^2 + (lineY[pt+1] - lineY[pt])^2)
		endfor 
		lineDivAVG /= 30
		
		ptsPerPix = round(avgDelta/lineDivAVG)
		Print "------ Line " + num2str(j) + " ------"
		Print "Points per pixel = " + num2str(ptsPerPix)
		
		for(i = 0; i < 8; i += 1) //for each direction
			
			if(thetaMode)
				dirLabel = "T"
			elseif(DSiMode)
				dirLabel = "D"
			else
				dirLabel = num2str(i)
			endif
			
			for(k = 0; k < numpnts(trials); k += 1) //will still run once if there are no trials (e.g. running on AVGs instead)
				//set suffix to trial number by default if no suffix (e.g. "AVG") specified
				if(ParamIsDefault(suffix))
					suffixStr = num2str(trials[k])
				else
					suffixStr = suffix
				endif		
				
				//conditional to allow using theta matrix as src
				if(thetaMode||DSiMode)
					WAVE Matrix = $(GetDataFolder(1) + prefix + suffixStr) //e.g. 'relThetaMatrix' + 'AVGraw'
				else
					WAVE Matrix = $(GetDataFolder(1) + prefix + dirLabel + "_" + suffixStr)
				endif
				
				if(!is3D)
					ImageLineProfile xWave=lineX, yWave=lineY, srcwave=Matrix, width=lineWidth
					WAVE profile = $(GetDataFolder(1) + "W_ImageLineProfile")
					SetScale/P x 0, avgDelta/ptsPerPix,"m", profile 
				elseif(is3D)
					ImageLineProfile/P=-2 xWave=lineX, yWave=lineY, srcwave=Matrix, width=lineWidth
					WAVE profile = $(GetDataFolder(1) + "M_ImageLineProfile")
					SetScale/P x 0, avgDelta/ptsPerPix,"m", profile 
				endif
				
				if(!is3D)
					//############### difference calculator ###############
					//Prepare output
					Make/O/N=(maxDist) $("lineDiffAVG" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					WAVE lineDiffAVG = $("lineDiffAVG" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					Make/O/N=(maxDist) $("lineDiffSEM" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					WAVE lineDiffSEM = $("lineDiffSEM" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					Make/O/N=(maxDist) $("lineDiffSD" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					WAVE lineDiffSD = $("lineDiffSD" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					Make/O/N=(maxDist) $("lineDiffMED" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					WAVE lineDiffMED = $("lineDiffMED" + dirLabel + "_" + suffixStr + "_line" + num2str(j))
					
					SetScale/P x .3333,.3333,"μm", lineDiffAVG, lineDiffSEM, lineDiffSD, lineDiffMED
					
					first = 1
					for(h = 1; h <= maxDist; h += 1)
						Make/O/N=1 $("lineDiff" + dirLabel + "_" + suffixStr + "_line" + num2str(j) + "_d" + num2str(h))
						WAVE lineDiffW = $("lineDiff" + dirLabel + "_" + suffixStr + "_line" + num2str(j) + "_d" + num2str(h))
						 
						for(x1 = 0; x1 < numpnts(profile)-1; x1 += ptsPerPix)
							if((x1+h*ptsPerPix) < numpnts(profile))
								
								if (first == 1)
									lineDiffW[0] = {abs(profile[x1+h*ptsPerPix] - profile[x1])}
									first = 0
								else
									lineDiffW[numpnts(lineDiffW)] = {abs(profile[x1+h*ptsPerPix] - profile[x1])}
								endif
								
							endif
						
						endfor
						WaveStats/Q/W lineDiffW 
						WAVE stats = M_WaveStats // W flag in WaveStats puts results in M_WaveStats
						lineDiffAVG[h-1] = stats[3] //mean
						lineDiffSEM[h-1] = stats[26] //standard error
						lineDiffSD[h-1] = stats[4] //standard deviation
				
						Duplicate/O lineDiffW lineDiffWmed
						Sort lineDiffWmed lineDiffWmed
						lineDiffMED[h-1] = lineDiffWmed[(numpnts(lineDiffWmed)-1)/2] //median
						KillWaves lineDiffWmed, stats
						
						KillWaves lineDiffW //comment out if raw diff waves needed
					endfor //end of dist loop
					//############### end of difference calculator ###############
					
					Duplicate/O profile $(GetDataFolder(1) + "prf" + num2str(j) + "_" + prefix + dirLabel + "_" + suffixStr)
				elseif(is3D)
					for(xPos = 0; xPos < DimSize(profile, 0); xPos += 1)
						MatrixOP/O timeWave = row(profile, xPos)
						Redimension/N=(DimSize(profile, 1))/D timeWave
						bsln = mean(timeWave,5,15)
						timeWave = (timeWave - bsln)/bsln
						Redimension/N=(1,numpnts(timeWave))/D timeWave
						profile[xPos][] = timeWave[0][q]
					endfor
					Duplicate/O profile $(GetDataFolder(1) + "spTmDF" + num2str(j) + "_" + prefix + dirLabel + "_" + suffixStr)
				endif
				
			endfor
			
			if(thetaMode||DSiMode)
				break //stop after first loop 
			endif
			
		endfor //dir loop (i)
		
	endfor
	
	KillWaves profile
End
