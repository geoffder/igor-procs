#pragma rtGlobals=3		// Use modern global access method and strict wave access.
Function stackFixer(prefix, mode)
	String prefix
	Variable mode //0: hotSpotter naming format, 1:Original Nidaq scans
	
	String pathToMatrix
	Variable i, j, k, count, xWidth, yHeight, frames, xPos, yPos, zPos, directions, totalXshift, totalYshift
	Variable lastXshift, lastYshift
	
	//prompt inputs
	Variable firstScan, lastScan, overwrite, background, x1, x2, y1, y2
	String REFlocation
	
	directions = 8 //number of directions per trial
	
	switch(mode)
		case 0:
			Prompt REFlocation, "Folder with master REF (or blank): "
			Prompt overwrite, "Overwrite? (0/1): "
			Prompt background, "Background (sub-dendrite) levels: "
			Prompt x1, "Set left x value for alignment window: "
			Prompt x2, "Set right x value for alignment window: "
			Prompt y1, "Set first y value for alignment window: "
			Prompt y2, "Set last y value for alignment window: "
			DoPrompt "Leave coordinates blank to not use crop", REFlocation, overwrite, background, x1, x2, y1, y2
			
			Print "REFlocation: " + REFlocation
			Print "Overwrite: " + num2str(overwrite) + ", Background: " + num2str(background)
			Print  "x1: " + num2str(x1) + ", x2: "+ num2str(x2) + ", y1: " + num2str(y1) + ", y2: " + num2str(y2)
			
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
			print trials
			//To get appropriate matrix dimensions
			//first trial of first direction as example
			WAVE Matrix = $(GetDataFolder(1) + prefix + "0_" + num2str(trials[0]))
			xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
			yHeight = DimSize(Matrix, 1)
			frames = DimSize(Matrix, 2)
	
			Make/O/N=(xWidth,yHeight) avgMatrixREF, avgMatrixTEST
	
			//create reference matrix for aligning all other stacks (first direction of first trial)
			Redimension/S Matrix
			if(!overwrite && stringmatch(REFlocation,""))
				Duplicate/O Matrix $(prefix + "Fix" + "0" + "_" + num2str(trials[0])) //to fit in with aligned output naming
			endif
			Duplicate/O Matrix workMatrix //temporary matrices for spatial pre-filtering before making reference image
			Make/O/N=(xWidth,yHeight) tempMatrix
			
			if(stringmatch(REFlocation,""))
				//apply spatial filter to movie stack to prevent abberant pixels from skewing alignment
				for (zPos = 0; zPos < frames; zPos +=1) 
					for (yPos = 0; yPos < yHeight; yPos += 1)
						for(xPos = 0; xPos < xWidth; xPos += 1) 
							tempMatrix[xPos][yPos] = Matrix[xPos][yPos][zPos] //create 2D matrix on each frame of the stack
							if (tempMatrix[xPos][yPos] < background) //try to wipe background specifically
								tempMatrix[xPos][yPos] = 0
							endif
							if(tempMatrix[xPos][yPos] > 1000) //quash the rando pixels BEFORE smudging
								tempMatrix[xPos][yPos] = 1000 //now when averaged, values low enough to be wiped later
							endif 
						endfor
					endfor
					MatrixFilter/N=(7) avg tempMatrix //set prefilter type here (types: avg, median, etc)
					workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
				endfor
	
				//average filtered stack to make REF image to align TEST image (created later)
				for (yPos = 0; yPos < yHeight; yPos += 1)
					for(xPos = 0; xPos < xWidth; xPos += 1)
						MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
						avgMatrixREF[xPos][yPos] = mean(zWave) //Plug the mean intensity for this pixel in to avgMatrixREF
						if(avgMatrixREF[xPos][yPos] > 250)
							avgMatrixREF[xPos][yPos] = 250
						endif
						if (avgMatrixREF[xPos][yPos] < background/2) //try to wipe background specifically
							avgMatrixREF[xPos][yPos] = 0
						endif
					endfor
				endfor
				
				if(!(x1 == 0 && x2 == 0 && y1 == 0 && y2 == 0)) //If no coordinates are entered in to prompt, operate on entire image
					Duplicate/O/R=[x1, x2][y1, y2] avgMatrixREF tempMatrix
					Duplicate/O tempMatrix avgMatrixREF
				endif
			
			elseif(exists("root:" + REFlocation + ":" + "avgMatrixREF")) //use REFlocation to grab avgMatrixREF from a previously operated on folder
				WAVE Matrix = $("root:" + REFlocation + ":" + "avgMatrixREF")
				Duplicate/O Matrix avgMatrixREF //copy target REF into current data folder
			else
				print("Specified folder '" + REFlocation + "' does not exist or does not contain avgMatrixREF.")
				abort "Aborting... please ensure name is correct and/or that you are ready for this step."	
			endif //end of 'if no master REF entered'
			
			for (j = 0; j < numpnts(trials); j += 1)
				for (i = 0; i < directions; i += 1) 
					if (stringmatch(REFlocation, "") && j == 0 && i == 0) //if using master REF, do align first stack
						continue //skip aligning reference stack
					endif
			
					//point to the next stack to align
					WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(i) + "_" + num2str(trials[j])) 
			
					//temporary matrices for spatial pre-filtering before making test image for
					Redimension/S Matrix
					Duplicate/O Matrix workMatrix
					Make/O/N=(xWidth,yHeight) tempMatrix, avgMatrixTEST
			
					//apply spatial filter to movie stack to prevent abberant pixels from skewing alignment
					for (zPos = 0; zPos < frames; zPos +=1)
						for (yPos = 0; yPos < yHeight; yPos += 1)
							for(xPos = 0; xPos < xWidth; xPos += 1) 
								tempMatrix[xPos][yPos] = Matrix[xPos][yPos][zPos] //create 2D matrix on each frame of the stack
								if (tempMatrix[xPos][yPos] < background) //try to wipe background specifically
									tempMatrix[xPos][yPos] = 0
								endif
								if(tempMatrix[xPos][yPos] > 1000) //quash the rando pixels BEFORE smudging
									tempMatrix[xPos][yPos] = 1000 //now when averaged, values low enough to be wiped later
								endif 
							endfor
						endfor
						MatrixFilter/N=(7) avg tempMatrix //set prefilter type here (types: avg, median, etc)
						workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
					endfor
			
					//average filtered stack to make TEST image for alignment to REF image
					for (yPos = 0; yPos < yHeight; yPos += 1)
						for(xPos = 0; xPos < xWidth; xPos += 1)
							MatrixOp/O zWave = beam(workMatrix, xPos, yPos)  //Create a zWave for current pixel over time
							avgMatrixTEST[xPos][yPos] = mean(zWave) //Plug the mean intensity for this pixel in to avgMatrixTEST
							if(avgMatrixTEST[xPos][yPos] > 250)
								avgMatrixTEST[xPos][yPos] = 250
							endif
							if (avgMatrixTEST[xPos][yPos] < background/2) //try to wipe background specifically
								avgMatrixTEST[xPos][yPos] = 0
							endif
						endfor
					endfor
					
					if(!(x1 == 0 && x2 == 0 && y1 == 0 && y2 == 0)) //If no coordinates are entered in to prompt, operate on entire image
						Duplicate/O/R=[x1, x2][y1, y2] avgMatrixTEST tempMatrix
						Duplicate/O tempMatrix avgMatrixTEST
					endif
					
					Redimension/S avgMatrixREF, avgMatrixTEST, Matrix
					
					//ImageTransform/Q/O/IOFF={lastXshift,lastYshift,0} offsetImage avgMatrixTEST //testing...
					
					ImageRegistration/Q/CONV=1/TRNS={1,1,0}/ROT={0,0,0}/SKEW={0,0,0}/ZMAN testWave=avgMatrixTEST, refWave=avgMatrixREF
					
					//wave with output from ImageRegistration
					WAVE RegParams = $(GetDataFolder(1) + "W_RegParams") //RegParams[0] = dx, RegParams[1] = dy
					//WAVE RegOut = $(GetDataFolder(1) + "M_RegOut")
					
					//totalXshift = RegParams[0]
					//totalYshift = RegParams[1]
					
					//ImageRegistration/Q/CONV=1/TRNS={1,1,0}/ROT={0,0,0}/SKEW={0,0,0} testWave=RegOut, refWave=avgMatrixREF
					
					//totalXshift += RegParams[0]
					//totalYshift += RegParams[1]
					//for(k = 0; k < 10; k +=1)
					//	ImageRegistration/Q/CONV=0/TRNS={1,1,0}/ROT={0,0,0}/SKEW={0,0,0} testWave=avgMatrixTEST, refWave=avgMatrixREF
					//	WAVE RegParams = $(GetDataFolder(1) + "W_RegParams") //RegParams[0] = dx, RegParams[1] = dy
					//	ImageTransform/Q/O/IOFF={RegParams[0],RegParams[1],0} offsetImage avgMatrixTEST
					//	//accumulate the final vaules for offset
					//	totalXshift += RegParams[0]
					//	totalYshift += RegParams[1]
					//endfor
					
					ImageTransform/Q/IOFF={RegParams[0],RegParams[1],0} offsetImage Matrix //offset and overwrite original stack
					//ImageTransform/Q/IOFF={(RegParams[0]+lastXshift),(RegParams[1]+lastYshift),0} offsetImage Matrix
					//ImageTransform/Q/IOFF={totalXshift,totalYshift,0} offsetImage Matrix //offset and overwrite original stack
					//lastXshift += RegParams[0]
					//lastYshift += RegParams[1]	
					WAVE Matrix = $(GetDataFolder(1) + "M_OffsetImage")
					if(!overwrite)
						Duplicate/O Matrix $(prefix + "Fix" + num2str(i) + "_" + num2str(trials[j]))
					elseif(overwrite)
						Duplicate/O Matrix $(prefix + num2str(i) + "_" + num2str(trials[j]))
					endif
				endfor
	
			endfor
			
		break
		
		case 1:
			Prompt firstScan, "First scan: "
			Prompt lastScan, "Last scan: "
			Prompt overwrite, "Overwrite? (0/1): "
			Prompt background, "Background (sub-dendrite) levels: "
			Prompt x1, "Set left x value for alignment window: "
			Prompt x2, "Set right x value for alignment window: "
			Prompt y1, "Set bottom y value for alignment window: "
			Prompt y2, "Set top y value for alignment window: "
			DoPrompt "Leave coordinates blank to not use crop", firstScan, lastScan, overwrite, background, x1, x2, y1, y2
			
			Print "First scan: " + num2str(firstScan) + ", Last scan: " + num2str(lastScan) + ", Overwrite: " + num2str(overwrite) + ", Background: " + num2str(background)
			Print "x1: " + num2str(x1) + ", x2: "+ num2str(x2) + ", y1: " + num2str(y1) + ", y2: " + num2str(y2)
			
			String zeroPads
			
			//example full path-- root:Nidaq_Scans:d1TTX_c100_000:d1TTX_c100_000_ch2
			//prefix = d1TTX_c100
			
			//need to use this if block for building matrix paths due to Nidaq scan naming convention
			if(firstScan < 10)
				zeroPads = "00"
			elseif(firstScan < 100)
				zeroPads = "0"
			else
				zeroPads = ""
			endif
			
			//To get appropriate matrix dimensions using first trial as example
			WAVE Matrix = $("root:Nidaq_Scans:" + prefix + "_" + zeroPads + num2str(firstScan) + ":" +  prefix + "_" +zeroPads + num2str(firstScan) + "_ch2")
			xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
			yHeight = DimSize(Matrix, 1)
			frames = DimSize(Matrix, 2)
			Make/O/N=(xWidth,yHeight) avgMatrixREF
	
			//create reference matrix for aligning all other stacks (first direction of first trial)
			Redimension/S Matrix
			if(!overwrite)
				Duplicate/O Matrix $(prefix + "Fix" + "_" + num2str(firstScan))
			endif
			Duplicate/O Matrix workMatrix //temporary matrices for spatial pre-filtering before making reference image
			Make/O/N=(xWidth,yHeight) tempMatrix
	
			//apply spatial filter to movie stack to prevent abberant pixels from skewing alignment
			for (zPos = 0; zPos < frames; zPos +=1)
				for (yPos = 0; yPos < yHeight; yPos += 1)
					for(xPos = 0; xPos < xWidth; xPos += 1) 
						tempMatrix[xPos][yPos] = Matrix[xPos][yPos][zPos] //create 2D matrix on each frame of the stack
						if (tempMatrix[xPos][yPos] < background) //try to wipe background specifically
							tempMatrix[xPos][yPos] = 0
						endif
						if(tempMatrix[xPos][yPos] > 1000) //quash the rando pixels BEFORE smudging
							tempMatrix[xPos][yPos] = 1000 //now when averaged, values low enough to be wiped later
						endif 
					endfor
				endfor
				MatrixFilter/N=(7) avg tempMatrix //set prefilter type here (types: avg, median, etc)
				workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
			endfor
	
			//average filtered stack to make REF image to align TEST image (created later)
			for (yPos = 0; yPos < yHeight; yPos += 1)
				for(xPos = 0; xPos < xWidth; xPos += 1)
					MatrixOp/O zWave = beam(workMatrix, xPos, yPos)   //Create a zWave for current pixel over time
					avgMatrixREF[xPos][yPos] = mean(zWave,3,30) //Plug the mean intensity for this pixel in to avgMatrixREF (3, 30 for baseline signal)
					if(avgMatrixREF[xPos][yPos] > 200)
						avgMatrixREF[xPos][yPos] = 200
					endif
					if (avgMatrixREF[xPos][yPos] < background/2) //try to wipe background specifically
						avgMatrixREF[xPos][yPos] = 0
					endif
				endfor
			endfor
			
			if(!(x1 == 0 && x2 == 0 && y1 == 0 && y2 == 0)) //If no coordinates are entered in to prompt, operate on entire image
				Duplicate/O/R=[x1, x2][y1, y2] avgMatrixREF tempMatrix
				Duplicate/O tempMatrix avgMatrixREF
			endif
			
			for (j = firstScan+1; j <= lastScan; j += 1)
				
				if(j < 10)
					zeroPads = "00"
				elseif(j < 100)
					zeroPads = "0"
				else
					zeroPads = ""
				endif	
				
				//point to the next stack to align
				WAVE Matrix = $("root:Nidaq_Scans:" + prefix + "_" + zeroPads + num2str(j) + ":" +  prefix +"_" + zeroPads + num2str(j) + "_ch2") 
				
				//temporary matrices for spatial pre-filtering before making test image for
				Redimension/S Matrix
				Duplicate/O Matrix workMatrix
				Make/O/N=(xWidth,yHeight) tempMatrix, avgMatrixTEST
		
				//apply spatial filter to movie stack to prevent abberant pixels from skewing alignment
				for (zPos = 0; zPos < frames; zPos +=1)
					for (yPos = 0; yPos < yHeight; yPos += 1)
						for(xPos = 0; xPos < xWidth; xPos += 1) 
							tempMatrix[xPos][yPos] = Matrix[xPos][yPos][zPos] //create 2D matrix on each frame of the stack
							if (tempMatrix[xPos][yPos] < background) //try to wipe background specifically
								tempMatrix[xPos][yPos] = 0
							endif 
							if(tempMatrix[xPos][yPos] > 1000) //quash the rando pixels BEFORE smudging 
								tempMatrix[xPos][yPos] = 1000 //now when averaged, values low enough to be wiped later
							endif 
						endfor
					endfor
					MatrixFilter/N=(7) avg tempMatrix //set prefilter type here (types: avg, median, etc)
					workMatrix[][][zPos] = tempMatrix[p][q] //then spatially filter it and plug it in to the corresponding frame
				endfor
			
				//average filtered stack to make TEST image for alignment to REF image
				for (yPos = 0; yPos < yHeight; yPos += 1)
					for(xPos = 0; xPos < xWidth; xPos += 1)
						MatrixOp/O zWave = beam(workMatrix, xPos, yPos)  //Create a zWave for current pixel over time
						avgMatrixTEST[xPos][yPos] = mean(zWave,3,30) //Plug the mean intensity for this pixel in to avgMatrixTEST
						if(avgMatrixTEST[xPos][yPos] > 200)
							avgMatrixTEST[xPos][yPos] = 200
						endif
						if (avgMatrixTEST[xPos][yPos] < background/2) //try to wipe background specifically
							avgMatrixTEST[xPos][yPos] = 0
						endif
					endfor
				endfor
							
				if(!(x1 == 0 && x2 == 0 && y1 == 0 && y2 == 0)) //If no coordinates are entered in to prompt, operate on entire image
					Duplicate/O/R=[x1, x2][y1, y2] avgMatrixTEST tempMatrix
					Duplicate/O tempMatrix avgMatrixTEST
				endif
			
				Redimension/S avgMatrixREF, avgMatrixTEST //ensure matrices are single precision FP32 for ImageRegistration
				ImageRegistration/Q/TRNS={1,1,0}/ROT={0,0,0}/SKEW={0,0,0} testWave=avgMatrixTEST, refWave=avgMatrixREF
				//wave with output from ImageRegistration
				WAVE RegParams = $(GetDataFolder(1) + "W_RegParams") //RegParams[0] = dx, RegParams[1] = dy
				
				ImageTransform/Q/IOFF={RegParams[0],RegParams[1],0} offsetImage Matrix //offset and overwrite original stack
				Redimension/W Matrix //return scan to original format, I16 (16bit integer)
				WAVE Matrix = $(GetDataFolder(1) + "M_OffsetImage")
				if(overwrite)
					//overwrite the original with transformed matrix output
					Redimension/W Matrix //make M_OffsetImage output I16 before overwriting original 
					Duplicate/O Matrix $("root:Nidaq_Scans:" + prefix + "_" + zeroPads + num2str(j) + ":" +  prefix +"_" + zeroPads + num2str(j) + "_ch2") 
				else
					Duplicate/O Matrix $(prefix + "Fix" + num2str(i) + "_" + num2str(j)) 
				endif
				
			endfor
			
			WAVE Matrix = $("root:Nidaq_Scans:" + prefix + "_" + zeroPads + num2str(firstScan) + ":" +  prefix + "_" +zeroPads + num2str(firstScan) + "_ch2")
			Redimension/W Matrix //return first (reference) scan to original format, I16 (16bit integer)
		break
		
		default:
			print("Please select a mode: ")
			print("0: hotSpotter format; 1: Nidaq scans")
		break
	
	endswitch
	
	KillWaves/Z workMatrix, tempMatrix, zWave, M_RegOut, M_RegMaskOut, W_RegParams, M_OffsetImage		
end

//Function/D Median(w) // Returns median value of wave w
//	Wave w
//	
//	Variable result
//	
//	Duplicate/O w, temp // Make a clone of wave
//	Sort temp, temp // Sort clone
//	result = temp[numpnts(temp)/2] //get middle point of sorted wave (median)
//	KillWaves temp // Kill clone
//	
//	return result
//end