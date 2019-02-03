#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function scanRenamer(channel, prefixOld, prefixNew, firstScan, numTrials)
	String prefixOld, prefixNew
	Variable channel, firstScan, numTrials
	
	Variable i, j, currentScan, zeroOffset, xDelta, yDelta
	String zeroPads
	
	
	//set to 1 or 0 to turn on offset zeroing (scaling origin of x and y) 
	zeroOffset = 1
	
	currentScan = firstScan
	
	for(i = 1; i <= numTrials; i += 1)
		
		for(j = 0; j < 8; j += 1)
			if(currentScan < 10)
				zeroPads = "00"
			elseif(currentScan < 100)
				zeroPads = "0"
			else
				zeroPads = ""
			endif
			
			WAVE Matrix = $("root:Nidaq_Scans:" + prefixOld + "_" + zeroPads + num2str(currentScan) + ":" +  prefixOld +"_" + zeroPads + num2str(currentScan) + "_ch" + num2str(channel))
			Duplicate Matrix $(prefixNew + num2str(j) + "_" + num2str(i)) 
			
			if(zeroOffset)
				WAVE Matrix = $(prefixNew + num2str(j) + "_" + num2str(i))
				
				xDelta = DimDelta(Matrix, 0)
				yDelta = DimDelta(Matrix, 1)
				
				SetScale/P x 0, xDelta, "m", Matrix
				SetScale/P y 0, yDelta, "m", Matrix
			endif
			
			currentScan += 1
		endfor
		
	endfor
	
End
