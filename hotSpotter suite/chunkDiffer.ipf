#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function chunkDiffer(prefix, diffChunkSize, maxDist, numLines)
	String prefix
	Variable diffChunkSize, maxDist, numLines
	
	String pathToMatrix
	Variable h, i, j, k, count, first, x1
	
	Make/O/N=1 trials
	count = 0
	for(i = 1; i < 12; i += 1)
		pathToMatrix = GetDataFolder(1) + prefix + "0_" + num2str(i)
		if(exists(pathToMatrix))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	for(j = 1; j <= numLines; j += 1)
		//Prepare trial average output
		Make/O/N=(maxDist) $("chunkDiffAVG_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		WAVE chunkDiffAVG_trAv = $("chunkDiffAVG_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		Make/O/N=(maxDist) $("chunkDiffSEM_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		WAVE chunkDiffSEM_trAv = $("chunkDiffSEM_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		Make/O/N=(maxDist) $("chunkDiffSD_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		WAVE chunkDiffSD_trAv = $("chunkDiffSD_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		Make/O/N=(maxDist) $("chunkDiffMED_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
		WAVE chunkDiffMED_trAv = $("chunkDiffMED_Sz" + num2str(diffChunkSize) + "_trAv_ln" + num2str(j))
			
		for(k = 0; k < numpnts(trials); k+= 1)
			
			//Point to chunkTheta wave of chosen size and current trial
			WAVE chunkThetas = $(GetDataFolder(1) + "thetas_chnkSz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			
			//Prepare trial by trial output
			Make/O/N=(maxDist) $("chunkDiffAVG_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			WAVE chunkDiffAVG = $("chunkDiffAVG_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			Make/O/N=(maxDist) $("chunkDiffSEM_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			WAVE chunkDiffSEM = $("chunkDiffSEM_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			Make/O/N=(maxDist) $("chunkDiffSD_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			WAVE chunkDiffSD = $("chunkDiffSD_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			Make/O/N=(maxDist) $("chunkDiffMED_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			WAVE chunkDiffMED = $("chunkDiffMED_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j))
			
			SetScale/P x diffChunkSize, diffChunkSize,"distance (pts)", chunkDiffAVG, chunkDiffSEM, chunkDiffSD,chunkDiffMED
			
			first = 1
			for(h = 1; h <= maxDist; h += 1)
				Make/O/N=1 $("chunkDiff_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j) + "_d" + num2str(h))
				WAVE chunkDiffW = $("chunkDiff_Sz" + num2str(diffChunkSize) + "_" + num2str(trials[k]) + "_ln" + num2str(j) + "_d" + num2str(h))
				 
				for(x1 = 0; x1 < numpnts(chunkThetas)-1; x1 += 1)
					if((x1+h) < numpnts(chunkThetas))
						
						if (first == 1)
							chunkDiffW[0] = {abs(chunkThetas[x1+h] - chunkThetas[x1])}
							first = 0
						else
							chunkDiffW[numpnts(chunkDiffW)] = {abs(chunkThetas[x1+h] - chunkThetas[x1])}
						endif
						
					endif
				endfor
				
				WaveStats/Q/W chunkDiffW 
				WAVE stats = M_WaveStats // W flag in WaveStats puts results in M_WaveStats
				chunkDiffAVG[h-1] = stats[3] //mean
				chunkDiffSEM[h-1] = stats[26] //standard error
				chunkDiffSD[h-1] = stats[4] //standard deviation
		
				Duplicate/O chunkDiffW chunkDiffWmed
				Sort chunkDiffWmed chunkDiffWmed
				chunkDiffMED[h-1] = chunkDiffWmed[(numpnts(chunkDiffWmed)-1)/2] //median
				KillWaves chunkDiffWmed
				
				KillWaves chunkDiffW //comment out if raw diff waves needed
			
			endfor
			
			//Accumulate trial data for averaging
			chunkDiffAVG_trAv += chunkDiffAVG
			chunkDiffSEM_trAv += chunkDiffSEM
			chunkDiffSD_trAv += chunkDiffSD
			chunkDiffMED_trAv += chunkDiffMED
		
		endfor // end of trials loop (k)
		
		//Finalize trial averages
		chunkDiffAVG_trAv /= numpnts(trials)
		chunkDiffSEM_trAv /= numpnts(trials)
		chunkDiffSD_trAv /= numpnts(trials)
		chunkDiffMED_trAv /= numpnts(trials)
	
	endfor //end of lines loop (j)
	
	KillWaves/Z M_WaveStats
End