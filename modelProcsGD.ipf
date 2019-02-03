#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <New Polar Graphs>

Function crfLoader(numTrials, basePath, areaMode, manyDirs)
	Variable numTrials, areaMode, manyDirs
	String basePath
	//"C:\Users\Geoff\Desktop\pref2_null17p8_nf5\pref2_null17p8_nf5_sens"
	
	Variable i, j
	String colInfo, path
	
	if(!manyDirs)
		Make/O dirs = {0, 45, 90, 135, 180}
	else
		Make/O dirs = {0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180}
	endif
	
	for(i = 0; i < numTrials; i += 1)
		for(j = 0; j < numpnts(dirs); j += 1)
			if(!areaMode)
				path = basePath + "\\t" + num2str(i) + "\dir_" + num2str(dirs[j]) + "\CRF.dat"
			else
				path = basePath + "\\t" + num2str(i) + "\dir_" + num2str(dirs[j]) + "\areaCRF.dat"
			endif
			colInfo = "C=1,N=APV"+num2str(j)+"_"+num2str(i)+";C=1,N=CTRL"+num2str(j)+"_"+num2str(i)+";"
			LoadWave/O/Q/D/H/G/K=1/ENCG={1,8}/B=(colInfo)/A path
		endfor
	endfor
	
	KillWaves/Z dirs
End

Function peakAreaLoader(name, numTrials, basePath, manyDirs)
	String name, basePath
	Variable numTrials, manyDirs
	//name ~= gNMDA / iNMDA / iEXC / etc
	//basePath ~= "C:\Users\Geoff\Desktop\modelFigure\pref2_null17p8_nf5\sens"
	
	Variable i, j
	String colInfo, path
	
	if(!manyDirs)
		Make/O dirs = {0, 45, 90, 135, 180}
	else
		Make/O dirs = {0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180}
	endif
	
	for(i = 0; i < numTrials; i += 1)
		for(j = 0; j < numpnts(dirs); j += 1)
			//peak
			path = basePath+"\\t"+num2str(i)+"\dir_"+num2str(dirs[j])+"\\"+name+"peak.dat"
			colInfo = "C=1,N="+name+"peak"+num2str(j)+"_"+num2str(i)+";"
			LoadWave/O/Q/D/H/G/K=1/ENCG={1,4}/B=(colInfo)/A path
			//area
			path = basePath+"\\t"+num2str(i)+"\dir_"+num2str(dirs[j])+"\\"+name+"area.dat"
			colInfo = "C=1,N="+name+"area"+num2str(j)+"_"+num2str(i)+";"
			LoadWave/O/Q/D/H/G/K=1/ENCG={1,8}/B=(colInfo)/A path
		endfor
	endfor
	
	KillWaves/Z dirs
End

Function ivLoader(name, basePath, neg, pos, step)
	String name, basePath
	Variable neg, pos, step
	//name ~= NMDA / Kv / etc
	//basePath ~= "C:\Users\Geoff\Desktop\contrast2comp_Kv_stuff\IV"
	//neg = negative potential, pos = positive potential
	//step = voltage step size
	
	Variable i
	String colInfo, path
	
	for(i = neg; i <= pos; i += step)
		//voltage-clamp current recording
		path = basePath+"\iv_"+name+"_"+num2str(i)+".dat"
		if(i>=0) //cannot have dash (-) character in name
			colInfo = "C=1,N=iv_"+name+"_p"+num2str(i)+";"
		else
			colInfo = "C=1,N=iv_"+name+"_n"+num2str(abs(i))+";"
		endif
		LoadWave/O/Q/D/H/G/K=1/ENCG={1,4}/B=(colInfo)/A path
		//direct channel current
		path = basePath+"\i_"+name+"_"+num2str(i)+".dat"
		if(i>=0)
			colInfo = "C=1,N=i_"+name+"_p"+num2str(i)+";"
		else
			colInfo = "C=1,N=i_"+name+"_n"+num2str(abs(i))+";"
		endif
		LoadWave/O/Q/D/H/G/K=1/ENCG={1,4}/B=(colInfo)/A path
		//direct channel conductance
		path = basePath+"\g_"+name+"_"+num2str(i)+".dat"
		if(i>=0)
			colInfo = "C=1,N=g_"+name+"_p"+num2str(i)+";"
		else
			colInfo = "C=1,N=g_"+name+"_n"+num2str(abs(i))+";"
		endif
		LoadWave/O/Q/D/H/G/K=1/ENCG={1,4}/B=(colInfo)/A path
	endfor
End

Function dirSlicer(condition,name,contrastName,contrastPt,suffix)
	String condition, name, contrastName, suffix
	Variable contrastPt
	
	Variable i, count
	String pathToWave

	//find how many directions there are and create an array
	//to loop over for this function
	Make/O/N=1 dirs
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = name + num2str(i) + suffix
		if(exists(pathToWave))
			dirs[count] = {i}
			count += 1
		endif
	endfor
	
	Make/o/n=(numpnts(dirs)) $(condition+name+contrastName+suffix)
	WAVE dirSlice = $(condition+name+contrastName+suffix)
	
	for(i = 0; i < numpnts(dirs); i += 1)
		WAVE currDir = $(name + num2str(i) + suffix)
		dirSlice[i] = currDir[contrastPt]
	endfor
End

Function contraDSi()
	Variable i, j, numContrasts, count, xsum, ysum
	
	Make/O dirLabel = {0, 45, 90, 135, 180}
	Make/O circle = {0, 45, 90, 135, 180, 225, 270, 315}
	
	WAVE APV = $(GetDataFolder(1) + "APV0_avg")
	numContrasts = numpnts(APV)
	Make/O/N=(numContrasts) DSiAPV, DSiCTRL
	Make/O/N=(numpnts(circle)) dirAPV, dirCTRL, xpts, ypts
	
	for(i = 0; i < numContrasts; i += 1)
		for(j = 0; j < numpnts(dirLabel); j += 1)
			WAVE APV = $(GetDataFolder(1) + "APV" + num2str(j) + "_avg")
			WAVE CTRL = $(GetDataFolder(1) + "CTRL" + num2str(j) + "_avg")
			dirAPV[j] = APV[i]
			dirCTRL[j] = CTRL[i]
		endfor
		count = 3
		for(j = 5; j < numpnts(circle); j +=1)
			dirAPV[j] = dirAPV[count]
			dirCTRL[j] = dirCTRL[count]
			count -= 1
		endfor
		xpts = dirAPV*cos(circle*pi/180)
		ypts = dirAPV*sin(circle*pi/180)
		xsum = sum(xpts)
		ysum = sum(ypts)
		if(sum(dirAPV)>0)
			DSiAPV[i] = sqrt(xsum^2+ysum^2)/sum(dirAPV)
		else
			DSiAPV[i] = 0
		endif
		xpts = dirCTRL*cos(circle*pi/180)
		ypts = dirCTRL*sin(circle*pi/180)
		xsum = sum(xpts)
		ysum = sum(ypts)
		if(sum(dirCTRL)>0)
			DSiCTRL[i] = sqrt(xsum^2+ysum^2)/sum(dirCTRL)
		else
			DSiCTRL[i] = 0
		endif
	endfor
	
	KillWaves/Z dirAPV, dirCTRL, xpts, ypts
End

Function averager(prefix)
	String prefix
	
	Variable count, i, j
	String pathToWave
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = prefix + "0_" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	//same for direction
	Make/O/N=1 dirs
	count = 0
	for(i = 0; i < 20; i += 1)
		pathToWave = prefix + num2str(i) + "_0"
		if(exists(pathToWave))
			dirs[count] = {i}
			count += 1
		endif
	endfor
	
	for(i = 0; i < numpnts(dirs); i += 1)
		for(j = 0; j < numpnts(trials); j += 1)
			WAVE curr = $(prefix + num2str(i) + "_" + num2str(j))
			if(j == 0)
				Duplicate/O curr $(prefix + num2str(i) + "_AVG")
				WAVE avg = $(prefix + num2str(i) + "_AVG")
			else
				avg += curr
			endif	
		endfor
		avg /= numpnts(trials)
	endfor
	
	KillWaves/Z trials, dirs
End

Function sampleFixer(name)
	String name
	
	WAVE target = $(name)
	Duplicate/O/R=[0,148] target front
	Duplicate/O/R=[149,*] target back
	
	Resample/UP=10 back
	Concatenate/O {front,back}, $(name + "_RS")
	WAVE new = $(name + "_RS")
	SetScale/P x 0.34,0.34,"", new
	
	Killwaves front, back
End

Function normalizer(prefix1, prefix2, suffix)
	String prefix1, prefix2, suffix //suffix = "_AVG"
	//prefix1 indicates the wave set to be normalized to
	
	Variable i, count, refMax
	String pathToWave
	
	//check number of directions
	Make/O/N=1 dirs
	count = 0
	for(i = 0; i < 20; i += 1)
		pathToWave = prefix1 + num2str(i) + suffix
		if(exists(pathToWave))
			dirs[count] = {i}
			count += 1
		endif
	endfor
	
	WAVE ref = $(prefix1 + "0" + suffix)
	refMax = wavemax(ref)
	for(i = 0; i < numpnts(dirs); i += 1)
		WAVE curr = $(prefix1 + num2str(i) + suffix)
		Duplicate/O curr $(prefix1 + num2str(i) + "_NORM")
		WAVE normd = $(prefix1 + num2str(i) + "_NORM")
		normd /= refMax
	endfor
	
	for(i = 0; i < numpnts(dirs); i += 1)
		WAVE curr = $(prefix2 + num2str(i) + suffix)
		Duplicate/O curr $(prefix2 + num2str(i) + "_NORM")
		WAVE normd = $(prefix2 + num2str(i) + "_NORM")
		normd /= refMax
	endfor
	
	KillWaves/Z dirs
End

Function dirSpkLoader(name, numRhos, basePath)
	String name, basePath
	Variable numRhos
	//name ~= dirSpks
	//basePath ~= "C:\Users\Geoff\Desktop\seed10k_noReset\"
	
	Variable i, j
	String colInfo, path, type, fold, rhoStr
	
	for(i = 0; i < 3; i += 1)//sp, tm, sptm
		if (!i)
			fold = "sp"
			type = "sp" 
		elseif (i == 1)
			fold = "\tm"
			type = "tm"
		else
			fold = "sptm"
			type = "sptm"
		endif
		for(j = 0; j < numRhos; j += 1)
			if (!j)
				rhoStr = "r0" 
			else
				rhoStr = "r0." + num2str(j)
			endif
			path = basePath+"\\"+type+"\\"+rhoStr+"\\"+name+".dat"
			colInfo = "C=1,N="+type+"_"+name+"_r"+num2str(j)+";"
			LoadWave/O/Q/M/D/H/G/K=1/ENCG={1,4}/B=(colInfo)/A path
		endfor
	endfor
End

Function modelPolar(type, name, rho, displayTrials)
	String type, name
	Variable rho, displayTrials
	
	Variable i, j, count, trials
	Variable xsum, ysum, radius, theta
	String pathToMatrix, graphname
	
	Make/O dirs = {225, 270, 315, 0, 45, 90,135, 180}
	
	//To get appropriate matrix dimensions
	WAVE Matrix = $(type+ "_" + name + "_r0") //first trial of first direction as example
	trials = DimSize(Matrix, 1) //returns the size of each dimension (x, y, z) 
	
	Make/O/N=(numpnts(dirs)) xpts, ypts
	
	WAVE Matrix = $(type + "_" + name + "_r" + num2str(rho))
	Make/O/N=(trials) $(type+"_thetas_r"+num2str(rho)), $(type+"_DSis_r"+num2str(rho))
	Make/O/N=(2,trials) $(type+"_DSiArrows_r" + num2str(rho)), $(type+"_thetaArrows_r" + num2str(rho))
	WAVE thetas = $(type+"_thetas_r"+num2str(rho))
	WAVE DSis = $(type+"_DSis_r"+num2str(rho))
	WAVE thetaArrows = $(type+"_thetaArrows_r" + num2str(rho))
	WAVE DSiArrows = $(type+"_DSiArrows_r" + num2str(rho))
	for(j = 0; j < trials; j += 1)
		MatrixOp/O col = col(Matrix, j)
		xpts = col*cos(dirs*pi/180) //convert each value to x and y  coordinates
		ypts = col*sin(dirs*pi/180)
		xsum = sum(xpts) //sum across all 8 directions
		ysum = sum(ypts) 
		radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
		theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
		if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
			thetas[j] = 360 + theta
		else
			thetas[j] = theta
		endif
		thetaArrows[0][j] = theta
		thetaArrows[1][j] = theta
		DSis[j] = radius/sum(col)
		if(displayTrials)
			DSiArrows[1] =  DSis[0] * 6//wavemax(Matrix)
		else
			DSiArrows[1] =  DSis[0] * 4
	endif
	endfor	
	
	//SEM of DSis
	//Wavestats/Q DSis
	//print "DSi SEM: " + num2str(V_sem)
	//Wavestats/Q thetas //can't do this on 360 degrees thetas (need +/-)
	//print "theta SEM: " + num2str(V_sem)
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~AVG BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Make/O/N=1 $(type+"_avgtheta_r"+num2str(rho)), $(type+"_avgDSi_r"+num2str(rho))
	WAVE avgTheta = $(type+"_avgtheta_r"+num2str(rho))
	WAVE avgDSi = $(type+"_avgDSi_r"+num2str(rho))
	Make/O/N=2 $(type+"_avgthetaArrow_r"+num2str(rho)), $(type+"_avgDSiArrow_r"+num2str(rho))
	WAVE avgThetaArrow = $(type+"_avgthetaArrow_r"+num2str(rho))
	WAVE avgDSiArrow = $(type+"_avgDSiArrow_r"+num2str(rho))
	Make/O/N=(numpnts(dirs)) $(type+"_dirAvg_r"+num2str(rho)), $(type+"_dirSDev_r"+num2str(rho))
	Make/O/N=(numpnts(dirs)) $(type+"_dirSEM_r"+num2str(rho))
	WAVE dirAvg = $(type+"_dirAvg_r"+num2str(rho))
	WAVE dirSDev = $(type+"_dirSDev_r"+num2str(rho))
	WAVE dirSEM = $(type+"_dirSEM_r"+num2str(rho))
	for(j = 0; j < trials; j += 1)
		MatrixOp/O col = col(Matrix,j)
		dirAvg += col
	endfor
	dirAvg /= trials
	
	xpts = dirAvg*cos(dirs*pi/180) //convert each value to x and y  coordinates
	ypts = dirAvg*sin(dirs*pi/180)
	xsum = sum(xpts) //sum across all 8 directions
	ysum = sum(ypts) 
	radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
	theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
	if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
		avgtheta[0] = 360 + theta
	else
		avgtheta[0] = theta
	endif
	avgthetaArrow[0] = theta
	avgthetaArrow[1]= theta
	avgDSi[0] = radius/sum(dirAvg)
	if(displayTrials)
		avgDSiArrow[1] =  avgDSi[0] * 6//wavemax(Matrix)
	else
		avgDSiArrow[1] =  avgDSi[0] * 4
	endif
	
	for(j = 0;j < numpnts(dirs); j+= 1)
		MatrixOp/O row = row(Matrix,j)
		Wavestats/Q row
		dirSDev[j] = V_sdev
		dirSEM[j] = V_sem
	endfor
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~END AVG BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//sort angle and corresponding responses to be clockwise
	Duplicate/O dirs dirSort 
	Duplicate/O dirAvg $(type+"_dirAvgSort_r"+num2str(rho))
	WAVE dirAvgSort = $(type+"_dirAvgSort_r"+num2str(rho))
	Duplicate/O dirSDev $(type+"_dirSDevSort_r"+num2str(rho))
	WAVE dirSDevSort = $(type+"_dirSDevSort_r"+num2str(rho))
	Duplicate/O dirSEM $(type+"_dirSEMSort_r"+num2str(rho))
	WAVE dirSEMSort = $(type+"_dirSEMSort_r"+num2str(rho))
	Duplicate/O Matrix $(type+ "_" + name + "_r" + num2str(rho) + "_srt")
	WAVE sortMatrix = $(type+ "_" + name + "_r" + num2str(rho) + "_srt")
	SortColumns keyWaves=dirSort, sortWaves=sortMatrix//sort angleSort ascending and sort dirDFs to angleSort
	Sort dirSort dirSort, dirAvgSort, dirSDevSort, dirSEMSort
	
	// SEM lines
	Duplicate/O dirAvgSort $(type+"_dirAvgSortMinSEM_r"+num2str(rho)), $(type+"_dirAvgSortPluSEM_r"+num2str(rho))
	WAVE dirAvgSortMinSEM = $(type+"_dirAvgSortMinSEM_r"+num2str(rho))
	WAVE dirAvgSortPluSEM = $(type+"_dirAvgSortPluSEM_r"+num2str(rho))
	dirAvgSortMinSEM -= dirSEMsort
	dirAvgSortPluSEM += dirSEMsort
	Duplicate/O dirAvgSort $(type+"_dirAvgSortMinSD_r"+num2str(rho)), $(type+"_dirAvgSortPluSD_r"+num2str(rho))
	WAVE dirAvgSortMinSD = $(type+"_dirAvgSortMinSD_r"+num2str(rho))
	WAVE dirAvgSortPluSD = $(type+"_dirAvgSortPluSD_r"+num2str(rho))
	dirAvgSortMinSD -= dirSDevSort
	dirAvgSortPluSD += dirSDevSort
	
	// close the loops
	dirSort[numpnts(dirs)] = {dirSort[0]}
	dirAvgSort[numpnts(dirs)] = {dirAvgSort[0]}	
	dirSDevSort[numpnts(dirs)] = {dirSDevSort[0]}
	dirSEMSort[numpnts(dirs)] = {dirSEMSort[0]}	
	dirAvgSortMinSEM[numpnts(dirs)] = {dirAvgSortMinSEM[0]}	
	dirAvgSortPluSEM[numpnts(dirs)] = {dirAvgSortPluSEM[0]}	
	dirAvgSortMinSD[numpnts(dirs)] = {dirAvgSortMinSD[0]}	
	dirAvgSortPluSD[numpnts(dirs)] = {dirAvgSortPluSD[0]}	
	InsertPoints/M=0 numpnts(dirs),1,sortMatrix
	Duplicate/O sortMatrix dirSortMat
	for(j = 0; j < trials;j+=1)
		sortMatrix[numpnts(dirs)][j] = sortMatrix[0][j]
		dirSortMat[][j] = dirSort[p] //also make matrix of sorted dirs
	endfor
	
	if(displayTrials)
		//Build polar plots
		graphname = type + "_" + name + "_r" + num2str(rho) + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, sortMatrix, dirSortMat, 360)	
		WMPolarAppendTrace(graphname, DSiArrows, thetaArrows, 360)
		ModifyGraph rgb = (52428,52428,52428)
		WMPolarAppendTrace(graphname, dirAvgSort, dirSort, 360)	
		WMPolarAppendTrace(graphname, avgDSiArrow, avgthetaArrow, 360)
		//WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		WMPolarSetManualRadiusRange(graphname, 0, 6)
	else
		graphname = type + "_" + name + "_r" + num2str(rho) + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, DSiArrows, thetaArrows, 360)
		WMPolarAppendTrace(graphname, dirAvgSortMinSEM, dirSort, 360)
		WMPolarAppendTrace(graphname, dirAvgSortPluSEM, dirSort, 360)
		//WMPolarAppendTrace(graphname, dirAvgSortMinSD, dirSort, 360)
		//WMPolarAppendTrace(graphname, dirAvgSortPluSD, dirSort, 360)
		ModifyGraph rgb = (52428,52428,52428)
		WMPolarAppendTrace(graphname, dirAvgSort, dirSort, 360)	
		WMPolarAppendTrace(graphname, avgDSiArrow, avgthetaArrow, 360)
		//WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		WMPolarSetManualRadiusRange(graphname, 0, 4)
	endif
	Killwaves/Z col,row, xpts, ypts
End

Function buildGrid(prefix, firstSyn,lastSyn,trial)
	String prefix
	Variable firstSyn, lastSyn, trial
	
	Variable integralMode = 0
	Variable areaMode = 1
	Variable thresholdMode = 1
	Variable threshold = -57
	Variable filterMode = 0
	
	Variable i, j, k, perSpX, perX, perY, perSpY
	String botAxName, leftAxName, rightAxName
	make/O/N=401 timeAx = x*.5
	//directions are [225 270 315 0 45 90 135 180]
	make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	if(!areaMode)
		perSpX = .02
		perX = (1-perSpx*7)/8
	else
		perSpX = .02
		perX = (1-perSpx*8)/9
	endif
	perSpY = .03
	perY = (1-perSpY*(lastSyn-firstSyn))/(lastSyn-firstSyn+1)
	
	for(i = 0; i <= lastSyn-firstSyn; i+=1)
		leftAxName = "s"+num2str(i)	
		rightAxName = "r"+num2str(i)
		if(!filterMode)
			Make/O/N=8 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areas")
			WAVE areas = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areas")
			if(thresholdMode)
				Make/O/N=8 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreas")
				WAVE thrAreas = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreas")
			endif
		else
			Make/O/N=8 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areasFltr")
			WAVE areas = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areasFltr")
		endif
		for(j = 0; j < 8; j+=1)
			botAxName = "d"+num2str(j)
			WAVE vm = $(prefix+"_"+num2str(trial)+"_"+num2str(firstSyn+i)+"_"+num2str(j)+"_1")
			if(filterMode)
				Duplicate/O vm $("fltr_"+prefix+"_"+num2str(trial)+"_"+num2str(firstSyn+i)+"_"+num2str(j)+"_1")
				WAVE vm = $("fltr_"+prefix+"_"+num2str(trial)+"_"+num2str(firstSyn+i)+"_"+num2str(j)+"_1")
				Make/O/D/N=0 coefs
				//FilterFIR/DIM=0/HI={0.0005,0.01,101}/WINF=KaiserBessel20/COEF coefs, vm
				FilterFIR/DIM=0/HI={0.001,0.02,101}/WINF=KaiserBessel20/COEF coefs, vm
			endif
			
			if(thresholdMode)
				Duplicate/O vm temp
				for(k = 0; k < numpnts(temp); k+=1)
					if(temp[k] < threshold)
						temp[k] = threshold
					endif
				endfor
				temp -= wavemin(temp)
				thrAreas[j] = sum(temp)
			endif
			Duplicate/O vm temp
			temp -= wavemin(temp)
			areas[j] = sum(temp)
			if(!i && !j)
				Display/L=$leftAxName/B=$botAxName vm vs timeAx
			else
				AppendToGraph/L=$leftAxName/B=$botAxName vm vs timeAx
			endif
			if(integralMode)
				Duplicate/O vm $("integ_"+num2str(trial)+"_"+num2str(firstSyn+i)+"_"+num2str(j)+"_1")
				WAVE integ = $("integ_"+num2str(trial)+"_"+num2str(firstSyn+i)+"_"+num2str(j)+"_1")
				//Smooth/B 55, integ
				integ -= wavemin(integ)
				Integrate/METH=0 integ
				AppendToGraph/R=$rightAxName/B=$botAxName integ vs timeAx
			endif
			ModifyGraph axisEnab($botAxName)={(j*(perX+perSpX)),(j*(perX+perSpX)+perX)}
			ModifyGraph freePos($botAxName)=0
			if(j == 7 && areaMode)
				Make/O W_coef = {5,5,5}
				FuncFit/Q vonMises W_coef areas /X=dirs /D
				if(thresholdMode)
					Make/O W_coef = {5,5,5}
					FuncFit/Q vonMises W_coef thrAreas /X=dirs /D
				endif
				if(!filterMode)
					if(thresholdMode)
						WAVE thrAreas_fit = $("fit_syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreas")
					endif
						WAVE areas_fit = $("fit_syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areas")
				else
					WAVE areas_fit = $("fit_syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areasFltr")
				endif
				botAxName = "d"+num2str(8)
				if(!thresholdMode)	
					AppendToGraph/R=$rightAxName/B=$botAxName areas vs dirs
					AppendToGraph/R=$rightAxName/B=$botAxName areas_fit
				else
					AppendToGraph/R=$rightAxName/B=$botAxName thrAreas vs dirs
					AppendToGraph/R=$rightAxName/B=$botAxName thrAreas_fit
				endif
				ModifyGraph axisEnab($botAxName)={(8*(perX+perSpX)),(8*(perX+perSpX)+perX)}
				ModifyGraph freePos($botAxName)=0
			endif
		endfor
		ModifyGraph axisEnab($leftAxName)={(i*(perY+perSpY)),(i*(perY+perSpY)+perY)}
		ModifyGraph freePos($leftAxName)=0
		if(integralMode || areaMode)
			ModifyGraph axisEnab($rightAxName)={(i*(perY+perSpY)),(i*(perY+perSpY)+perY)}
			ModifyGraph freePos($rightAxName)=0
		endif
	endfor
	Killwaves/Z temp
End

Function vmGridPolars(firstSyn,lastSyn,trial,rho)
	Variable firstSyn, lastSyn, trial, rho
	
	Variable i, xsum, ysum, radius, theta
	String graphname, polarPath
	WAVE dirs = $("dirs")
	
	Variable areaRange = 1000
	Variable thrAreaRange = 300
	
	for(i = 0; i <= lastSyn-firstSyn; i+=1)
		Duplicate/O dirs dirSort
		
		// grab areas generated by buildGrid()
		WAVE areas = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areas")
		Duplicate/O areas $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areasC")
		WAVE areasCirc = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areasC")
		// also threshold areas
		WAVE thrAreas = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreas")
		Duplicate/O thrAreas $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreasC")
		WAVE thrAreasCirc = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreasC")
		
		// vector calculations
		Make/O/N=1 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arTht")
		Make/O/N=1 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arDSi")
		Make/O/N=1 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArTht")
		Make/O/N=1 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArDSi")
		WAVE areaTheta = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arTht")
		WAVE areaDSi = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arDSi")
		WAVE thrAreaTheta = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArTht")
		WAVE thrAreaDSi = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArDSi")
		Make/O/N=2 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arThtAro")
		Make/O/N=2 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arDSiAro")
		Make/O/N=2 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArThtAro")
		Make/O/N=2 $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArDSiAro")
		WAVE areaThetaAro = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arThtAro")
		WAVE areaDSiAro = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_arDSiAro")
		WAVE thrAreaThetaAro = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArThtAro")
		WAVE thrAreaDSiAro = $("syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"_thrArDSiAro")
		
		Make/O/N=(numpnts(dirs)) xpts, ypts
		//first for plain area
		xpts = areas*cos(dirs*pi/180) //convert each value to x and y  coordinates
		ypts = areas*sin(dirs*pi/180)
		xsum = sum(xpts) //sum across all 8 directions
		ysum = sum(ypts) 
		radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
		theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
		if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
			areaTheta[0] = 360 + theta
		else
			areaTheta[0] = theta
		endif
		areaThetaAro[0] = theta
		areaThetaAro[1]= theta
		areaDSi[0] = radius/sum(areas)
		areaDSiAro[1] =  areaDSi[0] * areaRange//wavemax(areas)
		
		// then for threshold area
		xpts = thrAreas*cos(dirs*pi/180) //convert each value to x and y  coordinates
		ypts = thrAreas*sin(dirs*pi/180)
		xsum = sum(xpts) //sum across all 8 directions
		ysum = sum(ypts) 
		radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
		theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
		if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
			thrAreaTheta[0] = 360 + theta
		else
			thrAreaTheta[0] = theta
		endif
		thrAreaThetaAro[0] = theta
		thrAreaThetaAro[1]= theta
		thrAreaDSi[0] = radius/sum(thrAreas)
		thrAreaDSiAro[1] =  thrAreaDSi[0] * thrAreaRange//wavemax(thrAreas)
		
		//sort and close circles
		Sort dirSort dirSort, areasCirc, thrAreasCirc
		dirSort[numpnts(dirs)] = {dirSort[0]}
		areasCirc[numpnts(dirs)] = {areasCirc[0]}
		thrAreasCirc[numpnts(dirs)] = {thrAreasCirc[0]}
		
		// make polarplots
		// regular area
		graphname =  "rho"+num2str(rho)+"_"+"syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"areas"+"_Polar"
		WMNewPolarGraph("_default_", graphname)
		//polarPath = "root:Packages:WMPolarGraphs:"+graphname+":"
		//NVAR temp = $(polarPath+"doRadiusTickLabels")
		//temp = 0
		WMPolarGraphSetVar(graphname,"doAngleTickLabels",0)
		WMPolarGraphSetVar(graphname,"doMinorRadiusTicks",0)
		WMPolarGraphSetVar(graphname,"doMinorAngleTicks",0)
		WMPolarGraphSetVar(graphname,"majorAngleInc",90)
		WMPolarGraphSetVar(graphname,"innerRadius",0)
		WMPolarGraphSetVar(graphname,"outerRadius",thrAreaRange)
		WMPolarGraphSetVar(graphname,"majorRadiusInc",(areaRange/2))
		//put the vectors and areas on
		WMPolarAppendTrace(graphname, areasCirc, dirSort, 360)
		WMPolarAppendTrace(graphname, areaDSiAro, areaThetaAro, 360)
		WMPolarSetManualRadiusRange(graphname, 0, areaRange)
		// threshold area
		graphname = "rho"+num2str(rho)+"_"+"syn"+num2str(i+firstSyn)+"tr"+num2str(trial)+"thrAreas"+"_Polar"
		polarPath = "root:Packages:WMPolarGraphs:'_default_':"
		WMNewPolarGraph("_default_", graphname)
		polarPath = "root:Packages:WMPolarGraphs:"+graphname+":"
		//NVAR temp = $(polarPath+"doRadiusTickLabels")
		//temp = 0
		WMPolarGraphSetVar(graphname,"doAngleTickLabels",0)
		WMPolarGraphSetVar(graphname,"doMinorRadiusTicks",0)
		WMPolarGraphSetVar(graphname,"doMinorAngleTicks",0)
		WMPolarGraphSetVar(graphname,"majorAngleInc",90)
		WMPolarGraphSetVar(graphname,"innerRadius",0)
		WMPolarGraphSetVar(graphname,"outerRadius",thrAreaRange)
		WMPolarGraphSetVar(graphname,"majorRadiusInc",(thrAreaRange/2))
		//put the vectors and areas on
		WMPolarAppendTrace(graphname, thrAreasCirc, dirSort, 360)
		WMPolarAppendTrace(graphname, thrAreaDSiAro, thrAreaThetaAro, 360)
		WMPolarSetManualRadiusRange(graphname, 0, thrAreaRange)
		//WAVE/T tw=$WMPolarGraphTracesTW(graphname)
		//WMPolarUpdateAxes(tw,1)	
	endfor

End
Function currVmGrid(balTrial, unbalTrial)
	Variable balTrial, unbalTrial

	Variable i, j, k, perSpX, perX, perY, perSpY
	String botAxName, leftAxName, rightAxName
	Make/O/N=401 timeAx_Vm = x*.5
	Make/O/N=371 timeAx_I = x*.5 + 30*.5 //30 deleted points
	Variable numRows = 3
	//directions are [225 270 315 0 45 90 135 180]
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	Make/O balClr = {0,0,0}
	Make/O unbalClr = {65535,16385,16385}
	perSpX = .02
	perX = (1-perSpx*7)/8
	perSpY = .03
	perY = (1-perSpY*numRows)/(numRows)
	
	Make/O/N=(371,8) balExcDirs, balInhDirs, unbalExcDirs, unbalInhDirs
	Make/O/N=(401,8) balVmDirs,unbalVmDirs
	//currents
	leftAxName = "current"
	WAVE balExcMat = $("rho9_trialExc")
	WAVE balInhMat = $("rho9_trialInh")
	WAVE unbalExcMat = $("rho0_trialExc")
	WAVE unbalInhMat = $("rho0_trialInh")
	for(i = 0; i < 8; i+=1)
		botAxName = "d" + num2str(i)
		MatrixOp/O temp = col(unbalInhMat,(unbalTrial*8+i))
		unbalInhDirs[][i] = temp[p]
		MatrixOp/O temp = col(unbalExcMat,(unbalTrial*8+i))
		unbalExcDirs[][i] = temp[p]
		MatrixOp/O temp = col(balInhMat,(balTrial*8+i))
		balInhDirs[][i] = temp[p]
		MatrixOp/O temp = col(balExcMat,(balTrial*8+i))
		balExcDirs[][i] = temp[p]
		if(!i)
			Display
		endif
		AppendToGraph/L=$leftAxName/B=$botAxName/C=(unbalClr[0],unbalClr[1],unbalClr[2]) unbalInhDirs[][i] vs timeAx_I
		AppendToGraph/L=$leftAxName/B=$botAxName/C=(unbalClr[0],unbalClr[1],unbalClr[2]) unbalExcDirs[][i] vs timeAx_I
		AppendToGraph/L=$leftAxName/B=$botAxName/C=(balClr[0],balClr[1],balClr[2]) balInhDirs[][i] vs timeAx_I
		AppendToGraph/L=$leftAxName/B=$botAxName/C=(balClr[0],balClr[1],balClr[2]) balExcDirs[][i] vs timeAx_I
		//adjust bottom axis position
		ModifyGraph axisEnab($botAxName)={(i*(perX+perSpX)),(i*(perX+perSpX)+perX)}
		ModifyGraph freePos($botAxName)=0
	endfor
	//adjust left axis position
	ModifyGraph axisEnab($leftAxName)={(2*(perY+perSpY)),(2*(perY+perSpY)+perY)}
	ModifyGraph freePos($leftAxName)=0
		
	//balanced spiking Vm
	leftAxName = "balVm"
	WAVE spkMat = $("rho9_Vm")
	for(i = 0; i < 8; i+=1)
		botAxName = "d" + num2str(i)
		MatrixOp/O temp = col(spkMat,(balTrial*8+i))
		balVmDirs[][i] = temp[p]
		AppendToGraph/L=$leftAxName/B=$botAxName/C=(balClr[0],balClr[1],balClr[2]) balVmDirs[][i] vs timeAx_Vm
		//adjust bottom axis position
		ModifyGraph axisEnab($botAxName)={(i*(perX+perSpX)),(i*(perX+perSpX)+perX)}
		ModifyGraph freePos($botAxName)=0
	endfor
	//adjust left axis position
	ModifyGraph axisEnab($leftAxName)={(1*(perY+perSpY)),(1*(perY+perSpY)+perY)}
	ModifyGraph freePos($leftAxName)=0
	
	//unbalanced spiking Vm
	leftAxName = "unbalVm"
	WAVE spkMat = $("rho0_Vm")
	for(i = 0; i < 8; i+=1)
		botAxName = "d" + num2str(i)
		MatrixOp/O temp = col(spkMat,(unbalTrial*8+i))
		unbalVmDirs[][i] = temp[p]
		AppendToGraph/L=$leftAxName/B=$botAxName/C=(unbalClr[0],unbalClr[1],unbalClr[2]) unbalVmDirs[][i] vs timeAx_Vm
		//adjust bottom axis position
		ModifyGraph axisEnab($botAxName)={(i*(perX+perSpX)),(i*(perX+perSpX)+perX)}
		ModifyGraph freePos($botAxName)=0
	endfor
	//adjust left axis position
	ModifyGraph axisEnab($leftAxName)={(0*(perY+perSpY)),(0*(perY+perSpY)+perY)}
	ModifyGraph freePos($leftAxName)=0
End

Function tripletPuller(Matrix)
	WAVE Matrix
	//pull out waves and name in this format:
	//(prefix+"_"+num2str(trial)+"_"+num2str(firstSyn+i)+"_"+num2str(j)+"_1")
	//Vm_trialNum_synNum_dirNum_1 (because other func already written)
	//synNum starts from 0 (not same as true dendNumber or termNumber)
	//dirNum is 0 to 7
	Variable numSyns = 3
	Variable numDirs = 8
	Variable numTrials = DimSize(Matrix,1)/(numSyns*numDirs)
	Variable i, j, k
	
	for(i = 0; i < numTrials; i += 1)
		for(j = 0; j < numDirs; j += 1)
			for(k = 0; k < numSyns; k += 1)
				MatrixOp/O col = col(Matrix,(i*numDirs*numSyns+j*numSyns+k))
				Duplicate/O col $("Vm_"+num2str(i)+"_"+num2str(k)+"_"+num2str(j)+"_1")
			endfor
		endfor
	endfor
	Killwaves/Z col
End

Function treeRecLoader(path, numTrials)
	//e.g. "C:\Users\Geoff\Desktop\\treeRecs\\rho9\dirTreeRecs"
	//the format of these matrices will be like this:
	//|------dir0------||------dir1------|....
	//|seg0 seg1 seg2..||seg0 seg1 seg2..|
	//so, 8 blocks of 700 (dends*2) per trial matrix
	String path
	Variable numTrials
	
	Variable i
	String fullpath
	
	for(i = 0; i < numtrials; i += 1)
		fullpath = path + num2str(i) + ".dat"
		LoadWave/A=trial/O/Q/M/D/H/G/K=1/ENCG={1,4} fullpath
	endfor	
End

Function treeRecDS(mode)
	Variable mode //0:Vm; 1:iCa
	
	Variable i, j, k, m, count, wmin
	Variable xsum, ysum, theta, radius
	String pathToWave
	
	Variable numDirs = 8
	Variable numSegs = 700
	Variable thresholdVm = -50//-55//-50//-55//-52 //-57
	Variable thresholdCa = 0.00015//0.00015//.0003
	Variable CaPeak = 1 // use peak current instead of area
	Variable saveDir = 7
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = "trial" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	Make/O/N=(numDirs) thrAreas, xpts, ypts
	Make/O/N=(numSegs) thetas, absThetas, DSis, dirThrAreas
	for(i = 0; i < numpnts(trials); i += 1)
		WAVE Matrix = $("trial" + num2str(trials[i]))
		for(j = 0; j < numSegs; j+=1)
			for(k = 0; k < numDirs; k += 1)
				//pull a seg recording from current direction
				MatrixOp/O seg = col(Matrix,(k*numSegs+j))
				if(!mode) // Vm
					//get thresholded area of Vm
					for(m = 0; m < numpnts(seg); m+=1)
						if(seg[m] < thresholdVm)
							seg[m] = thresholdVm
						endif
					endfor
					wmin = wavemin(seg)
					seg -= wmin
					thrAreas[k] = sum(seg)
				else // iCa
					//get thresholded area of Ca
					seg = seg*-1 // flip the negative current
					if(!CaPeak)
						for(m = 0; m < numpnts(seg); m+=1)
							if(seg[m] < thresholdCa)
								seg[m] = thresholdCa
							endif
						endfor
						wmin = wavemin(seg)
						seg -= wmin
						thrAreas[k] = sum(seg)
					else
						thrAreas[k] = wavemax(seg)
					endif
				endif
				
				if(k==saveDir)
					dirThrAreas[j] = thrAreas[k]
				endif
			endfor
			xpts = thrAreas*cos(dirs*pi/180) //convert each value to x and y  coordinates
			ypts = thrAreas*sin(dirs*pi/180)
			xsum = sum(xpts) //sum across all 8 directions
			ysum = sum(ypts) 
			radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
			theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
			//if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
			//	thetas[j] = 360 + theta
			//else
			//	thetas[j] = theta
			//endif
			thetas[j] = theta
			absThetas[j] = abs(theta)
			DSis[j] = radius/sum(thrAreas)
			
		endfor
		Duplicate/O thetas $("thetas"+num2str(trials[i]))
		Duplicate/O absThetas $("absThetas"+num2str(trials[i]))
		Duplicate/O DSis $("DSis"+num2str(trials[i]))
		Duplicate/O dirThrAreas $("dirThrAreas"+num2str(trials[i]))
	endfor
	KillWaves/Z thetas, absThetas, DSis, seg, xpts, ypts
End

Function treeRecSingle(mode)
	Variable mode //0:Vm; 1:iCa
	
	Variable i, j, k, m, count, wmin
	String pathToWave
	
	Variable numSegs = 700
	Variable thresholdVm = -50//-55//-50//-55//-52 //-57
	Variable thresholdCa = 0.00015//0.00015//.0003
	Variable CaPeak = 1 // use peak current instead of area

	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 100; i += 1)
		pathToWave = "trial" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Make/O/N=(numSegs) dirThrAreas
	for(i = 0; i < numpnts(trials); i += 1)
		WAVE Matrix = $("trial" + num2str(trials[i]))
		for(j = 0; j < numSegs; j+=1)
			//pull a seg recording from current direction
			MatrixOp/O seg = col(Matrix,j)
			if(!mode) // Vm
				//get thresholded area of Vm
				for(m = 0; m < numpnts(seg); m+=1)
					if(seg[m] < thresholdVm)
						seg[m] = thresholdVm
					endif
				endfor
				wmin = wavemin(seg)
				seg -= wmin
				dirThrAreas[j] = sum(seg)
			else // iCa
				//get thresholded area of Ca
				seg = seg*-1 // flip the negative current
				if(!CaPeak)
					for(m = 0; m < numpnts(seg); m+=1)
						if(seg[m] < thresholdCa)
							seg[m] = thresholdCa
						endif
					endfor
					wmin = wavemin(seg)
					seg -= wmin
					dirThrAreas[j] = sum(seg)
				else
					dirThrAreas[j] = wavemax(seg)
				endif
			endif
			
		endfor
		Duplicate/O dirThrAreas $("dirThrAreas"+num2str(trials[i]))
	endfor
	KillWaves/Z seg, xpts, ypts, dirThrAreas

End

Function thetaSegPainter(trial,rangeMax)
	Variable trial, rangeMax
	
	Variable i, clr
	WAVE xLoc = root:xLoc
	WAVE yLoc = root:yLoc
	WAVE colorWave = root:MyColorWave
	
	Variable rangeMin = 0
	
	WAVE thetas = $("absThetas"+num2str(trial))
	for(i = 0; i < numpnts(yLoc); i += 1)
		clr = thetas[i]/rangeMax * DimSize(colorWave,0)
		if(clr > DimSize(colorWave,0)-1)
			clr = DimSize(colorWave,0)-1
		endif
		ModifyGraph rgb(yLoc[i])=(colorWave[clr][0],colorWave[clr][1],colorWave[clr][2])
	endfor
End

Function DSiSegPainter()
	
	Variable i, clr
	WAVE xLoc = root:xLoc
	WAVE yLoc = root:yLoc
	WAVE colorWave = root:MyColorWave
	
	Variable rangeMax = .8
	
	WAVE DSis = $("DSis0")
	for(i = 0; i < numpnts(yLoc); i += 1)
		clr = DSis[i]/rangeMax * DimSize(colorWave,0)
		if(clr > DimSize(colorWave,0)-1)
			clr = DimSize(colorWave,0)-1
		endif
		ModifyGraph rgb(yLoc[i])=(colorWave[clr][0],colorWave[clr][1],colorWave[clr][2])
	endfor
End

Function cableCalc(prefix,mode)
	String prefix
	Variable mode //0:regular; 1:only for sites above a threshold; 2:residual mode
	
	Variable i, k, j, count, numLocs, iSgn, jSgn, numCombos
	Variable lastBin, thisBin
	String pathToWave
	WAVE cableDists = root:cableDistances
	
	Variable threshold, thrAreaMode
	
	if(!CmpStr(prefix,"dirThrAreas"))
		thrAreaMode = 1
	else
		thrAreaMode = 0
	endif
	
	threshold = 30//.00002//10//0.0001 // area based (look at dirThrAreas waves to get idea)
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 15; i += 1)
		pathToWave = prefix + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Variable binSize = 3 //number of microns cable dist
	Variable XYbins = 60
	
	numLocs = DimSize(cableDists,0)
	KillWaves/Z diffs, dists
	Make/O/N=1 diffs, dists
	
	for(k = 0; k < numpnts(trials); k+=1)
		WAVE thetas = $(prefix + num2str(trials[k]))
		count = 0
		for(i = 0; i < numLocs; i += 1)
			//ADD IN: (if mode=1) check if synapse is above certain threshold
			//if not above, skip to next one. Only want comparisons from sites
			//where something happened
			if((thrAreaMode && mode==1 && thetas[i] > threshold) || mode != 1)
				iSgn = sign(thetas[i])
				for(j = i+1; j < numLocs; j+= 1)
					jSgn = sign(thetas[j])
					if((iSgn + jSgn)==0)
						diffs[count] = {abs(thetas[i])+abs(thetas[j])}
					else
						diffs[count] = {abs(thetas[i]-thetas[j])}
					endif
					dists[count] = {cableDists[i][j]}
					count += 1
				endfor
			endif
		endfor
		
		//now binning operations are needed
		Killwaves/Z diffSort,distSort
		Duplicate/O diffs diffSort
		Duplicate/O dists distSort
		Sort distSort, distSort, diffSort
		
		KillWaves/Z sumBins, avgBins, medBins
		Make/O/N=(floor(wavemax(cableDists)/binSize)+1) sumBins, avgBins, medBins
		SetScale/P x binSize,binSize,"", sumBins, avgBins, medBins
		lastBin = 0
		//find pts for beginning and end of bin and pool data
		for(i = 0; i < numpnts(sumBins); i += 1)
			
			FindLevel/Q/P distSort, ((i+1)*binSize)
			if(!V_flag)
				thisBin = V_LevelX
			else
				thisBin = numpnts(distSort)-1
			endif
			sumBins[i] = sum(diffSort,lastBin,thisBin)
			//avgBins[i] = sumBins[i]/(thisBin-lastBin+1)
			avgBins[i] = sum(diffSort,lastBin,thisBin)/(thisBin-lastBin+1)
			medBins[i] = median(diffSort,lastBin,thisBin)
			
			lastBin = thisBin
		endfor
		
		Duplicate/O diffSort $(prefix+num2str(trials[k])+"_diffs")
		Duplicate/O distSort $(prefix+num2str(trials[k])+"_dists")
		Duplicate/O sumBins $(prefix+num2str(trials[k])+"_sumBins")
		Duplicate/O avgBins $(prefix+num2str(trials[k])+"_avgBins")
		Duplicate/O medBins $(prefix+num2str(trials[k])+"_medBins")
		
		if(!k)
			Killwaves/Z $(prefix+"AVG_sumBins"), $(prefix+"AVG_avgBins")
			Duplicate/O sumBins $(prefix+"AVG_sumBins")
			Duplicate/O avgBins $(prefix+"AVG_avgBins")
			WAVE sumBinsAVG = $(prefix+"AVG_sumBins")
			WAVE avgBinsAVG = $(prefix+"AVG_avgBins")
		else
			sumBinsAVG += sumBins
			avgBinsAVG += avgBins
		endif
	endfor
	
	//calculate trial avgs/med
	sumBinsAVG /= numpnts(trials)
	avgBinsAVG /= numpnts(trials)
	Killwaves/Z $(prefix+"MED_medBins")
	Make/O/N=(numpnts(cableDists)) $(prefix+"MED_medBins")
	WAVE medBinsMED = $(prefix+"MED_medBins")
	SetScale/P x binSize,binSize,"",medBinsMED
	Duplicate/O $(prefix+num2str(trials[0])+"_medBins") medMatrix
	for(i = 1; i < numpnts(trials); i += 1)
		WAVE next = $(prefix+num2str(trials[i])+"_medBins")
		Concatenate {medMatrix,next}, temp
		Duplicate/O temp medMatrix
	endfor
	for(i = 0; i < numpnts(medBins); i += 1)
		MatrixOp/O bin = row(medMatrix,i)
		medBinsMED[i] = median(bin)
	endfor
	
	//create X (dist bin) and Y (diff) waves of all trials for scatters (fitting)
	Make/O/N=1 $(prefix+"X_medBins"), $(prefix+"Y_medBins")
	Make/O/N=1 $(prefix+"X_avgBins"), $(prefix+"Y_avgBins")
	WAVE avgBinsX = $(prefix+"X_avgBins")
	WAVE avgBinsY = $(prefix+"Y_avgBins")
	WAVE medBinsX = $(prefix+"X_medBins")
	WAVE medBinsY = $(prefix+"Y_medBins")
	count = 0
	for(i = 0; i < numpnts(trials); i += 1)
		WAVE nextAVG = $(prefix+num2str(trials[i])+"_avgBins")
		WAVE nextMED = $(prefix+num2str(trials[i])+"_medBins")
		for(k = 0; k < XYbins; k += 1)
			avgBinsX[count] = {binSize+binSize*k}
			medBinsX[count] = {binSize+binSize*k}
			avgBinsY[count] = {nextAVG[k]}
			medBinsY[count] = {nextMED[k]}
			count += 1
		endfor
	endfor
	
	//trial variance calculations
	KillWaves/Z corrs, dists
	Make/O/N=1 corrs, dists
	Make/O/N=(numpnts(trials)) iTemp,jTemp
	count = 0
	for(i = 0; i < numLocs-1; i += 1)
			
		for(j = i+1; j < numLocs; j+= 1)
			for(k = 0; k < numpnts(trials); k+=1)
				WAVE thetas = $(prefix + num2str(trials[k]))
				iTemp[k] = thetas[i]
				jTemp[k] = thetas[j]
			endfor
			corrs[count] = {StatsCorrelation(iTemp,jTemp)}
			dists[count] = {cableDists[i][j]}
			if(numtype(corrs[count]) == 2)
				corrs[count] = 0
			endif
			count += 1
		endfor
		
	endfor
	
	//now binning operations are needed, again.
	KillWaves/Z corrSort, distSort
	Duplicate/O corrs corrSort
	Duplicate/O dists distSort
	Sort distSort, distSort, corrSort
	
	Make/O/N=(floor(wavemax(dists)/binSize)+1) sumCorrBins, avgCorrBins, medCorrBins
	SetScale/P x binSize,binSize,"", avgCorrBins, medCorrBins
	lastBin = 0
	//find pts for beginning and end of bin and pool data
	for(i = 0; i < numpnts(sumCorrBins); i += 1)
		
		FindLevel/Q/P distSort, ((i+1)*binSize)
		if(!V_flag)
			thisBin = V_LevelX
		else
			thisBin = numpnts(distSort)-1
		endif
		sumCorrBins[i] = sum(corrSort,lastBin,thisBin)
		avgCorrBins[i] = sumCorrBins[i]/(thisBin-lastBin+1)
		medCorrBins[i] = median(corrSort,lastBin,thisBin)
		
		lastBin = thisBin
	endfor
	
	Duplicate/O avgCorrBins $(prefix+"_avgCorrBins")
	Duplicate/O medCorrBins $(prefix+"_medCorrBins")
	
	KillWaves/Z dists,diffs,distSort,diffSort,sumBins,avgBins,medBins
	KillWaves/Z sumCorrBins,avgCorrBins,medCorrBins,medMatrix,bin,temp
	KillWaves/Z iTemp,jTemp
End

Function linePuller(trial)
	Variable trial
	
	Variable i, j, k, count
	Variable numDirs = 8
	Variable numSegs = 700
	Variable segsPerSec = 2 // default recording pos .5 and 1 of each section
	Variable threshold = -55 //-57
	Variable saveDir = 0
	String pathToWave
	
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = "trial" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Make/O sections = {268, 158, 159, 160, 161, 162, 163, 164, 275}
	
	WAVE Matrix = $("trial" + num2str(trial))
	
	// for chosen trial, create a trimmed down matrix that just has the sections
	// in order in blocks of direction for simplified line plotting
	for(k = 0; k < numDirs; k += 1)
		
		for(i = 0; i < numpnts(sections); i += 1)
		
			for(j = 0; j < segsPerSec; j+=1)
				MatrixOp/O col = col(Matrix,(k*numSegs+sections[i]+j))
				
			endfor
		endfor
	endfor
End

Function cableResid(prefix)
	String prefix

	Variable i, k, j, count, numLocs, iSgn, jSgn, numCombos
	Variable lastBin, thisBin, h
	String pathToWave
	WAVE cableDists = root:cableDistances
	
	// {-135, -90 , -45 , 0, 45, 90, 135, 180}
	Variable numDirs = 8
	Make/O targetDirs = {7}
	
	//find how many trials there are and create an array to loop over for this function
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = prefix + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	Variable numTrials = numpnts(trials)
	
	Variable binSize = 3 //number of microns cable dist
	Variable XYbins = 40
	
	numLocs = DimSize(cableDists,0)
	KillWaves/Z corrs, dists
	Make/O/N=1 corrs, dists
	
	WAVE Matrix = $(prefix+num2str(trials[0]))
	Variable numRows = DimSize(Matrix,0)
	Variable numCols = DimSize(Matrix,1)
	
	KillWaves/Z avgMat
	Make/O/N=(numRows,numCols) avgMat
	for(k = 0; k < numpnts(trials); k += 1)
		WAVE Matrix = $(prefix+num2str(trials[k]))
		MatrixOp/O avgMat = avgMat + Matrix
	endfor
	MatrixOp/O avgMat = avgMat / numTrials
		
	count = 0
	for(k = 0; k < numpnts(trials); k+=1)
		WAVE Matrix = $(prefix + num2str(trials[k]))
		MatrixOp/O residuals = Matrix - avgMat
		for(j = 0; j < numpnts(targetDirs); j += 1)
			for(i = 0; i < numLocs; i += 1)
				MatrixOp/O iTemp = col(residuals,(targetDirs[j]*numLocs+i))
				//MatrixOp/O iTemp = col(Matrix,(targetDirs[j]*numLocs+i)) //raw (non-residual)
				for(h = i+1; h < numLocs; h += 1)
					MatrixOp/O hTemp = col(residuals,(targetDirs[j]*numLocs+h))
					//MatrixOp/O hTemp = col(Matrix,(targetDirs[j]*numLocs+h))//raw (non-residual)
					corrs[count] = {StatsCorrelation(iTemp,hTemp)} 
					dists[count] = {cableDists[i][h]}
					count += 1
				endfor
			endfor
		endfor
		print trials[k] //progress
	endfor	
	
	//now binning operations are needed
	Killwaves/Z corrSort,distSort
	Duplicate/O corrs corrSort
	Duplicate/O dists distSort
	Sort distSort, distSort, corrSort
	
	KillWaves/Z avgBins, medBins
	Make/O/N=(floor(wavemax(dists)/binSize)+1) avgBins, medBins
	SetScale/P x binSize,binSize,"", avgBins, medBins
	lastBin = 0
	//find pts for beginning and end of bin and pool data
	for(i = 0; i < numpnts(avgBins); i += 1)
		
		FindLevel/Q/P distSort, ((i+1)*binSize)
		if(!V_flag)
			thisBin = V_LevelX
		else
			thisBin = numpnts(distSort)-1
		endif
		
		avgBins[i] = sum(corrSort,lastBin,thisBin)/(thisBin-lastBin+1)
		medBins[i] = median(corrSort,lastBin,thisBin)
		
		lastBin = thisBin
	endfor
	
	Duplicate/O corrSort $(prefix+"_resCorrs")
	Duplicate/O distSort $(prefix+"_resDists")
	Duplicate/O avgBins $(prefix+"_resAvgBins")
	Duplicate/O medBins $(prefix+"_resMedBins")
	
	// cleanup
	KillWaves/Z dists,corrs,distSort,corrSort,avgBins,medBins
	KillWaves/Z bin,temp,iTemp,hTemp,residuals
End

// grab and use Vm from indicated segments (see locations mapped on to DSGC)
// middle segment (.5) of terminal sections will have synapses, since .5 and 1
// positions are recorded, dendNum will be /2 of xyLoc index on graph
// recordings used for polars etc are indicated in args {#,#,#,.. }
Function synPuller(prefix, segNums, mode)
	String prefix
	Wave segNums
	Variable mode
	
	Variable i,j,k,m,count
	Variable numSegs, numDirs, thresholdVm, thrAreaRange
	Variable thresholdCa, CaPeak, wmin
	String pathToWave, graphname, polarPath
	
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	numSegs = 700
	numDirs = numpnts(dirs)
	thresholdVm = -55
	CaPeak = 0 // use peak current instead of area
	thresholdCa = .00015//0
	thrAreaRange = 250 // set radius of polarplot (~max thrArea value)
	
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = "trial" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Variable xsum, ysum, radius,theta, DSi
	Make/O/N=(numDirs) thrAreas, xpts, ypts
	Make/O/N=2 thrAreaDSiAro, thrAreaThetaAro
	for(i = 0; i < numpnts(trials); i += 1)
		WAVE Matrix = $("trial" + num2str(trials[i]))
		for(j = 0; j < numpnts(segNums); j+=1)
			for(k = 0; k < numDirs; k += 1)
				//pull a seg recording from current direction
				MatrixOp/O seg = col(Matrix,(k*numSegs+segNums[j]))
				if(!mode) // Vm
					//get thresholded area of Vm
					for(m = 0; m < numpnts(seg); m+=1)
						if(seg[m] < thresholdVm)
							seg[m] = thresholdVm
						endif
					endfor
					wmin = wavemin(seg)
					seg -= wmin
					thrAreas[k] = sum(seg)
				else // iCa
					//get thresholded area of Ca
					seg = seg*-1 // flip the negative current
					if(!CaPeak)
						for(m = 0; m < numpnts(seg); m+=1)
							if(seg[m] < thresholdCa)
								seg[m] = thresholdCa
							endif
						endfor
						wmin = wavemin(seg)
						seg -= wmin
						thrAreas[k] = sum(seg)
					else
						thrAreas[k] = wavemax(seg)
					endif
				endif
//				//pull a seg recording from current direction
//				MatrixOp/O col = col(Matrix,(k*numSegs+segNums[j]))
//				//get thresholded area
//				for(m = 0; m < numpnts(col); m+=1)
//					if(col[m] < thresholdVm)
//						col[m] = thresholdVm
//					endif
//				endfor
//				col -= wavemin(col)
//				thrAreas[k] = sum(col)
			endfor
			xpts = thrAreas*cos(dirs*pi/180) //convert each value to x and y  coordinates
			ypts = thrAreas*sin(dirs*pi/180)
			xsum = sum(xpts) //sum across all 8 directions
			ysum = sum(ypts) 
			radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
			// metrics
			theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
			DSi= radius/sum(thrAreas)
			if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
				theta = 360 + theta
			endif
			// polar arrows
			thrAreaThetaAro[0] = theta
			thrAreaThetaAro[1]= theta
			thrAreaDSiAro[0] = 0
			thrAreaDSiAro[1] =  DSi * wavemax(thrAreas)//thrAreaRange
			
			//duplicate waves for sorting
			Duplicate/O dirs dirSort
			Duplicate/O thrAreas thrAreasCirc
			//sort and close circles
			Sort dirSort dirSort, thrAreasCirc
			dirSort[numpnts(dirs)] = {dirSort[0]}
			thrAreasCirc[numpnts(dirs)] = {thrAreasCirc[0]}
			
			// make polarplots
			graphname = "dend"+num2str(segNums[j])+"_tr"+num2str(trials[i])+"_thrAreas"+"_Polar"
			polarPath = "root:Packages:WMPolarGraphs:'_default_':"
			WMNewPolarGraph("_default_", graphname)
			polarPath = "root:Packages:WMPolarGraphs:"+graphname+":"
		
			// alter polar graph appearance
			WMPolarGraphSetVar(graphname,"doAngleTickLabels",0)
			WMPolarGraphSetVar(graphname,"doRadiusTickLabels",0)
			WMPolarGraphSetVar(graphname,"doMinorRadiusTicks",0)
			WMPolarGraphSetVar(graphname,"doMinorAngleTicks",0)
			WMPolarGraphSetVar(graphname,"majorAngleInc",90)//90
			WMPolarGraphSetVar(graphname,"innerRadius",0)
			WMPolarGraphSetVar(graphname,"outerRadius",wavemax(thrAreasCirc))//thrAreaRange
			WMPolarGraphSetVar(graphname,"majorRadiusInc",wavemax(thrAreasCirc)/2)//(thrAreaRange/2)
			WMPolarGraphSetVar(graphname,"radiusApproxTicks",1)
			WMPolarGraphSetVar(graphname,"angleApproxTicks",4)
			
			//save data to be plotted in a junk folder to prevent overwrite
			Duplicate/O thrAreasCirc $(GetDataFolder(1)+"junk:thrArCirc_dend"+num2str(segNums[j])+"_tr"+num2str(trials[i]))
			WAVE thrArCirc = $(GetDataFolder(1)+"junk:thrArCirc_dend"+num2str(segNums[j])+"_tr"+num2str(trials[i]))
			Duplicate/O thrAreaThetaAro $(GetDataFolder(1)+"junk:thrArThetaAro_dend"+num2str(segNums[j])+"_tr"+num2str(trials[i]))
			WAVE thrArThetaAro = $(GetDataFolder(1)+"junk:thrArThetaAro_dend"+num2str(segNums[j])+"_tr"+num2str(trials[i]))
			Duplicate/O thrAreaDSiAro $(GetDataFolder(1)+"junk:thrArDSiAro_dend"+num2str(segNums[j])+"_tr"+num2str(trials[i]))
			WAVE thrArDSiAro = $(GetDataFolder(1)+"junk:thrArDSiAro_dend"+num2str(segNums[j])+"_tr"+num2str(trials[i]))
			
			//put the vectors and areas on
			WMPolarAppendTrace(graphname, thrArCirc, dirSort, 360)
			WMPolarAppendTrace(graphname, thrArDSiAro, thrArThetaAro, 360)
			WMPolarSetManualRadiusRange(graphname, 0, wavemax(thrArCirc))//thrAreaRange
			//WAVE/T tw=$WMPolarGraphTracesTW(graphname)
			//WMPolarUpdateAxes(tw,1)	
		endfor
	endfor
End

Function buildTheGrid(prefix, trial,segNums)
	String prefix
	Variable trial
	Wave segNums
	
	Variable numSegs = 700
	Variable threshold = -55
	
	Variable i, j, k, perSpX, perX, perY, perSpY, count
	String botAxName, leftAxName, rightAxName, pathToWave
	
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = prefix + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Make/O/N=401 timeAx = x*.5
	//directions are [225 270 315 0 45 90 135 180]
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	perSpX = .02
	perX = (1-perSpx*7)/8
	perSpY = .03
	perY = (1-perSpY*(numpnts(segNums)-1))/numpnts(segNums)
	
	WAVE Matrix = $(prefix + num2str(trial))
	
	for(i = 0; i < numpnts(segNums); i+=1)
		leftAxName = "s"+num2str(i)	
		rightAxName = "r"+num2str(i)
		//WAVE Matrix = $(prefix + num2str(i))
		for(j = 0; j < 8; j+=1)
			botAxName = "d"+num2str(j)
			//pull a seg recording from current direction
			MatrixOp/O vm = col(Matrix,(j*numSegs+segNums[i]))
			
			if(!i && !j)
				//Display/L=$leftAxName/B=$botAxName vm vs timeAx
				Display/L=$leftAxName/B=$botAxName Matrix[][(j*numSegs+segNums[i])] vs timeAx
			else
				//AppendToGraph/L=$leftAxName/B=$botAxName vm vs timeAx
				AppendToGraph/L=$leftAxName/B=$botAxName Matrix[][(j*numSegs+segNums[i])] vs timeAx
			endif
			ModifyGraph axisEnab($botAxName)={(j*(perX+perSpX)),(j*(perX+perSpX)+perX)}
			ModifyGraph freePos($botAxName)=0
		endfor
		ModifyGraph axisEnab($leftAxName)={(i*(perY+perSpY)),(i*(perY+perSpY)+perY)}
		ModifyGraph freePos($leftAxName)=0
	endfor
	Killwaves/Z temp
End

Function powerTest(prefix, direction)
	String prefix
	Variable direction
	
	Variable i,j,pnts,count
	String pathToWave
	Variable numSegs = 700
	Variable numFreqs = 600
	Variable freqRes = .1
	
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = "trial" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	WAVE Matrix = $(prefix + num2str(trials[0]))
	if(mod(DimSize(Matrix,0),2))
		pnts = DimSize(Matrix,0)+1
	else
		pnts = DimSize(Matrix,0)
	endif
	Make/O/N=(pnts) timeWave = x*.0005
	
	count = 0
	for(i = 0; i < numpnts(trials); i += 1)
		WAVE Matrix = $(prefix + num2str(trials[i]))
		for(j = 0; j < numSegs; j+=1)
			MatrixOp/O seg = col(Matrix,(direction*numSegs+j))
			if(mod(numpnts(seg),2))
				seg[numpnts(seg)] = {seg[numpnts(seg)-1]}
			endif
			LombPeriodogram/NF=(numFreqs)/FR=(freqRes) timeWave, seg
			
			WAVE periodogram = W_LombPeriodogram
			if(!j)
				Duplicate/O periodogram periodAVG
			else
				periodAVG += periodogram
			endif
			count += 1
		endfor
	endfor
	periodAVG /= count
End

//matrix of 8 direction spike number data. no funny business
Function simplePolar(prefix, displayTrials)
	String prefix
	Variable displayTrials
	
	Variable i, j, count, trials
	Variable xsum, ysum, radius, theta
	String pathToMatrix, graphname
	
	Make/O dirs = {0, 45, 90, 135, 180, 225, 270, 315}
	
	//To get appropriate matrix dimensions
	WAVE Matrix = $(prefix) //first trial of first direction as example
	trials = DimSize(Matrix, 1) //returns the size of each dimension (x, y, z) 
	
	Make/O/N=(numpnts(dirs)) xpts, ypts
	
	WAVE Matrix = $(prefix)
	Make/O/N=(trials) $(prefix+"_theta"), $(prefix+"_DSis")
	Make/O/N=(2,trials) $(prefix+"_DSiArrows"), $(prefix+"_thetaArrows")
	WAVE thetas = $(prefix+"_thetas")
	WAVE DSis = $(prefix+"_DSis")
	WAVE thetaArrows = $(prefix+"_thetaArrows")
	WAVE DSiArrows = $(prefix+"_DSiArrows")
	for(j = 0; j < trials; j += 1)
		MatrixOp/O col = col(Matrix, j)
		xpts = col*cos(dirs*pi/180) //convert each value to x and y  coordinates
		ypts = col*sin(dirs*pi/180)
		xsum = sum(xpts) //sum across all 8 directions
		ysum = sum(ypts) 
		radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
		theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
		if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
			thetas[j] = 360 + theta
		else
			thetas[j] = theta
		endif
		thetaArrows[0][j] = theta
		thetaArrows[1][j] = theta
		DSis[j] = radius/sum(col)
		if(displayTrials)
			DSiArrows[1][j] =  DSis[j] * 60//wavemax(Matrix)
		else
			DSiArrows[1][j] =  DSis[j] * 60
	endif
	endfor	
	
	//SEM of DSis
	//Wavestats/Q DSis
	//print "DSi SEM: " + num2str(V_sem)
	//Wavestats/Q thetas //can't do this on 360 degrees thetas (need +/-)
	//print "theta SEM: " + num2str(V_sem)
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~AVG BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Make/O/N=1 $("avgtheta"), $("avgDSi")
	WAVE avgTheta = $("avgtheta")
	WAVE avgDSi = $("avgDSi")
	Make/O/N=2 $("avgthetaArrow"), $("avgDSiArrow")
	WAVE avgThetaArrow = $("avgthetaArrow")
	WAVE avgDSiArrow = $("avgDSiArrow")
	Make/O/N=(numpnts(dirs)) $("dirAvg"), $("dirSDev")
	Make/O/N=(numpnts(dirs)) $("dirSEM")
	WAVE dirAvg = $("dirAvg")
	WAVE dirSDev = $("dirSDev")
	WAVE dirSEM = $("dirSEM")
	for(j = 0; j < trials; j += 1)
		MatrixOp/O col = col(Matrix,j)
		dirAvg += col
	endfor
	dirAvg /= trials
	
	xpts = dirAvg*cos(dirs*pi/180) //convert each value to x and y  coordinates
	ypts = dirAvg*sin(dirs*pi/180)
	xsum = sum(xpts) //sum across all 8 directions
	ysum = sum(ypts) 
	radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
	theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
	if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
		avgtheta[0] = 360 + theta
	else
		avgtheta[0] = theta
	endif
	avgthetaArrow[0] = theta
	avgthetaArrow[1]= theta
	avgDSi[0] = radius/sum(dirAvg)
	if(displayTrials)
		avgDSiArrow[1] =  avgDSi[0] * 60//wavemax(Matrix)
	else
		avgDSiArrow[1] =  avgDSi[0] * 60
	endif
	
	for(j = 0;j < numpnts(dirs); j+= 1)
		MatrixOp/O row = row(Matrix,j)
		Wavestats/Q row
		dirSDev[j] = V_sdev
		dirSEM[j] = V_sem
	endfor
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~END AVG BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//sort angle and corresponding responses to be clockwise
	Duplicate/O dirs dirSort 
	Duplicate/O dirAvg $("dirAvgSort")
	WAVE dirAvgSort = $("dirAvgSort")
	Duplicate/O dirSDev $("dirSDevSort")
	WAVE dirSDevSort = $("dirSDevSort")
	Duplicate/O dirSEM $("dirSEMSort")
	WAVE dirSEMSort = $("dirSEMSort")
	Duplicate/O Matrix $(prefix + "_srt")
	WAVE sortMatrix = $(prefix + "_srt")
	SortColumns keyWaves=dirSort, sortWaves=sortMatrix//sort angleSort ascending and sort dirDFs to angleSort
	Sort dirSort dirSort, dirAvgSort, dirSDevSort, dirSEMSort
	
	Duplicate/O dirAvgSort $("dirAvgSortMinSEM"), $("dirAvgSortPluSEM")
	WAVE dirAvgSortMinSEM = $("dirAvgSortMinSEM")
	WAVE dirAvgSortPluSEM = $("dirAvgSortPluSEM")
	dirAvgSortMinSEM -= dirSEMsort
	dirAvgSortPluSEM += dirSEMsort
	Duplicate/O dirAvgSort $("dirAvgSortMinSD"), $("dirAvgSortPluSD")
	WAVE dirAvgSortMinSD = $("dirAvgSortMinSD")
	WAVE dirAvgSortPluSD = $("dirAvgSortPluSD")
	dirAvgSortMinSD -= dirSDevSort
	dirAvgSortPluSD += dirSDevSort
	
	// close the loops
	dirSort[numpnts(dirs)] = {dirSort[0]}
	dirAvgSort[numpnts(dirs)] = {dirAvgSort[0]}	
	dirSDevSort[numpnts(dirs)] = {dirSDevSort[0]}
	dirSEMSort[numpnts(dirs)] = {dirSEMSort[0]}	
	dirAvgSortMinSEM[numpnts(dirs)] = {dirAvgSortMinSEM[0]}	
	dirAvgSortPluSEM[numpnts(dirs)] = {dirAvgSortPluSEM[0]}	
	dirAvgSortMinSD[numpnts(dirs)] = {dirAvgSortMinSD[0]}	
	dirAvgSortPluSD[numpnts(dirs)] = {dirAvgSortPluSD[0]}	
	InsertPoints/M=0 numpnts(dirs),1,sortMatrix
	Duplicate/O sortMatrix dirSortMat
	for(j = 0; j < trials;j+=1)
		sortMatrix[numpnts(dirs)][j] = sortMatrix[0][j]
		dirSortMat[][j] = dirSort[p] //also make matrix of sorted dirs
	endfor
	
	if(displayTrials)
		//Build polar plots
		graphname = prefix + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, sortMatrix, dirSortMat, 360)	
		WMPolarAppendTrace(graphname, DSiArrows, thetaArrows, 360)
		ModifyGraph rgb = (52428,52428,52428)
		WMPolarAppendTrace(graphname, dirAvgSort, dirSort, 360)	
		WMPolarAppendTrace(graphname, avgDSiArrow, avgthetaArrow, 360)
		//WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		WMPolarSetManualRadiusRange(graphname, 0, 60)
	else
		graphname = prefix + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, DSiArrows, thetaArrows, 360)
		WMPolarAppendTrace(graphname, dirAvgSortMinSEM, dirSort, 360)
		WMPolarAppendTrace(graphname, dirAvgSortPluSEM, dirSort, 360)
		//WMPolarAppendTrace(graphname, dirAvgSortMinSD, dirSort, 360)
		//WMPolarAppendTrace(graphname, dirAvgSortPluSD, dirSort, 360)
		ModifyGraph rgb = (52428,52428,52428)
		WMPolarAppendTrace(graphname, dirAvgSort, dirSort, 360)	
		WMPolarAppendTrace(graphname, avgDSiArrow, avgthetaArrow, 360)
		//WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		WMPolarSetManualRadiusRange(graphname, 0, 60)
	endif
	Killwaves/Z col,row, xpts, ypts
End

Function rebuildMat(prefix)
	String prefix
	
	Variable i
	Killwaves/Z spikeTimes
	Make/N=(22,160)/O spikeTimes 
	
	for(i = 0; i < 160; i+=1)
		WAVE temp = $(prefix + num2str(i))
		spikeTimes[][i] = temp[p]
	endfor
End

Function alignment(matName)
	String matName
	
	WAVE Matrix = $matName
	Variable numDirs = 8
	variable numTrials = DimSize(Matrix ,1)/numDirs
	Variable i, j
	Make/O offsetDirs = {0,1,2,3,4,5,6,7}//{1,2,3,4,5}
	KillWaves/Z avgFirstSpk
	Make/O/N=(numpnts(offsetDirs)) avgFirstSpk
	Make/O/N=(numpnts(offsetDirs)) firstSpk
	
	for(i = 0; i < numTrials; i+=1)
		for (j = 0; j < numpnts(offsetDirs); j+=1)
	 		if (Matrix[0][i*numDirs+offsetDirs[j]] < firstSpk[j] || !i)
	 			firstSpk[j] = Matrix[0][i*numDirs+offsetDirs[j]]
	 		endif 
	 		avgFirstSpk[j] += Matrix[0][i*numDirs+offsetDirs[j]]
		endfor
	endfor
	avgFirstSpk /= numTrials
End

Function dsTimeEvolution(matName, binSize, numBins)
	String matName
	Variable binSize, numBins
	
	WAVE Matrix = $matName
	Variable numDirs = 8
	variable numTrials = DimSize(Matrix ,1)/numDirs
	Variable maxSpks = DimSize(Matrix, 0)
	Variable i, j, k, m, xsum, ysum, radius
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	//get proper null dir times by "SR" experiment
	//these are from firstSpike from each offset dir, the avg for the nulls
	//Make/O start = {78.86,63.6,89,98.3,82,61.4,78.86,78.86} //original
	//Make/O start = {27.5,49.2,85.7,98.3,82,58.1,37.8,50.7} //SR first spike (orig)
	Make/O start = {47,49.2,85.7,98.3,82,58.1,37.8,50.7} //SR first spike (mod dir0, ignore rogue spike)
	//Make/O start = {60.225,70.58,100.77,111.78,99.935,70.525,61.515,71.22} //SR AVG first spike (too late)
	
	KillWaves/Z bins
	Make/O/N=(numBins,numDirs) bins 
	//spike bins for each direction
	for(i = 0; i < numTrials; i+=1)
		for (j = 0; j < numDirs; j+=1)
	 		for(k = 0; k < numBins; k+=1)
	 			for(m = 0; m < maxSpks; m+=1)
	 				if (Matrix[m][i*numDirs+j] > start[j]-.1+binSize*k && Matrix[m][i*numDirs+j] < start[j]-.1+binSize*(k+1))
	 					bins[k][j] += 1
	 				endif
	 			endfor
	 		endfor
		endfor
	endfor
	
	Make/O/N=(numBins) binDSis
	Make/O/N=(numDirs) xpts, ypts
	for(i = 0; i < numBins; i+=1)
		MatrixOP/O temp = Row(bins,i)^t
		Duplicate/O temp $("bin"+num2str(binSize)+"_"+num2str(i+1))
		
		xpts = temp*cos(dirs*pi/180) //convert each value to x and y  coordinates
		ypts = temp*sin(dirs*pi/180)
		xsum = sum(xpts) //sum across all 8 directions
		ysum = sum(ypts) 
		radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
		
		binDSis[i] = radius/sum(temp)
	endfor
	Duplicate/O binDSis $("bin"+num2str(binSize)+"_DSis")
	
	KillWaves/Z temp, binDSis, xpts, ypts, dirSort
End

Function varianceEvolution(matName, binSize, numBins)
	String matName
	Variable binSize, numBins
	
	WAVE Matrix = $matName
	Variable numDirs = 8
	variable numTrials = DimSize(Matrix ,1)/numDirs
	Variable maxSpks = DimSize(Matrix, 0)
	Variable i, j, k, m, xsum, ysum, radius, theta
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	//get proper null dir times by "SR" experiment
	//these are from firstSpike from each offset dir, the avg for the nulls
	//Make/O start = {78.86,63.6,89,98.3,82,61.4,78.86,78.86} //original
	//Make/O start = {27.5,49.2,85.7,98.3,82,58.1,37.8,50.7} //SR first spike (orig)
	Make/O start = {47,49.2,85.7,98.3,82,58.1,37.8,50.7} //SR first spike (mod dir0, ignore rogue spike)
	//Make/O start = {60.225,70.58,100.77,111.78,99.935,70.525,61.515,71.22} //SR AVG first spike (too late)
	
	KillWaves/Z bins
	Make/O/N=(numBins,numDirs, numTrials) bins 
	//spike bins for each direction
	for(i = 0; i < numTrials; i+=1)
		for (j = 0; j < numDirs; j+=1)
	 		for(k = 0; k < numBins; k+=1)
	 			for(m = 0; m < maxSpks; m+=1)
	 				if (Matrix[m][i*numDirs+j] > start[j]-.1+binSize*k && Matrix[m][i*numDirs+j] < start[j]-.1+binSize*(k+1))
	 					bins[k][j][i] += 1
	 				endif
	 			endfor
	 		endfor
		endfor
	endfor
	
	Make/O/N=(numBins, numTrials) binDSis, binThetas
	Make/O/N=(numDirs) xpts, ypts
	for (i = 0; i < numTrials; i += 1)
		MatrixOP/O slice = layer(bins, i)
		for(j = 0; j < numBins; j+=1)
			MatrixOP/O temp = Row(slice,j)^t
			//Duplicate/O temp $("bin"+num2str(binSize)+"_"+num2str(j+1))
			xpts = temp*cos(dirs*pi/180) //convert each value to x and y  coordinates
			ypts = temp*sin(dirs*pi/180)
			xsum = sum(xpts) //sum across all 8 directions
			ysum = sum(ypts) 
			radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
			theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
	
			binThetas[j][i] = theta
			binDSis[j][i] = radius/sum(temp)
		endfor
	endfor
	Duplicate/O binDSis $("bin"+num2str(binSize)+"_DSis")
	Duplicate/O binThetas $("bin"+num2str(binSize)+"_Thetas")
	
	KillWaves/Z temp, slice, binDSis, binThetas, xpts, ypts
End

Function binTuningPolars(prefix,binSize,numBins)
	String prefix
	Variable binSize, numBins
	
	Variable i, j, k, m, xsum, ysum, radius, theta
	Variable manualMax, maxRadius
	manualMax = 1
	maxRadius = 120
	String graphname
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	Make/O/N=(numpnts(dirs)) xpts, ypts
	Make/O/N=(2,numBins) $(prefix+num2str(binSize)+"_DSiArrows"), $(prefix+num2str(binSize)+"_thetaArrows")
	WAVE thetaArrows = $(prefix+num2str(binSize)+"_thetaArrows")
	WAVE DSiArrows = $(prefix+num2str(binSize)+"_DSiArrows")
	
	for (i = 0; i < numBins; i += 1)
		WAVE binSpikes = $(prefix+num2str(binSize)+"_"+num2str(i+1))
		
		xpts = binSpikes*cos(dirs*pi/180) //convert each value to x and y  coordinates
		ypts = binSpikes*sin(dirs*pi/180)
		xsum = sum(xpts) //sum across all 8 directions
		ysum = sum(ypts) 
		radius = sqrt(xsum^2+ysum^2) //magnitude of the vector sum
		theta  = atan2(ysum,xsum)*180/pi //angle of the vector sum
		
		if (theta < 0) //above calculations come out as -180 to 180, convert these to values out of 360
			thetaArrows[0][i] = 360 + theta
			thetaArrows[1][i] = 360 + theta
		else
			thetaArrows[0][i] = theta
			thetaArrows[1][i] = theta
		endif
		
		if (manualMax)
			DSiArrows[1][i] = maxRadius * radius/sum(binSpikes)
		else
			DSiArrows[1][i] = wavemax(binSpikes)*radius/sum(binSpikes)
		endif
	endfor
	
	KillWaves/Z dirsCirc
	Duplicate/O dirs dirsCirc
	dirsCirc[numpnts(dirs)] = {dirs[0]}
	for (i = 0; i < numBins; i += 1)
		
		WAVE binSpikes = $(prefix+num2str(binSize)+"_"+num2str(i+1))
		Duplicate/O binSpikes $("bin"+num2str(binSize)+"SpikesCirc_"+num2str(i+1))
		WAVE binSpikesCirc = $("bin"+num2str(binSize)+"SpikesCirc_"+num2str(i+1))
		binSpikesCirc[numpnts(dirs)] = {binSpikesCirc[0]}
		
		MatrixOp/O dsiArr = Col(DSiArrows,i)
		MatrixOp/O thetaArr = Col(thetaArrows,i)
		Duplicate/O dsiArr $("dsiArr"+num2str(binSize)+"_"+num2str(i+1))
		Duplicate/O thetaArr $("thetaArr"+num2str(binSize)+"_"+num2str(i+1))
		WAVE dsiArrow = $("dsiArr"+num2str(binSize)+"_"+num2str(i+1))
		WAVE thetaArrow = $("thetaArr"+num2str(binSize)+"_"+num2str(i+1))
		
		graphname = prefix+num2str(binSize)+ "_" + num2str(i+1) + "_Polar"
		WMNewPolarGraph("_default_", graphname)
		WMPolarAppendTrace(graphname, dsiArrow, thetaArrow, 360)
		WMPolarAppendTrace(graphname,binSpikesCirc, dirsCirc, 360)
		//ModifyGraph rgb = (52428,52428,52428)
		//WMPolarSetZeroAngleWhere(graphname, "right", radiusOrigin = -.01) //so zeros do not distort shape
		if (manualMax)
			WMPolarSetManualRadiusRange(graphname, 0, maxRadius)
		endif
	endfor
	
	KillWaves/Z temp, binDSis, xpts, ypts, dirSort, dsiArr, thetaArr
End

Function trialSpikes(prefix, trials)
	String prefix
	Wave trials
		
	Variable i, j, k, perSpX, perX, perY, perSpY, count
	String botAxName, leftAxName, rightAxName, pathToWave
	
	Make/O/N=4501 timeAx = x*.1
	//directions are [225 270 315 0 45 90 135 180]
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	perSpX = .02
	perX = (1-perSpx*7)/8
	perSpY = .03
	perY = (1-perSpY*(numpnts(trials)-1))/numpnts(trials)
	
	WAVE Matrix = $(prefix)
	
	for(i = 0; i < numpnts(trials); i+=1)
		leftAxName = "t"+num2str(i)	
		for(j = 0; j < 8; j+=1)
			botAxName = "d"+num2str(j)
			//pull Vm recording of current direction
			//MatrixOp/O vm = col(Matrix,(trials[i]*numpnts(dirs)+j))
			
			if(!i && !j)
				//Display/L=$leftAxName/B=$botAxName vm vs timeAx
				Display/L=$leftAxName/B=$botAxName Matrix[][(trials[i]*numpnts(dirs)+j)] vs timeAx
			else
				//AppendToGraph/L=$leftAxName/B=$botAxName vm vs timeAx
				AppendToGraph/L=$leftAxName/B=$botAxName Matrix[][(trials[i]*numpnts(dirs)+j)] vs timeAx
			endif
			ModifyGraph axisEnab($botAxName)={(j*(perX+perSpX)),(j*(perX+perSpX)+perX)}
			ModifyGraph freePos($botAxName)=0
		endfor
		ModifyGraph axisEnab($leftAxName)={(i*(perY+perSpY)),(i*(perY+perSpY)+perY)}
		ModifyGraph freePos($leftAxName)=0
	endfor
	Killwaves/Z temp
End

Function chargePlot(numTrials)
	Variable numTrials
	Variable i, j, count
	
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	WAVE balExcMat = $("rho9_trialExc")
	WAVE balInhMat = $("rho9_trialInh")
	WAVE unbalExcMat = $("rho0_trialExc")
	WAVE unbalInhMat = $("rho0_trialInh")
	//Make/n=(371,8)/o rho9_avgExc, rho9_avgInh, rho0_avgExc, rho0_avgInh
	
	Make/n=(numTrials,8)/o rho9_excCharge, rho9_inhCharge, rho0_excCharge, rho0_inhCharge
	Make/n=(numTrials,8)/o rho9_ratioCharge, rho0_ratioCharge
	Make/n=(numTrials*8)/o rho9_ratioScatter, rho0_ratioScatter, dirScatter
	
	count = 0
	for(j = 0; j < numTrials; j += 1)
		for(i = 0; i < 8; i+=1)
			//unbalanced
			MatrixOp/O inh = col(unbalInhMat,(j*8+i))
			rho0_inhCharge[j][i] = sum(inh) // inh charge
			MatrixOp/O exc = col(unbalExcMat,(j*8+i))
			rho0_excCharge[j][i] = sum(exc) //exc charge
			rho0_ratioCharge[j][i] = abs(sum(exc)/sum(inh)) //ratio
			rho0_ratioScatter[count] = rho0_ratioCharge[j][i]
			//balanced
			MatrixOp/O inh = col(balInhMat,(j*8+i))
			rho9_inhCharge[j][i] = sum(inh) // inh charge
			MatrixOp/O exc = col(balExcMat,(j*8+i))
			rho9_excCharge[j][i] = sum(exc) //exc charge
			rho9_ratioCharge[j][i] = abs(sum(exc)/sum(inh))
			rho9_ratioScatter[count] = rho9_ratioCharge[j][i]
			
			dirScatter[count] = dirs[i]
			count += 1
		endfor
	endfor
	
	//use ratioCharge matrices to calc avg and SEM
	Killwaves/Z rho9_ratioAvg, rho0_ratioAvg, rho9_ratioError, rho0_ratioError
	Make/N=8/O rho9_ratioAvg, rho0_ratioAvg, rho9_ratioError, rho0_ratioError 
	for(j = 0; j < numTrials; j += 1)
		MatrixOp/O temp = row(rho0_ratioCharge, j)
		rho0_ratioAvg += temp
		MatrixOp/O temp = row(rho9_ratioCharge, j)
		rho9_ratioAvg += temp
	endfor
	rho0_ratioAvg /= numTrials
	rho9_ratioAvg /= numTrials
	//finish this with the sem
	for(i = 0; i < 8; i += 1)
		MatrixOp/O temp = col(rho0_ratioCharge, i)
		Wavestats/Q temp
		rho0_ratioError[i] = V_SEM
		MatrixOp/O temp = col(rho9_ratioCharge, i)
		Wavestats/Q temp
		rho9_ratioError[i] = V_SEM
	endfor
	
	Display
	for(j = 0; j < numTrials; j += 1)
		AppendToGraph rho0_ratioCharge[j][] vs dirs
		AppendToGraph rho9_ratioCharge[j][] vs dirs
	endfor
	
	Display rho0_ratioScatter vs dirScatter
	AppendToGraph rho9_ratioScatter vs dirScatter
	
	Killwaves/Z exc,inh
End

Function spikePlot(numTrials)
	Variable numTrials
	Variable i, j, count
	
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	WAVE balSpkMat = $("rho9_spks")
	WAVE unbalSpkMat = $("rho0_spks")
	
	Make/n=(numTrials*8)/O rho9_spkScatter, rho0_spkScatter, dirSpkScatter
	
	count = 0
	for(j = 0; j < numTrials; j += 1)
		for(i = 0; i < 8; i+=1)
			rho0_spkScatter[count] = unbalSpkMat[i][j]
			rho9_spkScatter[count] = balSpkMat[i][j]
			
			dirSpkScatter[count] = dirs[i]
			count += 1
		endfor
	endfor
	
	Killwaves/Z rho9_spkAvg, rho0_spkAvg
	Make/N=8/O rho9_spkAvg, rho0_spkAvg 
	for(j=0; j < numTrials; j += 1)
		// summing for average
		MatrixOp/O temp = col(unbalSpkMat,j)
		rho0_spkAvg += temp
		MatrixOp/O temp = col(balSpkMat,j)
		rho9_spkAvg += temp
	endfor
	rho9_spkAvg /= numTrials
	rho0_spkAvg /= numTrials
	
	Killwaves/Z rho9_spkError, rho0_spkError
	Make/N=8/O rho9_spkError, rho0_spkError
	for(i=0; i < 8; i += 1)
		// error bar calculations
		MatrixOp/O temp = row(unbalSpkMat,i)
		Wavestats/Q temp
		rho0_spkError[i] = V_SEM
		MatrixOp/O temp = row(balSpkMat,i)
		Wavestats/Q temp
		rho9_spkError[i] = V_SEM
	endfor
	
	Display
	for(j = 0; j < numTrials; j += 1)
		AppendToGraph unbalSpkMat[][j] vs dirs
		AppendToGraph balSpkMat[][j] vs dirs
	endfor
	
	Display rho0_spkScatter vs dirSpkScatter
	AppendToGraph rho9_spkScatter vs dirSpkScatter
	
	Killwaves/Z temp
End

Function treeRecAmpVar(mode)
	Variable mode //0:Vm; 1:iCa
	
	Variable i, j, k, m, count, wmin
	Variable xsum, ysum, theta, radius
	String pathToWave
	
	Variable numDirs = 8
	Variable numSegs = 700
	
	Make/O/N=1 trials
	count = 0
	for(i = 0; i < 30; i += 1)
		pathToWave = "trial" + num2str(i)
		if(exists(pathToWave))
			trials[count] = {i}
			count += 1
		endif
	endfor
	
	Make/O dirs = {-135, -90 , -45 , 0, 45, 90, 135, 180}
	
	// desired outputs: avg peak voltage of all sites, avg variance of peak voltage of all sites
	// avg variance of voltage traces of all sites  (these are all for all trials, can lump)
	// any benefit of dividing up on site by site basis, or just pool everything?
	
	KillWaves/Z peaks, peakPools, varPools, sitePeakVars, dirPeakVarAvgs
	KillWaves/Z dirPeakAvgs, dirVarAvgs, dirPeakVarSEMs, dirPeakSEMs, dirVarSEMs
	KillWaves/Z dirPeakVarMeanNormAvgs, dirPeakVarMeanNormSEMs
	Make/O/N=(numpnts(trials)) peaks
	Make/O/N=(numDirs, numSegs*numpnts(trials)) peakPools, varPools
	Make/O/N=(numDirs, numSegs) sitePeakVars, sitePeakVarMeanNorms
	Make/O/N=(numDirs) dirPeakVarAvgs, dirPeakAvgs, dirVarAvgs
	Make/O/N=(numDirs) dirPeakVarSEMs, dirPeakSEMs, dirVarSEMs
	Make/O/N=(numDirs) dirPeakVarMeanNormAvgs, dirPeakVarMeanNormSEMs
	for(k = 0; k < numDirs; k += 1)
		count = 0
		for(j = 0; j < numSegs; j+=1)
			for(i = 0; i < numpnts(trials); i += 1)
		
				WAVE Matrix = $("trial" + num2str(trials[i]))
		
				//pull a seg recording from current direction
				MatrixOp/O seg = col(Matrix,(k*numSegs+j))
				
				if(mode) // iCa
					seg = seg*-1 // flip the negative current
				endif //leave as is for Vm
				
				Wavestats/Q seg
				peaks[i] = V_max
				peakPools[k][count] = V_max
				varPools[k][count] = V_sdev
			
				count += 1
			endfor
			Wavestats/Q peaks
			sitePeakVars[k][j] = V_sdev
			sitePeakVarMeanNorms[k][j] = V_sdev / (V_avg + 63)
		endfor
		//peaks
		MatrixOp/O temp = row(peakPools,k)
		Wavestats/Q temp
		dirPeakAvgs[k] = V_avg
		dirPeakSEMs[k] = V_sem
		//variance
		MatrixOp/O temp = row(varPools,k)
		Wavestats/Q temp
		dirVarAvgs[k] = V_avg
		dirVarSEMs[k] = V_sem
		//peak variances
		MatrixOp/O temp = row(sitePeakVars,k)
		Wavestats/Q temp
		dirPeakVarAvgs[k] = V_avg
		dirPeakVarSEMs[k] = V_sem
		// mean normalized peak variances
		MatrixOp/O temp = row(sitePeakVarMeanNorms,k)
		Wavestats/Q temp
		dirPeakVarMeanNormAvgs[k] = V_avg
		dirPeakVarMeanNormSEMs[k] = V_sem
	endfor
	
	MatrixOp/O nullPeaks = row(peakPools,7)
	
	KillWaves/Z seg, temp, peaks
End