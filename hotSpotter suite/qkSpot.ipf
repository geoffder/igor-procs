#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//quick online analysis for picking out DS cells in Gcamp retinas
//Operates on matrices in the standard folder tree of 2P_LSM
//parent folder is now twoP_Scans (rather than Nidaq_scans)
//indicate channel (1), prefix ("Scan"), firstScan (0), and numTrials (1)
Function qkSpotOnline(channel, prefix, firstScan, numTrials)
	String prefix
	Variable channel, firstScan, numTrials
	
	Variable i, j, z, currentScan
	String zeroPads
	
	Variable xWidth, yHeight, frames
	Variable filterSize = 7 //n x n size of filter
	Variable preFilter = 1
	Variable useDarkSub = 1 //use dark subtracted values for dF calc
	
	//fudge baseline of pixels that survive darkF mask (cells)
	//this serves to fix dF/F0 issues arrising from near 0 F0s
	Variable bslnInflate = 1
	Variable inflation = 10 
	
	//Make/O angle = {0, 180, 45, 225, 90, 270, 135, 315}
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45}
	
	Variable bslnStart, bslnEnd, peakStart, peakEnd
	bslnStart = 10
	bslnEnd = 18
	peakStart = 23
	peakEnd = 35
	
	Variable first = 1
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
			
			WAVE Matrix = $("root:twoP_Scans:"+prefix+"_"+zeroPads+num2str(currentScan)+":"+prefix+"_"+zeroPads+num2str(currentScan)+"_ch"+num2str(channel)) 
			
			if(first)
				xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
				yHeight = DimSize(Matrix, 1)
				frames = DimSize(Matrix, 2)
				Make/O/N=(xWidth,yHeight) tempMatrix
				first = 0
			endif
			
			Duplicate/O/RMD=[][][bslnStart,bslnEnd] Matrix bsln
			Duplicate/O/RMD=[][][peakStart,peakEnd] Matrix peak
			
			if(preFilter)
				for(z = 0; z <= bslnEnd-bslnStart; z+=1)
					tempMatrix[][] = bsln[p][q][z]
					MatrixFilter/N=(filterSize) avg tempMatrix
					bsln[][][z] = tempMatrix[p][q]
				endfor
				for(z = 0; z <= peakEnd-peakStart; z+=1)
					tempMatrix[][] = peak[p][q][z]
					MatrixFilter/N=(filterSize) avg tempMatrix
					peak[][][z] = tempMatrix[p][q]
				endfor
			endif
			
			//subtract global F avg (dark and cells)
			MatrixOp/O avgF = mean(bsln)
			MatrixOp/O bslnSub = bsln - avgF[0][0][0]
			MatrixOp/O peakSub = peak - avgF[0][0][0]
			
			//collapse stacks of bsln and peak regions to 2D
			MatrixOp/O bsln2D =  sumBeams(bsln)/(bslnEnd-bslnStart)
			MatrixOp/O peak2D =  sumBeams(peak)/(peakEnd-peakStart)
			MatrixOp/O bsln2Dsub =  sumBeams(bslnSub)/(bslnEnd-bslnStart)
			MatrixOp/O peak2Dsub =  sumBeams(peakSub)/(peakEnd-peakStart)
			
			//non-cells will be negative from above subtraction
			//replace all negative values with 0 (acts as mask)
			MatrixOp/O negs = sgn(bsln2Dsub)
			MatrixOp/O negs = replace(negs,-1,0)
			MatrixOp/O bsln2D = bsln2D * negs
			MatrixOp/O bsln2Dsub = bsln2Dsub * negs
			if(bslnInflate) //increase base F for cells only
				MatrixOp/O negs = replace(negs,1,inflation)
				MatrixOp/O bsln2Dsub = bsln2Dsub + negs
			endif
			
			//ditto for peak matrix
			MatrixOp/O negs = sgn(peak2Dsub)
			MatrixOp/O negs = replace(negs,-1,0)
			MatrixOp/O peak2D = peak2D * negs
			MatrixOp/O peak2Dsub = peak2Dsub * negs
			if(bslnInflate) //increase base F for cells only
				MatrixOp/O negs = replace(negs,1,inflation)
				MatrixOp/O peak2Dsub = peak2Dsub + negs
			endif
			
			//dividing is creating problems because Fs for cells are too low
			//not much greater than background (so dF/F0 gives #s in 1000s)
			//bslnInflate added above to address this
			if(useDarkSub)
				MatrixOp/O dF = (peak2Dsub - bsln2Dsub)/(bsln2Dsub)
			else
				MatrixOp/O dF = (peak2D - bsln2D)/(bsln2D)
			endif	
			MatrixOp/O negs = sgn(df)
			MatrixOp/O negs = replace(negs,-1,0)
			MatrixOp/O dF = dF * negs
			
			//stacks to store 2D matrices for each direction
			if(!j)
				Make/O/N=(xWidth,yHeight,8) dFstack, xStack, yStack
				Make/O/N=(xWidth,yHeight,8) peak2Dstack, bsln2Dstack
			endif
			
			dFstack[][][j] = dF[p][q]
			//begin vector calculations
			MatrixOp/O tempMatrix = dF * cos(angle[j]*pi/180)
			xStack[][][j] = tempMatrix[p][q]
			MatrixOp/O tempMatrix = dF * sin(angle[j]*pi/180)
			yStack[][][j] = tempMatrix[p][q]
			
			//even if using non-sub matrices for math, store these
			//because the numbers are more intuitive to look at
			bsln2Dstack[][][j] = bsln2Dsub[p][q]
			peak2Dstack[][][j] = peak2Dsub[p][q]
			
			currentScan += 1
		endfor
		
		//calculate theta and DSi from X and Y vectors
		MatrixOp/O x2D =  sumBeams(xStack)
		MatrixOp/O y2D =  sumBeams(yStack)
		MatrixOp/O theta = atan2(y2D, x2D)*180/pi
		MatrixOp/O theta = replace(theta,0,NaN)
		x2D = x2D^2
		y2D = y2D^2
		MatrixOp/O radius = sqrt(x2D + y2D)
		MatrixOp/O DSi = radius / sumBeams(dFStack)
		
		//save for trial
		Duplicate/O bsln2Dstack $("bsln2Dstack_" + num2str(i))
		Duplicate/O peak2Dstack $("peak2Dstack_" + num2str(i))
		Duplicate/O dFstack $("dFstack_" + num2str(i))
		Duplicate/O theta $("theta_" + num2str(i))
		Duplicate/O DSi $("DSi_" + num2str(i))
		
		//cleanup
		Killwaves/Z bsln,peak,dF,x2D,y2D,radius,bsln2D,peak2D
		KillWaves/Z DSi, theta, peakSub,bslnSub,negs,avgF
		KillWaves/Z bsln2Dsub, peak2Dsub
	endfor
	//cleanup
	KillWaves/Z xStack,yStack,tempMatrix,workMat
	KillWaves/Z dFstack,bsln2Dstack,peak2Dstack
End

Function qkSpot() : ButtonControl
	
	Variable count, i, j, z, k
	String pathToMatrix
	
	//get values from panel
	CONTROLINFO/W=spot matPrefix
	String prefix = S_value
	CONTROLINFO/W=spot preFilter
	Variable preFon = V_value
	CONTROLINFO/W=spot postFilter
	Variable postFon = V_value
	CONTROLINFO/W=spot preFilterType
	String preFtp = S_value
	CONTROLINFO/W=spot postFilterType
	String postFtp = S_value
	CONTROLINFO/W=spot preFilterSz
	Variable preFsz = V_value
	CONTROLINFO/W=spot postFilterSz
	Variable postFsz = V_value
	//CONTROLINFO/W=spot nzThreshold
	//Variable nzThresh = V_value
	CONTROLINFO/W=spot bslnStart
	Variable b1 = V_value
	CONTROLINFO/W=spot bslnEnd
	Variable b2 = V_value
	CONTROLINFO/W=spot peakStart
	Variable p1 = V_value
	CONTROLINFO/W=spot peakEnd
	Variable p2 = V_value
	CONTROLINFO/W=spot prefTheta
	Variable prefTheta = V_value
	CONTROLINFO/W=spot subSlider
	Variable userSubtract = V_value
	CONTROLINFO/W=spot inflSlider
	Variable inflation = V_value
	
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
	
	Variable xWidth, yHeight, frames
	Variable useDarkSub = 1 //use dark subtracted values for dF calc
	Variable saveExtras = 0 //if 1 make bsln and peak matrixces as well
	
	//get dimensions and setup matrices
	WAVE Matrix = $(GetDataFolder(1) + prefix + "0_" + num2str(trials[0]))
	xWidth = DimSize(Matrix, 0) //DimSize returns the size of each dimension (x, y, z) 
	yHeight = DimSize(Matrix, 1)
	frames = DimSize(Matrix, 2)
	Make/O/N=(xWidth,yHeight) tempMatrix
	Make/O/N=(xWidth,yHeight,numpnts(trials)+1) DSi2Dstack, theta2Dstack
	//avg movies for each direction
	for(k = 0; k < 8; k += 1)
		Make/O/N=(xWidth,yHeight,frames) $(prefix+num2str(k)+"_AVG")
	endfor
				
	//Make/O angle = {0, 180, 45, 225, 90, 270, 135, 315}
	Make/O angle = {90, 270, 180, 0, 135, 315, 225, 45}
	
	for(i = 0; i <= numpnts(trials); i += 1)
		if(i == numpnts(trials))
			//avg movies for each direction
			if(numpnts(trials)>1)
				for(k = 0; k < 8; k += 1)
					WAVE dirAVG = $(prefix+num2str(k)+"_AVG")
					dirAVG /= numpnts(trials)
				endfor
			endif
		endif
		
		for(j = 0; j < 8; j += 1)
			
			if(i != numpnts(trials) && numpnts(trials) > 1)
				WAVE Matrix = $(GetDataFolder(1) + prefix + num2str(j) + "_" + num2str(trials[i]))
				WAVE dirAVG = $(prefix+num2str(j)+"_AVG")
				dirAVG += Matrix
			elseif(numpnts(trials) > 1)
				WAVE Matrix = $(prefix+num2str(j)+"_AVG")
			endif
			
			Duplicate/O/RMD=[][][b1,b2] Matrix bsln
			Duplicate/O/RMD=[][][p1,p2] Matrix peak
			
			if(preFon)
				for(z = 0; z <= b2-b1; z+=1)
					tempMatrix[][] = bsln[p][q][z]
					if(StringMatch(preFtp,"avg"))
						MatrixFilter/N=(preFsz) avg tempMatrix
					elseif(StringMatch(preFtp,"median"))
						MatrixFilter/N=(preFsz) median tempMatrix
					elseif(StringMatch(preFtp,"gauss"))
						MatrixFilter/N=(preFsz) gauss tempMatrix
					endif
					bsln[][][z] = tempMatrix[p][q]
				endfor
				for(z = 0; z <= p2-p1; z+=1)
					tempMatrix[][] = peak[p][q][z]
					if(StringMatch(preFtp,"avg"))
						MatrixFilter/N=(preFsz) avg tempMatrix
					elseif(StringMatch(preFtp,"median"))
						MatrixFilter/N=(preFsz) median tempMatrix
					elseif(StringMatch(preFtp,"gauss"))
						MatrixFilter/N=(preFsz) gauss tempMatrix
					endif
					peak[][][z] = tempMatrix[p][q]
				endfor
			endif
			
			//subtract global F avg (dark and cells)
			Duplicate/O/RMD=[xWidth*.25,xWidth*.75][][] bsln bslnCrop
			MatrixOp/O temp = mean(bslnCrop)
			MatrixOp/O avgF = sumbeams(temp)/(b2-b1)
			//print avgF
			NVAR bkgF = root:bkgF
			bkgF = avgF[0]
			MatrixOp/O bslnSub = bsln - avgF[0] + avgF[0]*userSubtract
			MatrixOp/O peakSub = peak - avgF[0] + avgF[0]*userSubtract
			
			//collapse stacks of bsln and peak regions to 2D
			MatrixOp/O bsln2D =  sumBeams(bsln)/(b2-b1)
			MatrixOp/O peak2D =  sumBeams(peak)/(p2-p1)
			MatrixOp/O bsln2Dsub =  sumBeams(bslnSub)/(b2-b1)
			MatrixOp/O peak2Dsub =  sumBeams(peakSub)/(p2-p1)
			
			//non-cells will be negative from above subtraction
			//replace all negative values with 0 (acts as mask)
			MatrixOp/O bsln2Dsub = replace(bsln2Dsub,0,-1)
			MatrixOp/O negs = sgn(bsln2Dsub)
			//make peaks corresponding to neg bslns neg also (zero later)
			MatrixOp/O peak2D = peak2D * negs
			MatrixOp/O peak2Dsub = peak2Dsub * negs
			//zero negs in bsln
			MatrixOp/O negs = replace(negs,-1,0)
			MatrixOp/O bsln2D = bsln2D * negs
			MatrixOp/O bsln2Dsub = bsln2Dsub * negs
			//NaN out the zeros
			MatrixOp/O bsln2D = replace(bsln2D,0,NaN)
			MatrixOp/O bsln2Dsub = replace(bsln2Dsub,0,NaN)
			//fudge baseline of pixels that survive darkF mask (cells)
			//this serves to fix dF/F0 issues arrising from near 0 F0s
			MatrixOp/O negs = replace(negs,1,inflation*avgF[0][0][0])
			MatrixOp/O bsln2Dsub = bsln2Dsub + negs
			
			//ditto for peak matrix
			MatrixOp/O peak2Dsub = replace(peak2Dsub,0,-1)
			MatrixOp/O negs = sgn(peak2Dsub)
			MatrixOp/O negs = replace(negs,-1,0)
			MatrixOp/O peak2D = peak2D * negs
			MatrixOp/O peak2Dsub = peak2Dsub * negs
			//NaN out the zeros
			MatrixOp/O peak2D = replace(peak2D,0,NaN)
			MatrixOp/O peak2Dsub = replace(peak2Dsub,0,NaN)
			//inflation
			MatrixOp/O negs = replace(negs,1,inflation*avgF[0][0][0])
			MatrixOp/O peak2Dsub = peak2Dsub + negs
			
			//dividing is creating problems because Fs for cells are too low
			//not much greater than background (so dF/F0 gives #s in 1000s)
			//bslnInflate added above to address this
			if(useDarkSub)
				MatrixOp/O dF = (peak2Dsub - bsln2Dsub)/(bsln2Dsub)
			else
				MatrixOp/O dF = (peak2D - bsln2D)/(bsln2D)
			endif	
			MatrixOp/O negs = sgn(df)
			MatrixOp/O negs = replace(negs,-1,0)
			MatrixOp/O dF = dF * negs
			
			//stacks to store 2D matrices for each direction
			if(!j)
				Make/O/N=(xWidth,yHeight,8) dFstack, xStack, yStack
				Make/O/N=(xWidth,yHeight,8) peak2Dstack, bsln2Dstack
			endif
			
			dFstack[][][j] = dF[p][q]
			//begin vector calculations
			MatrixOp/O tempMatrix = dF * cos(angle[j]*pi/180)
			xStack[][][j] = tempMatrix[p][q]
			MatrixOp/O tempMatrix = dF * sin(angle[j]*pi/180)
			yStack[][][j] = tempMatrix[p][q]
			
			//even if using non-sub matrices for math, store these
			//because the numbers are more intuitive to look at
			bsln2Dstack[][][j] = bsln2Dsub[p][q]
			peak2Dstack[][][j] = peak2Dsub[p][q]
		endfor
		
		//calculate theta and DSi from X and Y vectors
		MatrixOp/O x2D =  sumBeams(xStack)
		MatrixOp/O y2D =  sumBeams(yStack)
		MatrixOp/O theta = atan2(y2D, x2D)*180/pi
		MatrixOp/O theta = replace(theta,0,NaN)
		//make proper relative theta
		MatrixOp/O negs = sgn(theta)
		MatrixOp/O negs = replace(negs,-1,360)
		MatrixOp/O negs = replace(negs,1,0)
		MatrixOp/O temp = theta + negs
		MatrixOp/O theta = temp - prefTheta
		//now -360 if above 180 and +360 if below -180
		Duplicate/O theta theta2
		MatrixOp/O temp = theta - 180
		MatrixOp/O negs = sgn(temp)
		MatrixOp/O negs = replace(negs,1,-180)
		MatrixOp/O negs = replace(negs,-1,180)
		MatrixOp/O theta = temp + negs
		MatrixOp/O temp = theta2 + 180
		MatrixOp/O negs = sgn(temp)
		MatrixOp/O negs = replace(negs,-1,180)
		MatrixOp/O negs = replace(negs,1,-180)
		MatrixOp/O theta = temp + negs
		//Ben snippet to try later:
		//temp = (temp < 0) ? temp+180 : temp
		
		x2D = x2D^2
		y2D = y2D^2
		MatrixOp/O radius = sqrt(x2D + y2D)
		MatrixOp/O DSi = radius / sumBeams(dFStack)
		
		if(postFon)
			if(StringMatch(postFtp,"avg"))
				MatrixFilter/N=(postFsz) avg theta
				MatrixFilter/N=(postFsz) avg DSi
			elseif(StringMatch(postFtp,"median"))
				MatrixFilter/N=(postFsz) median theta
				MatrixFilter/N=(postFsz) median DSi
			elseif(StringMatch(postFtp,"gauss"))
				MatrixFilter/N=(postFsz) gauss theta
				MatrixFilter/N=(postFsz) gauss DSi
			endif
		endif
		//store theta and DSi in stacks for easy viewing
		//avg is last layer of stack (after trials)
		theta2Dstack[][][i] = theta[p][q]
		DSi2Dstack[][][i] = DSi[p][q]
		//save for trial
		if(i != numpnts(trials))
			if(saveExtras)
				Duplicate/O bsln2Dstack $("bsln2Dstack_" + num2str(trials[i]))
				Duplicate/O peak2Dstack $("peak2Dstack_" + num2str(trials[i]))
			endif
			Duplicate/O dFstack $("dFstack_" + num2str(trials[i]))
			Duplicate/O theta $("theta_" + num2str(trials[i]))
			Duplicate/O DSi $("DSi_" + num2str(trials[i]))
		elseif(numpnts(trials) > 1)
			if(saveExtras)
				Duplicate/O bsln2Dstack $("bsln2Dstack_AVG")
				Duplicate/O peak2Dstack $("peak2Dstack_AVG")
			endif
			Duplicate/O dFstack $("dFstack_AVG")
			Duplicate/O theta $("theta_AVG")
			Duplicate/O DSi $("DSi_AVG")
		endif
		//cleanup
		Killwaves/Z bsln,peak,dF,x2D,y2D,radius,bsln2D,peak2D
		KillWaves/Z DSi, theta, peakSub,bslnSub,negs//,avgF
		KillWaves/Z bsln2Dsub, peak2Dsub
	endfor

	//Display or Update theta and DSi plots
	//DoWindow relativeTheta
	if(!V_flag)
		//Display as "relativeTheta"
		//Slider preFilterSz win=relativeTheta,pos={35,52},size={150,35},value=0,ticks=numpnts(trials),vert=0,limits={0,numpnts(trials),1},proc=thetaSlider
		//AppendToGraph theta2DStack 
	endif
	//DoWindow DSi
	if(!V_flag)
		//Display DSi2DStack as "DSi"need to display as image check the code for this
		//Slider preFilterSz win=DSi,pos={35,52},size={150,35},value=0,ticks=numpnts(trials),vert=0,limits={0,numpnts(trials),1},proc=DSiSlider
	endif
	//cleanup
	KillWaves/Z xStack,yStack,tempMatrix,workMat
	KillWaves/Z dFstack,bsln2Dstack,peak2Dstack
End

Function spotPanel()
	Variable/G root:bkgF 
	//make panel unless it already exists
	DoWindow spot
	if(V_flag)
		KillWindow spot
	endif
	newPanel/N=spot/W=(200,200,400,500)
	//matrix prefix and preferred theta
	SetVariable matPrefix title="prefix",pos={14,9},value= _STR:"dir",size={80,14},limits={-inf,inf,0}//,proc=setvarProc
	SetVariable prefTheta title="prefTheta", pos={99,9}, value=_NUM:0,limits={0,360,.1},size={80,14}//,proc=setvarProc
	//pre-filter
	SetDrawEnv fsize= 10
	DrawText 7,44,"Pre-filter:"
	SetDrawEnv fsize= 10
	DrawText 7.5,77,"Size:\r(nxn)"
	PopupMenu preFilterType pos={85,32},value="avg;median;gauss"//,proc=popMenuProc
	CheckBox preFilter pos={57,33},size={61,44},side=1,title="on",value=1//,proc=CheckProc
	Slider preFilterSz win=spot, pos={35,52}, size={150,35},value=3,ticks=10,vert=0,limits={3,21,2}//,proc=sliderProc
	//post-filter
	SetDrawEnv fsize= 10
	DrawText 7,103.5,"Post-filter:"
	SetDrawEnv fsize= 10
	DrawText 8.5,132,"Size:\r(nxn)"
	PopupMenu postFilterType pos={85,92},value="avg;median;gauss "
	CheckBox postFilter pos={57,93},size={24,12},side=1,title="on",value=0//,proc=CheckProc,value=0
	Slider postFilterSz win=spot, pos={35,110}, size={150,35},value=3,ticks=10,vert=0,limits={3,21,2}//,proc=sliderProc
	//baseline and response times
	SetVariable bslnStart,pos={14,149}, title="bslnStart",size={79,14},value= _NUM:3,limits={0,inf,1}//,proc=setvarProc
	SetVariable bslnEnd title="bslnEnd ", pos={14,172}, size={79,14}, value= _NUM:10,limits={0,inf,1}//,proc=setvarProc
	SetVariable peakStart pos={101,146}, title="peakStart",size={79,14},value= _NUM:20,limits={0,inf,1}//,proc=setvarProc
	SetVariable peakEnd pos={101,173}, title="peakEnd ", size={79,14}, value= _NUM:25,limits={0,inf,1}//,proc=setvarProc
	//ssss
	SetDrawEnv fsize= 10
	DrawText 8,229,"adjust\rbkgr\rsbtrc"
	Slider subSlider win=spot, pos={44,196}, size={150,35},value=0,ticks=15,vert=0,limits={-7,7,.00005}
	SetDrawEnv fsize= 10
	DrawText 8,270,"incr\rdend\rFluor"
	Slider inflSlider win=spot, pos={44,238}, size={150,35},value=0,ticks=15,vert=0,limits={-7,7,.00005}
	//set number of deviations response must be above baseline
	//(removed for now, would slow down program to implement)
	//SetDrawEnv fsize= 10
	//DrawText 8,245,"Noise Threshold:"
	//CheckBox nzThreshCheck pos={65,249},size={24,12},side=1,title="on",value=0//,proc=CheckProc,value=0
	//SetDrawEnv fsize= 10
	//DrawText 8,229,"noise \rthresh\r(n*SD)"
	//Slider nzThreshold win=spot, pos={44,196}, size={150,35},value=1,ticks=10,vert=0,limits={0,5,.25}//,proc=sliderProc
	//run qkSpot with values from panel
	Button runQkSpot title="qkSpot( )",size={50,20},pos={145,275},proc=ButtonProc
	
End

Function spotPanelV()
	Variable/G root:bkgF 
	//make panel unless it already exists
	DoWindow spot
	if(V_flag)
		KillWindow spot
	endif
	newPanel/N=spot/W=(200,200,600,700)
	//matrix prefix and preferred theta
	SetVariable matPrefix title="prefix",pos={14,9},value= _STR:"dir",size={80,14},limits={-inf,inf,0}//,proc=setPrefixProc
	SetVariable prefTheta title="prefTheta", pos={99,9}, value=_NUM:0,limits={0,360,.1},size={120,20}//,proc=setvarProc
	//pre-filter
	SetDrawEnv fsize= 10
	DrawText 7,44,"Pre-filter:"
	SetDrawEnv fsize= 10
	DrawText 7.5,77,"Size:\r(nxn)"
	PopupMenu preFilterType pos={95,33},value="avg;median;gauss"//,proc=popMenuProc
	CheckBox preFilter pos={57,33},size={61,44},side=1,title="on",value=1//,proc=CheckProc
	Slider preFilterSz win=spot, pos={35,52}, size={150,35},value=3,ticks=10,vert=0,limits={3,21,2}//,proc=sliderProc
	//post-filter
	SetDrawEnv fsize= 10
	DrawText 7,120,"Post-filter:"
	SetDrawEnv fsize= 10
	DrawText 8.5,150,"Size:\r(nxn)"
	PopupMenu postFilterType pos={95,110},value="avg;median;gauss "
	CheckBox postFilter pos={57,110},size={24,12},side=1,title="on",value=0//,proc=CheckProc,value=0
	Slider postFilterSz win=spot, pos={35,125}, size={150,35},value=3,ticks=10,vert=0,limits={3,21,2}//,proc=sliderProc
	//baseline and response times
	SetVariable bslnStart,pos={14,180}, title="bslnStart",size={100,14},value= _NUM:3,limits={0,inf,1}//,proc=setvarProc
	SetVariable bslnEnd title="bslnEnd ", pos={14,200}, size={100,14}, value= _NUM:10,limits={0,inf,1}//,proc=setvarProc
	SetVariable peakStart pos={120,180}, title="peakStart",size={100,14},value= _NUM:20,limits={0,inf,1}//,proc=setvarProc
	SetVariable peakEnd pos={120,200}, title="peakEnd ", size={100,14}, value= _NUM:25,limits={0,inf,1}//,proc=setvarProc
	//ssss
	SetDrawEnv fsize= 10
	DrawText 8,280,"adjust\rbkgr\rsbtrc"
	Slider subSlider win=spot, pos={44,240}, size={150,35},value=0,ticks=15,vert=0,limits={-4,4,.0000005},proc=sliderProc
	ValDisplay subDisplay pos={199,243},title="actual\rvalue",size={80,20}
	SetDrawEnv fsize= 10
	DrawText 8,320,"incr\rdend\rFluor"
	Slider inflSlider win=spot, pos={44,300}, size={150,35},value=0,ticks=15,vert=0,limits={-4,4,.0000005},proc=sliderProc
	ValDisplay inflDisplay pos={199,303},title="actual\rvalue",size={80,20}
	//set number of deviations response must be above baseline
	//(removed for now, would slow down program to implement)
	//SetDrawEnv fsize= 10
	//DrawText 8,245,"Noise Threshold:"
	//CheckBox nzThreshCheck pos={65,249},size={24,12},side=1,title="on",value=0//,proc=CheckProc,value=0
	//SetDrawEnv fsize= 10
	//DrawText 8,229,"noise \rthresh\r(n*SD)"
	//Slider nzThreshold win=spot, pos={44,196}, size={150,35},value=1,ticks=10,vert=0,limits={0,5,.25}//,proc=sliderProc
	//run qkSpot with values from panel
	Button runQkSpot title="qkSpot( )",size={70,20},pos={100,350},proc=ButtonProc
End

Function sliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				//Variable curval = sa.curval
				NVAR bkgF = root:bkgF
				CONTROLINFO/W=spot subSlider
				Variable userSubtract = V_value
				CONTROLINFO/W=spot inflSlider
				Variable inflation = V_value
				ValDisplay subDisplay value=_NUM:(-bkgF+userSubtract*bkgF)
				ValDisplay inflDisplay value=_NUM:(inflation*bkgF)
			endif
			break
	endswitch
	return 0
End
Function DSiSlider(sa) : SliderControl
	STRUCT WMSliderAction &sa
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				//change displayed frame
			endif
			break
	endswitch
	return 0
End
Function thetaSlider(sa) : SliderControl
	STRUCT WMSliderAction &sa
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				//change displayed frame
			endif
			break
	endswitch
	return 0
End
Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
				//qkSpot()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End
Function SetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
				//qkSpot()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End
//Function SetPrefixProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
				
				String prefix = sa.curval
				WAVE Matrix = $(prefix+"0_1")
				CONTROLINFO/W=spot bslnStart
				Variable b1 = V_value
				CONTROLINFO/W=spot bslnEnd
				Variable b2 = V_value
				Duplicate/O/RMD=[][][b1,b2] Matrix bsln
				MatrixOp/O avgF = mean(bsln)
				KillWaves bsln
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End
Function CheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			//qkSpot()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			qkSpot()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End