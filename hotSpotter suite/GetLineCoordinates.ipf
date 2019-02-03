#pragma rtGlobals=3		// Use modern global access method and strict wave access.
 
 // Based on PrintImageROI (http://www.igorexchange.com/node/1217) written by harneit
 
Function GetLineCoordinates(lineNumber)
	Variable lineNumber
	
// 1. get the marquee in point units (axis units directly accessible only for 1D graphs!)
// sanity test: let there be exactly one image
	String TopGraphImages = ImageNameList("",";")
	if( ItemsInList(TopGraphImages) != 1 )
		Print "There is either no image or more than one on the graph"
		return -1
	endif
	
	//get first and last points from current both x and y waves for the current line of interest
	//e.g. line1_x and line1_y
	
// 2. get the plot area in point units
	GetWindow kwTopWin, psize
	Variable ppL = V_left, ppR = V_right, ppT = V_top, ppB = V_bottom
//	print "plot area in points (left, top, right, bottom) =", ppL, ppT, ppR, ppB
 
// 3. get the axis limits
// get limits of axes that image is displayed against
	String theImage = StringFromList(0, TopGraphImages)
	String XAxis = StringByKey( "XAXIS", ImageInfo("", theImage, 0) )
	String YAxis = StringByKey( "YAXIS", ImageInfo("", theImage, 0) )
	String XAxisLimits = StringFromList(2, StringByKey( "SETAXISCMD", axisinfo("", XAxis) ), " " )
	String YAxisLimits = StringFromList(2, StringByKey( "SETAXISCMD", axisinfo("", YAxis) ), " " )
	
// order left/right correctly, remember whether a swap occurred
	Variable axMin, axMax
	if ( strlen(XAxisLimits) == 0 )		// i.e. SetAxis/A bottom
		axMin = DimOffset($theImage,0) - ( DimDelta($theImage,0) / 2 )
		axMax = axMin + DimDelta($theImage,0)*DimSize($theImage,0)
	else
		axMin = str2num( StringFromList(0, XAxisLimits, ",") )
		axMax = str2num( StringFromList(1, XAxisLimits, ",") )
	endif
	Variable paL = min(axMin, axMax), paR = max(axMin, axMax), xSwap = (paL != axMin)
 
// order top/bottom correctly, remember whether a swap occurred
	if ( strlen(YAxisLimits) == 0 )		// i.e. SetAxis/A left
		axMin = DimOffset($theImage,1) - ( DimDelta($theImage,1) / 2 )
		axMax = axMin + DimDelta($theImage,1)*DimSize($theImage,1)
	else
		axMin = str2num( StringFromList(0, YAxisLimits, ",") )
		axMax = str2num( StringFromList(1, YAxisLimits, ",") )
	endif
	Variable paT = max(axMin, axMax), paB = min(axMin, axMax), ySwap = (paB != axMin)

//	print "plot area in axis units (left, top, right, bottom) =", paL, paT, paR, paB; print xSwap, ySwap
 
// 4. express the marquee coordinates in axis units
	Variable maL = paL + (mpL - ppL) * (paR - paL)/(ppR - ppL), maR = paL + (mpR - ppL) * (paR - paL)/(ppR - ppL), swap
	if( xSwap )
		swap = maL; maL = maR; maR = swap
	endif
	Variable maT = paT + (mpT - ppT) * (paB - paT)/(ppB - ppT), maB = paT + (mpB - ppT) * (paB - paT)/(ppB - ppT)
	if( ySwap )
		swap = maT; maT = maB; maB = swap
	endif
	print "marquee in axis units (left, top, right, bottom) =", maL, maT, maR, maB
End
