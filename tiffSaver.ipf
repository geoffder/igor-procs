#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function tiffSaver(first, numScans)
	Variable first, numScans
	
	Variable i
	String pad, toPath
	
	PathInfo home
	toPath = S_path
	
	for(i=first; i < first+numScans; i++)
		if (i < 10)
			pad = "00"
		elseif (i < 100)
			pad = "0"
		else
			pad = ""
		endif
		WAVE scan = $(GetDataFolder(1)+"Scan_"+pad+num2str(i)+":Scan_"+pad+num2str(i)+"_ch1")
		ImageSave/T="tiff"/U/S/O scan toPath+"Scan_"+pad+num2str(i)+"_ch1.tif"
	
	endfor
End
