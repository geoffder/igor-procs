#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//Continue button kills pause window
Function PauseForROI_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	
	DoWindow/K PauseForROI // Kill panel
End
