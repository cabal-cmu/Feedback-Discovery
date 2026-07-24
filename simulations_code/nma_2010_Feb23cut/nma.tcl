# NMA - FMRIB's DCM Tool
#
# Mark Woolrich, Saad Jbabdi and Matthew Webster, FMRIB Image Analysis Group
#
# Copyright (C) 2008 University of Oxford
#
# TCLCOPYRIGHT

source [ file dirname [ info script ] ]/fslstart.tcl

set MYEXEC    [ info nameofexecutable ]
set MYSHELL   [ file tail $MYEXEC ]
if { [  string match -nocase *wish* $MYSHELL ] } {
option add *LabelEntry.e.background grey95
}

proc nma { w } {
    global FSLDIR NMA
    # ---- Set up Frames ----
    toplevel $w
    wm title $w "NMA MINE"
    wm iconname $w "NMA"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    frame $w.f
    frame $w.f.data
    frame $w.roi
    frame $w.subjects
    TitleFrame $w.f.basic -relief groove -text "Main Settings" -font {Helvetica 12 bold} 
    TitleFrame $w.model -relief groove -text "Models" -font {Helvetica 12 bold}

    set NMA(nROI) 1
    LabelSpinBox $w.f.nROI -label "Number of ROIs   " -textvariable NMA(nROI) -range " 1 100 1 "  -command "$w.f.nROI.spin.e validate; nmaROI:update $w 1" -modifycmd  "nmaROI:update $w 1"
    set NMA(nEV) 1
    LabelSpinBox $w.f.nEV -label "Number of EVs    " -textvariable NMA(nEV) -range " 1 100 1 "  -command "$w.f.nEV.spin.e validate; nmaEV:update $w 1" -modifycmd  "nmaEV:update $w 1"
    set NMA(nModel) 1
    LabelSpinBox $w.f.nModel -label "Number of models" -textvariable NMA(nModel) -range " 1 100 1 "  -command "$w.f.nModel.spin.e validate; nmaModel:update $w" -modifycmd  "nmaModel:update $w"
    set NMA(nSubjects) 1
    LabelSpinBox $w.f.nSubject -label "Number of subjects" -textvariable NMA(nSubjects) -range " 1 100 1 "  -command "$w.f.nSubjects.spin.e validate; nmaModel:update $w" -modifycmd  "nmaSubject:update $w"

    set NMA(tr) 3
    LabelSpinBox $w.f.tr -label "TR (s) " -textvariable NMA(tr) -range {0.0001 200000 0.25 }
    pack  $w.f.nROI $w.f.nEV $w.f.nModel $w.f.nSubject $w.f.tr -in [ $w.f.basic getframe ] -anchor w
    set NMA(outputDirectory) [ exec sh -c " pwd " ]/output.nma
    FileEntry  $w.f.outputDirectory -textvariable NMA(outputDirectory) -label "Output Directory" -title "Select" -width 50 -filedialog directory -filetypes { }

    set NMA(dataIsSingle) 0
    checkbutton $w.f.dataIsSingle -text "Inputs are single timeseries" -variable NMA(dataIsSingle) -command "nma:updateDataMode $w"
    set NMA(inferROI) 1
    checkbutton $w.f.inferROI -text "infer ROIs" -variable NMA(inferROI) 
    set NMA(standardImage) "$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"
    FileEntry  $w.f.standardImage -textvariable NMA(standardImage) -label "Standard Image" -title "Select" -width 50 -filedialog directory -filetypes IMAGE

    pack $w.f.dataIsSingle $w.f.inferROI $w.f.standardImage -in $w.f.data -anchor w
    pack $w.f.outputDirectory $w.f.data -in [ $w.f.basic getframe ] -anchor w


    NoteBook $w.model.setupModel -side top -bd 2 -tabpady {5 10} -arcradius 3 
    nmaModel:update $w

    NoteBook $w.subjects.setupSubject -side top -bd 2 -tabpady {5 10} -arcradius 3
    nmaSubject:update $w


    pack $w.model.setupModel -in [ $w.model getframe ]
    pack $w.f.basic $w.subjects.setupSubject $w.f $w.model $w.subjects -expand yes -fill both

    frame $w.controlButtons

    button $w.controlButtons.go -command "nma:apply" -text "Go" -width 5
    button $w.controlButtons.save -command "feat_file:setup_dialog $w a a a [namespace current] *.nsf {Save NMA setup} {nma:save $w} {}" -text "Save" -width 5
    button $w.controlButtons.load -command "feat_file:setup_dialog $w a a a [namespace current] *.nsf {Load NMA setup} {nma:load $w} {}" -text "Load" -width 5
    button $w.controlButtons.cancel -command "destroy $w" -text "Exit" -width 5

    pack $w.controlButtons.go $w.controlButtons.cancel $w.controlButtons.save $w.controlButtons.load -side left
    pack $w.controlButtons

    collapsible frame $w.advanced -title "Advanced Options"
    set NMA(nBurnIn) 10000
    LabelSpinBox $w.nBurnIn -label "Number of discarded samples  " -textvariable NMA(nBurnIn) -range " 0 10000000 1 " 
    set NMA(nSamples) 5000
    LabelSpinBox $w.nSamples -label "Number of MCMC samples  " -textvariable NMA(nSamples) -range " 0 10000000 1 "
    set NMA(highPassSeconds) 100
    LabelSpinBox $w.nHighPass -label "High pass filter cutoff (s)" -textvariable NMA(highPassSeconds) -range " 0 10000000 1 "
    set NMA(nJumps) 1
    LabelSpinBox $w.nJumps -label "Jumps per Sample" -textvariable NMA(nJumps) -range " 0 10000000 1 "
    set NMA(resolution) 1
    LabelSpinBox $w.resolution -label "Resolution Factor" -textvariable NMA(resolution) -range " 0 100 0.1"

    pack  $w.nBurnIn $w.nSamples $w.nHighPass $w.nJumps  $w.resolution -in $w.advanced.b -anchor w -pady 2
    pack $w.advanced 
}

proc nma:updateDataMode { w } {
global NMA
    pack forget $w.f.standardImage
    for { set subject 1 } { $subject <= $NMA(nSubjects) } { incr subject 1 } {
	pack forget $w.subjects.inputData${subject} $w.subjects.maskData${subject} $w.subjects.transform${subject} -in $w.subjects.data -anchor w
    }
    if { ! $NMA(dataIsSingle) } {
	pack $w.f.standardImage -in $w.f.data -anchor w
	for { set subject 1 } { $subject <= $NMA(nSubjects) } { incr subject 1 } {
	    pack  $w.subjects.inputData${subject} $w.subjects.maskData${subject} $w.subjects.transform${subject} -in $w.subjects.data${subject} -anchor w
	}
    }
}

proc nmaSubject:update { w } {
global NMA
    #Delete any notebook pages beyond those wanted
    for { set i [ llength [ $w.subjects.setupSubject pages ] ]  } { $i > $NMA(nSubjects) } { incr i -1 } {
	nmaSubject:destroy $w $i
    }
    #Add notebook pages to fill number wanted 
    for { set i [ expr [ llength [ $w.subjects.setupSubject pages ] ] + 1 ] } { $i <= $NMA(nSubjects) } { incr i 1 } {
	nmaSubject:create $w $i
    }
    #Has the active page been deleted?
    if { [ $w.subjects.setupSubject raise ] == "" } {
	$w.subjects.setupSubject raise Subject$NMA(nSubjects) 
    } 
}

proc nmaSubject:create { w subject } {
global NMA
    $w.subjects.setupSubject insert $subject Subject${subject} -text "Subject ${subject}"

    frame $w.subjects.evSubject${subject}
    NoteBook $w.subjects.evSubject${subject}.setupEV -side top -bd 2 -tabpady {5 10} -arcradius 3
    nmaEV:update $w $subject
    pack $w.subjects.evSubject${subject}.setupEV 

    frame $w.subjects.roiSubject${subject}
    NoteBook $w.subjects.roiSubject${subject}.setupROI -side top -bd 2 -tabpady {5 10} -arcradius 3 
    nmaROI:update $w $subject
    pack $w.subjects.roiSubject${subject}.setupROI

    if { ! [ winfo exists $w.subjects.data${subject} ] } {
	frame $w.subjects.data${subject}
	FileEntry  $w.subjects.inputData${subject} -textvariable NMA(inputData,$subject) -label "Input Timeseries" -command "nma:updateOutputName" -title "Select" -width 50 -filedialog directory -filetypes IMAGE
	FileEntry  $w.subjects.maskData${subject} -textvariable NMA(maskData,$subject) -label "Subject Mask" -title "Select" -width 50 -filedialog directory -filetypes IMAGE
	FileEntry  $w.subjects.transform${subject} -textvariable NMA(transform,$subject) -label "Standard to example transform" -title "Select" -width 50 -filedialog directory -filetypes *
	pack  $w.subjects.inputData${subject} $w.subjects.maskData${subject} $w.subjects.transform${subject} -in $w.subjects.data${subject} -anchor w -pady 2

	FileEntry $w.subjects.confounds${subject} -textvariable NMA(confounds,$subject) -label "Confounds" -title "Select" -width 50 -filedialog directory -filetypes *
    }
    pack $w.subjects.evSubject${subject} $w.subjects.roiSubject${subject} $w.subjects.data${subject} $w.subjects.confounds${subject} -in [ $w.subjects.setupSubject getframe Subject${subject} ] -anchor w -pady 5

   $w.subjects.setupSubject compute_size
    
   if { $subject > 1 } {
       nmaSubject:copy $w [ expr $subject - 1 ] $subject
   }

}

proc nma:updateOutputName { { name output } } {
global FSLDIR NMA
    set NMA(outputDirectory)  [ fsl:exec "$FSLDIR/bin/remove_ext $name" -n ].nma
}

proc nmaSubject:destroy { w subject } {
    if { $subject <= [ llength [ $w.subjects.setupSubject pages ] ] } {
	$w.subjects.setupSubject delete Subject${subject}
	destroy $w.subjects.evSubject${subject} 
	destroy $w.subjects.roiSubject${subject}
    }
}

proc nmaSubject:copy { w sourceSubject destinationSubject } {
    global NMA
    set NMA(inputData,$destinationSubject) $NMA(inputData,$sourceSubject)
    set NMA(maskData,$destinationSubject) $NMA(maskData,$sourceSubject)
    set NMA(confounds,$destinationSubject) $NMA(confounds,$sourceSubject)
    set NMA(transform,$destinationSubject) $NMA(transform,$sourceSubject)
    for { set i 1 } { $i <= $NMA(nEV) } { incr i 1 } {
       set NMA(EVname,$i,$destinationSubject) $NMA(EVname,$i,$sourceSubject)
       set NMA(EVfile,$i,$destinationSubject) $NMA(EVfile,$i,$sourceSubject)
       set NMA(EVmodulation,$i,$destinationSubject) {$NMA(EVmodulation,$i,$sourceSubject)}
	if { $NMA(EVmodulation,$i,$sourceSubject) == 1 } {
	    $w.subjects.evSubject${destinationSubject}.amplitude$i select
	} else {
	    $w.subjects.evSubject${destinationSubject}.amplitude$i deselect
	}
    }
    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
       set NMA(ROIname,$i,$destinationSubject) $NMA(ROIname,$i,$sourceSubject)
       set NMA(ROIfile,$i,$destinationSubject) $NMA(ROIfile,$i,$sourceSubject)
    }
}

proc nmaModel:update { w } {
global NMA
    #Delete any notebook pages beyond those wanted
    for { set i [ llength [ $w.model.setupModel pages ] ]  } { $i > $NMA(nModel) } { incr i -1 } {
	nmaModel:destroy $w $i
    }
    #Add notebook pages to fill number wanted 
    for { set i [ expr [ llength [ $w.model.setupModel pages ] ] + 1 ] } { $i <= $NMA(nModel) } { incr i 1 } {
	nmaModel:create $w $i
    }
    #Has the active page been deleted?
    if { [ $w.model.setupModel raise ] == "" } {
	$w.model.setupModel raise model$NMA(nModel) 
    } 
}

proc nmaModel:destroy { w i } {
    $w.model.setupModel delete model$i
}

proc nmaModel:create { w i } {
global NMA
    $w.model.setupModel insert $i model$i -text "Model $i"

    if { ! [ winfo exists $w.model.name$i ] } {
	if { ! [ info exists NMA(modelName,$i) ] } { set NMA(modelName,$i) "model$i" } 
	LabelEntry  $w.model.name$i -textvariable NMA(modelName,$i) -label "Model Name" -width 20
	frame $w.model.matrixButtons$i
	button $w.model.matrixButtons$i.setupMatrixA -command "nma:setupMatrixA $w $i" -text "Setup A Matrix" -width 9
	button $w.model.matrixButtons$i.setupMatrixB -command "nma:setupMatrixB $w $i" -text "Setup B Matrix" -width 9
	button $w.model.matrixButtons$i.setupMatrixC -command "nma:setupMatrixC $w $i" -text "Setup C Matrix" -width 9
	button $w.model.matrixButtons$i.setupMatrixD -command "nma:setupMatrixD $w $i" -text "Setup D Matrix" -width 9
    }
    pack $w.model.matrixButtons$i.setupMatrixA $w.model.matrixButtons$i.setupMatrixB  $w.model.matrixButtons$i.setupMatrixC  $w.model.matrixButtons$i.setupMatrixD -side left
    pack $w.model.name$i -in  [ $w.model.setupModel getframe model$i ] -anchor w
    pack $w.model.matrixButtons$i -in  [ $w.model.setupModel getframe model$i ] -anchor w
    $w.model.setupModel compute_size
}

proc nmaROI:update { w subject } {
global NMA
    #Delete any notebook pages beyond those wanted
    for { set i [ llength [ $w.subjects.roiSubject${subject}.setupROI pages ] ]  } { $i > $NMA(nROI) } { incr i -1 } {
	nmaROI:destroy $w $i $subject
    }
    #Add notebook pages to fill number wanted 
    for { set i [ expr [ llength [ $w.subjects.roiSubject${subject}.setupROI pages ] ] + 1 ] } { $i <= $NMA(nROI) } { incr i 1 } {
	nmaROI:create $w $i $subject
    }
    #Has the active page been deleted?
    if { [ $w.subjects.roiSubject${subject}.setupROI raise ] == "" } {
	$w.subjects.roiSubject${subject}.setupROI raise roi$NMA(nROI) 
    } 
}

proc nmaROI:destroy { w i subject } {
    $w.subjects.roiSubject${subject}.setupROI delete roi$i
}

proc nmaROI:create { w i subject } {
global NMA
    $w.subjects.roiSubject${subject}.setupROI insert $i roi$i -text "ROI $i"

    if { ! [ winfo exists $w.subjects.roiSubject${subject}.file$i ] } {
	if { ! [ info exists NMA(ROIname,$i,$subject) ] } { set NMA(ROIname,$i,$subject) "roi$i" } 
	LabelEntry  $w.subjects.roiSubject${subject}.name$i -textvariable NMA(ROIname,$i,$subject) -label "ROI Name" -width 20 -state disabled -disabledbackground lightgrey -disabledforeground black
    }
    if { ! [ winfo exists $w.subjects.roiSubject${subject}.file$i ] } {
	FileEntry  $w.subjects.roiSubject${subject}.file$i -textvariable NMA(ROIfile,$i,$subject) -label "ROI Image/Data" -command "nmaROI:updateName $i $subject" -title "Select" -width 50 -filedialog directory -filetypes * 
    }
    pack $w.subjects.roiSubject${subject}.name$i -in  [ $w.subjects.roiSubject${subject}.setupROI getframe roi$i ] -anchor w
    pack $w.subjects.roiSubject${subject}.file$i -in  [ $w.subjects.roiSubject${subject}.setupROI getframe roi$i ] -anchor w
    $w.subjects.roiSubject${subject}.setupROI compute_size
}

proc nmaROI:updateName { i subject { dummy dummy } } {
global FSLDIR NMA
    if { ! $NMA(dataIsSingle) } {
	set NMA(ROIname,$i,$subject) [ fsl:exec "basename `$FSLDIR/bin/remove_ext $NMA(ROIfile,$i,$subject)` _mask" -n ]
    } else {
	set NMA(ROIname,$i,$subject) [ fsl:exec "basename $NMA(ROIfile,$i,$subject) .txt" -n ]
    }
}

proc nmaEV:update { w subject } {
global NMA 
    #Delete any notebook pages beyond those wanted
    for { set i [ llength [ $w.subjects.evSubject${subject}.setupEV pages ] ]  } { $i > $NMA(nEV) } { incr i -1 } {
	nmaEV:destroy $w $i $subject
    }
    #Add notebook pages to fill number wanted 
    for { set i [ expr [ llength [ $w.subjects.evSubject${subject}.setupEV pages ] ] + 1 ] } { $i <= $NMA(nEV) } { incr i 1 } {
	nmaEV:create $w $i $subject 
    }
    #Has the active page been deleted?
    if { [ $w.subjects.evSubject${subject}.setupEV raise ] == "" } {
	$w.subjects.evSubject${subject}.setupEV raise ev$NMA(nEV) 
    } 
}

proc nmaEV:destroy { w i subject }  {
    $w.subjects.evSubject${subject}.setupEV delete ev$i
}

proc nmaEV:create { w i subject } {
global NMA
    $w.subjects.evSubject${subject}.setupEV insert $i ev$i -text "EV $i"

    if { ! [ winfo exists $w.subjects.evSubject${subject}.file$i ] } {
        if { ! [ info exists NMA(EVname,$i,$subject) ] } { set NMA(EVname,$i,$subject) "ev$i" }
	LabelEntry  $w.subjects.evSubject${subject}.name$i -textvariable NMA(EVname,$i,$subject) -label "EV Name" -width 20 -state disabled -disabledbackground lightgrey -disabledforeground black
    }
    if { ! [ winfo exists $w.subjects.evSubject${subject}.file$i ] } {
	FileEntry  $w.subjects.evSubject${subject}.file$i -textvariable NMA(EVfile,$i,$subject) -label "EV Timeseries File"  -command "nmaEV:updateName $i $subject" -title "Select" -width 50 -filedialog directory -filetypes *
    }
    if { ! [ winfo exists $w.subjects.evSubject${subject}.amplitude$i ] } {
	checkbutton $w.subjects.evSubject${subject}.amplitude$i -text "Model between-trial/epoch variance" -variable NMA(EVmodulation,$i,$subject)
    }
    pack $w.subjects.evSubject${subject}.name$i -in  [ $w.subjects.evSubject${subject}.setupEV getframe ev$i ] -anchor w
    pack $w.subjects.evSubject${subject}.file$i -in  [ $w.subjects.evSubject${subject}.setupEV getframe ev$i ] -anchor w
    pack $w.subjects.evSubject${subject}.amplitude$i -in  [ $w.subjects.evSubject${subject}.setupEV getframe ev$i ] -anchor w
    $w.subjects.evSubject${subject}.setupEV compute_size
}

proc nmaEV:updateName { i subject { dummy dummy } } {
global NMA
    set NMA(EVname,$i,$subject) [ fsl:exec "basename $NMA(EVfile,$i,$subject)" ]
    set NMA(EVname,$i,$subject) [ file rootname $NMA(EVname,$i,$subject) ]
}

proc nma:save { w fileName } {
global NMA
    set outputFile [ open $fileName w ]
    puts $outputFile "set NMA(outputDirectory) {$NMA(outputDirectory)}"
    puts $outputFile "set NMA(nSubjects) $NMA(nSubjects)"
    puts $outputFile "set NMA(dataIsSingle) $NMA(dataIsSingle)"
    puts $outputFile "set NMA(inferROI) $NMA(inferROI)"
    puts $outputFile "set NMA(standardImage) {$NMA(standardImage)}"
    puts $outputFile "set NMA(tr) $NMA(tr)"

    puts $outputFile "set NMA(nROI) $NMA(nROI)"
    puts $outputFile "set NMA(nEV) $NMA(nEV)"
    for { set subject 1 } { $subject <= $NMA(nSubjects) } { incr subject 1 } {
    puts $outputFile "set NMA(inputData,$subject) {$NMA(inputData,$subject)}"
    puts $outputFile "set NMA(maskData,$subject) {$NMA(maskData,$subject)}"
    puts $outputFile "set NMA(confounds,$subject) {$NMA(confounds,$subject)}"
    puts $outputFile "set NMA(transform,$subject) {$NMA(transform,$subject)}"
	for { set i 1 } { $i <= $NMA(nEV) } { incr i 1 } {
	    puts $outputFile "set NMA(EVname,$i,$subject) {$NMA(EVname,$i,$subject)}"
	    puts $outputFile "set NMA(EVfile,$i,$subject) {$NMA(EVfile,$i,$subject)}"
	    puts $outputFile "set NMA(EVmodulation,$i,$subject) {$NMA(EVmodulation,$i,$subject)}"
	}
	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	    puts $outputFile "set NMA(ROIname,$i,$subject) {$NMA(ROIname,$i,$subject)}"
	    puts $outputFile "set NMA(ROIfile,$i,$subject) {$NMA(ROIfile,$i,$subject)}"
	}
    }

    puts $outputFile "set NMA(nModel) $NMA(nModel)"
    for { set model 1 } { $model <= $NMA(nModel) } { incr model 1 } {
	puts $outputFile "set NMA(modelName,$model) {$NMA(modelName,$model)}"
	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	    for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
		if { [ info exists NMA($model,A,$i,$j) ] } { puts $outputFile "set NMA($model,A,$i,$j) {$NMA($model,A,$i,$j)}" } else {
		    puts $outputFile "set NMA($model,A,$i,$j) 0"
		}
	    }
	}
	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	    for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
		for { set k 1 } { $k <= $NMA(nEV) } { incr k 1 } {
		    if { [ info exists NMA($model,B,$i,$j,$k) ] } {
			puts $outputFile "set NMA($model,B,$i,$j,$k) {$NMA($model,B,$i,$j,$k)}"
		    } else {
			puts $outputFile "set NMA($model,B,$i,$j,$k) 0"
		    }
		}
	    }
	}
	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	    for { set j 1 } { $j <= $NMA(nEV) } { incr j 1 } {
		if { [ info exists NMA($model,C,$i,$j) ] } { puts $outputFile "set NMA($model,C,$i,$j) {$NMA($model,C,$i,$j)}" } else {
		    puts $outputFile "set NMA($model,C,$i,$j) 0"
		}
	    }
	}
	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	    for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
		for { set k 1 } { $k <= $NMA(nROI) } { incr k 1 } {
		    if { [ info exists NMA($model,D,$i,$j,$k) ] } {
			puts $outputFile "set NMA($model,D,$i,$j,$k) {$NMA($model,D,$i,$j,$k)}"
		    } else {
			puts $outputFile "set NMA($model,D,$i,$j,$k) 0"
		    }
		}
	    }
	}
    }

    close $outputFile
}

proc nma:load { w fileName { GUIisActive 1 } } {
global NMA
    source $fileName
    if { $GUIisActive == 1 } {
	nmaModel:update $w
	for { set subject $NMA(nSubjects) } { $subject >= 1 } { incr subject -1 } {
	    nmaSubject:destroy $w $subject
	}
	nmaSubject:update $w
    }
}

proc nma:apply { } {
global NMA FSLDIR FSLDEVDIR
#Create output directory
    fsl:exec "mkdir -p $NMA(outputDirectory)"
    nma:save 0 ${NMA(outputDirectory)}/design.nsf
    fsl:exec "$FSLDIR/bin/imcp $NMA(standardImage) $NMA(outputDirectory)/standard_brain"
    for { set subject 1 } { $subject <= $NMA(nSubjects) } { incr subject 1 } {
	#End of common files
	#NMA(dataIs4D)
	fsl:exec "mkdir -p  $NMA(outputDirectory)/subject${subject}"
	fsl:exec "$FSLDIR/bin/imcp $NMA(inputData,$subject) $NMA(outputDirectory)/subject${subject}/data"
	fsl:exec "$FSLDIR/bin/imcp $NMA(maskData,$subject) $NMA(outputDirectory)/subject${subject}/mask"
	fsl:exec "cp $NMA(transform,$subject) $NMA(outputDirectory)/subject${subject}/standard2example_func.mat"
	fsl:exec "cp $NMA(confounds,$subject) $NMA(outputDirectory)/subject${subject}/confound_evs.txt"

	for { set i 1 } { $i <= $NMA(nEV) } { incr i 1 } {
	    fsl:exec "cp $NMA(EVfile,$i,$subject) $NMA(outputDirectory)/subject${subject}"
	}

	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {	
	    fsl:exec "cp $NMA(ROIfile,$i,$subject) $NMA(outputDirectory)/subject${subject}" 
	} 

#modify NMA source to allow ROI in subject directory

	set EVnames "--stim="
	for { set i 1 } { $i <= $NMA(nEV) } { incr i 1 } {
	    if { $i != $NMA(nEV) } { set EVnames "$EVnames$NMA(EVname,$i,$subject),"
	    } else { set EVnames "$EVnames$NMA(EVname,$i,$subject)" }
	}
	set ROInames "--nn="
	for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	    if { $i != $NMA(nROI) } { set ROInames "$ROInames$NMA(ROIname,$i,$subject)," 
	    } else { set ROInames "$ROInames$NMA(ROIname,$i,$subject)" }
	} 
	set amplitudeSettings "--sam="
	for { set i 1 } { $i <= $NMA(nEV) } { incr i 1 } {
	    if { $i != $NMA(nEV) } { set amplitudeSettings "$amplitudeSettings$NMA(EVmodulation,$i,$subject)," 
	    } else { set amplitudeSettings "$amplitudeSettings$NMA(EVmodulation,$i,$subject)" }
	} 

	for { set model 1 } { $model <= $NMA(nModel) } { incr model 1 } {
	    fsl:exec "mkdir -p $NMA(outputDirectory)/$NMA(modelName,$model)"
	    nma:writeMatrixA $NMA(outputDirectory)/$NMA(modelName,$model) $model
	    nma:writeMatrixB $NMA(outputDirectory)/$NMA(modelName,$model) $model
	    nma:writeMatrixC $NMA(outputDirectory)/$NMA(modelName,$model) $model
	    nma:writeMatrixD $NMA(outputDirectory)/$NMA(modelName,$model) $model

	    set NMA(dataDirectory) $NMA(outputDirectory)
	    set NMA(modelDirectory) $NMA(outputDirectory)/$NMA(modelName,$model)
	    set theCommand "--ld=$NMA(outputDirectory)/subject${subject}/$NMA(modelName,$model)_results --dd=$NMA(dataDirectory) --md=$NMA(modelDirectory) --sub=subject${subject} --tr=$NMA(tr) --resfactor=16"
	    if { $NMA(inferROI) == 1 } {
		set theCommand "$theCommand --dm=1"
	    } else {
		set theCommand "$theCommand --dm=2"
	    }
	    set theCommand "$theCommand $EVnames $ROInames $amplitudeSettings --bi=10000 --ns=5000 --se=1 --scc --pev --hm=balloon"
	    if { $NMA(dataIsSingle) } {
		set theCommand "$theCommand --sts"
	    }
	    fsl:exec "$FSLDIR/bin/nma $theCommand"
	}
    }
}
#Matrix functions
proc nma:setupMatrixA { w model } {
global FSLDIR NMA
    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0
    wm iconname $w0 "$model:A Matrix"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm
    wm title $w0 "A Matrix"

    frame $w0.gridFrame -padx 5 -pady 5
    
    label $w0.source -text "From" 
    pack $w0.source
    frame $w0.matrix
    label $w0.to -text "To" 

    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	label $w0.xname$i -text $NMA(ROIname,$i,1)
	label $w0.yname$i -text $NMA(ROIname,$i,1)
	grid  $w0.xname$i -in  $w0.gridFrame -row 0 -column $i
	grid  $w0.yname$i -in  $w0.gridFrame -row $i -column 0
	for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
	    if { $i != $j } {
		checkbutton $w0.button$i$j -variable NMA($model,A,$i,$j)
		grid $w0.button$i$j -in  $w0.gridFrame -row $i -column $j
	    }
	}
    }
    pack $w0.to $w0.gridFrame -side left -in $w0.matrix
    button $w0.done -command "destroy $w0" -text "Done" -width 5
    pack $w0.matrix $w0.done
}

proc nma:setupMatrixB { w model } {
global FSLDIR NMA
    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0
    wm iconname $w0 "B Matrix"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm
    wm title $w0 "B Matrix"

    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	frame $w0.matrix$i
	label $w0.source$i -text "Modulator" 
	pack $w0.source$i -in  $w0.matrix$i
	label $w0.to$i -text "From" 
	frame $w0.gridFrame$i -padx 5 -pady 5
        label $w0.xlabel$i -text "To $NMA(ROIname,$i,1)" -fg white -bg black
	grid  $w0.xlabel$i -in  $w0.gridFrame$i -row 0 -column 0
	for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
	    if {  $i != $j } {
		label $w0.yname$i$j -text $NMA(ROIname,$j,1)
	        grid  $w0.yname$i$j -in  $w0.gridFrame$i -row $j -column 0
		for { set k 1 } { $k <= $NMA(nEV) } { incr k 1 } {
		    label $w0.xname$i$j$k -text $NMA(EVname,$k,1)
		    grid  $w0.xname$i$j$k -in  $w0.gridFrame$i -row 0 -column $k
		    checkbutton $w0.button$i$j$k -variable NMA($model,B,$i,$j,$k)
		    grid $w0.button$i$j$k -in  $w0.gridFrame$i -row $j -column $k
		}		
	    }
	}
	pack $w0.to$i $w0.gridFrame$i -side left -in  $w0.matrix$i
	pack $w0.matrix$i
    }
    button $w0.done -command "destroy $w0" -text "Done" -width 5
    pack $w0.done
}

proc nma:setupMatrixC { w model } {
global FSLDIR NMA
    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0
    wm iconname $w0 "C Matrix"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm
    wm title $w0 "C Matrix"

    frame $w0.gridFrame -padx 5 -pady 5
    label $w0.source -text "From" 
    pack $w0.source
    frame $w0.matrix
    label $w0.to -text "To" 

    for { set i 1 } { $i <= $NMA(nEV) } { incr i 1 } {
	label $w0.yname$i -text $NMA(EVname,$i,1)
	grid  $w0.yname$i -in  $w0.gridFrame -row 0 -column $i
    }

    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	label $w0.xname$i -text $NMA(ROIname,$i,1)
	grid  $w0.xname$i -in  $w0.gridFrame -row $i -column 0
	for { set j 1 } { $j <= $NMA(nEV) } { incr j 1 } {
	    checkbutton $w0.button$i$j -variable NMA($model,C,$i,$j)
	    grid $w0.button$i$j -in  $w0.gridFrame -row $i -column $j
	}
    }
    pack $w0.to $w0.gridFrame -side left -in $w0.matrix
    button $w0.done -command "destroy $w0" -text "Done" -width 5
    pack $w0.matrix $w0.done
}

proc nma:setupMatrixD { w model } {
global FSLDIR NMA
    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0
    wm iconname $w0 "D Matrix"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm
    wm title $w0 "D Matrix"


    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	frame $w0.matrix$i
	label $w0.source$i -text "Modulator" 
	pack $w0.source$i -in  $w0.matrix$i
	label $w0.to$i -text "From" 
	frame $w0.gridFrame$i -padx 5 -pady 5
        label $w0.xlabel$i -text "To $NMA(ROIname,$i,1)" -fg white -bg black
	grid  $w0.xlabel$i -in  $w0.gridFrame$i -row 0 -column 0
	for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
	    if {  $i != $j } {
		label $w0.yname$i$j -text $NMA(ROIname,$j,1)
	        grid  $w0.yname$i$j -in  $w0.gridFrame$i -row $j -column 0
		for { set k 1 } { $k <= $NMA(nROI) } { incr k 1 } {
		    label $w0.xname$i$j$k -text $NMA(ROIname,$k,1)
		    grid  $w0.xname$i$j$k -in  $w0.gridFrame$i -row 0 -column $k
		    checkbutton $w0.button$i$j$k -variable NMA($model,D,$i,$j,$k)
		    grid $w0.button$i$j$k -in  $w0.gridFrame$i -row $j -column $k
		}		
	    }
	}

	pack $w0.to$i $w0.gridFrame$i -side left -in  $w0.matrix$i
	pack $w0.matrix$i
    }
    button $w0.done -command "destroy $w0" -text "Done" -width 5
    pack $w0.done
}

proc nma:writeMatrixA { pathName model } {
global NMA
    set outputFile [open "$pathName/matA.txt" w]
    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
	    if { [ info exists NMA($model,A,$i,$j) ] } {
		puts -nonewline $outputFile $NMA($model,A,$i,$j)
		puts -nonewline $outputFile { }
	    } else {
		puts -nonewline $outputFile 0
		puts -nonewline $outputFile { }
	    }
	}
	puts $outputFile {}
    }
    close $outputFile
}

proc nma:writeMatrixB { pathName model } {
global NMA

    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
        set outputFile [open "$pathName/matB_into_$NMA(ROIname,$i,1).txt" w]
	for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
	    for { set k 1 } { $k <= $NMA(nEV) } { incr k 1 } {
		if { [ info exists NMA($model,B,$i,$j,$k) ] } {
		    puts -nonewline $outputFile $NMA($model,B,$i,$j,$k)
		    puts -nonewline $outputFile { }
		} else {
		    puts -nonewline $outputFile 0
		    puts -nonewline $outputFile { }
		}
	    }
	    puts $outputFile {}
	}
	close $outputFile
    }
}

proc nma:writeMatrixC { pathName model } {
global NMA
    set outputFile [open "$pathName/matC.txt" w]
    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
	for { set j 1 } { $j <= $NMA(nEV) } { incr j 1 } {
	    if { [ info exists NMA($model,C,$i,$j) ] } {
		puts -nonewline $outputFile $NMA($model,C,$i,$j)
		puts -nonewline $outputFile { }
	    } else {
		puts -nonewline $outputFile 0
		puts -nonewline $outputFile { }
	    }
	}
	puts $outputFile {}
    }
    close $outputFile
}

proc nma:writeMatrixD { pathName model } {
global NMA

    for { set i 1 } { $i <= $NMA(nROI) } { incr i 1 } {
        set outputFile [open "$pathName/matD_into_$NMA(ROIname,$i,1).txt" w]
	for { set j 1 } { $j <= $NMA(nROI) } { incr j 1 } {
	    for { set k 1 } { $k <= $NMA(nROI) } { incr k 1 } {
		if { [ info exists NMA($model,D,$i,$j,$k) ] } {
		    puts -nonewline $outputFile $NMA($model,D,$i,$j,$k)
		    puts -nonewline $outputFile { }
		} else {
		    puts -nonewline $outputFile 0
		    puts -nonewline $outputFile { }
		}
	    }
	    puts $outputFile {}
	}
	close $outputFile
    }
}

if { [  string match -nocase *wish* $MYSHELL ] } {
    wm withdraw .
    nma .rename
    tkwait window .rename
}
