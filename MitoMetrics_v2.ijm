/*	MitoMetrics_v2.ijm
	***********************
	Author: 			Winnok H. De Vos & Marlies Verschuuren
	Modified by: 		Marlies Verschuuren
	Date Created: 		Nov 24th, 2012 
	Date Last Modified:	June 26th, 2023 
	
	References: 
	* De Vos WH, Van Neste L, Dieriks B, Joss GH, Van Oostveldt P. High content image cytometry in the context of subnuclear organization. Cytometry A. 2010 Jan;77(1):64-75. doi: 10.1002/cyto.a.20807. PMID: 19821512.
	* De Puysseleyr L, De Puysseleyr K, Vanrompay D, De Vos WH. Quantifying the growth of chlamydia suis in cell culture using high-content microscopy. Microsc Res Tech. 2017 Apr;80(4):350-356. doi: 10.1002/jemt.22799. Epub 2016 Nov 12. PMID: 27862609.
	* Sieprath T, Corne T, Robijns J, Koopman WJH, De Vos WH. Cellular Redox Profiling Using High-content Microscopy. J Vis Exp. 2017 May 14;(123):55449. doi: 10.3791/55449. PMID: 28570523; PMCID: PMC5607941.
 	
 	Description: 
 	Macro Set dedicated to cytometric analysis of n-dimensional images containing at least one channel of a nuclear dye, as well as a cellular counterstain and a mitochondria marker. 
 	It returns nuclear, cellular and mitochondia shape descriptors and intensity parameters.
 	Build on Cellblocks_v19.ijm (https://github.com/DeVosLab/CellBlocks).
 	
	Requires FeatureJ and imagescience plugins (E. Meijering): download at http://www.imagescience.org/meijering/software/featurej/
	Also requires Stardist and CSBDeep plugins (overlapping nuclei segmentation): https://imagej.net/StarDist
	Also requires Trained_Cell_Segmentation.bsh plugin and a trained model (*.model)
	Also requires PTBIOP plugin for cellpose and a working cellpose environment: https://github.com/BIOP/ijl-utilities-wrappers and https://github.com/MouseLand/cellpose#local-installation
	Also requires Cellpose_Script_Wrapper.groovy plugin for cellpose
	
	Change Log
	+ 
 	_________________________________________________________________
*/

/*
 	***********************

	Variable initiation

	***********************
*/


//	String variables
var cells_results					= "";										//	cell region results
var cells_roi_set 					= "";										//	cell region ROIsets
var cells_segmentation_method		= "Threshold";								//  Methods used for cell detection
var cells_threshold					= "Li";									//	threshold method for segmentation of cells 
var dir								= "";										//	directory
var log_path						= "";										//	path for the log file
var micron							= getInfo("micrometer.abbreviation");		// 	micro symbol
var mito_results					= "";										//	mitochondria results
var mito_roi_set					= "";										//	mitochondria roi set
var mito_threshold_method			= "Triangle";									// 	mitochondria threshold method
var nuclei_roi_set 					= "";										//	nuclear ROIsets
var nuclei_results 					= "";										//	nuclear results
var nuclei_segmentation_method		= "Gauss";									// 	Method used for nuclei detection
var nuclei_threshold				= "Li";										//	threshold method for segmentation of nuclei 
var order							= "xyczt(default)";							//	hyperstack dimension order
var output_dir						= "";										//	dir for analysis output
var results							= "";										//	summarized results	

var suffix							= ".tif";									//	suffix for specifying the file type

//	Number variables
var channels						= 3;										//	number of channels
var cells_channel					= 2;										//	channel used for segmentation of cells 
var cells_diameter					= 100;										//  approximate diameter of cells 
var cells_filter_scale				= 1;										//	radius for cell smoothing
var cells_fixed_threshold_value		= 200;										//	fixed maximum threshold for cell segmentation (if auto doesn't work well);
var fields							= 4;										//	number of xy positions
var image_height					= 1000;										//	image height
var image_width						= 1000;										//	image width
var iterations						= 25;										// 	iterations of dilations in region growing (for determining cell boundaries)
var slices							= 7;										//	number of z-slices
var mito_channel					= 1;										//	mito channel
var mito_fixed_threshold_value		= 5;										// 	fixed threshold
var mito_minsize					= 2;										// 	min size
var mito_rolling					= 50;										//	rolling ball radius (µm) subtract background
var mito_scale						= 5;										//	scale laplace
var nuclei_channel					= 3;										//	channel used for segmentation of nuclei 
var nuclei_filter_scale				= 3;										// 	radius for nuclei smoothing/laplace
var nuclei_fixed_threshold_value	= 200;										//	fixed maximum threshold for nuclei segmentation (if auto doesn't work well);
var nuclei_min_circularity			= 0.0;										//	min circularity
var nuclei_min_area					= 50;										//	calibrated min nuclear size (in µm2)
var nuclei_max_area					= 600;										//	calibrated max nuclear size (in µm2)
var nuclei_overlap 					= 0.3;										//	nuclei_overlap amount tolerated for Stardist nuclei detection
var nuclei_probability				= 0.15;										//	minimal nuclei_probability for Stardist nuclei detection
var pixel_size						= 0.0642373;								//	pixel size (µm)

//	Boolean Parameters
var exclude_nuclei					= true;										//	analyze cellular regions without nuclear area
var flatfield						= false;									//	perform flatfield correction	
var mito_background					= true;
var mito_despeckle					= true;
var mito_clahe						= false;									
var nuclei_background				= false;										//	subtract nuclei_background for nuclei segmentation
var nuclei_clahe					= false;									// 	local contrast enhancement for nuclei segmentation
var nuclei_watershed				= false;									//	use nuclei_watershed for nuclei segmentation
var segment_cells					= true;										//	analyze cell and cytoplasmic ROIs
var segment_mito					= true;										//	analyze mitochondria
var segment_nuclei					= true;										//	analyze nuclear ROIs (currently redundant)

//	Arrays
var cells_segmentation_methods		= newArray("Threshold","Dilation","Trained Model","Cellpose");
var cols 							= newArray("01","02","03","04","05","06","07","08","09","10","11","12");
var dimensions						= newArray("xyczt(default)","xyctz","xytcz","xytzc","xyztc","xyzct");		
var file_list						= newArray(0);								
var file_types 						= newArray(".tif",".tiff",".nd2",".ids",".jpg",".mvd2",".czi");		
var nuclei_segmentation_methods		= newArray("Gauss","Laplace","Stardist");
var objects 						= newArray("Nuclei","Cells","Mitochondria");
var prefixes 						= newArray(0);
var rows 							= newArray("A","B","C","D","E","F","G","H");
var threshold_methods				= getList("threshold.methods");	
var threshold_methods				= Array.concat(threshold_methods,"Fixed");	

/*
 	***********************

		Macros

	***********************
*/

macro "Autorun"
{
	erase(1);
}

macro "Setup Action Tool - C888 T5f16S"
{
	setup();
}

macro "Analyse Single Image Action Tool - C888 T5f161"
{
	erase(0);
	setBatchMode(true);
	dir = getInfo("image.directory");
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir))File.makeDirectory(output_dir);
	start = getTime();
	title = getTitle; 
	prefix = substring(title,0,lastIndexOf(title,suffix));
	setFileNames(prefix);
	id = getImageID;
	if(flatfield)id = flatfield_correct(id);
	masknAllid = segmentNuclei(id,nuclei_channel,1); 
	
	// keep only the non-excluded nuclei (pos nuclei)
	newImage("posnuclei", "8-bit Black", image_width, image_height, 1); 
	masknid = getImageID;
	selectImage(masknid); 
	setForegroundColor(255,255,255);
	roiManager("Deselect");
	roiManager("Fill");
	run("Convert to Mask");

	nuclei_nr = roiManager("count");
	if(nuclei_nr>0)roiManager("Save",nuclei_roi_set);
	if(nuclei_nr>0 && segment_cells)
	{
		maskcid = segmentRegions(id, masknAllid, masknid, cells_channel, iterations);
		selectImage(maskcid);
		getStatistics(area, mean, min, max, std, histogram);
		if(max>0)roiManager("Save",cells_roi_set);
		else {File.delete(nuclei_roi_set); nuclei_nr=0;} 
		roiManager("reset");
		selectImage(masknAllid); close();
		
		individual=1;
		maskmid=segmentMito(id,mito_channel,masknid,maskcid,individual);
		mito_nr = roiManager("count");
		if(mito_nr>0)roiManager("Save",mito_roi_set);
	}
	
	if(isOpen(masknid)){selectImage(masknid); close;}
	if(isOpen(maskcid)){selectImage(maskcid); close;}
	if(isOpen(maskmid)){selectImage(maskmid); close;}
	
	readout = analyzeRegions(id);
	if(readout)summarizeResults(maskmid);
	print((getTime()-start)/1000,"sec");
	print("Analysis Done");
	setBatchMode("exit and display");
}

macro "Batch Analysis Action Tool - C888 T5f16#"
{
	erase(1);
	setBatchMode(true);
	setDirectory();
	prefixes = scanFiles();
	fields = prefixes.length;
	setup();
	start = getTime();
	for(i=0;i<fields;i++)
	{
		prefix = prefixes[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id = getImageID;
		if(flatfield)id = flatfield_correct(id);
		masknAllid = segmentNuclei(id,nuclei_channel,1); 
		
		// keep only the non-excluded nuclei (pos nuclei)
		newImage("posnuclei", "8-bit Black", image_width, image_height, 1); 
		masknid = getImageID;
		selectImage(masknid); 
		setForegroundColor(255,255,255);
		roiManager("Deselect");
		roiManager("Fill");
		if(is("Inverting LUT"))run("Invert LUT");
	
		nuclei_nr = roiManager("count");
		if(nuclei_nr>0)roiManager("Save",nuclei_roi_set);
		if(nuclei_nr>0 && segment_cells)
		{
			maskcid = segmentRegions(id, masknAllid, masknid, cells_channel, iterations);
			selectImage(maskcid);
			getStatistics(area, mean, min, max, std, histogram);
			if(max>0)roiManager("Save",cells_roi_set);
			else {File.delete(nuclei_roi_set); nuclei_nr=0;} 
			roiManager("reset");
			selectImage(masknAllid); close();
			
			individual=1;
			maskmid=segmentMito(id,mito_channel,masknid,maskcid,individual);
			mito_nr = roiManager("count");
			if(mito_nr>0)roiManager("Save",mito_roi_set);
		}
		
		if(isOpen(masknid)){selectImage(masknid); close;}
		if(isOpen(maskcid)){selectImage(maskcid); close;}
		if(isOpen(maskmid)){selectImage(maskmid); close;}
		
		readout = analyzeRegions(id);
		if(readout)summarizeResults(maskmid);
		selectImage(id); close;
		erase(0);
	}
	print((getTime()-start)/1000,"sec");
	if(isOpen("Log")){selectWindow("Log");saveAs("txt",log_path);}
	print("Complete Analysis Done");
	setBatchMode("exit and display");
}


macro "Segment Nuclei Action Tool - C999 H11f5f8cf3f181100 C999 P11f5f8cf3f18110 Ceee V5558"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	c 		= getNumber("Nuclear Channel",nuclei_channel);
	mid 	= segmentNuclei(id,c,0);
	setBatchMode("exit and display");
	roiManager("Show None");
	roiManager("Show All");
}

macro "Segment Cells Action Tool - C999 H11f5f8cf3f181100 Ceee P11f5f8cf3f18110"
{
	erase(0);
	setBatchMode(true);
	id = getImageID;
	c1 = getNumber("Nuclear Channel",nuclei_channel);
	c2 = getNumber("Cellular Mask Channel",cells_channel); // 0 if no additional mask
	masknAllid = segmentNuclei(id,nuclei_channel,1); 
	// keep only the non-excluded nuclei (pos nuclei)
	newImage("posnuclei", "16-bit Black", image_width, image_height, 1); 
	masknid = getImageID;
	selectImage(masknid); 
	roiManager("Deselect");
	roiManager("Fill");
	run("Convert to Mask");
	nuclei_nr = roiManager("count");
	if(nuclei_nr>0 && segment_cells)
	{
		maskcid = segmentRegions(id, masknAllid, masknid, cells_channel, iterations);
		selectImage(maskcid);
		getStatistics(area, mean, min, max, std, histogram);
		selectImage(masknAllid); close();
		selectImage(masknid); close();
		selectImage(maskcid); close();
	}
	setBatchMode("exit and display");
	roiManager("Show None");
	roiManager("Show All");
}

macro "Segment Mitochondria Action Tool - C999 H11f5f8cf3f181100 C999 P11f5f8cf3f18110 Ceee V3329 V6519 V8519 Va827"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	c 		= getNumber("Mito Channel",mito_channel); 
	
	selectImage(id);
	getDimensions(width, height, channels, slices, frames);
	
	newImage("Untitled", "8-bit white", width, height, 1); masknid=getImageID();
	newImage("Untitled", "8-bit white", width, height, 1); maskcid=getImageID();
	
	individual=0;
	maskmid = segmentMito(id,c,masknid,maskcid,individual);
	selectImage(masknid); close();
	selectImage(maskcid); close();
	selectImage(maskmid); close();
	setBatchMode("exit and display");
	roiManager("Show None");
	roiManager("Show All");
}

macro "Toggle Overlay Action Tool - Caaa O11ee"
{
	toggleOverlay();
}

macro "[t] Toggle Overlay"
{
	toggleOverlay();
}

macro "Verification Stack Action Tool - C888 T5f16V"
{
	erase(1);
	setBatchMode(true);
	setDirectory();
	prefixes = scanFiles();
	subset = getBoolean(""+prefixes.length+" images, select subset?");
	if(subset)
	{
		//	Compose an array containing all 96 well labels
		allWells = newArray(0);	
		for(i = 0; i < rows.length; i++)
  		{
  			for(j = 0; j < cols.length; j++)
  			{
  				allWells = Array.concat(allWells, rows[i] + cols[j]);
  			}
  		}
  		selection = "";
  		defaults = newArray(allWells.length);
  		Array.fill(defaults,0);		
  		selWells = newArray(allWells.length);		
		Dialog.create("Verification Stack");
		Dialog.addMessage("Select wells to visualize...")
		Dialog.addCheckboxGroup(rows.length,cols.length,allWells,defaults);
  		Dialog.show;
  		for(i = 0; i < rows.length; i++)
  		{
  			for(j = 0; j < cols.length; j++)
  			{
  				v = Dialog.getCheckbox();
  				if(v>0)selection = selection + allWells[i*cols.length+j];
  			}
  		}
		names = newArray(0);
  		for(i = 0;i<prefixes.length;i++)
		{
			name = prefixes[i];
			well = substring(name,indexOf(name,"Well")+4,indexOf(name,"Well")+7);
			if(indexOf(selection,well)>-1)names = Array.concat(names,name);
		}
	}
	else names = prefixes;
	createOverlay(names);
	setBatchMode("exit and display");
	run("Channels Tool... ");
}

/*
 	***********************

		Functions

	***********************
*/
function projectImages()
{
	erase(1);
	Dialog.create("Project Images...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Image",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	dir 		= getDirectory("");
	file_list 	= getFileList(dir);
	destination_dir 	= dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			print(i+1);
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			ns = nSlices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+ns/channels+" frames=1 display=Color");
			id = getImageID;
			title = getTitle;
			run("Z Project...", "projection=[Max Intensity]");
			zid = getImageID;		
			selectImage(zid); saveAs(suffix,destination_dir+pre+title+suffix);
			selectImage(id); close;
			selectImage(zid); close;
		}
	}
	print("Done");
}

function selectImages()
{
	erase(1);
	Dialog.create("Z-Select Images...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Image",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	dir 		= getDirectory("");
	file_list 		= getFileList(dir);
	destination_dir 	= dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	for(i = 0;i < file_list.length; i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			print(i+1);
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			id 		= getImageID;
			title 	= getTitle;
			getDimensions(width, height, channels, slices, frames);
			run("Duplicate...","title=Div duplicate");
			did 	= getImageID;
			run("Find Edges", "stack");
			ss = newArray(channels);	//	sharpest slices
			for(c = 1; c <= channels; c++)
			{
				stdmax = 0;
				selectImage(did);
				Stack.setChannel(c);
				for(z = 1; z <= slices; z++)
				{
					Stack.setSlice(z);
					getRawStatistics(nPixels, mean, min, max, std);
					if(std>stdmax)
					{
						stdmax = std;
						ss[c-1] = z;
					}
				}
			}
			selectImage(did); close;		
			for(c = 1; c <= channels; c++)
			{
				selectImage(id);
				Stack.setChannel(c);
				Stack.setSlice(ss[c-1]);
				if(c==1)
				{
					run("Duplicate...","title=Sel duplicate channels="+c+" slices="+ss[c-1]);
					zid = getImageID;
				}
				else
				{
					run("Select All");
					run("Copy");
					selectImage(zid);
					run("Add Slice");
					run("Paste");
					run("Select None");
				}
			}					
			selectImage(zid); saveAs(suffix,destination_dir+pre+title+suffix);
			selectImage(id); close;
			selectImage(zid); close;
		}
	}
	print("Done");
}

function setOptions()
{
	run("Options...", "iterations=1 count=1");
	run("Colors...", "foreground=white nuclei_background=black selection=yellow");
	run("Overlay Options...", "stroke=red width=1 fill=none");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
}

function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}

function erase(all)
{
	if(all){
		print("\\Clear");
		run("Close All");
	}
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
}

function setDirectory()
{
	dir = getDirectory("Choose a Source Directory");
	file_list = getFileList(dir);
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir))File.makeDirectory(output_dir);
	log_path = output_dir+"Log.txt";
}

function setFileNames(prefix)
{
	nuclei_roi_set 		= output_dir+prefix+"_nuclei_roi_set.zip";
	nuclei_results 		= output_dir+prefix+"_nuclei_results.txt";
	cells_roi_set 		= output_dir+prefix+"_cells_roi_set.zip";
	cells_results		= output_dir+prefix+"_cells_results.txt";
	mito_roi_set 		= output_dir+prefix+"_mito_roi_set.zip";
	mito_results		= output_dir+prefix+"_mito_results.txt";
	results				= output_dir+prefix+"_summary.txt";
}

function scanFiles()
{
	prefixes = newArray(0);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,suffix) && indexOf(path,"flatfield")<0)
		{
			print(path);
			prefixes = Array.concat(prefixes,substring(file_list[i],0,lastIndexOf(file_list[i],suffix)));			
		}
	}
	return prefixes;
}

function setup()
{
	setOptions();
	Dialog.createNonBlocking("MitoMetrics_v2 Settings");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("------------------------------------   General parameters  ------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Define here which regions you wish to analyze and whether you want to perform a distance calculation between ROI sets\n(not all comparisons are currently available)\n", 12, "#999999");
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Image Type", file_types, suffix);
	Dialog.addNumber("Pixel Size", pixel_size, 3, 5, micron+"");
	Dialog.addToSameRow();
	Dialog.addNumber("Field Number",fields, 0, 5, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Channel Number", channels, 0, 5, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Detectable Objects:");
	labels = newArray(3);
	defaults = newArray(3);
	labels[0] = "Nuclei";		
	defaults[0] = segment_nuclei;
	labels[1] = "Cells";		
	defaults[1] = segment_cells;
	labels[2] = "Mitochondria";		
	defaults[2] = segment_mito;
	Dialog.setInsets(0,150,0);
	Dialog.addCheckboxGroup(1,3,labels,defaults);
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("---------------------------------------------     Nuclei parameters  --------- ------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Nuclei are segmented by classic smoothing or Laplacian enhancement, or using a trained classifier (Stardist).\nIf the latter is used, standard settings are not considered, but object filtering is always applied.\n", 12, "#999999");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Nuclei Channel", nuclei_channel, 0, 4, "");
	Dialog.addToSameRow();
	Dialog.addChoice("Segmentation Method", nuclei_segmentation_methods, nuclei_segmentation_method);
	Dialog.addToSameRow();
	Dialog.addNumber("Blur Radius", nuclei_filter_scale, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("--------------------------------------------------------------------------------------------------------------   ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Standard Settings:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Background Subtraction", nuclei_background);
	Dialog.addToSameRow();
	Dialog.addCheckbox("Contrast Enhancement", nuclei_clahe);
	Dialog.setInsets(0,0,0);	
	Dialog.addCheckbox("Watershed Separation", nuclei_watershed);
	Dialog.addToSameRow();
	Dialog.addChoice("Threshold Method", threshold_methods, nuclei_threshold);
	Dialog.addToSameRow();
	Dialog.addNumber("Fixed Threshold", nuclei_fixed_threshold_value, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("--------------------------------------------------------------------------------------------------------------   ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Stardist settings:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Probability", nuclei_probability, 2, 4, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Tolerated overlap", nuclei_overlap, 2, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("--------------------------------------------------------------------------------------------------------------   ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Object filters:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Min. Circularity", nuclei_min_circularity, 2, 5, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Min. Area", nuclei_min_area, 0, 5, micron+"2");
	Dialog.addToSameRow();
	Dialog.addNumber("Max. Area", nuclei_max_area, 0, 5, micron+"2");	
	Dialog.setInsets(20,0,0);
	Dialog.addMessage("---------------------------------------------     Cell parameters  ---------------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Cells are segmented by thresholding a selected channel of a cell stain or if absent, by dilating nuclear ROI seeds.\nRegion growing is performed for a defined number of iterations, or if iterations is set to 0, sheer Voronoi tesselation is applied. \n", 12, "#999999");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Cell Channel", cells_channel,0,4," ");
	Dialog.addToSameRow();
	Dialog.addChoice("Segmentation Method", cells_segmentation_methods, cells_segmentation_method);
	Dialog.addToSameRow();
	Dialog.addCheckbox("Exclude Nuclei",exclude_nuclei);
	Dialog.addToSameRow();
	Dialog.addNumber("Filter scale ",cells_filter_scale,0,3,"");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("--------------------------------------------------------------------------------------------------------------   ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Threshold settings:                                                                                                      Dilation settings:                Cellpose settings:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Threshold Method", threshold_methods, cells_threshold);
	Dialog.addToSameRow();
	Dialog.addNumber("Fixed Threshold", cells_fixed_threshold_value, 0, 4, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Grow cycles", iterations, 0, 5, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Expected diameter", cells_diameter, 0, 5, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("---------------------------------------------     Mirochondria parameters  --------- ------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Mitochondria channel", mito_channel, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Backgroud Subtraction",mito_background);
	Dialog.addToSameRow();
	Dialog.addNumber("Rolling ball", mito_rolling, 2, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Contrast enhancement",mito_clahe);
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Scale laplace", mito_scale, 2, 4, "");
	Dialog.addToSameRow();
	Dialog.addChoice("Threshold Method", threshold_methods, mito_threshold_method);
	Dialog.addToSameRow();
	Dialog.addNumber("Fixed Threshold", mito_fixed_threshold_value, 0, 4, "");
	Dialog.addToSameRow();
	Dialog.addCheckbox("Despeckle",mito_despeckle);
	Dialog.addToSameRow();
	Dialog.addNumber("Min size", mito_minsize, 2, 4, "µm");
	Dialog.show();

	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	
	suffix							= Dialog.getChoice();		print("Image Type:",suffix);
	pixel_size						= Dialog.getNumber(); 		print("Pixel Size:",pixel_size);
	fields 							= Dialog.getNumber();	 	print("Fields:",fields);
	channels 						= Dialog.getNumber();		print("Channels:",channels);
	segment_nuclei					= Dialog.getCheckbox(); 	print("Segment Nuclei:",segment_nuclei);
	segment_cells 					= Dialog.getCheckbox();		print("Segment Cells:",segment_cells);
	segment_mito					= Dialog.getCheckbox();		print("Segment Mitochondria:",segment_mito);
	nuclei_channel 					= Dialog.getNumber();		print("Nuclear Channel:",nuclei_channel);
	nuclei_segmentation_method		= Dialog.getChoice();		print("Nuclei Segmentation Method:",nuclei_segmentation_method);
	nuclei_filter_scale				= Dialog.getNumber();		print("Nuclei Filter Scale:",nuclei_filter_scale);
	nuclei_background				= Dialog.getCheckbox();		print("Background Subtraction:",nuclei_background);
	nuclei_clahe					= Dialog.getCheckbox();		print("Clahe:",nuclei_clahe);
	nuclei_watershed 				= Dialog.getCheckbox();		print("Watershed:",nuclei_watershed);
	nuclei_threshold				= Dialog.getChoice();		print("Nuclear Autothreshold:",nuclei_threshold);
	nuclei_fixed_threshold_value	= Dialog.getNumber();		print("Fixed Threshold Value:",nuclei_fixed_threshold_value);
	nuclei_probability 				= Dialog.getNumber();		print("Probability:",nuclei_probability);
	nuclei_overlap 					= Dialog.getNumber();		print("Overlap:",nuclei_overlap);
	nuclei_min_circularity			= Dialog.getNumber();		print("Min Circ:",nuclei_min_circularity);
	nuclei_min_area					= Dialog.getNumber();		print("Min Nuclear Size:",nuclei_min_area);
	nuclei_max_area					= Dialog.getNumber();		print("Max Nuclear Size:",nuclei_max_area);
	cells_channel 					= Dialog.getNumber();		print("Cell Channel:",cells_channel);
	cells_segmentation_method		= Dialog.getChoice();		print("Cell Segmentation Method:",cells_segmentation_method);
	if(cells_segmentation_method=="Trained Model")
	{
		model_dir = getDirectory("Where is the trained cell model?");
		list  = getFileList(model_dir);
		for(i=0;i<list.length;i++)
		{
			path = model_dir+list[i];
			if (endsWith(list[i], ".model"))
			{
				model = path;
				print("location of the cell classification model:",model);
			}			
		}
	}
	exclude_nuclei					= Dialog.getCheckbox();		print("Exclude Nuclear Area From Cell Analysis", exclude_nuclei);
	cells_filter_scale				= Dialog.getNumber();		print("Cell Filter Scale:", cells_filter_scale);
	cells_threshold					= Dialog.getChoice();		print("Cell Autothreshold:",cells_threshold);
	cells_fixed_threshold_value		= Dialog.getNumber();		print("Fixed Threshold Value:",cells_fixed_threshold_value);
	iterations						= Dialog.getNumber();		print("Region Growing Iterations:",iterations);
	cells_diameter					= Dialog.getNumber();		print("Cellpose diameter:",cells_diameter);
	mito_channel					= Dialog.getNumber();		print("Mitochondria channel:", mito_channel);
	mito_background					= Dialog.getCheckbox();   	print("Mitochondria background subtraction:",mito_background);
	mito_rolling					= Dialog.getNumber();		print("Mitochondria rolling ball radius background subtraction:",mito_rolling);
	mito_clahe						= Dialog.getCheckbox();		print("Mitochondria contrast enhancement:", mito_clahe);
	mito_scale						= Dialog.getNumber();		print("Mitochondria Laplace scale:",mito_scale);
	mito_threshold_method			= Dialog.getChoice();		print("Mitochondria Autothreshold:",mito_threshold_method);
	mito_fixed_threshold_value		= Dialog.getNumber();		print("Mitochondria Fixed threshold:", mito_fixed_threshold_value);
	mito_despeckle					= Dialog.getCheckbox();		print("Mitochondria Despeckle:",mito_despeckle);
	mito_minsize					= Dialog.getNumber();		print("Mitochondria min size:",mito_minsize);
	
}

function calibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!=micron)run("Properties...", " unit="+micron+" pixel_width="+pixel_size+" pixel_height="+pixel_size);
	else pixel_size = pixelWidth;
}

function decalibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!="pixel")run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
}

function flatfield_correct(id){
	//	flatField does a flatfield field correction, using an image of fluorescent plastic
	//	if flatfield contains different channels and input image only one, c directs to the correct channel to select in the flatfield stack
	selectImage(id);
	title = getTitle;
	flatpath = dir+"flatfield.tif";
	open(flatpath);
	fid = getImageID;
	selectImage(fid);
	getRawStatistics(n,flatmean);
	imageCalculator("Divide create 32-bit stack", id,fid); 
	selectImage(fid); close;
	selectImage(id); close;
	selectWindow("Result of "+title);
	id = getImageID;
	selectImage(id);
	rename(title);  
	run("Multiply...","value="+flatmean+" stack");
	return id;
}

function segmentNuclei(id,c,sel)
{
	// input = multichannel image, output = roiset of nuclear ROIs and if(sel==1) mask incl. border objects
	// output = an image (mid) that contains all ROIs (also touching borders) and roiset of full nuclei
	mid = 0;
	selectImage(id);
	image_width = getWidth;
	image_height = getHeight;
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID; // the nuclear channel image to be turned into a binary image
	calibrateImage(cid);
	// preprocess the image
	selectImage(cid);
	if(nuclei_clahe)run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	if(nuclei_background)run("Subtract Background...", "rolling="+round(30/pixel_size));
	if(nuclei_segmentation_method != "Laplace" && nuclei_filter_scale > 0)run("Gaussian Blur...", "sigma="+nuclei_filter_scale);
	else if(nuclei_segmentation_method == "Laplace")
	{
		run("FeatureJ Laplacian", "compute smoothing="+nuclei_filter_scale); // scale to be adapted depending on nuclei size and SNR
		selectImage(cid); close;
		selectImage("copy Laplacian");
		rename("copy");
		cid = getImageID;
		selectImage(cid);
	}
	if(nuclei_segmentation_method != "Stardist")
	{
		if(nuclei_threshold=="Fixed")
		{
			if(nuclei_segmentation_method == "Laplace")setAutoThreshold("Default ");
			else setAutoThreshold("Default dark");
			getThreshold(mit,mat); 
			setThreshold(nuclei_fixed_threshold_value,mat);
		}
		else {
			if(nuclei_segmentation_method == "Laplace")setAutoThreshold(nuclei_threshold); 
			else setAutoThreshold(nuclei_threshold+" dark"); 
		}
		getThreshold(mit,mat); print("Nuclei Threshold:",mit,mat);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Fill Holes");
		if(nuclei_watershed)run("Watershed");
	}
	else if(nuclei_segmentation_method == "Stardist")
	{
		selectImage(cid);
		run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
		+"args=['input':"+'copy'+", 'modelChoice':'Versatile (fluorescent nuclei)',"
		+"'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99',"
		+"'probThresh':'"+nuclei_probability+"', 'nmsThresh':'"+nuclei_overlap+"', 'outputType':'ROI Manager', 'nTiles':'1', "
		+"'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', "
		+"'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
		selectImage(cid); close;
		newImage("copy", "8-bit black", image_width, image_height, 1);
		calibrateImage(cid);
		cid = getImageID;
		selectImage(cid);
		selectImage("copy");
		nr_rois = roiManager("count");
		for(r = 0; r < nr_rois; r++)
		{
			roiManager("select",r);
			run("Enlarge...", "enlarge=1");
			run("Clear");
			run("Enlarge...", "enlarge=-1");
			run("Fill");
		}
		roiManager("Deselect");
		roiManager("reset");
		setThreshold(1,255);
		run("Convert to Mask");
	}
	run("Set Measurements...", "area centroid mean integrated redirect=None decimal=2");
	if(sel)
	{
		selectImage(cid);
		run("Analyze Particles...", "size="+nuclei_min_area+"-"+nuclei_max_area+" circularity="+nuclei_min_circularity+"-1.00 show=Masks clear include");
		if(isOpen("Mask of copy"))
		{
			selectWindow("Mask of copy"); 
			mid = getImageID;   // the full mask of all particles (for more accurate cell segmentation)
		}
	}
	selectImage(cid);
	run("Analyze Particles...", "size="+nuclei_min_area+"-"+nuclei_max_area+" circularity="+nuclei_min_circularity+"-1.00 show=Nothing exclude clear include add");
	rmc = roiManager("count"); print(rmc,"Nuc. ROI");
	if(rmc==0 && isOpen(mid)){selectImage(mid); close; mid=0;}
	for(i=0;i<rmc;i++)
	{
		roiManager("select",i);
		if(i<9)roiManager("Rename","000"+i+1);
		else if(i<99)roiManager("Rename","00"+i+1);	
		else if(i<999)roiManager("Rename","0"+i+1);	
		else roiManager("Rename",i+1);
	}
	run("Clear Results");
	roiManager("deselect"); 
	roiManager("Measure");
	selectImage(cid); close; 
	return mid;
}

function segmentRegions(id, mid,masknid, c, iterations){
	selectImage(mid); 
	run("Select None");
	name = getTitle;
	
	// generate voronoi regions from all detected nuclei (including edges) to have some rough boundaries between touching cells
	run("Duplicate...","title=voronoi");
	vid = getImageID;
	selectImage(vid);
	run("Voronoi");
	setThreshold(1, 255);
	run("Convert to Mask");
	run("Invert");
	
	// generate dilated nuclei (using more accurate EDM mapping (by x iterations) requires the biovioxxel package) - now just using the enlarge function
	if(cells_segmentation_method=="Dilation" && iterations != 0)
	{
		selectImage(mid);
		run("Duplicate...","title=dilate");
		did = getImageID;
		selectImage(did); 
		
		//run("EDM Binary Operations", "iterations="+iterations+" operation=dilate"); 
		run("Create Selection"); 
		run("Enlarge...", "enlarge="+iterations+" pixel");
		roiManager("Add");
		run("Select All");
		run("Fill");
		sel = roiManager("count")-1;
		roiManager("select", sel);
		run("Clear", "slice");
		roiManager("Delete");
		run("Select None");
		run("Invert LUT");
		imageCalculator("AND create", "voronoi","dilate");
		selectImage(did); close; 
		selectImage(vid); close;
		selectWindow("Result of voronoi");
		if(!is("Inverting LUT"))run("Invert LUT");
		//run("Invert");
		vid = getImageID;
		selectImage(vid);
		rename("Cell_ROI");
	}

	// use a cellular counterstain to make accurate cell ROIs
	if(cells_segmentation_method!="Dilation"){
		selectImage(id);
		run("Select None");
		if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
		else{setSlice(c);run("Duplicate...","title=copy ");}
		cid = getImageID;
		selectImage(cid);
		if(cells_filter_scale > 0)run("Gaussian Blur...", "sigma="+cells_filter_scale);
		if(cells_segmentation_method=="Trained Model")
		{
			print("Applying trained model, takes time!"); //model needs to be in the FIJI/scripts foldr along with the Trained_Glia_Segmentation.bsh script in teh FIJI/scripts folder
			run("Trained Cell Segmentation","model="+model+" image="+cid);
			selectImage(cid); close;
			selectWindow("Classification result");
			rename("copy");
			cid = getImageID();		
			setAutoThreshold("Default "); 	
		}
		else if(cells_segmentation_method=="Threshold")
		{
			if(cells_threshold=="Fixed"){
				setAutoThreshold("Default dark");
				getThreshold(mit,mat); 
				setThreshold(cells_fixed_threshold_value,mat);
			}
			else {
				setAutoThreshold(cells_threshold+" dark"); 
			}
		}
		else if(cells_segmentation_method=="Cellpose"){
			print("Cellpose Start");

			//Create RGB image for cellpose: R=1, G=2, B=3
			selectImage(id);
			run("Select None");
			if(Stack.isHyperstack)run("Duplicate...", "title=Nuclei_Red duplicate channels="+nuclei_channel);	
			else{ setSlice(nuclei_channel); run("Duplicate...","title=Nuclei_Red ");}
			redid=getImageID();
			selectImage(id);
			run("Select None");
			if(Stack.isHyperstack)run("Duplicate...", "title=Cells_Green duplicate channels="+cells_channel);	
			else{setSlice(cells_channel);run("Duplicate...","title=Cells_Green ");}
			greenid=getImageID();
			run("Merge Channels...", "c1=[Nuclei_Red] c2=[Cells_Green] create");
			compid=getImageID();
			run("RGB Color");
			rgbid=getImageID();
			selectImage(cid); close;

			//Run Cellpose and exclude roi borders from label map
			run("Cellpose Script Wrapper", "imp="+id+" diameter="+cells_diameter+" channel_nuc="+1+" channel_cyto="+2);
			
			selectWindow("Cellpose_Label");
			cid = getImageID();
			rename("copy");
			selectImage(rgbid); close;
			selectImage(compid); close;
			setForegroundColor(0, 0, 0);
			getRawStatistics(nPixels, mean, min, max, std, histogram);	
			maxMask=max;
			minPixels=(nuclei_min_area/(pixel_size*pixel_size))/10;  //min pixel number equal to 1/10 of min nuc area	
			maxPixels=nPixels;		
			for(i=1; i<=maxMask; i++)
			{
				selectImage(cid);
				setThreshold(i, i, "raw");
				run("Create Selection");
				getRawStatistics(nPixels, mean, min, max, std, histogram);	
				
				if(nPixels > minPixels && nPixels!=maxPixels){
					run("Draw", "slice");
				}
			}
			setThreshold(1, maxMask, "raw");
			setForegroundColor(255, 255, 255);
			print("Cellpose End");
		}
		selectImage(cid);
		getThreshold(mit,mat); 
		print("Cell Threshold:",mit,mat);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		//run("Fill Holes");	
		
		if(cells_segmentation_method!="Cellpose")
		{
			imageCalculator("AND create", "voronoi","copy"); // apply the voronoi bounderies to teh cell ROIs
			selectImage(cid); close; 
			selectImage(vid); close;
			selectWindow("Result of voronoi");
			vid = getImageID;
		}else { //cellpose
			selectImage(vid); close;
			selectImage(cid);
			vid=getImageID();
		}
		selectImage(vid);
		rename("Cell_ROI");
	}

	// fuse cell and nuclei image to avoid non-overlapping ROI (not with cellpose)
	if(cells_segmentation_method!="Cellpose"){
		imageCalculator("OR create", "Cell_ROI","posnuclei"); 
		selectImage(vid); close;
		selectWindow("Result of Cell_ROI");
		rename("Cells");
		vid = getImageID;
	}else{
		selectImage(vid);
		rename("Cells");
	}
	
	// make labeled nuclei centroid mask
	newImage("centroid", "16-bit Black", image_width, image_height, 1); 
	maskcentroidid=getImageID();
	calibrateImage(maskcentroidid);
	run("16-bit"); //Needed if > 255 nuclei detected
	rmc = roiManager("count"); // rmc = number of nuclei
	print(rmc,"retained nuclei");
	selectImage(maskcentroidid);
	run("Set Measurements...", "centroid redirect=None decimal=4");
	roiManager("measure");
	for(i=0;i<rmc;i++)
	{
		Mx=getResult("X", i)/pixel_size;
		My=getResult("Y", i)/pixel_size;
	 	makeOval(Mx-2, My-2, 4, 4);
		run("Set...", "value="+i+1);		
	}
	run("Select None");
	
	// Add cell rois and rename with matching nuclear index
	selectImage(vid);
	run("Analyze Particles...", "size="+nuclei_min_area+"-Infinity circularity=0.00-1.00 show=Nothing add");
	rmcb = roiManager("count"); // number of cell regions larger than a nucleus 
	print(rmcb-rmc,"Number of detected cell regions larger than a min. nuclear area"); 
	selectImage(maskcentroidid);
	for(i=rmcb-1;i>=rmc;i--)
	{
		roiManager("select",i);
		getRawStatistics(np,mean,min,max);
		if(max==0){roiManager("delete");} 				// no nuc so not retained
		else if(max<10)roiManager("Rename","000"+max);	// assigned to correct nuc
		else if(max<100)roiManager("Rename","00"+max);	
		else if(max<1000)roiManager("Rename","0"+max);	
		else roiManager("Rename",max);
	}	
	selectImage(maskcentroidid); close();
	
	rmcc = roiManager("count"); //all cell regions with one nucleus	
	print(rmcc-rmc,"Number of unique cell regions that overlap with a nucleus"); 
	
	// exclude nuclei
	if(rmcc>rmc)
	{
		if(exclude_nuclei) //define cytoplasmic regions (cells without nuclei)
		{
			for(i=0;i<rmc;i++)
			{
				roiManager("select",i); 
				roi_name_a = Roi.getName();
				getRawStatistics(np_a); //area
				for(j=rmc;j<rmcc;j++)
				{
					roiManager("select", j); 
					roi_name_b = Roi.getName(); 
					getRawStatistics(np_b);
					if (matches(roi_name_a, roi_name_b) && np_a!=np_b) //cell roi is not equally large as nuc roi
					{ 
						//print("Matched with Cell ROI",roi_name_b);
						couple = newArray(i,j);
						roiManager("select",couple); 
						roiManager("XOR"); // only retain cytoplasm
						roiManager("Add"); // add cytoplasm roi
						roiManager("select",roiManager("count")-1);
						roiManager("Rename",roi_name_a);
					}
					else if(matches(roi_name_a, roi_name_b) && np_a==np_b)
					{
						print("Matched with Cell ROI",roi_name_b,"But discarded due to equal size (no cytoplasm detected)");
					}
				}
			}
			roiManager("select",Array.getSequence(rmcc)); 
			roiManager("Delete"); 		//remove all prior ROIs	(nuclei and cells)
		} 
		else if(!exclude_nuclei)
		{
			roiManager("select",Array.getSequence(rmc)); 
			roiManager("Delete"); //remove all nuclei ROIs and retain cell ROIs
		} 
		run("Select None");
		cell_nr = roiManager("count");		
	}
	else 
	{
		cell_nr = 0;
	}
	print(cell_nr, "Number of unique cell regions retained with area larger than a nucleus");
	//selectImage(mid); close;
	selectImage(vid); close;
	//selectImage(pid); close;
	if(isOpen("copy")){selectWindow("copy");close;}
	
	newImage("maskCell", "8-bit Black", image_width, image_height, 1); 
	
	maskcid = getImageID;
	selectImage(maskcid); 
	setForegroundColor(255, 255, 255);
	setBackgroundColor(0, 0, 0);
	roiManager("Deselect");
	roiManager("Fill");
	if(!is("Inverting LUT"))run("Invert LUT");
	
	return maskcid;
}

function segmentMito(id,c,masknid,maskcid,individualMito)
{
	selectImage(id);
	run("Select None");
	if(Stack.isHyperstack)run("Duplicate...", "title=copyChannel duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copyChannel ");}
	cid = getImageID;
	selectImage(cid);
	title=getTitle();
	calibrateImage(cid);
	
	//Pre-processing
	if(mito_background)run("Subtract Background...", "rolling="+mito_rolling);
	if(mito_clahe) run("Enhance Contrast...", "saturated=0.1 normalize process_all");
	
	//multi-scale laplace
	e = 0;
	while(e<mito_scale)
	{			
		e++;
		selectImage(cid);
		run("FeatureJ Laplacian", "compute smoothing="+e);
		selectWindow(title+" Laplacian");
		run("Multiply...","value="+e*e);
		rename("scale "+e);
		eid = getImageID;
		if(e>1)
		{
			selectImage(eid);run("Select All");run("Copy");close;
			selectImage(fid);run("Add Slice");run("Paste");
		}
		else fid = getImageID;
	}
	selectImage(fid);
	run("Z Project...", "start=1 projection=[Sum Slices]");
	mAllid = getImageID;
	selectImage(fid); close;
	
	selectImage(mAllid);
	if(mito_threshold_method=="Fixed"){
		setAutoThreshold("Default ");
		getThreshold(mit,mat); 
		setThreshold(minOf(mit,-mito_fixed_threshold_value),-mito_fixed_threshold_value);
	}else{
		setAutoThreshold(""+mito_threshold_method);
	}
	run("Convert to Mask");
	rename("maskMitoAll");
	selectImage(mAllid);


	//Combine with nuclear and cell mask
	selectImage(masknid);
	rename("maskNuc");
	selectImage(maskcid);
	rename("maskCell");
	imageCalculator("OR create", "maskCell","maskNuc"); 
	rename("maskRegion");
	imageCalculator("AND create", "maskRegion","maskMitoAll"); 
	maskmid=getImageID();
	rename("maskMito");
	selectImage(maskmid);
	if(mito_despeckle){run("Despeckle");run("Despeckle");}
	run("Convert to Mask");
	
	run("Set Measurements...", "area mean standard perimeter shape feret's integrated redirect=["+title+"] decimal=4");
	run("Analyze Particles...", "size="+mito_minsize+"-Infinity circularity=0.00-1.00 show=Nothing add");
	nr=roiManager("count");
	
	//	A second pass on the inverted image captures the holes (of touching mitochondria) that were missed by 
	//	the first pass and reassigns them to the correct ROI
	if(individualMito && nr>0)
	{
		selectImage(maskmid);run("Invert");
		for(r=0;r<nr;r++)
		{
			roiManager("select",r);
			run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing add");
			nrn=roiManager("count");
			if(nrn>nr)
			{
				indices=newArray(0);
				for(ind=nr;ind<nrn;ind++){indices=Array.concat(indices,ind);}
				tindices=Array.concat(r,indices);
				roiManager("select",tindices);
				roiManager("XOR");roiManager("Update");	
				roiManager("Deselect");roiManager("select",indices);roiManager("Delete");
			}
		}
	}else {
		selectImage(maskmid);
		roiManager("Deselect");
		run("Select None");
		run("Create Selection");
		roiManager("Add");
	}
	
	roiManager("Deselect");
	run("Select None");
	selectImage(mAllid); close();
	selectImage(cid); close();
	selectWindow("maskRegion"); close();
	
	return maskmid;
}

function analyzeRegions(id)
{
	erase(0); 
	mask = 0;
	readout = 1;
	//	analyze cell rois
	selectImage(id);
	calibrateImage(id);
	if(File.exists(cells_roi_set))
	{
		run("Set Measurements...", "area mean standard modal min centroid center perimeter shape integrated median skewness kurtosis redirect=None decimal=4");
		roiManager("Open",cells_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		
		sortResults(); // organize results per channel
		
		if(!isOpen(mask))
		{
			newImage("Mask", "32-bit Black",image_width, image_height, 1); 	//	reference image for spot assignments
			mask = getImageID; 
		}
		selectImage(mask);
		for(j=0;j<rmc;j++)
		{
			roiManager("select",j);
			index = getInfo("roi.name");
			run("Set...", "value="+0-index);							//	negative values for cytoplasm, positive for nuclei
			setResult("Nuclei",j,index);
		}	
		saveAs("Measurements",cells_results);
		erase(0);
	}	
	//	analyze nuclear rois
	if(File.exists(nuclei_roi_set))
	{
		run("Set Measurements...", "area mean standard modal min centroid center perimeter shape feret's integrated median skewness kurtosis redirect=None decimal=4");
		roiManager("Open",nuclei_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
	
		sortResults(); //organise results per channel
				
		if(!isOpen(mask))
		{
			newImage("Mask", "32-bit Black",image_width, image_height, 1); //	reference image for spot assignments 
			mask = getImageID;
		}
		selectImage(mask);
		for(j=0;j<rmc;j++)
		{
			roiManager("select",j);
			index = getInfo("roi.name");
			run("Set...", "value="+index);					//	negative values for cytoplasm, positive for nuclei
			setResult("Cell",j,index);
		}	
		run("Select None");
		updateResults;
		saveAs("Measurements",nuclei_results);
		erase(0);
	}	
	
	//	analyze mito rois
	if(File.exists(mito_roi_set))
	{
		run("Set Measurements...", "area mean standard modal min centroid center perimeter shape integrated median skewness kurtosis redirect=None decimal=4");
		roiManager("Open",mito_roi_set);
		mito_nr = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		// determine the location of the spots (cell vs. nucleus)
		selectImage(mask); setSlice(1);
		roiManager("deselect");
		roiManager("Measure");
		nindices = newArray(mito_nr);
		cindices = newArray(mito_nr);	
		
		for(j=0;j<mito_nr;j++)
		{
			min = getResult("Min",j);
			max = getResult("Max",j);
			if(max>0){
				nindices[j] = max; 
				if(exclude_nuclei){
					if(min > 0){cindices[j] = 0;}
					else{cindices[j] = -min;}
				}
				else{cindices[j] = max;}
			}else if(min<0){nindices[j] = 0; cindices[j] = -min;}
		}	
		run("Clear Results");
		IJ.renameResults("Temp","Results");
		for(j=0;j<mito_nr;j++)
		{
			setResult("Nucleus",j,nindices[j]);
			setResult("Cell",j,cindices[j]);
		}
		updateResults;
		saveAs("Measurements",mito_results);
		erase(0);
	}
		
		/*
		selectImage(maskmid);
		selectWindow("maskMito");
		run("2D Analysis", "count total mean total_0 mean_0 mean_1 mean_2 branches total_1 mean_3 branch branch_0 mean_4 if=Count perform_0 area perimeter form aspect branches_0 total_2 mean_5 branch_1 branch_2 mean_6 longest mask=None mask_0=Mask second=None second_0=[Channel 2] third=None third_0=[Channel 3] =None to=None then=None");
		run("2D Analysis", "count total mean total_0 mean_0 mean_1 mean_2 branches total_1 mean_3 branch branch_0 mean_4 if=Count perform_0 area perimeter form aspect branches_0 total_2 mean_5 branch_1 branch_2 mean_6 longest mask=None mask_0=Mask second=None second_0=[Channel 2] third=None third_0=[Channel 3] =None to=None then=None");
		setBatchMode(true);
		
		roiManager("Open",mito_roi_set);
		mito_nr = roiManager("count");
		
		// determine the location of the mitochondria (cell vs. nucleus)
		run("Set Measurements...", "min max redirect=None decimal=4");
		selectImage(mask);
		roiManager("deselect");
		roiManager("Measure");
		nindices = newArray(mito_nr);
		cindices = newArray(mito_nr);	
		
		for(j=0;j<mito_nr;j++)
		{
			min = getResult("Min",j);
			max = getResult("Max",j);
			if(max>0){
				nindices[j] = max; 
				if(exclude_nuclei){
					if(min > 0){cindices[j] = 0;}
					else{cindices[j] = -min;}
				}
				else{cindices[j] = max;}
			}else if(min<0){nindices[j] = 0; cindices[j] = -min;}
		}	
		run("Clear Results");
		
		IJ.renameResults("2D Analysis Data - per Mito","Results");
		for(j=0;j<mito_nr;j++)
		{
			setResult("Nucleus",j,nindices[j]);
			setResult("Cell",j,cindices[j]);
		}
		updateResults;
		saveAs("Measurements",mito_results);
		selectImage(maskmid); close();
		erase(0);
		*/

	if(isOpen(mask)){selectImage(mask); close;}
	else readout = 0;	
	
	return readout;
}

function summarizeResults(maskmid)
{
	// 	open nuclei results
	run("Results... ", "open=["+nuclei_results+"]");
	nnr 			= nResults;
	nindices		= newArray(nnr);
	resultLabels 	= getResultLabels();
	matrix 			= results2matrix(resultLabels);
	selectWindow("Results"); 
	run("Close");
	for(r=0;r<nnr;r++)
	{
		for(s=0;s<resultLabels.length;s++)
		{
			selectImage(matrix);
			p = getPixel(s,r);
			if(resultLabels[s]!="Cell" && resultLabels[s]!="X_MC1" && resultLabels[s]!="Y_MC1")setResult("Nucl_SC"+nuclei_channel+"_"+resultLabels[s],r,p); // Label all nuclear measured parameters except for the cell or X and Y indices with a "Nucl" prefix
			else if(resultLabels[s]=="X_MC1")setResult("X",r,p);  //exception for X,Y coordinates for ease of tracing-back
			else if(resultLabels[s]=="Y_MC1")setResult("Y",r,p); 
			else setResult(resultLabels[s],r,p);
		}
	}
	updateResults;
	selectImage(matrix); close;
	
	//	append cellular results
	if(File.exists(cells_results))
	{	
		// once in a while a cell index is different from a nuclear index
		for(r=0;r<nnr;r++){nindices[r]=getResult("Cell",r)-1;} 
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+cells_results+"]");
		cnr				= nResults;
		cindices		= newArray(cnr);
		for(r=0;r<cnr;r++){cindices[r]=getResult("Nuclei",r)-1;}
		 
		//Append results based on roi index in cell results 
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(r=0;r<cnr;r++) 
		{
			for(s=0;s<resultLabels.length;s++)
			{
				selectImage(matrix);
				p = getPixel(s,r);
				setResult("Cell_SC"+cells_channel+"_"+resultLabels[s],cindices[r],p); // Label all cytoplasmic measured parameters with a "Cell" prefix
			}
		}
		updateResults;
		selectImage(matrix); close;
	}
	
	if(File.exists(mito_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+mito_results+"]");
		mnr 			= nResults;
		nindices 		= newArray(mnr);
		cindices 		= newArray(mnr);
		for(j=0;j<mnr;j++)
		{
			nindices[j] = getResult("Nucleus",j)-1;
			cindices[j] = getResult("Cell",j)-1;
		}	
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(s=0;s<resultLabels.length;s++)
		{
			if(resultLabels[s] != "Nucleus" && resultLabels[s] != "Cell")
			{
				nvalues 	= newArray(nnr);
				cvalues 	= newArray(nnr);
				nnumber 	= newArray(nnr);
				cnumber 	= newArray(nnr);
				for(r=0;r<mnr;r++)
				{
					selectImage(matrix);
					p = getPixel(s,r);
					if(nindices[r]>=0)
					{
						nvalues[nindices[r]] += p;  
						nnumber[nindices[r]] += 1;	
					}
					if(cindices[r]>=0)
					{
						cvalues[cindices[r]] += p;  
						cnumber[cindices[r]] += 1;	
					}
				}
				
				for(r=0;r<nnr;r++)
				{
					setResult("Mito_SC"+mito_channel+"_NrPerNuc",r,nnumber[r]);
					setResult("Mito_SC"+mito_channel+"_"+resultLabels[s]+"_SumPerNuc",r,nvalues[r]);              
					setResult("Mito_SC"+mito_channel+"_"+resultLabels[s]+"_MeanPerNuc",r,nvalues[r]/nnumber[r]);
					if(segment_cells)
					{
						setResult("Mito_SC"+mito_channel+"_NrPerCell",r,cnumber[r]);
						setResult("Mito_SC"+mito_channel+"_"+resultLabels[s]+"_SumPerCell",r,cvalues[r]);
						setResult("Mito_SC"+mito_channel+"_"+resultLabels[s]+"_MeanPerCell",r,cvalues[r]/cnumber[r]);
					}
				}
			}
		}
		selectImage(matrix); close;
		updateResults();
	}
	selectWindow("Results"); saveAs("Measurements",results);
}



function sortResults()
{
	resultLabels = getResultLabels();
	matrix = results2matrix(resultLabels);
	matrix2results(matrix,resultLabels,channels);
}

function getResultLabels()
{
	selectWindow("Results");
	ls 				= split(getInfo(),'\n');
	rr 				= split(ls[0],'\t'); 
	nparams 		= rr.length-1;			
	resultLabels 	= newArray(nparams);
	for(j=1;j<=nparams;j++){resultLabels[j-1]=rr[j];}
	return resultLabels;
}

function results2matrix(resultLabels)
{
	h = nResults;
	w = resultLabels.length;
	newImage("Matrix", "32-bit Black",w, h, 1);
	matrix = getImageID;
	for(j=0;j<w;j++)
	{
		for(r=0;r<h;r++)
		{
			v = getResult(resultLabels[j],r);
			selectImage(matrix);
			setPixel(j,r,v);
		}
	}
	run("Clear Results");
	return matrix;
}

function matrix2results(matrix,resultLabels,channels)
{
	selectImage(matrix);
	w = getWidth;
	h = getHeight;
	for(c=0;c<channels;c++)
	{
		start = c*h/channels;
		end = c*h/channels+h/channels;
		for(k=0;k<w;k++)
		{
			for(j=start;j<end;j++)
			{
				selectImage(matrix);
				p = getPixel(k,j);
				setResult(resultLabels[k]+"_MC"+c+1,j-start,p); // MC for measurement channel
			}
		}
	}
	selectImage(matrix); close;
	updateResults;
}

function toggleOverlay()
{	
	run("Select None"); 
	roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}
	
function createOverlay(names)
{
	setForegroundColor(25, 25, 25);
	fields = names.length;
	print(fields,"images");
	for(i=0;i<fields;i++)
	{
		prefix = names[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id = getImageID;
		Stack.getDimensions(w,h,channels,slices,frames); 
		if(!Stack.isHyperStack && channels == 1)
		{
			channels = slices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
		}
		id = getImageID;

		if(segment_nuclei)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(nuclei_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",nuclei_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
	
		if(segment_cells)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(cells_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",cells_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}

		if(segment_mito)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(mito_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",mito_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
	}
	run("Concatenate...", "all_open title=[Concatenated Stacks]");
	Stack.getDimensions(w,h,newchannels,slices,frames);
	for(c=1;c<=channels;c++){Stack.setChannel(c);Stack.setFrame(round(frames/2));resetMinAndMax;}
	range = pow(2,bitDepth);
	for(c=channels+1;c<=newchannels;c++){Stack.setChannel(c);setMinAndMax(0,range/2);}
	run("Make Composite");
}