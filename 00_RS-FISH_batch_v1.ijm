close("*");
// This macro script runs the RS (radial symmetry) FISH FIJI plug-in on a subset of .tif images in all the sub-directories of the defined parent dir (dialog window)
// After identifying the best parameters using the RS-FISH plugin GUI (interactive mode) on example images, this macro can be executed on larger sets of images.
// Parameters can be defined in the beginning of the script and will be used on all images in the specified parent dir.

// To run this macro in the ImageJ headless mode via the command line on the curry cluster use:
// <FIJI/DIR/PATH>/ImageJ-linux64 --headless --run </PATH/TO/THIS/SCRIPT>/RS_macro.ijm &> </PATH/TO/WHERE/YOU/WANT/YOUR/LOGFILE>.log
// NOTE: For headless mode the dir needs to be defined in this script or passed as an argument with dir=getArgument();
// ImageJ-linux64 --headless -macro Macro_TxtFolderList_IMS_maxZ_splitCh_saveTIFF_headless_LF.ijm dir

// The detection result table (.csv file with spot coordinates and intensities) will be saved in the dir of each input image.

// Select the parent dir containing the input images (also in subdirs)
dir = getDirectory("Choose the parent directory of your experiment (containing the images to process directly or in subfolders)");
outdir = dir + "/RSFISH/";
File.makeDirectory(outdir);

// Ask for a specific string at the image file end to subset the file list the script should run on, e.g. channel name like "637nm.tif"
// If nothing is specified "_maxZ.tif" will be used as default
var imend = "_maxZ.tif"
imgfilter = getString("Enter an endsWith-string for subsetting image files in parent dir:",imend);

//////// Define RS-FISH parameters //////////

anisotropyCoefficient = 1.0; 
ransac = "RANSAC";						// options are "RANSAC" (log value "SIMPLE") / "No RANSAC" / "MULTICONSENSU"
imMin = 0; 								// img min intensity
imMax = 65536; 							// max intensity
sigmaDoG = 0.85; 
thresholdDoG = 0.002;
supportRadius = 2;
inlierRatio = 0.1;						// meaning: min inlier ratio
maxError = 1.5; 						// meaning: max error range
intensityThreshold = 0;  				// meaning: spot intensity threshold
bsMethod = "No background subtraction";	// Log file 0 / 1 / 2 / 3 / 4 options correspond to "No background subtraction" / "Mean" / "Median" / "RANSAC on Mean" / "RANSAC on Median"
bsMaxError = 0.05;						// background subtraction param
bsInlierRatio = 0.0;					// background subtraction param
useMultithread = "";					// Not from log file (only in advanced mode)! If you wish to use multithreading "use_multithreading", else "" (empty string)
numThreads = 40;						// multithread param
blockSizX = 128;                     	// multithread param
blockSizY = 128;						// multithread param
blockSizZ = 16;							// multithread param

// RS-FISH parameters variable
RS_parameters = " mode=Advanced anisotropy=" + anisotropyCoefficient + " robust_fitting=[" + ransac + "] use_anisotropy" + 
	" image_min=" + imMin + " image_max=" + imMax + " sigma=" + sigmaDoG + " threshold=" + thresholdDoG + 
	" support=" + supportRadius + " min_inlier_ratio=" + inlierRatio + " max_error=" + maxError + " spot_intensity_threshold=" + intensityThreshold + 
	" background=[" + bsMethod + "] background_subtraction_max_error=" + bsMaxError + " background_subtraction_min_inlier_ratio=" + bsInlierRatio + 
	" " + useMultithread + " num_threads=" + numThreads + " block_size_x=" + blockSizX + " block_size_y=" + blockSizY + " block_size_z=" + blockSizZ;

ransac_sub = split(ransac, ' ');
ransac_sub = ransac_sub[0];

bsMethod_sub = split(bsMethod, ' ');
bsMethod_sub = bsMethod_sub[0];

// set some options
setBatchMode(true);

///////////////////////////////////////////////////

////////////// DEFINE FUNCTIONS //////////////////

// Functions for generating date-time string (by Anne)
function GetTimeString() {
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	TimeString = "Date: "+DayNames[dayOfWeek]+" ";
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

function GetShortTimeString() {
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	ShortTimeString = "" + year + "";
	month = month + 1;
	if (month < 10) {ShortTimeString = ShortTimeString + "0";}
	ShortTimeString = ShortTimeString + month + "";
	if (dayOfMonth<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+dayOfMonth+"_";
	if (hour<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+hour+"h";
	if (minute<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+minute+"m";
	if (second<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+second+"s";
	return ShortTimeString;
}

// Find all specified files in subdirs:
function walkFiles(dir) {
	list = getFileList(dir);
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], "/"))
		   walkFiles(""+dir+list[i]);

		// Select only files with specified endsWith-string and process
		else  if (endsWith(list[i], imgfilter)) 
		   processImage(dir, list[i]);
	}
}

// Processing of images by RS-FISH with the defined parameters
function processImage(dirPath, imName) {
	
	open("" + dirPath + imName);
	
	// set what to delete
	var imNwoend = imName.replace("_maxZ.tif", "");
	
	results_csv_path = outdir + imNwoend + "_spots.csv";
	
	RSparams =  "image=" + imName + 
	" mode=Advanced anisotropy=" + anisotropyCoefficient + " robust_fitting=[" + ransac + "] use_anisotropy" + 
	" image_min=" + imMin + " image_max=" + imMax + " sigma=" + sigmaDoG + " threshold=" + thresholdDoG + 
	" support=" + supportRadius + " min_inlier_ratio=" + inlierRatio + " max_error=" + maxError + " spot_intensity_threshold=" + intensityThreshold + 
	" background=[" + bsMethod + "] background_subtraction_max_error=" + bsMaxError + " background_subtraction_min_inlier_ratio=" + bsInlierRatio + 
	" results_file=[" + results_csv_path + "]" + 
	" " + useMultithread + " num_threads=" + numThreads + " block_size_x=" + blockSizX + " block_size_y=" + blockSizY + " block_size_z=" + blockSizZ;
	print(RSparams);
	startTime = getTime();
	run("RS-FISH", RSparams);
	exeTime = getTime() - startTime; //in miliseconds

	// Save RS-FISH parameters to a txt file 
	txt_RS_parameters = File.open(outdir + imNwoend + "_parameters.txt");
	print(txt_RS_parameters, RS_parameters);
	File.close(txt_RS_parameters);

	// Close all windows:
	run("Close All");	
	while (nImages>0) { 
		selectImage(nImages); 
		close(); 
	
    } 
} 

//////////////////// RUN ///////////////////////

// Get start date-time string for log
datetimestringS = GetShortTimeString();

walkFiles(dir);

// Get end date-time string for log
datetimestringE = GetShortTimeString();

// Create and open a log file
logfile = File.open(outdir + datetimestringS + "_log.txt");
print(logfile, "####################################");
print(logfile, "Run started:\n" + datetimestringS);
print(logfile, "Analyzing data in " + dir);
print(logfile, "Processing files containing the following terminal string: " + "*" + imgfilter);
print(logfile, "Finished run at:\n" + datetimestringE);
print(logfile, "####################################");
print(logfile, getInfo("log"));
File.close(logfile);
run("Close All");
run("Collect Garbage");