// Before starting, close everything that is open
close("*");

// Define which channels correspond to DAPI (405), Chromocenters (or other substructures (Substr)) and RNA spots (spots)
var channel405 = 5, channelSubstr = 4, channelSpots = 3;

// Create output directory
dir = getDirectory("Choose the directory of your experiment containing .tif files");
outdir = dir + "bg_and_ff_corrected/";
File.makeDirectory(outdir);

// Create log file
logfile = File.open(outdir + "log.txt");
print(logfile, "Analyzing data in " + dir);

// Get gain file and dark counts for flat-field correction
//FFCdir = getDirectory("Choose the directory which contains the files for flat-field correction");
FCCdir = "/Volumes/sd17B002/Microscopy-Dragonfly/Image-Analysis/Flat-Field-Correction/2022-05-16/";
FFCfiles = getFileList(FFCdir);

Dialog.create("Manual settings");
Dialog.addChoice("405 gain image:", FFCfiles);
Dialog.addChoice("dark counts:", FFCfiles);
Dialog.addSlider("405 channel", 1, 5, channel405);
Dialog.addSlider("substructure channel", 1, 5, channelSubstr);
Dialog.addSlider("spot channel", 1, 5, channelSpots);
Dialog.show();

gain405 = Dialog.getChoice();
dark = Dialog.getChoice();
channel405 = Dialog.getNumber();
channelSubstr = Dialog.getNumber();
channelSpots = Dialog.getNumber();

print(logfile, "Using " + gain405 + " for flatfield correction");
print(logfile, "Using " + dark + " for flatfield correction");
print(logfile, "Using channel " + channel405 + " as 405 channel");
print(logfile, "Using channel " + channelSubstr + " as substructure channel");
print(logfile, "Using channel " + channelSpots + " as spot channel");
print(logfile, "Flat-field uncorrected image for channel " + channel405 + " not saved separately");

// Create output directories
DAPIdir = outdir + "/405_ffc/";
substrdir = outdir + "/substructure/";
spotsdir = outdir + "/spot/";
otherdir = outdir + "/other/";
mergedir = outdir + "/merged/";
File.makeDirectory(DAPIdir);
File.makeDirectory(substrdir);
File.makeDirectory(spotsdir);
File.makeDirectory(otherdir);
File.makeDirectory(mergedir);

setBatchMode(true);
// load the images
for (i=0; i<filelist.length; i++) {
	if (endsWith(filelist[i], ".tif")) {
		image = filelist[i];
		im = dir + image;
		run("Bio-Formats Importer", "open=im autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		
		// convert minandmax for all channels
		Stack.getDimensions(width, height, channels, slices, frames);
			for(j=0; j<channels; j++){
				Stack.setChannel(j+1);
				setMinAndMax(0, 65535);}
		
		// split channels
		run("Split Channels");
		channelnbr = nImages;
		for (c = 1; c <= channelnbr; c++){
			selectImage(c);
			// run("Set Scale...", "distance=0 known=0 unit=pixel"); // for scale bars later: look up values for your objective
			rename("C" + c + "original");}
		
		// subtract dark counts from all channels
		open(FFCdir + dark);
		rename("dark");
		setMinAndMax(0, 65535);
		for (c = 1; c <= channelnbr; c++) {
			imageCalculator("Subtract create","C" + c + "original","dark");
			rename("C" + c);
			setMinAndMax(0, 65535);
			if (c == channel405) {
				
			} else if (c == channelsubstr) {
				saveAs("tiff", substrdir + "C" + c + "-" + image);
				rename("C" + c);
			} else if (c == channelspots) {
				saveAs("tiff", spotsdir + "C" + c + "-" + image);
				rename("C" + c);
			} else {
				saveAs("tiff", otherdir + "C" + c + "-" + image);
				rename("C" + c);
			}
		}
		
		// flat field correction for 405
		open(FFCdir + gain405);
		rename("gain");
		setMinAndMax(0, 65535);
		imageCalculator("Multiply create", "C" + ch_405, "gain");
		setMinAndMax(0, 65535);
		saveAs("tiff", DAPIdir + "C" + ch_405 + "-" + image);
		close("C" + ch_405);
		selectWindow("C" + ch_405 + "-" + image);
		rename("C" + ch_405);
		
		// save new merged image with the FFC DAPI
		if (channelnbr == 2) {
			run("Merge Channels...", "c1=C1 c2=C2 create keep");
		} else if (channelnbr == 3) {
			run("Merge Channels...", "c1=C1 c2=C2 c3=C3 create keep");
		} else if (channelnbr == 4) {
			run("Merge Channels...", "c1=C1 c2=C2 c3=C3 c4=C4 create keep");
		} else {
			run("Merge Channels...", "c1=C1 c2=C2 c3=C3 c4=C4 c5=C5 create keep");
		}
		Stack.setDisplayMode("grayscale");
		saveAs("tiff", mergedir + filelist[i]);
		close("*");
	}
}


setBatchMode(false);
