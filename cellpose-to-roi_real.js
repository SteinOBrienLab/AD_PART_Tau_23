importClass(Packages.ij.IJ);
importClass(Packages.ij.plugin.frame.RoiManager);
importClass(Packages.ij.io.DirectoryChooser);

// Get directory
var dc = new DirectoryChooser("Choose the directory containing cellpose_seg.npy");
var dir = dc.getDirectory();

// Get ROI Manager
var rm = RoiManager.getInstance();
if (rm == null) {
    rm = new RoiManager();
}
rm.reset();

// Import NPY file using Bio-Formats
var path = dir + "_seg.npy";
IJ.run("Bio-Formats Importer", "open=[" + path + "]");

// Get current image
var imp = IJ.getImage();

// Convert to 8-bit if needed
if (imp.getBitDepth() != 8) {
    IJ.run(imp, "8-bit", "");
}

// Create binary masks for each cell
IJ.setAutoThreshold(imp, "Default dark");
IJ.run(imp, "Convert to Mask", "");

// Find and analyze particles
IJ.run(imp, "Analyze Particles...", "size=10-Infinity pixel show=Nothing clear include add");

// Save ROIs
rm.runCommand("Show All");
rm.save(dir + "CellposeROIs.zip");

// Create labeled image
IJ.newImage("Labels", "8-bit black", imp.getWidth(), imp.getHeight(), 1);
rm.runCommand("Fill");
rm.runCommand("Show All");

// Clean up
IJ.run("Select None");

IJ.log("Conversion complete!");
IJ.log("Number of ROIs created: " + rm.getCount());
