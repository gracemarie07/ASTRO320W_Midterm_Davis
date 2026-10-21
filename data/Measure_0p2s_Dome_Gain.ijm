// ASTRO 320W: ten disjoint pairs of ORIGINAL 0.2-second dome flats.
// Run from the small AIJ toolbar: Plugins > Macros > Run...
// Choose the ORIGINAL Dome_Paired folder, then a folder for new outputs.
// Writes CSV/TXT only; never saves, renames, or overwrites any FITS image.
// Reopens originals; ignores any previously opened/edited image windows.
// No alignment, normalization, smoothing, or outlier replacement.
// ROI is zero-based, top-left (714,714), width=height=100 pixels.
// Every freshly opened flat is converted to 32-bit and has 100 DN subtracted.
// Gain is the photon-only estimator, NOT read-noise corrected; expect low-signal bias.

inputDir = getDirectory("Select ORIGINAL dome flats (Dome_Paired)");
outputDir = getDirectory("Choose folder for NEW gain results");
stamp = d2s(getTime(), 0);
outputCSV = outputDir + "dome_0p2s_gain_pairs_" + stamp + ".csv";
outputTXT = outputDir + "dome_0p2s_gain_summary_" + stamp + ".txt";
if (File.exists(outputCSV) || File.exists(outputTXT))
    exit("Output already exists. Run again to use a new timestamp.");

// Require exactly the known original sequence; unrelated files are ignored.
for (j=1; j<=20; j++) {
    path = inputDir + flatName(j);
    if (!File.exists(path))
        exit("Missing original file: " + flatName(j) + "\nChoose the original Dome_Paired folder.");
    header = File.openAsRawString(path, 28800);
    if (parseFloat(headerValue(header, "EXPTIME")) != 0.2 ||
        parseFloat(headerValue(header, "XBINNING")) != 2 ||
        parseFloat(headerValue(header, "YBINNING")) != 2 ||
        parseFloat(headerValue(header, "PEDESTAL")) != -100 ||
        indexOf(headerValue(header, "CALSTAT"), "D") < 0)
        exit("Unexpected or missing header setting in " + flatName(j) +
            "\nExpected 0.2 s, 2x2, PEDESTAL=-100, CALSTAT containing D.");
}

csv = "pair,file_1,file_2,roi_x,roi_y,roi_width,roi_height,pedestal_added_DN,F1_mean_DN,F1_SD_DN,F1_min_DN,F1_max_DN,F2_mean_DN,F2_SD_DN,F2_min_DN,F2_max_DN,difference_mean_DN,difference_SD_DN,difference_variance_DN2,photon_only_gain_e_per_DN,pair_level_difference_percent\n";
gains = newArray(10);
sumGain = 0;
setBatchMode(true);

for (p=0; p<10; p++) {
    n1 = 2*p+1;
    n2 = n1+1;
    id1 = openCorrected(inputDir + flatName(n1));
    getStatistics(area1, mean1, min1, max1, sd1);
    run("Select None");

    id2 = openCorrected(inputDir + flatName(n2));
    getStatistics(area2, mean2, min2, max2, sd2);
    run("Select None");

    // Image IDs ensure an older window with the same filename is not used.
    imageCalculator("Subtract create 32-bit", id1, id2);
    idDiff = getImageID();
    resetThreshold();
    makeRectangle(714, 714, 100, 100);
    getStatistics(areaD, meanD, minD, maxD, sdD);

    // Exact linearity check for the same ROI, allowing floating-point rounding.
    if (abs(meanD-(mean1-mean2)) > 0.005 || sdD <= 0)
        stopPair("Difference statistics failed consistency check.", id1, id2, idDiff);

    // Independently verified original-file reference pair. Stops if the reader/settings/data differ.
    if (p == 0) {
        if (abs(mean1-559.6255) > 0.01 ||
            abs(mean2-565.4219) > 0.01 ||
            abs(meanD+5.7964) > 0.01 ||
            abs(sdD-36.691320) > 0.01)
            stopPair("First pair does not match the verified original-file result.\n" +
                     "Measured F1="+d2s(mean1,6)+", F2="+d2s(mean2,6)+
                     ", difference SD="+d2s(sdD,6)+
                     "\nDo not apply any extra correction; send this message for review.",
                     id1, id2, idDiff);
    }

    variance = sdD*sdD;
    gain = (mean1+mean2)/variance;
    if (gain <= 0)
        stopPair("Nonpositive gain estimate; inspect the inputs.", id1, id2, idDiff);
    gains[p] = gain;
    sumGain = sumGain+gain;
    levelDifference = 100*abs(mean1-mean2)/((mean1+mean2)/2);

    csv = csv + n1 + "-" + n2 + "," + flatName(n1) + "," + flatName(n2) +
          ",714,714,100,100,-100," +
          d2s(mean1,9)+","+d2s(sd1,9)+","+d2s(min1,6)+","+d2s(max1,6)+","+
          d2s(mean2,9)+","+d2s(sd2,9)+","+d2s(min2,6)+","+d2s(max2,6)+","+
          d2s(meanD,9)+","+d2s(sdD,9)+","+d2s(variance,9)+","+
          d2s(gain,9)+","+d2s(levelDifference,6)+"\n";

    selectImage(idDiff); close();
    selectImage(id2); close();
    selectImage(id1); close();
    showProgress(p+1, 10);
}

meanGain = sumGain/10;
ss = 0;
for (p=0; p<10; p++) ss=ss+(gains[p]-meanGain)*(gains[p]-meanGain);
sampleSD = sqrt(ss/9);
sem = sampleSD/sqrt(10);
summary = "ASTRO 320W -- 0.2-second dome-flat gain check\n\n" +
          "Original input folder: " + inputDir + "\n" +
          "Frames: 20; disjoint pairs: 10\n" +
          "ROI: zero-based top-left x=714, y=714, width=100, height=100\n" +
          "Header checks: 0.2 seconds; 2x2; PEDESTAL=-100; dark-corrected\n" +
          "Processing: originals reopened; 32-bit conversion; subtract 100 DN once\n" +
          "No additional bias/dark subtraction, alignment, normalization or smoothing\n" +
          "First pair independently checked against original-file calculations\n\n" +
          "Estimator: g=(mean_F1+mean_F2)/(difference_SD^2)\n" +
          "Mean photon-only gain: " + d2s(meanGain,6) + " e-/DN\n" +
          "Pair-to-pair sample SD: " + d2s(sampleSD,6) + " e-/DN\n" +
          "Nominal SEM (SD/sqrt(10)): " + d2s(sem,6) + " e-/DN\n\n" +
          "LIMITS: This is a photon-only approximation with no read-noise\n" +
          "correction. Pair scatter/SEM exclude shared calibration errors,\n" +
          "illumination changes, detector nonlinearity and other systematics.\n" +
          "This is a 2x2 gain estimate; do not automatically transfer it to 1x1.\n" +
          "No image files or existing Results tables were overwritten.\n\n" +
          "CSV: " + outputCSV + "\n";
File.saveString(csv, outputCSV);
File.saveString(summary, outputTXT);
setBatchMode(false);
showMessage("Gain measurements finished", "Saved ten pair measurements and a summary:\n\n" +
            outputCSV + "\n\n" + outputTXT);

function flatName(n) {
    // Frame number occupies four digits, then '02' indicates the 0.2-second set.
    digits = "0000" + d2s(n,0);
    digits = substring(digits, lengthOf(digits)-4);
    return "DOME_FLAT_20_EXP-" + digits + "02.fit";
}

function headerValue(header, keyword) {
    for (offset=0; offset+80<=lengthOf(header); offset=offset+80) {
        card=substring(header, offset, offset+80);
        key=String.trim(substring(card,0,8));
        if (key == "END") return "";
        if (key == keyword) {
            value=substring(card,10);
            slash=indexOf(value,"/");
            if (slash>=0) value=substring(value,0,slash);
            value=replace(value,"'","");
            return String.trim(value);
        }
    }
    return "";
}

function openCorrected(path) {
    open(path);
    id=getImageID();
    getDimensions(w,h,channels,slices,frames);
    if (w!=1528 || h!=1528 || channels!=1 || slices!=1 || frames!=1 || bitDepth()!=16) {
        close();
        setBatchMode(false);
        exit("Expected one original 1528x1528, 16-bit image:\n"+path);
    }
    run("32-bit");
    run("Select None");
    resetThreshold();
    run("Subtract...", "value=100");
    makeRectangle(714,714,100,100);
    return id;
}

function stopPair(message, a, b, d) {
    selectImage(d); close();
    selectImage(b); close();
    selectImage(a); close();
    setBatchMode(false);
    exit(message);
}
