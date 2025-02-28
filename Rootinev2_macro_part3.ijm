////// Code for preprocessing, root segmentation and analysis - Rootine v.2

/////// 1. Define functions 
function ReadParameters(string) {
	filestring=File.openAsString(string);
	rows=split(filestring, "\n");
	param=newArray(rows.length);
	for(i=0; i<rows.length; i++){
	columns=split(rows[i],"\t");
	ColumnToRead=1;
	param[i]=parseFloat(columns[ColumnToRead]);}
	return param;
}

function GetExtendedRoiHist(image) {
	selectWindow("input");
	getDimensions(width, height, channels, slices, frames);
	nClas=256;
	ExtendedRoiHist=newArray(nClas);
	for(i=1;i<=slices;i+=10){
		setSlice(i);
		if(i<=slices/2)
		{
			frac=(i-1)/((slices/2)-1);
			x=frac*(ROI_param[4]-ROI_param[0])+ROI_param[0];
			y=frac*(ROI_param[5]-ROI_param[1])+ROI_param[1];
			w=frac*(ROI_param[6]-ROI_param[2])+ROI_param[2];
			h=frac*(ROI_param[7]-ROI_param[3])+ROI_param[3];
	
		}
		else
		{
			frac=(i-(slices/2))/(slices-(slices/2));
			x=frac*(ROI_param[8]-ROI_param[4])+ROI_param[4];
			y=frac*(ROI_param[9]-ROI_param[5])+ROI_param[5];
			w=frac*(ROI_param[10]-ROI_param[6])+ROI_param[6];		
		    h=frac*(ROI_param[11]-ROI_param[7])+ROI_param[7];		
		}
		w=w+extend_diameter; //extend the diameter to representatively capture the grey scale peak of the plastic wall of the pot
		h=h+extend_diameter; //extend the diameter to representatively capture the grey scale peak of the plastic wall of the pot 
		run("Specify...", "width="+w+" height="+h+" x="+x+" y="+y+" slice="+i+" oval centered");
		getHistogram(values, counts,nClas);
		for(j=0;j<nClas;j++){
				ExtendedRoiHist[j]+=counts[j];
		}
	}
	return ExtendedRoiHist;
}

function GetP1() {
	dummy = 0;
	for(i=1;i<histogram_cutoff;i++)
	{  
		if(ExtendedRoiHist[i]>dummy){
			p1=i;
			dummy=ExtendedRoiHist[i];
		}
	}
	return p1;
}
	
function GetP2() {
	dummy = 0;
	for(i=histogram_cutoff;i<255;i++)
	{	
		if(ExtendedRoiHist[i]>dummy){
			p2=i;
			dummy=ExtendedRoiHist[i];
		}
	}
	return p2;
}

function Createmask(image){
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	for(i=1;i<=slices;i++)
	{
		setSlice(i);
		if(i<=slices/2)
		{
			frac=(i-1)/((slices/2)-1);
			x=frac*(ROI_param[4]-ROI_param[0])+ROI_param[0];
			y=frac*(ROI_param[5]-ROI_param[1])+ROI_param[1];
			w=frac*(ROI_param[6]- ROI_param[2])+ ROI_param[2];
			h=frac*(ROI_param[7]-ROI_param[3])+ROI_param[3];
		}
		else
		{
			frac=(i-(slices/2))/(slices-(slices/2));
			x=frac*(ROI_param[8]-ROI_param[4])+ROI_param[4];
			y=frac*(ROI_param[9]-ROI_param[5])+ROI_param[5];
			w=frac*(ROI_param[10]-ROI_param[6])+ROI_param[6];		
		    h=frac*(ROI_param[11]-ROI_param[7])+ROI_param[7];		
		}
		run("Specify...", "width="+w+" height="+h+" x="+x+" y="+y+" slice="+i+" oval centered");
	    setBackgroundColor(0, 0, 0);
		run("Clear Outside", "slice");
	}
}

function computeRootDiameterDistribution_depth1(image) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	Stack.getStatistics(count, mean, min2, max2, std);
	nClas2=floor(max2+1)-floor(min2);
	cumul2=newArray(nClas2);
	for(i=1;i<=slices;i++){
		setSlice(i);
		getHistogram(values2, counts2,nClas2,floor(min2),floor(max2+1));
		for(j=0;j<nClas2;j++){
			cumul2[j]+=counts2[j];
				}
			}
	file2 = File.open(""+name+"_root_diameter_distribution_depth1.txt");
	for(j=0;j<nClas2;j++){
	print(file2,values2[j]+"	"+cumul2[j]);
	}
	File.close(file2);
}

function computeRootDiameterDistribution_depth2(image) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	Stack.getStatistics(count, mean, min2, max2, std);
	nClas2=floor(max2+1)-floor(min2);
	cumul2=newArray(nClas2);
	for(i=1;i<=slices;i++){
		setSlice(i);
		getHistogram(values2, counts2,nClas2,floor(min2),floor(max2+1));
		for(j=0;j<nClas2;j++){
			cumul2[j]+=counts2[j];
				}
			}
	file2 = File.open(""+name+"_root_diameter_distribution_depth2.txt");
	for(j=0;j<nClas2;j++){
	print(file2,values2[j]+"	"+cumul2[j]);
	}
	File.close(file2);
}

/////// 2. Get file directory
name=getArgument(); // Retrieves the "name" argument from the bash file
//name="/media/phalempi/MP_04/image"; The file path can alternatively be given as a string character (i.e. with quotation marks " "). 
// Note that under Windows, the "/" needs to be replaced with "\".

/////// 3. Read parameters 
//Read preprocessing and segmentation parameters from txt file and stores them in an array called "param"
param = ReadParameters(""+name+"_param.txt");

//Read parameters of the bounded ROI and stores them in an array called "ROI_param"
ROI_param = ReadParameters(""+name+"_ROI_param.txt");

//Read parameters from vol_info txt file and stores them in an array called "vol_info"
vol_info = ReadParameters(""+name+"_vol_info.txt");

/////// 4. Import and concatenate the denoised image parts
// Import 
for (i = 1; i <= vol_info[0]; i++) {
	run("MHD/MHA...","open="+name+"_norm_"+IJ.pad(i, 2)+"_nlm.mhd");  // here the "IJ.pad" pads 'i' with 1 leading zero (image name becomes "_norm_0i"). 
	rename(""+i+"");
}

// Concantenate
for (i = 1; i < vol_info[0]; i++) {
	run("Concatenate...", "  title=1 image1=1 image2="+(i+1)+" image3=[-- None --]");
}

// Here we subtract the greyvalue offset of 50000 that we created during the part 1 of the macro. 
// This is offset was added to circumvent the non-linear contrast stretching of low grey values generated by the 3DNLM we used.
// For more information, please refer to the paper. After subtraction of the constant offset, the image is converted to 8-bit again and saved. 

run("Subtract...", "value=50000 stack"); 
run("Conversions...", " ");
run("8-bit");
run("MHD/MHA ...", "save="+name+"_nlm_con.mhd");
rename("input");

/////// 5. Wall Detection
//Get extended Roi histogram
extend_diameter=50; // Number of voxels by which the bounded ROI will be extended
ExtendedRoiHist =  GetExtendedRoiHist("input");

// Retrieve P1 and P2
// In order to retrieve P1 and P2, we analyse the histogram in two different parts.
// Here we define the histogram cut-off as being the value right in the middle of the 8-bit grey value range (i.e. = 128). 
// With this value, we assume that both the grey value peaks of the pot wall (p1) and the soil matrix (p2) will have grey values < 128 and >= 128 respectively. 
// The functions "GetP1" and "GetP2" retrieve the corresponding grey values of the two peaks in the lowest part (i.e. GetP1) and in the highest part (i.e. GetP2) of the histogram;
histogram_cutoff=round(255/2); 
p1 = GetP1();
p2 = GetP2();

//Write P1 and P2 in a text file for future inspection
file = File.open(""+name+"_histogram_peaks.txt");
print(file, "First peak of the histogram (p1)	"+p1);
print(file, "Second peak of the histogram (p2)	"+p2);
File.close(file);

//Create the mask
selectWindow("input");
getDimensions(width, height, channels, slices, frames);
newImage("mask", "8-bit white", width, height, slices);
Createmask("mask"); // see the function in line 74 for more infos.

//Save the mask image for future inspection
run("MHD/MHA compressed ...", "save="+name+"_mask.mha");

//Remove the wall from the input image
imageCalculator("AND create stack", "mask","input");
rename("ROI");
close("original");
close("mask");
close("input");

/////// 6. Edge enhancement
selectWindow("ROI");
run("Unsharp Mask...", "radius="+param[0]+" mask="+param[1]+" stack");
rename("Unsharp");

/////// 7. Background removal via an absolute difference transform (ADT). 
//Calculate the root average grey value (v_r) and the threshold value of the ADT image (t_ADT);
v_r=round(param[2]*(p2-p1)+p1); 
t_ADT=255-(param[3]/2); 

//Computing the absolute difference transform of the unsharpened image. 
run("Macro...", "code=v=255-abs(v-"+v_r+") stack"); // ADT
rename("abs");
run("Duplicate...", "title=mask duplicate");
selectWindow("mask");
//Thresholding the ADT image.
setThreshold(0, t_ADT);
setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Dark black");
imageCalculator("Subtract create stack", "abs","mask");  // Note: Subtracting 255 from background voxels in abs will result in negative values and will thus be set to zero in 8-bit images, while the ADT values along roots remain unchanged (subtract zero)
rename("ADT");
close("abs");
close("mask");
close("Result of abs");

/////// 8. Root segmentation 
// Define root segmentation parameters
// Parameters common to both resolutions 
q=0.5; // the normalised smoothing strength of the tubeness filter. 
t_hyst_opt = floor(-194.55*pow(q,2)+534.78*q-139.27); // The optimal low threshold for hysteresis thresholding of the tubeness images. It is calculated according to the model regression in Fig. 4C in the manuscript.
t_hyst_high = 200 ; // The high threshold for hysteresis thresholding is kept constant and high enough for every tubeness resolutions and scales. 

// Root segmentation parameters at the original resolution
f_s_orig=1; // the resolution factor for the detection of fine roots at the original resolution.
sigma_orig = param[4]*q*f_s_orig; // sigma value used with the tubeness filter at the original resolution.

// Root segmentation parameters at the coarse resolution
f_s_coarse=0.5; // the resolution factor for the detection of coarse roots.
dr_inc = 4; // the increment at which root of increasing diameter will be detected.
d_r = newArray(((param[5]-(param[4] + dr_inc))/dr_inc)+1); // We create an array containing the targeted root diameter values.
sigma_coarse = newArray(d_r.length); // We create an array for the sigma values at the coarse scale. 
d_r[0] = param[4] + dr_inc; // smallest root diameter targeted at the coarse resolution
sigma_coarse[0] = d_r[0]*f_s_coarse*q; // smallest sigma value used at the coarse resolution
for (i = 1; i < d_r.length; i++) { // We calculate the root diameter target and the corresponding sigma values for each root diameter increment.
	d_r[i] = d_r[i-1] + dr_inc; // The root diameter target d_r is incremented by dr_inc at each iteration step
	sigma_coarse[i]= d_r[i]*f_s_coarse*q; // The sigma value is calculated for each targeted d_r value. 
}

/// 8.1 Detect fine roots at the original resolution
// First, we split the ADT image in several layers in order to not exceed the "Tubeness" plugin memory limitation (i.e. > 2Gb).
// This is crucial step and the number of layers should be defined accordingly so that (Stack_size / nr_layers) < 2 Gb. Otherwise, an error message will be thrown.
// Note that black slices are added at the bottom and at the top of each filtered layer. 
// Also note that the tubeness intensities decrease for slices close to the first and last slices in the stack. This decrease in tubeness intensities causes some
// some gap formation in the root system if it is not taken into account. To cope with that, we create an overlap in each filtered layer. Here we set this overlap equal to 8 slices 
// (e.g. the top of the first part will contain 8 of the first slices of the bottom of the second part and so forth). Note that we cannot create an overlap at the first and last slice of the concatenated stack.   

selectWindow("ADT");
getDimensions(width, height, channels, slices, frames); ; // We retrieve the number of slices in the ADT image
overlap=8; 
for(i=1;i<=vol_info[0];i++)
{
	if(i==1)bottom=1;
	else bottom=floor(slices*(i-1)/vol_info[0])+1; 
	
	if(i==vol_info[0])top=slices; 
	else top=floor(slices*i/vol_info[0])+overlap; 
	
	selectWindow("ADT");
	run("Duplicate...", "title=layer"+i+" duplicate range="+bottom+"-"+top);
}

//Run Tubeness plugin for each layer
for(i=1; i<=vol_info[0]; i++){
	selectWindow("layer"+i+"");
	run("Tubeness", "sigma="+sigma_orig+" use");
	close("layer"+i+"");
}

//After filtering with Tubeness, we remove the overlap prior to concatenation.
selectWindow("ADT");
getDimensions(width, height, channels, slices, frames); 
top=floor(slices/vol_info[0])+overlap; 
for(i=1;i<=vol_info[0];i++)
{
	selectWindow("tubeness of layer"+i);
	if(i!=vol_info[0]) run("Slice Remover", "first="+(top-overlap/2)+1+" last="+top+" increment=1");
	if(i!=1) run("Slice Remover", "first=1 last="+overlap/2+" increment=1");
}

//Concatenate the tubeness results of every layer
for (i = 1; i < vol_info[0]; i++) 
{run("Concatenate...", "  title=[tubeness of layer1] image1=[tubeness of layer1] image2=[tubeness of layer"+(i+1)+"] image3=[-- None --]");}

// Conversion of the 32-bit tubeness results to 8-bit 
run("Conversions...", "scale"); // We use the "scale" argument to normalize the grey value, i.e. the highest grey value in 32-bit is set to the grey value of 255 in 8-bit.
run("8-bit");
rename("Tubeness_orig");

//Hysterisis thresholding at the fine resolution
selectWindow("Tubeness_orig");
run("3D Hysteresis Thresholding", "high="+t_hyst_high+" low="+t_hyst_opt+" connectivity");
run("MHD/MHA compressed ...", "save="+name+"_fine_roots.mha"); // save the fine roots results for future inspection
rename("fine_roots"); 

/// 8.2 Detect coarse roots
//Downscale the image in order to reduce processing time without considerable loss of information for the coarse roots
selectWindow("ADT");
run("Scale...", "x="+f_s_coarse+" y="+f_s_coarse+" z="+f_s_coarse+" interpolation=None average process create");
rename("ADT_downscaled");
close("ADT");

//Segment the roots for every sigma value at the coarse resolution
for (i = 0; i < sigma_coarse.length; i++) {
	selectWindow("ADT_downscaled");
	run("Tubeness", "sigma="+sigma_coarse[i]+" use");
	run("Conversions...", "scale"); // here again, we use the "scale" argument to normalize the grey value, i.e. the highest grey value in 32-bit is set to the grey value of 255 in 8-bit.
	run("8-bit");
	run("3D Hysteresis Thresholding", "high="+t_hyst_high+" low="+t_hyst_opt+" connectivity");
	rename("roots_filtered_with_sigma"+sigma_coarse[i]+""); // Note that the sigma values are divided by 2 since rescaling is already accounted for in line 260.
    close("tubeness of ADT_downscaled");
}
close("ADT_downscaled");

// Merge roots segmented at every scale considered at the coarse resolution
opened_image = getList("image.titles"); // we retrieve the title of every opened image and create a list with those titles.
tubeness_coarse_scale=newArray(opened_image.length-1); // we create an array which will contain the name of the images created at every tubeness scale considered at the coarse resolution. 
j=0; // dummy variable
for (i = 0; i < opened_image.length; i++) 
{	print(opened_image[i]);
	if(startsWith(opened_image[i], "roots_filtered_with_sigma")){tubeness_coarse_scale[j]=opened_image[i];j++;} // we fill the "tubeness_coarse_scale" array with the name of the images starting with "roots_filtered_with_sigma" only.
}
for (i = 1; i < tubeness_coarse_scale.length; i++)
{
	imageCalculator("Max stack", ""+tubeness_coarse_scale[0]+"",""+tubeness_coarse_scale[i]+""); // All scales of the coarse resolution are stacked together using the "MAX" operator, i.e. a voxel is assigned to the root class if it is assigned to roots in at least one scale.
	close(""+tubeness_coarse_scale[i]+"");
}
rename("coarse_roots_downscaled");

//Upscale the image to its original resolution
selectWindow("coarse_roots_downscaled");
run("Scale...", "x="+(1/f_s_coarse)+" y="+(1/f_s_coarse)+" z="+(1/f_s_coarse)+" interpolation=None average process create");
rename("coarse_roots");
close("coarse_roots_downscaled");

//Combination of the fine and the coarse roots. Both resolutions results are stacked together using the "MAX" operator, i.e. a voxel is assigned to the root class if it is assigned to roots in at least one resolution.
imageCalculator("Max create stack", "coarse_roots","fine_roots");
rename("combined_roots");
run("MHD/MHA compressed ...", "save="+name+"_combined_roots.mha"); // save results for future inspection      
close("fine_roots");
close("coarse_roots");

/////// 9. False positives removal 
/// 9.1  Median filter. 
// This step smoothens the root surfaces and trims the oversegmentation voxels extending from the roots into the surroundings.
selectWindow("combined_roots");
run("Median 3D...", "x="+param[6]+" y="+param[6]+" z="+param[6]+"");                         
rename("trimmed_combined_roots");

/// 9.2  Keep connnected objects
// We ensure the connectivity of all root segments (e.g. the roots not directly emerging from the seed) by adding slices at the top of the stack. The number of added slices corresponds
// to the overlap that we created when splitting the layers during filtering with tubeness at the original (see line 269). The rationale is that the segmentation results for these slices are not 
// optimal due to the tubeness intensity decrease at the top of the image. 

// We duplicate the trimmed_combined_roots results as it will be used later on during the false negatives recovery step.
selectWindow("trimmed_combined_roots");
getDimensions(width, height, channels, slices, frames); 
run("Duplicate...", "title=trimmed_combined_roots_duplicate duplicate range=1-"+slices-(overlap*3)+""); // We already remove the overlap during duplication. Here the overlap value is multiplied by 3. See the explanation in line 400-402.

// We create a white image having the width and height of "trimmed_combined_roots" but a number of z slices equal to overlap.
newImage("whiteslices_for_keeplargest", "8-bit white", width, height, overlap); 
// Remove the top slices before adding the extra white slices
selectWindow("trimmed_combined_roots");
run("Slice Remover", "first="+slices-overlap+1+" last="+slices+" increment=1");
rename("trimmed_combined_roots_kept");
// We add the white slices at the top of "trimmed_combined_roots_kept"
run("Concatenate...", "  title=combined_root_whiteslice keep image1=trimmed_combined_roots_kept image2=whiteslices_for_keeplargest image3=[-- None --]");
close("trimmed_combined_roots_kept");
close("whiteslices_for_keeplargest");

//Keep Largest is applied in order to keep the connected roots only.
selectWindow("combined_root_whiteslice");
run("Connected Components Labeling", "connectivity=6 type=float"); // Connected Components Labeling with the argument "float" allows to have a number of labels equivalent to the range of 32-bit precision. 
rename("label");
run("Keep Largest Label"); // Keep largest identifies the largest label and keeps it (i.e. removes everything else).

// We remove the white slices at the top of the results of "keep largest".
// Important note. Some of the noise present in the image before Keep largest is connected to the white slices that we added at the top.
// Therefore, these artefacts are not removed during keep largest and are still present after removing the white slices at the top, if and when these artefacts extend in the z-direction. 
// This effect usually extends just a couple of slices downward. To get rid of this artefact, we remove an extra amount of slices by multiplying the overlap value by 3. 
run("Slice Remover", "first="+slices-(overlap*3)+1+" last="+slices+" increment=1");
rename("connected_roots");
close("label");
close("combined_root_whiteslice");

/////// 10. False negatives recovery 
// This step labels and evaluates every object unconnected to the root system and test whether it fulfills shape criteria which evoke typical shape of roots.
// We evaluate every unconnected object based on two criteria, i.e. its “Vesselness” and its size. 
imageCalculator("Subtract stack", "trimmed_combined_roots_duplicate","connected_roots"); // We subtract the connected roots in order to only evaluate the unconnected objetcs in the "trimmed_combined_roots" image.

// To evaluate the vesselness, we derive a “vesselness score” based on the analysis of the length of the semi-axes of fitting ellipsoids.
// The fitting ellipsoids are computed for every object with "Analyze Regions 3D". visit https://imagej.net/MorphoLibJ.html#Analyze_Regions_3D for more information.
// Note that the computed radii are sorted in decreasing order. We adpoted another notation in the manuscript. Corresponding radii are R1 = lambda3, R2 = lambda2 and R3 = lambda1.
selectWindow("trimmed_combined_roots_duplicate");
run("Connected Components Labeling", "connectivity=6 type=float");  // Connected Components Labeling with the argument "float" allows to have a number of labels equivalent to the range of 32-bit precision. 
run("Analyze Regions 3D", "equivalent_ellipsoid ellipsoid_elongations surface_area_method=[Crofton (13 dirs.)] euler_connectivity=C26"); // We compute the fitting ellipsoids for every object. visit https://imagej.net/MorphoLibJ.html#Analyze_Regions_3D for more information.
Table.rename("trimmed_combined_roots_duplicate-lbl-morpho", "Results"); //Command needed to convert the "trimmed_combined_roots_duplicate-lbl-morpho" table into a Results window from which getResult will work
for(i=0;i<nResults;i++)
{
	Rb=getResult("Elli.R3",i)/sqrt(getResult("Elli.R2",i)*getResult("Elli.R1",i)); // Calculation the geometrical ratio Rb to distinguish "blob-like" structures
	Ra=getResult("Elli.R2",i)/getResult("Elli.R1",i); // Calculation the geometrical ratio Ra to distinguish between a plate and a line
	setResult("Vesselness", i, exp(-Rb*Rb)*exp(-Ra*Ra)); // Calculation of the vesselness score for every object.
	flag=0;
	if(getResult("Vesselness", i)>param[7] && getResult("Elli.R1", i)>param[8])
	{
		flag=255; // If both the vesselness and the size score of an object are greater than the defined threshold, a flag (i.e. of value of 255) is assigned to this object.
	}
	setResult("Flag", i, flag);
}

run("Assign Measure to Label", "results="+"Results"+" column="+"Flag"+" min=0 max=255"); // combines a label image with a results table, and creates a new image for which each pixel/voxel is assigned the measurement value corresponding to the label it belongs to.
run("Conversions...", "scale");
run("8-bit");		
rename("false_negatives"); 
imageCalculator("Max create stack", "connected_roots","false_negatives"); // The false negatives are added to the connected roots.
rename("root_system"); 
close("trimmed_combined_roots_duplicate");
close("trimmed_combined_roots_duplicate-lbl");
close("false_negatives");
close("connected_roots");

//Saving the root segmentation results 
run("MHD/MHA compressed ...", "save="+name+"_root_system.mha");

// Duplicate the result for futher RSA analysis
selectWindow("root_system"); 
run("Duplicate...", "title=root_system duplicate");

//Saving the maximum projection in z-direction for a quick assessment of segmentation quality in a single 2D image 
run("Z Project...", "projection=[Max Intensity]");
saveAs("Jpeg", ""+name+"_Z_projection_root_system.jpg");

/////// 11. Determine root length
//Skeletonize root system
selectWindow("root_system");
run("Skeletonize (2D/3D)");
rename("skeleton");
getDimensions(width, height, channels, slices, frames); 
run("MHD/MHA compressed ...", "save="+name+"_skeleton.mha");

//Create a txt file in which the root length results will be stored.
file = File.open(""+name+"_root_length.txt");
print(file,"rootlength[pix]");
selectWindow("skeleton");
getDimensions(width, height, channels, slices, frames); 

//Analyse the skeleton for the two considered depth.
for (i = 3; i < vol_info.length; i=i+2) {
	bottom=vol_info[i]; 
	top=vol_info[i+1];
	selectWindow("skeleton");
	run("Duplicate...", "title=copy duplicate range="+bottom+"-"+top);
	total_branch_length=newArray(vol_info.length);
	run("Analyze Skeleton (2D/3D)", "prune=none show");
	selectWindow("Branch information");
	IJ.renameResults("Branch information","Results"); //Command needed to convert the Branch information window into a Results window from which getResult will work
	selectWindow("Results");
	sum=0;
	for (j = 0; j < nResults; j++) {
		sum+=getResult("Branch length",j);
		}
	total_branch_length[i]=sum;
	close("copy");
	close("Tagged skeleton");
	print(file,""+total_branch_length[i]);
}
File.close(file);

/////// 12. Determine root diameter distribution 
//Compute Root Thickness
selectWindow("root_system");
run("Geometry to Distance Map", "threshold=128");
close("root_system");
run("Distance Map to Distance Ridge");
close("root_system_EDT");
run("Distance Ridge to Local Thickness");
close("root_system_EDT_DR");
run("Local Thickness to Cleaned-Up Local Thickness");
close("root_system_EDT_DR_LT");
run("Subtract...", "value=1 stack"); // accouting for the hidden dilation step during LocThick. More information under : https://www.optinav.info/LocalThicknessEd.pdf
rename("root_thickness");

//Merge Skeleton with Thickness
selectWindow("skeleton");
run("Subtract...", "value=254 stack"); // accouting for the hidden dilation step during LocThick by trimming the results with the original root system
imageCalculator("Multiply create 32-bit stack", "root_thickness","skeleton");
rename("root_thickness_skel");
close("skeleton");
close("root_thickness");

//Retrieve the histogram for the two analysed depths
selectWindow("root_thickness_skel");
run("Duplicate...", "title=depth_1 duplicate range="+vol_info[3]+"-"+vol_info[4]+"");
computeRootDiameterDistribution_depth1("depth_1");
close("depth_1");

selectWindow("root_thickness_skel");
run("Duplicate...", "title=depth_2 duplicate range="+vol_info[5]+"-"+vol_info[6]+"");
computeRootDiameterDistribution_depth2("depth_2");
close("depth_2");

/////// 13. Quit ImageJ/Fiji
run("Quit");