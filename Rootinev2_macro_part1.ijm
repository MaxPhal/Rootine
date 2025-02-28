///////////// Code for image stitching, greyscale normalization and preparation for filtering
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

/////// 2. Get file directory
name=getArgument(); // Retrieves the "name" argument from the bash file
//name="/media/phalempi/MP_04/image"; The file path can alternatively be given as string character (i.e. with quotation marks " "). 
// Note that with Windows, the "/" needs to be replaced with "\".

// Read parameters from vol_info txt file 
vol_info = ReadParameters(""+name+"_vol_info.txt");

/////// 3. Open the two raw bottom and top CT images
run("MHD/MHA...","open="+name+"_bottom.mhd");
rename("bottom");

run("MHD/MHA...","open="+name+"_top.mhd");
rename("top");

/////// 4. Image stitching, normalisation and saving the image in several parts prior filtering
// Remove overlapping slices in the two scans "bottom" and "top"
selectWindow("bottom");
run("Slice Keeper", "first=1 last="+vol_info[1]+" increment=1");
rename("bottom_nooverlap");
selectWindow("top");
getDimensions(width, height, channels, slices, frames);
run("Slice Keeper", "first="+vol_info[2]+" last="+slices+" increment=1");
rename("top_nooverlap");
//Image stitching of the two scans "bottom" and "top" without their overlap
run("Concatenate...", "  title=[Concatenated] image1=bottom_nooverlap image2=top_nooverlap image3=[-- None --]");
close("bottom");
close("top");

// Adapting the number of slices in the concatenated stack so that it can be split according to the number of layers defined.
selectWindow("Concatenated")
getDimensions(width, height, channels, slices, frames)
while (slices%vol_info[0]!=0){ // removes the top slices until the number of slices in the stack is divisable by nr_layers
run("Slice Remover", "first="+slices+" last="+slices+" increment=1");
selectWindow("Concatenated");
getDimensions(width, height, channels, slices, frames);
}

// Correcting for the grey value discontinuity present at the boundary of the two stitched images. 
// see https://imagejdocu.tudor.lu/plugin/stacks/attenuation_correction/start for more infos.
selectWindow("Concatenated");
getDimensions(width, height, channels, slices, frames);
run("Attenuation Correction", "opening=1 reference="+slices/2); // We use the slice in the middle of the stack as the reference slice.
selectWindow("Correction of Concatenated");
rename("corrected");
close("Concatenated");
close("Background of Concatenated");

// Preparation for filtering with the 3DNLM filter implemented in ITK. 
// This filter was originally implemented for MRI images with a Rician noise instead of Gaussian noise model.
// This noise model leads to a non-linear contrast stretching of low grey values in the filtered image.
// We circumvent this by adding a large constant offset in 16-bit (i.e. here 50000). 
// If a 3D NLM filter in another software (e.g. Avizo, Dragonfly, skicit-image) is used, this step can be commented out.

selectWindow("corrected");
run("Conversions...", " ");
run("16-bit");
run("Add...", "value=50000 stack");
getDimensions(width, height, channels, slices, frames);

// Saving the image in several parts. The number of parts created is equal to nr_layers as defined in the vol_info.txt file.
for(i=1;i<=vol_info[0];i++)
{
	if(i==1)bottom=1;
	else bottom=floor(slices*((i-1)/vol_info[0]))+1;

	if(i==vol_info[0]) top=slices;
	else top=floor(slices*(i/vol_info[0]));
	selectWindow("corrected");
	run("Duplicate...", "title=dummy duplicate range="+bottom+"-"+top);
	run("MHD/MHA ...", "save="+name+"_norm_"+IJ.pad(i, 2)+".mhd"); // here the "IJ.pad" pads 'i' with 1 a leading zero (image name becomes "_norm_0i"). 
}

/////// 5. Quit FIJI
run("Quit");

