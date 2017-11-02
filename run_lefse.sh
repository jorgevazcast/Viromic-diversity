### Software
LEfSe_dir=~/software/nsegata-lefse-e3cabe93a0d1
### infiles
mkdir LEfSe
mv LEfSe.infile ./LEfSe
cd LEfSe
## The infile was created in R
infile=LEfSe.infile
# The current directory
indir=`pwd`
outir=$indir/out_lefse_LDA
mkdir $outir
# First row contain the class name for each sample "C" or "H"
class=1
# No subclasses
subclass=0
# Samples
samples=2

# Print the help and program options
python2.7 $LEfSe_dir/format_input.py -h
# format_input.py, create the LEfSe infile
# -o set the normalization value (default -1.0 meaning no normalization)
python2.7 $LEfSe_dir/format_input.py $indir/$infile $outir/$infile.format_input -c $class -s $subclass -u $samples -o 1000000

# Print the help and program options
python2.7 $LEfSe_dir/run_lefse.py -h 
# Performs the LEfSe analysis
# -w set the alpha value for the Wilcoxon test (default 0.05). In this case is going to be 0.00001
# -l set the threshold on the absolute value of the logarithmic LDA score (default 2.0). Change value to 3
# -s 1 set the multiple testing correction options
python2.7 $LEfSe_dir/run_lefse.py $outir/$infile.format_input -w 0.00001 -s 1 -l 3 $outir/$infile.res 

# Plot the barplot in a pdf format with a resolution of 600 dpi
python2.7 $LEfSe_dir/plot_res.py $outir/$infile.res $outir/$infile.pdf --dpi 600 --format pdf

# Plot the hierachical cladogram in a pdf format with a resolution of 600 dpi
python2.7 $LEfSe_dir/plot_cladogram.py $outir/$infile.res $outir/$infile.cladogram.pdf --dpi 600 --format pdf

# Optional: plot indiviudal barplots or those taxa that were statisticall diferent between conditions
# mkdir $outir/biomarkers_raw_images
# python2.7 $LEfSe_dir/plot_features.py $outir/$infile.format_input  $outir/$infile.res  $outir/biomarkers_raw_images/


