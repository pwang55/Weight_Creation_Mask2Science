#!/bin/bash

file=$1		# input file
path0=$2	# optional path to where the file is
BACKSIZE=64

usage='
Usage:

	In the same directory as the script:
	$ ./create_weight.sh file_full_path/filename

	In the target file directory:
	$ script_path/./create_weight.sh filename

	If running using a FOR loop in the script directory, path to files may need to be provided:
	$ for file in `cat filelist.txt`; do
	> ./create_weight $file file_full_path
	> done

	If full_file_path/filename is given, second argument path will be ignored. Otherwise, 2nd argument path will be used
	to locate the filename given in 1st argument.


This script takes corrected_zwicky1953_4ne_0001.fits as input (example filename), creates
corrected_zwicky1953_4ne_0001.wt.fits
Then mask corrected_*fits with optical filtermask to create
science_zwicky1953_4ne_0001.fits and .wt.fits


'

# If no argument was given at all, show docs and end the script.
if [ -z "$file" ] && [ -z "$path0" ]
then
        echo "$usage"
        exit 1
fi


# Get the absolute directory of this script so that it can find extra/files
script_dir=$(cd `dirname $0` && pwd)

# If argument 2 is provided and end is not /, add it; otherwise just use it as path. If argument 2 is not given, path would just be empty
if [ -n "$path0" ] && [ "${path0: -1}" != "/" ]; then
        path=$path0/
else
        path=$path0
fi
# If argument 1 is full path, reset path to empty
len_file=`echo $file | awk '{n=split($1,a,"/"); print n}'`
if (( $len_file > 1 )); then
        path=""
fi


# If path/file doesn't exist and path is not given, end the script
if [ ! -e $path$file ] && [ -z $path0 ] && (( $len_file == 1 )); then
        echo -e "File: \t $file \t can't be found, path argument may be needed, or filename is incorrect. Script end."
        exit 1
# If path/file doesn't exist and path is given, filename or path might be incorrect, end the script
elif [ ! -e $path$file ]; then
        echo -e "File: \t $path$file \t can't be found, check the file name or path, Script end."
        exit 1
fi


##############################################



# Get file name before .fits (including full path if full path to file is given)
base=`echo $path$file | sed -e 's/\.fits//g'`

sex $path$file -c ${script_dir}/extra/weight.config -PARAMETERS_NAME ${script_dir}/extra/test.param -FILTER_NAME ${script_dir}/extra/gauss_4.0_7x7.conv -WEIGHT_TYPE NONE -BACK_SIZE $BACKSIZE -CHECKIMAGE_TYPE BACKGROUND_RMS,SEGMENTATION -CHECKIMAGE_NAME ${base}.bgrms.fits,${base}.seg.fits


python ${script_dir}/rms2weight.py ${base}.bgrms.fits ${base}.seg.fits

rm ${base}.bgrms.fits
rm ${base}.seg.fits


python ${script_dir}/mask2science.py ${path}$file



