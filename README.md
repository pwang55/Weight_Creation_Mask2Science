

Usage:

In Script folder:
`$ python photutils_weight_science.py path_to_file/filename`

In Script folder, if running command line for loop, path to the files need to be provided so that the script can find the files:
```
$ for files in `cat correctedlist.txt`; do
> python photutils_weight_science.py $files path_to_files
> done
```

In files folder:
`$ python photutils_weight_science.py filename`

To change backsize, use option: backsize=64 or -backsize=64 after filename.

This script takes corrected fits image as input, use photutils to calculate background RMS map, convert it to weight map, 
and mask corrected image and weight image with mask, output to science frames. Lightstreaks will be detected and masked in the weight image.


------------------------------------------

Script Usage (this is the older version):

In the same directory as the script:
`$ ./create_weight.sh file_full_path/filename`

In the target file directory:
`$ script_path/./create_weight.sh filename`

If running using a FOR loop in the script directory, path to files may need to be provided:
```
$ for file in `cat filelist.txt`; do
> ./create_weight $file file_full_path
> done
```

If full_file_path/filename is given, second argument path will be ignored. Otherwise, 2nd argument path will be used to locate the filename given in 1st argument.


This script takes corrected_zwicky1953_4ne_0001.fits as input (example filename), creates corrected_zwicky1953_4ne_0001.wt.fits, then mask corrected_*fits with optical filtermask to create science_zwicky1953_4ne_0001.fits and .wt.fits


