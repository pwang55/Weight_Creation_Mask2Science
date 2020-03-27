

Usage:

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


