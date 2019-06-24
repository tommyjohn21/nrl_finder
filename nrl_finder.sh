#!/bin/bash

##########################
### Default parameters ###
##########################
out_format=maxes
offset=75
freq_cutoff=40
count_thr=10

#####################
### Input parsing ###
#####################
while getopts ":hb:o:y:f:t:" opt; do
  case ${opt} in
    h ) # help flag
			echo ''
			printf "usage: \tnrl_finder.sh -b </path/to/bedfile.bed> [-h] [-o <out_format>] \n\t\t[-y <offset>] [-f <freq_cutoff>] [-t <count_thr>]\n"
		  echo ''
			printf "The nrl_finder.sh bash script wrapper passes a .bed file name and\n"
			printf "additional parameters to the nrl.finder function located in bin/nrl.py.\n"
			printf "Outputs are given in tab-delimited form for use in piping.\n"
			printf "\nDependencies:\n\tpython, numpy, scipy, matplotlib\n"
    	printf "\nRequired Args:\n\t-b\tbedfile: string denoting filepath to desired .bed file,\n\t\te.g. \"./demo/wt.bed\"\n"
			printf "\nOptional Args:\n\t-o\tout_format: declares desired output, select from one\n\t\tof the strings below:"
			printf "\n\t\t\"maxes\"  => (default) returns the read lengths of\n\t\t\t    the computed NRL values (i.e. maxes)"
      printf "\n\t\t\"mins\"   => returns the read lengths of the computed\n\t\t\t    local minima"
      printf "\n\t\t\"auc\"    => returns the computed area under the curve"
      printf "\n\t\t\"figure\" => prints a .png figure to the current\n\t\t\t    directory that grahically shows each computed\\n\t\t\t    quantity"
	    printf "\n\t-y\toffset: first read length to include for filtering\n\t\t(default 75)"			
			printf "\n\t-f\tfreq_cutoff: parameter that determines the low-pass cutoff\n\t\tfrequency for filtering (default 40)"
			printf "\n\t-t\tcount_thr: number of histogram counts needed to declare\n\t\ta true extrema (default: 10), note: to include all\n\t\tcomputed extrema, set count_thr to 0"
			printf "\n\nThis code is provided free of charge and as-is. Reproduce, redistribute\nor modify it as you see fit. Use at your own risk."
	    printf "\n\nDev: Tommy Wilson"
			echo -e "\n"
			exit 0
      ;;
		b ) # bed file input
			bedfile=$OPTARG
			;;
		o ) # output format
			out_format=$OPTARG
			;;
		y ) # offset
			offset=$OPTARG
			;;
		f ) # frequency cutoff for Butterworth filter
			freq_cutoff=$OPTARG
			;;
		t ) # count threshold
			count_thr=$OPTARG
			;;
    \? ) printf "usage: \tnrl_finder.sh -b </path/to/bedfile.bed> [-h] [-o <out_format>] \n\t\t[-y <offset>] [-f <freq_cutoff>] [-t <count_thr>]\n"
			exit 1
			;;
  esac
done

####################
### Input checks ###
####################
# Return usage if no input is supplied
if [ $OPTIND -eq 1 ]
then 
	printf "usage: \tnrl_finder.sh -b </path/to/bedfile.bed> [-h] [-o <out_format>] \n\t\t[-y <offset>] [-f <freq_cutoff>] [-t <count_thr>]\n"
	printf "\nUse -h (help) flag for additional information.\n"
  exit 1
fi

# Ensure bed file is supplied
if [ ! "$bedfile" ]
then
    echo -e "Error! Bedfile must be supplied with option -b.\nUse -h (help) flag for additional information."
		printf "usage: \tnrl_finder.sh -b </path/to/bedfile.bed> [-h] [-o <out_format>] \n\t\t[-y <offset>] [-f <freq_cutoff>] [-t <count_thr>]\n"
    exit 1
fi
# Ensure that a valid output option is chosen
if [ "$out_format" != "maxes" ] && [ "$out_format" != "mins" ] && [ "$out_format" != "auc" ] && [ "$out_format" != "figure" ] 
then
    echo "Error! Output format in option -o must be one of: maxes, mins, auc, figure."
    exit 1
fi
# Assert offset is numeric
if ! [[ $offset =~ ^[0-9]+$ ]]
then
  echo "Error! Offset parameter specified by -y should be numeric."
  exit 1
fi
# Assert frequency cutoff is numeric
if ! [[ $freq_cutoff =~ ^[0-9]+$ ]]
then
  echo "Error! Frequency cutoff parameter specified by -f should be numeric."
  exit 1
fi
# Assert count threshold is numeric
if ! [[ $count_thr =~ ^[0-9]+$ ]]
then
  echo "Error! Count threshold parameter specified by -t should be numeric."
  exit 1
fi

#########################
### Callout to python ###
#########################
# Use the above parameters
# Output data in tab-delimited format
output=$(python -c "import bin.nrl as nrl; print(nrl.finder(\"$bedfile\",\"$out_format\",$offset,$freq_cutoff,$count_thr))")

#####################
### Format output ###
#####################
# Remove None return from Python output if returned
if [ "$output" == "None" ]
then
	output=
else
	# Format as tab-delimited string without brackets
	output=$(echo "$output" | tr " " \\t | tr -d [])
fi 


#####################
### Return output ###
#####################
# Print to terminal
if [ ! -z "$output" ]
then
	echo -e "$output"
fi
exit 0




