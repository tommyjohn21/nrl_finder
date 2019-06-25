import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def finder(bedfile,out_format="maxes",offset=75,freq_cutoff=40,count_thr=10):
    '''
    The finder functon reads in a .bed file, extracts read lengths and computes NRLs from a filtered histogram of read lengths.
    
    Dependencies:
        numpy, scipy, matplotlib
    
    Recommended usage:
        nrls = nrl.finder("path/to/bedfile.bed")
    
    Required Args:
        bedfile: string denoting filepath to desired .bed file, e.g. "./demo/wt.bed"
    Optional Args:
        out_format: declares desired output, select from one of the strings below
          "maxes"  => (default) returns the read lengths of the computed NRL values (i.e. maxes)
          "mins"   => returns the read lengths of the computed local minima
          "auc"    => returns the computed area under the curve
          "figure" => prints a .png figure to the current directory that grahically shows each computed quantity
        offset: first read length to include for filtering (default 75)
        freq_cutoff: parameter that determines the low-pass cutoff frequency for filtering (default 40)
        count_thr: number of histogram counts needed to declare a true extrema (default: 10)
          note: to include all computed extrema, set count_thr to 0
    
    This code is provided free of charge and as-is. Reproduce, redistribute or modify it as you see fit. Use at your own risk.
     
    Dev: Tommy Wilson
    '''
    
    # Extract read lengths
    lengths=read_bed(bedfile)
    # Extract quantities of interest (nrls/maxes, mins, auc, etc)
    output=hunt_nrls(lengths,offset,freq_cutoff,count_thr)
    # Construct plot
    output=construct_plot(bedfile,output)
    # Format output for returning and return it
    output=format_output(bedfile,out_format,output)
    return output

def read_bed(bedfile):
    '''
    Read_bed extracts read lengths from a .bed file and ignores any headers
    '''
    
    # Initialize output
    lengths=[]
    # Read bedfile, cycle through lines
    with open(bedfile) as f:
        for num,line in enumerate(f,0):
            # Exclude headers                    
            if line.startswith("track") or line.startswith("browser"):
                continue
            array=line.split() # Parse line
            for i in range(1,3): # Assert that 2nd, 3rd columns are digits
                assert(array[i].isdigit()), "Value {} of line {} has a non-digit".format(i+1,num+1)
            lengths.append(int(array[2])-int(array[1])) # Collect length                
    # Output lengths
    return lengths

def hunt_nrls(lengths,offset,freq_cutoff,count_thr):
    '''
    Hunt_nrls computes a histogram of read lengths and low-pass filters a piece of that histogram (excluding offset) with a 6th-order Butterworth filter at the frequency defined by freq_cutoff. Extrema (minima and maxima) are then extracted, excluding those whose count number falls below a pre-defined threshold (given by count_thr).
    '''
    
    # Create histogram, extract counts
    bin_edges=np.arange(0.5,np.max(lengths)+1) # Define bin edges
    hist=np.histogram(lengths,bin_edges) # Compute counts per bin
    bin_edges=bin_edges[1:]-0.5 # Book-keeping
    counts=hist[0] # Take counts exclusively
    
    # Cut out offset
    counts_tf=counts[offset:] # Counts to filter
    bin_edges_tf=bin_edges[offset:] # Bin edges to filter
    
    # Generate and implement Butterworth filter
    b,a=signal.butter(6,freq_cutoff*1.0/1000) # Polynomials of IIR filter (force floating point for Python 2 compatibility)
    counts_f=signal.filtfilt(b,a,counts_tf) # Filtered counts
    bin_edges_f=bin_edges_tf # Book-keeping: 
    
    # Find minima and maxima of filtered data via second derivative test
    dcounts=np.diff(counts_f) # Numerical first derivative
    ddcounts=np.diff(dcounts) # Numerical second derivative
    crt_pt=np.where(np.append(np.array(False),np.multiply(dcounts,np.roll(dcounts,-1))<0))[0] # Identify critical points
    sdt=ddcounts[crt_pt-2] # Second derivative test
    crt_pt=crt_pt+offset # Book-keeping: index critical points to include offset
    
    # Extract extrema, area under the curve
    maxes=crt_pt[sdt<0] # Local maxima by second derivative test
    mins=crt_pt[sdt>0] # Local minima by second derivative test
    auc=np.sum(counts[mins[0]:]) # Area under the curve

    # Use only mins/maxes that have an adequate number of counts
    maxes=maxes[counts[maxes]>count_thr]
    mins=mins[counts[mins]>count_thr]

    # Organize outputs into dictionary and return
    output=dict()
    output["bin_edges"]=bin_edges
    output["counts"]=counts
    output["bin_edges_f"]=bin_edges_f
    output["counts_f"]=counts_f
    output["maxes"]=maxes
    output["mins"]=mins
    output["auc"]=auc
    return output
    
def construct_plot(bedfile,output):
    '''
    Construct_plot constructs a useful plot for visualization of computed quantities.
    '''
    
    f=plt.figure() # Initialize figure
    
    # Baseline plot images
    plt.fill_between(output["bin_edges"][output["mins"][0]:],output["counts"][output["mins"][0]:],color="#a8a8a8") # Plot AUC region
    plt.semilogy(output["bin_edges"],output["counts"],"k") # Plot counted histogram data
    plt.semilogy(output["bin_edges_f"],output["counts_f"],color="#d62728") # Plot filtered histogram data
    
    # Plot extrema
    find = lambda x,y=output["bin_edges_f"]: np.where(y==x)[0].tolist()[0] # Book-keeping: lambda to convert indices from non-filtered to filtered data
    plt.semilogy(output["bin_edges"][output["maxes"]],output["counts_f"][np.array(list(map(find,output["maxes"].tolist())))],"o",color="#1f77b4") # Maxes
    plt.semilogy(output["bin_edges"][output["mins"]],output["counts_f"][np.array(list(map(find,output["mins"].tolist())))],"o",color="#8F2EFF") # Mins

    # Adjust axes for viewability
    plt.xlim(0,round(output["bin_edges"][np.where(output["counts_f"]>1)[0][-1]]*1.1)) # Adjust x-axis for scaling
    plt.ylim(1,np.exp(np.log(np.max(output["counts"]))*1.1)) # Adjust y-axis for scaling

    # Plot labels
    plt.xlabel("Fragment length")
    plt.ylabel("Counts")
    plt.title("NRL Finder: "+bedfile)
    plt.legend(("Raw counts","Filtered counts","Maxes (NRLs)","Mins","Area under the curve"),loc="upper right")

    # Return with figure
    output["figure"]=f
    return output

def format_output(bedfile,out_format,output):
    '''
    Format_output determines which piece of information was queried by the user (option: out_format).
    '''
    
    # Call output by dictionary keys
    if out_format is "maxes" or \
      out_format is "mins" or \
      out_format is "auc":
        output=output[out_format]
    # If figure desired, print it and escape
    elif out_format is "figure":
        # Save figure in current directory with the bed file name
        output[out_format].savefig(bedfile[bedfile.rfind('/')+1:].replace('.bed','.png'))
        output=None
    # Throw error if out_format key is not one of the above
    else:
        raise KeyError("\"{}\" is not a valid out_format".format(out_format))
        
    return output
    