import SortedNoDupeBedOverlap as bdO
import os
import yaml
import pandas as pd
from scipy.io import mmread

def read_bed_file(fileName, delimiter='\t', skiprows=0, header=None):
    ## default setting assumed input: tab separated, no header, no index
    ## assumed 0, 1, 2 columns are chromosome, start, end
    ## output: a list of peaks (chr:str, start:int, end:int)     
    
    data = pd.read_csv(fileName, delimiter = delimiter, 
                       skiprows = skiprows, header = header)
    output = data.values.tolist()
    out = list((row[0], int(row[1]), int(row[2])) for row in output)
    return out

def read_yaml(file):
    with open(file) as reffile: 
        info = yaml.load(reffile,Loader=yaml.FullLoader)
    return info

def read_reference_data(ref_folder, verbose = False):
    ## reading reference data
    
    refinfo = read_yaml(os.path.join(ref_folder, 'info.yaml'))
    bkg_file = os.path.join(ref_folder, refinfo['background'])
    fn_target = [(c, os.path.join(ref_folder, refinfo['cell_type'][c])) for c in refinfo['cell_type']]

    if (verbose):
        print('Using background file: ',refinfo['background'])
        print('With cell types:')
        for cell in refinfo['cell_type']:
            print (cell,':',refinfo['cell_type'][cell])    
    
    ref_bk_peak = read_bed_file(bkg_file)
    
    ref_subtype_peaks = []
    for c, f1 in fn_target:
        scan = set(read_bed_file(f1))    
        ref_subtype_peaks.append((c, scan))
        
    return ref_subtype_peaks, ref_bk_peak

def read_atac_count_csv(input_filename, verbose = False):
    ## read a count matrix in csv format (peak x cells)
    ## assumed the count matrix has a header and a index columm
    ## the header stores the cell barcodes
    ## the index columns stores open chromatin regions encoded as strings "chr_start_end"
    
    if(verbose): print("Reading Input File")

    data = pd.read_csv(input_filename, index_col = 0)
    #data.rename( columns={'Unnamed: 0':'id2peak'}, inplace=True)
    
    id2peak = []
    for c in data.index:
        temp = c.split('_')
        id2peak.append((temp[0], int(temp[1]), int(temp[2])))

    id2peak.sort( key = lambda x: (x[0], x[1]))    
    return data, id2peak

def read_atac_count_10x(input_folder, verbose = False):
    ## read a count matrix in 10x format 
    ## assumed three files in the input folder
    ## matrix.mtx: the count matrix (peaks by cells)is stored in the Matrix Market File Format
    ## barcodes.tsv: a tsv file stores the cell barcodes
    ## peaks.bed: stores open chromatin regions 
    ## return 
    ## mat: a coo sparse array (cells, peaks)
    ## id2bc: vector of cell barcodes (cells)
    ## id2peak: a list of intervals. Length = #peaks
    
    if(verbose): print("Reading Input File")

    mat = mmread(os.path.join(input_folder, 'matrix.mtx')).T
    id2bc = pd.read_csv(os.path.join(input_folder, 'barcodes.tsv'), header = None)[0].values
    id2peak = read_bed_file(os.path.join(input_folder, 'peaks.bed') )
    return mat, id2bc, id2peak

def read_input(input_filename, verbose = False):
    return read_atac_count_10x(input_folder, verbose = verbose)

def make_dataframe(values, index, columns):
    df = pd.DataFrame(values, 
                      columns = columns, 
                      index = index
                      )
    return df
    