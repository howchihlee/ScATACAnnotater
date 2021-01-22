import csv
import numpy as np
import pandas as pd

def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

def _check_bed_list(bed_list):
    return not(all((isinstance(x[0], str)) and (isinstance(x[1], int)) and (isinstance(x[2], int))) for x in bed_list)
    
def BedOverlap(bed1, bed2, overlapPercent):
    ## inputs: bed1, bed2, list of (chr, start, end)
    ## assumed bed1 and bed2 are sorted
    ## overlapPercent: float ranges from 0 - 1.
    ## output: a list of bed1 peaks that has over overlap higer than overlapPercent
    
    ### TODO change to raise error message
    if _check_bed_list(bed1): return "Invalid Input in Bed1"
    if _check_bed_list(bed2): return "Invalid Input in Bed2"    
    p0, p1 = 0, 0
    bed1Perc = [0.] * len(bed1)
    overlapPercent -= 1e-7 ## https://stackoverflow.com/questions/41690800/greater-than-or-equal-for-floats-failed
    
    while p0 < len(bed1) and p1 < len(bed2):
        while bed1[p0][0] != bed2[p1][0]: 
            if bed1[p0][0] < bed2[p1][0]:
                p0 += 1
            else:
                p1 += 1
            if(p0 >= len(bed1)):break
            if(p1 >= len(bed2)):break
            
        if(p0 >= len(bed1)):break
        if(p1 >= len(bed2)):break
            
        a0, b0, a1, b1 = int(bed1[p0][1]), int(bed1[p0][2]), int(bed2[p1][1]), int(bed2[p1][2])
        
        if (a0<b0<a1<b1):
            p0 += 1
        elif (a1<b1<a0<b0):
            p1 += 1
        elif (a0<=a1<=b0<=b1):
            bed1Perc[p0] += overlap(a0, b0, a1, b1) / (b0 - a0) # determine the percentage of overlap with bed1
            p0 += 1
        elif a1<=a0<=b0<=b1:
            bed1Perc[p0] += 1.
            p0 += 1    
        elif (a0<=a1<=b1<=b0) or (a0<=a1<=b1<=b0) or (a1<=a0<=b1<=b0):
            bed1Perc[p0] += overlap(a0, b0, a1, b1) / (b0 - a0) # determine the percentage of overlap with bed1
            p1 += 1
            
    return [p1 for p1, overlap in zip(bed1, bed1Perc) if overlap >= overlapPercent]
    #return [(p1, overlap) for p1, overlap in zip(bed1, bed1Perc)]

def CompBeds():
    delim = "\t"
    print("Enter file name for bed 1: ")
    bed1Name = input()
    print("Enter file name for bed 2: ")
    bed2Name = input()
    bed1 = BedScan(bed1Name, delim)
    bed2 = BedScan(bed2Name, delim)
    print("enter minimum overlap: ")
    minOverlap = int(input())
    result = BedOverlap(bed1, bed2, minOverlap)
    print(result)
    #This assumes that the bed files are sorted, perhaps using this command sort -k2 -n file > out
