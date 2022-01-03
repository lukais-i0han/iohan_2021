#!/usr/bin/env python
# coding: utf-8

# Necessary libraries
from IPython.display import clear_output

import synapseclient
import pandas as pd
import subprocess
import os.path

# Set pathway
local = "/data4/lukais.iohan/5XFAD"

tablePath = local + "/metadata/fastq.csv"
df = pd.read_csv(tablePath)
df = df.sort_values(by=["name"])
samples = df.specimenID.unique()


fqPath = local + "/fastq/"
klPath = local + "/kallisto/"
os.popen('mkdir '+fqPath) # Create fastq Path
os.popen('mkdir '+klPath) # Create kallisto output path

# Directory to builded kallisto index
kalRef = "/data4/lukais.iohan/5XFAD/refs/Mus_musculus.GRCm38.re99.idx"

#### Synapse login
syn = synapseclient.Synapse()

syn.login('lukais.iohan', '338486319ccil')

def runKallisto(samples,path):
    
    # Create kallisto output dir
    kalDir = klPath + samples
    # Create directory
    os.popen('mkdir '+kalDir)
    
    # Import ids
    
    print("\nSamples:")
    name1, name2 = df.name[df.specimenID==samples]
    r1, r2 = df.id[df.specimenID==samples]
    print("Sample.1: "+name1)
    print("r1: "+r1)
    print("Sample.2: "+name2)
    print("r2: "+r2)
    
    if os.path.isfile(kalDir+'/abundance.tsv'):
        print("Sample " + samples + " was already processed!")
        pass
    
    else:
        print("Downloading " + samples + " ...")
        file1 = syn.get(r1, downloadLocation=fqPath) # Download file1 to kallisto
        file2 = syn.get(r2, downloadLocation=fqPath) # Download file2 to kallisto

        # Run Kallisto
        print("Running kallisto on " + samples + " ...")
        subprocess.call('kallisto quant -i '+
                 kalRef + ' -t 20 -o '+ kalDir+ ' '+
                file1.path+' '+file2.path, shell=True)

        # Remove files
        os.popen('rm '+file1.path)
        os.popen('rm '+file2.path)
    
for i in samples:
    runKallisto(samples = i, path = local)


