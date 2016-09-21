import glob
import pandas as pd

def get_r1(sample, seq_dir):
    return(glob.glob(os.path.join(seq_dir, sample + "_*_R1_*.fastq.gz")))

def get_r2(sample, seq_dir):
    return(glob.glob(os.path.join(seq_dir, sample + "_*_R2_*.fastq.gz")))

def open_sample_sheet(sample_sheet_fp, lanes=False):
    sample_sheet = pd.read_excel(sample_sheet_fp, skiprows = 19, sheetname=1, header=19)

    if lanes:
        sample_sheet = sample_sheet.loc[df['Lane'].isin(lanes)]

    return(sample_sheet)



