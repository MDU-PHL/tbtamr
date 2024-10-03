
from .CustomLog import logger
import logging
from mutamr import Fastq2vcf


def generatevcf(read1,read2,threads,ram,seq_id,keep,mindepth,minfrac,force,mtb,tmp) -> str:
    logger.info(f"You have provided reads as input - will now use mutAMR to generate vcf - please be patient.")
    V = Fastq2vcf.Fastq2Vcf(read1 = read1,
                read2= read2,
                threads=threads,
                ram = ram,
                seq_id= seq_id,
                keep = keep,
                mtb = mtb,
                mindepth = mindepth,
                minfrac =minfrac,
                force = force,
                tmp = tmp)
    vcf = V.run()
    

    return vcf