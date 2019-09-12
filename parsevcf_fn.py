import vcf
import argparse
from collections import defaultdict
import pandas as pd

# example run:
# python CommonUniqueVariants.py -v vcf1.vcf -v vcf2.vcf -v vcf3.vcf -v vcf4.vcf -c cosmic.vcf 
vcf_variants = defaultdict()
vcf_positions = defaultdict()
def parseVCF(vcffile):
    pos = []
    var = []
    vcfread = vcf.Reader(open(vcffile))
    samples = vcfread.samples
    sampleName = vcffile
    if len(samples) == 0:
        print("no sample found, consolidating it as known mutations")
        sampleName += "_known_mutations"

    
    for record in vcfread:
        alts = record.ALT
        chrm = record.CHROM.replace("chr","")
        chrm = chrm.replace("Chr","")
        pos.append(chrm + "_" + str(record.POS))
        vcf_positions[sampleName] = set(pos)

        if "_known_mutations" in sampleName:
            alts = record.ALT
            for alt in alts:
                var.append("_".join([chrm, str(record.POS),record.REF, str(alt)]))
                vcf_variants[sampleName] = set(var)
        
        else:
            for sample in samples:
                sampleName += sample
                alts = record.genotype(sample).gt_bases
                try:
                    alts = alts.split("|")
                    if len(alts) == 1:
                        alts = record.genotype(sample).gt_bases.split("/")
                except AttributeError:
                    continue
                for alt in alts:
                    var.append("_".join([chrm, str(record.POS),record.REF, str(alt)]))
                    vcf_variants[sampleName] = set(var)
    print(vcf_variants)

parseVCF("sample1_final_filtered_snps.vcf")



