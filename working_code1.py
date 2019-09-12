import vcf
import argparse
from collections import defaultdict
import pandas as pd

# example run:
# python CommonUniqueVariants.py -v vcf1.vcf -v vcf2.vcf -v vcf3.vcf -v vcf4.vcf -c cosmic.vcf 
vcf_variants = {}
vcf_positions = {}
def parseVCF(vcffile):
    pos = []
    var = []
    vcfread = vcf.Reader(open(vcffile))
    samples = vcfread.samples
    sampleName = vcffile
    if len(samples) == 0:
        print("no sample found, consolidating it as known mutations")
        sampleName += "_known_mutations"
    else:
        for sample in samples:
            sampleName += sample
    
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

def parseCosmic(cosmicVCF):
	vcfreader = vcf.Reader(open(cosmicVCF))
	vcfreader._parse_info
	cosmicMutations = set()
	cosmicMutationScores = dict()
	for record in vcfreader: 
		for alt in record.ALT:
			chrm = record.CHROM.replace('chr', '')
			chrm = chrm.replace('Chr', '')
			identifier = "_".join([chrm, str(record.POS), record.REF, str(alt)])
			cosmicMutations.add(identifier)
			# Fill in fathmm score consolidation here
			# cosmicMutationScores[identifier] = vcfreader._parse_info(record.INFO) 
	return cosmicMutations, cosmicMutationScores


def common_mutations(vcf_positions, vcf_variants, cosmicMutations):
	samples = vcf_positions.keys()
	print("SAMPLES being compared:", samples)
	
	allCommon = set()
	commonVariants = {}
	cosmicOverlap = {}

	## Common across all samples:
	for sample in samples:
		if allCommon == set():
			allCommon = vcf_positions[sample]
		else:
			allCommon = allCommon.intersection(vcf_positions[sample])
	print("Common Across all samples:", allCommon)

	## Common variants between every 2 samples and with cosmic:
	for sample1 in samples:
		commonVariants[sample1] = {}
		for sample2 in samples:
			if sample1 == sample2:
				continue
			else:
				commonVariants[sample1][sample2] = len(vcf_positions[sample1].intersection(vcf_positions[sample2]))
		cosmicOverlap[sample1] = len(vcf_variants[sample1].intersection(cosmicMutations))

	print("COMMON VARIANTS matix:")
	print(pd.DataFrame(commonVariants).to_string())
	print("COSMIC OVERLAP:", cosmicOverlap)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcfs", action="append", help="list vcfs")
    parser.add_argument("-c", "--cosmic", help="cosmic vcf")
    parser.add_argument("-o", "--output_path", help="output path")

    args = parser.parse_args()
    for eachVcf in args.vcfs:
    	print("PARSING", eachVcf)
    	parseVCF(eachVcf)

    cosmicMutations, cosmicMutationScores = parseCosmic(args.cosmic)
    cvcfs_with_cosmic = common_mutations(vcf_positions, vcf_variants, cosmicMutations)

if __name__ == "__main__":
    main()
