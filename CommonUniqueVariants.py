import vcf
import argparse
from collections import defaultdict
import pandas as pd

# example run:
# python CommonUniqueVariants.py -v vcf1.vcf -v vcf2.vcf -v vcf3.vcf -v vcf4.vcf -c cosmic.vcf 

vcf_variants = defaultdict(set)
vcf_positions = defaultdict(set)

def parseVCF(vcffile):
	vcfreader = vcf.Reader(open(vcffile))
	samples = vcfreader.samples
	sampleName = None
	if len(samples) == 0:
		print("No sample found, consolidating it as known_mutations")
		sampleName = vcffile + "_known_mutations"

	for record in vcfreader:
		alts = record.ALT
		chrm = record.CHROM.replace('chr', '')
		chrm = chrm.replace('Chr', '')

		if sampleName:
			if "_known_mutations" in sampleName:
				vcf_positions[sampleName].add(chrm + "_" + str(record.POS))
				alts = record.ALT
				for alt in alts:
					vcf_variants[sampleName].add("_".join([chrm, str(record.POS), record.REF, str(alt)]))

		else:
			for sample in samples:
				sampleName = vcffile + "_" + sample
				vcf_positions[sampleName].add(chrm + "_" + str(record.POS))
				alts = record.genotype(sample).gt_bases.split("|")
				if len(alts) == 1:
					alts = record.genotype(sample).gt_bases.split("/")
				for alt in alts:
					vcf_variants[sampleName].add("_".join([chrm, str(record.POS), record.REF, str(alt)]))

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
