import vcf
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

# example run:
# python CommonUniqueVariants.py -v vcf1.vcf -v vcf2.vcf -v vcf3.vcf -v vcf4.vcf -c cosmic.vcf 

vcf_variants = defaultdict(set)
vcf_positions = defaultdict(set)
#vcf_features = defaultdict(set)
#feature = {}

def parseVCF(vcffile):
	vcfreader = vcf.Reader(open(vcffile))
	samples = vcfreader.samples
	sampleName = None
	if len(samples) == 0:
		print("No sample found, consolidating it as known_mutations")
		sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_known_mutations"

	for record in vcfreader:
		alts = record.ALT
		chrm = record.CHROM.replace('chr', '')
		chrm = chrm.replace('Chr', '')
		info = record.INFO
		pos = chrm + "_" + str(record.POS)
		if sampleName and "_known_mutations" in sampleName:
				vcf_positions[sampleName].add(pos)
				alts = record.ALT
				#feature[sampleName]
				for alt in alts:
					vcf_variants[sampleName].add("_".join([pos, record.REF, str(alt)]))

		else:
			for sample in samples:
				sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_" + sample
				vcf_positions[sampleName].add(pos)
				try:
				    alts = record.genotype(sample).gt_bases.split("|")
				except:
					continue
				if len(alts) == 1:
					alts = record.genotype(sample).gt_bases.split("/")
				for alt in alts:
					vcf_variants[sampleName].add("_".join([pos, record.REF, str(alt)]))

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
		cosmicOverlap[sample1] = (vcf_variants[sample1].intersection(cosmicMutations))

	fathmm = {}
	df = pd.read_csv("fat1.csv",dtype=object)
	cos = df.to_dict("records")
	for val in cosmicOverlap.values():
		for value in val:
			values = value.split("_")
			for dicts in cos:
				if values[0] == dicts["chr"] and int(values[1]) >= float(dicts["pos_1"]) and int(values[1]) <= float(dicts["pos_2"]):
					fathmm[value] = [dicts["FATHMM_score"],dicts["FATHMM_prediction"],dicts["gene_name"]]

	print("COMMON VARIANTS matix:")
	sample_matrix = pd.DataFrame(commonVariants)
	print (sample_matrix.to_string())
	sample_matrix.to_csv("sample_matrix.csv")
	print("COSMIC OVERLAP:", cosmicOverlap)
	#cos_overlap = pd.DataFrame.from_dict(cosmicOverlap, orient='index')
	with open("cosmic_overlap.txt", "w") as txtfile:
		txtfile.write(str(cosmicOverlap))
	print("FATHMM SCORES: ", fathmm)
	fathmm_df = pd.DataFrame(fathmm)
	fathmm_df = fathmm_df.T
	fathmm_df.columns =["Fathmm_score", "Fathmm_prediction","Gene_Name"]
	fathmm_df["Fathmm_score"] = fathmm_df["Fathmm_score"].astype(float)
	fathmm_df.to_csv("fathmm_consolidation.csv")

	return (commonVariants, cosmicOverlap)

def unique_mutations(vcf_variants):
	samples = vcf_variants.keys()
	vcf_exclusives = {}

	## Unique variants in every sample
	for sample1 in samples:
		vcf_exclusives[sample1] = vcf_variants[sample1]
		for sample2 in samples:
			if sample1 == sample2:
				continue
			else:
				vcf_exclusives[sample1] = vcf_exclusives[sample1] - vcf_variants[sample2]

		fout = open(sample1 + "_uniqueVariants.txt", 'w')
		for variant in vcf_exclusives[sample1]:
			fout.writelines(variant + "\n")
		fout.close()

	return vcf_exclusives

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
    (commonVariants, cosmicOverlap) = common_mutations(vcf_positions, vcf_variants, cosmicMutations)
    vcf_exclusives = unique_mutations(vcf_variants)

if __name__ == "__main__":
    main()
