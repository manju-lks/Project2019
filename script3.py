import vcf
from collections import Counter 
vcf_sample1_snp = vcf.Reader(open("sample1_final_filtered_snps.vcf","r"))
vcf_sample2_snp = vcf.Reader(open("sample2_final_filtered_snps.vcf","r"))
vcf_sample3_snp = vcf.Reader(open("sample3_final_filtered_snps.vcf","r"))
vcf_sample4_snp = vcf.Reader(open("sample4_final_filtered_snps.vcf","r"))
vcf_sample5_snp = vcf.Reader(open("sample5_final_filtered_snps.vcf","r"))
vcf_sample6_snp = vcf.Reader(open("sample6_final_filtered_snps.vcf","r"))
vcf_sample7_snp = vcf.Reader(open("sample7_final_filtered_snps.vcf","r"))


lis_vcfs=[vcf_sample1_snp, vcf_sample2_snp, vcf_sample3_snp, vcf_sample4_snp, vcf_sample5_snp,
          vcf_sample6_snp, vcf_sample7_snp]

#A list of lists with chromosome, positions, reference and alternate alleles

res_list=[]
for vcff in lis_vcfs:
	for record in vcff:
		res_list.append([record.CHROM, record.POS, record.REF,
                         ' '.join([str(elem) for elem in (record.ALT)])])

#count of occurence of each list
#if the count is equal to the number of input files they are taken as common ones 

n=len(lis_vcfs)       
counts=Counter(tuple(x) for x in res_list)
val=[]
for k in counts.keys():
	if counts[k] == n:
		val.append(k)
        
print("number of snps common in all ", n, " files is ", len(val), " and they are, \n", val)
