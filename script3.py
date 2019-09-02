import vcf
from collections import Counter 
vcf_sample1_snp = vcf.Reader(open("sample1_final_filtered_snps.vcf","r"))
vcf_sample2_snp = vcf.Reader(open("sample2_final_filtered_snps.vcf","r"))
vcf_sample3_snp = vcf.Reader(open("sample3_final_filtered_snps.vcf","r"))
vcf_sample4_snp = vcf.Reader(open("sample4_final_filtered_snps.vcf","r"))
vcf_sample5_snp = vcf.Reader(open("sample5_final_filtered_snps.vcf","r"))
vcf_sample6_snp = vcf.Reader(open("sample6_final_filtered_snps.vcf","r"))
vcf_sample7_snp = vcf.Reader(open("sample7_final_filtered_snps.vcf","r"))

lis_vcfs=[vcf_sample1_snp, vcf_sample2_snp, vcf_sample3_snp, vcf_sample4_snp, vcf_sample5_snp, vcf_sample6_snp, vcf_sample7_snp]

res_list=[]
for vcff in lis_vcfs:
	for record in vcff:
		res_list.append([record.CHROM, record.POS])
counts=Counter(tuple(x) for x in res_list)
val=[]
for k in counts.keys():
	if counts[k] == 7:
		val.append(k)
print(val)
