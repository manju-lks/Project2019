import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pylab import savefig
import img2pdf 
from PIL import Image 
import os

def heatmap(sample_matrix):
	heat_df = pd.read_csv(sample_matrix, index_col=0)
	heat_df = heat_df.fillna(0)
	a4_dims = (10, 10)
	fig, ax = plt.subplots(figsize = a4_dims)
	heat_plot = sns.heatmap(ax = ax, data = heat_df, annot = True)
	heat_plot.figure.savefig("heat_plot.png",dpi=400.0)


def cosmic_frequency(cosmicOverlap):
	with open("cosmic_overlap.txt") as fh:
		for line in fh:
			cos_overlap = eval(line)
	cos_dict={}
	for key, val in cos_overlap.items():
    		cos_dict[key] = len(val)
	plt.figure(figsize=(10,10))
	plt.bar(list(cos_dict.keys()), cos_dict.values(), color= "c")
	plt.savefig("cosmic_overlap.png",dpi=400.0)


def fathmm_plot(fathmm):
	df = pd.read_csv(fathmm)
	df1 = df[df.Fathmm_prediction == 'PATHOGENIC']
	sns.distplot(df1['Fathmm_score'], hist = True, kde = False, label='pathogenic')
	df1 = df[df.Fathmm_prediction == 'NEUTRAL']
	sns.distplot(df1['Fathmm_score'], hist = True, kde = False, label='neutral')

	# Plot formatting

	plt.legend(prop={'size': 12})
	plt.title('Frequency of Variants w.r.t Fathmm Prediction and Score')
	plt.xlabel('Fathmm Score')
	plt.ylabel('Frequency of Variants')
	plt.savefig("fathmm_frequency.png",dpi=400.0)

def pdf(im_list):
	a4_list = []
	count = 0
	for img in im_list:
		count += 1
		img = Image.open(img)
		img = img.convert('RGB')
		img = img.resize((1000, 1000))
		a4im=Image.new("RGB",(1500, 1500),(255, 255, 255))
		a4im.paste(img, (50,50))
		if count < len(im_list):
			a4_list.append(a4im)
    
	a4im.save("output.pdf", save_all=True, append_images=a4_list)    
	

def main():
	
	fathmm_plot("fathmm_consolidation.csv")
	heatmap("sample_matrix.csv")
	cosmic_frequency("cosmic_overlap.csv")
	pdf(["heat_plot.png", "cosmic_overlap.png","fathmm_frequency.png"])
if __name__ == "__main__":
	main()
	

