#!/usr/bin/python3
'''

Author	  : Sonal Rashmi / Sanket
Date		: 17/07/2020 / 08/01/2020
Description : It takes the IPD directory post successful completion of IPD run, generates a html report with Basic Alignment stats(using picard), Coverage of SARS-CoV2(using samtools), abundance plots (Overall and SARS-CoV2), Novel variant and Clade assessment
version : 2
'''

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
from globals import *
from gisaidvcfannotator import *
import subprocess
from markdown import markdown
from localblastn import *
from gisaidmetadataparser import *
from variantcladeassessment import *

'''
Input to the object: Output Directory of IPD
'''

class CustomError(Exception):
	def __int__(self,value):
		self.value=value
	def __str__(self):
		return repr(self.value)
class CoV2ReportGenerator(object):
	def __init__(self, outdir=None):
		GlobalVar.initialize()
		if outdir:
			self.outdir=outdir
			self.vcf_map={}
			self.countfile_map={}
			self.total_fragment_list=[]
			self.coverage_map={}
			self.summaryout_map={}
			self.contig_map={}
			self.cov2outdir=outdir+"cov2output/"
			if not os.path.exists(self.cov2outdir):
				os.mkdir(self.cov2outdir)
			for file in os.listdir(self.outdir):
				if file.endswith("_count.bam"):
					bam=os.path.join(self.outdir, file)
					bam_sample=os.path.basename(bam)
					sample=bam.replace("_count.bam","")
					summary_file_out=self.bam_summary_generation(bam)
					if summary_file_out:
						self.summaryout_map[sample]=summary_file_out
					coverage_file_out=self.bam_coverage(bam)
					if coverage_file_out:
						self.coverage_map[sample]=coverage_file_out
				if file.endswith("_finalcounts.tsv"):
					countfile=os.path.join(self.outdir, file)
					countfile_sample=os.path.basename(countfile)
					sample=countfile_sample.replace("_finalcounts.tsv","")
					self.countfile_map[sample]=countfile

				if file.endswith("_featurecounts.tsv.summary"):
					featurecountsummary=os.path.join(self.outdir, file)
					data=pd.read_csv(featurecountsummary,sep="\t",index_col = "Status")
					data_t=pd.DataFrame.transpose(data)
					total=data_t['Assigned'].to_list()[0]
					self.total_fragment_list.append(total)
					#print(self.total_fragment_list)

				if file.endswith("_final_annotated.vcf"):
					vcf=os.path.join(self.outdir, file)
					vcf_sample=os.path.basename(vcf)
					sample=vcf_sample.replace("_final_annotated.vcf","")
					self.vcf_map[sample]=vcf
			#where clade assignment happens
			print(self.vcf_map)
			self.varcladeassessmentobj_=VariantCladeAssessment(self.vcf_map)
			#self.clade_assessment_obj=GisaidVcfAnnotator(self.vcf_map)
		else:
			print("Input directory %s not found!" %(outdir))
			sys.exit(0)

	def bam_summary_generation(self, bam=None):
		summary_file=None
		if bam:
			if os.path.isfile(bam):
				#sort bam
				sortedbam=bam.replace(".bam","_sorted.bam")
				cmd=GlobalVar.picard_+" SortSam INPUT="+ bam +" OUTPUT="+sortedbam+" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR="+os.path.join(self.outdir, "tmp")
				cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
				cprocess.check_returncode()

				#bam indexing
				cmd=GlobalVar.picard_+" BuildBamIndex INPUT="+ sortedbam +" OUTPUT="+sortedbam+".bai VALIDATION_STRINGENCY=SILENT  TMP_DIR="+os.path.join(self.outdir, "tmp")
				cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
				cprocess.check_returncode()

				#picard summary
				summaryout=sortedbam.replace("_count_sorted.bam","_summary.tsv")
				cmd=GlobalVar.picard_+" CollectAlignmentSummaryMetrics R="+ GlobalVar.hspathofa_ +" I="+ sortedbam +" O="+summaryout
				cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
				cprocess.check_returncode()

		if os.path.isfile(summaryout):
			summary_file=summaryout
		else:
			print("Error in generating summary file "+bam)
		return summary_file

	def bam_coverage(self, bam=None):
		coverage_out=None
		if bam:
			if os.path.isfile(bam):
				#samtools depth
				sortedbam=bam.replace(".bam","_sorted.bam")
				coverage=sortedbam.replace("_count_sorted.bam","_coverage.tsv")
				cmd=GlobalVar.samtools_+" depth -r NC_045512.2:1-29903 -a "+ sortedbam + " -o "+coverage
				cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
				cprocess.check_returncode()
		if os.path.isfile(coverage):
			coverage_out=coverage
		else:
			print("Error in generating coverage for "+bam)
		return coverage_out

	def basic_stats_tabulation_for_file(self,file_name,sample):
		data=open(file_name,'r')
		Pandasobj=pd.read_csv(file_name,skiprows=5,nrows=3,sep="\t",engine='python',header=1)
		PandasDf=pd.DataFrame(data=Pandasobj)
		col = ['TOTAL_READS','PF_READS_ALIGNED','PCT_PF_READS_ALIGNED','MEAN_READ_LENGTH']
		if PandasDf["CATEGORY"].shape[0] == 3:
			PandasDf=PandasDf.loc[PandasDf["CATEGORY"] == "FIRST_OF_PAIR"]
			PandasDf=PandasDf.rename(index={2:0})
		elif PandasDf['CATEGORY'].shape[0] == 1:
			PandasDf=PandasDf.loc[PandasDf["CATEGORY"] == "UNPAIRED"]
		df = pd.DataFrame(PandasDf, columns=col)
		df['Percent_Aligned_Reads']=df['PCT_PF_READS_ALIGNED']*100
		df=df.drop(['PCT_PF_READS_ALIGNED'],axis=1)
		df=df.rename(columns={"TOTAL_READS":"Total_Reads","PF_READS_ALIGNED":"Aligned_Reads","Percent_Aligned_Reads":"Percent_Aligned_Reads","MEAN_READ_LENGTH":"Mean_Read_Length"})
		df=df[["Total_Reads","Aligned_Reads","Percent_Aligned_Reads","Mean_Read_Length"]]
		df_t = pd.DataFrame.transpose(df).rename(columns={0:os.path.basename(sample)})
		return df_t

	def get_basic_stats_tabulation_for_the_batch(self):
		output_file=None
		all_sample_df = pd.DataFrame()
		cnt=0
		for sample in self.summaryout_map:
			if self.summaryout_map[sample]:
				file=self.summaryout_map[sample]
				data=open(file,'r')
				temp_df=self.basic_stats_tabulation_for_file(file,sample)
				if cnt == 0:
					all_sample_df = temp_df
				else:
					all_sample_df = all_sample_df.merge(temp_df, left_index=True, right_index=True)
				cnt=cnt+1
		output_file=self.cov2outdir+"Basic_stats_summary.csv"
		all_sample_df.to_csv(output_file)
		if os.path.isfile(output_file):
			print("Basic Stats summary tabulation generated!")
		else:
			print("Error in basic Stats summary output generation.")
		return output_file

	def get_coverage_plot(self):
		image_file_list=[]
		coverage_out=None
		#coverage_out=open(os.path.join(self.cov2outdir,"cov2_coverage_compilation.csv"),'a+',encoding='UTF-8')
		coverage_out=open(self.cov2outdir+"cov2_coverage_compilation.csv",'a+',encoding='UTF-8')
		for sample in self.coverage_map:
			if self.coverage_map[sample]:
				x = []
				y = []
				table = pd.read_csv(self.coverage_map[sample] , sep='\t',header=None)
				alignment = pd.DataFrame(data=table)
				x=alignment[1]
				y=alignment[2]
				axes = plt.gca()
				axes.set_ylim([0,500])

				median=np.median(y)
				ss=sample.split("/")
				sample=ss[len(ss)-1]
				#base name obtained - change Nov 5,2020 - Sanket
				title_var= sample+" (Median Coverage = "+str(median)+" )"
				plt.plot(x,y,color="grey")
				plt.xlabel('Position in SARS CoV2 Genome')
				plt.ylabel('Depth of Coverage')
				plt.title(title_var, loc='center')
				plt.savefig(self.cov2outdir+sample+"_coverage.png"))
				plt.close(fig=None)
				image_file_list.append(self.cov2outdir+sample+"_coverage.png")
				print(str(sample)+"\t"+str(median),file=coverage_out)
		print("Coverage output generated!!!")
		return image_file_list

	def get_abundance_plot(self):
		human=[]
		CoV2=[]
		pathogen=[]
		CoV2_fpkm=[]
		pathogen_all=[]
		sample_list=[]
		image_file_list=[]
		cov2_fpkm_out=open(self.cov2outdir+"cov2_fpkm_compilation.csv",'a+',encoding='UTF-8')
		for sample in self.countfile_map:
			file=self.countfile_map[sample]
			data=open(file,'r')
			table = pd.read_csv(file, sep='\t')
			count = pd.DataFrame(data=table)
			human_df=count[count['Feature'].str[0:4].isin(['ENSG'])].round(2)
			human.append(round(human_df['Fragments'].sum(),2))
			CoV2_df=count[count.ID == "NC_045512.2"]
			CoV2.append(round(CoV2_df["Fragments"],2).to_string(index=False))
			CoV2_fpkm.append(round(CoV2_df["FPKM"],2).to_string(index=False))
			pathogen_df=count[~count['Feature'].str[0:4].isin(['ENSG'])]
			pathogen_all.append(round(pathogen_df['Fragments'].sum(),2))
			sample_list.append(sample)
			print(str(sample)+"\t"+str(round(CoV2_df["FPKM"],2).to_string(index=False)),file=cov2_fpkm_out)

		print("Abundance Summary generated!")

		CoV2=list(map(float, CoV2))
		CoV2_fpkm=list(map(float, CoV2_fpkm))
		Total=self.total_fragment_list
		pathogen=[x1 - x2 for (x1, x2) in zip(pathogen_all, CoV2)]
		Unaligned=[round(x4 - (x1 + x2 + x3),2) for (x1, x2, x3, x4) in zip(human, pathogen, CoV2, Total)]

		#Plot FPKM
		CoV2_fpkm_log2=[]
		for i in CoV2_fpkm:
			if i > 0:
				x=np.log2(i)
			else:
				x=i
			CoV2_fpkm_log2.append(x)
		sample=sample_list
		log2FPKM=CoV2_fpkm_log2
		#print(sample,log2FPKM)
		sample_pos = [i for i, _ in enumerate(sample)]
		fpkm_plot_file=self.cov2outdir+'CoV2_FPKM.png'
		plt.bar(sample_pos, log2FPKM, color='red', width=1)
		plt.xlabel('')
		plt.ylabel('CoV2 log2 FPKM')
		plt.title("SARS-CoV-2 quantification (FPKM)")
		plt.xticks(sample_pos, sample, rotation=90)
		plt.savefig(fpkm_plot_file, bbox_inches = 'tight', pad_inches = 0.1)
		#plt.savefig(fpkm_plot_file, bbox_inches = 'tight')
		plt.close(fig=None)

		print("FPKM plot is generated")

		#Plot Stack bar plot
		Sample = sample_list
		CoV2= np.array([(x/y)*100 for x, y in zip(map(float, CoV2), map(int, Total))])
		Human= np.array([(x/y)*100 for x, y in zip(map(float, human), map(int, Total))])
		Pathogen= np.array([(x/y)*100 for x, y in zip(map(float, pathogen), map(int, Total))])
		Unaligned= np.array([(x/y)*100 for x, y in zip(map(float, Unaligned), map(int, Total))])
		ind = [x for x, _ in enumerate(Sample)]
		#print(Sample)
		#print(CoV2)
		#print(Human)
		#print(Pathogen)
		#print(Unaligned)
		plt.bar(ind, Unaligned, width=0.5, label='Unaligned', color='blue', bottom=Human + Pathogen+ CoV2)
		plt.bar(ind, Human, width=0.5, label='Human', color='green', bottom=Pathogen + CoV2)
		plt.bar(ind, Pathogen, width=0.5, label='Pathogen', color='yellow', bottom=CoV2)
		plt.bar(ind, CoV2, width=0.5, label='CoV2', color='red')

		stack_plot_file=self.cov2outdir+'Samples_stackbar.png'
		plt.xticks(ind, Sample, rotation = 90)
		plt.ylabel("Relative Composition")
		plt.xlabel("")
		plt.legend(loc="upper right")
		plt.title("Sample Composition")
		plt.savefig(stack_plot_file, bbox_inches = 'tight', pad_inches = 0.1)
		plt.close(fig=None)

		print("Abundance stack plot is generated!")
		#print("Abundance stack plot is generated!")

		image_file_list=[fpkm_plot_file,stack_plot_file]
		return image_file_list

	def get_clade_assessment(self):
		clade_df=pd.DataFrame()
		try:
			#clade_df=self.clade_assessment_obj.get_pandas_df_clade_assessment()
			clade_df=self.varcladeassessmentobj_.get_clade_assessment_data_frame()
		except Exception as e:
			print(e)
			print("Error in Clade Assessment")
			sys.exit(0)
		cov2_clade_out=self.cov2outdir+"Variant_based_clade_assessment.csv"
		clade_df.to_csv(cov2_clade_out)

		return clade_df

	def get_novel_variant(self):
		novel_var_df=pd.DataFrame()
		try:
			#novel_var_df=self.clade_assessment_obj.get_pandas_df_novel_variant()
			novel_var_df=self.varcladeassessmentobj_.get_novel_variant_data_frame()
		except:
			print("Error in Novel Variant")
		novel_var_out=self.cov2outdir+"Novel_variant.csv")
		novel_var_df.to_csv(novel_var_out)

		return novel_var_df


'''
input: output directory of IPD post successful run
'''
import argparse

def main():
	parser = argparse.ArgumentParser(description='Generate an automated report for IPD analysed, SARS-CoV-2 sequenced samples')
	parser.add_argument('-dir',const=None, help="IPD output directory location", dest='inputdir')
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(0)
	args=parser.parse_args()
	outdir=args.inputdir
	if not outdir.endswith("/"):
		outdir = outdir +"/"
	if os.path.exists(outdir):
		obj=CoV2ReportGenerator(outdir)
		out_file=obj.get_basic_stats_tabulation_for_the_batch()
		cov_plot_list=obj.get_coverage_plot()
		ab_plot_list=obj.get_abundance_plot()
		clade_df=obj.get_clade_assessment()
		novel_variant_df=obj.get_novel_variant()
		cov2outdir=outdir+"cov2output/"
		output_file = open(cov2outdir+"Output.html", "w", encoding="utf-8", errors="xmlcharrefreplace")

		header=markdown("#**IPD REPORT**")
		output_file.write(header)

		header1 = markdown("###**Sequence Statistics**")
		output_file.write(header1)

		df=pd.read_table(out_file, sep=',',index_col=0)
		html = df.round(2).to_html()
		output_file.write(html)

		header2 = markdown("###**Coverage Plot**")
		output_file.write(header2)

		if len(cov_plot_list) > 2:
			for i,k in zip(cov_plot_list[0::2], cov_plot_list[1::2]):
				image="![]("+i+") ![]("+k+")"
				html_image = markdown(image)
				output_file.write(html_image)
		elif len(cov_plot_list) == 1:
			image="![]("+cov_plot_list[0]+")"
			html_image = markdown(image)
			output_file.write(html_image)

		#Commented to reason the image placement error
		#image2="![]("+ab_plot_list[0]+")
        #image2=![]("+ab_plot_list[1]+")"
		image2="![]("+ab_plot_list[0]+")"
		image0="![]("+ab_plot_list[1]+")"

		header3 = markdown("###**Relative Abundance**")
		output_file.write(header3)

		html_image2 = markdown(image2)
		output_file.write(html_image2)
		#Added for testing
		html_image0 = markdown(image0)
		output_file.write(html_image0)

		header6 = markdown("###**Novel Variants**")
		output_file.write(header6)

		if not novel_variant_df.empty:
			html2 = novel_variant_df.round(2).to_html(index=False,justify='center')
		else:
			html2 = markdown("No novel variants in the sample.")
		output_file.write(html2)

		header5 = markdown("###**Variant Based Related Strains**")
		output_file.write(header5)

		if not clade_df.empty:
			html3 = clade_df.round(2).to_html(index=False,justify='center')
		else:
			html3 = markdown("No intersecting variants found. Clade cannot be determined.")
		output_file.write(html3)
		output_file.close()
		#HTML to PDF report
		#import pdfkit
		#pdfkit.from_file(os.path.join(cov2outdir,"Output.html"), os.path.join(cov2outdir,"OutputReport.pdf"))
		#2.0 uses wkhtmltopdf
		cmd= "wkhtmltopdf " + cov2outdir+"Output.html " +cov2outdir+"OutputReport.pdf"
		cprocess=subprocess.run(cmd, shell=True)
		cprocess.check_returncode()
if __name__ =="__main__":
		main()
