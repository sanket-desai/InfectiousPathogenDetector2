#!/usr/bin/python3
'''
Author:         Sanket Desai
Date:           08/01/2021
Version:        1.0
IPD Version:    2.0
Description:    Given a VCF file these set of classes will return basic Assessment
                number of variants per file, novel variant list and annotation, closest isolate

Read the clade and subclade files in the cov2reference directory


'''
from globals import *
import pysam
import sys
import os
import pandas as pd

class TabvarRecord(object):
    def __init__(self, l):
        si=l.strip().split("\t")
        if len(si) != 17:
            print(si)
            print("Format issue! Please check the tabvar file")
            sys.exit(0)
        self.chr_=si[0]
        self.pos_=si[1]
        self.type_=si[2]
        self.ref_=si[3]
        self.alt_=si[4]
        self.evidence_=si[5]
        self.feature=si[6]
        self.strand_=si[7]
        self.nt_pos_=si[8]
        self.aa_pos_=si[9]
        self.locus_=si[10]
        self.gene_=si[11]
        self.product_=si[12]
        self.effect_=si[13]
        self.samplename_=si[14]
        self.gisaid_id_=si[15]
        self.date_=si[16]
    def key(self):
        return self.ref_+"_"+self.pos_+"_"+self.alt_

class CladeTSVRecord(object):
    def __init__(self, rec):
        sl=rec.strip().split('\t')
        try:
            self.clade_=sl[0]
            self.gene_=sl[1]
            self.site_=int(sl[2])
            self.alt_=sl[3]
        except:
            raise GlobalFormatException("Clade TSV file is in incorrect format! Please check.")
    def record_as_list(self):
        return [self.clade_, self.gene_, self.site_, self.alt_]
    def key(self):
        return str(self.site_)+"_"+self.alt_

class CladesTSVParser(object):
    def __init__(self, fname, hasheader):
        self.header_=''
        self.filename_=fname
        self.site_record_map_={} # site_alt -> [arr of record] # order clade gene site alt
        self.clade_sitealt_map_={} # clade -> [arr of site_alt keys]
        fi=open(self.filename_)
        if hasheader:
            self.header_ = fi.readline().strip().split('\t')
        for l in fi:
            if len(l)>1:
                cr=CladeTSVRecord(l)
                self.site_record_map_[cr.key()]=cr
                if not cr.clade_ in self.clade_sitealt_map_:
                    self.clade_sitealt_map_[cr.clade_]=[ cr.key() ]
                else:
                    csa=self.clade_sitealt_map_[cr.clade_]
                    csa.append(cr.key())
                    self.clade_sitealt_map_[cr.clade_]=csa
        fi.close()
    def get_keys_for_clade(self,sclade):
        return self.clade_sitealt_map_[sclade]
    def get_number_of_clade_defining_variants(self, sclade):
        return len(self.clade_sitealt_map_[sclade])
    def get_number_of_subclade_defining_variants(self, sclade):
        return len(self.subclade_sitealt_map_[sclade])
    def get_clades_matching_variant(self, pysamvariant):
        cla=[] #multiple clades could be mappng same variant
        for va in pysamvariant.alts:
            vk=str(pysamvariant.pos)+"_"+va
            if vk in self.site_record_map_:
                cla.append(self.site_record_map_[vk].clade_)
        return cla

class CladeAssignment(object):
    def __init__(self, vcffile):
        self.clade_score_map_={}
        self.subclade_score_map_={}
        #stores the clade -> number of clade supporting variants found in vcf
        #finally the one with highest support will be returned as an assigned clade
        self.cladetsvparser_obj_=CladesTSVParser(GlobalVar.cladestsv_, True)
        #primary clade assignment
        #ctsvp=CladesTSVParser(GlobalVar.cladestsv_, True)
        try:
            vf=pysam.VariantFile(vcffile)
            while(1):
                try:
                    vnext=next(vf)
                    arrclades=self.cladetsvparser_obj_.get_clades_matching_variant(vnext)
                    for c in arrclades:
                        if not c in self.clade_score_map_:
                            self.clade_score_map_[c]=1
                        else:
                            self.clade_score_map_[c]=self.clade_score_map_[c]+1
                except:
                    break
        except Exception as e:
            print("Clade assignment failed due to followig exception: ")
            print(e)
            sys.exit(0)
        self.subcladetsvparser_obj_=CladesTSVParser(GlobalVar.subcladestsv_, False)
        #sctsvp=CladesTSVParser(GlobalVar.subcladestsv_, False)
        try:
            vf=pysam.VariantFile(vcffile)
            while(1):
                try:
                    vnext=next(vf)
                    arrclades=self.subcladetsvparser_obj_.get_clades_matching_variant(vnext)
                    for c in arrclades:
                        if not c in self.subclade_score_map_:
                            self.subclade_score_map_[c]=1
                        else:
                            self.subclade_score_map_[c]=self.subclade_score_map_[c]+1
                except:
                    break
        except Exception as e:
            print("Sub-clade assignment failed due to followig exception: ")
            print(e)
            sys.exit(0)
    def assigned_clade(self):
        maxscoreclade=''
        maxscore=0
        for c in self.clade_score_map_:
            if self.clade_score_map_[c] > maxscore:
                maxscoreclade=c
                maxscore=self.clade_score_map_[c]
        #All clade defining variants match
        if not maxscoreclade == '':
            if not maxscore >= self.cladetsvparser_obj_.get_number_of_clade_defining_variants(maxscoreclade):
                maxscoreclade=''
        return maxscoreclade
    def assigned_subclade(self):
        maxscoreclade=''
        maxscore=0
        for c in self.subclade_score_map_:
            if self.subclade_score_map_[c] > maxscore:
                maxscoreclade=c
                maxscore=self.subclade_score_map_[c]
        #All clade defining variants match
        if not maxscoreclade == '':
            if not maxscore >= self.cladetsvparser_obj_.get_number_of_clade_defining_variants(maxscoreclade):
                maxscoreclade=''
        return maxscoreclade

class VariantAssessment(object):
    def __init__(self, vcffile):
        self.number_of_variants_=0
        self.novel_variant_list_=[] # list of pysam variant objects which are not found in the database
        pysamvcf=pysam.VariantFile(vcffile)
        self.closest_gisaid_sample_map_={} # gisaid / EPI -> number of matching variants
        tabfile=pysam.TabixFile(GlobalVar.tabvardatabase_)
        while(1):
            try:
                isnovel=True
                vrec=next(pysamvcf)
                chrom=vrec.chrom
                if chrom.startswith("NC_045512"):
                    if chrom.find('.')>0:
                        chrom=chrom[:chrom.find('.')]
                    self.number_of_variants_+=1
                    query=chrom+":"+str(vrec.pos)+"-"+str(vrec.pos)
                    matchinggisaid=[]
                    for trec in tabfile.fetch(query):
                        rtabvar=TabvarRecord(trec)
                        #check if variants are same
                        dkey=rtabvar.key()
                        for al in vrec.alts:
                            vkey=vrec.ref+"_"+str(vrec.pos)+"_"+al
                            if vkey == dkey:
                                isnovel=False
                                if not rtabvar.gisaid_id_ in matchinggisaid:
                                    matchinggisaid.append(rtabvar.gisaid_id_)
                else:
                    print("Ignored variant %s" %(chrom))
                for m in matchinggisaid:
                    if not m in self.closest_gisaid_sample_map_:
                        self.closest_gisaid_sample_map_[m]=1
                    else:
                        self.closest_gisaid_sample_map_[m]=self.closest_gisaid_sample_map_[m]+1
                if isnovel:
                    if vrec.chrom.startswith('NC_045512'):
                        self.novel_variant_list_.append(vrec)
            except Exception as se:
                print(se)
                break
        self.closest_gisaid_number_of_variant_overlap_=0
        self.closest_gisaid_id_=''
        for c in self.closest_gisaid_sample_map_:
            if self.closest_gisaid_sample_map_[c] > self.closest_gisaid_number_of_variant_overlap_:
                self.closest_gisaid_number_of_variant_overlap_=self.closest_gisaid_sample_map_[c]
                self.closest_gisaid_id_=str(c)
    def novel_variant_list(self):
        return self.novel_variant_list_
    def closest_gisaid_sample(self):
        return self.closest_gisaid_id_
    def closest_gisaid_number_of_variant_overlap(self):
        return self.closest_gisaid_number_of_variant_overlap_
    def get_number_of_variants(self):
        return self.number_of_variants_

#A wrapper class to perform both variant check and clade assignment and tabulate in panda frame
class VariantCladeAssessment(object):
    def __init__(self, vcfmap): # , outdir): #{ samplename->vcffile path }
    #self.cov2outdir=os.path.join(self.outdir, "cov2output")
        data=[] #array of array each line. Fill and then give to dataFrame function
        novdfcol=["Sample","Genome","Position","Reference","Altered","Consequence","Gene","Transcript","Protein Change"]
        cladedfcol=["Sample","Number of variants", "Related GISAID genome", "Clade", "Subclade"]
        novvararr=[]
        cladesarr=[]
        self.novelvardf_=pd.DataFrame()
        self.cladedf_=pd.DataFrame()
        for sample in vcfmap:
            vcffile=vcfmap[sample]
            #print("Processing sample %s" %(sample))
            va=VariantAssessment(vcffile)
            ca=CladeAssignment(vcffile)
            if len(va.novel_variant_list_) > 0:
                for n in va.novel_variant_list_:
                    ann=n.info['ANN'][0].split('|')
                    alt=ann[0]
                    conseq=ann[1]
                    gene=ann[3]
                    transcript=ann[4]
                    protein_change=ann[10]
                    nvrec=[ sample, n.chrom, str(n.pos), n.ref, alt, conseq, gene, transcript, protein_change]
                    novvararr.append(nvrec)
            else:
                print("No novel variants found!")
            #start filling clade line
            cnumvar=va.get_number_of_variants()
            print("Number of variants %d " %(cnumvar))
            crelatedgisaid=va.closest_gisaid_sample()+ " (" + str(va.closest_gisaid_number_of_variant_overlap()) +")"
            cclade=ca.assigned_clade()
            csubclade=ca.assigned_subclade()
            cladesarr.append([sample, str(cnumvar), crelatedgisaid, cclade, csubclade])
        if len(novvararr)>0:
            self.novelvardf_ = pd.DataFrame(novvararr, columns = novdfcol)
        else:
            self.novelvardf_= pd.DataFrame([ [str(sample),"","","","","","","",""] ], columns= novdfcol)
        self.cladedf_ = pd.DataFrame(cladesarr, columns = cladedfcol)
        #xxxxxxxxxx
    def get_novel_variant_data_frame(self):
        return self.novelvardf_
    def get_clade_assessment_data_frame(self):
        return self.cladedf_
		#cov2_clade_out=os.path.join(self.cov2outdir,"Variant_based_clade_assessment.csv")
		#clade_df.to_csv(cov2_clade_out)
