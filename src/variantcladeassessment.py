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
        #self.site_record_map_={} # site_alt -> [arr of record] # order clade gene site alt
        self.site_clades_map_={}
        self.clades_=[]
        #self.clade_probscore_map_={} # clade -> probability score for a sample
        fi=open(self.filename_)
        if hasheader:
            self.header_ = fi.readline().strip().split('\t')
        for l in fi:
            if len(l)>1:
                cr=CladeTSVRecord(l)
                if not cr.clade_ in self.site_clades_map_:
                    if not cr.clade_ in self.clades_:
                        self.clades_.append(cr.clade_)
                    self.site_clades_map_[cr.key()]=[cr.clade_]
                else:
                    tempclades=self.site_clades_map_[cr.key()]
                    if not cr.clade_ in tempclades:
                        tempclades.append(cr.clade_)
                    self.site_clades_map_[cr.key()]=tempclades
        fi.close()
    def get_clades_for_variant(self,skey):
        return self.site_clades_map_[skey]
    def is_informative_variant(self, pysamvariant):
        skey=str(pysamvariant.pos)+"_"+pysamvariant.alts[0]
        return skey in self.site_clades_map_
    def get_clade_probability_score(self, varfile):
        sample_clade_score_map={}
        for sc in self.clades_:
            sample_clade_score_map[sc]=0
        while(1):
            try:
                vnext=next(varfile)
                if self.is_informative_variant(vnext):
                    skey=str(vnext.pos)+"_"+vnext.alts[0]
                    cladesfound=self.site_clades_map_[skey]
                    numclades=len(cladesfound)
                    for c in self.clades_:
                        if c in cladesfound:
                            sample_clade_score_map[c]=sample_clade_score_map[c]+(1/numclades)
                        else:
                            sample_clade_score_map[c]=sample_clade_score_map[c]+0
            except:
                break
        return sample_clade_score_map
    def get_max_clade_probability(self,vfile):
        maxprobclade=''
        a=self.get_clade_probability_score(vfile)
        max=0
        for c in a:
            if a[c] > max:
                max=a[c]
                maxprobclade=c
        return maxprobclade
class CladeAssignment(object):
    def __init__(self, vcffile):
        vf=pysam.VariantFile(vcffile)
        ctp=CladesTSVParser(GlobalVar.cladestsv_, True)
        self.maxscoreclade_=ctp.get_max_clade_probability(vf)
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

class GisaidCladeTSVRecord(object):
    def __init__(self, rec):
        srec=rec.split("\t")
        self.parentclade_=""
        self.clade_=""
        self.id_=""
        self.name_=""
        self.gisaidnumber_=0
        self.date_=""
        if len(srec) < 2:
            raise GlobalFormatException("GISAID Clade TSV format error for record: %s" %(rec.strip()))
        else:
            self.clade_ = srec[1]
            ssrec=srec[0].split("|")
            if len(ssrec) == 3:
                self.name_=ssrec[0]
                self.id_=ssrec[1]
                self.date_=ssrec[2]
            else:
                raise GlobalFormatException("GISAID Clade TSV format error for record: %s" %(rec.strip()))
            self.gisaidnumber_=int(self.id_.split('_')[2])

class GisaidCladeTSVParser(object):
    def __init__(self):
        gctsv=open(GlobalVar.gisaidcladestsv_)
        self.header_=gctsv.readline().strip()
        self.gisaidnumber_record_map_={} #store GISAID number; ex 736498 for ID EPI_ISL_736498 -> GisaidCladeTSVRecord
        for g in gctsv:
            grec=GisaidCladeTSVRecord(g.strip())
            self.gisaidnumber_record_map_[grec.gisaidnumber_]=grec
    def get_clade(self, ginum):
        gr=self.gisaidnumber_record_map_[ginum]
        return gr.clade_
    def get_parentclade(self, ginum):
        gr=self.gisaidnumber_record_map_[ginum]
        return gr.parentclade_

#A wrapper class to perform both variant check and clade assignment and tabulate in panda frame
class VariantCladeAssessment(object):
    def __init__(self, vcfmap): # , outdir): #{ samplename->vcffile path }
    #self.cov2outdir=os.path.join(self.outdir, "cov2output")
        data=[] #array of array each line. Fill and then give to dataFrame function
        novdfcol=["Sample","Genome","Position","Reference","Altered","Consequence","Gene","Transcript","Protein Change"]
        cladedfcol=["Sample","Number of variants", "Related GISAID genome", "Clade"]
        novvararr=[]
        cladesarr=[]
        self.novelvardf_=pd.DataFrame()
        self.cladedf_=pd.DataFrame()
        #gisaidcladetsvpar=GisaidCladeTSVParser()
        for sample in vcfmap:
            vcffile=vcfmap[sample]
            #print("Processing sample %s" %(sample))
            va=VariantAssessment(vcffile)
            #ca=CladeAssignment(vcffile)
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
            #cclade=ca.assigned_clade()
            #csubclade=ca.assigned_subclade()
            #ccladeprob= '%.2f' % ca.assigned_clade_probability()
            #csubcladeprob= '%.2f' % ca.assigned_subclade_probablility()
            #cladesarr.append([sample, str(cnumvar), crelatedgisaid, cclade+"("+ccladeprob+")" ])
            cclade=""
            ca=CladeAssignment(vcffile)
            if ca.maxscoreclade_ != '':
                cclade=ca.maxscoreclade_
            '''
            if va.closest_gisaid_sample() != "":
                cg=va.closest_gisaid_sample().split('-')
                cginum=int(cg[2])
                cclade=gisaidcladetsvpar.get_clade(cginum)
            '''
            cladesarr.append([sample, str(cnumvar), crelatedgisaid, cclade ])
            print("######## Sample : %s ; Related GISAID sample: %s ; Clade : %s" %(sample, crelatedgisaid, cclade))
        if len(novvararr)>0:
            self.novelvardf_ = pd.DataFrame(novvararr, columns = novdfcol)
        else:
            self.novelvardf_= pd.DataFrame([ ["","","","","","","","",""] ], columns= novdfcol)
        self.cladedf_ = pd.DataFrame(cladesarr, columns = cladedfcol)
        #xxxxxxxxxx
    def get_novel_variant_data_frame(self):
        return self.novelvardf_
    def get_clade_assessment_data_frame(self):
        return self.cladedf_
