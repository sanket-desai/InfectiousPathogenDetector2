#!/usr/bin/python3
'''

Author	    : Sanket Desai / Aishwarya Rane
Date		: 10/03/2021
Description : Automatic updater for the SARS-CoV-2 module variant database

'''
import argparse
import sys
import os
import subprocess
import pysam
#from globals import *
import ntpath
import multiprocessing
import datetime
from commandlineparser import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pathlib
import glob
import pandas as pd

class TabvarRecord(object):
    def __init__(self, l):
        si=l.strip().split(",")
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
    def mutation_key(self):
        return self.ref_+"_"+self.pos_+"_"+self.alt_
def is_equal(l1, l2): #compare if elements of two arrays are same, irrespective of order
    #algorithm: 1. check lengths are same, if yes 2 sort both the arrays ascending or descending and use ==
    if len(l1) == len(l2):
        l1.sort()
        l2.sort()
        if l1 == l2:
            return True
        else:
            return False
    else:
        return False

def file_choices(choices,fname):
    ext = ''.join(pathlib.Path(fname).suffixes)
    fnamewrong=True
    for c in choices:
        if fname.endswith(c):
            fnamewrong=False
    if ext.lower() not in choices and fnamewrong:
        print(ext)
        raise argparse.ArgumentTypeError('File must have proper extension')
    return fname

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise NotADirectoryError(path)
#Read original fasta file, filter sequences (trim N, length filter), create metadata file, split fasta into individual files and save
class GISAIDPreprocessing(object):
    def __init__(self, gisaidfasta, minlength, maxN, outdir):
        self.fastadoutdir_=outdir+"splitfasta/"
        if not os.path.exists(self.fastadoutdir_):
            os.mkdir(self.fastadoutdir_)
        self.gisaid_fasta_=gisaidfasta
        self.metadatafile_=outdir+"metadata.tsv"
        ofmetadata=open(self.metadatafile_, 'w')
        for rec in SeqIO.parse(self.gisaid_fasta_,"fasta"):
            cleanedseq=self.sequence_cleanup(rec)
            if len(cleanedseq) > minlength and cleanedseq.seq.count("N") <= maxN : #if sequence pass two criteria
                #change header
                #originalheader=cleanedseq.id #saved for metadata
                cleanedseq.id=cleanedseq.id.replace("/","_")
                cleanedseq.id=cleanedseq.id.replace("|","_")
                cleanedseq.id=cleanedseq.id.replace("-","_")
                f=cleanedseq.id
                estart=f.find("EPI_ISL_")
                if estart != -1: #if EPI ID not found in header, ignore the sequence
                    seqn=f[:estart-1]
                    temp=f[estart:]
                    stemp=temp.split("_")
                    try:
                        eslid=stemp[0]+"_"+stemp[1]+"_"+stemp[2]
                        date="-".join(stemp[3:])
                        ofmetadata.write("%s\t%s\t%s\t%s\n" %(cleanedseq.id, seqn, eslid, date) )
                        ofasta=open(self.fastadoutdir_+cleanedseq.id+".fa", 'w')
                        s=SeqIO.write(cleanedseq, ofasta, 'fasta')
                        ofasta.close()
                    except Exception as e:
                        pass #unknown exception
        ofmetadata.close()

    def sequence_cleanup(self, r): #input -> fasta seqIO record
        nr=r
        nsq=r.seq
        clean=True
        while nsq.endswith("N"):
            nsq=nsq[:-1]
            clean=False
        while nsq.endswith("-"):
            nsq=nsq[:-1]
            clean=False
        while nsq.endswith("N"):
            nsq=nsq[:-1]
            clean=False
        while nsq.startswith("-"):
            nsq=nsq[1:]
            clean=False
        while nsq.startswith("N"):
            nsq=nsq[1:]
            clean=False
        while nsq.startswith("-"):
            nsq=nsq[1:]
            clean=False
        nr.seq=nsq
        nr.description=""
        return nr

#define parser and subparsers
def main():
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-f',dest='ifasta',type=str, required=True,help="GISAID fasta sequence file")
    parent_parser.add_argument('-l',dest='mingenlen',type=int,default=29000, required=False,help="Minimum SARS-CoV-2 genome length to consider for database creation")
    parent_parser.add_argument('-t',dest='threads',type=int,default=2,required=False,help="Set number of threads (default = 4)")
    parent_parser.add_argument('-N',dest='maxnumn',type=int,default=10000,required=False,help="Maximum number of 'N' allowed in the genome sequence")
    parent_parser.add_argument('-o',dest='outdir',type= lambda path:dir_path(path),required=True,default=os.getcwd(),help="Output directory to process the intermediate files")
    #parser = argparse.ArgumentParser(add_help=False)
    sarscov2gbk="../data/cov2moduleref/NC_045512.gb"
    repfileheader="chrom\tpos\ttype\tref\talt\tevidence\tfeattype\tstrand\tnt_pos\taa_poseffect\tlocus_tag\tgene\tproduct\tsamplename\tid\tdate"
    args = parent_parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    else:
        if not args.outdir.endswith("/"):
            args.outdir=args.outdir+"/"
        metadatafile=args.outdir+"metadata.tsv"
        splitfastadir=args.outdir+"splitfasta/"
        snippyoutdir=args.outdir+"snippyoutdir/"
        combinedsnippyfname=args.outdir+"combinedsnippyout.tab"
        samplemutationprofilefname=args.outdir+"sammutprofile.tab"
        repmutationprofilefname=args.outdir+"representative_mutation_profiles.tab"
        repsamplelistfname=args.outdir+"representative_samplelist.tab"
        ofcombinedsnippyout=open(combinedsnippyfname,'w')
        ofcombinedsnippyout.write("%s\n" %(repfileheader))
        gp=GISAIDPreprocessing(args.ifasta, args.mingenlen, args.maxnumn, args.outdir)
        #Run snippy for each fasta
        
	    for f in glob.glob(splitfastadir+"*.fa"):
            pre=os.path.splitext(os.path.basename(f))[0]
            cmd="snippy --ctgs " +f+ " --ref " +sarscov2gbk+ " --outdir " +snippyoutdir+ " --prefix " +pre+ " --force --cpus "+str(args.threads)
            cprocess=subprocess.run(cmd, shell=True)
            cprocess.check_returncode()
        #combine variants from tab files
        ind=1
        uniqrepnames=[]
        for f in glob.glob(snippyoutdir+"*.tab"):
            print("%d %s" %(ind, f))
            ind=ind+1
            estart=f.find("EPI_ISL_")
            fi=open(f)
            if estart != -1:
                seqn=f[:estart-1].replace(".tab","")
                temp=f[estart:]
                stemp=temp.split("_")
                try:
                    eslid=stemp[0]+"_"+stemp[1]+"_"+stemp[2]
                    date="-".join(stemp[3:])
                    header=fi.readline()
                    for i in fi:
                        si=i.strip().split("\t")
                        if len(si) < 14:
                            while len(si) < 14:
                                si.append("")
                            #print(si)
                            #print("Exiting inside first if")
                            #sys.exit(0)
                        #ofo.write("%s\t%s\t%s\t%s\n" %(i.strip(), seqn, eslid, date.replace(".tab","")) )
                        #ii=",".join( i.strip().split("\t") )
                        ii=",".join( si )
                        if len(seqn)==0 or len(eslid) == 0 or len(ii) == 0:
                            print(si)
                            sys.exit(0)
                        sseqn=seqn.split("/")
                        seqn=sseqn[len(sseqn)-1]    
                        ofcombinedsnippyout.write("%s,%s,%s,%s\n" %(ii, seqn, eslid, date.replace(".tab","")))
                except Exception as e:
                    print(e)
                    print("Error in the file parsing!")
            fi.close()
        ofcombinedsnippyout.close()
        print("Combined tabvar created")
        fi2=open(combinedsnippyfname)
        sammutprofilefo=open(samplemutationprofilefname,'w')
        sample_mutation_map={} # samples -> mutations
        header=fi2.readline()
        for i in fi2:
            tb=TabvarRecord(i)
            if not tb.samplename_ in sample_mutation_map:
                sample_mutation_map[tb.samplename_]=[tb.mutation_key()]
            else:
                arr=sample_mutation_map[tb.samplename_]
                arr.append(tb.mutation_key())
                sample_mutation_map[tb.samplename_]=arr
        fi2.close()
        numvar_samples_map={} # number of variants -> [ sample names ]
        for i in sample_mutation_map:
            a=sample_mutation_map[i]
            #a.sort()
            #sammutprofilefo.write("%s\t%s\n" %(i,",".join(a)))
            if len(a) in numvar_samples_map:
                arr2=numvar_samples_map[len(a)]
                arr2.append(i)
                numvar_samples_map[len(a)]=arr2
            else:
                numvar_samples_map[len(a)]=[i]
        #print("%s\t%d" %(i,len(a)))
        repmutfo=open(repmutationprofilefname,'w')
        repsamfo=open(repsamplelistfname,'w')
        for n in numvar_samples_map:
            rep_mutation_map={}
            repsam_samples_map={}
            samples=numvar_samples_map[n]
            print("variantcount:%d;samplecount:%d" %(n, len(samples)))
            tobeskipped=[]
            j=1
            for s in samples:
                if not s in tobeskipped:
                    b=sample_mutation_map[s]
                    b.sort()
                    rep_mutation_map[s]=b
                    repsam_samples_map[s]=[]
                    nkeys=samples[j:]
                    for nn in nkeys:
                        if not nn in tobeskipped:
                            if is_equal(sample_mutation_map[s], sample_mutation_map[nn]):
                                tobeskipped.append(nn)
                                arr3=repsam_samples_map[s]
                                arr3.append(nn)
                                repsam_samples_map[s]=arr3
                ss=s.split("/")
                s=ss[len(ss)-1]
                sammutprofilefo.write("%s\t%s\n" %(s,",".join(b)))
                j+=1
            for m in rep_mutation_map:
                ss=m.split("/")
                mm=ss[len(ss)-1]
                repmutfo.write("%s\t%s\n" %(mm,",".join(rep_mutation_map[m])))
            for o in rep_mutation_map:
                ss=o.split("/")
                oo=ss[len(ss)-1]
                if not oo in uniqrepnames:
                    uniqrepnames.append(oo)
                repsamfo.write("%s\t%s\n" %(oo,",".join(repsam_samples_map[o])))
        repmutfo.close()
        repsamfo.close()
        sammutprofilefo.close()
        fi2=open(combinedsnippyfname)
        header=fi2.readline().strip()
        ###take keys from repfile, read from combined snippyfname and create pandas, sort and write back into file
        outdf=pd.DataFrame()
        print("Extracting representative variants!")
        for l in fi2:
            sl=l.strip().split(",")
            if sl[14] in uniqrepnames:
                outdf=outdf.append([sl])
        print("Writing representative tabvar to file!")
        #outdf=outdf.sort_values(by=1, ascending=True)
        reptabvar=args.outdir + "representative_tabvar.tsv"
        outdf.to_csv(reptabvar, sep="\t", index=False, header=False)
        fi2.close()
        print("Sorted representative database created!")
        cmd="sort -k2,2n "+reptabvar +" | bgzip > "+args.outdir+"representative_tabvar_sorted.tsv.gz"
        cprocess=subprocess.run(cmd, shell=True)
        cprocess.check_returncode()
        cmd="tabix -b 2 -e 2 -f "+args.outdir+"representative_tabvar_sorted.tsv.gz"
        cprocess=subprocess.run(cmd, shell=True)
        cprocess.check_returncode()
if __name__ == "__main__":
    main()
