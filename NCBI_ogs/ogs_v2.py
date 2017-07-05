import pandas as pd
import re
import os
import subprocess
from sys import argv
import logging

'''
usage: python ogs_v1.py biosample_ID file.gff3
example: python ogs_v1.py GCA_000696205.1 clec_OGS_v1_2_with_pep_CDS.gff3
'''
Biosample=argv[1]
file_in=argv[2]
#file_in='clec_OGS_v1_2_with_pep_CDS.gff3'
print Biosample
print file_in

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S =%p',
                    filename='ogs.log',
                    filemode='w',
                    level=logging.INFO)
logging.info('start program')
#gff_data = pd.read_csv('clec_OGS_v1_2_with_pep_CDS.gff3', sep="\t", header = None,comment='#')
'''
remove NOTES attribute
name attributes: gene, mRNA, pseudogenic_transcript

'''

f_out=open('file_out.txt', 'w')
#1. remove %09
#preprocess= subprocess.Popen("sed s/%09//g "+ file_in +">temp.gff", stdout=subprocess.PIPE,shell=True)
#gff_data = pd.read_csv('temp.gff', sep="\t", header = None,comment='#')
gff_data = pd.read_csv('clec_OGS_v1_2_with_pep_CDS.2.gff3', sep="\t", header = None,comment='#')
#2. remodel pseudogene
def Remodel_pseudogenes():
    pseudo_list=[]


    df0 = gff_data.where(gff_data[2].str.contains("pseudogene"))
    idx0 = df0.dropna().index.tolist()
    for index in idx0:
        gff_data.iloc[index, 2] = "gene"
       # gff_data.iloc[index, 8] = gff_data.loc[index, 8] + ";pseudogene=unknown"
        # get gene_ID
        gene_ID = "".join(re.match("^ID=([^;]+);.+$", gff_data.iloc[index, 8]).groups())
        # print gene_ID
        pseudo_list.append(gene_ID)

    #pseudo_list=['CLEC026002', 'CLEC026001','CLEC0260214']

    #for each child of pseudogene, alter col9
    for pi in range(len(pseudo_list)):
        df1 = gff_data.where(gff_data[8].str.contains(pseudo_list[pi]))
        idx2 = df1.dropna().index.tolist()
        for i in idx2:
            gff_data.iloc[i, 8] = str(gff_data.loc[i, 8]).rstrip(';') + ";pseudogene=unknown"
        ## check it's child col3 whether has multi-name , then throw out warning !
           # if (gff_data.loc[i,2]!="gene") and (gff_data.loc[i,2]!="pseudogenic_transcript" or gff_data.loc[i,2]!= "pseudogenic_exon"):
            if (gff_data.loc[i, 2] != "gene") and not ( re.search( "pseudogenic_transcript",gff_data.loc[i, 2]) or re.search( "pseudogenic_exon",gff_data.loc[i, 2])):
                child_ID = "".join(re.match("^ID=([^;]+);.+$", gff_data.iloc[i, 8]).groups())
                feature=gff_data.loc[i,2]
                logging.warning(" multiple types of child features : "+"ID= "+child_ID+";  feature= "+feature)

   # gff_data.to_csv('clec_OGS_v1_2.step2.gff3', header=False, index=False, sep='\t')
#3-1 util check for locus tag in col9
def check_locus():
    df0 = gff_data.where(gff_data[2].str.contains("gene"))
    idx0 = df0.dropna().index.tolist()
    df1 = df0.where(df0[8].str.contains("locus_tag"))
    idx1 = df1.dropna().index.tolist()
    if len(idx1) == len(idx0):
        f_out.write("All gene features contain a locus_tag attribute")
    elif 0 < len(idx1) < len(idx0):
        f_out.write("Some gene features contain a locus_tag attribute")
    elif len(idx1) == 0:
        f_out.write("No gene features contain a locus_tag attribute")
#3-2. util get locus tag
def locus_tag_util(biosample):
    import subprocess
    import logging
    process = subprocess.Popen("esearch -db assembly -query " + biosample ,stdout=subprocess.PIPE, shell=True)
    process2= subprocess.Popen("efetch -format docsum" ,stdin=process.stdout , stdout=subprocess.PIPE,shell=True)
    process3=subprocess.Popen("xtract -pattern DocumentSummary -element BioprojectId", stdin=process2.stdout,stdout=subprocess.PIPE, shell=True)
    proc_stdout_bioproject = process3.communicate()[0].strip()

    process31=subprocess.Popen("efetch -db bioproject -id "+proc_stdout_bioproject+" -format xml",stdout=subprocess.PIPE, shell=True)
    #process32=subprocess.Popen("xtract -pattern DocumentSummary -element LocusTagPrefix", stdin=process31.stdout,stdout=subprocess.PIPE, shell=True)

    process32=subprocess.Popen("grep LocusTagPrefix -", stdin=process31.stdout,stdout=subprocess.PIPE, shell=True)

    proc_stdout_locus = process32.communicate()[0].strip()
    ##print "bioproject\n"+ proc_stdout_bioproject
    ##print "locus\n"+proc_stdout_locus


    biosample_1 = subprocess.Popen("esearch -db assembly -query GCA_000696205.1" ,stdout=subprocess.PIPE, shell=True)
    biosample_2= subprocess.Popen("efetch -format docsum" ,stdin=biosample_1.stdout , stdout=subprocess.PIPE,shell=True)
    #biosample_3=subprocess.Popen("xtract -pattern DocumentSummary -element BiosampleID", stdin=biosample_2.stdout,stdout=subprocess.PIPE, shell=True)
    biosample_3=subprocess.Popen("xtract -pattern DocumentSummary -element BioSampleAccn", stdin=biosample_2.stdout,stdout=subprocess.PIPE, shell=True)
    proc_stdout_biosample = biosample_3.communicate()[0].strip()

    #print "biosample: "+ proc_stdout_biosample #SAMN02645558

    # efetch -db bioproject -id 229125 -format xml | grep "LocusTagPrefix"

    import re
    locus_list=proc_stdout_locus.lstrip("\t").split("\n")

    locus_dict={}
    for i in range(len(locus_list)):
        biosample_id="".join(re.match(".*<LocusTagPrefix biosample_id=\"(\w+)\".+>", locus_list[i]).groups())
        locustag = "".join(re.match(".*<LocusTagPrefix biosample_id=\"\w+\">\s*([\w]+).+>", locus_list[i]).groups())
    #    print biosample_id,locustag
        locus_dict.update({biosample_id:locustag})
   # print locus_dict
    try:
        LOCUS_TAG=locus_dict.get(proc_stdout_biosample)
    #    print LOCUS_TAG
        return LOCUS_TAG
    except:
        logging.warning("biosample ID not in locus prefix list")
#print locus_tag_util()
locus_tag= locus_tag_util(Biosample)
#locus_tag='W904'
#4. ID attribute  mRNAS , CDS transcript_id, protein_id
def ID_attribute():
    df1 = gff_data.where(gff_data[2].str.contains("mRNA"))
    idx1 = df1.dropna().index.tolist()
    for i in idx1:
        try:
            ID="".join(re.match("^.*Parent=([^;-]+).*$", gff_data.iloc[i, 8]).groups())
            gff_data.iloc[i, 8] = str(gff_data.loc[i, 8]).rstrip(";") + ";transcript_id=gnl|"+locus_tag+"|"+ID+"-RA"
            #print gff_data.iloc[i,8]
        except:
            print "mRNA_ID ERROR"
            print  gff_data.iloc[i,8]
def cds_id_attribute():
#need to throw a warning if there is no polypetide available
#CDS,
    df2 = gff_data.where(gff_data[2].str.contains("CDS"))
    idx2 = df2.dropna().index.tolist()
    for i in idx2:
            try:
                ID = "".join(re.match("^.*Parent=([^;\-]+).*$", gff_data.iloc[i, 8]).groups())
                #print ID
                ## search for polypeptide first
                pep=ID+"-PA"
                test=gff_data.where(gff_data[8].str.contains("ID="+pep))
                if(   len(test.dropna().index.tolist())<1):
                    logging.warning("CDS without polypeptide attribute!"+ID)

                gff_data.iloc[i, 8] = str(gff_data.loc[i, 8]).rstrip(";") + ";protein_id=gnl|" + locus_tag + "|" + ID+"-PA"
               # print gff_data.iloc[i,8]
            except:
                print "CDS_ID ERROR"
                print gff_data.iloc[i,8]
#5. remove note attribute
def remove_Note():
        df1 = gff_data.where(gff_data[8].str.contains("Note="))
        idx2 = df1.dropna().index.tolist()
        ## check if the Note is the last term!
        for i in idx2:
                Note_term = "".join(re.match("^.*Note=([^;]+).+$", gff_data.iloc[i, 8]).groups())
              #  print Note_term
                Note_new= re.sub(Note_term, "Additional notes are available on the i5k Workspace gene page.", gff_data.iloc[i,8])
              #  print Note_new
                gff_data.iloc[i, 8] = Note_new
              #  print gff_data.iloc[i,8]

#6. Name attributes(gene)
def name_Attributes_gene():
        df1 = gff_data.where(gff_data[2].str.contains("gene"))
        idx1 = df1.dropna().index.tolist()
        #print len(idx1) # 14087
        df2 = df1.where(df1[8].str.contains("Name="))
        idx2 = df2.dropna().index.tolist()
        #print idx2,len(idx2) #696
        for i in idx2:
                try:
                        name_attribute = "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups()).rstrip("\s")
                        temp=name_attribute+";gene="+name_attribute
                        name_new = re.sub(name_attribute, (temp), gff_data.iloc[i, 8])
                        gff_data.iloc[i, 8] = name_new
                except:
                        logging.warning( "ERROR gene name attribute:"+gff_data.iloc[i,8])
                        name_attribute= "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups()).rstrip("\s")
                        new = name_attribute+";gene="+name_attribute
                        #a= re.sub(name_attribute, new, gff_data.iloc[i, 8])
                        print new

#7. Name_attribute for mRNA, pseudo_transcript
# multi conditions for features:
def name_Attribute_mRNA_pseudo():
        feature_type = ['mRNA', 'pseudogenic_transcript']  #total: 14214, only 2 pseudogenic_transcript
        df1 = gff_data.where(gff_data[2].str.contains('|'.join(feature_type)))
        idx1 = df1.dropna().index.tolist()
        df2 = df1.where(df1[8].str.contains("Name="))
        idx2 = df2.dropna().index.tolist()
        #print df2, len(idx2) # 1481
        for i in idx2:
                 try:
                        name_attribute = "".join(re.match("^.*Name=([^;]+).+$", gff_data.iloc[i, 8]).groups())
                        temp=name_attribute+";product="+name_attribute
                        name_new = re.sub(name_attribute, (temp), gff_data.iloc[i, 8])
                        #print name_new
                        gff_data.iloc[i, 8] = name_new
                 except:
                        logging.warning( "ERROR mRNA, pseudogenic_transcript name attribute:"+gff_data.iloc[i,8])
                        print "".join(re.match("^.*Name=([^;]+).+$", gff_data.iloc[i, 8]).groups())

#8. name attribute (tRNAs, rRNAs and ncRNAs)
def name_Attribute_trncRNA():
    feature_type = ['tRNA', 'rRNA','ncRNA']  #total: 14214, only 2 pseudogenic_transcript
    df1 = gff_data.where(gff_data[2].str.contains('|'.join(feature_type)))
    idx1 = df1.dropna().index.tolist()
    #print df1 #3
    df2 = df1.where(df1[8].str.contains("Name="))
    idx2 = df2.dropna().index.tolist()
    for i in idx2:
        try:
            name_attribute = "".join(re.match("^.*Name=([^;]+).+$", gff_data.iloc[i, 8]).groups())
            name_new = re.sub("Name=", "product=" , gff_data.iloc[i, 8])
           # print name_new
            gff_data.iloc[i, 8] = name_new
        except:

            logging.warning("ERROR rRNA,tRNA,ncRNA  name attribute:" + gff_data.iloc[i, 8])
            print "".join(re.match("^.+Name=([^;]+).+$", gff_data.iloc[i, 8]).groups())

#9. stop_codon_read_through, stop_codon_readthrough, deletion, insertion, substitution
def Scan_for_genome_alterations():
    dict={}
    alteration_list=["stop_codon_read_through", "stop_codon_readthrough", "deletion", "insertion", "substitution"]
    for i in range(len(alteration_list)):
        df1 = gff_data.where(gff_data[2].str.contains(alteration_list[i]))
        dx1 = df1.dropna().index.tolist()
        dict.update({alteration_list[i]:len(dx1)})
  #  print "Program Summary:"
    f_out.write("\nProgram Summary:")
    sum= " ".join([str(v)+" of type "+k+","  for k,v in dict.items()])
   # print "Found the following features in file [input gff3 file]: "+sum.rstrip(",")+"."
    f_out.write( "\nFound the following features in file [input gff3 file]: "+sum.rstrip(",")+".")
    #print "You may need to manually add a transl_except attribute to the corresponding CDS feature lines. See https://www.ncbi.nlm.nih.gov/sites/genbank/genomes_gff/ for details."
    f_out.write("\nYou may need to manually add a transl_except attribute to the corresponding CDS feature lines. See https://www.ncbi.nlm.nih.gov/sites/genbank/genomes_gff/ for details.")

#10. Description attributes
def description_attribute():
    df1 = gff_data.where(gff_data[8].str.contains("description="))
    dx1 = df1.dropna().index.tolist() #1254
    df2= gff_data.where(gff_data[2].str.contains("gene"))
    df21=df2.where(df2[8].str.contains("description="))
    dx2 = df21.dropna().index.tolist() #425
    #print len(dx1), len(dx2)
    #print str(len(dx1)-len(dx2))+"  description attributes found outside of gene features."
    f_out.write("\n"+str(len(dx1)-len(dx2))+"  description attributes found outside of gene features.")

#11. Symbol
#Program will scan for symbol attributes in features other than gene features.
def symbol_attribute():
    df1 = gff_data.where(gff_data[8].str.contains("symbol="))
    dx1 = df1.dropna().index.tolist() #1254

    df2= gff_data.where(gff_data[2].str.contains("gene"))
    df21=df2.where(df2[8].str.contains("description="))
    dx2 = df21.dropna().index.tolist() #425
    #print len(dx1), len(dx2)
    print str(len(dx1)-len(dx2))+"  symbol attributes found outside of gene features."
    f_out.write("\n"+str(len(dx1)-len(dx2))+"  symbol attributes found outside of gene features.")

#12. Dbxrefs
## first for gene feature, copy ID to Dbxref
def Dbxref():
    df1=gff_data.where(gff_data[2].str.contains("gene"))
    dx1 = df1.dropna().index.tolist() #1254
    for i in dx1:
        try:
            gene_ID = "".join(re.match("^ID=([^;]+).*$", gff_data.iloc[i, 8]).groups())
            print gene_ID

            # print name_new
            gff_data.iloc[i, 8] =  gff_data.iloc[i, 8] +"; Dbxref: I5KNAL:"+gene_ID

        except:
            print "error"
    ## second, for non gene feature, remove Dbxref tag
    df2=gff_data.where(gff_data[2].str.contains("gene"))
    dx2 = df2.dropna().index.tolist() #1254
    total=list(range(len(gff_data)))
    ret_list = list(set(dx2)^set(total)) #non gene feature
    for i in ret_list:
       # try:
            if re.search("Dbxref",gff_data.iloc[i,8]):

                if (re.match("^.*(Dbxref=[^;]+);.*$", gff_data.iloc[i, 8]) ):
                    Dbxref= "".join(re.match("^.*(Dbxref=[^;]+).*$", gff_data.iloc[i, 8]).groups())+";"
                    new_dbxref=re.sub(Dbxref,"",gff_data.iloc[i,8])
                    gff_data.iloc[i,8]=new_dbxref
                   # print gff_data.iloc[i,8]
                   # print new_dbxref
                elif (re.match("^.*(Dbxref=[^;]);$", gff_data.iloc[i, 8])):
                    Dbxref = "".join(re.match("^.*(Dbxref=[^;]+);$", gff_data.iloc[i, 8]).groups()) + ";"
                    new_dbxref = re.sub(Dbxref, "", gff_data.iloc[i, 8])
                    gff_data.iloc[i, 8] = new_dbxref
                elif(re.match("^.*(Dbxref=[^;])$", gff_data.iloc[i, 8])):
                    Dbxref = "".join(re.match("^.*(Dbxref=[^;]+)$", gff_data.iloc[i, 8]).groups())
                    new_dbxref = re.sub(Dbxref, "", gff_data.iloc[i, 8])
                    gff_data.iloc[i, 8] = new_dbxref

Remodel_pseudogenes() #2
check_locus()#3
locus_tag_util(Biosample) #3
remove_Note() #5
name_Attributes_gene()#6
name_Attribute_mRNA_pseudo() #7
name_Attribute_trncRNA() #8
Scan_for_genome_alterations() #9
description_attribute() #10
symbol_attribute() #11
ID_attribute() #4
cds_id_attribute() #4
Dbxref() #12

logging.info('Finish program, start generating new gff3 file.')
gff_data.to_csv('clec_OGS_v1_3.pep.gff3', header=False, index=False, sep='\t')

