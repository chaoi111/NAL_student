import re
import os
import subprocess
import sys
from sys import argv
import shlex
import logging

logging.basicConfig(level=logging.INFO,filename='logging.txt')
#args= shlex.split(sys.argv)


file2=sys.argv[2]
file1=sys.argv[1]
#ps0= subprocess.call("touch {}_result.txt".format(file1),shell=True)
open('{}_result.txt'.format(file1),'w').close() 
ps1 = subprocess.Popen('samtools view -H {}'.format(file1), shell=True,stdout=subprocess.PIPE)
ps2 = subprocess.Popen("awk '{print $2\"\t\"$3}' -",shell=True,stdin=ps1.stdout,stdout=subprocess.PIPE)
ps1.stdout.close() # matters?
ps3= subprocess.Popen("sed 's/SN://g' -", shell=True,stdin= ps2.stdout,stdout=subprocess.PIPE)
ps2.stdout.close() # why to do this ?
ps4 = subprocess.Popen("sed 's/LN://g' - > {}_result.txt".format(file1), shell=True, stdin=ps3.stdout,stdout=subprocess.PIPE)
ps3.stdout.close()
output = ps4.communicate()[0] 

def bam_dic(bam_out):
   dicts1={}
   dict1={}
   for line in bam_out:

       try:
        a =line.rstrip("\n").split("\t")
      
        dict1={a[0]:a[1]}
        dicts1.update(dict1)
       except:
        print "can't read the line"
   return dicts1

def fasta_dic(fasta_file):
    ccount =0
    nowID=""
    preID=""
    dicts2={}
    dict2={}

    for line in fasta_file:
        line =line.strip("\n")
        if line[0]==">":
            if preID=="": # initialize only once
                preID = "".join(re.match(">(\S+)", line).groups())
                nowID=preID
            preID=nowID
            nowID="".join(re.match(">(\S+)",line).groups())
            dict2={preID:ccount}
            dicts2.update(dict2)
            ccount = 0
            continue

        if re.search("[aAtTcCgGNn]",line[0]):
            ccount = ccount+len(line)
    dicts2.update( {nowID: ccount})
    return dicts2


def run(file1,file2):
  dicts1={}
  dicts2={}
  dicts1= bam_dic(file1)
  dicts2=fasta_dic(file2)
  for key in dicts1:
      if key in dicts2:
          if int(dicts1.get(key))!= int(dicts2.get(key)): # error type1 both id exist
               logging.info('error1_incorrect_length: {}\t{}\t{}'.format(key,dicts1.get(key),dicts2.get(key)))
      elif key not in dicts2:
          logging.info('error2_nonexist_seq: {}\t{}'.format(key,dicts1.get(key)))


if  __name__  == '__main__':
     file2=sys.argv[2]
     with open("{}_result.txt".format(argv[1]))as mid_file:
       with open(file2) as f2:

           run(mid_file,f2)
