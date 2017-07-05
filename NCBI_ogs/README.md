## NCBI_OGS
*https://gitlab.com/i5k_Workspace/workspace_roadmap/issues/406
*clean Official Gene Sets prior to NCBI submission
## Background
-We need to submit the Official Gene Sets generated at the i5k Workspace to NCBI/GenBank. To do so, we need to 1) clean up of the original gff, and 2) run NCBI's new table2asn program. This issue regards step 1 - cleaning up the gff. It would be good to write a program to do as much of the cleaning up as we can. Below I list steps for this program, but these steps are likely to change as we work on our input files and get feedback from NCBI. We will also have to do some manual fixes. Nonetheless, the output of this issue should be a mostly functioning program to clean gffs, as described below. 

-More information on NCBI's requirements for the gff3 file: https://www.ncbi.nlm.nih.gov/sites/genbank/genomes_gff/.
## Function
2. Remodel_pseudogenes() 
3. check_locus()
3. locus_tag_util(Biosample) 
4. ID_attribute() 
4. cds_id_attribute() 
5. remove_Note() 
6. name_Attributes_gene()
7.name_Attribute_mRNA_pseudo() 
8. name_Attribute_trncRNA() 
9. Scan_for_genome_alterations() 
10. description_attribute() 
11. symbol_attribute() 
12. Dbxref() 
