#!/bin/python

"""This script annotates 'gene_exp.diff', 'promoters.diff', 'splicing.diff', 'cds.diff', and 'isoform_exp.diff' cuffdiff tables. It generates 1 file for all results, 1 file for p<0.05, and 1 file for q<0.05. For significant values (i.e. q<0.05) it also generates tables containg all pair-wise comparisons in different sheets as well as gene ontology enrichment files for biological processes (BP), cellular component (CC), and moleculr function (MF)."""

import pandas as pd
import numpy as np
import os
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


############### Check available BioMart data ######################
biomaRt = importr("biomaRt") # do not quote
#print(biomaRt.listMarts())
ensemblMart=biomaRt.useMart("ensembl") # do not quote 
#print(biomaRt.listDatasets(ensemblMart))
ensembl=biomaRt.useDataset("celegans_gene_ensembl", mart=ensemblMart) # do not quote
#print(biomaRt.listFilters(ensembl))
#print(biomaRt.listAttributes(ensembl))
###################################################################


############### Set Biomart resources ######################
dataset = 'celegans_gene_ensembl'
data_filter = 'ensembl_gene_id' # filter for gene_id in GTF file, 'link_ensembl_gene_id'
data_output_biotypes = ['ensembl_gene_id','gene_biotype'] # order needs to be kept i.e. id, biotype
data_output_goterms = ['ensembl_gene_id','go_id', 'name_1006'] # order needs to be kept i.e. id, go_id, go_name
############################################################


############### Check DAVID's ID for your set ######################
# http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html#input_list
DAVID_id='WORMBASE_GENE_ID'
####################################################################


################### paths to files #########################
diff_out = '/Users/jboucas/Desktop/Hanjie/cuffdiff_output'
original_gtf = '/Users/jboucas/Desktop/Hanjie/WBcel235.78.gtf'
merged_fixed_gtf = '/Users/jboucas/Desktop/Hanjie/cuffcmp.combined.gtf'
python_output = 'python_out'
os.chdir(diff_out)

if not os.path.exists(python_output):
    os.makedirs(python_output)


########## Get list of gene names and respective ids present in the data set

genes = pd.DataFrame()
for file in ['gene_exp.diff', 'promoters.diff', 'splicing.diff', 'cds.diff', 'isoform_exp.diff']:
    df = pd.read_table(file)
    df = df[['gene']]
    genes = pd.concat([genes,df]).drop_duplicates()
genes = genes.astype(str)
genes = pd.DataFrame(genes.gene.str.split(',').tolist())[0]
genes = genes.drop_duplicates()
genes = genes.tolist()

gtf = pd.read_table(original_gtf, sep='\t', skiprows=6, header=None)
gene_name = pd.DataFrame(gtf[8].str.split('gene_name').tolist())[1]
gene_name = pd.DataFrame(gene_name.str.split(';',1).tolist())
gene_name = pd.DataFrame(gene_name[0].str.split('"').tolist())[1]
            
gene_id = pd.DataFrame(gtf[8].str.split('gene_id').tolist())[1]
gene_id = gene_id.astype(str)
gene_id = pd.DataFrame(gene_id.str.split(';',1).tolist())
gene_id = pd.DataFrame(gene_id[0].str.split('"').tolist())[1]
            
name_id = pd.concat([gene_name, gene_id], axis=1).drop_duplicates()
name_id.columns = ['g_name','g_id']
name_id = name_id[name_id['g_name'].isin(genes)]

genes = name_id['g_id'].tolist()

name_id.to_csv(python_output+'/genes_table.txt', sep="\t",index=False)

del gtf, gene_name, gene_id



# Use R and BioMart to retrieve biotypes and gene ontoloty information

biotypes=biomaRt.getBM(attributes=data_output_biotypes, filters=data_filter, values=genes, mart=ensembl)
goterms=biomaRt.getBM(attributes=data_output_goterms, filters=data_filter, values=genes, mart=ensembl)
bio_go=robjects.r.merge(biotypes, goterms,by=1, all='TRUE')
bio_go.to_csvfile(python_output+'/biotypes_go_raw.txt', quote=False, sep='\t', row_names=False)



# generate biotypes and go terms table using R/biomart output table.

ontology = pd.read_table(python_output+"/biotypes_go_raw.txt")
ontology.columns = ['g_id','gene_biotype','GO_id','GO_term']
ontology = pd.merge(name_id, ontology, how='outer', on='g_id')
ontology = ontology[['g_name','gene_biotype','GO_id','GO_term']]
ontology.columns = ['gene_name','gene_biotype','GO_id','GO_term']

final = pd.DataFrame(columns = ['gene_name','gene_biotype','GO_id','GO_term'])

genes = ontology[['gene_name']]
genes = genes.drop_duplicates()
for gene in list(genes.gene_name):
    ontology_gene = ontology[ontology['gene_name'] == gene]
    ontology_gene_go = ontology_gene[['GO_id','GO_term']]
    if len(ontology_gene_go.index) >= 1:
        ontology_gene_go = ontology_gene_go.transpose()
        ontology_gene_go.to_csv("tmp.txt", sep=";", header=False, index=False)
        ontology_gene_go_number = pd.read_table("tmp.txt", sep ="\t", header=None, nrows = 1)
        ontology_gene_go_name = pd.read_table("tmp.txt", sep ="\t", header=None, skiprows=1, nrows = 1)
        ontology_gene_go = pd.concat([ontology_gene_go_number, ontology_gene_go_name], axis=1)
        ontology_gene_go.columns = ['GO_id','GO_term']
    else:
        ontology_gene_go = pd.DataFrame(columns = ['GO_id','GO_term'])
        
    ontology_gene = ontology_gene[['gene_name','gene_biotype']].drop_duplicates()
    ontology_gene_go = ontology_gene_go.reset_index()
    ontology_gene = ontology_gene.reset_index()
    ontology_gene = pd.concat([ontology_gene, ontology_gene_go], axis = 1)
    ontology_gene = ontology_gene[['gene_name','gene_biotype','GO_id','GO_term']]
    final = pd.concat([final, ontology_gene])
os.remove("tmp.txt")
final.reset_index()
final=final[['gene_name','gene_biotype','GO_id','GO_term']]
final.to_csv(python_output+"/biotypes_go.txt", sep= "\t")

del ontology, genes, ontology_gene, ontology_gene_go, ontology_gene_go_number, ontology_gene_go_name, final



##################### Functions required for getting analysis from DAVID

def DAVIDenrich(listF, idType, bgF='', resF='', bgName = 'Background1',listName='List1', category = '', thd=0.1, ct=2):
    from suds.client import Client
    import os
   
    if len(listF) > 0 and os.path.exists(listF):
        inputListIds = ','.join(open(listF).read().split('\n'))
        print 'List loaded.'        
    else:
        print 'No list loaded.'
        raise

    flagBg = False
    if len(bgF) > 0 and os.path.exists(bgF):
        inputBgIds = ','.join(open(bgF).read().split('\n'))
        flagBg = True
        print 'Use file background.'
    else:
        print 'Use default background.'

    client = Client('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl')
    print 'User Authentication:',client.service.authenticate('Jorge.Boucas@age.mpg.de')

    listType = 0
    print 'Percentage mapped(list):', client.service.addList(inputListIds,idType,listName,listType)
    if flagBg:
        listType = 1
        print 'Percentage mapped(background):', client.service.addList(inputBgIds,idType,bgName,listType)

    print 'Use categories:', client.service.setCategories(category)
    chartReport = client.service.getChartReport(thd,ct)
    chartRow = len(chartReport)
    print 'Total chart records:',chartRow
    
    if len(resF) == 0 or not os.path.exists(resF):
        if flagBg:
            resF = listF + '.withBG.chartReport'
        else:
            resF = listF + '.chartReport'
    with open(resF, 'w') as fOut:
        fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
        for row in chartReport:
            rowDict = dict(row)
            categoryName = str(rowDict['categoryName'])
            termName = str(rowDict['termName'])
            listHits = str(rowDict['listHits'])
            percent = str(rowDict['percent'])
            ease = str(rowDict['ease'])
            Genes = str(rowDict['geneIds'])
            listTotals = str(rowDict['listTotals'])
            popHits = str(rowDict['popHits'])
            popTotals = str(rowDict['popTotals'])
            foldEnrichment = str(rowDict['foldEnrichment'])
            bonferroni = str(rowDict['bonferroni'])
            benjamini = str(rowDict['benjamini'])
            FDR = str(rowDict['afdr'])
            rowList = [categoryName,termName,listHits,percent,ease,Genes,listTotals,popHits,popTotals,foldEnrichment,bonferroni,benjamini,FDR]
            fOut.write('\t'.join(rowList)+'\n')
        print 'write file:', resF, 'finished!'
        
def DAVID_get(cat, filtered_table, all_genes_table):
    IDs_table = pd.merge(filtered_table, all_genes_table, how='left', left_on='identifier', right_on='g_name')
    IDs_table = IDs_table[['g_id']].dropna()
    IDs_table.to_csv('targets_tmp.txt',sep='\t',header=False,index=False)
        
    background = all_genes_table[['g_id']].dropna()
    background.to_csv('background_tmp.txt',sep='\t',header=False,index=False)
    
    DAVIDenrich(listF = './targets_tmp.txt', bgF = './background_tmp.txt', idType = DAVID_id, bgName = 'all_RNAseq_genes', listName = 'changed_genes', category = cat)
    enrich=pd.read_csv('targets_tmp.txt.withBG.chartReport',sep='\t')
    
    os.remove('targets_tmp.txt')
    os.remove('background_tmp.txt')  
    os.remove('targets_tmp.txt.withBG.chartReport')
    
    terms=enrich['Term'].tolist()
    enrichN=pd.DataFrame()
    for term in terms:
        tmp=enrich[enrich['Term']==term]
        tmp=tmp.reset_index(drop=True)
        ids=tmp.xs(0)['Genes']
        ids=pd.DataFrame(data=ids.split(", "))
        ids.columns=['g_id']
        ids['g_id']=ids['g_id'].map(str.lower)
        all_genes_table['g_id']=all_genes_table['g_id'].map(str.lower)
        ids=pd.merge(ids, all_genes_table, how='left', left_on='g_id', right_on='g_id')
        names=ids['g_name'].tolist()
        names = ', '.join(names)
        tmp=tmp.replace(to_replace=tmp.xs(0)['Genes'], value=names)
        enrichN=pd.concat([enrichN, tmp])
    enrichN=enrichN.reset_index(drop=True)    
    
    return enrichN



# create excel report tables

bio_go = pd.read_table(python_output+"/biotypes_go.txt", sep= "\t")

for sig, label in zip([0.05, 2, 'yes'],['diff_p.05.xlsx','diff_all.xlsx','diff_sig.xlsx']):
    writer = pd.ExcelWriter(python_output+'/'+label)
    if sig == 'yes':
        writer_bp = pd.ExcelWriter(python_output+'/bio_process_'+label)
        writer_cc = pd.ExcelWriter(python_output+'/cell_component_'+label)
        writer_mf = pd.ExcelWriter(python_output+'/mol_function_'+label)
    for imp, out, outshort in zip(['gene_exp.diff', 'promoters.diff', 'splicing.diff', 'cds.diff', 'isoform_exp.diff'], ['diff_gene_expression','diff_promoter','diff_splicing','diff_cds','diff_isoforms'], ['geneexp','prom','splic','cds','iso']):
        df = pd.read_table(imp)
        df = df.sort('p_value')
        df = df.sort('q_value')
        if sig == 'yes':
            df = df[df['significant'] == 'yes']
        else:
            df = df[df['p_value'] < sig]
        df = df.reset_index()
        df['gene'] = df['gene'].astype(str)
        tmp = pd.DataFrame(df.gene.str.split(',',1).tolist())
        tmp = pd.DataFrame(tmp.ix[:,0])
        tmp.columns = ['identifier']
        df = pd.concat([df,tmp], axis=1)
        df = pd.merge(df, bio_go, how='left', left_on='identifier', right_on='gene_name')
                     
        if imp == 'isoform_exp.diff': # for isoform_exp.diff we want to have the transcript references
            gtf = pd.read_table(merged_fixed_gtf, sep='\t', skiprows=6, header=None)
            t_id = pd.DataFrame(gtf[8].str.split('transcript_id').tolist())[1]
            t_id = pd.DataFrame(t_id.str.split(';',1).tolist())
            t_id = pd.DataFrame(t_id[0].str.split('"').tolist())[1]
            
            n_ref = pd.DataFrame(gtf[8].str.split('nearest_ref').tolist())[1]
            n_ref = n_ref.astype(str)
            n_ref = pd.DataFrame(n_ref.str.split(';',1).tolist())
            n_ref = pd.DataFrame(n_ref[0].str.split('"').tolist())[1]
            
            id_ref = pd.concat([t_id, n_ref], axis=1).drop_duplicates()
            id_ref.columns = ['transcript_id','nearest_ref']
            
            df = pd.merge(id_ref, df, how='right', left_on='transcript_id', right_on='test_id')
              
        """for significant changes also report overlaps between the days, pair-wise, as well as go ontology enrichemnt for each table from DAVID"""
        if sig == 'yes': 
            
            def DAVID_write(table_of_interest, name_vs_id_table, sheet_name):
                        enr = DAVID_get('GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT', table_of_interest, name_vs_id_table)
                        enr[enr['Category'] == 'GOTERM_BP_FAT'].to_excel(writer_bp, outshort+'_'+sheet_name, index=False)
                        enr[enr['Category'] == 'GOTERM_CC_FAT'].to_excel(writer_cc, outshort+'_'+sheet_name, index=False)
                        enr[enr['Category'] == 'GOTERM_MF_FAT'].to_excel(writer_mf, outshort+'_'+sheet_name, index=False)
            
            sample1 = df[['sample_1']]
            sample1.columns=['samples']
            sample2 = df[['sample_2']]
            sample2.columns=['samples']
            samples=pd.concat([sample1,sample2])
            samples=samples.drop_duplicates()
            samples=samples['samples'].tolist()

            for sample1 in samples:
                for sample2 in samples:
                    if sample1 != sample2:
                        if not os.path.exists(sample1+sample2): 
                            if not os.path.exists(sample2+sample1):

                                os.makedirs(sample1+sample2)
                                df_pair = df[df['sample_1'].isin([sample1,sample2])][df['sample_2'].isin([sample1,sample2])]
            
                                DAVID_write(df_pair, name_id, sample1+'_vs_'+sample2)
            
                                df_pair.drop(['test_id','index','gene_id','Unnamed: 0','identifier','gene_name'], axis=1, inplace=True)
            
                                if imp not in ['gene_exp.diff','isoform_exp.diff']:
                                    df_pair.drop(['value_1','value_2','test_stat'], axis=1, inplace=True)
            
                                df_pair.to_excel(writer, outshort+'_'+sample1+'_vs_'+sample2, index=False)
        
            for sample1 in samples:
                for sample2 in samples:
                    if sample1 != sample2:
                        if os.path.exists(sample1+sample2):
                            os.rmdir(sample1+sample2)
                        elif os.path.exists(sample2+sample1):
                            os.rmdir(sample2+sample1)
                            
        df.drop(['test_id','index','gene_id','Unnamed: 0','identifier','gene_name'], axis=1, inplace=True)
          
        if imp not in ['gene_exp.diff','isoform_exp.diff']:
            df.drop(['value_1','value_2','test_stat'], axis=1, inplace=True)
                
        df.to_excel(writer, outshort+'_'+'ALL', index=False)

    writer.save()
    if sig == 'yes':
        writer_bp.save()
        writer_cc.save()
        writer_mf.save()


exit



