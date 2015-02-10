#!/bin/bash

#  all_in_one.sh
#  
#
#  Created by JBoucas on 30/01/15.
#

#rm -r ../tmp
#rm -r ../slurm_logs
#rm -r ../fastqc
#rm -r ../raw_trimmed
#rm -r ../tophat_output
#rm -r ../cufflinks_output
#rm -r ../cuffdiff_output
#rm -r ../cuffmerge_output
rm ../assemblies.txt

mkdir ../tmp
mkdir ../slurm_logs
mkdir ../fastqc
mkdir ../raw_trimmed
mkdir ../tophat_output
mkdir ../cufflinks_output
mkdir ../cuffdiff_output

top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw)/
rawt=$(readlink -f ../raw_trimmed)/
ann=/data/refs/caenorhabditis_elegans/WBcel235_78

# 1st step: fastqc

cd $raw

for file in $(ls *.fq.gz); do echo "#!/bin/bash
cp $raw$file $tmp
cd $tmp
gunzip $file
fastqc -t 2 -o ../fastqc ${file::(-3)}
rm ${file::(-3)}
rm ${file::(-15)}_fastqc.sh" > $tmp${file::(-15)}_fastqc.sh

cd $tmp
chmod 755 ${file::(-15)}_fastqc.sh
rm ${file::(-15)}_fastqc.id
rm ../slurm_logs/fastqc.${file::(-15)}.*.out
sbatch -o ../slurm_logs/fastqc.${file::(-15)}.%j.out ${file::(-15)}_fastqc.sh 2>&1 | tee ${file::(-15)}_fastqc.id

done

pause 3

# capture fastqc job ids

cd $tmp
rm all_fastqc.id
for file in $(ls *_fastqc.id); do
id=$(cat $file | grep 'Submitted batch job')
echo -n :${id:20} >> all_fastqc.id
done
ids=$(cat all_fastqc.id)


# 2nd step: flexbar

cd $raw

for file in $(ls *1_sequence.fq.gz); do echo "#!/bin/bash

flexbar -r $raw${file::(-16)}1_sequence.fq.gz -p $raw${file::(-16)}2_sequence.fq.gz -t ../raw_trimmed/${file::(-17)} -n 32 -a ~/documents/TruSeqAdapters.txt -ao 10 -u 5 -q 20 -m 20 -f i1.8 -ae ANY

cd $top
gzip raw_trimmed/${file::(-17)}_1.fastq
gzip raw_trimmed/${file::(-17)}_2.fastq
rm $tmp${file::(-17)}_flexbar.sh" > $tmp${file::(-17)}_flexbar.sh

cd $tmp
chmod 755 ${file::(-17)}_flexbar.sh
rm ../slurm_logs/flexbar.${file::(-17)}.*.out
rm ${file::(-17)}_flexbar.id
sbatch -d afterok$ids --cpus-per-task=40 --mem=64 -o ../slurm_logs/flexbar.${file::(-17)}.%j.out ${file::(-17)}_flexbar.sh 2>&1 | tee ${file::(-17)}_flexbar.id

done

pause 3

# capture flexbar job ids

cd $tmp
rm all_flexbar.id
for file in $(ls *_flexbar.id); do
id=$(cat $file | grep 'Submitted batch job')
echo -n :${id:20} >> all_flexbar.id
done
ids=$(cat all_flexbar.id)

# 3rd step: tophat and cufflinks

cd $raw

for file in $(ls *1_sequence.fq.gz); do echo "#!/bin/bash

tophat2 -p 40 --library-type fr-firststrand --transcriptome-index $ann/cuffcmp_GTF_index/cuffcmp_GTF.WBcel235.78 -o ../tophat_output/${file::(-17)} $ann/WBcel235.dna.toplevel $rawt${file::(-16)}1.fastq.gz $rawt${file::(-16)}2.fastq.gz

cufflinks -p 40 --library-type fr-firststrand -g $ann/cuffcmp_GTF.WBcel235.78.gtf --no-faux-reads -o ../cufflinks_output/${file::(-17)} ../tophat_output/${file::(-17)}/accepted_hits.bam

rm ${file::(-17)}_th_cl.sh" > $tmp${file::(-17)}_th_cl.sh

cd $tmp
chmod 755 ${file::(-17)}_th_cl.sh
rm ../slurm_logs/th_cl.${file::(-17)}.*.out
rm ${file::(-17)}_th_cl.id
sbatch -d afterok$ids --cpus-per-task=40 -o ../slurm_logs/th_cl.${file::(-17)}.%j.out ${file::(-17)}_th_cl.sh 2>&1 | tee ${file::(-17)}_th_cl.id

done

pause 3

# capture tophat_cufflins job ids

cd $tmp
rm all_th_cl.id
for file in $(ls *_th_cl.id); do
id=$(cat $file | grep 'Submitted batch job')
echo -n :${id:20} >> all_th_cl.id
done
ids=$(cat all_th_cl.id)

# 4th step: select transcripts with full read coverage, cuffmerge, cuffcompare merged GTF with ensembl GTF, cuffdiff

echo "#!/bin/bash

cd $raw
for file in $(ls *1_sequence.fq.gz); do
cd $top
cd tophat_output
mv ${file::(-17)} ${file:16:(-17)}
cd ../cufflinks_output
mv ${file::(-17)} ${file:16:(-17)}
cat ${file:16:(-17)}/transcripts.gtf | grep yes > ${file:16:(-17)}/transcripts_full_read.gtf
echo 'cufflinks_output/${file:16:(-17)}/transcripts_full_read.gtf' >> ../assemblies.txt
done

cd $top

cuffmerge -g $ann/cuffcmp_GTF.WBcel235.78.gtf -s $ann/WBcel235.dna.toplevel.fa -o cuffmerge_output -p 40 assemblies.txt

rm assemblies.txt
mkdir cuffmerge_output/cuffcompare_output
cd cuffmerge_output/cuffcompare_outputt

cuffcompare -V -CG -r $ann/WBcel235.78.gtf -s $ann/chromosomes ../merged.gtf

cd $top
cd tophat_output

cuffdiff -o ../cuffdiff_output -p 40 -L day_0,day_7,day_14,day_21 -u ../cuffmerge_output/cuffcompare_output/cuffcmp.combined.gtf BR1_0/accepted_hits.bam,BR2_0/accepted_hits.bam,BR3_0/accepted_hits.bam BR1_7/accepted_hits.bam,BR2_7/accepted_hits.bam,BR3_7/accepted_hits.bam BR1_14/accepted_hits.bam,BR2_14/accepted_hits.bam,BR3_14/accepted_hits.bam BR1_21/accepted_hits.bam,BR2_21/accepted_hits.bam,BR3_21/accepted_hits.bam

rm cuffmerge_diff.sh
rm -r ../tmp" > cuffmerge_diff.sh

cd $top
cd tophat_output
chmod 755 cuffmerge_diff.sh
rm ../slurm_logs/cuffmerge_diff.*.out
sbatch -d afterok$ids --cpus-per-task=40 -o ../slurm_logs/cuffmerge_diff.%j.out cuffmerge_diff.sh

exit