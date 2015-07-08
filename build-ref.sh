#!/bin/bash

# This script builds the references and annotations folder for each selected organism and selected release. 
# eg. caenorhabditis_elegans release: WBcell235.78
# Level1 folder = caenorhabditis_elegans
# Level2 folder = WBcell235_78
# reference genome (fasta and index) = WBcel235.dna.toplevel
# reference GTF = WBcel235.78.gtf
# reference GTF index = caenorhabditis_elegans/WBcell235_78/GTF_index/
# cuffcompare fixed GTF = cuffcmp_GTF.WBcel235.78.gtf
# other cuffcompare output files = cuffcmp.results.tar.bz2
# cuffcompare fixed GTF index = caenorhabditis_elegans/WBcell235_78/GTF_cuffcmp_GTF_index/
# chromosomes fasta = caenorhabditis_elegans/WBcell235_78/fasta

# example tophat options: 
# --transcriptome-index $references_directory/caenorhabditis_elegans/WBcell235_78/GTF_cuffcmp_GTF_index 
# for the reference indexed genome use $references_directory/caenorhabditis_elegans/WBcell235_78/bowtie2/WBcel235.dna.toplevel

# example cufflinks options:
# -g $references_directory/caenorhabditis_elegans/WBcell235_78/cuffcmp_GTF.WBcel235.78.gtf

# example cuffcompare options:
# -s $references_directory/caenorhabditis_elegans/WBcell235_78/chromosomes  

curl -l ftp://ftp.ensembl.org/pub/current_fasta/

printf "Paste your organism from the list above: "
read organism

if [ ! -d "organism" ]; then mkdir $organism; fi

cd $organism

# check the latest release
gtf_import=$(curl -l ftp://ftp.ensembl.org/pub/current_gtf/$organism/ | grep gtf.gz)
releaseA=$(echo $gtf_import | cut -f2 -d.)
releaseB=$(echo $gtf_import | cut -f3 -d.)
sep=_
release=$releaseA$sep$releaseB 


# ask user if latest relases should be intalled or and older one

printf "The latest realse is the $release, do you wish to add the (l)atest release or an (o)lder release? (l/o) "
read answer1

if [ $answer1 == 'l' ]; then

printf "Setting up the latest release"

elif [ $answer1 == 'o' ]; then

curl -l ftp://ftp.ensembl.org/pub/ | grep release

printf "Paste a release from above (eg. release-66): "
read answer2

old_gtf_path=$answer2/gtf/
old_fasta_path=$answer2/fasta/

printf "ftp://ftp.ensembl.org/pub/$old_gtf_path$organism/"

gtf_import=$(curl -l ftp://ftp.ensembl.org/pub/$old_gtf_path$organism/ | grep gtf.gz)
releaseA=$(echo $gtf_import | cut -f2 -d.)
releaseB=$(echo $gtf_import | cut -f3 -d.)
sep=_
release=$releaseA$sep$releaseB 

else 

printf "
Exiting... you are only allowed to give in l or o
"

exit

fi

# if the release does not exist on the organism folder add it. Otherwise, ask if the user wants to add more components

if [ ! -d "$release" ]; then

mkdir $release
cd $release

if [ $answer1 == 'l' ]; then

wget ftp://ftp.ensembl.org/pub/current_fasta/$organism/dna/*.dna.*
wget ftp://ftp.ensembl.org/pub/current_gtf/$organism/*.gtf.gz

else if [ $answer1 == 'o' ]; then

wget ftp://ftp.ensembl.org/pub/$old_fasta_path$organism/dna/*.dna.*
wget ftp://ftp.ensembl.org/pub/$old_gtf_path$organism/*.gtf.gz

fi; fi

gunzip *.gz

namep1=$(echo $gtf_import | cut -f1 -d.)

for file in $(ls $namep1.*); do
nname=${file#"$namep1."}
nname=${nname%"$namep1."}
mv $file $nname; done

mkdir chromosomes
for file in $(ls *chromosome*); do
if [ ${file} != chromosomes: ]; then
chrp1=$(echo $file | cut -f4 -d.)
chrp2=$(echo $file | cut -f5 -d.)
sep=.
chr=$chrp1$sep$chrp2

mv ${file} chromosomes/${chr}; fi; done

mv *nonchromosomal* chromosomes/nonchromosomal.fa

original=$(ls *dna.toplevel.fa)
toplevel=${original#".fa"}
toplevel=${toplevel%".fa"}


# BOWTIE2 index
mkdir bowtie2
cd bowtie2
ln -s ../${original} ${original}
module load Bowtie2
bowtie2-build-s $original $toplevel 2>&1 | tee bowtie.log
which bowtie2 >> bowtie.log

cd ..
gtf=$(ls *.gtf)

# Fix GTF with cuffcompare

module load Cufflinks
cuffcompare -V -CG -s chromosomes -r $gtf $gtf 2>&1 | tee cuffcompare.log
which cuffcompare >> cuffcompare.log

mv cuffcmp.combined.gtf cuffcmp_GTF.$gtf
tar -jcvf cuffcmp.results.tar.bz2 cuffcmp.* --remove-files

# Generate TOPHAT transcriptome indexes


printf "
Indexing cuffcompare GTF
"

module load TopHat
mkdir cuffcmp_GTF_index
tophat2 -G cuffcmp_GTF.$gtf --transcriptome-index cuffcmp_GTF_index bowtie2/$toplevel 2>&1 | tee index.log
which tophat2 >> index.log
mv index.log cuffcmp_GTF_index/index.log
rm -r tophat_out


printf "
Indexing GTF
"

mkdir GTF_index
tophat2 -G $gtf --transcriptome-index GTF_index bowtie2/$toplevel 2>&1 | tee index.log
which tophat2 >> index.log
mv index.log GTF_index/index.log
rm -r tophat_out


# BWA index creation 

printf "
Generating BWA index
"

module load BWA
mkdir bwa
cd bwa
full_path=$(pwd)
ln -s ../${original} ${original}  
echo "#!/bin/bash
bwa index -a bwtsw -p ${full_path}/${original::(-3)} ${original}
which bwa" > bwa.sh
chmod 770 bwa.sh; sbatch -p himem,hugemem,blade -o bwa.log bwa.sh
cd ..


# STAR index creation

printf "
Generating STAR index
"

module load STAR
mkdir star
cd star
full_path=$(pwd)
ln -s ../${original} ${original}
ln -s ../${gtf} ${gtf}
echo "#!/bin/bash
STAR --runMode genomeGenerate --genomeDir ${full_path} --genomeFastaFiles ${original} --runThreadN 20 --sjdbGTFfile ${gtf} --sjdbOverhang 100 --limitGenomeGenerateRAM 240000000000
which STAR
" > star.sh; chmod 770 star.sh; sbatch --cpus-per-task=20 -p himem,hugemem,blade --mem=256GB -o star.log star.sh
cd ..

else

cd $release

printf "
You already have this release in:
"
pwd
printf "
Press enter to list components.
"
read
ls
printf " 
Do you wish to add more components? (y/n) 
"
read answer

# If the user wants to add more components then they should be implemented into the script.

if [ $answer == 'n' ]; then exit;

elif [ $answer == 'y' ]; then

printf "
Please contact me: Jorge.Boucas@age.mpg.de | +49 (0)221 37970 312
"; 

else

printf "
Exiting.. you can only type y or n
"
exit

fi; fi

exit
