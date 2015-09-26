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

usage="USAGE:
  ./$(basename $0) /genome/base/directory/
"

script=$(readlink -f $0)

[[ -z $1 ]] && echo -e "$usage\nError: missing folder name" && exit 65

if [[ ! -d ${1} ]]; then mkdir -p ${1}; fi
dest=$(readlink -f ${1}) 
cd ${dest} 

echo "* fetching organism list"
available_organisms=$(curl -#l ftp://ftp.ensembl.org/pub/current_fasta/)
echo "* available organisms:"
printf '%s\n' "$available_organisms" | column

# atlernative: ncftpls -l ftp://ftp.ensembl.org/pub/current_fasta/


printf "* Paste your organism from the list above (eg. saccharomyces_cerevisiae): "
read organism

if [[ -z $(printf ':%s:' $available_organisms | grep -o ":$organism:") ]]; then
  echo "Error: provided organism '$organism' is not in available list"
  exit 65
fi

# check the latest release
echo "* checking latest release"
gtf_import=$(curl -#l ftp://ftp.ensembl.org/pub/current_gtf/$organism/ | grep gtf.gz)
releaseA=$(echo $gtf_import | cut -f2 -d.)
releaseB=$(echo $gtf_import | cut -f3 -d.)
sep=_
release=$releaseA$sep$releaseB 


# ask user if latest relases should be intalled or and older one

printf "* The latest realse is the $release, do you wish to add the (l)atest release or an (o)lder release? (l/o) "
read answer1

if [ $answer1 == 'l' ]; then

  printf "* Setting up the latest release\n"
  path_curl_dna="ftp://ftp.ensembl.org/pub/current_fasta/$organism/dna/"
  path_curl_gtf="ftp://ftp.ensembl.org/pub/current_gtf/$organism/"

elif [ $answer1 == 'o' ]; then

  echo "* fetching release list"
  available_releases=$(curl -#l ftp://ftp.ensembl.org/pub/ | grep release | sort)
  echo "* available releases:"
  printf '%s\n' "$available_releases" | column

  printf "* Paste a release from above (eg. release-66): "
  read answer2
  
  if [[ -z $(printf ':%s:' $available_releases | grep -o ":$answer2:") ]]; then
    echo "Error: provided release '$answer2' is not in available list"
    exit 65
  fi
  old_gtf_path=$answer2/gtf/
  old_fasta_path=$answer2/fasta/
  
  printf "* using FTP server address:\nftp://ftp.ensembl.org/pub/$old_gtf_path$organism/\n"
  
  gtf_import=$(curl -sl ftp://ftp.ensembl.org/pub/$old_gtf_path$organism/ | grep gtf.gz)
  releaseA=$(echo $gtf_import | cut -f2 -d.)
  releaseB=$(echo $gtf_import | cut -f3 -d.)
  sep=_
  release=$releaseA$sep$releaseB 

  path_curl_dna="ftp://ftp.ensembl.org/pub/$old_fasta_path$organism/dna/"
  path_curl_gtf="ftp://ftp.ensembl.org/pub/$old_gtf_path$organism/"

else 

  printf "\nExiting... you are only allowed to give in l or o\n"
  exit

fi

# if the release does not exist on the organism folder add it. Otherwise, ask if the user wants to add more components

if [[ -d $organism/$release ]]
then
  
  cd $organism/$release

  printf "\nYou already have this release in:\n"
  pwd
  printf "\nPress enter to list components.\n"
  read
  ls
  printf "\nDo you wish to add more components? (y/n)\n"
  read answer

  # If the user wants to add more components then they should be implemented into the script.

  if [ $answer == 'y' ]; then

    printf "\nPlease contact me: Jorge.Boucas@age.mpg.de | +49 (0)221 37970 312\n"

  fi
  
  exit

fi

# create directories

mkdir -p $organism/$release
cd $organism/$release

# download reference data

echo "* fetching list of dna files"
files_curl_dna=$(curl -#l $path_curl_dna | grep "\.dna\.")
echo "* fetching list of gtf files"
files_curl_gtf=$(curl -#l $path_curl_gtf | grep '\.gtf\.gz$')
for f in $files_curl_dna; do
  echo "* downloading $f:"
  curl -#O $path_curl_dna/$f
done
for f in $files_curl_gtf; do
  echo "* downloading $f:"
  curl -#O $path_curl_gtf/$f
done

echo "* extracting"
gunzip *.gz

echo "* renaming files"
namep1=$(echo $gtf_import | cut -f1 -d.)

for file in $(ls $namep1.*); do
  nname=${file#"$namep1."}
  nname=${nname%"$namep1."}
  mv $file $nname
done

mkdir chromosomes
for file in $(ls *chromosome*); do
  if [ ${file} != chromosomes: ]; then
    chrp1=$(echo $file | cut -f4 -d.)
    chrp2=$(echo $file | cut -f5 -d.)
    sep=.
    chr=$chrp1$sep$chrp2

    mv ${file} chromosomes/${chr}
  fi 
done

if [[ -n $(ls *nonchromosomal* 2> /dev/null) ]]; then
  mv *nonchromosomal* chromosomes/nonchromosomal.fa
fi

echo "* checking gtf"
for g in $(ls *.gtf); do
  ab=abinitio
  if [ "${g/$ab}" != "${g}" ]; then
    printf "\n${g} is not the correct GTF, looking for correct GTF..\n"
  else
    printf "\nFound ${g}\n"
    gtf=${g}
  fi
done


ori=$(ls *dna.toplevel.fa)
rel=${ori#".dna.toplevel.fa"}
rel=${rel%".dna.toplevel.fa"}
pri=${rel}.dna.primary_assembly.fa

mkdir logs

for g in ${ori} ${pri}; do
  
  if [ "${g}" == "${pri}" ]; then
    original=${pri}
    toplevel=${original#".fa"}
    toplevel=${toplevel%".fa"}
    gLabel="primary_"
  else
    original=${ori}
    toplevel=${original#".fa"}
    toplevel=${toplevel%".fa"}
    gLabel="toplevel_"
  fi
  
  if [ -e ${original} ]; then
    
    echo
    echo "*** running alignments for $g"
    echo

    # BOWTIE2 index
    log="logs/${gLabel}bowtie.log"
    echo "* creating bowtie2 index, logging to $log"
    mkdir ${gLabel}bowtie2
    cd ${gLabel}bowtie2
    ln -s ../${original} ${original}
    module load Bowtie2
    #bowtie2-build-s $original $toplevel 2>&1 | tee ${gLabel}bowtie.log
    bowtie2-build-s $original $toplevel &> ../$log
    echo "$(date)" >> ../$log
    which bowtie2 >> ../$log
    cd ..
    
    # Fix GTF with cuffcompare
    log="logs/cuffcompare.log"
    echo "* fixing GTF file with cuffcompare, logging to $log"
    module load Cufflinks
    #cuffcompare -V -CG -s chromosomes -r $gtf $gtf 2>&1 | tee cuffcompare.log
    cuffcompare -V -CG -s chromosomes -r $gtf $gtf &> $log
    echo "$(date)" >> $log
    which cuffcompare >> $log
    
    mv cuffcmp.combined.gtf cuffcmp_GTF.$gtf
    tar -jcvf cuffcmp.results.tar.bz2 cuffcmp.* --remove-files
    
    # Generate TOPHAT transcriptome indexes
    log="logs/${gLabel}tophat.cuffcomp_GTF.index.log"
    echo "* indexing cuffcompare GTF, logging to $log"
    module load TopHat
    mkdir ${gLabel}tophat_cuffcmp_GTF_index
    #tophat2 -G cuffcmp_GTF.$gtf --transcriptome-index ${gLabel}tophat_cuffcmp_GTF_index ${gLabel}bowtie2/$toplevel 2>&1 | tee ${gLabel}index.log
    tophat2 -G cuffcmp_GTF.$gtf --transcriptome-index ${gLabel}tophat_cuffcmp_GTF_index ${gLabel}bowtie2/$toplevel &> $log
    echo "$(date)" >> $log
    which tophat2 >> $log
    rm -r tophat_out
    
    
    log="logs/tophat.GTF.index.log"
    echo "* indexing GTF, logging to $log"
    mkdir ${gLabel}tophat_GTF_index
    #tophat2 -G $gtf --transcriptome-index ${gLabel}tophat_GTF_index ${gLabel}bowtie2/$toplevel 2>&1 | tee ${gLabel}index.log
    tophat2 -G $gtf --transcriptome-index ${gLabel}tophat_GTF_index ${gLabel}bowtie2/$toplevel &> $log
    echo "$(date)" >> $log
    which tophat2 >> $log
    rm -r tophat_out
    
    
    # BWA index creation 
    
    echo "* Generating BWA index (via slurm)"
    module load BWA
    mkdir ${gLabel}bwa
    cd ${gLabel}bwa
    full_path=$(pwd)
    ln -s ../${original} ${original}  
    echo "#!/bin/bash
bwa index -a bwtsw -p ${full_path}/${original::(-3)} ${original}
echo '$(date)'
which bwa
cp ${gLabel}bwa.log ../logs
" > ${gLabel}bwa.sh
    chmod 770 ${gLabel}bwa.sh; sbatch -p himem,hugemem,blade -o ${gLabel}bwa.log ${gLabel}bwa.sh
    cd ..
    
    
    # STAR index creation
    
    echo "* generating STAR index (via slurm)"
    module load STAR
    mkdir ${gLabel}star
    cd ${gLabel}star
    wSTAR=$(which STARshort)
    IFS='/' read -ra VER <<< "$wSTAR"
    VER=${VER[3]}
    mkdir ${VER}
    cd ${VER}
    full_path=$(pwd)
    ln -s ../../${original} ${original}
    ln -s ../../${gtf} ${gtf}
    echo "#!/bin/bash
STARshort --runMode genomeGenerate --genomeDir ${full_path} --genomeFastaFiles ${original} --runThreadN 18 --sjdbGTFfile ${gtf} --sjdbOverhang 100 --limitGenomeGenerateRAM 90543555797
echo '$(date)'
which STARshort
cp ${gLabel}star.log ../../logs/
" > ${gLabel}star.sh; chmod 770 ${gLabel}star.sh; sbatch --cpus-per-task=18 -p himem,hugemem,blade --mem=100GB -o ${gLabel}star.log ${gLabel}star.sh
    cd ../../
    
    
    # HiSat index creation
    echo "Generating HiSat index (via slurm)"
    mkdir ${gLabel}hisat
    cd ${gLabel}hisat
    full_path=$(pwd)
    ln -s ../${original} ${original}
    echo "#!/bin/bash
module load HISAT
hisat-build ${original} ${original}
echo '$(date)'
which hisat
cp ${gLabel}hisat.log ../logs/
" > ${gLabel}hisat.sh; chmod 770 ${gLabel}hisat.sh; sbatch -p himem,hugemem,blade -o ${gLabel}hisat.log --cpus-per-task=4 ${gLabel}hisat.sh
    cd ..
    
  else
    
    echo "* file '${g}' not found, alignment skipped"
    
  fi
    
done

cp ${script} logs/script.log

exit
