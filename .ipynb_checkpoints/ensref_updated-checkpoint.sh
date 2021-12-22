#!/bin/bash

HOMESOURCE="source ~/.bashrc"
SLURMPARTITION="hooli"
SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.9.sif /bin/bash"

### legal

# copyright (c) 2016, the Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing, Cologne
# <bioinformatics at age dot mpg dot de>

### about

# ensref - build and index genome data bases from ensembl
#
# base folder: 
# - /<base directory>/<organism name>/<release number>/
# base files:
# - BUILD_XX_RELEASE_YY = build and release version info, empty file
# - original.gtf = original reference gtf file
# - original.*.gtf = varying reference gtf file
# - original.*.fa = dna sequences fasta files (e.g. toplevel, ...)
# - cuffcompare.gtf = original.gtf annotated with tophat/cuffcompare
# index files:
# - <dna reference level>_<aligner>, e.g.
# - toplevel_tophat2/index.fa = tophat2 index base name
# - toplevel_tophat2_cuffcompare/index.fa = tophat2/cuffcompare index base name
# log files:
# - log/*.log = stdout and stderr log of index build scripts
# - log/*.sh = index build scripts


### variables and options

unset LIST
unset RELEASE
unset ORGANISM

BIN=$(basename $0)
BASE="/beegfs/common/genomes"
VER_STAR="2.7.3a"
USAGE="USAGE:
  $BIN [OPTIONS]
OPTIONS:
  -h, --help                : show this help and exist
  -l, --list                : return available organisms/release-ids
  -b, --base, --prefix PATH : select the base directory, (default: $BASE)
  -o, --organism ORGANISM   : define organism, e.g. 'caenorhabditis_elegans'
  -r, --release RELEASE     : define release version, e.g. '81'
                              (default: automatic detection of current version)
  --ver-star VERSION        : version of STAR (default $VER_STAR)
  -n, --no-slurm            : use bash to run scripts, not sbatch from slurm
                              (not yet integrated)
  -c, --cpus NUMBER         : define number of CPUs (for STAR, default: 12)
                              (not yet integrated)
"

while [[ $# -gt 0 ]]
do
  case $1 in
    -h|--help)
      echo "$USAGE"
      exit
      ;;
    -l|--list)
      LIST=true
      ;;
    -b|--base|--prefix)
      BASE=$2
      shift
      ;;
    -o|--organism)
      ORGANISM=$2
      shift
      ;;
    -r|--release)
      RELEASE=$2
      shift
      ;;
    --version-star)
      VER_STAR=$2
      shift
      ;;
    -n|--no-slurm)
      unset SLURM
      echo "* warning: option '-n' not yet integrated"
      ;;
    -c|--cpus)
      CPUS=$2
      shift
      echo "* warning: option '-c' not yet integrated"
      ;;
    *)
      echo "error: unsupported option '$1'. See '$(basename $0) --help'."
      exit 65
      ;;
  esac
  shift
done

### functions

get_organism () {
  curl -#l ftp://ftp.ensembl.org/pub/current_fasta/
  # alternative: ncftpls -l ftp://ftp.ensembl.org/pub/current_fasta/
}

get_releases () {
  curl -#l ftp://ftp.ensembl.org/pub/ | grep release | sed 's|release-||' | sort
}

get_current () {
  echo "* getting latest release"
  gtf=$(curl -#l ftp://ftp.ensembl.org/pub/current_gtf/$1/ | grep gtf.gz | uniq)
  relA=$(echo $gtf | cut -f2 -d.)
  relB=$(echo $gtf | cut -f3 -d.)
  sep=_
  rel=$relA$sep$relB 
  echo $rel
}

### check base folder

if [[ ! -d $BASE || ! -x $BASE || ! -w $BASE ]]
then
  echo "error: base directory '$BASE' not found/accessable/writeable"
  exit 65
fi

echo "* fetching organism list"
ORGANISMS=($(get_organism))
echo "* fetching release list"
RELEASES=($(get_releases))
RELEASES=($(for REL in ${RELEASES[@]}; do echo $((REL)); done)) 
CURRENT=$(for REL in ${RELEASES[@]}; do echo $((REL)); done | sort -n | tail -n 1)
echo ${CURRENT}

## print available organisms/release-ids
if [[ $LIST ]]
then
  echo
  echo "* available organisms:"
  printf '%s\n' "${ORGANISMS[@]}" | column
  echo
  echo "* available releases:"
  printf '%s\n' "${RELEASES[@]}" | column
  echo
  exit
fi

### check organism

if [[ -z $ORGANISM || -z $(echo ${ORGANISMS[@]} | grep $ORGANISM) ]]
then
  echo "* error: organism '$ORGANISM' not available"
  exit 65
fi

### check release

if [[ ! $RELEASE ]]
then
  RELEASE=$CURRENT
  echo $RELEASE
fi

if [[ -z $RELEASE || -z $(echo ${RELEASES[@]} | grep -o $RELEASE) ]]
then
  echo "Error: release '$RELEASE' not available"
  exit 65
fi

## Defining paths
URL_GTF="ftp://ftp.ensembl.org/pub/release-$RELEASE/gtf/$ORGANISM/"
URL_DNA="ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/$ORGANISM/dna/"
TARGET=$BASE/$ORGANISM/$RELEASE

### notify

echo "* settings:
  ORG : $ORGANISM
  REL : $RELEASE (current: $CURRENT)
  URL : $URL_GTF
  DIR : $TARGET
"

### generate target directory

if [[ -d $TARGET ]]
then
  echo "* target directory '$TARGET' exists"
fi

mkdir -p $TARGET

if [[ ! -x $TARGET || ! -w $TARGET ]]
then
  echo "error: target directory '$TARGET' not accessable/writeable"
  exit 65
fi

cd $TARGET

### download and extract

echo "* fetching file lists"
FILES_GTF=($(curl -#l $URL_GTF | grep "\.gtf\.gz$"))
FILES_DNA=($(curl -#l $URL_DNA | grep "\.dna\."))
BUILD=$(printf '%s\n' "${FILES_GTF[@]}" | grep "${RELEASE}.gtf.gz")
BUILD=${BUILD#*.}
BUILD=${BUILD%".${RELEASE}.gtf.gz"}
echo ${BUILD}

touch BUILD_${BUILD}_RELEASE_${RELEASE}

echo "* downloading and extracting"

echo "  README (gtf)"
curl -# $URL_GTF/README > README_gtf

shopt -s nocasematch

for FILE in ${FILES_GTF[@]}
do
  FILE_UNZIP=${FILE%%.gz}
  SUB_STR=$(printf '%s\n' "${FILE_UNZIP//${ORGANISM}.${BUILD}.${RELEASE}./}")
  FILE_REBASED=original.${SUB_STR}
  echo "  ${FILE%%.gz} : "
  [[ -f $FILE_REBASED ]] && echo "exists" && continue
  curl -#O $URL_GTF/$FILE
  srun -p ${SLURMPARTITION} gunzip $FILE
  mv $FILE_UNZIP $FILE_REBASED
done

echo "  README (fasta)"
curl -# $URL_DNA/README > README_fa

mkdir -p chromosomes
for FILE in ${FILES_DNA[@]}
do
  FILE_UNZIP=${FILE%%.gz}
  FILE_CHROM=chromosomes/${FILE_UNZIP##*chromosome.}
  SUB_STR=$(printf '%s\n' "${FILE_UNZIP//${ORGANISM}.${BUILD}.dna./}")
  FILE_REBASED=original.${SUB_STR}
  echo "  ${FILE_UNZIP} : "
  [[ -f ${FILE_REBASED} || -f ${FILE_CHROM} ]] && echo "exists" && continue
  [[ $FILE =~ nonchromosomal && -f chromosomes/nonchromosomal.fa ]] && echo "exists" && continue
  curl -#O $URL_DNA/$FILE
  srun -p ${SLURMPARTITION} gunzip $FILE
  if [[ $FILE =~ chromosome ]]
  then
    mv $FILE_UNZIP $FILE_CHROM
  elif [[ $FILE =~ nonchromosomal ]]
  then
    mv $FILE_UNZIP chromosomes/nonchromosomal.fa
  else
    mv $FILE_UNZIP $FILE_REBASED
  fi
done

###############################################################

echo "Alignment for index"

###############################################################

mkdir -p log

gtf=original.gtf

for fa in original.toplevel.fa original.primary_assembly.fa; do

label_fa=${fa#original.}
label_fa=${label_fa%.fa}

#######
### bwa
#######
  
name=${label_fa}_bwa
if [[ -d $name ]]
then
    echo "* indexing $name: skipped (remove folder to rebuild)"
else
    echo "* indexing $name"
    jobname=ensref_${ORGANISM}_${RELEASE}_$name
    cmd="
module purge
module load bwa
#which bwa

date

mkdir $name
cd $name

ln -s ../$fa index.fa
bwa index -a bwtsw -p index.fa index.fa
"
    echo "${cmd}" > log/${name}.sh
    chmod u+x log/${name}.sh

    id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${name}.log
#SBATCH -c 4
#SBATCH --job-name=${jobname}

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

log/${name}.sh

echo "Done bilding bwa index."

SIN
EOF
)
    
fi

#######
### star
#######

if [ "${ORGANISM}" == "homo_sapiens" ]
then
    mem=170GB
    mem_star=170000000000
else
    mem=130GB
    mem_star=125000000000
fi

name=${label_fa}_star_$VER_STAR
if [[ -d $name ]]
then
    echo "* indexing $name: skipped (remove folder to rebuild)"
else
    echo "* indexing $name"
    jobname=ensref_${ORGANISM}_${RELEASE}_$name
    cmd="
module purge
module load star
#which star

date

mkdir $name
cd $name

ln -s ../$fa index.fa
ln -s ../$gtf index.gtf

STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles index.fa \
  --runThreadN 12 --sjdbGTFfile index.gtf --sjdbOverhang 100 \
  --limitGenomeGenerateRAM ${mem_star}
"
    echo "${cmd}" > log/${name}.sh
    chmod u+x log/${name}.sh
    
    id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${name}.log
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --job-name=${jobname}
#SBATCH --mem=${mem}

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

log/${name}.sh

echo "Done bilding STAR index."

SIN
EOF
)

fi

#######
### hisat2
#######
  
name=${label_fa}_hisat2
if [[ -d $name ]]
then
    echo "* indexing $name: skipped (remove folder to rebuild)"
else
    echo "* indexing $name"
    jobname=ensref_${ORGANISM}_${RELEASE}_$name
    cmd="
module purge
module load hisat
#which hisat2

date

mkdir $name
cd $name

ln -s ../$fa index.fa
hisat2-build index.fa index.fa
"
    echo "${cmd}" > log/${name}.sh
    chmod u+x log/${name}.sh
    
    id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${name}.log
#SBATCH -c 4
#SBATCH --job-name=${jobname}

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

log/${name}.sh

echo "Done bilding hisat2 index."

SIN
EOF
)

fi

#######
### bowtie2
#######
  
name=${label_fa}_bowtie2
if [[ -d $name ]]
then
    echo "* indexing $name: skipped (remove folder to rebuild)"
else
    echo "* indexing $name"
    jobname=ensref_${ORGANISM}_${RELEASE}_$name
    cmd="
module purge
module load bowtie

date

mkdir $name
cd $name

ln -s ../$fa index.fa
bowtie2-build index.fa index.fa
"
    echo "${cmd}" > log/${name}.sh
    chmod u+x log/${name}.sh

    id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${name}.log
#SBATCH -c 4
#SBATCH --job-name=${jobname}

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

log/${name}.sh

echo "Done bilding bowtie2 index."

SIN
EOF
)

    srun -d afterok:${id} echo "Done bilding bowtie2 index."
    
fi

#######
### tophat2_cuffcompare
#######
  
# NOTE:
# cuffcompare error 'Warning: could not parse ID or Parent from GFF line'
# can be ignored/bypassed with dropping 'gene' features from the gtf file
# e.g. 
# awk -F '\t' '{if ($3 != "gene") print $0}' $gtf > ${gtf}.tmp
# see the question
# https://www.biostars.org/p/126423/
# and solution
# https://groups.google.com/forum/#!msg/tuxedo-tools-users/FTKA4qozJIc/p47AwnCXxvwJ

name=${label_fa}_tophat2_cuffcompare
if [[ -d $name ]]
then
  echo "* indexing $name: skipped (remove folder to rebuild)"
else
  echo "* indexing $name"
  jobname=ensref_${ORGANISM}_${RELEASE}_$name
  cmd="
module purge
module load tophat
module load cufflinks
module load bowtie
#which cuffliks
#which tophat2

date

mkdir $name
cuffcompare -V -CG -s chromosomes -r $gtf $gtf 2>&1
mv cuffcmp.combined.gtf cuffcompare.gtf

tar -jcvf cuffcompare.results.tar.bz2 cuffcmp.* --remove-files
tophat2 -o ${name}.tmp -G cuffcompare.gtf --transcriptome-index ${name}/index \
  ${label_fa}_bowtie2/index.fa 2>&1

rm -r ${name}.tmp
cd $name

for f in *
do
  [[ \$f =~ ^index\.fa ]] && continue
  mv \$f \${f/#index/index.fa}
done
"
      echo "${cmd}" > log/${name}.sh
      chmod u+x log/${name}.sh
      
      id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${name}.log
#SBATCH -c 4
#SBATCH --job-name=${jobname}

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

log/${name}.sh

echo "Done bilding tophat2_cuffcompare index."

SIN
EOF
)

fi


#######
### tophat2
#######

name=${label_fa}_tophat2
if [[ -d $name ]]
then
  echo "* indexing $name: skipped (remove folder to rebuild)"
else
  echo "* indexing $name"
  jobname=ensref_${ORGANISM}_${RELEASE}_$name
  cmd="
module purge
module load bowtie
module load tophat
#which bowtie2
#which tophat2

date

mkdir $name
tophat2 -G $gtf -o $name --transcriptome-index $name/index ${label_fa}_bowtie2/index.fa 2>&1

cd $name

for f in *
do
  [[ \$f =~ ^index\.fa ]] && continue
  [[ -d \$f ]] && continue
  mv \$f \${f/#index/index.fa}
done
"
    echo "${cmd}" > log/${name}.sh
    chmod u+x log/${name}.sh
    
    id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${name}.log
#SBATCH -c 4
#SBATCH --job-name=${jobname}

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

log/${name}.sh

echo "Done bilding tophat2 index."

SIN
EOF
)

fi

srun -d afterok:${id} echo "** Done with ${fa} **"

#######
### index samtools
#######

cd $TARGET

awk '{print $1 "\t" $2}' original.toplevel.fa.fai > original.toplevel.genome

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output log/${fa%.fa}_samtoolsIndex.log
#SBATCH -c 4
#SBATCH --job-name=${fa%.fa}_samtoolsIndex

echo \${HOSTNAME}
date

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load samtools

samtools faidx ${fa}

SIN
EOF
)

echo "Done bilding samtools index."

done

chmod -R 775 $TARGET

echo "**** Done ****"

exit