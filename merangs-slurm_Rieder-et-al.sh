#!/bin/bash

# merangs-slurm.sh - a meRanTK pipeline for slurm HPC systems
#
# This version implements the alignment against the ensembl genome
# using the meRanGs (aligner = star) program.
#
# Author: Sven E. Templer <sven.templer@gmail.com>
#
# Copyright (c) 2016 Bioinformatics Core Facility of the
# Max Planck Institute for Biology of Ageing, Cologne, Germany

unset bin dir pid
declare -A bin  # software binaries
declare -A dir  # output directories
declare -A pid  # slurm process ids

################################################################################
# USER INPUT
################################################################################

### root project directory to put output folder into
# all current paths are relative to 'home' or absolute
home="/beegfs/group_bit/data/projects/departments/Xiangru_Xu/XX_RNAmeth"
### notification mail address
email="$USER@age.mpg.de"
### suffix
# added to output folders, example:
# suffix '-testing' results in ./merangs-testing_align/ instead ./merangs_align
suffix="Khoddami"
### igenome
# run `curl ftp.ncbi.nlm.nih.gov/refseq/` to check available transcriptomes
#genome="ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz"
# ensembl genome
genome_gtf="/beegfs/common/genomes/mus_musculus/83/original.gtf"
genome_fa="/beegfs/common/genomes/mus_musculus/83/original.toplevel.fa"
# adapters to trim
adapter="/beegfs/common/genomes/adapters/All.fasta"
### sample information
# space separated strings, comments '#.*' and empty lines are ignored,
# columns are:
# * 'name' for sample names
# * 'group' for grouping the samples (2 levels required)
# * 'file1' single end reads fastq file
# * 'file2' paired reads, optional
# the template comes from the public data used in Rieder et al. 2015 on
# the data of Khoddami & Cairns 2013
samples="
# name group file1 (file2)
s1 wt reads_public_GSE44359/SRR727515_1.fastq # used in Rieder et al.
#s2 wt reads_public_GSE44359/SRR727516_1.fastq
#s3 wt reads_public_GSE44359/SRR727517_1.fastq
#s4 wt reads_public_GSE44359/SRR727518_1.fastq
s5 dnmt2 reads_public_GSE44359/SRR727519_1.fastq # used in Rieder et al.
#s6 dnmt2 reads_public_GSE44359/SRR727520_1.fastq
#s7 dnmt2 reads_public_GSE44359/SRR727521_1.fastq
#s8 dnmt2 reads_public_GSE44359/SRR727522_1.fastq
"
# read size (used for -sj0 in align/index, and for -rl in call)
readsize=100
### software paths
# absolute paths, or relative to home
# path to meRanTK software v 1.1.1a
meran="software/meRanTK_1.1.1a"
bin[meranmap]="$meran/util/mkRefSeq2GeneMap.pl"
bin[merangs]="$meran/meRanGs"
bin[merancall]="$meran/meRanCall"
bin[merancompare]="$meran/meRanCompare"
bin[sam]="$meran/extutil/x86_64-linux-thread-multi/tophat2/samtools_0.1.18"
bin[flexbar]="/software/Flexbar/2.5/flexbar"
bin[fqdump]="/software/sratoolkit/2.4.5-2/bin/fastq-dump"
bin[htseq]="htseq-count"

################################################################################
# functions
################################################################################

quit () {
  echo -e "error: $@"
  exit 65
}
file_check () {
  [[ -z $1 ]] && return
  [[ ! -f $1 ]] && quit "missing file '$1'"
}
cat_array () {
  local IFS=','
  echo "$*"
}
in_array () {
  local IFS=':'
  local v="$1"
  shift
  echo ":$*:" | grep -o ":$v:"
}
submit_sbatch () {
  local sout=slurm/${1}.out
  local spid=slurm/${1}.pid
  local sbat=slurm/${1}.sbatch
  local pid=""
  [[ -z $1 ]] && quit "missing arguments for submit_sbatch"
  mkdir -p $(dirname $sout)
  echo -n "* $sbat " >&2
  if [[ -f $spid ]]
  then
    echo -n "already submitted, skipped" >&2
  else
    echo '#!/bin/bash' > $sbat
    echo -e "#SBATCH -o $sout\n#SBATCH -J $1\n$2" >> $sbat
    echo test > $spid
    sbatch --parsable $sbat > $spid
    [[ $? -gt 0 ]] && quit "sbatch error"
    echo -n "submitted" >&2
  fi
  pid=$(cat $spid)
  echo " ($pid)" >&2
  echo $pid
  return
}
printsamples () {
  local what=$1
  local line=$2
  local col=""
  case $what
  in
    name) col=1 ;;
    group) col=2 ;;
    file1) col=3 ;;
    file2) col=4 ;;
  esac
  echo "$samples" | \
    sed 's/#.*//' | \
    egrep -v "^\ *$" | \
    awk -v col=$col '{print $col}' | \
    sed -n "${line}p"
}

################################################################################
# check
################################################################################

# directories
dir[index]="merangs${suffix}_index"
dir[align]="merangs${suffix}_align"
dir[call]="merangs${suffix}_call"
dir[compare]="merangs${suffix}_compare"
dir[trim]="merangs${suffix}_trim"
dir[htseq]="merangs${suffix}_htseqcount"
[[ ! -w $home ]] && quit "'$home' is not writeable"
cd $home
mkdir -p ${dir[@]}

# samples
[[ -z $(printsamples) ]] && quit "samples empty/not defined"
for f in $(printsamples file1) $(printsamples file2)
do
  file_check $f
done
groups=($(echo "$(printsamples group)" | sort | uniq))
[[ ${#groups[@]} -ne 2 ]] && quit "number of groups is not 2"
si=$(seq $(printsamples | wc -l))
comparison=${groups[0]}_vs_${groups[1]}

# user settings
[[ -z $email ]] && quit "provide an email address"

# software
sbin="sbatch gunzip wget realpath awk"
for s in ${sbin}
do
  [[ -z $(which $s) ]] && quit "system software '$s' missing in \$PATH"
done
for i in ${!bin[@]}
do
  s=${bin[$i]}
  if [[ -f $s ]]
  then
    b=$(realpath $s)
  else
    b=$(which $s)
  fi
  [[ ! -x $b ]] && quit "software '$s' missing or not executable"
  bin[$i]="$b"
done

echo -e "\n* samples\n"
echo -e "$(printsamples)\n" | sed 's/^/\ \ /'
echo
echo "* checking jobs"
echo

################################################################################
# public data from Khoddami & Cairns 2013
################################################################################

# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44359
dir[data]="reads_public_GSE44359"
dir[dataftp]="ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP018/SRP018704"
mkdir -p ${dir[data]}

### download

echo "* retrieving public data"
if [[ -d ${dir[dataftp]} ]]
then
  echo -e "  directory '${dir[dataftp]}' exists, skipped download"
else
  echo "  downloading to '${dir[dataftp]}'"
  wget -r ftp://${dir[dataftp]}
fi

### dump

sra=$(find ${dir[dataftp]} -type f -name "*sra")
for s in $sra
do
  n=$(basename $s)
  job=${dir[data]}/$n
  sbatch="${bin[fqdump]} --split-files -O ${dir[data]} $s"
  pid[data]+=":$(submit_sbatch $job "$sbatch")"
done

################################################################################
# index
################################################################################

# limitGenomeGenerateRAM not less than124544990592

### download

# used local mirror of ensembl














### prepare

job=${dir[index]}/prepare
sbatch="
cd ${dir[index]}
echo collecting ids from 
echo . $genome_fa
awk '/^>/ {print \$1}' $genome_fa | sort | uniq > ids.fa
echo . $genome_gtf
awk '/^[^#]/ {print \$1}' $genome_gtf | sort | uniq | sed 's/^/>/' > ids.gtf
echo merging ids
comm -12 ids.gtf ids.fa > ids.intersect
echo subsetting 
echo . $genome_fa
awk 'FNR==NR {a[\$1]++; next} /^>/ {flag=0} \$1 in a {flag=1} flag {print}' \
  ids.intersect $genome_fa > index.fa
echo . $genome_gtf
grep -f <(sed 's/^>/^/' ids.gtf) $genome_gtf > index.gtf
echo OK
"
pid[index_prepare]=":$(submit_sbatch $job "$sbatch")"

### index

job=${dir[index]}/index
sbatch="
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
#SBATCH -d afterok${pid[index_prepare]}
#SBATCH -n 18

${bin[merangs]} mkbsidx \
  -sjO $((readsize - 1)) \
  -GTFtagEPT transcript_id \
  -GTFtagEPG gene_id \
  -star_limitGenomeGenerateRAM 350000000000 \
  -t 18 \
  -id ${dir[index]} \
  -fa ${dir[index]}/index.fa \
  -GTF ${dir[index]}/index.gtf
"
pid[index]=":$(submit_sbatch $job "$sbatch")"


################################################################################
# trimming/filtering
################################################################################

cat $adapter > ${dir[trim]}/adapter.fa
# add the adapter Rieder et al. trimmed
echo -e ">Rieder et. al\nTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGC" >> ${dir[trim]}/adapter.fa

for i in $si
do
  name=$(printsamples name $i)
  f1=$(printsamples file1 $i)
  f2=$(printsamples file2 $i)
  fq="-r $f1"
  [[ -n $f2 ]] && fq+=" -p $f2"
  job=${dir[trim]}/$name
  sbatch="
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
#SBATCH -d afterok${pid[data]}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 18

mkdir $job

${bin[flexbar]} \
  -t $job \
  -n 18 \
  -a ${dir[trim]}/adapter.fa \
  -ao 10 \
  -u 5 \
  -q 30 \
  -m 40 \
  -f i1.8 \
  -ae ANY \
  $fq
"
  pid[trim,$name]=":$(submit_sbatch $job "$sbatch")"
done

################################################################################
# align
################################################################################

for i in $si
do
  name=$(printsamples name $i)
  f1=$(printsamples file1 $i)
  f2=$(printsamples file2 $i)
  if [[ -z $f2 ]]
  then
    f1=${dir[trim]}/${name}.fastq
    fq="-f $f1"
  else
    f1=${dir[trim]}/${name}_1.fastq
    f2=${dir[trim]}/${name}_2.fastq
    fq="-f $f1 -r $f2"
  fi
  job=${dir[align]}/$name
  sbatch="
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
#SBATCH -d afterok${pid[index]}${pid[data]}${pid[trim,$name]}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 18

mkdir -p $job

${bin[merangs]} align \
  -t 18 \
  --star_outFilterScoreMinOverLread 0.9 \
  --star_outFilterMatchNminOverLread 0.9 \
  --star_outFilterMismatchNoverLmax 0.05 \
  --star_outFilterMismatchNmax 5 \
  -MM \
  -star_outFilterMultimapNmax 20 \
  -bg \
  -mbgc 10 \
  -mbp \
  -sjO $((readsize - 1)) \
  -o $job \
  -S aligned.sam \
  -un -ud $job/unaligned \
  -id ${dir[index]} \
  $fq

if [[ ! -f $job/aligned_cut.sam ]]; then
  cut -f1-11 $job/aligned.sam > $job/aligned_cut.sam
fi
if [[ ! -f $job/aligned.bam ]]; then
  ${bin[sam]} view -bS -o $job/aligned.bam $job/aligned_cut.sam
fi
if [[ ! -f $job/aligned_sorted.bam ]]; then
  ${bin[sam]} sort $job/aligned.bam $job/aligned_sorted
fi
if [[ ! -f $job/aligned_sorted.bai ]]; then
  ${bin[sam]} index $job/aligned_sorted.bam $job/aligned_sorted.bai
fi
"
  pid[align,$name]=":$(submit_sbatch $job "$sbatch")"
done

################################################################################
# htseq count
################################################################################

mkdir -p ${dir[htseq]}
for i in $si
do
  name=$(printsamples name $i)
  job=${dir[htseq]}/$name
  bam=${dir[align]}/$name/aligned_sorted.bam
  sbatch="
#SBATCH -d afterok${pid[align,$name]}
${bin[htseq]} -f bam -m intersection-strict -a 10 -t exon -i gene_id $bam $genome_gtf > $job
"
  pid[htseq]+=":$(submit_sbatch $job "$sbatch")"
done

################################################################################
# call
################################################################################

for i in $si
do
  name=$(printsamples name $i)
  job=${dir[call]}/$name
  sbatch="
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
#SBATCH -d afterok${pid[align,$name]}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
${bin[merancall]} \
  -p 8 \
  -fs5 0 \
  -fs3 1 \
  -rs5 0 \
  -rs3 0 \
  -tref \
  -rl $readsize \
  -np \
  -bed63 \
  -sc 5 \
  -cr 1 \
  -fdr 0.01 \
  -ei 0.1 \
  -md 10 \
  -o $job \
  -bam ${dir[align]}/$name/aligned_sorted.bam \
  -f ${dir[index]}/index.fa
"
  #-fs5 5 -rs5 5 -rl 50 -sc 10 -gref -fdr 0.05 \
  pid[call]+=":$(submit_sbatch $job "$sbatch")"
done

################################################################################
# compare
################################################################################

job=${dir[compare]}
src=${dir[call]}
pth=$(realpath ${bin[merancompare]})
fa=($(printsamples | awk -v g=${groups[0]} '$2==g {print $1}'))
fa=($(printf "$src/%s_FDR_0.01.txt " ${fa[@]}))
fa=($(realpath ${fa[@]}))
fa=$(cat_array ${fa[@]})
fb=($(printsamples | awk -v g=${groups[1]} '$2==g {print $1}'))
fb=($(printf "$src/%s_FDR_0.01.txt " ${fb[@]}))
fb=($(realpath ${fb[@]}))
fb=$(cat_array ${fb[@]})

sbatch="
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
#SBATCH -d afterok${pid[call]}

cd $job
cat > samples.$comparison << EOF
$(printsamples)
EOF
$pth \
  -mr 3 -s 0.01 -fdr 0.01 \
  -na ${groups[0]} \
  -nb ${groups[1]} \
  -fa $fa \
  -fb $fb
"
pid[compare]=":$(submit_sbatch $job "$sbatch")"

################################################################################
# notify
################################################################################

pid[notification]=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -J merangs${suffix}_done
#SBATCH --mail-type=END
#SBATCH --mail-user=$email
#SBATCH -d afterany${pid[compare]}${pid[htseq]}
echo done
EOF
)

echo
echo "* done - notification email will be send to '$email'"
echo "  (after last job ended and if any job failed)"
echo

