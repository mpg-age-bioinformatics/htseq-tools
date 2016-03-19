#!/bin/bash

oldFasta=$1
newFasta=$2

module load BLAT
module load kentutils
module load kentTools

mkdir ../tmp
mkdir ../slurm_logs

top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
logs=$(readlink -f ../slurm_logs)/

#rm -rf ../chains_output

mkdir ../chains_output
mkdir ../chains_output/oldFasta
mkdir ../chains_output/newFasta
mkdir ../chains_output/bit

cp ${oldFasta} ../chains_output/oldFasta
cp ${newFasta} ../chains_output/newFasta

oldFasta=$(readlink -f ../chains_output/oldFasta/$(basename ${oldFasta}))
newFasta=$(readlink -f ../chains_output/newFasta/$(basename ${newFasta}))

old=$(readlink -f ../chains_output/oldFasta)/
new=$(readlink -f ../chains_output/newFasta)/

ChunkSize="3000"

cd ${old}
#srun faSplit about -lift ${oldFasta} ${ChunkSize} fragment
srun faToTwoBit ${oldFasta} $(basename ${oldFasta} .fa).2bit
srun twoBitInfo $(basename ${oldFasta} .fa).2bit chrom.sizes # stdout | sort -k2nr >

cd ${new}
#srun faSplit size -lift=$(basename ${newFasta} .fa).lft ${newFasta} ${ChunkSize} fragment
srun faSplit size -lift=$(basename ${newFasta} .fa).lft -oneFile ${newFasta} ${ChunkSize} fragment
srun faToTwoBit ${newFasta} $(basename ${newFasta} .fa).2bit
srun twoBitInfo $(basename ${newFasta} .fa).2bit chrom.sizes

cd ${old}
srun blat $(basename ${oldFasta} .fa).2bit /dev/null /dev/null -tileSize=11 -makeOoc=11.ooc -repMatch=100
ooc=$(readlink -f 11.ooc)

cd ${new}
IDS=
for i in $(ls fragment*.*);
    do echo "#!/bin/bash
    cd ${new}
    echo 'starting blat'
    blat ${old}$(basename ${oldFasta} .fa).2bit ${i} OLD.$(basename ${i} .fa).psl -tileSize=11 -ooc=${ooc} -minScore=100 -minIdentity=98 -fastMap
    echo 'starting liftUp'
    liftUp -pslQ $(basename ${i} .fa).psl ${new}$(basename ${newFasta} .fa).lft warn OLD.$(basename ${i} .fa).psl
    echo 'starting axtChain'
    axtChain -linearGap=medium -psl $(basename ${i} .fa).psl ${old}$(basename ${oldFasta} .fa).2bit ${new}$(basename ${newFasta} .fa).2bit $(basename ${i} .fa).chain
    rm ${tmp}$(basename ${i} .fa).sh
    exit
    " > ${tmp}$(basename ${i} .fa).sh
    chmod 755 ${tmp}$(basename ${i} .fa).sh
    rm ${logs}$(basename ${i} .fa).*.out > /dev/null 2&>1
    ID=$(sbatch -o ${logs}$(basename ${i} .fa).%j.out ${tmp}$(basename ${i} .fa).sh)
    IDS=${IDS}:${ID:20}
    echo ${ID}
done
echo "Blat jobs running: ${IDS}"
echo "#!/bin/bash
cd ${new}
mkdir chains
mv *.chain chains
DIR=chains
BATCH_SIZE=100
SUBFOLDER_NAME=chains.part
COUNTER=1

# Merge chain files
while [ `find $DIR -maxdepth 1 -type f| wc -l` -gt $BATCH_SIZE ] ; do
    NEW_DIR=$DIR/${SUBFOLDER_NAME}${COUNTER}
    mkdir $NEW_DIR
    find $DIR -maxdepth 1 -type f | head -n $BATCH_SIZE | xargs -I {} mv {} $NEW_DIR
    cd $DIR/${SUBFOLDER_NAME}${COUNTER}
    chainMergeSort *.chain > ${SUBFOLDER_NAME}${COUNTER}.chain
    cd ${new}
    let COUNTER++
done
NEW_DIR=$DIR/${SUBFOLDER_NAME}${COUNTER}
mkdir ${NEW_DIR}
mv chains/*.chain ${NEW_DIR}
cd ${NEW_DIR}
chainMergeSort *.chain > ${SUBFOLDER_NAME}${COUNTER}.chain
cd ${new}${DIR} 
for i in $(ls -d ${SUBFOLDER_NAME}*);
    do mv ${i}/${i}.chain ${new}${DIR}
done

chainMergeSort ${SUBFOLDER_NAME}*.chain | chainSplit -lump=100 chainSplit stdin
endsInLf chainSplit/*.chain

mkdir netSplit
mkdir overSplit
for f in $(ls chainSplit/*.chain);
    do chainNet ${f} ${old}chrom.sizes ${new}chrom.sizes netSplit/$(basename ${f} .chain).net /dev/null
    netChainSubset netSplit/$(basename ${f} .chain).net ${f} stdout | chainStitchId stdin overSplit/$(basename ${f} .chain).chain
done
endsInLf netSplit/*.net
endsInLf overSplit/*.chain

cat overSplit/*.chain | gzip -c > ${top}chains_output/$(basename ${oldFasta} .fa)_To_$(basename ${newFasta} .fa).chain.gz
exit
" > ${tmp}final.chain.sh
chmod 755 ${tmp}final.chain.sh
rm ${logs}final.chain.*.out > /dev/null 2&>1
sbatch -d afterok${IDS} -o ${logs}final.chain.%j.out ${tmp}final.chain.sh
exit
