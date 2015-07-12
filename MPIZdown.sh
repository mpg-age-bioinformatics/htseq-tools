#!/bin/bash

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -l|--links_file)
    f=${2}
    shift # past argument
    ;;
    -o|--output_folder)
    folder=${2}
    shift # past argument
    ;;
    *)
    ;;
esac
shift # past argument or value
done

f=$(readlink -f ${f})

if [[ ! -e ${folder} ]]; then
mkdir -p ${folder}; fi
folder=$(readlink -f ${folder}) 

cd ${folder}

while read line; do down=${line:0:4}; if [[ ${down} == "http" ]]; then wget ${line}; fi; done < ${f}
while read line; do md5=${line:0:6}; if [[ ${md5} == "md5sum" ]]; then echo ${line:7} >> md5.tmp; fi; done < ${f}
md5sum -c md5.tmp 2>&1 | tee md5.log

echo $(date) >> downloads.log
cat ${f} >> downloads.log
cat md5.log >> downloads.log

rm md5.log md5.tmp

exit

