#!/bin/bash
echo -e "    ________  _______                                                     \n   / ____/ / / /  _( )_____                                               
  / /_  / /_/ // / |// ___/                                               \n / __/ / __  // /   (__  )                                                
/_/___/_/_/_/___/__/____/__       ______    _    __     ___               \n  / ___//   |  / __ \/ ___/      / ____/___| |  / /    |__ \              
  \__ \/ /| | / /_/ /\__ \______/ /   / __ \ | / /_______/ /              \n ___/ / ___ |/ _, _/___/ /_____/ /___/ /_/ / |/ /_____/ __/               
/____/_/ _|_/_/ |_|/____/______\____/\____/|___/ ____/____/______         \n| |     / /   | / ___/_  __/ ____/ |     / /   |/_  __/ ____/ __ \        
| | /| / / /| | \__ \ / / / __/  | | /| / / /| | / / / __/ / /_/ /        \n| |/ |/ / ___ |___/ // / / /___  | |/ |/ / ___ |/ / / /___/ _, _/         
|__/|__/_/_ |_/____//_/ /_____/__|__/|__/_/__|_/_/_/_____/_/_|_|__________\n  / ___// / / / __ \ |  / / ____/  _/ /   / /   /   |  / | / / ____/ ____/
  \__ \/ / / / /_/ / | / / __/  / // /   / /   / /| | /  |/ / /   / __/   \n ___/ / /_/ / _, _/| |/ / /____/ // /___/ /___/ ___ |/ /|  / /___/ /___   
/____/\____/_/ |_| |___/_____/___/_____/_____/_/  |_/_/ |_/\____/_____/   \n                                                                          "

echo""
echo "Extracting lineage "${1}" from fastq reads"
echo "Stay tuned for updates!" 
echo "Visit us at https://github.com/folkehelseinstituttet/Wastewater_SARS-CoV-2"
sleep 5s
echo ""
echo ""

mkdir /Data/results/ExtractedLineages

for dir in $(ls -d */ | grep -v 'results')
do
    
    numberoffiles=$(ls ${dir}*.fastq.gz | wc -l)
    SKIP="FALSE"

    if (( ${numberoffiles} == 1 ))
    then
    FILESIZE=$(stat -c%s ${dir}/*.fastq.gz)

    if (( ${FILESIZE} < 1000000))
    then
    SKIP="TRUE"
    echo "Skipping "${dir}
    fi
    fi

    if [[ ${numberoffiles} > 1  || ${SKIP} == "FALSE" ]] 
    then
    echo ""
    echo "Processing "${numberoffiles}" fastq.gz files in "${dir}
    echo ""
      

    cd ${dir}
  	Reads=$(ls *.fastq.gz)
    cat ${Reads} > ${dir%/}.fastq.gz
    gzip -d ${dir%/}.fastq.gz
    seqkit seq ${dir%/}.fastq -M ${6} -m ${5} -Q ${1} > ${dir%/}.filtered.fastq
    rm ${dir%/}.fastq
    seqkit fq2fa ${dir%/}.filtered.fastq -o ${dir%/}.uncompressed.fasta
    
    cp /home/docker/CommonFiles/reference/reference.msh ./reference.msh
    cp /home/docker/CommonFiles/Tools/mash ./mash
    Rscript /home/docker/CommonFiles/Tools/MashXtractor.R ${1}
    cp *Xtracted* /Data/results/ExtractedLineages
    cd ..
done