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

echo "Quality:"${1}"/Noise Cutoff:"${2}
echo "Analyzing Spike gene from "${3}" to "${4}
echo "Read size between"${5}" and "${6}
echo -e "Quality, noise, read size and region to analyze can be set using these flags: \n -e qual=Q, -e noise=N, -e m=min, -e M=max -e start=S, -e end=E"
sleep 3s

RefBowtie2=/home/docker/CommonFiles/reference/SpikeSars-CoV-2
RefSpike=/home/docker/CommonFiles/reference/SpikeRef.fa
Tools=/home/docker/CommonFiles/Tools

basedir=/Data

mkdir -p /home/docker/results/

source activate nextclade
nextclade dataset get --name 'sars-cov-2' --output-dir '/home/docker/nc_sars-cov-2'
conda deactivate

for dir in $(ls -d */)
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
    echo "Processing "${numberoffiles}"fastq.gz files in"${dir}

    cd ${dir}
	  Reads=$(ls *.fastq.gz)
    cat ${Reads} > ${dir%/}.fastq.gz
    gzip -d ${dir%/}.fastq.gz
    seqkit seq ${dir%/}.fastq -M ${6} -m ${5} -Q ${1} > ${dir%/}.filtered.fastq
    rm ${dir%/}.fastq
    #(bowtie2 -p 8 -x ${RefBowtie2} -U ${dir%/}.fastq -S ${dir%/}.sam) 2> ${dir%/}_Bowtie2summary.txt
    #(tanoti -r ${RefSpike} -i ${dir%/}.fastq -o ${dir%/}.sam -u) 2> Tanoti_${dir%/}.summary.txt
    minimap2 --secondary=no -ax map-ont ${RefSpike} ${dir%/}.filtered.fastq > ${dir%/}.sam  
    samtools view -F 1024 -F 256 -F4 -F 2048 -bS ${dir%/}.sam | samtools sort -o ${dir%/}.sorted.bam
    samtools index ${dir%/}.sorted.bam
    samtools mpileup -aa -A -d 0 -Q 0 --reference ${RefSpike} ${dir%/}.sorted.bam | ivar consensus -t 0 -n N -m 20 -p ${dir%/}_consensus 
    rm ${dir%/}.sam
    rm ${dir%/}.filtered.fastq
    ${Tools}/FINex2 -f ${dir%/}.sorted.bam > ${dir%/}.noise.tsv 
    cp ${Tools}/bbasereaderHC ./bbasereader
    cat ${dir%/}_consensus.fa ${RefSpike} > spike.cons.fa
    
    source activate nextclade
    nextclade --input-fasta spike.cons.fa --input-dataset /home/docker/nc_sars-cov-2 --output-csv dummy.csv --output-fasta spike.cons.aligned.fa
    conda deactivate

    Rscript ${Tools}/AnalysisWW.R $(pwd)/ ${2} ${3} ${4}
    #Rscript ${Tools}/AnalysisFixedWW.R $(pwd)/
    mv Variants.fa /home/docker/results/${dir%/}.variants.fa
    mv VariantResults.xlsx /home/docker/results/${dir%/}.results.xlsx
    mv Coverage.pdf /home/docker/results/${dir%/}.coverage.pdf
    mv Noise.pdf /home/docker/results/${dir%/}.noise.pdf
    mv ${dir%/}.noise.tsv  /home/docker/results/${dir%/}.noise.tsv 
    mv ${dir%/}.sorted.bam /home/docker/results/${dir%/}.sorted.bam
    mv ${dir%/}.sorted.bam.bai /home/docker/results/${dir%/}.sorted.bam.bai
    mv *_consensus.qual.txt /home/docker/results/${dir%/}_consensus.qual.txt 
    mv ${dir%/}_consensus.fa /home/docker/results/${dir%/}_consensus.fa
    rm dummy.csv
    rm bbasereader
    rm Rplots.pdf
    rm spike.cons*
    cd ${basedir}  
    fi
    #rm *.fasta
done

cp -R /home/docker/results /Data/results

#post analysis Surveillence
cd /Data/results
mkdir sequences
cat *.variants.fa > sequences/merged_variants.fa
source activate nextclade
nextclade --input-fasta sequences/merged_variants.fa --input-dataset /home/docker/nc_sars-cov-2 --output-csv Nextclade.results.csv --output-fasta merged_variants.aligned.fa
conda deactivate
Rscript ${Tools}/postanalysisWW.R
rm merged_variants*

#post analysis Fixed
mkdir bam analysis QC
mkdir analysis/SurveillenceMode
mkdir analysis/FixedMode
mv /Data/results/*.bam /Data/results/bam
mv /Data/results/*.bai /Data/results/bam
mv /Data/results/*.fa /Data/results/sequences
mv /Data/results/*.qual.txt /Data/results/QC
mv /Data/results/*noise.tsv /Data/results/QC
mv /Data/results/*coverage.pdf /Data/results/QC
mv /Data/results/*noise.pdf /Data/results/QC
mv /Data/results/*results.xlsx /Data/results/analysis/SurveillenceMode
mv /Data/results/*Barplot.pdf /Data/results/analysis/SurveillenceMode
mv /Data/results/*Sankeyplot*.pdf /Data/results/analysis/SurveillenceMode
mv /Data/results/Results.Aggregated.xlsx /Data/results/analysis/SurveillenceMode
mv /Data/results/Nextclade.results.csv /Data/results/analysis/SurveillenceMode/Variants.nextclade.csv
rm /Data/results/Rplots.pdf
#To be changed after adding fixed mode
mv /Data/results/analysis/SurveillenceMode/* /Data/results/analysis/
rm -rf /Data/results/analysis/SurveillenceMode
rm -rf /Data/results/analysis/FixedMode
#mv *ResultsFixedPos.xlsx /analysis/FixedMode



