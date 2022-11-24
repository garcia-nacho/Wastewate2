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
echo "Quality:"${1}"/Noise Cutoff:"${2}
echo "Analyzing Spike gene from nt "${3}" to "${4}
echo "Read size between "${5}" and "${6} "nt"
echo "Mode "${7}
echo "Trimming "${8}"nt"
echo ${9}" Technology"
echo ""
echo -e "Quality, noise, read size, region to analyze and trimming can be set using these flags: \n -e qual=Q, -e noise=N, -e m=min, -e M=max -e start=S, -e end=E -e trim=20"
echo "Modes independent or dependent can be set with -e mode=i or -e mode=d"
echo "Stay tuned for updates!" 
echo "Visit us at https://github.com/folkehelseinstituttet/Wastewater_SARS-CoV-2"
sleep 5s
echo ""
echo ""

RefBowtie2=/home/docker/CommonFiles/reference/SpikeSars-CoV-2
RefSpike=/home/docker/CommonFiles/reference/SpikeRef.fa
Tools=/home/docker/CommonFiles/Tools

basedir=/Data

mkdir -p /home/docker/results/

source activate nextclade
echo "Downloading Nextclade database"
nextclade dataset get --name 'sars-cov-2' --output-dir '/home/docker/nc_sars-cov-2'
conda deactivate

echo ""
echo "Preparing Spike References"
Rscript /home/docker/CommonFiles/Tools/Spike_Extractor.R

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
    echo ""
    echo "Processing "${numberoffiles}" fastq.gz files in "${dir}
    echo ""
    
    

    cd ${dir}
	  Reads=$(ls *.fastq.gz)
    cat ${Reads} > ${dir%/}.fastq.gz
    gzip -d ${dir%/}.fastq.gz
    seqkit seq ${dir%/}.fastq -M ${6} -m ${5} -Q ${1} > ${dir%/}.filtered.fastq

    if (( ${8} != 0))
    then
    echo "Trimming "${8}"nt"
    source activate cutadaptenv
    cutadapt --cut ${8} -o ./trimmed1.fastq ./${dir%/}.filtered.fastq 
    rm ./${dir%/}.filtered.fastq
    cutadapt --cut -${8} -o ./trimmed2.fastq ./trimmed1.fastq
    rm ./trimmed1.fastq
    mv ./trimmed2.fastq ./${dir%/}.filtered.fastq
    conda deactivate
    fi

    rm ${dir%/}.fastq
    #(bowtie2 -p 8 -x ${RefBowtie2} -U ${dir%/}.fastq -S ${dir%/}.sam) 2> ${dir%/}_Bowtie2summary.txt
    #(tanoti -r ${RefSpike} -i ${dir%/}.fastq -o ${dir%/}.sam -u) 2> Tanoti_${dir%/}.summary.txt
    minimap2 --secondary=no -ax map-ont ${RefSpike} ${dir%/}.filtered.fastq > ${dir%/}.sam  
    samtools view -F 1024 -F 256 -F4 -F 2048 -bS ${dir%/}.sam | samtools sort -o ${dir%/}.sorted.bam
    samtools index ${dir%/}.sorted.bam

    minimap2 --secondary=no -ax map-ont /Data/VariantSpike.fasta ${dir%/}.filtered.fastq > ${dir%/}_voc.sam  
  
    Rscript /home/docker/CommonFiles/Tools/SamParser.R

    #samtools depth -a ${dir%/}_voc.sorted.bam > ${dir%/}_voc_depth.tsv
    #samtools mpileup -aa -A -d 0 -Q 0 --reference ${RefSpike} ./${dir%/}.sorted.bam | ivar consensus -t 0 -n N -m 20 -p ${dir%/}_consensus 
    ls
    samtools mpileup -aa -A -d 0 -q 0 --reference ${RefSpike} ./${dir%/}.sorted.bam | ivar consensus -t 0 -n N -m 20 -p ${dir%/}_consensus

    rm ${dir%/}.sam
    rm ${dir%/}_voc.sam
    #rm ${dir%/}_voc.sorted.bam
    #rm ${dir%/}_voc.sorted.bam.bai
      
    rm ${dir%/}.filtered.fastq
    ${Tools}/FINex2 -f ${dir%/}.sorted.bam > ${dir%/}.noise.tsv 
    cp ${Tools}/bbasereaderHC ./bbasereader
    cat ${dir%/}_consensus.fa ${RefSpike} > spike.cons.fa
    
    source activate nextclade
    nextclade --input-fasta spike.cons.fa --input-dataset /home/docker/nc_sars-cov-2 --output-csv dummy.csv --output-fasta spike.cons.aligned.fa
    conda deactivate
    #nextalign -i spike.cons.fa -o spike.cons.aligned.fa -r /home/docker/CommonFiles/reference/SpikeRef.fa

    if [ ${7} == d ]
    then
      Rscript ${Tools}/AnalysisWW.R $(pwd)/ ${2} ${3} ${4}
      mkdir /home/docker/results/${dir%/}
      mv Variants.fa /home/docker/results/${dir%/}/${dir%/}.variants.fa
      mv VariantResults.xlsx /home/docker/results/${dir%/}/${dir%/}.results.xlsx
      mv Coverage.pdf /home/docker/results/${dir%/}/${dir%/}.coverage.pdf
      mv Noise.pdf /home/docker/results/${dir%/}/${dir%/}.noise.pdf
      mv ${dir%/}.noise.tsv  /home/docker/results/${dir%/}/${dir%/}.noise.tsv 
      mv ${dir%/}.sorted.bam /home/docker/results/${dir%/}/${dir%/}.sorted.bam
      mv ${dir%/}.sorted.bam.bai /home/docker/results/${dir%/}/${dir%/}.sorted.bam.bai
      mv ${dir%/}_voc_depth.tsv /home/docker/results/${dir%/}/${dir%/}_voc_depth.tsv     
      mv *_consensus.qual.txt /home/docker/results/${dir%/}/${dir%/}_consensus.qual.txt 
      mv ${dir%/}_consensus.fa /home/docker/results/${dir%/}/${dir%/}_consensus.fa
      mv spike.cons.aligned.fa /home/docker/results/${dir%/}/spike.cons.aligned.fa
      mv ${dir%/}_voc_depth.csv /home/docker/results/${dir%/}_voc_depth.csv

    else
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
      mv ${dir%/}_voc_depth.tsv /home/docker/results/${dir%/}_voc_depth.tsv 
    fi

    rm dummy.csv
    rm bbasereader
    rm Rplots.pdf
    rm spike.cons*
    cd ${basedir}  
    fi
    #rm *.fasta
done

cp -R /home/docker/results /Data/results

if [ ${7} == d ]
then
  cd /Data/results
  Rscript ${Tools}/POIgenerator.R $(pwd)/ ${2} ${3} ${4}
  for dir2 in $(ls -d */)
  do
    echo "Dependent mode on" ${dir2}
    cd ${dir2}
    cp ${Tools}/bbasereaderHC ./bbasereader
    cp ../poi.temp ./poi.temp 

    Rscript ${Tools}/AnalysisWWDep.R $(pwd)/ ${2} ${3} ${4}

    mv ${dir2%/}.variants.fa /Data/results/${dir2%/}.variants.fa
    mv *pdf /Data/results/
    mv ${dir2%/}.results.xlsx /Data/results/${dir2%/}.results.xlsx
    mv ${dir2%/}.noise.tsv  /Data/results/${dir2%/}.noise.tsv 
    mv ${dir2%/}.sorted.bam /Data/results/${dir2%/}.sorted.bam
    mv ${dir2%/}.sorted.bam.bai /Data/results/${dir2%/}.sorted.bam.bai
    mv *_consensus.qual.txt /Data/results/${dir2%/}_consensus.qual.txt 
    mv ${dir2%/}_consensus.fa /Data/results/${dir2%/}_consensus.fa
    rm spike.cons.aligned.fa
    mv VariantsDependent.fa /Data/results/${dir2%/}.variants.dependent.fa
    mv VariantResultsDependent.xlsx /Data/results/${dir2%/}.results.dependent.xlsx
    rm poi.temp
    rm bbasereader
    cd ..
    rm -rf ${dir2}
  done
fi

rm poi.temp

if [ ${7} == d ]
then
cd /Data/results
mkdir sequences
cat *.variants.dependent.fa > sequences/merged_variants.fa

else
cd /Data/results
mkdir sequences
cat *.variants.fa > sequences/merged_variants.fa
fi

source activate nextclade
nextclade --input-fasta sequences/merged_variants.fa --input-dataset /home/docker/nc_sars-cov-2 --output-csv Nextclade.results.csv --output-fasta merged_variants.aligned.fa
conda deactivate
Rscript ${Tools}/postanalysisWW.R


#post analysis Fixed
mkdir bam analysis QC
 
mkdir analysis/FixedMode
mv /Data/results/*.bam /Data/results/bam
mv /Data/results/*.bai /Data/results/bam
mv /Data/results/*.fa /Data/results/sequences
mv /Data/results/*.qual.txt /Data/results/QC
mv /Data/results/*noise.tsv /Data/results/QC
mv /Data/results/*coverage.pdf /Data/results/QC
mv /Data/results/*noise.pdf /Data/results/QC
mv /Data/results/*.xlsx /Data/results/analysis
mv /Data/results/*Barplot.pdf /Data/results/analysis
mv /Data/results/*Sankeyplot*.pdf /Data/results/analysis
mv /Data/results/Results.Aggregated.xlsx /Data/results/analysis/SurveillenceMode
mv /Data/results/Nextclade.results.csv /Data/results/analysis/Variants.nextclade.csv
rm /Data/results/Rplots.pdf
rm /Data/results/merged_variants*.fasta
rm /Data/results/merged_variants*.csv
#To be changed after adding fixed mode


rm -rf /Data/results/analysis/SurveillenceMode
rm -rf /Data/results/analysis/FixedMode

cd /Data/results
cp ${Tools}/bbasereaderHC ./bbasereader
Rscript /home/docker/CommonFiles/Tools/AnalysisSingleMuts.R $(pwd)/ ${2}
Rscript /home/docker/CommonFiles/Tools/NoiseMosaic.R $(pwd)/ ${2}
rm Rplots.pdf
Rscript /home/docker/CommonFiles/Tools/MappedVariantPlotter.R
mv /Data/results/*.pdf /Data/results/analysis
mv /Data/results/*.xlsx /Data/results/analysis
rm ./bbasereader
mv /Data/VariantSpike.fasta /Data/results/analysis/ReferencesSpike.fasta
mv /Data/results/*_voc_depth.csv /Data/results/QC