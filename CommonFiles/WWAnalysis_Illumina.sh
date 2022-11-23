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
echo "Noise Cutoff:"${2}
echo "Analyzing Spike entire Spike gene"

echo "Trimming "${8}"nt"
echo ${9}" Technology"
echo ""
echo -e "Noise and trimming can be set using these flags: \n  -e noise=N -e trim=20"
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
Rscript SpikeExtractor.R

for dir in $(ls -d */)
do
    
    numberoffiles=$(ls ${dir}*.fastq.gz | wc -l)

    cd ${dir}
	Reads1=$(ls *_R1_*.fastq.gz)
    Reads2=$(ls *_R2_*.fastq.gz)
    cat ${Reads1} > ${dir%/}.R1.fastq.gz
    cat ${Reads2} > ${dir%/}.R2.fastq.gz

    if (( ${8} != 0))
    then
    echo "Trimming "${8}"nt"
    source activate cutadaptenv
    cutadapt --cut ${8} -o ./${dir%/}.R1.filtered.fastq.gz ${dir%/}.R1.fastq.gz  
    cutadapt --cut -${8} -o ./${dir%/}.R2.filtered.fastq.gz ${dir%/}.R2.fastq.gz
    rm ${dir%/}.R1.fastq.gz
    rm ${dir%/}.R2.fastq.gz
    
    conda deactivate
    else
    mv ${dir%/}.R1.fastq.gz ${dir%/}.R1.filtered.fastq.gz
    mv ${dir%/}.R2.fastq.gz ${dir%/}.R2.filtered.fastq.gz
    fi

    rm ${dir%/}.fastq
    (bowtie2 -p 8 -x ${RefBowtie2} -1 ${dir%/}.R1.filtered.fastq.gz -2 ${dir%/}.R1.filtered.fastq.gz -S ${dir%/}.sam) 2> ${dir%/}_Bowtie2summary.txt
    #(tanoti -r ${RefSpike} -i ${dir%/}.fastq -o ${dir%/}.sam -u) 2> Tanoti_${dir%/}.summary.txt
    minimap2 --secondary=no -ax map-ont /Data/VariantSpike.fasta ${dir%/}.filtered.fastq > ${dir%/}_voc.sam  
    samtools view  -F 256 -F4 -F 2048 -bS ${dir%/}.sam | samtools sort -o ${dir%/}.sorted.bam
    samtools index ${dir%/}.sorted.bam

    samtools view -F 256 -F4 -F 2048 -bS ${dir%/}_voc.sam | samtools sort -o ${dir%/}_voc.sorted.bam
    samtools index ${dir%/}_voc.sorted.bam

    samtools depth -a ${dir%/}_voc.sorted.bam > ${dir%/}_voc_depth.tsv


    #samtools mpileup -aa -A -d 0 -Q 0 --reference ${RefSpike} ./${dir%/}.sorted.bam | ivar consensus -t 0 -n N -m 20 -p ${dir%/}_consensus 
    ls
    samtools mpileup -aa -A -d 0 -q 0 --reference ${RefSpike} ./${dir%/}.sorted.bam | ivar consensus -t 0 -n N -m 20 -p ${dir%/}_consensus

    rm ${dir%/}.sam
    rm ${dir%/}_voc.sam

    rm *.filtered.fastq.gz
    ${Tools}/FINex2 -f ${dir%/}.sorted.bam > ${dir%/}.noise.tsv 
    cp ${Tools}/bbasereaderHC ./bbasereader
    cat ${dir%/}_consensus.fa ${RefSpike} > spike.cons.fa
    
    source activate nextclade
    nextclade --input-fasta spike.cons.fa --input-dataset /home/docker/nc_sars-cov-2 --output-csv dummy.csv --output-fasta spike.cons.aligned.fa
    conda deactivate
    #nextalign -i spike.cons.fa -o spike.cons.aligned.fa -r /home/docker/CommonFiles/reference/SpikeRef.fa

    Rscript ${Tools}/AnalysisWW_Illumina.R $(pwd)/ ${2}
    #Rscript ${Tools}/AnalysisFixedWW.R $(pwd)/
    mv Coverage.pdf /home/docker/results/${dir%/}.coverage.pdf
    mv Noise.pdf /home/docker/results/${dir%/}.noise.pdf
    mv ${dir%/}.noise.tsv  /home/docker/results/${dir%/}.noise.tsv 
    mv ${dir%/}.sorted.bam /home/docker/results/${dir%/}.sorted.bam
    mv ${dir%/}_voc.sorted.bam /home/docker/results/${dir%/}_voc.sorted.bam
    mv ${dir%/}.sorted.bam.bai /home/docker/results/${dir%/}.sorted.bam.bai
    mv ${dir%/}_voc.sorted.bam.bai /home/docker/results/${dir%/}_voc.sorted.bam.bai
    mv ${dir%/}_voc_depth.tsv /home/docker/results/${dir%/}_voc_depth.tsv
    mv *_consensus.qual.txt /home/docker/results/${dir%/}_consensus.qual.txt 
    mv ${dir%/}_consensus.fa /home/docker/results/${dir%/}_consensus.fa
    mv *Bowtie2summary.txt /home/docker/results/${dir%/}_bowtie2_summary.txt
    rm dummy.csv
    rm Rplots.pdf
    rm spike.cons*
    rm bbasereader
    cd ${basedir}  
    
 
    #rm *.fasta
done

cp -R /home/docker/results /Data/results

cd /Data/results
#post analysis Fixed
mkdir bam analysis QC sequences
 
mkdir analysis/FixedMode
mv /Data/results/*.bam /Data/results/bam
mv /Data/results/*.bai /Data/results/bam
mv /Data/results/*.fa /Data/results/sequences
mv /Data/results/*.qual.txt /Data/results/QC
mv /Data/results/*noise.tsv /Data/results/QC
mv /Data/results/*bowtie2_summary.txt /Data/results/QC
mv /Data/results/*coverage.pdf /Data/results/QC
mv /Data/results/*noise.pdf /Data/results/QC
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
Rscript /home/docker/CommonFiles/Tools/MappedVariantPlotter.R
mv /Data/results/*.pdf /Data/results/analysis
mv /Data/results/*.xlsx /Data/results/analysis
rm ./bbasereader
