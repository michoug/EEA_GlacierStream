# Pipeline for metagenome analysis of EEA

Grégoire Michoud

November 2021

## Trimming and quality

### Installation

`mamba create -n trim_galore -c bioconda -c conda-forge trim-galore`

### Run

All paired `fastqc.gz` files should be placed in the `Raw` folder and should follow the convention `sampleName_R1.fastq.gz` and `sampleName_R2.fastq.gz`

``` bash
conda activate trim_galore
for i in Raw/*R1.fastq.gz;
do
    trim_galore -j 4 --fastqc --paired $i ${i%R1.fastq.gz}R2.fastq.gz;
done

mkdir Trim
mv Raw/*val_*.fastq.gz Trim
```

## Assembly

### Installation

`mamba create -n megahit -c bioconda -c conda-forge megahit`

### Script

``` bash
conda activate megahit
for i in Trim/*R1_val_1.fastq.gz;
do
    megahit -1 $i -2 ${i%R1_val_1.fastq.gz}R2_val_2.fastq.gz -m 2.4e+11 -t 28 --kmin-1pass --min-contig-len 1000 -o ${i%_mtg_sed_R1_val_1.fastq.gz}
done

mkdir Contigs
for i in */final.contigs.fa
do
    mv $i Contigs/${i%/final.contigs.fa}.fa
done
```

## Annotation

### Installation

`mamba create -n annotation -c bioconda -c conda-forge prodigal seqkit eggnog-mapper`

Had some issues with the conda installation of `eggnog-mapper`  so downloaded the release (https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/2.1.2.tar.gz)
and used this version inside the conda environment

### Script

``` bash
conda activate annotation

cd Contigs

# Rename contigs to avoid duplicates

for i in *fa
    do perl -pe "s/>/>${i%.fa}\_/g" $i > ${i%.fa}_rename.fa
done

# Split fasta files into 28 parts as prodigal is single threaded

for i in *_rename.fa
do
    seqkit split2 -p 28 -j 28 $i;
done

# Run in parallel as you wish to greatly speed up the process

for i in *split/*fa
    do prodigal -a $i\a -d ${i%.fa}.ffn -f gff -i $i -p meta -q -o ${i%.fa}.gff;
done

for i in *split
do
    cat $i\*.faa > ${i%.fa.split}.faa
    cat $i\*.ffn > ${i%.fa.split}.ffn
    cat $i\*.gff > ${i%.fa.split}.gff
done

cd ..

mkdir Annotations

mv Contigs/*faa Annotations
mv Contigs/*ffn Annotations
mv Contigs/*faa Annotations

# The --dbmem parameter should be removed if the memory of the computer used is a little low but its addition strongly reduce the run time

for i in Annotations/*faa
do
    /work/sber/Databases/eggnog-mapper-2.1.2/emapper.py --dbmem --resume --cpu 28 -i $i --itype proteins -m diamond --sensmode very-sensitive -o ${i%.faa}_egg.txt
done
```

## Coverage and filtering

The enzymes and genes of interest are :

- AG = 3.2.1.20
- BG = 3.2.1.21
- NAG = 3.2.1.14
- LAP = 3.4.11.1
- AP = 3.1.3.1
- recA = K03553

### Installation

`mamba create -n EEAdiverse -c bioconda -c conda-forge r-tidyverse samtools coverm cd-hit kraken2`

### Script

``` bash
conda activate R

mkdir EggNogAnalysis
cd EggNogAnalysis

cat ../Annotations/*_egg.txt > allEggNogAnnotations.txt

Rscript extractGenesOfInterest.R allEggNogAnnotations.txt EggNoglist.txt EggNogAnnotationsSubset.txt

# Concactenate all ffn files to facilitate the extraction of the genes of interest

cat ../Annotations/*ffn > temp.ffn

samtools faidx temp.ffn -r EggNoglist.txt -o EEA_genes.ffn

rm temp.ffn

# Cluster the genes with cd-hit to avoid duplicates and to facilitate the coverage calculation

cd-hit-est -i EEA_genes.ffn -o EEA_genes_cd -c 0.95 -T 28 -G 0 -aS 0.9 -g 1 -r 1 -d 0

# list reads files in the Trim folder

ls ../Trim > listreads.txt

# Calculate coverage of each of the genes of interest

coverm contig -r EEA_genes_cd.ffn -m trimmed_mean -o EEA_genes_cd_cov.txt -t 28 -c `< listreads.txt`
```

## Taxonomy and Final calculations

### Script

conda activate EEAdiverse

kraken2 --threads 28 --db maxikraken2_1903_140GB/ --use-names --confidence 0.5  --output EEA_genes_krakenOut.txt --report EEA_genes_krakenReport.txt EEA_genes_cd.ffn

