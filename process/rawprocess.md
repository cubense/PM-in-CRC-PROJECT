```{cellranger,velocyto}
cellranger count --id=crc1output \
 --transcriptome=../cellranger/refdata-gex-GRCh38-2020-A \
 --fastqs=../pm/crc \
 --sample=samople-01,samople-02,samople-03,samople-04 \
 --expect-cells=20000 \
 --localcores=16 \
 --localmem=64

velocyto run10x -m ../velocyto/GRCh38_rmsk.gtf ../singlecellrawdata/HT2020-17780-1/crc6output ../cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf

```
```{CNMF}
python ../cNMF/cnmf.py consensus --output-dir ../cnmf/mc/ --name pm_crc1 --components 28 --local-density-threshold 2 --local-neighborhood-size 0.2 --show-clustering
```

```{wes:gatk4,cnvkit,vep}
cat sample.txt | while read id
do
 arr=(${id})
 echo "start bqsr for ${id}" `date`
bwa mem -M -t 16 -R "@RG\tID:${arr[1]}\tSM:${arr[1]}\tLB:WES\tPL:Illumina" /hg38/Homo_sapiens_assembly38.fasta ./${arr[0]}_1_val_1.fq.gz ./${arr[0]}_2_val_2.fq.gz | samtools sort -@ 10 -o /wes/bqsr/${arr[1]}.bam -
gatk --java-options -Xmx20G MarkDuplicates \
  -I ../wes/bqsr/${arr[1]}.bam \
  --REMOVE_DUPLICATES true \
  -O ../wes/bqsr/${arr[1]}_marked.bam \
  -M ../wes/bqsr/${arr[1]}.metrics

gatk --java-options -Xmx20G FixMateInformation \
  -I ../wes/bqsr/${arr[1]}_marked.bam \
  -O ../wes/bqsr/${arr[1]}_fixed.bam \
  -ADD_MATE_CIGAR true
samtools index ../wes/bqsr/${arr[1]}_fixed.bam
gatk --java-options -Xmx20G BaseRecalibrator \
    -R $ref  \
    -I ../wes/bqsr/${arr[1]}_fixed.bam \
    --known-sites $snp \
    --known-sites $indel \
    -O ../wes/bqsr/${arr[1]}_recal.table
gatk --java-options -Xmx20G ApplyBQSR \
    -R $ref  \
    -I ../wes/bqsr/${arr[1]}_fixed.bam \
    -bqsr ../wes/bqsr/${arr[1]}_recal.table \
    -O ../wes/bqsr2/${arr[1]}_bqsr.bam
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" Mutect2 -R ${ref} \
 -I ${T} -tumor $(basename "$T" _bqsr.bam) \
 -I ${N} -normal $(basename "$N" _bqsr.bam) \
 --germline-resource $gnomad \
 --panel-of-normals ../wes/bqsr2/pon/pon.vcf.gz \
 -L ${bed} \
 -O ../wes/bqsr2/vcf/${sample}_mutect2.vcf

 gatk FilterMutectCalls \
  -R ${ref} \
  -V ../wes/bqsr2/vcf/${sample}_mutect2.vcf \
  -O ..wes/bqsr2/vcf/${sample}_somatic.vcf

cat pair.txt |while read id
do
 arr=(${id})
 vcf2maf.pl --input-vcf ${arr[1]}_vep_tmp.vcf \
 --output-maf ../wes/bqsr2/maf/${arr[1]}_vep.maf \
 --ref-fasta ../WES/hg38/Homo_sapiens_assembly38.fasta \
 --vep-data ../vep/vep \
 --vep-path ../envs/bin/ \
 --tumor-id ${arr[1]} \
 --normal-id ${arr[0]} \
 --ncbi-build GRCh38 
done
cnvkit.py batch *C_bqsr.bam *M_bqsr.bam --normal *BD_bqsr.bam \
    --access ../access-5kb-mappable.hg38.bed \
    --targets $bed --fasta $ref \
    --drop-low-coverage --scatter --diagram -p 10 \
    --output-reference my_reference.cnn --output-dir ./cnv/
cnvkit.py heatmap *CRC_bqsr.cns *PM_bqsr.cns *LM_bqsr.cns -o cm.pdf
```

```{qiime2}
docker run -t -i -v $(pwd):/data quay.io/qiime2/amplicon:2023.9 /bin/bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --input-format PairedEndFastqManifestPhred33V2 --output-path ./16s/pc.qza
qiime demux summarize --i-data pc.qza --o-visualization pc.qzv
qiime dada2 denoise-paired --i-demultiplexed-seqs pc.qza --p-trim-left-f 10 --p-trim-left-r 20 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 20 --o-table table-raw0.qza --o-representative-sequences rep-raw0.qza --o-denoising-stats denoising-raw0.qza
qiime feature-table summarize --i-table table-raw0.qza --o-visualization table-raw0.qzv --m-sample-metadata-file metadata.txt
qiime feature-table tabulate-seqs --i-data rep-raw0.qza --o-visualization rep-raw0.qzv
qiime metadata tabulate --m-input-file denoising-raw0.qza --o-visualization denoising-raw0.qzv
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-nb-classifier.qza --i-reads rep-raw0.qza --o-classification taxonomy-raw0.qza
qiime metadata tabulate --m-input-file taxonomy-raw0.qza --o-visualization taxonomy-raw0.qzv
qiime taxa barplot --i-table table-raw0.qza --i-taxonomy taxonomy-raw0.qza --m-metadata-file metadata.txt --o-visualization taxa-bar-plots.qzv

#filter
qiime feature-table filter-features --i-table table-raw0.qza --p-min-frequency 10 --o-filtered-table table_f1.qza
qiime taxa filter-table --i-table table_f1.qza --i-taxonomy taxonomy-raw0.qza --p-exclude mitochondria,chloroplast --o-filtered-table table_f2.qza
qiime feature-table filter-features --i-table table_f2.qza --p-min-samples 2 --o-filtered-table table_f3.qza
qiime taxa filter-table --i-table table_f3.qza --i-taxonomy taxonomy-raw0.qza --p-include p_ --o-filtered-table table_f4.qza
qiime taxa barplot --i-table table_f4.qza --i-taxonomy taxonomy_f.qza --m-metadata-file metadata.txt --o-visualization taxa-bar-plots_f.qzv

qiime feature-table summarize --i-table table_f4.qza --o-visualization table_f4.qzv --m-sample-metadata-file metadata.txt

qiime taxa filter-table --i-table table_f2.qza --i-taxonomy taxonomy-raw0.qza --p-include p_ --o-filtered-table table_f5.qza
qiime feature-table summarize --i-table table_f5.qza --o-visualization table_f5.qzv --m-sample-metadata-file metadata.txt

qiime feature-table filter-seqs --i-data rep-raw0.qza --i-table table_f4.qza --o-filtered-data rep_f4.qza
qiime feature-table tabulate-seqs --i-data rep_f4.qza --o-visualization rep_f4.qzv

qiime feature-classifier classify-sklearn --i-classifier silva-138-99-nb-classifier.qza --i-reads rep_f4.qza --o-classification taxonomy_f4.qza
qiime metadata tabulate --m-input-file taxonomy_f4.qza --o-visualization taxonomy_f4.qzv
qiime taxa barplot --i-table table_f4.qza --i-taxonomy taxonomy_f4.qza --m-metadata-file metadata.txt --o-visualization taxa-bar-plots_f4.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep_f4.qza --o-alignment aligned-rep_f4.qza --o-masked-alignment masked-aligned-rep_f4.qza --o-tree unrooted-tree_f4.qza --o-rooted-tree rooted-tree_f4.qza

qiime tools export --input-path table_f4.qza --output-path exported
qiime tools export --input-path rep_f4.qza --output-path exported
qiime tools export --input-path unrooted-tree_f4.qza --output-path exported
qiime tools export --input-path rooted-tree_f4.qza --output-path exported
qiime tools export --input-path taxonomy_f4.qza --output-path exported

biom convert -i ./exported/feature-table.biom -o ./exported/table.txt --to-tsv
```
 
