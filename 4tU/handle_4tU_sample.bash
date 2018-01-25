#! /bin/bash

# inputs:
# $1 - sample name
# $2 - annotation windows, a stranded bed file with annotations to count reads in
# $3 - remove reads aligned to this spike-in genome, if given

IDX="/cs/wetlab/genomics/scer/bowtie/sacCer3";
TLG="/cs/wetlab/genomics/scer/bowtie/only-ACG-sacCer3";
CL="/cs/wetlab/genomics/scer/genome/sacCer3_ordered.sizes";
W2BW="/cs/bd/tools/wigToBigWig";
B2BW="/cs/bd/tools/bedGraphToBigWig";

RMD_3P="/cs/bd/tools/4tU_scripts/3p_rm_dups.awk";
RMD_3PH="/cs/bd/tools/4tU_scripts/rm_dups_with-hdr-seq.awk";
F2A="/cs/bd/tools/4tU_scripts/fq2acg_fq.awk";
AS2S="/cs/bd/tools/4tU_scripts/acg_sam2sam.awk";
PTS="/cs/bd/tools/4tU_scripts/parse_tlg_sam.awk";
FTB="/cs/bd/tools/4tU_scripts/fit_truncated_binomial.py";
#FTB="/cs/phd/alonap/Dropbox/workspace/seqtools/4tU/fit_truncated_binomial.py";

annots=${2:-0};
SPK_IDX=${3:-0};
s=$1;
declare -A bam_list;
bed_list=();
wig_list=();
declare -A stts_files;

echo "handling $s"
echo "  building fastq..."
fq="fastq/"$s".fastq.gz";
cat tmp/hd*/$s* | awk -F "\t" '{print "@umi:"substr($4,8,8)"\n"$3"\n+\n"$7}' | gzip > "$fq"

if [ "$SPK_IDX" != "0" ]; then
  echo "count and remove spike-in reads..."
  spk_align="STATS/"$s".spk.align.stats";
  prelim_align="STATS/"$s".prelim.align.stats";
  tmp_bam="tmp/"$s".tmp.bam";
  stts_files["spk-align"]=$spk_align;
  stts_files["prelim-align"]=$prelim_align;
  fq_nspk="tmp/"$s".no.spk.fastq.gz"
  bowtie2 -p 8 -U $fq -x $SPK_IDX 2> $prelim_align > $tmp_bam
  # pass spike-in unaligned reads(f=4) to the pipeline
  samtools view $tmp_bam -b | samtools view -f 4 -b | samtools fastq  - 2> /dev/null | gzip > $fq_nspk
  # count spike-in alignd reads (f!=4) that do not align (f=4) to regular genome, since they may (probably) come from cerevisae..
  samtools view $tmp_bam -b | samtools view -F 4 -b | samtools fastq  - 2> /dev/null |\
  bowtie2 -p 8 -U - -x $IDX 2> $spk_align > /dev/null
else
  fq_nspk=$fq
fi

echo "  align reads to normal genome..."
bam_nrm="BAM/"$s"_mod-nrm.bam";
bam_list["nrm"]=$bam_nrm;
nrm_dhist="STATS/"$s".nrm.dhist";
nrm_align="STATS/"$s".nrm.align.stats";
stts_files["nrm-align"]=$nrm_align;
stts_files["nrm-dhist"]=$nrm_dhist;

zcat $fq_nspk | bowtie2 -p 8 -U - -x $IDX 2> $nrm_align |\
samtools sort | samtools view | awk -f $RMD_3P -v dhist_file=$nrm_dhist |\
samtools view -t $CL -h -b > $bam_nrm
samtools index $bam_nrm

echo "  align to 3 letter genome and split to rejected/unconverted/converted BAMs..."
bam_cnv="BAM/"$s"_mod-cnv.bam"; sam_cnv="tmp/"$s"_cnv.sam";
bam_unc="BAM/"$s"_mod-unc.bam"; sam_unc="tmp/"$s"_unc.sam";
bam_rej="BAM/"$s"_mod-rej.bam"; sam_rej="tmp/"$s"_rej.sam";
bam_list["cnv"]=$bam_cnv;
bam_list["unc"]=$bam_unc;
bam_list["rej"]=$bam_rej;

ow="tmp/"$s"_t-obs";
cw="tmp/"$s"_t-cnv";
wig_list+=($ow".c.wig")
wig_list+=($ow".w.wig")
wig_list+=($cw".c.wig")
wig_list+=($cw".w.wig")

bam_nrm="BAM/"$s"_mod-nrm.bam";
mutf="STATS/"$s".mut";
rt="STATS/"$s".read.t.hist";
gt="STATS/"$s".genomic.t.hist";
cl="STATS/"$s".cnv.stats";
ah="STATS/"$s".alen.hist";
sc="STATS/"$s".strippedc.cnt";
stts_files["mut"]=$mutf;
stts_files["gt"]=$gt;
stts_files["rt"]=$rt;
stts_files["cl"]=$cl;
stts_files["ah"]=$ah;
stts_files["sc"]=$sc;

tlg_dhist="STATS/"$s".tlg.dhist";
tlg_align="STATS/"$s".tlg.align.stats";
stts_files["tlg-align"]=$tlg_align;
stts_files["tlg-dhist"]=$tlg_dhist;

zcat $fq_nspk | paste - - - - | awk -F "\t" -f $F2A  -v stripc_file=$sc -v ahist_file=$ah |\
bowtie2 -p 8 -U - -x $TLG 2> $tlg_align | samtools sort | samtools view |\
awk -f $AS2S | awk -f $RMD_3PH -v dhist_file=$tlg_dhist |\
awk -f $PTS -v rej_sam=$sam_rej -v cnv_sam=$sam_cnv -v unc_sam=$sam_unc \
    -v to_wig=$ow -v tc_wig=$cw -v mut_tab=$mutf -v rt_hist=$rt -v gt_hist=$gt -v cl_stats=$cl;

samtools view -t $CL -h $sam_rej -b | samtools sort | samtools view -h | grep -v "^$" | samtools view -b > $bam_rej
samtools view -t $CL -h $sam_unc -b | samtools sort | samtools view -h | grep -v "^$" | samtools view -b > $bam_unc
samtools view -t $CL -h $sam_cnv -b | samtools sort | samtools view -h | grep -v "^$" | samtools view -b > $bam_cnv

samtools index $bam_rej;
samtools index $bam_unc;
samtools index $bam_cnv;

# fit data to (multiple) truncated binomial(s) to estimate conversion efficiency
tfit="STATS/"$s".p-cnv.fit";
echo "  fit a truncated binomial to estimate conversion probability -> $tfit"
stts_files["tbino_fit"]=$tfit;
python $FTB $rt -o $tfit

#if given, count reads in annotation windows
if [ $annots != "0" ]; then
  acnt="STATS/"$s".annot.cnt";
  echo "  count reads in annotation windows -> $acnt"
  stts_files["annot_cnt"]=$acnt;
  echo "annot	rej	unc	cnv" > $acnt
  bedtools multicov -bams $bam_rej $bam_unc $bam_cnv -bed $annots -S |\
  awk 'BEGIN {OFS="\t"} {if ($6+$7+$8>0) print NR,$6,$7,$8;}' >> $acnt
fi

echo "  generating BIGWIGs (coverage/T tracks)..."
# generating tmp bed files
for modtype in "${!bam_list[@]}"
do
  cbed="tmp/"$s"_mod-"$modtype"_str-c.bed.tmp";
  wbed="tmp/"$s"_mod-"$modtype"_str-w.bed.tmp";
  genomeCoverageBed -ibam ${bam_list[$modtype]} -g $CL -dz -strand + > $wbed
  genomeCoverageBed -ibam ${bam_list[$modtype]} -g $CL -dz -strand - > $cbed
  bed_list+=($cbed);
  bed_list+=($wbed);
done

# generating coverage bigwig files
for bedf in "${bed_list[@]}"; do
  n=`basename $bedf`
  n=${n%.bed.tmp}
  tmpf="tmp/"$n".tmp"; #loop temp file
  bwf="BIGWIG/"$n".bw"; # final bigwig
  cat $bedf | awk '{print $1,$2,$2+1,$3;}' | sort -k1,1 -k2,2n > $tmpf
  $B2BW $tmpf $CL $bwf
  rm $tmpf
done

# generating T count bigwigs
for wigf in "${wig_list[@]}"; do
  n=`basename $wigf`
  n=${n%.wig}
  bwf="BIGWIG/"$n".bw"; # final bigwig
  $W2BW $wigf $CL $bwf
done


#make a statatistics files index for easy merging:
index="STATS/"$s".index";
rm $index 2> /dev/null
touch $index
for stat in "${!stts_files[@]}"
do
  echo "$stat:${stts_files[$stat]}" >> $index
done

# all done! now need to:
# - collect 2 bigwig hubs (coverage and Ts)
# - collect sample statistics (# reads, aligned to normal genome, converted, unconverted, rejected, uniques)
# - collect duplicate histogrmas per sample
# - collect read T histograms into a matlab struct
# - collect genomic T histograms into a matlab struct
# - collect mutation rates into a matlab struct

