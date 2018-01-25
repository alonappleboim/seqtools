#! /bin/bash
# execute after all samples were handled by handle_4tU_sample.bash
#
# Arguments:
#  $1 - final matlab file name and struct name
#  $2 - whether tmp files should be removed or not, default=1
#
# now need to:
# - collect 2 bigwig hubs (coverage and Ts)
# - collect sample statistics (# reads, aligned to normal genome, converted, unconverted, rejected, uniques)
# - collect duplicate histogrmas per sample
# - collect read T histograms into a matlab struct
# - collect genomic T histograms into a matlab struct
# - collect mutation rates into a matlab struct

rmv_tmp=${2:-1};

B2M="/cs/bd/tools/scripts/bwhub2mat.py"
MRG="/cs/bd/tools/4tU/merge_stats.py"
CL="/cs/wetlab/genomics/scer/genome/sacCer3_ordered.sizes"

#collect sample names
samples=();
while read line; do
  read -ra arr <<<"$line"
  samples+=(${arr[1]})
done < "/dev/stdin"

# generate bigwig hubs and mat files
# handle general coverage hub - a typical sample: well-A1_convert-1_rep-1_stress-0_tU-0_mod-cnv.bam
covmat="tmp/cov.mat";
cov_ve="well-\w+_convert-(?P<cnv>\d)_rep-(?P<rep>\d)_stress-(?P<stress>\d+)_tU-(?P<tU>\d+)_mod-(?P<mod>\w+)_str-(?P<str>w|c)\.bw"
cov_ord="mod:str;tU:num;cnv:num;stress:num;rep:num;str:str";
echo "generating coverage mat file --> $covmat"
python $B2M "BIGWIG/" "cov4tU" -ve $cov_ve -cm $CL -ob $cov_ord 2> /dev/null > $covmat

tmat="tmp/t.mat";
t_ve="well-\w+_convert-(?P<cnv>\d)_rep-(?P<rep>\d)_stress-(?P<stress>\d+)_tU-(?P<tU>\d+)_t-(?P<t>\w+)\.(?P<str>c|w)\.bw"
t_ord="t:str;tU:num;cnv:num;stress:num;rep:num;str:str";
echo "generating genomic t mat file --> $tmat"
python $B2M "BIGWIG/" "Ts4tU" -ve $t_ve -cm $CL -ob $t_ord 2> /dev/null > $tmat

# merge all into one matlab file
matf=$1".mat"
s_ord="tU:num;convert:num;stress:num;rep:num";
python $MRG -on $1 -ob $s_ord > $matf
cp $matf ~/Dropbox/

# generate track hubs (add to hub index manually):
# stress time course
rm -Rf /cs/bd/track_hubs/4tU-stress #delete if exists
python /cs/bd/tools/seqtools/future/build_hub.py BIGWIG/ 4tU-stress \
"well-\w+_convert-1_rep-(?P<rep>\d+)_stress-(?P<stress>\d+)_tU-10_mod-(?P<mod>cnv|unc|nrm)_str-(?P<str>w|c)\.bw" \
-g "mod,rep,str" -c "mod" -o "stress:num;mod:str;rep:num;str:str" -tp "str=c:negateValues=on"

# remove temporary files
echo -n "delete tmp folder (Y/n)? "
read answer
if echo "$answer" | grep -iq "^n" ;then
  echo "then remember to remove it yourself!"
else
  rm tmp -Rf
  echo "deleted!"
fi

