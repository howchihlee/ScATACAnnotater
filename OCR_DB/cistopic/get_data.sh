## please download files following the instruction in readme before runing this script
for f in CD34-Bone-Marrow CD34-Cord-Blood CD4Tcell CD8Tcell Bcell Mono NKcell
do
fn0=$f._peaks.narrowPeak
fn1=$f.bed
cat $fn0 | awk -F '\t' '{ print $1"\t"$2"\t"$3 }' | sort -k1,1 -k2,2n > $fn1
gzip $fn1
done

zcat *.gz | sort -k1,1 -k2,2n > tmp.bed
bedtools merge -i tmp.bed > background.bed
gzip background.bed
rm tmp.bed