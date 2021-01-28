wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147672/suppl/GSE147672_scATAC_idr_peaks.tar.gz
tar xvfz GSE147672_scATAC_idr_peaks.tar.gz
rm GSE147672_scATAC_idr_peaks.tar.gz

for i in {1..24}
do
fn0=Cluster$i.idr.optimal.narrowPeak.gz
fn1=Cluster$i.bed
zcat $fn0 | awk -F '\t' '{ print $1"\t"$2"\t"$3 }' | sort -k1,1 -k2,2n > $fn1
gzip $fn1
rm $fn0
done

zcat *.gz | sort -k1,1 -k2,2n > tmp.bed
bedtools merge -i tmp.bed > background.bed
gzip background.bed
rm tmp.bed
