rm *.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_peaks_Occ1_Fcx1_Cbh1_joint_spp_RMMMcombined_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.Ast.diffPeaks.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.End.diffPeaks.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.Ex.diffPeaks.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.Mic.diffPeaks.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.In.diffPeaks.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.Opc.diffPeaks.bed.gz
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97887/suppl/GSE97887_occ1_MAINSPLITS.Oli.diffPeaks.bed.gz

gunzip *.bed.gz
for f in Ast End Ex Mic In Opc Oli
do
fn0=GSE97887_occ1_MAINSPLITS.$f.diffPeaks.bed
fn1=$f.diffPeaks.bed 
cat $fn0 | awk -F '\t' '{ print $1"\t"$2"\t"$3 }' | tail -n +2 | sort -k1,1 -k2,2n > $fn1
rm $fn0
done

fn0=GSE97887_peaks_Occ1_Fcx1_Cbh1_joint_spp_RMMMcombined_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.bed
fn1=background.bed
cat $fn0 | awk -F '\t' '{ print $1"\t"$2"\t"$3 }' | sort -k1,1 -k2,2n > $fn1
rm $fn0
gzip *.bed