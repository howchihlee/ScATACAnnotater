wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785_scATAC-Hematopoiesis-All.mtx.gz
mv GSE129785_scATAC-Hematopoiesis-All.mtx.gz GSE129785_scATAC-Hematopoiesis-All.mtx
python get_common_peaks.py

# generate background
zcat GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz | tail -n+2 | awk -F"_" '$1=$1' OFS="\t" > background.bed
gzip *.bed

# clean up
rm GSE129785_scATAC-Hematopoiesis-All*
