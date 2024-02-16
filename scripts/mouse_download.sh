#!/bin/bash

echo "[mouse-download]: Starting ..."

DATA_DIR="../data/external"

echo "[mouse-download]: -- Deng et al., 2014"
mkdir $DATA_DIR/Deng_et_al_2014
# ffq -l 1 GSE45719 > $DATA_DIR/Deng_et_al_2014/metadata.json
wget -q "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45719/suppl/GSE45719_RAW.tar" \
  -O $DATA_DIR/Deng_et_al_2014/GSE45719_RAW.tar && \
  tar -xf $DATA_DIR/Deng_et_al_2014/GSE45719_RAW.tar -C $DATA_DIR/Deng_et_al_2014

echo "[mouse-download]: -- Nowotschin et al., 2019"
mkdir $DATA_DIR/Nowotschin_et_al_2019
wget -q "https://s3.amazonaws.com/dp-lab-data-public/mouse_endoderm/sc_endoderm_all_cells.h5ad" \
  -O $DATA_DIR/Nowotschin_et_al_2019/sc_endoderm_all_cells.h5ad

echo "[mouse-download]: -- Biase et al., 2014"
# ffq -l 1 GSE57249 > Biase_et_al_2014/metadata.json
mkdir $DATA_DIR/Biase_et_al_2014
wget -q -O - "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57249/suppl/GSE57249_fpkm.txt.gz" | \
  gunzip -c > $DATA_DIR/Biase_et_al_2014/GSE57249_fpkm.txt
wget -q "https://scrnaseq-public-datasets.s3.amazonaws.com/manual-data/biase/biase_cell_types.txt" \
  -O $DATA_DIR/Biase_et_al_2014/biase_cell_types.txt

echo "[mouse-download]: -- Posfai et al., 2017"
# ffq -l 1 GSE84892 > Posfai_et_al_2017/metadata.json
mkdir $DATA_DIR/Posfai_et_al_2017
wget -q -O - "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84892/suppl/GSE84892_rpkm.txt.gz" | \
  gunzip -c > $DATA_DIR/Posfai_et_al_2017/GSE84892_rpkm.txt

echo "[mouse-download]: -- Goolam et al., 2016"
mkdir $DATA_DIR/Goolam_et_al_2016
wget -q "https://www.ebi.ac.uk/biostudies/files/E-MTAB-3321/Goolam_et_al_2015_count_table.tsv" \
  -O $DATA_DIR/Goolam_et_al_2016/Goolam_et_al_2015_count_table.tsv
wget -q "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3321/E-MTAB-3321.sdrf.txt" \
  -O $DATA_DIR/Goolam_et_al_2016/E-MTAB-3321.sdrf.txt

echo "[mouse-download]: -- Boroviak et al., 2015"
mkdir $DATA_DIR/Boroviak_et_al_2015
wget -q "https://ars.els-cdn.com/content/image/1-s2.0-S1534580715006589-mmc2.xlsx" \
  -O $DATA_DIR/Boroviak_et_al_2015/1-s2.0-S1534580715006589-mmc2.xlsx

echo "[mouse-download]: -- Chen et al., 2016"
mkdir $DATA_DIR/Chen_et_al_2016
wget -q -O - "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74155&format=file&file=GSE74155%5FBlastocyst%5Fgenes%5Ffpkm%2Etxt%2Egz" | \
  gunzip -c > $DATA_DIR/Chen_et_al_2016/GSE74155_Blastocyst_genes_fpkm.txt
wget -q -O - "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74155/matrix/GSE74155_series_matrix.txt.gz" | \
  gunzip -c > $DATA_DIR/Chen_et_al_2016/GSE74155_series_matrix.txt

echo "[mouse-download]: Job done ..."