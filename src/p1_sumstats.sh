cd ~/Cancer_GenomicSEM
echo "== build mini reference from mtag_in3 (SNP A1 A2 MAF) =="
head -1 mtag_in3/LUAD.txt
awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","MAF"; next} {print $1,$2,$3,$7}' \
  mtag_in3/LUAD.txt > mini_ref_mtag.txt
wc -l < mini_ref_mtag.txt
head -3 mini_ref_mtag.txt
