set -e
cd ~/Cancer_GenomicSEM
echo "== disk =="
df -h . | tail -2
echo "== loci273 lines =="
wc -l < loci273.txt

echo "== expanding reference (8.55M rows) =="
if [ ! -s ref_a1a2.txt ]; then
  zcat reference.1000G.SNP.A1A2.MAF.txt.gz > ref_a1a2.txt
fi
wc -l < ref_a1a2.txt

echo "== building mini reference: 273 loci + every 200th SNP as background =="
awk -v OFS="\t" '
  NR==FNR { k[$1]=1; next }
  FNR==1 { print "SNP","A1","A2","MAF"; next }
  {
    if ($1 in k) { print $1,$2,$3,$4; next }
    c++
    if (c % 200 == 0) print $1,$2,$3,$4
  }
' loci273.txt ref_a1a2.txt > mini_ref.txt

echo "== mini reference =="
wc -l < mini_ref.txt
head -3 mini_ref.txt

echo "== how many of the 273 loci made it in =="
awk 'NR==FNR{k[$1]=1;next} (FNR>1 && ($1 in k)){n++} END{print n+0}' loci273.txt mini_ref.txt
