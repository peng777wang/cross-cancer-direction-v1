cd ~/Cancer_GenomicSEM
echo "== header of clean/LUAD.txt.gz =="
zcat clean/LUAD.txt.gz | head -3
echo "== header of clean/GC.txt.gz =="
zcat clean/GC.txt.gz | head -3
echo "== header of munged LUAD.sumstats.gz =="
zcat LUAD.sumstats.gz | head -3
echo "== loci present in clean/LUAD.txt.gz =="
zcat clean/LUAD.txt.gz | awk 'NR==FNR{k[$1]=1;next} ($1 in k){n++} END{print n+0}' loci273.txt -
echo "== loci present in clean/GC.txt.gz =="
zcat clean/GC.txt.gz | awk 'NR==FNR{k[$1]=1;next} ($1 in k){n++} END{print n+0}' loci273.txt -
echo "== loci present in LUAD.sumstats.gz =="
zcat LUAD.sumstats.gz | awk 'NR==FNR{k[$1]=1;next} ($1 in k){n++} END{print n+0}' loci273.txt -
echo "== row counts =="
zcat clean/LUAD.txt.gz | wc -l
zcat LUAD.sumstats.gz | wc -l
