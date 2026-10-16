cd ~/Cancer_GenomicSEM
echo "== mtag_in* dirs =="
ls -l mtag_in mtag_in2 mtag_in3 2>/dev/null | head -40
echo "== mtag_out head =="
ls -l mtag_out | head -15
echo "== header of mtag_out/mtag_all_trait_1.txt =="
head -3 mtag_out/mtag_all_trait_1.txt
echo "== rows in mtag_out/mtag_all_trait_1.txt =="
wc -l < mtag_out/mtag_all_trait_1.txt
echo "== header of first mtag_in file =="
for d in mtag_in mtag_in2 mtag_in3; do
  f=$(ls $d 2>/dev/null | head -1)
  if [ -n "$f" ]; then echo "--- $d/$f"; head -3 "$d/$f"; wc -l < "$d/$f"; fi
done
echo "== mtag_all.log tail =="
tail -20 mtag_out/mtag_all.log 2>/dev/null
