raxmlHPC-AVX -s temp.fasta -m GTRCAT -p 12345 -n out &
sleep 2
kill $!
