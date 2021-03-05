for file in *.fasta
do
    muscle -in $file -out muscle/$file
done
