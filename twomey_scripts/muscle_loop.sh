for file in *.fa
do
    muscle -in $file -out muscle/$file
done
