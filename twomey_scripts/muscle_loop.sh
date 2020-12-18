for file in *.fa
do
    muscle $file -o muscle/$file
done
