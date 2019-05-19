echo "/Bash script run successfully/"
mask=$1
home=`pwd`

cd ./Files
in=*.fasta
for x in $in
do
 echo "/Clustalo run/"
 clustalw -INFILE="$x" -OUTFILE="$x.aligned.fasta" -TYPE=DNA -OUTPUT=FASTA
 echo "/Done/"
done

a=*.fasta.aligned.fasta
for x in $a
do
 echo "/Removing directory/"
 rm -R ${x%.fasta.aligned.fasta}"_"$mask
 echo "/Python script run/" 
 python "$home/alignment_reader_4.0.py" $x ${x%.fasta.aligned.fasta} $mask
 echo "/Done/"
 cd ${x%.fasta.aligned.fasta}"_"$mask
 pwd
 echo "/R script run/"
 Rscript "$home/R_plot.R" --no-save --no-restore --args `pwd` $mask ${x%.fasta.aligned.fasta}
 echo "/Done/"
 echo "/Cleaning the derictory/"
 mv ../${x%.fasta.aligned.fasta}.dnd ./
 mv ../$x ./
 echo "/Done/"
 echo "/Drawing trees/"
 fdrawtree -intreefile ${x%.fasta.aligned.fasta}.dnd -plotfile ${x%.fasta.aligned.fasta}"_"$mask"_"$wind.fdrawtree -auto
 fdrawgram -intreefile ${x%.fasta.aligned.fasta}.dnd -plotfile ${x%.fasta.aligned.fasta}"_"$mask"_"$wind.c.fdrawgram -auto -style c -previewer n
 cd ..
 echo "/Done/"
done
echo "/Bash script end successfully/"



