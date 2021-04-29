# -S /bin/sh

#Intersect Exons from CE-GTF of primates with corresponding circRNAs coordinates

#circRNAs files path
path_circ="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ForLiftOver/PSI_05"

#Exons CE Whippet path
path_exon="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10"


###For circRNAs

cd $path_circ

#files_circ=(*.sort)
files_circ=(LemurCirc_Coord_p05.sort)

###For Exons

cd $path_exon

#files_exon=(*.filt.bed)
files_exon=(Lemur_Exons_GTF-CE10.filt.bed)

length=${#files_exon[@]}

for ((i=0; i <= $length; i++)); do

        int="bedtools intersect -loj -a $path_exon/"${files_exon[$i]}" -b $path_circ/"${files_circ[$i]}"  > $path_exon/${files_exon[$i]%%_*}_ExonsCircRNA_GTF-CE10.bed"

        echo "$int" | qsub -l mem_requested=20G,tmp_requested=20G -N Int_${files_circ[$i]%%.*}_circ_Exon -q short.q -V

done

