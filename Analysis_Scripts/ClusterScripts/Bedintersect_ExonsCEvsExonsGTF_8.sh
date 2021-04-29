# -S /bin/sh

#Sort Exons from Whippet CE and from GTF of primates

#Exons GTF files path
path_gtf="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords"

#Exons CE Whippet path
path_whipp="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10"


###For GTF

cd $path_gtf

#files_gtf=(*Exons.sort)
files_gtf=(Lemur_Exons.sort)

###For CE Whippet

cd $path_whipp

#files_CE=(*.merge)
files_CE=(Coord_Lemur10.merge)

length=${#files_CE[@]}

for ((i=0; i <= $length; i++)); do

        int="bedtools intersect -wa -a $path_gtf/"${files_gtf[$i]}" -b $path_whipp/"${files_CE[$i]}"  > $path_whipp/${files_gtf[$i]%%_*}_Exons_GTF-CE10.bed"

        echo "$int" | qsub -l mem_requested=20G,tmp_requested=20G -N Int_${files[$i]%%.*}_GTF_CE -q short.q -V

done

