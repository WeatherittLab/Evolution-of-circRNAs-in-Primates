# -S /bin/sh

#Sort Exons from Whippet CE and from GTF of primates

#Exons GTF files path
path_gtf="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords"

#Exons CE Whippet path
path_whipp="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10"


###For GTF

#cd $path_gtf

#files=(*Exons.txt)

#length=${#files[@]}

#for ((i=0; i <= $length; i++)); do

	#sort="bedtools sort -i $path_gtf/"${files[$i]}" > $path_gtf/${files[$i]%%.*}.sort"

	#echo "$sort" | qsub -l mem_requested=20G,tmp_requested=20G -N Sort_${files[$i]%%.*}_ExonsGTFCoord -q short.q -V

#done



###For CE Whippet

cd $path_whipp

#files=(Coord*)
files=(Coord_Lemur10)

length=${#files[@]}

for ((i=0; i <= $length; i++)); do

        sort="bedtools sort -i $path_whipp/"${files[$i]}" > $path_whipp/${files[$i]%%.*}.sort"

        echo "$sort" | qsub -l mem_requested=20G,tmp_requested=20G -N Sort_${files[$i]%%.*}_ExonsCECoord -q short.q -V

done

