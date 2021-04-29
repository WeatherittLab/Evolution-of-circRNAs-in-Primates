# -S /bin/sh

#Merge Exons from Whippet

#Exons CE Whippet path
path_whipp="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10"


###For CE Whippet

cd $path_whipp

files=(Coord*.sort)
files=(Coord_Lemur10.sort)

length=${#files[@]}

for ((i=0; i <= $length; i++)); do

        merge="bedtools merge -i $path_whipp/"${files[$i]}" > $path_whipp/${files[$i]%%.*}.merge"

        echo "$merge" | qsub -l mem_requested=20G,tmp_requested=20G -N Merge_${files[$i]%%.*}_ExonsCECoord -q short.q -V

done

