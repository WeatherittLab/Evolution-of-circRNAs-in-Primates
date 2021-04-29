# -S /bin/sh

#Sort circRNAs coordinates of primates

#circRNAs files path
path="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ForLiftOver/PSI_05"

cd $path

#files=(*Coord_p05.txt)

files=(LemurCirc_Coord_p05.txt)

length=${#files[@]}

for ((i=0; i <= $length; i++)); do

	sort="bedtools sort -i $path/"${files[$i]}" > $path/${files[$i]%%.*}.sort"

	echo "$sort" | qsub -l mem_requested=20G,tmp_requested=20G -N Sort_${files[$i]%%.*}_CircCoord -q short.q -V

done
