#$ -S /bin/sh

#Script to run human  samples from Chimp Samples
#to get Whippet quant output

path="/share/ScratchGeneral/gabrod/circRNASamples_Primates/Chimp/Cutadapt/"

out_path="/share/ScratchGeneral/gabrod/circRNASamples_Primates/WhippetOut_circRNAs/"

Chimp_indx="/share/ClusterShare/biodata/contrib/rweather/species_index/species_index/Pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.95.jls"

cd $path

Ones=(*_1_cut.fastq.gz)
Twos=(*_2_cut.fastq.gz)

length=${#Ones[@]}

for ((i=0; i<= $length; i++)); do

        mkdir -p $out_path/${Ones[$i]%_*}

        Whippet_ChimpSamples="/share/ClusterShare/biodata/contrib/rweather/julia/julia-9d11f62bcb/bin/julia /home/robwea/bin/whippet/Whippet.jl/bin/whippet-quant.jl --circ $path/"${Ones[$i]}" $path/"${Twos[$i]}" -o $out_path/${Ones[$i]%_*}/${Ones[$i]%_*} -x $Chimp_indx"

        echo "$Whippet_ChimpSamples" | qsub -l mem_requested=40G,tmp_requested=40G -N Whippet_${Ones[$i]%_*} -q long.q -V

done

