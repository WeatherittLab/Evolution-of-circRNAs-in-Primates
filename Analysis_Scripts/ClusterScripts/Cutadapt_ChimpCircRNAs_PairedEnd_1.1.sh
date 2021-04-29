# -S /bin/sh

#Quality trimming of paired en data of Paired End of Chimp Samples

#For missing samples
path="/share/ScratchGeneral/gabrod/circRNASamples_Primates/Chimp"

#Out path for missing samples
out_path="/share/ScratchGeneral/gabrod/circRNASamples_Primates/Chimp/Cutadapt"

Path_cutadapt="/home/gabrod/.local/bin/"

cd $path

Ones=(*_1.fastq.gz)
Twos=(*_2.fastq.gz)

length=${#Ones[@]}

for ((i=0; i <= $length; i++)); do

	cut="$Path_cutadapt/./cutadapt --pair-filter=any -a AGATCGGAAGAG -A AGATCGGAAGAG -q 20,20 --minimum-length 35:35 --trim-n -o $out_path/${Ones[$i]%%.*}_cut.fastq.gz -p $out_path/${Twos[$i]%%.*}_cut.fastq.gz $path/"${Ones[$i]}" $path/"${Twos[$i]}"  "

	echo "$cut" | qsub -l mem_requested=20G,tmp_requested=20G -N Cutadapt_${Ones[$i]%%.*}_cutadapt -q short.q -V

done

