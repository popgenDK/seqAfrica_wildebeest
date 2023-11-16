#!/bin/bash

num_procs=6
num_jobs="\j"

inputSFS="Wildebeest_Black-Brindled4-WwhiteBeard.folded.sfs"
inputTPL="wildebeest.model0.folded.tpl"
inputEST="wildebeest.model0.folded.est"


for i in {1..100}; do

   while (( ${num_jobs@P} >= num_procs )); do
      #echo "${num_jobs@P} "
    wait -n
  done

ln -s  $inputEST model0_folded_${i}.est
ln -s  $inputTPL model0_folded_${i}.tpl
ln -s  $inputSFS model0_folded_${i}_MSFS.obs

nice -n10 /home/users/xi/software/fsc27_linux64/fsc2702 -t model0_folded_${i}.tpl -n500000 \
-e model0_folded_${i}.est \
-M -L100 -u -m -c10 -C 100 -r ${i} 2> model0_folded_${i}.log &

done
