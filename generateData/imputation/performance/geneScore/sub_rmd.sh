#INITIALISE FOLDERS
my_folder=/nfs/research1/marioni/alsu
out_folder=${my_folder}/clust_out/spatial
err_folder=${my_folder}/clust_err/spatial

#SELECT SCRIPT
#If you change this, you MUST update the wrapper's grep
script_name=score

#CHOOSE PARAMETERS
#RAM in megabytes
memory=90000
r_command="rusage[mem=${memory}]"
#num_processors
nproc=3

smg=/nfs/research1/marioni/alsu/singularity/R1.simg
script=/nfs/research1/marioni/alsu/spatial/mouse_embryo/amScripts_mouseEmbryoFISH/generateData/imputation/performance/geneScore/run_rmd.R

bsub -q research-rh74 -e ${err_folder}/${script_name} \
-o ${out_folder}/${script_name} \
-M $memory -R $r_command -n $nproc -J ${script_name} \
"singularity exec $smg Rscript $script"
