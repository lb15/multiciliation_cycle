#any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o ~/log/                        #-- output directory (fill in)
#$ -e ~/log/                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=50G
#$ -l scratch=10G
#$ -l h_rt=72:00:00
#$ -m ea                           #--email when done
#$ -M Lauren.Byrnes@ucsf.edu        #--email

module load CBI
module load r/4.0.5

dir=/wynton/group/reiter/lauren/

seur5=/wynton/group/reiter/lauren/E2f7_em3/WT_A/v5/basic_analysis/WT_A_v5_basic.rds
seur6=/wynton/group/reiter/lauren/E2f7_em3/WT_B/v5/basic_analysis/WT_B_v5_basic.rds
seur7=/wynton/group/reiter/lauren/E2f7_em3/Hom_A/v5/basic_analysis/Hom_A_v5_basic.rds
seur8=/wynton/group/reiter/lauren/E2f7_em3/Hom_B/v5/basic_analysis/Hom_B_v5_basic.rds

id=E2f7_em3_merge
version=v2

Rscript /wynton/group/reiter/lauren/seurat_integration/seurat_integration.R $seur5 $seur6 $seur7 $seur8 $id $version

