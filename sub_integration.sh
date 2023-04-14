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

#seur1=$dir/Ad1_agg_soupx/v4/basic_analysis/Ad1_agg_soupx_v4_basic_subset.rds
#seur2=$dir/mcc_paper/v2/Ad3_soupx_v5.rds
#seur3=$dir/Ad9_soupx/v5/basic_analysis/Ad9_soupx_v5_basic.rds
#seur4=$dir/Ad36_soupx/v5/basic_analysis/Ad36_soupx_v5_basic.rds

seur5=/wynton/group/reiter/lauren/E2f7_em3/WT_A/v5/basic_analysis/WT_A_v5_basic.rds
seur6=/wynton/group/reiter/lauren/E2f7_em3/WT_B/v5/basic_analysis/WT_B_v5_basic.rds
seur7=/wynton/group/reiter/lauren/E2f7_em3/Hom_A/v5/basic_analysis/Hom_A_v5_basic.rds
seur8=/wynton/group/reiter/lauren/E2f7_em3/Hom_B/v5/basic_analysis/Hom_B_v5_basic.rds

id=E2f7_em3_merge
version=v2

Rscript /wynton/group/reiter/lauren/seurat_integration/seurat_integration.R $seur5 $seur6 $seur7 $seur8 $id $version

#Rscript seurat_integration.R ~/10X_082519/Ad3_v3_scran.rds ~/10X_082519/Ad9_v3_scran.rds ~/10X_082519/Ad36_v5_scran.rds ~/10X_070219/WtAd3_v16_scran.rds ~/10X_042519/Tdsort_v6_scran.rds ~/10X_042519/NS_v17_scran.rds ~/10X_042519/Ad1_agg_v7_scran.rds merge_04_07_08_2019 v1

#Rscript seurat_integration.R ~/10X_082519/Ad3_v3_scran.rds ~/10X_082519/Ad9_v3_scran.rds ~/10X_082519/Ad36_v5_scran.rds ~/10X_070219/WtAd3_v16_scran.rds ~/10X_042519/NS_v17_scran.rds ~/10X_042519/Ad1_agg_v7_scran.rds merge_04_07_08_2019_NoTdSort v1

#Rscript seurat_integration.R ~/10X_082519/Ad3_v3_scran.rds ~/10X_082519/Ad9_v3_scran.rds ~/10X_082519/Ad36_v5_scran.rds ~/10X_070219/WtAd3_v16_scran.rds ~/10X_042519/Ad1_agg_v7_scran.rds merge_04_07_08_2019_WTonly v1

#Rscript seurat_integration.R ~/10X_082519/Ad3_v3_scran.rds ~/10X_082519/Ad9_v3_scran.rds ~/10X_082519/Ad36_v5_scran.rds ~/10X_070219/WtAd3_v16_scran.rds ~/10X_042519/Ad1_agg_v7_scran.rds ~/10X_082519/CT_d4_v2_scran.rds ~/10X_082519/CT_Ad3_v5_scran.rds merge_04_07_08_2019_WT_CT v1

#Rscript seurat_integration.R ~/10X_082519/Ad3_v3_scran.rds ~/10X_070219/WtAd3_v16_scran.rds ~/10X_042519/Ad1_agg_v7_scran.rds ~/10X_082519/CT_d4_v2_scran.rds ~/10X_082519/CT_Ad3_v5_scran.rds merge_04_07_08_2019_WT_CT_d4_Ad1_Ad3 v1

#Rscript seurat_integration.R ~/10X_082519/Ad3_v3_scran.rds ~/10X_070219/WtAd3_v16_scran.rds ~/10X_042519/Ad1_agg_v7_scran.rds ~/10X_082519/CT_d4_v3_scran.rds ~/10X_082519/CT_Ad3_v5_scran.rds merge_04_07_08_2019_WT_CT_d4_Ad1_Ad3 v2

#Rscript seurat_integration.R ~/well1_filtered.rds ~/well2_filtered.rds ~/well3_filtered.rds ~/well4_filtered.rds all_wells v1

#Rscript seurat_integration.R ~/well1_filtered.rds ~/well2_filtered.rds ~/well3_filtered_v2.rds ~/well4_filtered_v2.rds all_wells v2_lr

