### Jinliang Yang
### Sept. 21th, 2016

#### genotypes
gt0 <- read.csv("data/AllZeaGBSv2.7_publicSamples_metadata20140411.csv", header=T)
gt <- subset(gt0, Project %in% c("2010 Ames Lines", "AMES Inbreds", "Ames282") )
gt$Pedigree <- as.character(gt$Pedigree)

### 
sam <- read.csv("data/SAM_cellcount.csv")
sam <- subset(sam, !is.na(Count_Cells))
id <- as.character(unique(sam$Genotype))

sub <- subset(gt, toupper(DNASample) %in% toupper(id) | Pedigree %in% toupper(id))
as.character(unique(sub$DNASample))
length(unique(sub$DNASample, sub$Pedigree))

sub <- sub[order(sub$DNASample),]
write.table(sub, "cache/cellnum_GBS_sampleid.csv", sep=",", row.names=FALSE, quote=FALSE)
###>>> Manually curated the ids

#### extract subset of the data: biallelic SNPs only and no-update
ids <- read.csv("cache/cellnum_GBS_sampleid_curated.csv")
write.table(ids[, 1], "cache/id14.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

tabix <- "tabix -p vcf AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3_sorted.vcf.gz"
cmd <- paste("bcftools view -S /home/jolyang/Documents/Github/zmHapMap/cache/id14.txt",
             "--no-update -m2 -M2 -v snps", 
             "AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3_sorted.vcf.gz",
             "-Oz -o AllZeaGBSv2.7_id14_imputedV3b.vcf.gz")


### extract id and convert to PLINK
cmd1 <- c("cd /home/jolyang/dbcenter/HapMap/HapMap3")
cmd2 <- paste("plink -vcf AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3_sorted.vcf.gz",
              "--biallelic-only --snps-only --set-missing-var-ids @_# --out plink_chr1", 
              "--allow-extra-chr --freq")



set_farm_job(slurmsh = "slurm-script/bcf2plink.sh",
             shcode = c(cmd1, cmd2), wd = NULL, jobid = "maf",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemh", "3", "23000"))


#########################################
for(i in 1:10){
  shid <- paste0("slurm-script/runplink_", i, ".sh")
  cmd1 <- c("cd /home/jolyang/dbcenter/HapMap/HapMap3")
  cmd2 <- paste0("plink -vcf merged_flt_c", i, ".vcf.gz --biallelic-only --snps-only", 
                 " --set-missing-var-ids @_# --out plink_chr", i, 
                 " --allow-extra-chr --freq")
  cat(c(cmd1, cmd2), file=shid, sep="\n", append=FALSE)
}
shcode <- "sh slurm-script/runplink_$SLURM_ARRAY_TASK_ID.sh"

set_array_job(shid="slurm-script/runplink.sh", shcode=shcode,
              arrayjobs="1-10", wd=NULL, jobid="plink", email="yangjl0930@gmail.com",
              run = c(TRUE, "bigmemh", 4))

