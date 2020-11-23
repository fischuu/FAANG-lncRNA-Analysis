# Preparations:
# wget https://raw.githubusercontent.com/CSCfi/allas-cli-utils/master/allas_conf
# sudo apt install s3cmd python-openstackclient
# source allas_conf --mode s3cmd --user fischerd --project project_2001289

a-list lncrna_bos_taurus/bos_taurus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > bos_taurus_gtf.txt
a-list lncrna_capra_hircus/capra_hircus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > capra_hircus_gtf.txt
a-list lncrna_equus_caballus/equus_caballus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > equus_caballus_gtf.txt
a-list lncrna_gallus_gallus/gallus_gallus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > gallus_gallus_gtf.txt
a-list lncrna_ovis_aries/ovis_aries/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > ovis_aries_gtf.txt
a-list lncrna_sus_scrofa/sus_scrofa/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > sus_scrofa_gtf.txt

a-list lncrna_sus_scrofa/sus_scrofa/BAM | grep '.bam' | cut -d'/' -f4  | cut -d'_' -f1 | cut -d'.' -f1 | sort | uniq > sus_scrofa_bam.csv
a-list lncrna_bos_taurus/bos_taurus/BAM | grep '.bam' | cut -d'/' -f4  | cut -d'_' -f1 | cut -d'.' -f1 | sort | uniq > bos_taurus_bam.csv
a-list lncrna_ovis_aries/ovis_aries/BAM | grep '.bam' | cut -d'/' -f4  | cut -d'_' -f1 | cut -d'.' -f1 | sort | uniq > ovis_aries_bam.csv
a-list lncrna_gallus_gallus/gallus_gallus/BAM | grep '.bam' | cut -d'/' -f4  | cut -d'_' -f1 | cut -d'.' -f1 | sort | uniq > gallus_gallus_bam.csv
a-list lncrna_equus_caballus/equus_caballus/BAM | grep '.bam' | cut -d'/' -f4  | cut -d'_' -f1 | cut -d'.' -f1 | sort | uniq > equus_caballus_bam.csv
a-list lncrna_capra_hircus/capra_hircus/BAM | grep '.bam' | cut -d'/' -f4  | cut -d'_' -f1 | cut -d'.' -f1 | sort | uniq > capra_hircus_bam.csv

# a-list lncrna_sus_scrofa/sus_scrofa/FASTQ/ | tr -s ' ' | cut -d ' ' -f 4 | cut -d '/' -f 6 | grep gz.md5 | cut -d '_' -f1 | sort | uniq > sus_scrofa_fastq.csv
# a-list lncrna_bos_taurus/bos_taurus/FASTQ/ | tr -s ' ' | cut -d ' ' -f 4 | cut -d '/' -f 6 | grep gz.md5 | cut -d '_' -f1 | sort | uniq > bos_taurus_fastq.csv
# a-list lncrna_ovis_aries/ovis_aries/FASTQ/ | tr -s ' ' | cut -d ' ' -f 4 | cut -d '/' -f 6 | grep gz.md5 | cut -d '_' -f1 | sort | uniq > ovis_aries_fastq.csv
# a-list lncrna_gallus_gallus/gallus_gallus/FASTQ/ | tr -s ' ' | cut -d ' ' -f 4 | cut -d '/' -f 6 | grep gz.md5 | cut -d '_' -f1 | sort | uniq > gallus_gallus_fastq.csv
# a-list lncrna_equus_caballus/equus_caballus/FASTQ/ | tr -s ' ' | cut -d ' ' -f 4 | cut -d '/' -f 6 | grep gz.md5 | cut -d '_' -f1 | sort | uniq > equus_caballus_fastq.csv
# a-list lncrna_capra_hircus/capra_hircus/FASTQ/ | tr -s ' ' | cut -d ' ' -f 4 | cut -d '/' -f 6 | grep gz.md5 | cut -d '_' -f1 | sort | uniq > capra_hircus_fastq.csv

# Import the sample list from the google folder

# Import required libaries
  library("googlesheets4")

# Project defaults
  projFolder <- "/home/fischuu/ownCloud/Luke/Projects/FAANG-lncRNA/Analysis"
  options(stringsAsFactors=FALSE)

# Import the latest currated sample list from google drive
  sampleList <- read_sheet("https://docs.google.com/spreadsheets/d/1iQNOeGk2d00-Zjnk5e56pjCABBgqqpUmOSOTUOdYwAM/edit#gid=507446984", col_names=FALSE)

# Import pipeline results
  susSamples <- as.vector(as.matrix(read.table(file.path(projFolder, "sus_scrofa_gtf.txt"), header=FALSE)))
  bosSamples <- as.vector(as.matrix(read.table(file.path(projFolder, "bos_taurus_gtf.txt"), header=FALSE)))
  ovisSamples <- read.table(file.path(projFolder, "ovis_aries_gtf.txt"), header=FALSE)
  ovisSamples <- as.vector(as.matrix(as.data.frame(t(ovisSamples)[!grepl("temp", t(ovisSamples))], ncol=1)))  # Something went wrong with the pipe and it created some artifacts...,
  gallusSamples <- as.vector(as.matrix(read.table(file.path(projFolder, "gallus_gallus_gtf.txt"), header=FALSE)))
  equusSamples <- as.vector(as.matrix(read.table(file.path(projFolder, "equus_caballus_gtf.txt"), header=FALSE)))
  capraSamples <- as.vector(as.matrix(read.table(file.path(projFolder, "capra_hircus_gtf.txt"), header=FALSE)))

  totalSamples <- sum(length(susSamples),
                      length(bosSamples),
                      length(ovisSamples),
                      length(gallusSamples),
                      length(equusSamples),
                      length(capraSamples))

  sum(is.element(as.vector(as.matrix(sampleList[,1])),susSamples))+
  sum(is.element(as.vector(as.matrix(sampleList[,1])),bosSamples))+
  sum(is.element(as.vector(as.matrix(sampleList[,1])),ovisSamples))+
  sum(is.element(as.vector(as.matrix(sampleList[,1])),gallusSamples))+
  sum(is.element(as.vector(as.matrix(sampleList[,1])),equusSamples))+
  sum(is.element(as.vector(as.matrix(sampleList[,1])),capraSamples))

  