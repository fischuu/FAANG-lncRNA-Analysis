# Get the list of processed samples per species

a-list lncrna_bos_taurus/bos_taurus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > bos_taurus_processedSamples.txt
a-list lncrna_capra_hircus/capra_hircus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > capra_hircus_processedSamples.txt
a-list lncrna_equus_caballus/equus_caballus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > equus_caballus_processedSamples.txt
a-list lncrna_gallus_gallus/gallus_gallus/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > gallus_gallus_processedSamples.txt
a-list lncrna_ovis_aries/ovis_aries/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > ovis_aries_processedSamples.txt
a-list lncrna_sus_scrofa/sus_scrofa/GTF/Ref_fc | cut -d '/' -f5 | cut -d'_' -f1 | sort | uniq > sus_scrofa_processedSamples.txt