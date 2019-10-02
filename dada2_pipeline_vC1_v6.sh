#!/bin/sh

#species classification is possible but currently disabled
#taxonomy assignment occurs twice in this worklow; using DADA2 methods and using RDP classifier

conda_base=$(conda info --base)
export scripts='/path/to/dada2_pipeline_scripts'
#need better solution for this one:
export qiime='/path/to/qiime1/bin'
#custom species classifier can go here
#export v6species='/path/to/V6_species'
export database='/path/to/dada2_pipeline_databases'
export HDF5_USE_FILE_LOCKING='FALSE'

usage() {
	echo "
--------------------------------------------------------------------------------------------------------------------------------------

	USAGE: available options
		mandatory options:
				-i or --input : provide the name of your input folder (containing gzipped raw FASTQ files)
		optional options:
				-r or --run : provide the name of your sequencing run e.g. run17
				-t or --threads : give the number of threads to use (default is 2)
				-dt or --taxonomydb : give the name of your taxonomy database
				-ds or --speciesdb : give the name of your species taxonomy database
	  Help options:
				-h or --help : pop up the usage guide

		e.g. :
				sbatch dada2_pipeline_vC1_v6.sh -i raw_reads -r run17 -t 2 -dt rdp_train_set_16.fa.gz -ds rdp_species_assignment_16.fa.gz
				raw_reads is the name of the input folder
				run17 is the run name
				2 is the total number of threads
				rdp_train_set_16.fa.gz is the name of your taxonomy database
				rdp_species_assignment_16.fa.gz is the name of your species taxonomy database

--------------------------------------------------------------------------------------------------------------------------------------

	"
	}


#Give at least two variables
if [ $# -lt 2 ]
then
	echo "--------------give at least two arguments ( -i input_folder ) to represent input_name--------------"
	usage
	exit 1
fi


#default settings
threads=2
run_name=run_unknown
db_taxonomy=rdp_train_set_16.fa.gz
db_species=rdp_species_assignment_16.fa.gz


#Flag options
nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case $1 in
    -r | --run)
      run_name=$2
			echo "-------------------------------This is run: $2--------------------------------------------"
			shift 2
			;;

		-i | --input)
      input_folder=$2
			echo "-------------------------------input_folder: $2----------------------------------------"
			shift 2
			;;

    -t | --threads)
      threads=$2
  		echo "----------------------------The number of our threads: $2---------------------------------"
  		shift 2
  		;;

		-dt | --taxonomydb)
			db_taxonomy=$2
			echo "------------------------------databases you are using: $2---------------------------------"
			shift 2
			;;

		-ds | --speciesdb)
			db_species=$2
			echo "-------------------------------databases you are using: $2----------------------------------"
			shift 2
			;;

		-h | --help)
			usage
			exit 1
			;;
	esac
done

#
if [ -z "$input_folder" ]
then
	echo "---------------------------Can not go on without your input -------------------------------------"
	echo "--------------give at least two arguments ( -i or --input ) to represent input_folder---------------"
	usage
	exit 1
fi

source $conda_base/bin/activate dada2_flow_env

my_output=my_output_$run_name
mkdir $my_output
mkdir $my_output/End_user_results
mkdir $my_output/RDS_files
mkdir trimmed_reads
mkdir primer_trimmed_reads


#counting raw reads from unfiltered .gz files
echo "------------------------------STEP 1: Original raw data counting---------------------------------"
#use zcat for Linux
for i in $input_folder/*_R1_001.fastq.gz; do echo -n "$i" $'\t'; echo $(gzcat $i | wc -l) / 4 | bc; done > $my_output/LP_${run_name}_readCount.txt
#text editing: delete prefix and suffix
awk '{split($1,a,"/"); $1=a[length(a)]"\t";}1' $my_output/LP_${run_name}_readCount.txt > $my_output/LP_${run_name}_readCount.intermediate1.txt
awk '{split($1,a,"_"); $1=a[1]"\t";}1' $my_output/LP_${run_name}_readCount.intermediate1.txt > $my_output/LP_${run_name}_readCount.intermediate2.txt
#text editing: add headers to each column
echo -e "Sample\tRaw_reads" | cat - $my_output/LP_${run_name}_readCount.intermediate2.txt > $my_output/LP_${run_name}_readCount.final.txt
rm $my_output/LP_${run_name}_readCount.intermediate1.txt
rm $my_output/LP_${run_name}_readCount.intermediate2.txt


#We run Trimmomatic to remove low quality sequences from raw data (input_file_name)
echo "----------------------------STEP 2: Trim Low Quality Sequences------------------------------------------"
perl $scripts/run_trimmomatic_amplicon_mode.pl $input_folder/*


#move files that have "paired_trimmed" in their file name to trimmed_reads
echo "-----------------------------------STEP 3: Move Files---------------------------------------------------"
mv $input_folder/*paired_trimmed* trimmed_reads
rm $input_folder/*single_trimmed*


#primer sequences need to be removed so that DADA2 won't generate inaccurate data or chokes on ambigious bases
echo "-----------------------------STEP 4: Remove Primer Sequences--------------------------------------------"
parallel --link --jobs $threads 'cutadapt --pair-filter any --no-indels -e 0.15 --discard-untrimmed -g file:$database/amplicon_primers.fasta -G file:$database/amplicon_primers.fasta -o primer_trimmed_reads/{1/}.gz -p primer_trimmed_reads/{2/}.gz {1} {2} > primer_trimmed_reads/{1/}_cutadapt_log.txt' ::: trimmed_reads/*_R1_*.fastq ::: trimmed_reads/*_R2_*.fastq

#filter out reads with ambigious base calls (Ns)
echo "---------------------------------STEP 5: filter out reads with ambigious base calls (Ns)----------------"
$scripts/dada2_filter.R -f primer_trimmed_reads --maxN 0 --truncQ 0 --threads $threads --truncLen 0 --maxEE 2,2 --f_match _R1_.*fastq.gz --r_match _R2_.*fastq.gz -o $my_output/${run_name}_filtered_fastqs


#dada2 ASV calling - creates seqtab.rds file
echo "------------------------STEP 6: creates seqtab.rds file, hold on----------------------------------------"
$scripts/dada2_inference.R -f $my_output/${run_name}_filtered_fastqs/ --seed 4124 -t $threads --verbose --plot_errors --minOverlap 12 -o $my_output/${run_name}_seqtab.rds


#remove spurious features and assign taxonomy - creates seqtab_final.rds
echo "-----------------STEP 7: now we are filtering seqtab.rds, this step takes the longest--------------------"
$scripts/dada2_chimera_taxa.R -i $my_output/${run_name}_seqtab.rds -t 2 -r $database/$db_taxonomy -s $database/$db_species --allowMultiple --minBoot 80 --count_out $my_output/${run_name}_seqtab_final.rds --sp_out $my_output/${run_name}_species_final.rds --tax_out $my_output/${run_name}_tax_final.rds


#convert RDS objects to TSV (in classic biom format) and FASTA - creates file.biom.tsv and file_species.fasta
echo "---------------------------STEP 8: creating my_output files....almost done----------------------------------"
$scripts/convert_dada2_out.R -i $my_output/${run_name}_seqtab_final.rds -b $my_output/${run_name}_biom.tsv -f $my_output/${run_name}_species.fasta --taxa_in $my_output/${run_name}_species_final.rds --taxa_out $my_output/${run_name}_taxa_metadata.txt
#get out of the dada2 environment.... while you still can

echo "----------------------------------------STEP 9: convert TSV to BIOM file---------------------------------"
biom convert -i $my_output/${run_name}_biom.tsv -o $my_output/${run_name}_json.biom --to-json

#create the summary of your biom file
echo "--------------------------------STEP 10: create the summary of your biom file-----------------------------"
biom summarize-table -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_json_summary.txt


#create log file, text editing: extract columns we need, calculate the summation
echo "--------------------------------------------STEP 11: create log file---------------------------------------"
# paste $my_output/LP_${run_name}_readCount.final.txt $my_output/dada2_filter_read_counts.txt $my_output/dada2_inferred_read_counts.txt $my_output/dada2_nonchimera_counts.txt | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$12"\t"$15}' > $my_output/${run_name}_dada2_final1.txt
paste $my_output/LP_${run_name}_readCount.final.txt dada2_filter_read_counts.txt dada2_inferred_read_counts.txt dada2_nonchimera_counts.txt | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$12"\t"$15}' > $my_output/${run_name}_dada2_final1.txt
sed -i '1d' $my_output/${run_name}_dada2_final1.txt
echo -e "Sample\tRaw_reads\tQuality_trimmed_reads\tReads_without_'N'_bases\tDenoised_dereplicated_reads\tNonchimera_reads" | cat - $my_output/${run_name}_dada2_final1.txt > $my_output/${run_name}_dada2_final2.txt
awk 'NR==1{$7="Proportion_of_saved_reads"}NR>1{$7 = (($6)/($2))*100"\%"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $my_output/${run_name}_dada2_final2.txt > $my_output/${run_name}_dada2_final3.txt
awk '{print $0}''NR>1{sum2+=$2; sum3+=$3; sum4+=$4; sum5+=$5; sum6+=$6}END{print "Total_reads\t"sum2"\t"sum3"\t"sum4"\t"sum5"\t"sum6"\t"(sum6/sum2)*100"\%"}' $my_output/${run_name}_dada2_final3.txt > $my_output/${run_name}_dada2_QC_stats.txt

#remove intermediate files
rm $my_output/${run_name}_dada2_final1.txt
rm $my_output/${run_name}_dada2_final2.txt
rm $my_output/${run_name}_dada2_final3.txt
rm $my_output/LP_${run_name}_readCount.final.txt
rm $my_output/LP_${run_name}_readCount.txt

#move QC files to main output folder
mv dada2_filter_read_counts.txt $my_output
mv dada2_inferred_read_counts.txt $my_output
mv dada2_nonchimera_counts.txt $my_output
mv estimated_forward_err.pdf $my_output/${run_name}_estimated_forward_err.pdf
mv estimated_reverse_err.pdf $my_output/${run_name}_estimated_reverse_err.pdf
mv $my_output/*.rds $my_output/RDS_files

#taxonomy time
echo "--------------------------------------------STEP 12: taxonomy assignment ---------------------------------------"
#genus level classification
rdp_classifier classify -o $my_output/RDP_GENUS_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -f fixrank

#species level classification using a custom classifier
#rdp_classifier classify -o $my_output/RDP_V6species_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -t $v6species/rRNAClassifier.properties


perl $scripts/convert_default_settings_rdpclassifier_genus_taxa_to_biom_metadata.pl $my_output/RDP_GENUS_taxa_metadata.txt
#perl $scripts/convert_default_settings_rdpclassifier_speciesV6_taxa_to_biom_metadata.pl $my_output/RDP_V6species_taxa_metadata.txt

biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP_GENUS_c80_tax.biom --observation-metadata-fp $my_output/RDP_GENUS_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

#biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP_V6-SPECIES_c80_tax.biom --observation-metadata-fp $my_output/RDP_V6species_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

conda deactivate

$qiime/summarize_taxa.py -i $my_output/${run_name}_RDP_GENUS_c80_tax.biom -L 2,3,4,5,6 -o $my_output/${run_name}_RDP_GENUS_c80_tax_summarizeTaxa --absolute_abundance
biom convert -i $my_output/${run_name}_RDP_GENUS_c80_tax.biom -o $my_output/${run_name}_RDP_GENUS_c80_tax.tsv --to-tsv --header-key taxonomy

#$qiime/summarize_taxa.py -i $my_output/${run_name}_RDP_V6-SPECIES_c80_tax.biom -L 2,3,4,5,6,7 -o $my_output/${run_name}_RDP_V6-SPECIES_c80_tax_summarizeTaxa --absolute_abundance
#biom convert -i $my_output/${run_name}_RDP_V6-SPECIES_c80_tax.biom -o $my_output/${run_name}_RDP_V6-SPECIES_c80_tax.tsv --to-tsv --header-key taxonomy


#aureus merging
#echo "--------------------------------------------STEP 13: aureus merging---------------------------------------"
#data summation of bacteria that has "Staphylococcus_aureus"
#awk -F '\t' '$1 ~ /Staphylococcus_aureus/ {x=NF; for (i=1;i<x;i++) {sum[i] += $(i+1)} } END { for (i in sum) {print sum[i]|"sort"} }' $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7.txt > L7_intermediate1.txt
#change my_output record separator to tab instead of defalut new line
#awk 'BEGIN {ORS="\t"}; {print $0}' L7_intermediate1.txt > L7_intermediate2.txt
#adding "Staphylococcus_aureus<0.8" before the data
#awk '{print "Staphylococcus_aureus<0.8\t"$0}' L7_intermediate2.txt > L7_intermediate3.txt
#delete bacteria that has "Staphylococcus_aureus" from the original V6-SPECIES file
#awk '!/Staphylococcus_aureus/{print $0}' $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/run17.RDP212_V6-SPECIES_c80_tax_L7.txt > L7_intermediate4.txt
#add the generated Staphylococcus_aureus<0.8 line to the end of the original V6-SPECIES file
#cat L7_intermediate4.txt L7_intermediate3.txt > $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7_final.txt

#rm L7_intermediate1.txt
#rm L7_intermediate2.txt
#rm L7_intermediate3.txt
#rm L7_intermediate4.txt

#copy most important files to separate folder
cp $my_output/${run_name}_RDP_GENUS_c80_tax_summarizeTaxa/${run_name}_RDP_GENUS_c80_tax_L6.txt $my_output/End_user_results
#cp $my_output/${run_name}_RDP_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP_V6-SPECIES_c80_tax_L7.txt $my_output/End_user_results
cp $my_output/${run_name}_json_summary.txt $my_output/End_user_results
cp $my_output/${run_name}_dada2_QC_stats.txt $my_output/End_user_results

#echo "----------------------------------------FINISHED       CONGRATS--------------------------------------------"
