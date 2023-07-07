#!/bin/bash -login

cat << "EOF"
   ______          _ ___      
  / ____/__  _____(_) (_)___ _
 / /   / _ \/ ___/ / / / __ `/
/ /___/  __/ /__/ / / / /_/ / 
\____/\___/\___/_/_/_/\__,_/  
                                    
EOF

echo -e "
Hello there, I am Cecilia, your usearCh basEd ampliCon pIpeLine for Illumina dAta.
\n
Cecilia v.1.0 by Gian M. N. Benucci, P.hD.
email: benucci[at]msu[dot]edu
August 22, 2022
\n"

echo -e "This pipeline is based upon work supported by the Great Lakes Bioenergy Research
Center, U.S. Department of Energy, Office of Science, Office of Biological and Environmental
Research under award DE-SC0018409\n"

source ./config.yaml

cd $project_dir/rawdata/

echo -e "\n========== Comparing md5sum codes... ==========\n"

md5=$project_dir/rawdata/md5.txt

if [[ -f "$md5" ]]; then
		echo -e "\nAn md5 file exist. Now checking for matching codes.\n"
		md5sum md5* --check > tested_md5.results
		cat tested_md5.results
		resmd5=$(cat tested_md5.results | cut -f 2 -d" " | uniq)

		if [[ "$resmd5" == "OK" ]]; then
				echo -e "\n Good news! Files are identical. \n"
		else
				echo -e "\n Oh oh! You are in trouble. Your files are different from those at the source! \n"
				echo -e "\nSomething went wrong during files download. Try again, please.\n"
				exit
		fi
else
	echo "No md5 file was found. You should look for it, and start over again!"
fi

cd $project_dir/code/

echo -e "\n========== What I will do for you? Please see below... ==========\n"

echo -e "\n========== Prefiltering ==========\n" 

jid1=`sbatch 01_decompress-bash.sb | cut -d" " -f 4`
echo "$jid1: For starting, I will decompress your raw reads files."

jid2=`sbatch --dependency=afterok:$jid1 02_qualityCheck.sb | cut -d" " -f 4`
echo "$jid2: I will check the quality and generate statistics."

jid3=`sbatch --dependency=afterok:$jid1:$jid2 03_removePhix-bowtie2.sb | cut -d" " -f 4`
echo "$jid3: I will remove Phix reads with bowtie2."

# Assembling reads
if [[ "$assemble" == "yes" ]]; then
		jid4=`sbatch --dependency=afterok:$jid3 04_readAssembly-pear.sb | cut -d" " -f 4`
		echo "$jid4: I will assemble your reads with pear."
else 
		echo -e "\n You chose to NOT assemble the reads. \n";
fi

# Renaming reads 
if [[ "$assemble" == "yes" ]] && [[ "$rename" == "yes" ]]; then
		jid5=`sbatch --dependency=afterok:$jid4 05_readRename.sb | cut -d" " -f 4`
		echo "$jid5: I will rename your reads using "Sample..." as a basename."
elif [[ "$assemble" == "no" ]] && [[ "$rename" == "yes" ]]; then 
		jid5=`sbatch --dependency=afterok:$jid3 05_readRename.sb | cut -d" " -f 4`
		echo "$jid5: I will rename your reads using "Sample..." as a basename."
elif [[ "$assemble" == "no" ]] && [[ "$rename" == "no" ]] || [[ "$assemble" == "yes" ]] && [[ "$rename" == "no" ]]; then 
		echo -e "\n You chose NOT to rename the reads! I will use the original names. \n";
fi


# Stripping primers and adapters
if [[ "$assemble" == "no" ]] && [[ "$rename" == "no" ]]; then			
		jid6=`sbatch --dependency=afterok:$jid3 06_primerStrip-cutadapt.sb | cut -d" " -f 4`
		echo "$jid6: I will remove primers and adapters with cutadapt. A subset of reads is generated for double checking." 
elif  [[ "$assemble" == "yes" ]] && [[ "$rename" == "no" ]]; then		
		jid6=`sbatch --dependency=afterok:$jid4 06_primerStrip-cutadapt.sb | cut -d" " -f 4`
		echo "$jid6: I will remove primers and adapters with cutadapt. A subset of reads is generated for double checking." 
else		
		jid6=`sbatch --dependency=afterok:$jid5 06_primerStrip-cutadapt.sb | cut -d" " -f 4`
		echo "$jid6: I will remove primers and adapters with cutadapt. A subset of reads is generated for double checking." 
fi

jid7=`sbatch --dependency=afterok:$jid6 07_readTrimDemux-usearch.sb | cut -d" " -f 4`
echo "$jid7: I will calculate the max expected errors of your reads in usearch11."


if [[ "$test_length" == "yes" ]]; then
		jid8=`sbatch --dependency=afterok:$jid7 08_testLength-usearch.sb | cut -d" " -f 4`
		echo "$jid8: I will check for optimal length for trimming. This step is informative, but takes time."
else
		echo -e "\n You chose to NOT for optimal length! That's ok, I will go forward. \n";
fi


echo -e "\n========== Clustering ==========\n" 
if [[ "$test_length" == "yes" ]]; then
		jid9=`sbatch --dependency=afterok:$jid8 09_filterMaxee.sb | cut -d" " -f 4`
		echo "$jid9: I will filter your reads using -fastq_maxee $max_Eerr."
else
		jid9=`sbatch --dependency=afterok:$jid7 09_filterMaxee.sb | cut -d" " -f 4`
		echo "$jid9: I will filter your reads using -fastq_maxee $max_Eerr."
fi

jid10=`sbatch --dependency=afterok:$jid9 10_dereplicateReads-usearch.sb | cut -d" " -f 4`
echo -e "$jid10: I will dereplicate (i.e. collapse) identical reads."

# Clustering OTUs, ASVs, and/or SWARMs 
if [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jid11=`sbatch --dependency=afterok:$jid10 11_clusterUPARSE.sb | cut -d" " -f 4`
		echo "$jid11: I will cluster OTUs with the UPARSE algorithm in usearch11."
		jid12=`sbatch --dependency=afterok:$jid10 12_clusterUNOISE.sb | cut -d" " -f 4`
		echo "$jid12: I will generate ASVs with the unoise3 algorithm in usearch11."
		jid13=`sbatch --dependency=afterok:$jid10 13_clusterSWARM.sb | cut -d" " -f 4`
		echo "$jid13: I will generate swarms with the SWARM algorithm in swarm3."
elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jid11=`sbatch --dependency=afterok:$jid10 11_clusterUPARSE.sb | cut -d" " -f 4`
		echo "$jid11: I will cluster OTUs with the UPARSE algorithm in usearch11."
		jid13=`sbatch --dependency=afterok:$jid10 13_clusterSWARM.sb | cut -d" " -f 4`
		echo "$jid13: I will generate swarms with the SWARM algorithm in swarm3."
elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jid12=`sbatch --dependency=afterok:$jid10 12_clusterUNOISE.sb | cut -d" " -f 4`
		echo "$jid12: I will generate ASVs with the unoise3 algorithm in usearch11."		
		jid13=`sbatch --dependency=afterok:$jid10 13_clusterSWARM.sb | cut -d" " -f 4`
		echo "$jid13: I will generate swarms with the SWARM algorithm in swarm3."
elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jid13=`sbatch --dependency=afterok:$jid10 13_clusterSWARM.sb | cut -d" " -f 4`
		echo "$jid13: I will generate swarms with the SWARM algorithm in swarm3."
elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "no" ]]; then
		jid11=`sbatch --dependency=afterok:$jid10 11_clusterUPARSE.sb | cut -d" " -f 4`
		echo "$jid11: I will cluster OTUs with the UPARSE algorithm in usearch11."
		jid12=`sbatch --dependency=afterok:$jid10 12_clusterUNOISE.sb | cut -d" " -f 4`
		echo "$jid12: I will generate ASVs with the unoise3 algorithm in usearch11."
elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "no" ]]; then
		jid12=`sbatch --dependency=afterok:$jid10 12_clusterUNOISE.sb | cut -d" " -f 4`
		echo "$jid12: I will generate ASVs with the unoise3 algorithm in usearch11."	
elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "no" ]]; then
		jid11=`sbatch --dependency=afterok:$jid10 11_clusterUPARSE.sb | cut -d" " -f 4`
		echo "$jid11: I will cluster OTUs with the UPARSE algorithm in usearch11."
fi


echo -e "\n========== Assigning taxonomies ==========\n" 
if [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "yes" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC1=`sbatch --dependency=afterok:$jid11 14.1_AssTaxEuk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using the UNITE eukaryote database with CONSTAX."
		jidC2=`sbatch --dependency=afterok:$jid12 14.2_AssTaxEuk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using the UNITE eukaryote database with CONSTAX."
		jidC3=`sbatch --dependency=afterok:$jid13 14.3_AssTaxEuk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using the UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC1=`sbatch --dependency=afterok:$jid11 15.1_AssTaxProk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using SILVA prokaryotic database with CONSTAX."
		jidC2=`sbatch --dependency=afterok:$jid12 15.2_AssTaxProk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using SILVA prokaryotic database with CONSTAX."
		jidC3=`sbatch --dependency=afterok:$jid13 15.3_AssTaxProk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using SILVA prokaryotic database with CONSTAX."
	fi
elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "yes" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC1=`sbatch --dependency=afterok:$jid11 14.1_AssTaxEuk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using the UNITE eukaryote database with CONSTAX."
		jidC3=`sbatch --dependency=afterok:$jid13 14.3_AssTaxEuk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using the UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC1=`sbatch --dependency=afterok:$jid11 15.1_AssTaxProk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using SILVA prokaryotic database with CONSTAX."
		jidC3=`sbatch --dependency=afterok:$jid13 15.3_AssTaxProk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using SILVA prokaryotic database with CONSTAX."
	fi	
elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "yes" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC2=`sbatch --dependency=afterok:$jid12 14.2_AssTaxEuk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using the UNITE eukaryote database with CONSTAX."
		jidC3=`sbatch --dependency=afterok:$jid13 14.3_AssTaxEuk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using the UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC2=`sbatch --dependency=afterok:$jid12 15.2_AssTaxProk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using SILVA prokaryotic database with CONSTAX."
		jidC3=`sbatch --dependency=afterok:$jid13 15.3_AssTaxProk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using SILVA prokaryotic database with CONSTAX."
	fi
elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "yes" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC3=`sbatch --dependency=afterok:$jid13 14.3_AssTaxEuk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using the UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC3=`sbatch --dependency=afterok:$jid13 15.3_AssTaxProk-SWARM-constax.sb | cut -d" " -f 4`
		echo "$jidC3: I will assign taxonomy to SWARMs representative sequences using SILVA prokaryotic database with CONSTAX."
	fi
elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "no" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC1=`sbatch --dependency=afterok:$jid11 14.1_AssTaxEuk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using the UNITE eukaryote database with CONSTAX."
		jidC2=`sbatch --dependency=afterok:$jid12 14.2_AssTaxEuk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using the UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC1=`sbatch --dependency=afterok:$jid11 15.1_AssTaxProk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using SILVA prokaryotic database with CONSTAX."
		jidC2=`sbatch --dependency=afterok:$jid12 15.2_AssTaxProk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using SILVA prokaryotic database with CONSTAX."
	fi
elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "no" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC2=`sbatch --dependency=afterok:$jid12 14.2_AssTaxEuk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC2=`sbatch --dependency=afterok:$jid12 15.2_AssTaxProk-ASV-constax.sb | cut -d" " -f 4`
		echo "$jidC2: I will assign taxonomy to ASV sequences using SILVA prokaryotic database with CONSTAX."
	fi	
elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "no" ]]; then
	if [[ "$DNAmarker" == "ITS" ]] || [[ "$DNAmarker" == "18S" ]]; then	
		jidC1=`sbatch --dependency=afterok:$jid11 14.1_AssTaxEuk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using UNITE eukaryote database with CONSTAX."
	elif [[ "$DNAmarker" == "16S" ]]; then
		jidC1=`sbatch --dependency=afterok:$jid11 15.1_AssTaxProk-OTU-constax.sb | cut -d" " -f 4`
		echo "$jidC1: I will assign taxonomy to OTU sequences using SILVA prokaryotic database with CONSTAX."
	fi
fi


echo -e "\n========== Save and zip all results and reports ==========\n" 
if [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC1:$jidC2:$jidC3 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."

elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC1:$jidC3 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."

elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC2:$jidC3 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."

elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "yes" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC3 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."

elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "no" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC1:$jidC2 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."

elif [[ "$cluster_otu" == "no" ]] && [[ "$cluster_asv" == "yes" ]] && [[ "$cluster_swarm" == "no" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC2 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."

elif [[ "$cluster_otu" == "yes" ]] && [[ "$cluster_asv" == "no" ]] && [[ "$cluster_swarm" == "no" ]]; then
		jidC4=`sbatch --dependency=afterok:$jidC1 16_getResults-bash.sb | cut -d" " -f 4`
		echo "$jidC4: This will be fast. I will organize all the results for you."
fi

echo -e "\n========== These below are the submitted sabtch... ==========\n" 
echo -e "\n `sq` \n"

echo -e "\n========== 'This is the end, my friend'... Now, be patient, you have to wait a bit... ==========\n"
