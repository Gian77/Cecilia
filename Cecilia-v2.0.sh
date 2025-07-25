#!/bin/bash

cat << "EOF"
   ______          _ ___      
  / ____/__  _____(_) (_)___ _
 / /   / _ \/ ___/ / / / __ `/
/ /___/  __/ /__/ / / / /_/ / 
\____/\___/\___/_/_/_/\__,_/  
                                    
EOF

echo -e "\nHello there, I am Cecilia, your usearCh basEd ampliCon pIpeLine for Illumina dAta."
echo -e "Cecilia v.1.0 by Gian M. N. Benucci, Ph.D."
echo -e "email: benucci[at]msu[dot]edu\nAugust 22, 2022\n"
echo -e "This pipeline is based upon work supported by the Great Lakes Bioenergy Research\nCenter, U.S. Department of Energy, Office of Science, Office of Biological and Environmental\nResearch under award DE-SC0018409\n"

# Load configuration
source ./config.yaml

# Validate critical variables
[[ -z "$project_dir" ]] && { echo "Error: project_dir not set in config.yaml."; exit 1; }

cd "$project_dir/rawdata/"

echo -e "\n========== Comparing md5sum codes... ==========\n"
if [[ -f "md5.txt" ]]; then
    echo -e "\nAn md5 file exists. Checking file integrity...\n"
    md5sum md5* --check > tested_md5.results
    cat tested_md5.results
    resmd5=$(cut -f 2 -d" " tested_md5.results | uniq)

    if [[ "$resmd5" == "OK" ]]; then
        echo -e "\nGood news! Files match the source.\n"
    else
        echo -e "\nError: Files differ from the original source. Please re-download.\n"
        exit 1
    fi
else
    echo "No md5 file found. You should provide one and start over."
    exit 1
fi

cd "$project_dir/code/"
echo -e "\n========== What I will do for you? Please see below... ==========\n"

echo -e "\n========== Prefiltering ==========\n"
jid1=$(sbatch 01_decompress-bash.sb | cut -d" " -f 4)
echo "$jid1: Decompressing raw reads."

jid2=$(sbatch --dependency=afterok:$jid1 02_qualityCheck-fastqc.sb | cut -d" " -f 4)
echo "$jid2: Checking quality and generating stats."

jid3=$(sbatch --dependency=afterok:$jid1:$jid2 03_removePhix-usearch.sb | cut -d" " -f 4)
echo "$jid3: Removing PhiX reads in USEARCH."

if [[ "$assemble" == "yes" ]]; then
    jid4=$(sbatch --dependency=afterok:$jid3 04_readAssembly-usearch.sb | cut -d" " -f 4)
    echo "$jid4: Assembling reads in USEARCH."
fi

if [[ "$rename" == "yes" ]]; then
    if [[ "$assemble" == "yes" ]]; then
        jid5=$(sbatch --dependency=afterok:$jid4 05_readRename.sb | cut -d" " -f 4)
    else
        jid5=$(sbatch --dependency=afterok:$jid3 05_readRename.sb | cut -d" " -f 4)
    fi
    echo "$jid5: Renaming reads with 'Sample...' prefix."
fi

if [[ "$assemble" == "no" && "$rename" == "no" ]]; then
    jid6=$(sbatch --dependency=afterok:$jid3 06_primerStrip-cutadapt.sb | cut -d" " -f 4)
elif [[ "$assemble" == "yes" && "$rename" == "no" ]]; then
    jid6=$(sbatch --dependency=afterok:$jid4 06_primerStrip-cutadapt.sb | cut -d" " -f 4)
else
    jid6=$(sbatch --dependency=afterok:$jid5 06_primerStrip-cutadapt.sb | cut -d" " -f 4)
fi
echo "$jid6: Removing primers/adapters with Cutadapt."

jid7=$(sbatch --dependency=afterok:$jid6 07_readTrimDemux-usearch.sb | cut -d" " -f 4)
echo "$jid7: Trimming and demultiplexing reads."

if [[ "$test_length" == "yes" ]]; then
    jid8=$(sbatch --dependency=afterok:$jid7 08_testLength-usearch.sb | cut -d" " -f 4)
    echo "$jid8: Testing optimal read length."
fi

echo -e "\n========== Clustering ==========\n"
if [[ "$test_length" == "yes" ]]; then
    jid9=$(sbatch --dependency=afterok:$jid8 09_filterMaxee.sb | cut -d" " -f 4)
else
    jid9=$(sbatch --dependency=afterok:$jid7 09_filterMaxee.sb | cut -d" " -f 4)
fi
echo "$jid9: Filtering reads with maxee=$max_Eerr."

jid10=$(sbatch --dependency=afterok:$jid9 10_dereplicateReads-usearch.sb | cut -d" " -f 4)
echo "$jid10: Dereplicating reads."

if [[ "$cluster_otu" == "yes" ]]; then
    jid11=$(sbatch --dependency=afterok:$jid10 11_clusterUPARSE.sb | cut -d" " -f 4)
    echo "$jid11: Clustering OTUs with UPARSE."
fi

if [[ "$cluster_asv" == "yes" ]]; then
    jid12=$(sbatch --dependency=afterok:$jid10 12_clusterUNOISE.sb | cut -d" " -f 4)
    echo "$jid12: Generating ASVs with UNOISE."
    if [[ "$cluster_asv_to_otu" == "yes" ]]; then
        jid14=$(sbatch --dependency=afterok:$jid12 13_clusterASV.sb | cut -d" " -f 4)
        echo "$jid14: Clustering ASVs into OTUs."
    fi
fi

if [[ "$cluster_swarm" == "yes" ]]; then
    jid13=$(sbatch --dependency=afterok:$jid10 15_clusterSWARM.sb | cut -d" " -f 4)
    echo "$jid13: Clustering with SWARM."
fi

if [[ "$closed_ref_otu" == "yes" ]]; then
    jid15=$(sbatch --dependency=afterok:$jid10 14_closedRef_OTU.sb | cut -d" " -f 4)
    echo "$jid15: Performing closed-reference OTU clustering."
fi

echo -e "\n========== Assigning taxonomies ==========\n"

deps=()
[[ "$constax_taxonomy" == "yes" ]] && deps+=("$jid11" "$jid12" "$jid13" "$jid14")
if (( ${#deps[@]} )); then
    dep_string=$(IFS=":"; echo "${deps[*]}")
    jid16=$(sbatch --dependency=afterok:$dep_string 16_taxonomyConstax.sb | cut -d" " -f 4)
    echo "$jid16: Assigning taxonomy using CONSTAX2."
fi

deps=()
[[ "$sintax_taxonomy" == "yes" ]] && deps+=("$jid11" "$jid12" "$jid13" "$jid14")
if (( ${#deps[@]} )); then
    dep_string=$(IFS=":"; echo "${deps[*]}")
    jid17=$(sbatch --dependency=afterok:$dep_string 17_taxonomySintax.sb | cut -d" " -f 4)
    echo "$jid17: Assigning taxonomy using USEARCH SINTAX."
fi

echo -e "\n========== Packaging results ==========\n"
all_jids=("${jid1:-}" "${jid2:-}" "${jid3:-}" "${jid4:-}" "${jid5:-}" \
          "${jid6:-}" "${jid7:-}" "${jid8:-}" "${jid9:-}" "${jid10:-}" \
          "${jid11:-}" "${jid12:-}" "${jid13:-}" "${jid14:-}" "${jid15:-}" \
          "${jid16:-}" "${jid17:-}")

final_deps=()
for jid in "${all_jids[@]}"; do
    [[ -n "$jid" ]] && final_deps+=("$jid")
done

if (( ${#final_deps[@]} > 0 )); then
    dep_string=$(IFS=:; echo "${final_deps[*]}")
    jid18=$(sbatch --dependency=afterok:$dep_string 18_getResults-bash.sb | cut -d" " -f 4)
    echo "$jid18: Packaging all results after jobs: $dep_string"
else
    jid18=$(sbatch 18_getResults-bash.sb | cut -d" " -f 4)
    echo "$jid18: No prior jobs found. Packaging results now."
fi

echo -e "\n========== Job summary =========="
for i in {1..17}; do
    var="jid$i"
    [[ -n "${!var:-}" ]] && echo "$var: ${!var}"
done
[[ -n "$jid18" ]] && echo "jid18: $jid18"

echo -e "\n========== 'This is the end, my friend'... Now, be patient. =========="
