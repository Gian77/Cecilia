#!/bin/bash

cat << "EOF"
   ______          _ ___      
  / ____/__  _____(_) (_)___ _
 / /   / _ \/ ___/ / / / __ `/
/ /___/  __/ /__/ / / / /_/ / 
\____/\___/\___/_/_/_/\__,_/  
                                    
EOF

echo -e "\nHi! I am Cecilia, your usearCh basEd ampliCon pIpeLine for Illumina dAta."
echo -e "Cecilia v.2.0 by Gian M. N. Benucci, Ph.D."
echo -e "email: benucci[at]msu[dot]edu\nJuly 25, 2025\n"
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
# -----------------------------------------------------------------------------
# DECOMPRESS fastq
# -----------------------------------------------------------------------------
echo -e "\n========== Prefiltering ==========\n"
jid1=$(sbatch 01_decompress-bash.sb | cut -d" " -f 4)
echo "$jid1: Decompressing raw reads."

# -----------------------------------------------------------------------------
# QUALITY REPORTs
# -----------------------------------------------------------------------------
jid2=$(sbatch --dependency=afterok:$jid1 02_qualityCheck-fastqc.sb | cut -d" " -f 4)
echo "$jid2: Checking quality and generating stats."
# -----------------------------------------------------------------------------
# REMOVE PHIX
# -----------------------------------------------------------------------------
jid3=$(sbatch --dependency=afterok:$jid1:$jid2 03_removePhix-usearch.sb | cut -d" " -f 4)
echo "$jid3: Removing PhiX reads in USEARCH."

# -----------------------------------------------------------------------------
# ASSEMBLE
# -----------------------------------------------------------------------------
if [[ "$assemble" == "yes" ]]; then
    jid4=$(sbatch --dependency=afterok:$jid3 04_readAssembly-usearch.sb | cut -d" " -f 4)
    echo "$jid4: Assembling reads in USEARCH."
fi

# -----------------------------------------------------------------------------
# ASVs and OTUs
# -----------------------------------------------------------------------------
if [[ "$rename" == "yes" ]]; then
    if [[ "$assemble" == "yes" ]]; then
        jid5=$(sbatch --dependency=afterok:$jid4 05_readRename.sb | cut -d" " -f 4)
    else
        jid5=$(sbatch --dependency=afterok:$jid3 05_readRename.sb | cut -d" " -f 4)
    fi
    echo "$jid5: Renaming reads with 'Sample...' prefix."
fi

# -----------------------------------------------------------------------------
# ASVs and OTUs
# -----------------------------------------------------------------------------
if [[ "$assemble" == "no" && "$rename" == "no" ]]; then
    jid6=$(sbatch --dependency=afterok:$jid3 06_primerStrip-cutadapt.sb | cut -d" " -f 4)
elif [[ "$assemble" == "yes" && "$rename" == "no" ]]; then
    jid6=$(sbatch --dependency=afterok:$jid4 06_primerStrip-cutadapt.sb | cut -d" " -f 4)
else
    jid6=$(sbatch --dependency=afterok:$jid5 06_primerStrip-cutadapt.sb | cut -d" " -f 4)
fi
echo "$jid6: Removing primers/adapters with Cutadapt."

# -----------------------------------------------------------------------------
# ASVs and OTUs
# -----------------------------------------------------------------------------
jid7=$(sbatch --dependency=afterok:$jid6 07_readEEstats-usearch.sb | cut -d" " -f 4)
echo "$jid7: Generating max-expetced-erros statsitics."

# -----------------------------------------------------------------------------
# TEST for BEST LENGTH
# -----------------------------------------------------------------------------
if [[ "$test_length" == "yes" ]]; then
    jid8=$(sbatch --dependency=afterok:$jid7 08_testLength-usearch.sb | cut -d" " -f 4)
    echo "$jid8: Testing optimal read length."
fi

# -----------------------------------------------------------------------------
# FILTERING BASED ON MAX EXPECTED ERRORs
# -----------------------------------------------------------------------------
echo -e "\n========== FIltering based upon Max Expected Errors ==========\n"
if [[ "$test_length" == "yes" ]]; then
    jid9=$(sbatch --dependency=afterok:$jid8 09_filterMaxee.sb | cut -d" " -f 4)
else
    jid9=$(sbatch --dependency=afterok:$jid7 09_filterMaxee.sb | cut -d" " -f 4)
fi
echo "$jid9: Filtering reads with maxee=$max_Eerr."

# -----------------------------------------------------------------------------
# DEREPLICATING
# -----------------------------------------------------------------------------
jid10=$(sbatch --dependency=afterok:$jid9 10_dereplicateReads-usearch.sb | cut -d" " -f 4)
echo "$jid10: Dereplicating reads."

# -----------------------------------------------------------------------------
# ASVs and OTUs
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# TAXONOMY ASSIGNMENTS
# -----------------------------------------------------------------------------
echo -e "\n========== Assigning taxonomies ==========\n"

# Build dependency string
deps=""
for jid in "$jid11" "$jid12" "$jid13" "$jid14"; do
    if [[ -n "$jid" ]]; then
        deps="${deps:+$deps:}$jid"
    fi
done

if [[ -z "$deps" ]]; then
    echo "No upstream clustering jobs found."
else

    if [[ "$constax_taxonomy" == "yes" ]]; then
        jid16=$(sbatch --dependency=afterok:$deps 16_taxonomyConstax.sb | awk '{print $4}')
        echo "$jid16: Assigning taxonomy using CONSTAX2 (after $deps)."
    else
        echo "No taxonomy CONSTAX assignment selected."
    fi

    if [[ "$sintax_taxonomy" == "yes" ]]; then
        jid17=$(sbatch --dependency=afterok:$deps 17_taxonomySintax.sb | awk '{print $4}')
        echo "$jid17: Assigning taxonomy using USEARCH SINTAX (after $deps)."
    else
        echo "No taxonomy SINTAX assignment selected."
    fi

fi

# -----------------------------------------------------------------------------
# PACKING UP RESUTLS
# -----------------------------------------------------------------------------
echo -e "\n========== Packaging results ==========\n"

# Collect all job dependencies
all_jids=()
for i in {1..17}; do
    var="jid$i"
    [[ -n "${!var:-}" ]] && all_jids+=("${!var}")
done

# Check if any of the result-producing steps were run (10,11,13–17)
result_jobs=("${jid10:-}" "${jid11:-}" "${jid13:-}" "${jid14:-}" "${jid15:-}" "${jid16:-}" "${jid17:-}")
has_results=false
for jid in "${result_jobs[@]}"; do
    [[ -n "$jid" ]] && has_results=true && break
done

# Submit packaging job with appropriate message
if (( ${#all_jids[@]} > 0 )); then
    dep_string=$(IFS=:; echo "${all_jids[*]}")
    jid18=$(sbatch --dependency=afterok:$dep_string 18_getResults-bash.sb | cut -d" " -f 4)
    if ! $has_results; then
        echo "Warning: Packaging triggered, but no result-generating steps (10,11,13–17) detected."
    fi
    echo "$jid18: Packaging all results after jobs: $dep_string"
else
    jid18=$(sbatch 18_getResults-bash.sb | cut -d" " -f 4)
    echo "No prior jobs found. Packaging results immediately — they may be incomplete."
    echo "jid18: $jid18"
fi

echo -e "\n========== 'This is the end, my friend'... Now, be patient. =========="
