#!/usr/bin/env bash
# =============================================================================
# summarizeCounts.sh
#
# Joins the per-stage read counts that the pipeline already produces into two
# tidy TSVs:
#
#   read_tracking.tsv  one row per sample — read attrition through the
#                      per-sample stages (raw -> nophix -> assembled ->
#                      stripped). These stages are unambiguous: one value
#                      per sample.
#
#   otu_summary.tsv    long format, one row per (sample, length, method) —
#                      filtered reads, observed richness, and reads mapped
#                      back to the reference set. These stages fan out over
#                      truncation length and clustering method, so they
#                      cannot share a row with the above.
#
# Raw counts are recomputed here from the raw FASTQs rather than emitted by
# DECOMPRESS on purpose: adding an output to DECOMPRESS would change its
# script hash and invalidate every downstream task in the resume cache.
#
# Usage:
#   summarizeCounts.sh <rawdir> <outprefix>
# Expects in the working directory (staged by the process):
#   file_mapping.txt          (optional; absent when rename = false)
#   all.nophix.counts         (optional)
#   all.assembled.counts      (optional; absent when assemble = false)
#   pooled.counts
#   filtered_*bp.fastq        (any number)
#   otutable_*.txt / asvtable_*.txt  (any number)
# =============================================================================
set -euo pipefail

rawdir="${1:?usage: summarizeCounts.sh <rawdir> <outprefix>}"
prefix="${2:-.}"

track="${prefix}/read_tracking.tsv"
otusum="${prefix}/otu_summary.tsv"

# -----------------------------------------------------------------------------
# Count reads in a FASTQ regardless of compression. Prints the record count.
# -----------------------------------------------------------------------------
count_fastq() {
    local f="$1" n=0
    case "$f" in
        *.gz)  n=$(gzip -cd  -- "$f" | wc -l) ;;
        *.bz2) n=$(bzip2 -cd -- "$f" | wc -l) ;;
        *)     n=$(wc -l < "$f") ;;
    esac
    echo $(( n / 4 ))
}

# -----------------------------------------------------------------------------
# Strip the pipeline's filename decorations down to the original sample stem.
#   004_S4_L001_R1_001_nophix.fastq -> 004_S4_L001
#   004_S4_L001_assembled.fastq     -> 004_S4_L001
#   004_S4_L001_R1_001.fastq.gz     -> 004_S4_L001
# -----------------------------------------------------------------------------
stem() {
    printf '%s\n' "$1" \
        | sed -E 's/\.(fastq|fq)(\.gz|\.bz2)?$//' \
        | sed -E 's/_nophix$//' \
        | sed -E 's/_assembled$//' \
        | sed -E 's/_dc$//' \
        | sed -E 's/_R[12](_[0-9]+)?$//' \
        | sed -E 's/_[12]_[12]$//'
}

# =============================================================================
# read_tracking.tsv
# =============================================================================

# --- sample <-> original stem, from file_mapping.txt when renaming is on -----
: > .map_stem_to_sample
if [[ -s file_mapping.txt ]]; then
    while IFS=$'\t' read -r orig new; do
        [[ -n "${orig:-}" ]] || continue
        s=$(stem "$orig")
        samp="${new%.fastq}"
        printf '%s\t%s\n' "$s" "$samp" >> .map_stem_to_sample
    done < file_mapping.txt
fi

# --- raw counts, one entry per original stem (R1 only = read pairs) ---------
: > .raw_counts
shopt -s nullglob
raw_seen=""
for f in "$rawdir"/*; do
    b=$(basename "$f")
    case "$b" in
        *.fastq|*.fq|*.fastq.gz|*.fq.gz|*.fastq.bz2|*.fq.bz2) ;;
        *) continue ;;
    esac
    # For paired data count R1 only; R2 mirrors it and would double the total.
    case "$b" in
        *_R2_*|*_R2.*|*_2_2.*) continue ;;
    esac
    s=$(stem "$b")
    case " $raw_seen " in *" $s "*) continue ;; esac
    raw_seen="$raw_seen $s"
    printf '%s\t%s\n' "$s" "$(count_fastq "$f")" >> .raw_counts
done

# --- nophix / assembled counts, keyed by stem ------------------------------
parse_counts() {   # <file> -> "stem <tab> count", R1/first entry per stem
    local src="$1"
    [[ -s "$src" ]] || return 0
    sed -E 's/[[:space:]]*:[[:space:]]*/\t/' "$src" \
        | awk -F'\t' 'NF==2 {print $1"\t"$2}' \
        | while IFS=$'\t' read -r fname n; do
              case "$fname" in *_R2_*|*_R2.*) continue ;; esac
              printf '%s\t%s\n' "$(stem "$fname")" "$n"
          done \
        | awk -F'\t' '!seen[$1]++'
}
parse_counts all.nophix.counts    > .nophix_counts    || true
parse_counts all.assembled.counts > .assembled_counts || true

# --- write the table -------------------------------------------------------
{
    printf 'sample\torig_file\traw\tnophix\tassembled\tstripped\n'

    awk -F'\t' '
        FNR==NR && FILENAME==map      { samp[$1]=$2; order[++n]=$1; next }
        FILENAME==rawf                { raw[$1]=$2;      next }
        FILENAME==nopf                { nop[$1]=$2;      next }
        FILENAME==asmf                { asm[$1]=$2;      next }
        FILENAME==poolf               { split($0,a,/[[:space:]]*:[[:space:]]*/);
                                        pooled[a[1]]=a[2]; next }
        END {
            if (n == 0) {
                # rename = false: fall back to the stems seen in raw counts
                for (s in raw) order[++n]=s
            }
            for (i=1; i<=n; i++) {
                s = order[i]
                sm = (s in samp) ? samp[s] : s
                printf "%s\t%s\t%s\t%s\t%s\t%s\n", sm, s,
                    (s  in raw)    ? raw[s]     : "NA",
                    (s  in nop)    ? nop[s]     : "NA",
                    (s  in asm)    ? asm[s]     : "NA",
                    (sm in pooled) ? pooled[sm] : "NA"
            }
        }
    ' map=.map_stem_to_sample rawf=.raw_counts nopf=.nophix_counts \
      asmf=.assembled_counts poolf=pooled.counts \
      .map_stem_to_sample .raw_counts .nophix_counts .assembled_counts pooled.counts
} > "$track"

# =============================================================================
# otu_summary.tsv
# =============================================================================

# Authoritative sample roster, in pipeline order, taken from the table above.
tail -n +2 "$track" | cut -f1 > .sample_list

# Per-sample pooled (primer-stripped) counts, keyed by sampleNNN. This is the
# denominator for every method that runs `usearch -otutab <pooled.fastq>`:
# those map the whole pooled read set back to the reference sequences, not the
# length-filtered subset. Using the filtered count there yields >100%.
sed -E 's/[[:space:]]*:[[:space:]]*/\t/' pooled.counts > .pooled_counts

# --- per-sample filtered counts, per truncation length ---------------------
: > .filtered_counts
for f in filtered_*bp.fastq; do
    len=$(printf '%s\n' "$f" | sed -E 's/^filtered_([0-9]+)bp\.fastq$/\1/')
    awk 'NR%4==1' "$f" \
        | sed 's/^@//' \
        | awk -F'.' -v L="$len" '{c[$1]++} END {for (s in c) printf "%s\t%s\t%s\n", s, L, c[s]}' \
        >> .filtered_counts
done

# --- per-sample richness and mapped reads, per table -----------------------
{
    printf 'sample\tlen\tmethod\tfiltered\totus_obs\treads_mapped\tbasis\tbasis_reads\tpct_mapped\n'

    for t in otutable_*.txt asvtable_*.txt; do
        b=$(basename "$t" .txt)
        # otutable_UPARSE_250bp / asvtable_UNOISE_250bp
        # otutable_250bp_closedRef / otutable_asv_250bp_to97otus
        len=$(printf '%s\n' "$b" | grep -oE '[0-9]+bp' | head -1 | tr -d 'bp')
        case "$b" in
            *_closedRef)   method="closedRef" ;;
            *_to97otus)    method="ASV_TO_OTU" ;;
            otutable_*|asvtable_*)
                method=$(printf '%s\n' "$b" | sed -E 's/^(otu|asv)table_//; s/_[0-9]+bp.*$//') ;;
            *) method="$b" ;;
        esac
        [[ -n "$len" ]] || continue

        # Every method now builds its table with `usearch -otutab pooled.fastq`,
        # so each read carries its own sampleNNN label and the per-sample
        # columns are real. Denominator = that sample's pooled (stripped) count.
        basis="pooled"

        # usearch -otutab drops a sample's column entirely when nothing maps.
        # Iterate over the full sample roster instead of the table's columns so
        # every (sample, len, method) cell exists, with an explicit 0 rather
        # than a silently absent row.
        awk -F'\t' -v L="$len" -v M="$method" -v B="$basis" -v FC=.filtered_counts \
                   -v PC=.pooled_counts -v SL=.sample_list '
            BEGIN {
                while ((getline line < FC) > 0) {
                    split(line, p, "\t")
                    if (p[2] == L) filt[p[1]] = p[3]
                }
                while ((getline line < PC) > 0) {
                    split(line, q, "\t")
                    if (q[1] != "") pool[q[1]] = q[2]
                }
                while ((getline s < SL) > 0) if (s != "") roster[++nr] = s
            }
            NR==1 { for (j=2; j<=NF; j++) { name[j]=$j; col[$j]=j } ncol=NF; next }
            {
                for (j=2; j<=ncol; j++) {
                    sum[name[j]] += $j
                    if ($j+0 > 0) nz[name[j]]++
                }
            }
            END {
                for (i=1; i<=nr; i++) {
                    s = roster[i]
                    f = (s in filt) ? filt[s] : 0
                    d = (s in pool) ? pool[s] : 0
                    m = sum[s]+0
                    pct = (d > 0) ? sprintf("%.2f", 100*m/d) : "NA"
                    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                        s, L, M, (f>0 ? f : "NA"), nz[s]+0, m,
                        B, (d>0 ? d : "NA"), pct
                }
            }
        ' "$t"
    done
} | { read -r hdr; printf '%s\n' "$hdr"; sort -k1,1 -k2,2n -k3,3; } > "$otusum"

rm -f .map_stem_to_sample .raw_counts .nophix_counts .assembled_counts \
      .filtered_counts .sample_list .pooled_counts

echo "Wrote $track and $otusum"
wc -l "$track" "$otusum"
