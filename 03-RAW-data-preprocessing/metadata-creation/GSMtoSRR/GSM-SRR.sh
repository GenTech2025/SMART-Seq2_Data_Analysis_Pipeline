#!/bin/bash

# Function to get SRR ID from GSM ID
get_srr_id() {
    local gsm_id=$1
    esearch -db sra -query $gsm_id | efetch -format runinfo | grep $gsm_id
}

# List of hardcoded GSM IDs
gsm_ids=("GSM1833279" "GSM1833280" "GSM1833281" "GSM1833282" "GSM1833283" "GSM1833284" "GSM1833285" "GSM1833286"
         "GSM1833287" "GSM1833288" "GSM1833289" "GSM1833290" "GSM1833291" "GSM1833292" "GSM1833294" "GSM1833295"
         "GSM1833296" "GSM1833297" "GSM1833298" "GSM1833300" "GSM1833302" "GSM1833304" "GSM1833307" "GSM1833309"
         "GSM1833312" "GSM1833314" "GSM1833316" "GSM1833318" "GSM1833321" "GSM1833323" "GSM1833324" "GSM1833325"
         "GSM1833326" "GSM1833327" "GSM1833328" "GSM1833329" "GSM1833330" "GSM1833331" "GSM1833332" "GSM1833333"
         "GSM1833334" "GSM1985628" "GSM1985629" "GSM1985630" "GSM1985631" "GSM2151860" "GSM2151861" "GSM2151862")

# Create a CSV file to store the mappings
output_file="GSE71318_gsm_to_srr_mapping.csv"
echo "GSM_ID,SRR_ID" > $output_file

# Log file for errors
error_log="gsm_to_srr_errors.log"
echo "Errors:" > $error_log

# Loop through each hardcoded GSM ID and get the corresponding SRR ID
for gsm_id in "${gsm_ids[@]}"; do
    if [[ -z "$gsm_id" ]]; then
        continue
    fi
    srr_info=$(get_srr_id $gsm_id)
    if [[ -n "$srr_info" ]]; then
        srr_id=$(echo "$srr_info" | cut -d ',' -f 1)
        echo "$gsm_id,$srr_id" >> $output_file
    else
        echo "$gsm_id,NA" >> $output_file
        echo "Failed to map $gsm_id" >> $error_log
    fi
done

echo "Mapping completed. Check the file $output_file for the results and $error_log for any errors."

