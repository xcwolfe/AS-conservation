#!/bin/bash

# Input CSV file containing coordinates
input_csv="extracted_coords_IR_events.csv"

# Output file
output_file="phyloP_scores_combined.txt"

# Temporary file for bedGraph output
temp_bedGraph="temp.bedGraph"

# BigWig file
bigwig_file="ce11.phyloP26way.bw"

# Initialize output file with header
echo -e "AS_event_ID\tChr\tStart\tEnd\tAverage_Score" > "$output_file"

# Read each row of the CSV (excluding the header)
tail -n +2 "$input_csv" | while IFS=',' read -r AS_event_ID chr start end; do
    # Remove any surrounding quotes from chr, start and end
    chr=$(echo "$chr" | tr -d '"')
    start=$(echo "$start" | tr -d '"')
    end=$(echo "$end" | tr -d '"')

    # Prefix the chromosome name with "chr"
    chr="chr$chr"

    # Debugging to check inputs
    echo "Processing: AS_event_ID=$AS_event_ID, Chr=$chr, Start=$start, End=$end"

    # Check if start and end are valid integers
    if [[ "$start" =~ ^[0-9]+$ ]] && [[ "$end" =~ ^[0-9]+$ ]]; then
        # Calculate region length
        region_length=$((end - start + 1))

        # Generate bedGraph data for the region
        ./bigWigToBedGraph "$bigwig_file" "$temp_bedGraph" -chrom="$chr" -start="$start" -end="$end"

        # Extract the phyloP scores for the region
        ./bigWigSummary "$bigwig_file" "$chr" "$start" "$end" "$region_length" > temp_summary.txt

        # Format output (single-line scores)
        scores=$(tr '\n' '\t' < temp_summary.txt)

        # Calculate the average score
        scores_array=($scores)
        total=0
        count=0

        # Sum the scores and count the number of valid scores
        for score in "${scores_array[@]}"; do
            if [[ "$score" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                total=$(echo "$total + $score" | bc)
                ((count++))
            fi
        done

        # Calculate the average (if there are valid scores)
        if ((count > 0)); then
            average=$(echo "scale=5; $total / $count" | bc)
        else
            average="No Data"
        fi

        # Append only the average score to the output file
        echo -e "$AS_event_ID\t$chr\t$start\t$end\t$average" >> "$output_file"

        # Clean up temp files
        rm -f "$temp_bedGraph" temp_summary.txt
    else
        echo "Skipping invalid coordinates: $AS_event_ID ($chr, $start, $end)"
    fi
done
