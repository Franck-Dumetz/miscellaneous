#!/usr/bin/env python

#Summarising many Hisat2 --new-summary output

import os
import csv

# Prompt user for the input directory and output file name
input_dir = input("Enter the directory containing HISAT2 summary files: ").strip()
output_file = input("Enter the name of the output CSV file: ").strip()

# Prepare the CSV header
header = [
    "File", "Total Pairs", "Aligned 0 Times (Concordant/Discordant)",
    "Aligned Concordantly 1 Time", "Aligned Concordantly >1 Times",
    "Aligned Discordantly 1 Time", "Total Unpaired Reads",
    "Aligned 0 Times (Unpaired)", "Aligned 1 Time (Unpaired)",
    "Aligned >1 Times (Unpaired)", "Overall Alignment Rate"
]

# Function to extract data from a summary file
def parse_hisat2_summary(file_path):
    data = {}
    with open(file_path, 'r') as f:
        for line in f:
            if "Total pairs:" in line:
                data["Total Pairs"] = line.split(":")[1].strip()
            elif "Aligned concordantly or discordantly 0 time:" in line:
                data["Aligned 0 Times (Concordant/Discordant)"] = line.split(":")[1].split("(")[0].strip()
            elif "Aligned concordantly 1 time:" in line:
                data["Aligned Concordantly 1 Time"] = line.split(":")[1].split("(")[0].strip()
            elif "Aligned concordantly >1 times:" in line:
                data["Aligned Concordantly >1 Times"] = line.split(":")[1].split("(")[0].strip()
            elif "Aligned discordantly 1 time:" in line:
                data["Aligned Discordantly 1 Time"] = line.split(":")[1].split("(")[0].strip()
            elif "Total unpaired reads:" in line:
                data["Total Unpaired Reads"] = line.split(":")[1].strip()
            elif "Aligned 0 time:" in line and "(Unpaired)" not in data:
                data["Aligned 0 Times (Unpaired)"] = line.split(":")[1].split("(")[0].strip()
            elif "Aligned 1 time:" in line and "(Unpaired)" not in data:
                data["Aligned 1 Time (Unpaired)"] = line.split(":")[1].split("(")[0].strip()
            elif "Aligned >1 times:" in line and "(Unpaired)" not in data:
                data["Aligned >1 Times (Unpaired)"] = line.split(":")[1].split("(")[0].strip()
            elif "Overall alignment rate:" in line:
                data["Overall Alignment Rate"] = line.split(":")[1].strip().replace('%', '')
    return data

# Create the summary CSV
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)

    # Process each summary file in the directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):
            file_path = os.path.join(input_dir, filename)
            summary_data = parse_hisat2_summary(file_path)
            row = [
                filename,
                summary_data.get("Total Pairs", "N/A"),
                summary_data.get("Aligned 0 Times (Concordant/Discordant)", "N/A"),
                summary_data.get("Aligned Concordantly 1 Time", "N/A"),
                summary_data.get("Aligned Concordantly >1 Times", "N/A"),
                summary_data.get("Aligned Discordantly 1 Time", "N/A"),
                summary_data.get("Total Unpaired Reads", "N/A"),
                summary_data.get("Aligned 0 Times (Unpaired)", "N/A"),
                summary_data.get("Aligned 1 Time (Unpaired)", "N/A"),
                summary_data.get("Aligned >1 Times (Unpaired)", "N/A"),
                summary_data.get("Overall Alignment Rate", "N/A")
            ]
            writer.writerow(row)

print(f"Summary report saved to {output_file}")
