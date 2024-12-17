import os
import sys

def extract_metrics(file_path):
    """
    Extract the second line after the header from a Picard metrics file.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("## METRICS CLASS"):
                next(file)  # Skip header line
                return next(file).strip()  # Extract the metrics line
    return None

def summarize_metrics(directory, output_file):
    """
    Summarize metrics from all files in a directory into a tab-delimited table.
    """
    with open(output_file, 'w') as out:
        out.write("File\tLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\t"
                  "UNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\t"
                  "PERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n")
        
        for file_name in os.listdir(directory):
            if file_name.endswith("_metrics.txt"):
                file_path = os.path.join(directory, file_name)
                metrics_line = extract_metrics(file_path)
                if metrics_line:
                    out.write(f"{file_name}\t{metrics_line}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python summarize_metrics.py <metrics_directory> <output_file>")
        sys.exit(1)

    metrics_dir = sys.argv[1]
    output_summary = sys.argv[2]

    summarize_metrics(metrics_dir, output_summary)
    print(f"Summary saved to {output_summary}")
