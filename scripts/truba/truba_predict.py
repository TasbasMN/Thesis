import argparse
import time
import json
import pandas as pd
from scripts.utils import *

# cli arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", type=int, help="Start index")
parser.add_argument("-e", type=int, help="End index")
parser.add_argument("--vcf", help="VCF file to analyze")
parser.add_argument("--output_dir", default="./results", help="Output directory for runtime file, defaults to ./results")
args = parser.parse_args()

start = args.s
end = args.e
vcf_file = args.vcf
output_dir = args.output_dir

# Measure the start time
start_time = time.time()

# import df
df = pd.read_csv(vcf_file, sep="\t", header=None)
df.columns = ["chr", "pos", "unknown", "ref", "alt"]
df = df[start:end]

# augment vcf 
start_augment_vcf = time.time()
df = augment_cancer_vcf(df)
end_augment_vcf = time.time()

# finding matches
start_find_matches = time.time()
results = find_matches_for_vcfs(df)
results_mutated = find_matches_for_vcfs(df, mutated=True)
end_find_matches = time.time()

# adding labels
results["is_mutated"] = 0
results_mutated["is_mutated"] = 1
final = pd.concat([results, results_mutated], axis=0)

# Generate a unique identifier for the output file
timestamp = time.strftime("%Y%m%d-%H%M%S")
output_file = f"{output_dir}/{timestamp}_runtime_report.json"

final.to_csv(f"results/7_results_{start}_{end}.csv", index=False)

# Calculate the runtimes
augment_vcf_runtime = end_augment_vcf - start_augment_vcf
find_matches_runtime = end_find_matches - start_find_matches
total_runtime = time.time() - start_time

# Convert start time and end time to human-readable format
start_time_readable = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
end_time_readable = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# Write the runtimes and additional job details to the output file
data = {
    "Job Details": {
        "VCF file": vcf_file,
        "Start Row": start,
        "End Row": end,
        "Start Time": start_time_readable,
        "End Time": end_time_readable
    },
    "Runtimes": {
        "Augment VCF": augment_vcf_runtime,
        "Find Matches": find_matches_runtime,
        "Total Runtime": total_runtime
    }
}

with open(output_file, "w") as f:
    json.dump(data, f, indent=4)


