import sys
import yaml
import pandas as pd
from pysradb.sraweb import SRAweb

def fetch_sra_metadata(geo_id):
    """Fetch SRA metadata for a given GEO dataset."""
    db = SRAweb()
    
    print(f"Fetching metadata for GEO dataset {geo_id}...")
    metadata = db.gsm_to_srr(geo_id)

    if metadata.empty:
        raise ValueError(f"No SRA records found for {geo_id}")

    # Extract relevant columns
    metadata = metadata[["experiment_alias", "run_accession"]]

    # Map sample names to SRR IDs (list per GSM to handle multiple runs per sample)
    samples = metadata.groupby("experiment_alias")["run_accession"].apply(list).to_dict()

    return samples

def create_config(geo_id, output_file="config.yaml"):
    """Generate a config.yaml file for the Snakemake pipeline."""
    print(f"Generating config file for {geo_id}...")

    # Fetch metadata
    samples = fetch_sra_metadata(geo_id)

    # Reference genome (modify as needed)
    reference_genome = "references/refdata-gex-GRCh38-2024-A"

    # Build config structure
    config_data = {
        "datasets": {
            geo_id: {
                "samples": samples,
                "reference": reference_genome
            }
        }
    }

    # Write to YAML file
    with open(output_file, "w") as yaml_file:
        yaml.dump(config_data, yaml_file, default_flow_style=False)

    print(f"Config file saved: {output_file}")

# Run the script with user input
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_config.py GSEXXXXX")
        sys.exit(1)

    geo_id = sys.argv[1]
    create_config(geo_id)
