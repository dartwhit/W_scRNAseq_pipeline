"""
Rename already-downloaded FASTQ files from the old SRR-based naming scheme
to the new GSM+lane naming scheme expected by the updated Snakefile.

Old format: {SRR}_S1_L001_R{1,2}_001.fastq.gz
New format: {GSM}_S1_L{lane}_R{1,2}_001.fastq.gz

Run from the pipeline root directory:
    python rename_existing_fastqs.py GSE264508
"""

import sys
import os
import yaml

def main(dataset):
    with open("config.yaml") as f:
        config_data = yaml.safe_load(f)

    samples_dict = config_data["datasets"][dataset]["samples"]

    # Normalize to lists
    for gsm in samples_dict:
        if isinstance(samples_dict[gsm], str):
            samples_dict[gsm] = [samples_dict[gsm]]

    fastq_dir = f"datasets/{dataset}/fastq"

    renames = []   # (src, dst) pairs for valid files
    deletes = []   # files to remove (bad/incomplete downloads)

    for gsm, srr_list in samples_dict.items():
        for lane_idx, sra_id in enumerate(srr_list, start=1):
            lane = f"{lane_idx:03d}"
            for read in ["1", "2"]:
                old_name = f"{sra_id}_S1_L001_R{read}_001.fastq.gz"
                new_name = f"{gsm}_S1_L{lane}_R{read}_001.fastq.gz"
                old_path = os.path.join(fastq_dir, old_name)
                new_path = os.path.join(fastq_dir, new_name)

                if os.path.exists(new_path):
                    print(f"  SKIP (already exists): {new_name}")
                    continue

                if os.path.exists(old_path):
                    size = os.path.getsize(old_path)
                    if size == 0:
                        deletes.append(old_path)
                    else:
                        renames.append((old_path, new_path))

        # Also remove any partial single-file downloads for SRRs in this GSM
        for sra_id in srr_list:
            single_file = os.path.join(fastq_dir, f"{sra_id}.fastq.gz")
            if os.path.exists(single_file):
                deletes.append(single_file)

    print(f"\nFiles to rename: {len(renames)}")
    for src, dst in renames:
        print(f"  {os.path.basename(src)} -> {os.path.basename(dst)}")

    print(f"\nFiles to delete (bad/partial downloads): {len(deletes)}")
    for path in deletes:
        print(f"  {os.path.basename(path)}")

    if not renames and not deletes:
        print("\nNothing to do.")
        return

    confirm = input("\nProceed? [y/N] ").strip().lower()
    if confirm != "y":
        print("Aborted.")
        return

    for src, dst in renames:
        os.rename(src, dst)
        print(f"  Renamed: {os.path.basename(src)} -> {os.path.basename(dst)}")

    for path in deletes:
        os.remove(path)
        print(f"  Deleted: {os.path.basename(path)}")

    print("\nDone. Run snakemake to download remaining files and continue the pipeline.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python rename_existing_fastqs.py GSEXXXXX")
        sys.exit(1)
    main(sys.argv[1])
