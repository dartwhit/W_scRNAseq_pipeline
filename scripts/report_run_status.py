#!/usr/bin/env python3
"""Generate a reusable per-sample stage status report for this Snakemake pipeline.

The report uses:
1) config.yaml sample/lane definitions
2) expected pipeline outputs on disk
3) optional failure hints parsed from .snakemake/log
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import yaml


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def select_dataset(config_data: dict, dataset_arg: str | None) -> str:
    if dataset_arg:
        dataset = dataset_arg
    else:
        dataset = os.getenv("DATASET")

    if not dataset:
        datasets = sorted(config_data.get("datasets", {}).keys())
        if len(datasets) == 1:
            return datasets[0]
        raise ValueError(
            "Dataset is required. Use --dataset or set DATASET in the environment."
        )

    if dataset not in config_data.get("datasets", {}):
        raise ValueError(f"Dataset '{dataset}' not found in config.yaml")
    return dataset


def normalize_samples(dataset_cfg: dict) -> Dict[str, List[str]]:
    samples = dataset_cfg.get("samples", {})
    norm = {}
    for gsm, srrs in samples.items():
        if isinstance(srrs, str):
            norm[gsm] = [srrs]
        elif isinstance(srrs, list):
            norm[gsm] = [str(x) for x in srrs]
        else:
            raise ValueError(f"Sample '{gsm}' has unsupported value type: {type(srrs)}")
    return norm


def stage_status(exists_count: int, expected_count: int) -> str:
    if expected_count == 0:
        return "not_applicable"
    if exists_count == expected_count:
        return "ok"
    if exists_count == 0:
        return "missing"
    return "partial"


def parse_failure_hints(log_dir: Path, dataset: str) -> Tuple[Dict[str, List[str]], List[str]]:
    """Parse Snakemake logs and return per-sample + dataset-level failure hints."""
    sample_hints: Dict[str, List[str]] = defaultdict(list)
    dataset_hints: List[str] = []

    if not log_dir.exists():
        return sample_hints, dataset_hints

    log_files = [p for p in log_dir.rglob("*") if p.is_file()]
    log_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)

    # Keep this bounded but useful for previous runs.
    for log_path in log_files[:20]:
        try:
            text = log_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue

        rule_matches = list(re.finditer(r"Error in rule ([A-Za-z0-9_]+):", text))
        if not rule_matches:
            continue

        for idx, match in enumerate(rule_matches):
            rule = match.group(1)
            block_start = match.end()
            block_end = rule_matches[idx + 1].start() if idx + 1 < len(rule_matches) else len(text)
            block = text[block_start:block_end]

            wildcards_match = re.search(r"wildcards:\s*([^\n]+)", block)
            wildcards = {}
            if wildcards_match:
                for part in wildcards_match.group(1).split(","):
                    if "=" in part:
                        key, val = part.split("=", 1)
                        wildcards[key.strip()] = val.strip()

            if wildcards.get("dataset") not in (None, dataset):
                continue

            msg = f"{rule} in {log_path.name}"
            gsm = wildcards.get("gsm")
            if gsm:
                sample_hints[gsm].append(msg)
            else:
                dataset_hints.append(msg)

    return sample_hints, dataset_hints


def main() -> None:
    parser = argparse.ArgumentParser(description="Report completion/failure status from prior pipeline runs")
    parser.add_argument("--dataset", help="Dataset key from config.yaml (for example: GSE264508)")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument("--base-dir", default=".", help="Repository root or pipeline working directory")
    parser.add_argument(
        "--snakemake-log-dir",
        default=".snakemake/log",
        help="Path to Snakemake log directory (relative to base-dir unless absolute)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output TSV path. Defaults to datasets/<dataset>/status/<dataset>_run_status.tsv",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = base_dir / config_path

    config_data = load_config(config_path)
    dataset = select_dataset(config_data, args.dataset)
    samples = normalize_samples(config_data["datasets"][dataset])

    snakemake_log_dir = Path(args.snakemake_log_dir)
    if not snakemake_log_dir.is_absolute():
        snakemake_log_dir = base_dir / snakemake_log_dir
    sample_hints, dataset_hints = parse_failure_hints(snakemake_log_dir, dataset)

    seurat_rds = base_dir / "datasets" / dataset / "seurat" / f"{dataset}_annotated.rds"
    seurat_pdf = base_dir / "datasets" / dataset / "seurat" / f"{dataset}_plots.pdf"
    seurat_exists = int(seurat_rds.exists()) + int(seurat_pdf.exists())
    seurat_status = stage_status(seurat_exists, 2)

    if args.output:
        out_path = Path(args.output)
        if not out_path.is_absolute():
            out_path = base_dir / out_path
    else:
        out_path = base_dir / "datasets" / dataset / "status" / f"{dataset}_run_status.tsv"

    out_path.parent.mkdir(parents=True, exist_ok=True)

    headers = [
        "dataset",
        "sample",
        "num_sra_accessions",
        "download_fastq_status",
        "missing_fastq_files",
        "cellranger_status",
        "process_seurat_status",
        "failure_hints",
    ]

    rows = []
    for gsm in sorted(samples.keys()):
        srrs = samples[gsm]
        lanes = [f"{i:03d}" for i in range(1, len(srrs) + 1)]

        expected_fastqs = []
        for lane in lanes:
            expected_fastqs.append(base_dir / "datasets" / dataset / "fastq" / f"{gsm}_S1_L{lane}_R1_001.fastq.gz")
            expected_fastqs.append(base_dir / "datasets" / dataset / "fastq" / f"{gsm}_S1_L{lane}_R2_001.fastq.gz")

        existing_fastq_count = sum(path.exists() for path in expected_fastqs)
        download_status = stage_status(existing_fastq_count, len(expected_fastqs))
        missing_fastqs = [str(p.relative_to(base_dir)).replace("\\", "/") for p in expected_fastqs if not p.exists()]

        cellranger_matrix = (
            base_dir / "datasets" / dataset / "cellranger_output" / gsm / "outs" / "filtered_feature_bc_matrix"
        )
        cellranger_status = "ok" if cellranger_matrix.is_dir() else "missing"

        hints = sample_hints.get(gsm, [])
        rows.append(
            [
                dataset,
                gsm,
                str(len(srrs)),
                download_status,
                "; ".join(missing_fastqs) if missing_fastqs else "",
                cellranger_status,
                seurat_status,
                " | ".join(sorted(set(hints))),
            ]
        )

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)

    print(f"Wrote status report: {out_path}")
    if dataset_hints:
        print("Dataset-level failure hints from Snakemake logs:")
        for hint in sorted(set(dataset_hints)):
            print(f"- {hint}")


if __name__ == "__main__":
    main()
