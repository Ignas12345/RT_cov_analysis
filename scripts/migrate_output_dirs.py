"""
Migrate output files from the old directory layout to the new modular layout.

Old layout                                          New layout
--------------------------------------------------  --------------------------------------------------
{out_dir}/{sample}/reads_with_tso/                  {out_dir}/detect_tso/{sample}/
  {sample}_cdna_001_trimmed.fastq                     {sample}_cdna_001_trimmed.fastq
  {sample}_bc_001.fastq                               {sample}_bc_001.fastq

{out_dir}/{sample}/star/                            {out_dir}/star_mapping/{sample}/
  {sample}_Aligned.sortedByCoord.out.bam              {sample}_Aligned.sortedByCoord.out.bam
  (all other STAR output files)                       (all other STAR output files)

{out_dir}/genome_and_annotations/                   {out_dir}/filter_single_isoform_genes/
  single_isoform_genes.txt                            single_isoform_genes.txt

After moving all files, the now-empty old directories are removed.

Usage
-----
    # Dry run (default) — only prints what would happen:
    python scripts/migrate_output_dirs.py

    # Actually perform the moves:
    python scripts/migrate_output_dirs.py --execute
"""

import argparse
import os
import shutil
import sys
import yaml


def load_config(path="config.yaml"):
    with open(path) as fh:
        return yaml.safe_load(fh)


def move_file(src, dst, dry_run):
    if not os.path.exists(src):
        print(f"  [SKIP]   {src}  (not found)")
        return
    if os.path.exists(dst):
        print(f"  [SKIP]   {src}  (destination already exists: {dst})")
        return
    print(f"  [MOVE]   {src}")
    print(f"        -> {dst}")
    if not dry_run:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.move(src, dst)


def move_dir_contents(src_dir, dst_dir, dry_run):
    """Move every file/subdirectory inside src_dir into dst_dir."""
    if not os.path.isdir(src_dir):
        print(f"  [SKIP]   {src_dir}/  (directory not found)")
        return
    entries = os.listdir(src_dir)
    if not entries:
        print(f"  [SKIP]   {src_dir}/  (already empty)")
        return
    for entry in entries:
        src = os.path.join(src_dir, entry)
        dst = os.path.join(dst_dir, entry)
        if os.path.exists(dst):
            print(f"  [SKIP]   {src}  (destination already exists: {dst})")
            continue
        print(f"  [MOVE]   {src}")
        print(f"        -> {dst}")
        if not dry_run:
            os.makedirs(dst_dir, exist_ok=True)
            shutil.move(src, dst)


def remove_if_empty(path, dry_run):
    """Remove a directory (and empty parents up to out_dir) if it is empty."""
    if not os.path.isdir(path):
        return
    if os.listdir(path):
        print(f"  [KEEP]   {path}/  (not empty, skipping removal)")
        return
    print(f"  [RMDIR]  {path}/")
    if not dry_run:
        os.rmdir(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--execute", action="store_true",
                        help="Actually perform the moves (default: dry run)")
    parser.add_argument("--config", default="config.yaml",
                        help="Path to config.yaml (default: config.yaml)")
    args = parser.parse_args()

    dry_run = not args.execute
    if dry_run:
        print("=== DRY RUN — pass --execute to apply changes ===\n")

    cfg = load_config(args.config)
    out_dir = cfg["out_dir"]
    samples = cfg["samples"]

    # -----------------------------------------------------------------------
    # 1. detect_tso: {out_dir}/{sample}/reads_with_tso/ -> detect_tso/{sample}/
    # -----------------------------------------------------------------------
    print("--- detect_tso outputs ---")
    for sample in samples:
        src_dir = os.path.join(out_dir, sample, "reads_with_tso")
        dst_dir = os.path.join(out_dir, "detect_tso", sample)
        for fname in [
            f"{sample}_cdna_001_trimmed.fastq",
            f"{sample}_bc_001.fastq",
        ]:
            move_file(os.path.join(src_dir, fname),
                      os.path.join(dst_dir, fname),
                      dry_run)

    # -----------------------------------------------------------------------
    # 2. star_mapping: {out_dir}/{sample}/star/ -> star_mapping/{sample}/
    # -----------------------------------------------------------------------
    print("\n--- star_mapping outputs ---")
    for sample in samples:
        src_dir = os.path.join(out_dir, sample, "star")
        dst_dir = os.path.join(out_dir, "star_mapping", sample)
        move_dir_contents(src_dir, dst_dir, dry_run)

    # -----------------------------------------------------------------------
    # 3. filter_single_isoform_genes: genome_and_annotations/ -> filter_single_isoform_genes/
    # -----------------------------------------------------------------------
    print("\n--- filter_single_isoform_genes outputs ---")
    move_file(
        os.path.join(out_dir, "genome_and_annotations", "single_isoform_genes.txt"),
        os.path.join(out_dir, "filter_single_isoform_genes", "single_isoform_genes.txt"),
        dry_run,
    )

    # -----------------------------------------------------------------------
    # 4. Remove now-empty old directories
    # -----------------------------------------------------------------------
    print("\n--- removing old directories (if empty) ---")
    for sample in samples:
        remove_if_empty(os.path.join(out_dir, sample, "reads_with_tso"), dry_run)
        remove_if_empty(os.path.join(out_dir, sample, "star"), dry_run)
        remove_if_empty(os.path.join(out_dir, sample), dry_run)
    remove_if_empty(os.path.join(out_dir, "genome_and_annotations"), dry_run)

    print()
    if dry_run:
        print("=== Dry run complete. Run with --execute to apply the above changes. ===")
    else:
        print("=== Migration complete. ===")


if __name__ == "__main__":
    main()
