#!/usr/bin/env python3
"""Prepare public annotation inputs used by the APABenchmark workflow."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path


ANNOTATION_INPUTS = {
    "gencode.vM25.polyAs.gtf": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.polyAs.gtf.gz",
        "md5": "5062c08e11a7b8611d441d89cf72a65b",
        "kind": "gunzip",
    },
    "gencode.vM25.annotation.bed": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
        "md5": "0c38fc4ccbc731a2708fc91e7f1c2efd",
        "kind": "gtf_to_bed",
    },
    "atlas.clusters.3.0.GRCm38.GENCODE_M25.bed.gz": {
        "url": "https://polyasite.unibas.ch/download/atlas/3.0/GRCm38.GENCODE_M25/atlas.clusters.3.0.GRCm38.GENCODE_M25.bed.gz",
        "sha256": "f3ec3f26cb4f6d3389cbe1c72c18565a8570e104f699bb8f708a33a1932f86d5",
        "kind": "copy",
    },
    "MousePas/mm10.PAS.main.tsv": {
        "url": "https://exon.apps.wistar.org/polya_db/v4/download/4.1/MousePas.zip",
        "sha256": "fda7c18dd2ad67f816f1508d85beeee11bca02af349e46badf096c34bcb5cc9f",
        "member": "mm10.PAS.main.tsv",
        "kind": "zip_member",
    },
    "gencode.v40.polyAs.gtf": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.polyAs.gtf.gz",
        "md5": "54ff2f0a32a640ad1f0b1e8a89dc423b",
        "kind": "gunzip",
    },
    "gencode.v40.annotation.bed": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz",
        "md5": "14a867b82917c8c3006838c3a5053a3e",
        "kind": "gtf_to_bed",
    },
    "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz": {
        "url": "https://polyasite.unibas.ch/download/atlas/3.0/GRCh38.GENCODE_42/atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz",
        "sha256": "212173f4bf8fddf83b93fe37f9265aea27741faf03ed47a0435a099af06d6f6c",
        "kind": "copy",
    },
    "HumanPas/hg38.PAS.main.tsv": {
        "url": "https://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip",
        "sha256": "d00e774642ad2acdad857843ef1d5dc48e77b3f077389b6923d09450b941e9c2",
        "member": "hg38.PAS.main.tsv",
        "kind": "zip_member",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download and prepare public annotation inputs under data/raw_data/annotations."
    )
    parser.add_argument(
        "--data-root",
        default=None,
        help="Repository/project root containing data/. Defaults to the repository root.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Download cache. Defaults to data/raw_data/annotations/.cache.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing prepared files.")
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, cache_dir: Path, force: bool) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    filename = url.rsplit("/", 1)[-1]
    target = cache_dir / filename
    if target.exists() and not force:
        return target
    print(f"[download] {url}", file=sys.stderr)
    tmp = target.with_suffix(target.suffix + ".tmp")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(target)
    return target


def verify(source: Path, spec: dict[str, str]) -> None:
    if "md5" in spec:
        observed = digest(source, "md5")
        if observed != spec["md5"]:
            raise RuntimeError(f"MD5 mismatch for {source}: expected {spec['md5']}, observed {observed}")
    if "sha256" in spec:
        observed = digest(source, "sha256")
        if observed != spec["sha256"]:
            raise RuntimeError(
                f"SHA256 mismatch for {source}: expected {spec['sha256']}, observed {observed}"
            )


def atomic_path(target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=target.name + ".",
        suffix=".tmp",
        dir=target.parent,
        delete=False,
    )
    handle.close()
    return Path(handle.name)


def gunzip_to(source: Path, target: Path) -> None:
    tmp = atomic_path(target)
    with gzip.open(source, "rb") as src, tmp.open("wb") as out:
        shutil.copyfileobj(src, out)
    tmp.replace(target)


def copy_to(source: Path, target: Path) -> None:
    tmp = atomic_path(target)
    shutil.copyfile(source, tmp)
    tmp.replace(target)


def extract_zip_member(source: Path, member_name: str, target: Path) -> None:
    tmp = atomic_path(target)
    with zipfile.ZipFile(source) as archive:
        matched = [name for name in archive.namelist() if name.endswith(member_name)]
        if not matched:
            raise RuntimeError(f"Cannot find {member_name} in {source}")
        with archive.open(matched[0]) as src, tmp.open("wb") as out:
            shutil.copyfileobj(src, out)
    tmp.replace(target)


def gtf_to_bed(source_gz: Path, target: Path) -> None:
    tmp = atomic_path(target)
    with gzip.open(source_gz, "rt") as src, tmp.open("w") as out:
        for line in src:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            bed_start = max(int(start) - 1, 0)
            out.write(
                "\t".join(
                    [
                        chrom,
                        str(bed_start),
                        end,
                        ".",
                        score,
                        strand,
                        source,
                        feature,
                        frame,
                        attrs,
                    ]
                )
                + "\n"
            )
    tmp.replace(target)


def main() -> int:
    args = parse_args()
    root = Path(args.data_root).resolve() if args.data_root else repo_root()
    annotation_dir = root / "data/raw_data/annotations"
    cache_dir = Path(args.cache_dir).resolve() if args.cache_dir else annotation_dir / ".cache"
    annotation_dir.mkdir(parents=True, exist_ok=True)

    for rel_target, spec in ANNOTATION_INPUTS.items():
        target = annotation_dir / rel_target
        if target.exists() and not args.force:
            print(f"[exists] {target}", file=sys.stderr)
            continue
        source = download(spec["url"], cache_dir, args.force)
        verify(source, spec)
        print(f"[prepare] {target}", file=sys.stderr)
        if spec["kind"] == "gunzip":
            gunzip_to(source, target)
        elif spec["kind"] == "copy":
            copy_to(source, target)
        elif spec["kind"] == "zip_member":
            extract_zip_member(source, spec["member"], target)
        elif spec["kind"] == "gtf_to_bed":
            gtf_to_bed(source, target)
        else:
            raise RuntimeError(f"Unsupported input kind: {spec['kind']}")

    print("[prepare-public-inputs] done", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
