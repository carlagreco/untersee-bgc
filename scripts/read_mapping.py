import argparse
import csv
import json
import os
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

RESULT_COLUMNS = [
    "sra_accession",
    "mag_id",
    "bgc_id",
    "containment",
    "bioproject",
    "bgc_contig",
    "bgc_type",
    "bgc_length",
    "breadth",
    "mean_depth",
    "status",
]


def time_now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def run_cmd(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{exc.stderr.strip()}") from exc


def ensure_tools():
    required = ["curl", "minimap2", "samtools", "fasterq-dump"]
    missing = [x for x in required if shutil.which(x) is None]
    if missing:
        raise SystemExit(f"Missing required tools: {', '.join(missing)}")


def log_failure(log_path, accession, stage, detail):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as handle:
        handle.write(f"{time_now()}\t{accession}\t{stage}\t{detail}\n")


def remove_file(path):
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def load_bgcs_by_mag(regions_path, contig_bin_table):
    contig_to_mag = {}
    with open(contig_bin_table) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                contig_to_mag[parts[0]] = parts[1]

    bgcs_by_mag = defaultdict(list)
    with open(regions_path, newline="") as handle:
        for row in csv.DictReader(handle):
            contig_id = row["contig_id"]
            mag_id = contig_to_mag.get(contig_id)
            if not mag_id:
                continue
            bgc = {
                "bgc_id": f"{contig_id}.region{int(row['name']):03d}",
                "contig_id": contig_id,
                "bgc_type": row["product"],
                "start_0based": int(row["start"]),
                "end_1based": int(row["end"]),
                "length": int(row["length"]),
            }
            bgcs_by_mag[mag_id].append(bgc)
    return dict(bgcs_by_mag)


def load_jobs(hit_table, bgcs_by_mag, min_containment, only_set):
    jobs = defaultdict(list)
    with open(hit_table, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["sra_accession"]
            if only_set and accession not in only_set:
                continue
            try:
                containment = float(row["containment"])
            except ValueError:
                continue
            if containment <= min_containment:
                continue

            mag_id = row["mag_id"]
            bgcs = bgcs_by_mag.get(mag_id)
            if not bgcs:
                continue

            jobs[accession].append(
                {
                    "sra_accession": accession,
                    "mag_id": mag_id,
                    "containment": containment,
                    "bioproject": row["bioproject"],
                    "bgcs": bgcs,
                }
            )
    return dict(jobs)


def get_ena_fastq_urls(accession):
    api = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={accession}&result=read_run&fields=fastq_ftp&format=json"
    )
    try:
        with urlopen(Request(api), timeout=30) as response:
            payload = json.loads(response.read())
        if not payload:
            return []
        fastq_field = payload[0].get("fastq_ftp") or ""
        return [f"https://{part}" for part in fastq_field.split(";") if part]
    except (OSError, TimeoutError, HTTPError, URLError, json.JSONDecodeError) as exc:
        print(f"    ENA lookup failed for {accession}: {exc}", flush=True)
        return []


def download_reads(accession, reads_dir, threads):
    reads_dir.mkdir(parents=True, exist_ok=True)

    existing = sorted(reads_dir.glob(f"{accession}*.fastq*"))
    if existing:
        print(f"    Reusing existing reads for {accession}", flush=True)
        return existing, ""

    for attempt in [1, 2]:
        temp_dir = reads_dir / f"fasterq_tmp_{accession}_try{attempt}"
        if temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
        temp_dir.mkdir(parents=True, exist_ok=True)

        print(f"    Downloading {accession} with fasterq-dump (try {attempt}/2)...", flush=True)
        try:
            run_cmd(
                [
                    "fasterq-dump",
                    "--size-check",
                    "off",
                    "--split-files",
                    "--threads",
                    str(threads),
                    "--outdir",
                    str(reads_dir),
                    "--temp",
                    str(temp_dir),
                    accession,
                ]
            )
            download_files_exist = sorted(reads_dir.glob(f"{accession}*.fastq*"))
            if download_files_exist:
                shutil.rmtree(temp_dir, ignore_errors=True)
                return download_files_exist, ""
        except RuntimeError as exc:
            print(f"    fasterq-dump failed for {accession} on try {attempt}: {exc}", flush=True)
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    print(f"    Trying ENA fallback for {accession}...", flush=True)

    urls = get_ena_fastq_urls(accession)
    if not urls:
        print(f"    {accession}: URL not found on ENA", flush=True)
        return [], "ENA URL not found"

    downloaded = []
    for url in urls:
        target = reads_dir / url.rsplit("/", 1)[-1]
        print(f"    Downloading {target.name} from ENA...", flush=True)
        try:
            run_cmd(
                [
                    "curl",
                    "-sS",
                    "-L",
                    "--retry",
                    "3",
                    "--retry-delay",
                    "5",
                    "-o",
                    str(target),
                    url,
                ]
            )
        except RuntimeError:
            for path in downloaded:
                remove_file(path)
            remove_file(target)
            return [], "ENA download failed"

        if target.exists() and target.stat().st_size > 0:
            downloaded.append(target)
            continue

        for path in downloaded:
            remove_file(path)
        remove_file(target)
        return [], "ENA download failed"
    return sorted(downloaded), ""

def map_reads(accession, mag_id, ref_path, reads, bam_dir, threads):
    bam_dir.mkdir(parents=True, exist_ok=True)
    safe_mag = mag_id.replace("/", "_")
    sam_path = bam_dir / f"{accession}.{safe_mag}.sam"
    bam_path = bam_dir / f"{accession}.{safe_mag}.sorted.bam"

    print(f"    Mapping {accession} -> {ref_path.name} (-x sr)", flush=True)
    try:
        run_cmd(
            [
                "minimap2",
                "-a",
                "-x",
                "sr",
                "-t",
                str(threads),
                "-o",
                str(sam_path),
                str(ref_path),
                *[str(x) for x in reads],
            ]
        )
    except RuntimeError:
        remove_file(sam_path)
        return None, "minimap2 failed"
    if not sam_path.exists() or sam_path.stat().st_size == 0:
        remove_file(sam_path)
        return None, "minimap2 failed"

    try:
        run_cmd(["samtools", "sort", "-@", str(threads), "-o", str(bam_path), str(sam_path)])
    except RuntimeError:
        remove_file(sam_path)
        remove_file(bam_path)
        return None, "samtools sort failed"
    remove_file(sam_path)
    if not bam_path.exists() or bam_path.stat().st_size == 0:
        remove_file(bam_path)
        return None, "samtools sort failed"

    try:
        run_cmd(["samtools", "index", str(bam_path)])
    except RuntimeError:
        return None, "samtools index failed"

    return bam_path, ""


def calc_coverage(bam_path, bgc):
    region = f"{bgc['contig_id']}:{bgc['start_0based'] + 1}-{bgc['end_1based']}"
    cov_path = bam_path.with_suffix(".coverage.tsv")
    try:
        run_cmd(["samtools", "coverage", "-H", "-r", region, "-o", str(cov_path), str(bam_path)])
    except RuntimeError:
        remove_file(cov_path)
        return 0.0, 0.0, "absent"
    if not cov_path.exists():
        return 0.0, 0.0, "absent"

    try:
        fields = cov_path.read_text().strip().split("\t")
        breadth = float(fields[6]) / 100.0
        mean_depth = float(fields[7])
    except (OSError, ValueError, IndexError):
        remove_file(cov_path)
        return 0.0, 0.0, "absent"

    remove_file(cov_path)
    status = "present" if breadth >= 0.8 else "partial" if breadth >= 0.3 else "absent"
    return breadth, mean_depth, status


def write_sra_results(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def cleanup_reads(reads_dir, accession):
    for path in reads_dir.glob(f"{accession}*"):
        remove_file(path)


def process_accession(accession, groups, refs_dir, reads_dir, bam_dir, results_file, failure_log, threads, keep_bam):
    print(f"{accession}: {len(groups)} hit(s)", flush=True)
    reads, download_error = download_reads(accession, reads_dir, threads)
    if not reads:
        reason = download_error if download_error else "no FASTQ files"
        log_failure(failure_log, accession, "download_failed", reason)
        print(f"  Download failed for {accession}: {reason}", flush=True)
        return

    rows = []
    for hit in groups:
        mag_id = hit["mag_id"]
        ref_path = refs_dir / f"{mag_id}.bgc_contigs.fasta"
        if not ref_path.exists():
            log_failure(failure_log, accession, "missing_reference", ref_path.name)
            print(f"  Missing reference: {ref_path.name}", flush=True)
            continue

        bam_path, error = map_reads(accession, mag_id, ref_path, reads, bam_dir, threads)
        if bam_path is None:
            log_failure(failure_log, accession, f"mapping_failed:{mag_id}", error)
            print(f"  Mapping failed for {mag_id}", flush=True)
            continue

        for bgc in hit["bgcs"]:
            breadth, mean_depth, status = calc_coverage(bam_path, bgc)
            rows.append(
                {
                    "sra_accession": accession,
                    "mag_id": mag_id,
                    "bgc_id": bgc["bgc_id"],
                    "containment": f"{hit['containment']:.4f}",
                    "bioproject": hit["bioproject"],
                    "bgc_contig": bgc["contig_id"],
                    "bgc_type": bgc["bgc_type"],
                    "bgc_length": str(bgc["length"]),
                    "breadth": f"{breadth:.4f}",
                    "mean_depth": f"{mean_depth:.2f}",
                    "status": status,
                }
            )
            print(
                f"  {mag_id} | {bgc['bgc_id']}: breadth={breadth:.1%}, depth={mean_depth:.1f}x, {status}",
                flush=True,
            )

        if not keep_bam:
            remove_file(bam_path)
            remove_file(Path(str(bam_path) + ".bai"))

    if rows:
        write_sra_results(results_file, rows)
        print(f"  Wrote {len(rows)} row(s) to {results_file.name}", flush=True)
    else:
        print(f"  No rows written for {accession}", flush=True)

    cleanup_reads(reads_dir, accession)
    print(flush=True)


def parse_args():
    base = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description="Simple read mapping per SRA")
    parser.add_argument("--hit-table", type=Path, required=True)
    parser.add_argument("--regions", type=Path, default=base / "antismash" / "regions_summary.csv")
    parser.add_argument("--contig-bin-table", type=Path, default=base / "maxbin" / "contig_bin_table.tsv")
    parser.add_argument("--refs-dir", type=Path, default=base / "bgc_biogeography" / "refs")
    parser.add_argument("--outdir", type=Path, default=base / "bgc_biogeography")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-containment", type=float, default=0.7)
    parser.add_argument("--keep-bam", action="store_true")
    parser.add_argument("--only", nargs="+", metavar="SRR")
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_tools()

    outdir = args.outdir
    results_dir = outdir / "results"
    failure_log = outdir / "download_failures.log"
    reads_dir = outdir / "reads"
    bam_dir = outdir / "bams"

    bgcs_by_mag = load_bgcs_by_mag(args.regions, args.contig_bin_table)
    only_set = set(args.only) if args.only else None
    jobs = load_jobs(args.hit_table, bgcs_by_mag, args.min_containment, only_set)

    print(f"{len(jobs)} accession(s) in input", flush=True)
    for accession in sorted(jobs):
        results_file = results_dir / f"{accession}.coverage.tsv"
        if results_file.exists() and results_file.stat().st_size > 0:
            print(f"Skipping {accession}: {results_file.name} already exists", flush=True)
            continue
        process_accession(
            accession,
            jobs[accession],
            args.refs_dir,
            reads_dir,
            bam_dir,
            results_file,
            failure_log,
            args.threads,
            args.keep_bam,
        )

    print(f"Results dir: {results_dir}", flush=True)
    if failure_log.exists():
        print(f"Failures: {failure_log}", flush=True)


if __name__ == "__main__":
    main()
