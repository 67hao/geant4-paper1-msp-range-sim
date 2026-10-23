"""
run_worker.py - Specialized worker for MSP & Range Simulation
"""
import argparse
import os
import sys
import time
import json
import random
import shutil
import subprocess
from pathlib import Path
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)

import gdrive_helper
from googleapiclient.http import MediaFileUpload

def parse_args():
    parser = argparse.ArgumentParser(description="MSP & Range Swarm Worker")
    parser.add_argument("--runner-id", type=int, required=True)
    parser.add_argument("--num-runners", type=int, default=10)
    parser.add_argument("--primaries", type=int, default=3000)
    parser.add_argument("--gdrive-folder", type=str, required=True)
    parser.add_argument("--plan-file", type=str, default="plan_swarm.json")
    parser.add_argument("--manifest-file", type=str, default="manifest.csv")
    parser.add_argument("--github-user", type=str, default=None)
    parser.add_argument("--sim-bin", type=str, default="./build/sim_msp_range")
    return parser.parse_args()

def upload_file_to_drive(service, file_path: Path, parent_folder_id: str):
    q = f"'{parent_folder_id}' in parents and name = '{file_path.name}' and trashed = false"
    for attempt in range(5):
        try:
            res = service.files().list(q=q, fields="files(id)").execute().get("files", [])
            media = MediaFileUpload(str(file_path), resumable=True)
            if res:
                f_id = res[0]["id"]
                service.files().update(fileId=f_id, media_body=media).execute()
                return f_id
            body = {"name": file_path.name, "parents": [parent_folder_id]}
            new_f = service.files().create(body=body, media_body=media, fields="id").execute()
            return new_f["id"]
        except Exception:
            time.sleep(1.0 + random.uniform(0.5, 2.0))
    raise RuntimeError(f"Failed to upload {file_path.name}")

def main():
    args = parse_args()
    acc_name = args.github_user or os.environ.get("GITHUB_REPOSITORY_OWNER") or "runner"
    print("=" * 80)
    print(f"STARTING MSP & RANGE RUNNER #{args.runner_id} (Account: {acc_name})")
    print("=" * 80)

    sim_bin = Path(args.sim_bin)
    if not sim_bin.exists():
        raise FileNotFoundError(f"Executable not found: {sim_bin}")

    service = gdrive_helper.get_drive_service()
    with open("drive_folders.json", "r", encoding="utf-8") as f:
        cf = json.load(f)
    flux_sub_map = cf.get("flux_subfolders", {})
    outputs_sub_map = cf.get("outputs_subfolders", {})

    with open(args.plan_file, "r", encoding="utf-8") as f:
        full_plan = json.load(f)
    my_jobs = full_plan.get(acc_name, {}).get(str(args.runner_id), [])

    print(f"Assigned {len(my_jobs)} jobs to runner #{args.runner_id}:")
    for j in my_jobs:
        print(f"  * {j['job_name']}")

    if not my_jobs:
        print("No jobs for this runner. Exiting cleanly.")
        return

    for idx, job in enumerate(my_jobs, 1):
        job_name = job["job_name"]
        particle = job["particle"]
        mat = job["material"]
        sample = job["sample"]
        csv_msp = job["csv_msp"]
        csv_range = job["csv_range"]
        dest_flux_folder = flux_sub_map[sample]
        dest_outputs_folder = outputs_sub_map[sample]

        # Check if already on Drive
        is_on_drive = False
        try:
            q_chk = f"'{dest_flux_folder}' in parents and name = '{csv_msp}' and trashed = false"
            res = service.files().list(q=q_chk, fields="files(id)").execute().get("files", [])
            if res:
                is_on_drive = True
        except Exception:
            pass

        if is_on_drive:
            print(f"[{idx}/{len(my_jobs)}] SKIP: {job_name} already on Drive.")
            continue

        print(f"[{idx}/{len(my_jobs)}] Running: {sim_bin} {particle} {mat} energies_kev.txt sim_{particle}_{mat} 3000 200")
        t0 = time.time()
        res_run = subprocess.run([
            str(sim_bin.resolve()),
            particle,
            mat,
            "energies_kev.txt",
            f"sim_{particle}_{mat}",
            "3000",
            "200"
        ], capture_output=True, text=True)
        t_run = time.time() - t0

        if res_run.returncode != 0:
            print(f"  [ERROR] Execution failed for {job_name} (code {res_run.returncode})")
            print(res_run.stderr[-1000:])
            continue

        print(f"  * Completed in {t_run:.2f} s")

        p_msp = Path(csv_msp)
        p_range = Path(csv_range)
        if p_msp.exists():
            upload_file_to_drive(service, p_msp, dest_flux_folder)
            print(f"  * Uploaded -> 02_Flux_Data/{sample}/{csv_msp}")
        if p_range.exists():
            upload_file_to_drive(service, p_range, dest_flux_folder)
            print(f"  * Uploaded -> 02_Flux_Data/{sample}/{csv_range}")

        for prof in Path(".").glob(f"sim_{particle}_{mat}_profile_*.csv"):
            upload_file_to_drive(service, prof, dest_outputs_folder)
            prof.unlink()

        if p_msp.exists(): p_msp.unlink()
        if p_range.exists(): p_range.unlink()

    print(f"RUNNER #{args.runner_id} COMPLETED ALL ASSIGNED JOBS!")

if __name__ == "__main__":
    main()
