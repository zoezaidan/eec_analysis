#!/usr/bin/env python
import os
import ROOT

# ----------------------- CONFIG -----------------------
tree_name       = "hiEvtAnalyzer/HiTree"   # or "dir/subdir/tree"
n_jobs          = 100
job_dir         = "exec_jobs"
macro_name      = "/grid_mnt/vol_home/llr/cms/zaidan/CMSSW_10_6_48/src/EEC_Zoe/create_trees_eec.cpp"

# Inputs of the macro
input_root_file = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"
dataType = 1
pT_low   = 80
pT_high  = 140
n        = 1
btag        = True
aggregated  = True
matching    = True

# ------------------------------------------------------

# Create output directory
if not os.path.exists(job_dir):
    os.makedirs(job_dir)

# Open ROOT file
f = ROOT.TFile.Open(input_root_file)
if not f or f.IsZombie():
    raise RuntimeError("Could not open ROOT file: %s" % input_root_file)

# Support "dir/subdir/tree" paths
if "/" in tree_name:
    parts = tree_name.split("/")
    current = f
    for p in parts[:-1]:
        current = current.Get(p)
        if not current:
            raise RuntimeError("Directory '%s' not found in file" % p)
    tree = current.Get(parts[-1])
else:
    tree = f.Get(tree_name)

if not tree:
    raise RuntimeError("TTree '%s' not found in file" % tree_name)

total_events = tree.GetEntries()
print("Total events in TTree '%s': %d" % (tree_name, total_events))

# Compute chunk size
events_per_job = (total_events + n_jobs - 1) // n_jobs
print("Events per job: %d" % events_per_job)

# ----------------------- JOB CREATION -----------------------

for job_id in range(n_jobs):
    start = job_id * events_per_job
    end   = min(start + events_per_job, total_events)

    if start >= total_events:
        break

    jobfile = "%s/job_%d.sh" % (job_dir, job_id)

    with open(jobfile, "w") as fsh:
        fsh.write("#!/bin/bash\n")
        fsh.write('echo "Running on $(hostname)"\n\n')

        # Setup environment
        fsh.write("cd /home/llr/cms/zaidan/CMSSW_10_6_48/src/EEC_Zoe/\n")
        fsh.write("eval `cmsenv`\n\n")

        # ROOT macro call
        # Note: escaping quotes for shell
        macro_call = (
            'root -b -q "%s(%d, %d, %d, %d, %d, %d, %s, %s, %s, \\"job_%d.root\\")"\n'
            % (
                macro_name,
                dataType,
                pT_low,
                pT_high,
                n,
                start,
                end,
                "true" if btag else "false",
                "true" if aggregated else "false",
                "true" if matching else "false",
                job_id,
            )
        )

        fsh.write(macro_call)

    os.chmod(jobfile, 0o755)

print("Generated job scripts in: %s/" % job_dir)
