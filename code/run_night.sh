#!/bin/bash
# Stages run in sequence by construction -- no polling, no pgrep.
cd /Users/benoitjeanson/vsCode/TUD/IJEPES2026-OTSD/code
echo "=== highs start $(date) ===" >> logs/stages.log
julia --project=. experiments/run_ablation.jl highs > logs/highs_campaign.log 2>&1
echo "=== highs done $(date); clean gurobi re-run start ===" >> logs/stages.log
rm -rf results/*_gurobi_*
julia --project=. experiments/run_ablation.jl gurobi > logs/gurobi_clean.log 2>&1
echo "=== all done $(date) ===" >> logs/stages.log
