#!/bin/bash
#SBATCH --job-name=CoNbS_1Q3Q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=slurm_pipeline_%j.out
#SBATCH --error=slurm_pipeline_%j.err

echo "==================================="
echo "Starting 1Q / 3Q pipeline"
echo "==================================="
# -----------------------------------------------------------------------------
# [병렬 처리 프로세스 설명 (Parallel computing breakdown)]
# 1. julia -t auto: 
#    Julia가 SLURM 배치를 통해 할당된 16개의 CPU 쓰레드 자원(--cpus-per-task=16)을 
#    자동으로 인식하여 멀티쓰레딩 환경(-t auto)을 구성합니다.
#
# 2. calc_1Q / 3Q 계산 원리: 
#    calc 스크립트는 내부에서 src/CoreSimulation.jl 의 run_parallel_simulation! 
#    함수를 호출하며, 이곳에서 Threads.@threads 블록을 실행합니다.
#
# 3. 병렬 자원 활용법: 
#    설정된 샘플 개수(n_samples=5) 만큼 독립적인 열역학적 궤적(trajectory)이 
#    각각 독립된 쓰레드에서 병렬 계산됩니다. 즉 5개의 쓰레드가 각각 1개의 샘플을 
#    동시에(parallel) 전담하여 thermalization 및 correlation 측정을 수행하므로
#    계산 소모 시간을 최대 1/5 가량으로 대폭 단축시킵니다. 
# -----------------------------------------------------------------------------

echo ">> Running 1Q domain generation"
# 1Q: 5개의 SampledCorrelation 샘플을 병렬 쓰레드로 분석 및 결과 저장
julia -t auto calc_1Q_correlations.jl

echo ">> Running 3Q domain generation"
# 3Q: 5개의 SampledCorrelation 샘플을 병렬 쓰레드로 분석 및 결과 저장
julia -t auto calc_3Q_correlations.jl

echo ">> Running visualization"
julia plot_correlations.jl

echo "==================================="
echo "Pipeline finished"
echo "==================================="
