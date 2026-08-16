#!/usr/bin/env bash
# wsl-bootstrap.sh — LocalColabFold (AF2/AF2-Multimer) + LightDock inside WSL2
# Invoked by D:\bioai\bin\wsl-setup.ps1 (or manually):
#   wsl -d Ubuntu -- bash /mnt/d/bioai/bin/wsl-bootstrap.sh
set -euo pipefail
export HOME=/root
cd /root
unset PYTHONPATH

echo "== 1) GPU passthrough check =="
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total --format=csv
else
    echo "WARN: nvidia-smi not visible inside WSL (driver passthrough missing?)"
fi

echo "== 2) apt deps (TUNA mirror) =="
# mirrored networking has no IPv6 route; disable v6 so Python/conda use IPv4
sudo sysctl -w net.ipv6.conf.all.disable_ipv6=1 >/dev/null 2>&1 || true
sudo cp /etc/apt/sources.list /etc/apt/sources.list.bak 2>/dev/null || true
sudo tee /etc/apt/sources.list >/dev/null <<'EOF'
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-updates main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-backports main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-security main restricted universe multiverse
EOF
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -qq
sudo apt-get install -y -qq wget curl bzip2 git >/dev/null 2>&1 || true

echo "== 3) Miniforge (local file from Windows side, or TUNA mirror) =="
if [ ! -d "$HOME/miniforge3" ]; then
    if [ -f /mnt/d/bioai/wsl/Miniforge3-Linux-x86_64.sh ]; then
        bash /mnt/d/bioai/wsl/Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
    else
        wget -q https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
        bash /tmp/miniconda.sh -b -p "$HOME/miniforge3"
    fi
fi
export PATH="$HOME/miniforge3/bin:$PATH"

cat > "$HOME/.condarc" <<'EOF'
channels:
  - conda-forge
  - bioconda
  - defaults
default_channels:
  - https://mirrors.ustc.edu.cn/anaconda/pkgs/main
  - https://mirrors.ustc.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.ustc.edu.cn/anaconda/cloud
  bioconda: https://mirrors.ustc.edu.cn/anaconda/cloud
EOF

echo "== 4) colabfold env (AF2 / AF2-Multimer) =="
source "$HOME/miniforge3/etc/profile.d/conda.sh"
if ! conda env list | grep -q '^colabfold '; then
    conda create -y -n colabfold -c conda-forge -c bioconda python=3.10 colabfold
fi
conda activate colabfold
colabfold_batch --help | head -5

echo "== 5) LightDock (pip) =="
pip install -q lightdock || true
python -c "import lightdock; print('lightdock import ok')" || true

echo "== 6) shared param dir (Windows side: D:\\bioai\\models) =="
mkdir -p /mnt/d/bioai/models/colabfold /mnt/d/bioai/models/cache
grep -q COLABFOLDDIR "$HOME/.bashrc" || echo 'export COLABFOLDDIR=/mnt/d/bioai/models/colabfold' >> "$HOME/.bashrc"
grep -q XDG_CACHE_HOME "$HOME/.bashrc" || echo 'export XDG_CACHE_HOME=/mnt/d/bioai/models/cache' >> "$HOME/.bashrc"

echo "== 7) params present? =="
ls -lh /mnt/d/bioai/models/colabfold/ 2>/dev/null || echo "params missing -> run D:\\bioai\\bin\\download-af2-params.ps1 on Windows side"

echo "BOOTSTRAP DONE"
