#!/usr/bin/env bash
# wsl-download-params.sh v3 — AF2 params via GCS, reboot-safe
# Chunks live on /mnt/d (persist across reboots); resume with curl -C -.
# Layout matches colabfold: COLABFOLDDIR/params/<tars>
set -u
export HOME=/root
cd /root
unset PYTHONPATH

DEST=/mnt/d/bioai/models/colabfold/params
CKDIR=/mnt/d/bioai/models/colabfold/.chunks
mkdir -p "$DEST" "$CKDIR"
CHUNKS=8

dl() {
    local url="$1" out="$2" base="$3"
    local len
    len=$(curl -4 -sI --max-time 30 "$url" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1)
    if [ -z "$len" ] || [ "$len" -le 0 ]; then echo "no content-length for $url"; return 1; fi
    if [ -f "$out" ]; then
        local have
        have=$(stat -c%s "$out" 2>/dev/null || echo 0)
        if [ "$have" = "$len" ]; then echo "already complete: $out"; return 0; fi
        rm -f "$out"
    fi
    echo "== $base total=$((len/1048576))MB chunks=$CHUNKS start $(date +%H:%M:%S)"
    local size=$(( (len + CHUNKS - 1) / CHUNKS ))
    local i start end want
    for ((i=0;i<CHUNKS;i++)); do
        start=$((i*size)); end=$((start+size-1))
        [ $end -ge $len ] && end=$((len-1))
        [ $start -ge $len ] && continue
        want=$((end-start+1))
        (
            local got=0 a
            for ((a=1;a<=60;a++)); do
                got=$(stat -c%s "$CKDIR/$base.chunk$i" 2>/dev/null || echo 0)
                if [ "$got" = "$want" ]; then break; fi
                if [ "$got" -gt "$want" ]; then rm -f "$CKDIR/$base.chunk$i"; got=0; fi
                curl -4 -s -r "$((start+got))-$end" -o "$CKDIR/$base.chunk$i" --max-time 1800 "$url"
                got=$(stat -c%s "$CKDIR/$base.chunk$i" 2>/dev/null || echo 0)
                [ "$got" = "$want" ] && break
                sleep 3
            done
            if [ "$got" != "$want" ]; then echo "CHUNK $base.$i incomplete want=$want got=$got" >&2; fi
        ) &
    done
    wait
    : > "$out.tmp"
    for ((i=0;i<CHUNKS;i++)); do
        [ -f "$CKDIR/$base.chunk$i" ] && cat "$CKDIR/$base.chunk$i" >> "$out.tmp"
    done
    local got
    got=$(stat -c%s "$out.tmp")
    if [ "$got" = "$len" ]; then
        mv "$out.tmp" "$out"
        echo "VERIFIED OK: $out ($((got/1048576))MB)"
        rm -f "$CKDIR"/$base.chunk*
    else
        echo "ASSEMBLY MISMATCH got=$got want=$len (chunks kept for resume)"
    fi
}

rm -f "$DEST/alphafold_params_colab_2022-12-06.tar" "$DEST/alphafold_params_colab_2022-03-02.tar" "$DEST/alphafold_params_2021-07-14.tar"
rm -f "$CKDIR"/monomer.chunk* "$CKDIR"/multimer.chunk* 2>/dev/null
# colabfold v1.5.5 mapping:
#   2022-12-06.tar = alphafold2_multimer_v3 (marker download_complexes_multimer_v3_finished.txt)
#   2022-03-02.tar = alphafold2_multimer_v2 (marker download_complexes_multimer_v2_finished.txt)
#   2021-07-14.tar = AlphaFold2-ptm monomer   (marker download_finished.txt)
# Wave 1: multimer_v3 alone (needed for acceptance) — full 8-connection speed
dl "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar" "$DEST/alphafold_params_colab_2022-12-06.tar" "multimer_v3"
# Wave 2: v2 + monomer
dl "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-03-02.tar" "$DEST/alphafold_params_colab_2022-03-02.tar" "multimer_v2" &
dl "https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar" "$DEST/alphafold_params_2021-07-14.tar" "monomer_ptm" &
wait
echo "ALL PARAMS DONE"
