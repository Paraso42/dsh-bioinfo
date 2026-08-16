#!/usr/bin/env python3
# LICENSE NOTE: this script wraps prodigy-prot / prodigy-lig (Bonvin lab,
# Utrecht University). PRODIGY is free for ACADEMIC use only; a commercial
# license must be obtained from Utrecht University. Cite: Xue LC, et al.
# "PRODIGY: a web server for predicting the binding affinity of protein-protein
# complexes." Bioinformatics. 2016;32(23):3676-3678.
r"""prodigy_affinity.py — PRODIGY 结合亲和力预测包装(ΔG / Kd)

调用 Bonvin 实验室 PRODIGY(本机 venv 内 prodigy-prot / prodigy-lig),解析输出为 JSON。

蛋白-蛋白(prodigy-prot):
  & 'D:\bioai\venv\Scripts\python.exe' prodigy_affinity.py --complex complex.pdb --chains A B \
      --temperature 25 --out affinity.json
蛋白-小分子(prodigy-lig,配合 vina_dock/virtual_screen.py 使用):
  & 'D:\bioai\venv\Scripts\python.exe' prodigy_affinity.py --complex rec.pdb --ligand lig.sdf \
      --out affinity.json

输出示例:
  {"dg_kcal_mol": -15.2, "kd_M": 7.3e-12, "temperature_C": 25.0, ...}

解析基于 PRODIGY 2.x 标准输出("Predicted binding affinity (kcal.mol-1)" /
"Predicted dissociation constant (M)");若上游改格式则原样保存 raw 供排查。

编程调用:
  from prodigy_affinity import predict_pp, predict_lig
"""
import argparse
import json
import os
import re
import subprocess
import sys

PRODIGY_BIN = None  # 自动探测


def _find_prodigy():
    global PRODIGY_BIN
    if PRODIGY_BIN:
        return PRODIGY_BIN
    candidates = [
        os.path.join(sys.prefix, "Scripts", "prodigy.exe"),
        os.path.join(sys.prefix, "Scripts", "prodigy_lig.exe"),
    ]
    for c in candidates:
        if os.path.exists(c):
            PRODIGY_BIN = c
            return c
    # 已安装但 console script 名不同:问 pip
    try:
        r = subprocess.run([sys.executable, "-m", "pip", "show", "prodigy-prot"],
                           capture_output=True, text=True, timeout=60)
        if r.returncode == 0:
            raise FileNotFoundError(
                "prodigy-prot is installed but no prodigy.exe found in %s\\Scripts — "
                "inspect `pip show -f prodigy-prot`" % sys.prefix)
    except FileNotFoundError:
        raise
    raise FileNotFoundError(
        "PRODIGY not found: install with  pip install prodigy-prot  (protein-protein) "
        "and  pip install prodigy-lig  (protein-ligand) into D:\\bioai\\venv")


def _run(cmd, timeout=600):
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"   # PRODIGY 输出含 ° 等字符,GBK 控制台会崩溃
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env,
                       encoding="utf-8", errors="replace")
    if r.returncode != 0:
        raise RuntimeError("prodigy failed (rc=%s): %s" % (r.returncode, (r.stderr or r.stdout)[-800:]))
    return (r.stdout or "") + "\n" + (r.stderr or "")


def _parse(out):
    """从 PRODIGY 文本输出提取 ΔG / Kd / 温度。"""
    dg = kd = temp = None
    m = re.search(r"binding affinity \(kcal\.mol-1\):\s*([-+]?\d+(?:\.\d+)?)", out)
    if m:
        dg = float(m.group(1))
    m = re.search(r"dissociation constant \(M\)[^:]*:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", out)
    if m:
        kd = float(m.group(1))
    m = re.search(r"at\s+([-+]?\d+(?:\.\d+)?)\s*[°˚]?C", out)
    if m:
        temp = float(m.group(1))
    return dg, kd, temp


def predict_pp(complex_pdb, chains=("A", "B"), temperature=25.0, timeout=600):
    """蛋白-蛋白亲和力。返回 dict {dg_kcal_mol, kd_M, temperature_C, raw, n_ic_contacts...}"""
    bin_path = _find_prodigy()
    cmd = [bin_path, complex_pdb, "--selection", chains[0], chains[1],
           "--temperature", str(temperature)]
    out = _run(cmd, timeout=timeout)
    dg, kd, temp = _parse(out)
    report = {
        "mode": "protein-protein",
        "complex": complex_pdb, "chains": list(chains),
        "dg_kcal_mol": dg, "kd_M": kd,
        "temperature_C": temp if temp is not None else temperature,
    }
    # IC 统计(输出常带 "IC charged-charged / charged-polar / ...")
    m = re.search(r"IC\s+(\d+)\s+(\d+)\s+(\d+)", out)
    if m:
        report["ic_contacts"] = {"charged_charged": int(m.group(1)),
                                 "charged_polar": int(m.group(2)),
                                 "charged_apolar": int(m.group(3))}
    report["raw"] = out.strip()
    return report


def predict_lig(complex_pdb, ligand, temperature=25.0, timeout=600):
    """蛋白-小分子亲和力(prodigy-lig)。ligand 为 .sdf/.mol/.mol2 路径。"""
    bin_path = _find_prodigy()
    import glob
    lig_bin = os.path.join(sys.prefix, "Scripts", "prodigy_lig.exe")
    if os.path.exists(lig_bin):
        bin_path = lig_bin
    elif not bin_path.lower().endswith("prodigy_lig"):
        # prodigy-prot 的 exe 名存在但 prodigy-lig 未装时给出明确提示
        try:
            r = subprocess.run([sys.executable, "-m", "pip", "show", "prodigy-lig"],
                               capture_output=True, text=True, timeout=60)
            if r.returncode != 0:
                raise FileNotFoundError("prodigy-lig not installed: pip install prodigy-lig")
        except FileNotFoundError:
            raise
    cmd = [bin_path, "--receptor", complex_pdb, "--ligand", ligand,
           "--temperature", str(temperature)]
    out = _run(cmd, timeout=timeout)
    dg, kd, temp = _parse(out)
    return {
        "mode": "protein-ligand", "complex": complex_pdb, "ligand": ligand,
        "dg_kcal_mol": dg, "kd_M": kd,
        "temperature_C": temp if temp is not None else temperature,
        "raw": out.strip(),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="PRODIGY binding affinity wrapper (prodigy-prot / prodigy-lig)")
    ap.add_argument("--complex", required=True, help="complex PDB (protein-protein) or receptor PDB (with --ligand)")
    ap.add_argument("--chains", nargs="+", default=["A", "B"], help="two chain ids (default A B)")
    ap.add_argument("--ligand", help="ligand file (.sdf/.mol/.mol2) — switches to prodigy-lig mode")
    ap.add_argument("--temperature", type=float, default=25.0)
    ap.add_argument("--out", help="write JSON report to file")
    args = ap.parse_args(argv)
    try:
        if args.ligand:
            report = predict_lig(args.complex, args.ligand, temperature=args.temperature)
        else:
            if len(args.chains) != 2:
                print("ERROR: --chains needs exactly two chain ids", file=sys.stderr)
                return 2
            report = predict_pp(args.complex, chains=tuple(args.chains), temperature=args.temperature)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    print("mode: %s" % report["mode"])
    if report["dg_kcal_mol"] is not None:
        print("dG_bind = %.2f kcal/mol" % report["dg_kcal_mol"])
    if report["kd_M"] is not None:
        print("Kd      = %.3g M" % report["kd_M"])
    if "ic_contacts" in report:
        print("IC contacts: %s" % report["ic_contacts"])
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
        print("JSON written: %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
