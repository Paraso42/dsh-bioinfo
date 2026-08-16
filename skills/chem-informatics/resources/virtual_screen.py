#!/usr/bin/env python3
r"""virtual_screen.py — AutoDock Vina 批量虚拟筛选管道(化合物库 → 对接 → 打分排序 → 报告)

流程:受体 PDB(自动刚性 PDBQT,排除水与指定 HETATM)→ 逐配体 meeko 准备 →
Vina 对接 → 结果增量写 results.csv → 排序汇总 + 前 N 姿势转 PDB。
支持断点续跑(--resume:跳过已有结果行);可选 PRODIGY-Lig 对最优姿势预测结合亲和力。

用法:
  & 'D:\bioai\venv\Scripts\python.exe' virtual_screen.py \
      --receptor 1yph_protein.pdb --ligands library.csv --smiles-col smiles \
      --ref-ligand 1yph_ligand.pdb --exclude-res BEN \
      --exhaustiveness 8 --num-modes 3 --top 5 --outdir D:\bioai\jobs\vscreen1
  # 或显式盒子: --center 10 10 10 --size 20 20 20
  # 附加亲和力: --prodigy

编程调用:
  from virtual_screen import screen
  report = screen("rec.pdb", [("benzamidine", "NC(=N)c1ccccc1"), ...], center=..., size=...)
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys

DEFAULT_BIN_DIR = r"D:\bioai\bin"
AD_TYPE = {"C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P", "F": "F",
           "CL": "Cl", "BR": "Br", "I": "I", "H": "HD", "FE": "Fe", "ZN": "Zn",
           "MG": "Mg", "CA": "CA", "MN": "Mn", "CO": "Co", "NA": "Na", "K": "K"}


# ── 受体刚性准备(与 pdb_to_pdbqt.py 同算法,内联以保持脚本独立)──────────────
def pdb_to_pdbqt_text(pdb_text, exclude_res=("HOH", "WAT")):
    def _elem(line):
        if len(line) >= 78:
            e = line[76:78].strip().upper()
            if e:
                return e
        name = line[12:16].strip()
        if len(name) >= 2 and name[0].isdigit():
            name = name[1:]
        if name and name[0].isalpha():
            e = name[0]
            if len(name) >= 2 and name[1].isalpha() and name[1].islower():
                e += name[1].lower()
            return e.upper()
        return ""
    lines, serial, n = [], 0, 0
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        resname = line[17:20].strip()
        if resname in exclude_res:
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        try:
            occ = float(line[54:60])
        except ValueError:
            occ = 1.0
        try:
            bf = float(line[60:66])
        except ValueError:
            bf = 0.0
        atype = AD_TYPE.get(_elem(line)) or "C"
        serial += 1
        lines.append("%s%5d %-4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f      0.00 %-2s"
                     % (line[:6], serial, line[12:16], resname, line[21:22],
                        line[22:26].strip(), x, y, z, occ, bf, atype))
        n += 1
    if n == 0:
        raise ValueError("no ATOM/HETATM records parsed from receptor PDB")
    return "REMARK    rigid receptor prepared by virtual_screen.py\n" + "\n".join(lines) + "\n"


# ── 配体准备(meeko,同 vina_dock.py 配方)─────────────────────────────────────
def smiles_to_pdbqt(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("cannot parse SMILES: %s" % smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
    AllChem.MMFFOptimizeMolecule(mol)
    out = MoleculePreparation().prepare(mol)
    if isinstance(out, tuple):
        setups = out[0] if isinstance(out[0], list) else [out[0]]
    else:
        setups = out if isinstance(out, list) else [out]
    if not setups:
        raise ValueError("meeko produced no setup for %s" % smiles)
    result = PDBQTWriterLegacy.write_string(setups[0])
    pdbqt, is_ok, err = (result + (None,) * 3)[:3] if isinstance(result, tuple) else (result, True, None)
    if not is_ok or not pdbqt:
        raise ValueError("meeko write failed: %s" % err)
    return pdbqt


# ── Vina ─────────────────────────────────────────────────────────────────────
def find_vina_bin():
    hits = sorted(glob.glob(os.path.join(DEFAULT_BIN_DIR, "vina*.exe")))
    return hits[0] if hits else None


def run_vina(vina_bin, receptor_pdbqt, ligand_pdbqt, center, size, out_pdbqt,
             exhaustiveness=8, n_poses=3, seed=None):
    cmd = [vina_bin, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exhaustiveness), "--num_modes", str(n_poses),
           "--out", out_pdbqt]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError("vina failed (rc=%s): %s" % (r.returncode, (r.stdout + r.stderr)[-500:]))
    return r.stdout


def best_score_from_pdbqt(path):
    best = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                try:
                    s = float(line.split()[3])
                except (IndexError, ValueError):
                    continue
                best = s if best is None or s < best else best
    return best


def write_pose_pdb(model_lines, out_path):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("REMARK   converted from Vina PDBQT output\nMODEL        1\n")
        for line in model_lines:
            if line.startswith(("ATOM", "HETATM")):
                f.write(line[:54].rstrip() + "\n")
        f.write("ENDMDL\nEND\n")


def parse_first_model_lines(pdbqt_path):
    """取第一个(最优)MODEL 的 ATOM/HETATM 行。"""
    lines, in_model, seen = [], False, False
    with open(pdbqt_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model, seen = True, True
            elif line.startswith("ENDMDL"):
                if seen:
                    break
            elif in_model:
                lines.append(line)
    return lines


# ── 盒子推断(参考配体)───────────────────────────────────────────────────────
def box_from_ref_ligand(pdb_path, padding=8.0, min_side=20.0):
    xyz = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() not in ("HOH", "WAT"):
                try:
                    xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    continue
    if not xyz:
        raise ValueError("no coordinates in ref ligand %s" % pdb_path)
    xs = [p[0] for p in xyz]; ys = [p[1] for p in xyz]; zs = [p[2] for p in xyz]
    center = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2)
    size = tuple(max(min_side, (max(v) - min(v)) + 2 * padding)
                 for v in ((min(xs), max(xs)), (min(ys), max(ys)), (min(zs), max(zs))))
    return center, size


def _load_library(path, smiles_col=None):
    rows = []
    with open(path, encoding="utf-8-sig", newline="") as f:
        for rec in csv.DictReader(f):
            col = smiles_col or next((k for k in rec if "smiles" in k.lower()), None)
            if col is None:
                raise ValueError("--smiles-col required (no 'smiles' column found)")
            smi = (rec.get(col) or "").strip()
            if not smi:
                continue
            rid = rec.get("id") or rec.get("name") or rec.get("compound") or smi
            rows.append((str(rid), smi))
    if not rows:
        raise ValueError("no ligands parsed from %s" % path)
    return rows


# ── 主流程 ───────────────────────────────────────────────────────────────────
def screen(receptor_pdb, ligands, outdir, center=None, size=None,
           exhaustiveness=8, n_poses=3, seed=None, top=5, vina_bin=None,
           exclude_res=("HOH", "WAT"), resume=True, prodigy=False):
    vina_bin = vina_bin or find_vina_bin()
    if not vina_bin or not os.path.exists(vina_bin):
        raise FileNotFoundError("vina binary not found (D:\\bioai\\bin)")
    os.makedirs(outdir, exist_ok=True)

    # 受体准备(跳过已存在的)
    rec_pdbqt = os.path.join(outdir, "receptor.pdbqt")
    if not os.path.exists(rec_pdbqt) or os.path.getsize(rec_pdbqt) == 0:
        with open(receptor_pdb, encoding="utf-8") as f:
            text = pdb_to_pdbqt_text(f.read(), exclude_res=tuple(exclude_res))
        with open(rec_pdbqt, "w", encoding="utf-8") as f:
            f.write(text)

    results_csv = os.path.join(outdir, "results.csv")
    done = {}
    if resume and os.path.exists(results_csv):
        with open(results_csv, encoding="utf-8", newline="") as f:
            for rec in csv.DictReader(f):
                if rec.get("best_score"):
                    done[rec["id"]] = rec

    fieldnames = ["rank", "id", "smiles", "best_score", "ligand_pdbqt", "out_pdbqt",
                  "dg_kcal_mol", "kd_M", "error"]
    csv_file = open(results_csv, "a", encoding="utf-8", newline="")
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    if not resume or os.path.getsize(results_csv) == 0:
        writer.writeheader()

    rows = []
    try:
        for i, (rid, smi) in enumerate(ligands):
            if rid in done and not done[rid].get("error"):
                print("[%d/%d] %s: skip (done, score=%s)" % (i + 1, len(ligands), rid, done[rid]["best_score"]))
                rows.append(done[rid])
                continue
            safe = "".join(c if c.isalnum() else "_" for c in rid)[:40] or ("lig_%03d" % i)
            lig_pdbqt = os.path.join(outdir, "%s.pdbqt" % safe)
            out_pdbqt = os.path.join(outdir, "%s_out.pdbqt" % safe)
            rec = {"rank": None, "id": rid, "smiles": smi, "best_score": None,
                   "ligand_pdbqt": lig_pdbqt, "out_pdbqt": out_pdbqt,
                   "dg_kcal_mol": None, "kd_M": None, "error": None}
            try:
                if not os.path.exists(lig_pdbqt) or os.path.getsize(lig_pdbqt) == 0:
                    with open(lig_pdbqt, "w", encoding="utf-8") as f:
                        f.write(smiles_to_pdbqt(smi))
                run_vina(vina_bin, rec_pdbqt, lig_pdbqt, center, size, out_pdbqt,
                         exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed)
                score = best_score_from_pdbqt(out_pdbqt)
                rec["best_score"] = score if score is not None else None
                if rec["best_score"] is None:
                    rec["error"] = "no VINA RESULT in output"
                print("[%d/%d] %s: best score = %s" % (i + 1, len(ligands), rid, rec["best_score"]))
            except Exception as e:
                rec["error"] = str(e)
                print("[%d/%d] %s: ERROR %s" % (i + 1, len(ligands), rid, e))
            writer.writerow(rec)
            csv_file.flush()
            rows.append(rec)
    finally:
        csv_file.close()

    ok = [r for r in rows if r["error"] is None and r["best_score"] is not None]
    ok.sort(key=lambda r: r["best_score"])
    for rank, r in enumerate(ok):
        r["rank"] = rank + 1

    # 重写排序后的 CSV
    with open(results_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in ok + [r for r in rows if r not in ok]:
            w.writerow(r)

    # 前 N 姿势转 PDB + 可选 PRODIGY-Lig 亲和力
    top_poses = []
    for r in ok[:top]:
        pose_pdb = os.path.join(outdir, "top_%02d_%s.pdb" % (r["rank"], _safe(r["id"])))
        write_pose_pdb(parse_first_model_lines(r["out_pdbqt"]), pose_pdb)
        item = {"rank": r["rank"], "id": r["id"], "smiles": r["smiles"],
                "score": r["best_score"], "pose_pdb": pose_pdb}
        if prodigy:
            try:
                from prodigy_affinity import predict_lig
                aff = predict_lig(receptor_pdb, pose_pdb)
                r["dg_kcal_mol"] = aff.get("dg_kcal_mol")
                r["kd_M"] = aff.get("kd_M")
                item["dg_kcal_mol"], item["kd_M"] = aff.get("dg_kcal_mol"), aff.get("kd_M")
            except Exception as e:
                item["prodigy_error"] = str(e)
        top_poses.append(item)

    report = {
        "outdir": outdir, "receptor": receptor_pdb, "receptor_pdbqt": rec_pdbqt,
        "n_ligands": len(ligands), "n_ok": len(ok), "n_failed": len(rows) - len(ok),
        "center": list(center), "size": list(size),
        "exhaustiveness": exhaustiveness, "results_csv": results_csv,
        "top": top_poses,
        "hits": [{"rank": r["rank"], "id": r["id"], "smiles": r["smiles"], "score": r["best_score"],
                  "dg_kcal_mol": r.get("dg_kcal_mol"), "kd_M": r.get("kd_M")} for r in ok[:top]],
    }
    return report


def _safe(s):
    return "".join(c if c.isalnum() else "_" for c in str(s))[:40]


def _print_report(r):
    print("\n=== virtual screening summary ===")
    print("receptor: %s | ligands: %d (ok %d, failed %d)" % (
        r["receptor"], r["n_ligands"], r["n_ok"], r["n_failed"]))
    print("center %s size %s | exhaustiveness %d" % (r["center"], r["size"], r["exhaustiveness"]))
    print("rank  id                       score      dG(kcal/mol)  Kd(M)")
    for h in r["hits"]:
        print("%4d  %-24s %8.2f   %-12s  %s" % (
            h["rank"], h["id"][:24], h["score"],
            ("%.2f" % h["dg_kcal_mol"]) if h.get("dg_kcal_mol") is not None else "-",
            ("%.3g" % h["kd_M"]) if h.get("kd_M") is not None else "-"))
    print("results: %s" % r["results_csv"])


def main(argv=None):
    ap = argparse.ArgumentParser(description="Vina batch virtual screening pipeline")
    ap.add_argument("--receptor", required=True, help="receptor PDB (auto rigid PDBQT)")
    ap.add_argument("--ligands", required=True, help="library CSV (one SMILES per row)")
    ap.add_argument("--smiles-col")
    ap.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"))
    ap.add_argument("--size", type=float, nargs=3, metavar=("X", "Y", "Z"))
    ap.add_argument("--ref-ligand", help="PDB of a co-crystallized ligand; box auto-derived")
    ap.add_argument("--exclude-res", default="HOH,WAT",
                    help="comma list of residue names excluded from receptor (add co-ligand resnames)")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=3)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--no-resume", action="store_true")
    ap.add_argument("--prodigy", action="store_true", help="run PRODIGY-Lig on top poses")
    ap.add_argument("--out", help="write JSON report to file")
    args = ap.parse_args(argv)

    if args.center and args.size:
        center, size = tuple(args.center), tuple(args.size)
    elif args.ref_ligand:
        center, size = box_from_ref_ligand(args.ref_ligand)
        print("box from %s: center %s size %s" % (args.ref_ligand, center, size))
    else:
        print("ERROR: provide --center/--size or --ref-ligand", file=sys.stderr)
        return 2
    try:
        report = screen(args.receptor, _load_library(args.ligands, args.smiles_col),
                        outdir=args.outdir, center=center, size=size,
                        exhaustiveness=args.exhaustiveness, n_poses=args.num_modes,
                        seed=args.seed, top=args.top,
                        exclude_res=tuple(x.strip() for x in args.exclude_res.split(",") if x.strip()),
                        resume=not args.no_resume, prodigy=args.prodigy)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    _print_report(report)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
        print("JSON written: %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
