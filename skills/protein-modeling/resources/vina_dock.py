#!/usr/bin/env python3
r"""vina_dock.py — AutoDock Vina 对接包装(配体准备 → 对接 → 姿势解析排序 → 报告)

运行环境(需要 meeko/rdkit):
  & 'D:\bioai\venv\Scripts\python.exe' vina_dock.py \
      --receptor-pdbqt rec.pdbqt --smiles "CCO" \
      --center 10 10 10 --size 20 20 20 --outdir D:\bioai\jobs\dock1
可选 Bio 接触摘要:先 $env:PYTHONPATH='D:\biopython' 再运行(venv 内也能 import Bio)。

组件:
  1. 配体准备:--smiles / --ligand(.sdf|.mol|.mol2) → rdkit → meeko → PDBQT
  2. 受体:--receptor-pdbqt 直接给 PDBQT(meeko rigid / ADFR 预先生成);
     --receptor 给 PDB 时仅用于 Bio.PDB 接触摘要(不自动转换)
  3. 对接:外部 vina 二进制(默认自动探测 D:\bioai\bin\vina*.exe,可用 --vina-bin 覆盖)
  4. 解析:REMARK VINA RESULT 打分 → 排序 → 导出前 --top 个姿势为 PDB
  5. 报告:JSON(分数/路径)+ 人类可读摘要;有 Biopython 时附受体-配体接触对(4 Å)

编程调用:
  from vina_dock import dock
  report = dock(receptor_pdbqt="rec.pdbqt", smiles="CCO", center=(10,10,10), size=(20,20,20))
"""
import argparse
import glob
import json
import os
import subprocess
import sys

DEFAULT_BIN_DIR = r"D:\bioai\bin"
DEFAULT_OUTDIR = r"D:\bioai\jobs"


# ── 1. 配体准备(rdkit + meeko)─────────────────────────────────────────────
def _meeko_write(mol):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    preparator = MoleculePreparation()
    out = preparator.prepare(mol)
    if isinstance(out, tuple):
        setups = out[0] if isinstance(out[0], list) else [out[0]]
    else:
        setups = out if isinstance(out, list) else [out]
    if not setups:
        raise ValueError("meeko produced no setup")
    result = PDBQTWriterLegacy.write_string(setups[0])
    if isinstance(result, tuple):
        pdbqt, is_ok, err = (result + (None,) * 3)[:3]
    else:
        pdbqt, is_ok, err = result, True, None
    if not is_ok or not pdbqt:
        raise ValueError("meeko write failed: %s" % err)
    return pdbqt


def smiles_to_pdbqt(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("cannot parse SMILES: %s" % smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
    AllChem.MMFFOptimizeMolecule(mol)
    return _meeko_write(mol)


def ligand_file_to_pdbqt(path):
    from rdkit import Chem
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".mol"):
        mols = Chem.SDMolSupplier(path, removeHs=False)
        mol = next((m for m in mols if m is not None), None)
    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(path, removeHs=False)
    else:
        raise ValueError("unsupported ligand format: %s (use .sdf/.mol/.mol2 or --smiles)" % path)
    if mol is None:
        raise ValueError("no molecule parsed from " + path)
    mol = Chem.AddHs(mol)
    return _meeko_write(mol)


# ── 2. 对接(外部 vina 二进制)──────────────────────────────────────────────
def find_vina_bin(bin_dir=DEFAULT_BIN_DIR):
    if os.path.isdir(bin_dir):
        hits = sorted(glob.glob(os.path.join(bin_dir, "vina*.exe")))
        if hits:
            return hits[0]
    return None


def run_vina(vina_bin, receptor_pdbqt, ligand_pdbqt, center, size, out_pdbqt,
             exhaustiveness=16, n_poses=9, seed=None, cpu=None):
    cmd = [vina_bin, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exhaustiveness), "--num_modes", str(n_poses),
           "--out", out_pdbqt]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        tail = ((r.stdout or "") + (r.stderr or ""))[-600:]
        raise RuntimeError("vina failed (rc=%s): %s" % (r.returncode, tail))
    return r.stdout


# ── 3. PDBQT 解析与姿势导出 ────────────────────────────────────────────────
def parse_pdbqt(pdbqt_path):
    models, cur, cur_score = [], None, None
    with open(pdbqt_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("MODEL"):
                cur, cur_score = [], None
            elif line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        cur_score = float(parts[3])
                    except ValueError:
                        pass
            elif line.startswith("ENDMDL") and cur is not None:
                models.append({"model": len(models) + 1, "score": cur_score, "lines": cur})
                cur = None
            elif cur is not None:
                cur.append(line)
    return models


def write_pose_pdb(model_lines, out_path):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("REMARK   converted from Vina PDBQT output\n")
        f.write("MODEL        1\n")
        for line in model_lines:
            if line.startswith(("ATOM", "HETATM")):
                f.write(line[:54].rstrip() + "\n")
        f.write("ENDMDL\n")
        f.write("END\n")


# ── 4. 可选 Bio.PDB 接触摘要 ──────────────────────────────────────────────
def bio_contact_summary(receptor_pdb, pose_pdb, cutoff=4.0):
    try:
        from Bio.PDB import PDBParser, NeighborSearch
    except ImportError:
        return None
    p = PDBParser(QUIET=True)
    try:
        rec = p.get_structure("rec", receptor_pdb)[0]
        lig = p.get_structure("lig", pose_pdb)[0]
    except Exception as e:
        return {"error": str(e)}
    rec_atoms = [a for a in rec.get_atoms() if a.get_parent().get_resname() != "HOH"]
    lig_atoms = [a for a in lig.get_atoms() if a.get_parent().get_resname() != "HOH"]
    if not rec_atoms or not lig_atoms:
        return None
    ns = NeighborSearch(rec_atoms)
    pairs = set()
    for a in lig_atoms:
        for b in ns.search(a.coord, cutoff):
            key = (b.get_parent().id[1], a.get_parent().id[1])
            pairs.add(key)
    return {"cutoff": cutoff, "n_contacts": len(pairs),
            "receptor_residues": sorted({p[0] for p in pairs})[:100]}


# ── 5. 主流程 ─────────────────────────────────────────────────────────────
def dock(receptor_pdbqt=None, receptor_pdb=None, smiles=None, ligand=None,
         ligand_pdbqt=None, center=None, size=None, outdir=DEFAULT_OUTDIR,
         exhaustiveness=16, n_poses=9, seed=None, cpu=None, top=5, vina_bin=None,
         name="dock"):
    if not center or len(center) != 3:
        raise ValueError("--center x y z is required")
    if not size or len(size) != 3:
        raise ValueError("--size x y z is required")
    if vina_bin is None:
        vina_bin = find_vina_bin()
    if not vina_bin or not os.path.exists(vina_bin):
        raise FileNotFoundError("vina binary not found; pass --vina-bin or put vina*.exe in %s" % DEFAULT_BIN_DIR)

    os.makedirs(outdir, exist_ok=True)
    if receptor_pdbqt is None:
        # 自动刚性受体准备(PDB -> PDBQT,元素映射 + 0.00 电荷)
        if receptor_pdb and os.path.exists(receptor_pdb):
            try:
                from pdb_to_pdbqt import pdb_to_pdbqt_text
            except ImportError:
                raise ValueError("--receptor PDB auto-prep needs pdb_to_pdbqt.py beside this script; "
                                 "or prepare the PDBQT yourself and pass --receptor-pdbqt")
            with open(receptor_pdb, encoding="utf-8") as f:
                pdbqt_text = pdb_to_pdbqt_text(f.read())
            receptor_pdbqt = os.path.join(outdir, name + "_receptor.pdbqt")
            with open(receptor_pdbqt, "w", encoding="utf-8") as f:
                f.write(pdbqt_text)
        else:
            raise ValueError("provide --receptor-pdbqt or a --receptor PDB path")
    if ligand_pdbqt is None:
        if smiles:
            ligand_pdbqt = _tmp_pdbqt(smiles_to_pdbqt(smiles), outdir, name)
        elif ligand:
            ligand_pdbqt = _tmp_pdbqt(ligand_file_to_pdbqt(ligand), outdir, name)
        else:
            raise ValueError("provide --smiles, --ligand, or --ligand-pdbqt")

    out_pdbqt = os.path.join(outdir, name + "_out.pdbqt")
    run_vina(vina_bin, receptor_pdbqt, ligand_pdbqt, center, size, out_pdbqt,
             exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed, cpu=cpu)

    models = parse_pdbqt(out_pdbqt)
    models.sort(key=lambda m: m["score"] if m["score"] is not None else float("inf"))
    poses = []
    for i, m in enumerate(models[:top]):
        pose_pdb = os.path.join(outdir, "%s_pose%d.pdb" % (name, i + 1))
        write_pose_pdb(m["lines"], pose_pdb)
        item = {"model": m["model"], "score": m["score"], "pdb": pose_pdb}
        if receptor_pdb and os.path.exists(receptor_pdb):
            item["contacts"] = bio_contact_summary(receptor_pdb, pose_pdb)
        poses.append(item)

    report = {
        "name": name,
        "vina_bin": vina_bin,
        "receptor": receptor_pdbqt,
        "ligand": ligand_pdbqt,
        "center": center, "size": size,
        "exhaustiveness": exhaustiveness,
        "n_poses_found": len(models),
        "best_score": poses[0]["score"] if poses else None,
        "poses": poses,
        "output_pdbqt": out_pdbqt,
    }
    return report


def _tmp_pdbqt(text, outdir, name):
    path = os.path.join(outdir, name + "_ligand.pdbqt")
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return path


def _print_report(r):
    print("vina: %s" % r["vina_bin"])
    print("receptor: %s | center %s | size %s" % (r["receptor"], r["center"], r["size"]))
    print("poses found: %d | best score: %s" % (r["n_poses_found"], r["best_score"]))
    for p in r["poses"]:
        extra = ""
        if p.get("contacts"):
            c = p["contacts"]
            if "error" in c:
                extra = " (contacts: %s)" % c["error"]
            else:
                extra = " (contacts: %d, rec residues: %s)" % (c["n_contacts"], ",".join(map(str, c["receptor_residues"][:8])))
        print("  model %2d  score %7.2f  %s%s" % (p["model"], p["score"], os.path.basename(p["pdb"]), extra))


def main(argv=None):
    ap = argparse.ArgumentParser(description="AutoDock Vina docking wrapper (meeko prep + ranking)")
    ap.add_argument("--receptor-pdbqt", help="prepared receptor PDBQT")
    ap.add_argument("--receptor", help="receptor PDB (for Bio.PDB contact summary only)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--smiles", help="ligand SMILES")
    src.add_argument("--ligand", help="ligand file .sdf/.mol/.mol2")
    src.add_argument("--ligand-pdbqt", help="pre-made ligand PDBQT")
    ap.add_argument("--center", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"))
    ap.add_argument("--size", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"))
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    ap.add_argument("--name", default="dock")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--n-poses", type=int, default=9)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--cpu", type=int)
    ap.add_argument("--vina-bin")
    ap.add_argument("--out", help="write JSON report to file")
    args = ap.parse_args(argv)

    try:
        report = dock(receptor_pdbqt=args.receptor_pdbqt, receptor_pdb=args.receptor,
                      smiles=args.smiles, ligand=args.ligand, ligand_pdbqt=args.ligand_pdbqt,
                      center=tuple(args.center), size=tuple(args.size), outdir=args.outdir,
                      exhaustiveness=args.exhaustiveness, n_poses=args.n_poses, top=args.top,
                      seed=args.seed, cpu=args.cpu, vina_bin=args.vina_bin, name=args.name)
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
