#!/usr/bin/env python3
r"""mol_tools.py — RDKit 化学信息学套件(描述 / 标准化 / 相似性 / 子结构 / 构象 / 2D 图)

运行环境(需要 rdkit,本机 venv 已装):
  & 'D:\bioai\venv\Scripts\python.exe' mol_tools.py <子命令> ...

子命令:
  describe     SMILES → 物化性质表(Lipinski/QED/TPSA/可旋转键...),CSV 或 JSON
               mol_tools.py describe --smiles "CC(=O)Oc1ccccc1C(=O)O" --json
               mol_tools.py describe --library lib.csv --smiles-col smiles --out props.csv
  canonical    标准化 SMILES(去盐、规范互变异构、规范形式)
  depict       SMILES → 2D PNG(可 --highlight SMARTS)
  similarity   查询分子 vs 化合物库,输出按 Tanimoto 排序的 CSV(Morgan2/AtomPair)
  substructure SMARTS 子结构匹配(计数 + 高亮 PNG)
  conformers   ETKDGv3 生成 N 个构象 + MMFF 优化 → SDF

编程调用:
  from mol_tools import describe_mol, similarity_search, gen_conformers
"""
import argparse
import csv
import json
import os
import sys


# ── 解析与性质 ───────────────────────────────────────────────────────────────
def _rdmol(smiles, add_h=False):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("cannot parse SMILES: %s" % smiles)
    if add_h:
        mol = Chem.AddHs(mol)
    return mol


def canonical_smiles(smiles):
    """去盐 + 规范互变异构 + 规范 SMILES。"""
    from rdkit import Chem
    try:
        from rdkit.Chem.SaltRemover import SaltRemover   # RDKit >= 2023.09
    except ImportError:
        from rdkit.Chem import SaltRemover               # 旧版
    remover = SaltRemover()
    mol = _rdmol(smiles)
    stripped = remover.StripMol(mol)
    if stripped is None or stripped.GetNumAtoms() == 0:
        # 裸金属离子(如 [Na+])不在盐库:保留最大有机片段
        frags = Chem.GetMolFrags(mol, asMols=True)
        frags = [f for f in frags if f.GetNumAtoms() > 0]
        mol = max(frags, key=lambda f: sum(1 for a in f.GetAtoms() if a.GetAtomicNum() == 6)) if frags else mol
    else:
        mol = stripped
    Chem.RemoveStereochemistry(mol)
    from rdkit.Chem.MolStandardize import rdMolStandardize
    te = rdMolStandardize.TautomerEnumerator()
    mol = te.Canonicalize(mol)
    return Chem.MolToSmiles(mol)


def _lipinski_violations(mw, logp, hbd, hba):
    v = 0
    if mw > 500:
        v += 1
    if logp > 5:
        v += 1
    if hbd > 5:
        v += 1
    if hba > 10:
        v += 1
    return v


def describe_mol(smiles):
    from rdkit.Chem import Descriptors, QED, Lipinski
    from rdkit.Chem import Crippen
    mol = _rdmol(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    rings = Descriptors.RingCount(mol)
    aro = sum(1 for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
    return {
        "smiles": canonical_smiles(smiles),
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "hbd": int(hbd),
        "hba": int(hba),
        "rotatable_bonds": int(rot),
        "rings": int(rings),
        "aromatic_rings": int(aro),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "lipinski_violations": _lipinski_violations(mw, logp, hbd, hba),
        "lipinski_pass": _lipinski_violations(mw, logp, hbd, hba) <= 1,
        "qed": round(float(QED.qed(mol)), 3),
        "formula": Descriptors.MolecularFormula(mol) if hasattr(Descriptors, "MolecularFormula") else None,
    }


def _fingerprint(mol, ftype="morgan2", nbits=2048):
    from rdkit.Chem import rdFingerprintGenerator, AllChem
    if ftype == "morgan2":
        try:
            return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=nbits).GetFingerprint(mol)
        except Exception:
            return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits)
    if ftype == "morgan3":
        try:
            return rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=nbits).GetFingerprint(mol)
        except Exception:
            return AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=nbits)
    if ftype == "atompair":
        try:
            return rdFingerprintGenerator.GetAtomPairGenerator(fpSize=nbits).GetFingerprint(mol)
        except Exception:
            return AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=nbits)
    raise ValueError("unknown fptype: %s" % ftype)


def similarity_search(query_smiles, library, smiles_col=None, ftype="morgan2", nbits=2048):
    """library: CSV 路径或 [(id, smiles)]。返回按 Tanimoto 降序的列表。"""
    from rdkit import DataStructs
    q = _rdmol(query_smiles)
    qfp = _fingerprint(q, ftype, nbits)
    if isinstance(library, str):
        rows = []
        with open(library, encoding="utf-8-sig", newline="") as f:
            for rec in csv.DictReader(f):
                if smiles_col is None:
                    keys = [k for k in rec if "smiles" in k.lower()]
                    smiles_col = keys[0] if keys else None
                if smiles_col is None:
                    raise ValueError("--smiles-col required (no 'smiles' column found)")
                rows.append((rec.get("id") or rec.get("name") or rec.get(smiles_col), rec.get(smiles_col, "")))
    else:
        rows = list(library)
    out = []
    for rid, smi in rows:
        try:
            m = _rdmol(smi)
            sim = DataStructs.TanimotoSimilarity(qfp, _fingerprint(m, ftype, nbits))
            out.append({"id": rid, "smiles": canonical_smiles(smi), "tanimoto": round(sim, 4)})
        except Exception as e:
            out.append({"id": rid, "smiles": smi, "tanimoto": None, "error": str(e)})
    out.sort(key=lambda r: -(r["tanimoto"] if r["tanimoto"] is not None else -1))
    return out


def substructure(smiles, smarts, out_png=None):
    from rdkit import Chem
    mol = _rdmol(smiles)
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise ValueError("cannot parse SMARTS: %s" % smarts)
    matches = mol.GetSubstructMatches(patt)
    flat = sorted({i for m in matches for i in m})
    if out_png:
        from rdkit.Chem.Draw import MolToFile
        if matches:
            MolToFile(mol, out_png, size=(500, 400), highlightAtoms=flat)
        else:
            MolToFile(mol, out_png, size=(500, 400))
    return {"smiles": canonical_smiles(smiles), "smarts": smarts,
            "n_matches": len(matches), "matched_atom_indices": [list(m) for m in matches[:20]],
            "png": out_png}


def gen_conformers(smiles, n=10, out_sdf=None, seed=0xF00D):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("cannot parse SMILES: %s" % smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n, params=params)
    if not cids:
        raise ValueError("conformer embedding failed for %s" % smiles)
    props = AllChem.MMFFGetMoleculeProperties(mol)
    energies = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        if ff is None:
            energies.append(None)
            continue
        ff.Minimize()
        energies.append(ff.CalcEnergy())
    # 按能量排序并写出
    order = sorted(range(len(cids)), key=lambda i: energies[i] if energies[i] is not None else 1e18)
    if out_sdf:
        w = Chem.SDWriter(out_sdf)
        for rank, i in enumerate(order):
            mol.SetProp("_Name", "conf_%d" % (rank + 1))
            mol.SetProp("MMFF_energy", "%.4f" % (energies[i] if energies[i] is not None else float("nan")))
            w.write(mol, confId=cids[i])
        w.close()
    return {"smiles": canonical_smiles(smiles), "n_confs": len(cids),
            "mmff_energies_kcal": [round(e, 3) if e is not None else None for i, e in zip(order, [energies[i] for i in order])],
            "sdf": out_sdf}


def depict(smiles, out_png, highlight_smarts=None, size=(500, 400)):
    from rdkit.Chem.Draw import MolToFile
    mol = _rdmol(smiles)
    atoms = []
    if highlight_smarts:
        from rdkit import Chem
        patt = Chem.MolFromSmarts(highlight_smarts)
        if patt is not None:
            atoms = sorted({i for m in mol.GetSubstructMatches(patt) for i in m})
    MolToFile(mol, out_png, size=size, highlightAtoms=atoms or None)
    return out_png


# ── 主程序 ──────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(description="RDKit cheminformatics toolbox")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_desc = sub.add_parser("describe", help="molecular properties")
    p_desc.add_argument("--smiles")
    p_desc.add_argument("--library", help="CSV file of SMILES")
    p_desc.add_argument("--smiles-col", help="column name for SMILES")
    p_desc.add_argument("--json", action="store_true")
    p_desc.add_argument("--out")

    p_can = sub.add_parser("canonical", help="canonical SMILES")
    p_can.add_argument("smiles")

    p_dep = sub.add_parser("depict", help="2D depiction PNG")
    p_dep.add_argument("--smiles", required=True)
    p_dep.add_argument("--out", required=True)
    p_dep.add_argument("--highlight")

    p_sim = sub.add_parser("similarity", help="Tanimoto similarity search")
    p_sim.add_argument("--query", required=True)
    p_sim.add_argument("--library", required=True)
    p_sim.add_argument("--smiles-col")
    p_sim.add_argument("--fptype", default="morgan2", choices=["morgan2", "morgan3", "atompair"])
    p_sim.add_argument("--top", type=int, default=20)
    p_sim.add_argument("--out")

    p_sub = sub.add_parser("substructure", help="SMARTS substructure match")
    p_sub.add_argument("--smiles", required=True)
    p_sub.add_argument("--smarts", required=True)
    p_sub.add_argument("--png")

    p_conf = sub.add_parser("conformers", help="ETKDG conformer generation")
    p_conf.add_argument("--smiles", required=True)
    p_conf.add_argument("--n", type=int, default=10)
    p_conf.add_argument("--out", required=True)

    args = ap.parse_args(argv)
    try:
        if args.cmd == "describe":
            if args.smiles:
                rows = [describe_mol(args.smiles)]
            elif args.library:
                rows = []
                with open(args.library, encoding="utf-8-sig", newline="") as f:
                    for rec in csv.DictReader(f):
                        col = args.smiles_col or next((k for k in rec if "smiles" in k.lower()), None)
                        if col is None:
                            raise ValueError("--smiles-col required")
                        try:
                            r = describe_mol(rec[col])
                            r["id"] = rec.get("id") or rec.get("name") or rec[col]
                            rows.append(r)
                        except Exception as e:
                            rows.append({"id": rec.get("id", "?"), "smiles": rec.get(col, ""), "error": str(e)})
            else:
                raise ValueError("--smiles or --library required")
            if args.json:
                print(json.dumps(rows, ensure_ascii=False, indent=2))
            else:
                keys = ["id", "smiles", "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds",
                        "rings", "aromatic_rings", "lipinski_violations", "lipinski_pass", "qed"]
                print("\t".join(keys))
                for r in rows:
                    print("\t".join(str(r.get(k, "")) for k in keys))
            if args.out:
                with open(args.out, "w", encoding="utf-8", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                    w.writeheader()
                    w.writerows(rows)
                print("CSV written: %s" % args.out)
            return 0

        if args.cmd == "canonical":
            print(canonical_smiles(args.smiles))
            return 0

        if args.cmd == "depict":
            depict(args.smiles, args.out, highlight_smarts=args.highlight)
            print("PNG written: %s" % args.out)
            return 0

        if args.cmd == "similarity":
            rows = similarity_search(args.query, args.library, smiles_col=args.smiles_col,
                                     ftype=args.fptype)
            shown = [r for r in rows if r["tanimoto"] is not None][:args.top]
            for r in shown:
                print("%-6.4f  %s  %s" % (r["tanimoto"], str(r["id"]), r["smiles"]))
            n_err = sum(1 for r in rows if r["tanimoto"] is None)
            if n_err:
                print("(%d rows failed to parse)" % n_err)
            if args.out:
                with open(args.out, "w", encoding="utf-8", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=["rank", "id", "smiles", "tanimoto"])
                    w.writeheader()
                    for i, r in enumerate(rows):
                        w.writerow({"rank": i + 1, "id": r["id"], "smiles": r["smiles"],
                                    "tanimoto": r["tanimoto"] if r["tanimoto"] is not None else ""})
                print("CSV written: %s" % args.out)
            return 0

        if args.cmd == "substructure":
            r = substructure(args.smiles, args.smarts, out_png=args.png)
            print(json.dumps({k: v for k, v in r.items() if k != "png"}, ensure_ascii=False, indent=2))
            if args.png:
                print("PNG written: %s" % args.png)
            return 0

        if args.cmd == "conformers":
            r = gen_conformers(args.smiles, n=args.n, out_sdf=args.out)
            print(json.dumps(r, ensure_ascii=False, indent=2))
            return 0
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
