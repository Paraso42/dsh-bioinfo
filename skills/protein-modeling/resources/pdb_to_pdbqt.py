#!/usr/bin/env python3
r"""pdb_to_pdbqt.py — 刚性受体 PDB -> PDBQT 转换器(纯标准库)

AutoDock Vina 受体准备的最小实现:读 PDB 的 ATOM/HETATM 行,输出 PDBQT
(坐标 + 占据率 + 温度因子 + 电荷 0.00 + AutoDock4 原子类型)。
排除水(HOH/WAT)。刚性对接(rigid receptor)不依赖精确电荷,
Vina 用自己的打分函数,0.00 电荷 + 元素映射类型足够。

用法:
  & 'D:\bioai\venv\Scripts\python.exe' pdb_to_pdbqt.py receptor.pdb -o receptor.pdbqt
编程:
  from pdb_to_pdbqt import pdb_to_pdbqt_text
  text = pdb_to_pdbqt_text(open("receptor.pdb", encoding="utf-8").read())
"""
import argparse
import sys

# 元素 -> AutoDock4 原子类型(刚性受体最简映射)
AD_TYPE = {
    "C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "H": "HD", "FE": "Fe", "ZN": "Zn", "MG": "Mg", "CA": "CA",
    "MN": "Mn", "CO": "Co", "NA": "Na", "K": "K",
}
SKIP_RES = {"HOH", "WAT", "H2O"}


def _elem_of(line):
    """从 PDB ATOM/HETATM 行取元素符号(第 77-78 列,回退原子名推断)。"""
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


def pdb_to_pdbqt_text(pdb_text, model=0):
    lines = []
    cur_model = 0
    serial = 0
    n_atoms = 0
    for line in pdb_text.splitlines():
        if line.startswith("MODEL"):
            cur_model = int(line[10:14] or 1) - 1
            continue
        if cur_model != model:
            continue
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        rec = line[:6]
        name = line[12:16]
        resname = line[17:20].strip()
        if resname in SKIP_RES:
            continue
        chain = line[21:22]
        resseq = line[22:26].strip()
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
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
        elem = _elem_of(line)
        atype = AD_TYPE.get(elem) or "C"   # 未知元素(如 Se)回退 C,保证 Vina 可解析
        serial += 1
        lines.append(
            "%s%5d %-4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f      0.00 %-2s"
            % (rec, serial, name, resname, chain, resseq, x, y, z, occ, bf, atype)
        )
        n_atoms += 1
    if n_atoms == 0:
        raise ValueError("no ATOM/HETATM records parsed from PDB")
    header = (
        "REMARK    receptor prepared by pdb_to_pdbqt.py (rigid, charges 0.00)\n"
        "REMARK    AD4 atom types by element mapping\n"
    )
    return header + "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description="rigid receptor PDB -> PDBQT (stdlib)")
    ap.add_argument("pdb")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args(argv)
    with open(args.pdb, encoding="utf-8") as f:
        text = f.read()
    out = pdb_to_pdbqt_text(text)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(out)
    n = sum(1 for l in out.splitlines() if l.startswith(("ATOM", "HETATM")))
    print("wrote %s (%d atoms)" % (args.out, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
