#!/usr/bin/env python3
r"""md_mmgbsa.py — OpenMM MM-GBSA 结合自由能 + 显式溶剂 MD 协议(加热/平衡/采样 + RMSD/RMSF)

运行环境(D:\bioai\venv:openmm/numpy/matplotlib):
  # MM-GBSA(隐式溶剂,分钟级,蛋白-蛋白或蛋白单链稳定性):
  & 'D:\bioai\venv\Scripts\python.exe' md_mmgbsa.py --mode gb --complex complex.pdb \
      --rec-chains A --lig-chains B --out D:\bioai\jobs\mmgbsa1.json
  # 显式溶剂 MD(完整协议;--steps 控制采样长度):
  & 'D:\bioai\venv\Scripts\python.exe' md_mmgbsa.py --mode md --complex complex.pdb \
      --steps 10000 --outdir D:\bioai\jobs\md1 --platform CPU

输出:
  gb 模式:ΔG_bind = E(complex) - E(receptor) - E(ligand),分解为 internal(键/角/二面)
          与非键(范德华+GB+SA)两项;JSON 报告。绝对值偏"过稳定",用于排序/对比而非绝对 ΔG。
  md 模式:trajectory.dcd + state.csv(温度/势能/体积)+ rmsd.csv/rmsf.csv + PNG 图。

注意事项(本机实测配方,见 protein-modeling 技能):
  - MSE 无模板:脚本内部自动 MSE→MET、 SE → SD
  - 结晶水无模板:自动 deleteWater;缺氢自动 addHydrogens
  - 显式溶剂默认 TIP3P 10 Å 盒 + 0.15 M NaCl 中和
  - 8GB 笔记本 GPU 上 MD 建议 CPU 平台(确定性好、无显存抖动);--platform OpenCL 可用
"""
import argparse
import json
import os
import sys
import time

import numpy as np


# 标准氨基酸完整重原子集(缺侧链原子的晶体残基会被过滤,避免 OpenMM 模板错配)
AA_HEAVY = {
    "GLY": {"N", "CA", "C", "O"},
    "ALA": {"N", "CA", "C", "O", "CB"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
}


def load_pdb(path):
    import openmm.app as app
    with open(path, encoding="utf-8") as f:
        raw = f.read()
    fixed = raw.replace("MSE", "MET").replace(" SE ", " SD ")
    # Bio.PDB 清洗:去 HETATM(配体/离子/水)、无序 altloc、原子集不完整的残基;
    # 断口处自动补 OXT 封端(被剔除残基的邻位成为人工末端时 OpenMM 需要 OXT)
    import io
    import numpy as _np
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", io.StringIO(fixed))

    def _keep(res):
        if res.id[0] != " ":
            return False
        expected = AA_HEAVY.get(res.get_resname().strip().upper())
        if expected is None:
            return False
        names = {a.get_name().strip() for a in res}
        return expected <= names

    lines, serial = [], 0
    for model in struct:
        for chain in model:
            kept = [r for r in chain if _keep(r)]
            # 相邻保留残基的骨架连续性(C-N 距离)
            bonded = []
            for i in range(len(kept) - 1):
                a, b = kept[i], kept[i + 1]
                if "C" in a and "N" in b:
                    bonded.append(float(_np.linalg.norm(_np.asarray(a["C"].coord) -
                                                        _np.asarray(b["N"].coord))) < 2.2)
                else:
                    bonded.append(False)
            # 丢弃孤立残基(两侧均断:无法同时匹配 N 端与 C 端模板)
            final = []
            for i, res in enumerate(kept):
                has_prev = (i > 0 and bonded[i - 1]) or i == 0
                has_next = (i < len(kept) - 1 and bonded[i]) or i == len(kept) - 1
                if has_prev or has_next:
                    final.append(res)
            if not final:
                continue
            # 重算 final 内部连续性
            f_bonded = []
            for i in range(len(final) - 1):
                a, b = final[i], final[i + 1]
                if "C" in a and "N" in b:
                    f_bonded.append(float(_np.linalg.norm(_np.asarray(a["C"].coord) -
                                                          _np.asarray(b["N"].coord))) < 2.2)
                else:
                    f_bonded.append(False)
            prev_bonded = True
            for idx, res in enumerate(final):
                bonded_next = f_bonded[idx] if idx < len(f_bonded) else False
                if not prev_bonded:
                    lines.append("TER")   # 断口处 TER,防止 PDBFile 跨断口建肽键
                prev_bonded = bonded_next
                real_end = not bonded_next
                for a in res.get_atoms():
                    if a.element == "H":
                        continue
                    if a.get_name().strip() == "OXT" and not real_end:
                        continue   # 晶体常在中链残基带 OXT,会导致模板误判为 C 端
                    serial += 1
                    x, y, z = a.coord
                    lines.append(
                        "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
                        % (serial, a.get_name(), res.get_resname(), chain.id,
                           res.id[1], x, y, z, 1.0, 0.0, a.element))
                if real_end and "OXT" not in {a.get_name().strip() for a in res}:
                    # 近似 OXT 位置(C 四面体对侧),后续最小化会修正
                    ca = res["CA"].coord
                    c = res["C"].coord
                    o = res["O"].coord
                    v = _np.asarray(ca) - _np.asarray(c) + _np.asarray(o) - _np.asarray(c)
                    n = v / max(_np.linalg.norm(v), 1e-8)
                    oxt = _np.asarray(c) - 1.25 * n
                    serial += 1
                    lines.append(
                        "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           O"
                        % (serial, "OXT", res.get_resname(), chain.id,
                           res.id[1], oxt[0], oxt[1], oxt[2], 1.0, 0.0))
    pdb_text = "\n".join(lines) + "\nEND\n"
    return app.PDBFile(io.StringIO(pdb_text))


def build_modeller(pdb, ff, keep_chains=None):
    import openmm.app as app
    import openmm.unit as u
    m = app.Modeller(pdb.topology, pdb.positions)
    m.deleteWater()
    if keep_chains is not None:
        keep = set(keep_chains)
        dels = [c for c in m.topology.chains() if c.id not in keep]
        m.delete(dels)
    m.addHydrogens(ff)
    return m


# ── MM-GBSA(隐式 OBC2)───────────────────────────────────────────────────────
def _gb_energy(pdb_path, chains, platform_name="CPU"):
    import openmm
    import openmm.app as app
    import openmm.unit as u
    pdb = load_pdb(pdb_path)
    ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
    m = build_modeller(pdb, ff, keep_chains=chains)
    system = ff.createSystem(m.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
    integrator = openmm.LangevinMiddleIntegrator(300 * u.kelvin, 1 / u.picosecond,
                                                 0.002 * u.picoseconds)
    platform = openmm.Platform.getPlatformByName(platform_name)
    sim = app.Simulation(m.topology, system, integrator, platform)
    sim.context.setPositions(m.positions)
    sim.minimizeEnergy(maxIterations=500)
    state = sim.context.getState(getEnergy=True)
    total = state.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    # 分组:internal = 键合力;nonbonded = GBSAOBC/Nonbonded 力
    nb_idx = [i for i, f in enumerate(system.getForces())
              if "OBC" in f.__class__.__name__ or "Nonbonded" in f.__class__.__name__]
    internal = 0.0
    if nb_idx:
        s_int = sim.context.getState(getEnergy=True, groups=set(range(len(system.getForces()))) - set(nb_idx))
        internal = s_int.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    return {"total": total, "internal": internal, "nonbonded": total - internal,
            "n_atoms": m.topology.getNumAtoms()}


def mmgbsa(complex_pdb, rec_chains, lig_chains, platform="CPU"):
    """ΔG_bind = E_complex - E_rec - E_lig(OBC2 隐式)。返回报告 dict。"""
    t0 = time.time()
    e_complex = _gb_energy(complex_pdb, chains=None, platform_name=platform)
    e_rec = _gb_energy(complex_pdb, chains=rec_chains, platform_name=platform)
    e_lig = _gb_energy(complex_pdb, chains=lig_chains, platform_name=platform)
    report = {
        "mode": "mm-gbsa (OBC2 implicit)",
        "complex": complex_pdb, "rec_chains": list(rec_chains), "lig_chains": list(lig_chains),
        "platform": platform,
        "e_complex_kcal": {k: round(v, 2) for k, v in e_complex.items() if k != "n_atoms"},
        "e_receptor_kcal": {k: round(v, 2) for k, v in e_rec.items() if k != "n_atoms"},
        "e_ligand_kcal": {k: round(v, 2) for k, v in e_lig.items() if k != "n_atoms"},
        "n_atoms": {"complex": e_complex["n_atoms"], "receptor": e_rec["n_atoms"],
                    "ligand": e_lig["n_atoms"]},
        "dg_bind_kcal_mol": round(e_complex["total"] - e_rec["total"] - e_lig["total"], 2),
        "dg_internal_kcal_mol": round(e_complex["internal"] - e_rec["internal"] - e_lig["internal"], 2),
        "dg_nonbonded_kcal_mol": round(e_complex["nonbonded"] - e_rec["nonbonded"] - e_lig["nonbonded"], 2),
        "elapsed_s": round(time.time() - t0, 1),
        "note": "implicit-solvent single-point energies; absolute values over-stabilize — use for ranking/comparison, not absolute dG",
    }
    return report


# ── 显式溶剂 MD 协议 ─────────────────────────────────────────────────────────
def run_md(complex_pdb, outdir, steps=100000, platform="CPU", dt_ps=0.002,
           report_interval=1000, heating_steps=10000, equil_steps=50000):
    import openmm
    import openmm.app as app
    import openmm.unit as u
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    pdb = load_pdb(complex_pdb)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")   # tip3p 文件含 Na+/Cl- 离子模板
    m = app.Modeller(pdb.topology, pdb.positions)
    m.deleteWater()
    m.addHydrogens(ff)
    m.addSolvent(ff, padding=1.0 * u.nanometer, ionicStrength=0.15 * u.molar)
    system = ff.createSystem(m.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * u.nanometer, constraints=app.HBonds)
    integrator = openmm.LangevinMiddleIntegrator(300 * u.kelvin, 1 / u.picosecond,
                                                 dt_ps * u.picoseconds)
    platform_obj = openmm.Platform.getPlatformByName(platform)
    sim = app.Simulation(m.topology, system, integrator, platform_obj)
    sim.context.setPositions(m.positions)
    sim.minimizeEnergy(maxIterations=1000)

    log_path = os.path.join(outdir, "md.log")
    log = open(log_path, "w", encoding="utf-8")
    def _log(msg):
        print(msg)
        log.write(msg + "\n")
        log.flush()

    _log("minimized; heating %d steps 0->300K ..." % heating_steps)
    for i in range(heating_steps):
        integrator.setTemperature(300.0 * (i + 1) / heating_steps * u.kelvin)
        sim.step(1)
    _log("equilibrating %d steps (NVT, 300K) ..." % equil_steps)
    sim.step(equil_steps)

    dcd = os.path.join(outdir, "trajectory.dcd")
    csv_path = os.path.join(outdir, "state.csv")
    report_interval = min(report_interval, max(1, steps // 10))   # 短运行也保证 >=10 帧
    sim.reporters.append(app.DCDReporter(dcd, report_interval))
    sim.reporters.append(app.StateDataReporter(csv_path, report_interval, step=True,
                                               temperature=True, potentialEnergy=True,
                                               volume=True, density=True))
    _log("production: %d steps (%g ps) ..." % (steps, steps * dt_ps))
    for i in range(0, steps, report_interval):
        sim.step(min(report_interval, steps - i))
        if (i // report_interval) % 10 == 0:
            _log("  step %d/%d (%.0fs)" % (i, steps, time.time() - t0))
    _log("done in %.0f s" % (time.time() - t0))
    sim.reporters.clear()   # 强制报告器析构/落盘,避免残尾帧

    report = {"mode": "md-explicit", "complex": complex_pdb, "platform": platform,
              "steps": steps, "dt_ps": dt_ps, "dcd": dcd, "log": log_path,
              "state_csv": csv_path,
              "n_atoms": m.topology.getNumAtoms(), "elapsed_s": round(time.time() - t0, 1)}
    try:
        ana = analyze_trajectory(m.topology, dcd, outdir, dt_ps=dt_ps,
                                 report_interval=report_interval)
        report["analysis"] = ana
    except Exception as e:
        report["analysis_error"] = str(e)
        _log("analysis failed: %s" % e)
    return report


def _sanitize_dcd(path):
    """按 CHARMM DCD 布局解析头部,截掉末帧残缺尾(OpenMM DCDReporter 析构前
    可能把部分缓冲写入文件,导致 mdtraj 'premature end of file')。返回完整帧数。"""
    import struct
    with open(path, "rb") as f:
        raw = f.read()
    if len(raw) < 92:
        return 0
    o = 4 + 4 + 80 + 4                    # block1: marker + 'CORD'+80B + end
    m2 = struct.unpack("<i", raw[o:o + 4])[0]
    o += 4 + m2 + 4                       # block2 整块跳过
    o += 4                                 # block3 marker
    natom = struct.unpack("<i", raw[o:o + 4])[0]
    o += 8                                 # natom + end
    m4 = struct.unpack("<i", raw[o:o + 4])[0]
    o += 4 + m4 + 4                        # block4(晶胞信息)
    unit = natom * 12 + 8
    frames = (len(raw) - o) // unit
    tail = (len(raw) - o) % unit
    if tail and frames > 0:
        with open(path, "r+b") as f:
            f.truncate(o + frames * unit)
    return frames


def analyze_trajectory(topology, dcd, outdir, dt_ps=0.002, report_interval=1000):
    """mdtraj 读取 DCD + 对齐 + RMSD/RMSF(OpenMM DCDFile 只是写入器,parmed 4.3.1 读不了 DCD)。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mdtraj as md

    _sanitize_dcd(dcd)
    top = md.Topology.from_openmm(topology)
    traj = md.load(dcd, top=top)
    ca = traj.topology.select("name CA")
    if len(ca) < 3:
        raise RuntimeError("no CA atoms in topology")
    n_frames = traj.n_frames
    if n_frames < 2:
        raise RuntimeError("trajectory has <2 frames")
    traj.superpose(traj, 0, atom_indices=ca)
    rmsd = md.rmsd(traj, traj, 0, atom_indices=ca) * 10.0   # nm -> Å
    rmsf = md.rmsf(traj, traj, atom_indices=ca) * 10.0
    reskeys = [str(r.resSeq) + r.name for r in traj.topology.residues
               if any(a.name == "CA" for a in r.atoms)]

    rmsd_csv = os.path.join(outdir, "rmsd.csv")
    with open(rmsd_csv, "w", encoding="utf-8") as f:
        f.write("frame,rmsd_ca_A\n")
        for i, v in enumerate(rmsd):
            f.write("%d,%.3f\n" % (i, float(v)))
    rmsf_csv = os.path.join(outdir, "rmsf.csv")
    with open(rmsf_csv, "w", encoding="utf-8") as f:
        f.write("residue,rmsf_A\n")
        for k, v in zip(reskeys, rmsf):
            f.write("%s,%.3f\n" % (k, float(v)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(rmsd)
    axes[0].set_xlabel("frame"); axes[0].set_ylabel("RMSD (A)"); axes[0].set_title("CA RMSD")
    axes[1].plot(rmsf)
    axes[1].set_xlabel("residue index"); axes[1].set_ylabel("RMSF (A)"); axes[1].set_title("CA RMSF")
    fig.tight_layout()
    png = os.path.join(outdir, "rmsd_rmsf.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return {"n_frames": n_frames, "rmsd_csv": rmsd_csv, "rmsf_csv": rmsf_csv,
            "plot": png, "rmsd_mean": round(float(np.mean(rmsd)), 3),
            "rmsd_last": round(float(rmsd[-1]), 3),
            "rmsf_max": round(float(np.max(rmsf)), 3),
            "rmsf_max_residue": reskeys[int(np.argmax(rmsf))]}


def main(argv=None):
    ap = argparse.ArgumentParser(description="OpenMM MM-GBSA / explicit MD protocol")
    ap.add_argument("--mode", choices=["gb", "md"], required=True)
    ap.add_argument("--complex", required=True)
    ap.add_argument("--rec-chains", nargs="+", default=["A"])
    ap.add_argument("--lig-chains", nargs="+", default=["B"])
    ap.add_argument("--platform", default="CPU", choices=["CPU", "OpenCL", "CUDA"])
    ap.add_argument("--steps", type=int, default=100000)
    ap.add_argument("--dt-ps", type=float, default=0.002)
    ap.add_argument("--heating-steps", type=int, default=10000)
    ap.add_argument("--equil-steps", type=int, default=50000)
    ap.add_argument("--out", help="gb: JSON report path")
    ap.add_argument("--outdir", help="md: output directory")
    args = ap.parse_args(argv)
    try:
        if args.mode == "gb":
            report = mmgbsa(args.complex, tuple(args.rec_chains), tuple(args.lig_chains),
                            platform=args.platform)
            print("dG_bind = %.2f kcal/mol (internal %.2f | nonbonded %.2f)"
                  % (report["dg_bind_kcal_mol"], report["dg_internal_kcal_mol"],
                     report["dg_nonbonded_kcal_mol"]))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(report, f, ensure_ascii=False, indent=2)
                print("JSON written: %s" % args.out)
        else:
            outdir = args.outdir or os.path.join(os.path.dirname(args.complex), "md_out")
            report = run_md(args.complex, outdir, steps=args.steps, platform=args.platform,
                            dt_ps=args.dt_ps, heating_steps=args.heating_steps,
                            equil_steps=args.equil_steps)
            if report.get("analysis"):
                a = report["analysis"]
                print("MD done: RMSD mean=%.2f last=%.2f | RMSF max=%.2f (%s) | plot %s"
                      % (a["rmsd_mean"], a["rmsd_last"], a["rmsf_max"], a["rmsf_max_residue"], a["plot"]))
            with open(os.path.join(outdir, "md_report.json"), "w", encoding="utf-8") as f:
                json.dump(report, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
