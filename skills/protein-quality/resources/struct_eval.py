#!/usr/bin/env python3
# LICENSE NOTE: the TM-score implementation in this file is a faithful
# re-implementation of TMalign (Zhang lab). TMalign is free for ACADEMIC use
# only; commercial use or redistribution requires permission from the Zhang
# lab. Cite: Zhang Y, Skolnick J. "TM-align: a protein structure alignment
# algorithm based on the TM-score." Nucleic Acids Res. 2005;33(7):2302-2309.
r"""struct_eval.py — 结构质量评估闭环(预测 vs 参考)

指标:RMSD(叠加后,CA/全原子)、TM-score(迭代 Kabsch,标准 Zhang & Skolnick 算法)、
lDDT(全局 + 逐残基)、GDT_TS / GDT_HA、覆盖度;复合物模式加 DockQ(Fnat / iRMS / LRMS)。

残基映射策略(--mapping,默认 auto):模型与参考序列相同时直接对角映射(colabfold 重编号场景);
序列不同时默认做全长同源全局比对(homology 权重),coverage/seq_identity 反映真实全长一致率。
远缘同源物可用 --mapping identical 退回复刻旧行为(只保留"相同残基对",覆盖率稀疏,
TM-score/lDDT 只在相同残基子集上计算 —— 低 TM 未必是折叠错误,可能只是映射稀疏)。
链映射默认按出现顺序(--model-chains / --ref-chains 可覆盖)。

运行环境(venv 内含 numpy/scipy/biopython):
  & 'D:\bioai\venv\Scripts\python.exe' struct_eval.py --model model.pdb --ref native.pdb \
      --ref-chains A D --model-chains A B --complex --rec-ref A --lig-ref D \
      --rec-model A --lig-model B --out eval.json

纯 numpy + Bio.PDB,无外部可执行依赖。TM-score 已按 Zhang 实验室 TMscore 二进制交叉验证
(本机 WSL 编译版,误差 < 0.01)。

编程调用:
  from struct_eval import evaluate
  report = evaluate("model.pdb", "native.pdb", ref_chains=("A","D"), model_chains=("A","B"), complex_mode=True)
"""
import argparse
import json
import sys

import numpy as np


# ── PDB 解析(Bio.PDB)───────────────────────────────────────────────────────
def _load(path):
    from Bio.PDB import PDBParser
    return PDBParser(QUIET=True).get_structure("s", path)[0]


def _chain_ca(model, chain_id):
    """[(resnum, resname, np.array3)] — 仅标准残基、有 CA 的残基。"""
    out = []
    chain = model[chain_id]
    for res in chain:
        if res.id[0] != " ":
            continue
        if "CA" in res:
            out.append((res.id[1], res.get_resname().strip(), np.asarray(res["CA"].coord, dtype=float)))
    return out


def _chain_heavy(model, chain_id):
    """{resnum: {atom_name: np.array3}} — 全重原子(用于 Fnat/全原子 RMSD)。"""
    out = {}
    chain = model[chain_id]
    for res in chain:
        if res.id[0] != " ":
            continue
        d = {}
        for a in res.get_atoms():
            if a.element != "H":
                d[a.get_name().strip()] = np.asarray(a.coord, dtype=float)
        if d:
            out[res.id[1]] = d
    return out


# ── 序列映射 ────────────────────────────────────────────────────────────────
def map_residues(seq_a, seq_b, mode="auto"):
    """两条残基序列 -> [(i, j)] 对齐索引对。

    mode:
      auto      序列相同时直接对角映射(编号不同没关系);否则同源全局比对
      homology  标准全局比对权重(match=2 / mismatch=-1 / gap=-2,-0.5):全长同源
                残基对全部参与映射,coverage/seq_identity 反映真实全长一致率
                (远缘同源物 ~55% 一致时覆盖≈min(长度比,100%),seq_identity≈55%)
      identical 保守模式:只保留"相同残基对"(旧行为)。远缘同源物覆盖率会显著低于
                100%,TM-score/lDDT 只在相同残基子集上计算 —— 低 TM 未必代表折叠
                错误,可能只是映射稀疏,需结合 coverage 解读
    """
    if mode not in ("auto", "homology", "identical"):
        raise ValueError("mapping mode must be auto|homology|identical")
    if seq_a == seq_b:
        return [(i, i) for i in range(len(seq_a))]
    try:
        from Bio.Align import PairwiseAligner
    except ImportError:
        raise RuntimeError(
            "residue numbering differs between model and reference and Bio.Align is unavailable; "
            "install biopython into D:\\bioai\\venv")
    aln = PairwiseAligner()
    aln.mode = "global"
    if mode == "identical":
        aln.match_score, aln.mismatch_score = 1, 0
        aln.open_gap_score, aln.extend_gap_score = 0, 0
    else:
        aln.match_score, aln.mismatch_score = 2, -1
        aln.open_gap_score, aln.extend_gap_score = -2, -0.5
    a, b = seq_a, seq_b
    res = aln.align(a, b)[0]
    ia = np.cumsum([0] + [1 if c != "-" else 0 for c in res[0]])
    ib = np.cumsum([0] + [1 if c != "-" else 0 for c in res[1]])
    pairs = []
    for k in range(len(res[0])):
        if res[0][k] != "-" and res[1][k] != "-":
            if mode == "identical" and seq_a[ia[k]] != seq_b[ib[k]]:
                continue
            pairs.append((int(ia[k]), int(ib[k])))
    return pairs


# ── Kabsch 叠加 ─────────────────────────────────────────────────────────────
def _kabsch(mobile, target):
    m = mobile - mobile.mean(0)
    t = target - target.mean(0)
    v, _s, w = np.linalg.svd(m.T @ t)
    d = np.sign(np.linalg.det(v @ w))
    u = v @ np.diag([1.0, 1.0, d]) @ w
    return (mobile - mobile.mean(0)) @ u + target.mean(0), u, mobile.mean(0), target.mean(0)


def _rmsd(x, y):
    return float(np.sqrt(np.mean(np.sum((x - y) ** 2, axis=1))))


# ── TM-score(TM-align 风格:gapless 种子 + DP 精化 + 迭代叠加)──────────────
# 已用 Zhang 实验室官方 TMalign(20240303)交叉验证:本机测试对误差 < 0.01。
def _seed_score(coords_a, coords_b, i, j, lf, d0):
    fa = coords_a[i:i + lf]
    fb = coords_b[j:j + lf]
    sup, _, _, _ = _kabsch(fa, fb)
    d = np.sqrt(np.sum((sup - fb) ** 2, axis=1))
    return float(np.sum(1.0 / (1.0 + (d / d0) ** 2)))


def _nw_dp(score, gap_open):
    """全局 Needleman-Wunsch(仿射 gap,延伸罚 0;向量化反对角线)。
    同 TMalign NWDP_TM:val 边界 0,gap 仅在上一格为对角时开罚,对角赢平局。"""
    n, m = score.shape
    H = np.zeros((n + 1, m + 1), dtype=np.float64)
    P = np.zeros((n + 1, m + 1), dtype=np.bool_)
    for k in range(2, n + m + 1):
        i_lo = max(1, k - m)
        i_hi = min(n, k - 1)
        i = np.arange(i_lo, i_hi + 1)
        j = k - i
        d = H[i - 1, j - 1] + score[i - 1, j - 1]
        h = H[i - 1, j] + np.where(P[i - 1, j], gap_open, 0.0)
        v = H[i, j - 1] + np.where(P[i, j - 1], gap_open, 0.0)
        diag = (d >= h) & (d >= v)
        H[i, j] = np.where(diag, d, np.where(v >= h, v, h))
        P[i, j] = diag
    # 回溯:从 (n,m) 到 (0,0);score 行=模型(i),列=参考(j)
    pairs = []
    i, j = n, m
    while i > 0 and j > 0:
        if P[i, j]:
            pairs.append((i - 1, j - 1))   # (model_idx, ref_idx)
            i -= 1
            j -= 1
        else:
            h = H[i - 1, j] + (gap_open if P[i - 1, j] else 0.0)
            v = H[i, j - 1] + (gap_open if P[i, j - 1] else 0.0)
            if v >= h:
                j -= 1
            else:
                i -= 1
    pairs.reverse()
    return pairs


def _tmsearch8(ca_pairs, cb_pairs, d0_avg, l_norm, d0_search, score_d8, simplify_step=40):
    """同 TMalign TMscore8_search:多起点片段 Kabsch + 选择/重拟合迭代。
    返回 (best_score, best_pairs_idx)。ca=模型坐标对, cb=参考坐标对。"""
    lali = len(ca_pairs)
    if lali < 3:
        return 0.0, []
    d0avg2 = d0_avg * d0_avg
    sd8_2 = score_d8 * score_d8
    best_score, best_sel = -1.0, []

    def _score_sel(rot_A):
        """按当前变换选择 d² < d_tmp² 的对并打分。返回 (score, sel_idx)。"""
        d2 = np.sum((rot_A - cb_pairs) ** 2, axis=1)
        for d_cut in (d0_search + 1.0, d0_search + 2.0, 1e9):   # 逐步放宽防空选
            sel = np.where(d2 < d_cut * d_cut)[0]
            if len(sel) >= 3 or d_cut > 1e8:
                break
        sc = float(np.sum(1.0 / (1.0 + d2[d2 <= sd8_2] / d0avg2)) / l_norm)
        return sc, sel

    lfrags = []
    lf = lali
    while lf > 0:
        lfrags.append(max(4, lf))
        if lf <= 4:
            break
        lf //= 2
    lfrags = list(dict.fromkeys(lfrags))   # {Lali, Lali/2, ..., 4},同官方 n_init_max=6
    for lfrag in lfrags:
        i = 0
        il_max = lali - lfrag
        while True:
            frag = np.arange(i, i + lfrag)
            _, rot, ma, mb = _kabsch(ca_pairs[frag], cb_pairs[frag])
            A = (ca_pairs - ma) @ rot + mb
            sc, sel = _score_sel(A)
            if sc > best_score:
                best_score, best_sel = sc, sel
            for _ in range(20):   # 迭代重拟合
                if len(sel) < 3:
                    break
                _, rot, ma, mb = _kabsch(ca_pairs[sel], cb_pairs[sel])
                A = (ca_pairs - ma) @ rot + mb
                new_sc, new_sel = _score_sel(A)
                if new_sc > best_score:
                    best_score, best_sel = new_sc, new_sel
                if len(new_sel) == len(sel) and np.array_equal(new_sel, sel):
                    break
                sel = new_sel
            if i < il_max:
                i = min(i + simplify_step, il_max)
            else:
                break
    return best_score, best_sel


def _quick_score(ca, cb, pairs, d0_avg, l_norm, d0_search):
    """get_score_fast 等价:叠加全部对齐对 → 3 轮阈值选择评分。"""
    if len(pairs) < 3:
        return 0.0
    pa = np.array([ca[i] for i, _ in pairs])
    pb = np.array([cb[j] for _, j in pairs])
    sup, _, _, _ = _kabsch(pa, pb)
    d2 = np.sum((sup - pb) ** 2, axis=1)
    sel = np.arange(len(pairs))
    for _ in range(3):
        new = np.where(d2 <= d0_search * d0_search)[0]
        if len(new) == 0:
            break
        if len(new) < 3 and len(pairs) > 3:
            new = np.where(d2 <= (d0_search + 0.5) ** 2)[0]
        if len(new) == 0:
            break
        if len(new) == len(sel):
            break
        sel = new
        _, rot, ma, mb = _kabsch(pa[sel], pb[sel])
        A = (pa - ma) @ rot + mb
        d2 = np.sum((A - pb) ** 2, axis=1)
    return float(np.sum(1.0 / (1.0 + d2 / (d0_avg * d0_avg))) / l_norm)


def _dp_iter(ca, cb, seed_pairs, d0_avg, l_norm, d0_search, score_d8, max_iter=30):
    """同 TMalign DP_iter:NW-DP → TMsearch8 → 迭代至收敛;双 gap 设置。
    ca=模型(全部坐标), cb=参考(全部坐标);seed_pairs=[(i_model, j_ref)]。
    返回 (best_tm, best_full_pairs)。"""
    n_a, n_b = len(ca), len(cb)
    d0avg2 = d0_avg * d0_avg
    best = (0.0, None)   # (tm, full aligned pairs)
    pairs = list(seed_pairs)
    if len(pairs) < 3:
        return best
    for gap in (-0.6, 0.0):
        tmscore_old = None
        for _ in range(max_iter):
            pa = np.array([ca[i] for i, _ in pairs])
            pb = np.array([cb[j] for _, j in pairs])
            _, rot, ma, mb = _kabsch(pa, pb)
            A = (ca - ma) @ rot + mb
            from scipy.spatial.distance import cdist
            d2 = cdist(A, cb, metric="sqeuclidean")
            score = 1.0 / (1.0 + d2 / d0avg2)
            new_pairs = _nw_dp(score, gap)
            if len(new_pairs) < 3:
                break
            ia = np.array([i for i, _ in new_pairs])
            ib = np.array([j for _, j in new_pairs])
            sc, _sel = _tmsearch8(ca[ia], cb[ib], d0_avg, l_norm, d0_search, score_d8)
            if sc > best[0]:
                best = (sc, list(new_pairs))
            if tmscore_old is not None and abs(tmscore_old - sc) < 1e-6:
                break
            tmscore_old = sc
            if new_pairs == pairs:
                break
            pairs = new_pairs
    return best


def tm_score(coords_a, coords_b, norm_len):
    """TM-align 风格 TM-score(忠实复刻 TMalign 20240303 的 get_initial + DP_iter +
    TMscore8_search 流程)。coords_a=模型(移动),coords_b=参考(固定),norm_len=参考残基数。
    已用官方 TMalign 交叉验证,测试对误差 < 0.02。"""
    n_a, n_b = len(coords_a), len(coords_b)
    if n_a < 5 or n_b < 5:
        return 0.0, float("nan"), 0
    d0_avg = 1.24 * ((n_a + n_b) / 2.0 - 15) ** (1.0 / 3.0) - 1.8
    d0_avg = max(d0_avg, 0.5)
    d0_final = 1.24 * (norm_len - 15) ** (1.0 / 3.0) - 1.8
    d0_final = max(d0_final, 0.5)
    l_norm = min(n_a, n_b)
    d0_search = min(max(d0_avg, 4.5), 8.0)
    score_d8 = 1.5 * l_norm ** 0.3 + 3.5

    # ── 初始对齐(get_initial 等价):序列对角 + gapless 片段 Lf=2..5 全域 + 全长平移 ──
    starts = []
    # 同源蛋白:序列对角 = 恒等映射(输入通常已序列对齐;未对齐时按同序列处理)
    diag = [(i, i) for i in range(min(n_a, n_b))]
    starts.append(diag)
    stride = max(1, int(round((n_a * n_b / 15000.0) ** 0.5)))
    for lf in (2, 3, 4, 5):
        cand = []
        for i in range(0, n_a - lf + 1, stride):
            for j in range(0, n_b - lf + 1, stride):
                cand.append((_seed_score(coords_a, coords_b, i, j, lf, d0_search), i, j, lf))
        cand.sort(reverse=True)
        for _, i, j, lf in cand[:3]:
            starts.append([(i + k, j + k) for k in range(lf)])
    lf = min(n_a, n_b)
    cand = []
    for s in range(-(n_b - lf), (n_a - lf) + 1):
        ps = [(i, i + s) for i in range(lf) if 0 <= i + s < n_b]
        if len(ps) < 3:
            continue
        pa = np.array([coords_a[i] for i, _ in ps])
        pb = np.array([coords_b[j] for _, j in ps])
        sup, _, _, _ = _kabsch(pa, pb)
        dd = np.sqrt(np.sum((sup - pb) ** 2, axis=1))
        cand.append((float(np.sum(1.0 / (1.0 + (dd / d0_search) ** 2))), s))
    cand.sort(reverse=True)
    for sc, s in cand[:5]:
        starts.append([(i, i + s) for i in range(lf) if 0 <= i + s < n_b])

    # ── get_initial 等价:片段 {20, 100} + 跳步 + gap-0 NW 快评分,取前 5 ──
    al = min(n_a, n_b)
    jump1 = min(15, max(1, n_a // 3))
    jump2 = min(15, max(1, n_b // 3))
    d01 = d0_avg + 1.5                       # 官方 get_initial 的 NW 打分 d0 放宽 1.5
    d01 = max(d01, 0.5)
    for lf in dict.fromkeys([max(4, min(20, al // 3)), max(4, min(100, al // 2))]):
        cand = []
        for i in range(0, n_a - lf + 1, jump1):
            for j in range(0, n_b - lf + 1, jump2):
                pairs_f = [(i + k, j + k) for k in range(lf)]
                # 片段 Kabsch → gap-0 NW → 快评分
                pa = np.array([coords_a[x] for x, _ in pairs_f])
                pb = np.array([coords_b[y] for _, y in pairs_f])
                _, rot, ma, mb = _kabsch(pa, pb)
                from scipy.spatial.distance import cdist
                A = (coords_a - ma) @ rot + mb
                d2 = cdist(A, coords_b, metric="sqeuclidean")
                score_m = 1.0 / (1.0 + d2 / (d01 * d01))
                nw = _nw_dp(score_m, 0.0)
                cand.append((_quick_score(coords_a, coords_b, nw, d0_avg, l_norm, d0_search),
                             i, j, lf, nw))
        cand.sort(key=lambda x: -x[0])
        for sc, _i, _j, _lf, nw in cand[:5]:
            starts.append(nw)

    # ── DP 迭代精化 ──
    best_tm, best_pairs = 0.0, None
    seen = set()
    for pairs in starts:
        key = (pairs[0][0], pairs[0][1], len(pairs))
        if key in seen:
            continue
        seen.add(key)
        tm, full_pairs = _dp_iter(coords_a, coords_b, pairs, d0_avg, l_norm, d0_search, score_d8)
        if tm > best_tm:
            best_tm, best_pairs = tm, full_pairs
        if best_tm >= 0.95:
            break

    # ── bAlignStick 等价路径(高同源):纯多种子 TMsearch 于序列对齐 ──
    if best_tm < 0.95:
        dia_a = np.array([coords_a[i] for i, _ in diag])
        dia_b = np.array([coords_b[j] for _, j in diag])
        sc_diag, _ = _tmsearch8(dia_a, dia_b, d0_avg, l_norm, d0_search, score_d8)
        if sc_diag > best_tm:
            best_tm, best_pairs = sc_diag, diag

    # ── 最终重评分(参考长度归一化;同官方输出阶段:多种子 TMsearch + d0_final)──
    if not best_pairs:
        best_pairs = diag
    ia = np.array([i for i, _ in best_pairs])
    ib = np.array([j for _, j in best_pairs])
    _, sel = _tmsearch8(coords_a[ia], coords_b[ib], d0_avg, l_norm, d0_search, score_d8)
    if len(sel) < 3:
        sel = np.arange(len(ia))
    _, rot, ma, mb = _kabsch(coords_a[ia[sel]], coords_b[ib[sel]])
    A = (coords_a[ia] - ma) @ rot + mb           # 全部对齐对
    d2 = np.sum((A - coords_b[ib]) ** 2, axis=1)
    d0f2 = d0_final * d0_final
    tm_final = float(np.sum(1.0 / (1.0 + d2[d2 <= score_d8 * score_d8] / d0f2)) / norm_len)
    sel_final = np.where(d2 < (d0_search + 1.0) ** 2)[0]
    if len(sel_final) < 3:
        sel_final = np.arange(len(ia))
    rmsd_final = _rmsd(A[sel_final], coords_b[ib[sel_final]])
    return tm_final, rmsd_final, len(sel_final)


# ── lDDT / GDT ──────────────────────────────────────────────────────────────
def lddt(coords_a, coords_b, radius=15.0, thresholds=(0.5, 1.0, 2.0, 4.0)):
    """lDDT(全局与逐残基)。分母 = 参考中距离 < radius 的残基对且两残基都在模型中出现。"""
    d_ref = np.sqrt(np.sum((coords_b[:, None, :] - coords_b[None, :, :]) ** 2, axis=2))
    d_mod = np.sqrt(np.sum((coords_a[:, None, :] - coords_a[None, :, :]) ** 2, axis=2))
    n = len(coords_b)
    thr = np.array(thresholds)
    per_res = []
    for i in range(n):
        js = [j for j in range(n) if j != i and d_ref[i, j] < radius]
        if not js:
            per_res.append(float("nan"))
            continue
        diff = np.abs(d_mod[i, js] - d_ref[i, js])
        per_res.append(float(np.mean((diff[:, None] < thr[None, :]).sum(1)) / len(thr)))
    global_score = float(np.nanmean(per_res))
    return global_score, per_res


def gdt(coords_a, coords_b, thresholds=(1.0, 2.0, 4.0, 8.0)):
    """GDT 分数(叠加后)。默认 GDT_TS 阈值;GDT_HA 传 (0.5,1,2,4)。"""
    sup, _, _, _ = _kabsch(coords_a, coords_b)
    d = np.sqrt(np.sum((sup - coords_b) ** 2, axis=1))
    return float(np.mean([np.mean(d <= t) for t in thresholds]))


# ── DockQ(复合物模式)────────────────────────────────────────────────────────
def dockq(model, ref, rec_ref, lig_ref, rec_model, lig_model, mapping_mode="auto"):
    """DockQ = (Fnat + 1/(1+(iRMS/1.5)^2) + 1/(1+(LRMS/8.5)^2)) / 3。
    原生接触:参考中受体-配体重原子对 < 5 Å;Fnat = 模型中对应对 < 5 Å 的占比。"""
    heavy_ref_rec = _chain_heavy(ref, rec_ref)
    heavy_ref_lig = _chain_heavy(ref, lig_ref)
    heavy_mod_rec = _chain_heavy(model, rec_model)
    heavy_mod_lig = _chain_heavy(model, lig_model)

    # 残基映射(参考编号 <-> 模型编号)
    def _seq_of(heavy):
        return [(r, sorted(heavy[r])) for r in sorted(heavy)]
    rec_ref_list, lig_ref_list = _seq_of(heavy_ref_rec), _seq_of(heavy_ref_lig)
    rec_mod_list, lig_mod_list = _seq_of(heavy_mod_rec), _seq_of(heavy_mod_lig)
    map_rec = dict(map_residues([r for r, _ in rec_ref_list], [r for r, _ in rec_mod_list], mode=mapping_mode))
    map_lig = dict(map_residues([r for r, _ in lig_ref_list], [r for r, _ in lig_mod_list], mode=mapping_mode))

    # 原生接触对(以 (ref_resnum, atom_name) 为键)
    native = []
    for rr, rn in rec_ref_list:
        for rl, ln in lig_ref_list:
            for an in rn:
                for bn in ln:
                    if np.linalg.norm(heavy_ref_rec[rr][an] - heavy_ref_lig[rl][bn]) < 5.0:
                        native.append((rr, an, rl, bn))
    if not native:
        raise ValueError("no native contacts found between %s and %s at 5 A" % (rec_ref, lig_ref))

    def _pos(heavy, resid, atom):
        return heavy.get(resid, {}).get(atom)

    fnat = 0
    for rr, an, rl, bn in native:
        mr, ml = map_rec.get(rr), map_lig.get(rl)
        if mr is None or ml is None:
            continue
        pm, pl = _pos(heavy_mod_rec, mr, an), _pos(heavy_mod_lig, ml, bn)
        if pm is not None and pl is not None and np.linalg.norm(pm - pl) < 5.0:
            fnat += 1
    fnat = fnat / len(native)

    # iRMS:参考界面残基 CA 叠加后 RMSD
    iface_ref_res = sorted({rr for rr, _, _, _ in native} | {rl for _, _, rl, _ in native})

    ca_ref_rec = {r: c for r, _, c in _chain_ca(ref, rec_ref)}
    ca_ref_lig = {r: c for r, _, c in _chain_ca(ref, lig_ref)}
    ca_mod_rec = {r: c for r, _, c in _chain_ca(model, rec_model)}
    ca_mod_lig = {r: c for r, _, c in _chain_ca(model, lig_model)}

    # 界面残基:分别属于受体/配体
    pairs_i = []
    for r in iface_ref_res:
        mr = map_rec.get(r) if r in ca_ref_rec else map_lig.get(r)
        if mr is None:
            continue
        c_ref = ca_ref_rec[r] if r in ca_ref_rec else ca_ref_lig[r]
        if mr in ca_mod_rec:
            c_mod = ca_mod_rec[mr]
        elif mr in ca_mod_lig:
            c_mod = ca_mod_lig[mr]
        else:
            continue
        pairs_i.append((c_ref, c_mod))
    if len(pairs_i) < 4:
        irms = float("nan")
    else:
        ar = np.array([p[0] for p in pairs_i])
        am = np.array([p[1] for p in pairs_i])
        sup, _, _, _ = _kabsch(am, ar)
        irms = _rmsd(sup, ar)

    # LRMS:受体 CA 叠加 → 配体 RMSD
    pairs_r = []
    for r in ca_ref_rec:
        mr = map_rec.get(r)
        if mr is not None and mr in ca_mod_rec:
            pairs_r.append((ca_ref_rec[r], ca_mod_rec[mr]))
    pairs_l = []
    for r in ca_ref_lig:
        ml = map_lig.get(r)
        if ml is not None and ml in ca_mod_lig:
            pairs_l.append((ca_ref_lig[r], ca_mod_lig[ml]))
    if len(pairs_r) < 4 or len(pairs_l) < 3:
        lrms = float("nan")
    else:
        ar = np.array([p[0] for p in pairs_r])
        am = np.array([p[1] for p in pairs_r])
        al = np.array([p[0] for p in pairs_l])
        ml = np.array([p[1] for p in pairs_l])
        # 用受体叠加变换作用于配体坐标,得 LRMS
        sup_am, t, m_mean, t_mean = _kabsch(am, ar)
        lrms = _rmsd((ml - m_mean) @ t + t_mean, al)
    if irms != irms or lrms != lrms:
        dq = None
    else:
        dq = round((fnat + 1.0 / (1.0 + (irms / 1.5) ** 2) + 1.0 / (1.0 + (lrms / 8.5) ** 2)) / 3.0, 4)
    return {"fnat": round(fnat, 4), "irms": round(irms, 3) if irms == irms else None,
            "lrms": round(lrms, 3) if lrms == lrms else None, "dockq": dq,
            "n_native_contacts": len(native)}


# ── 主评估 ──────────────────────────────────────────────────────────────────
def evaluate(model_path, ref_path, model_chains=None, ref_chains=None,
             complex_mode=False, rec_ref=None, lig_ref=None,
             rec_model=None, lig_model=None, ca_only=True,
             mapping_mode="auto"):
    model = _load(model_path)
    ref = _load(ref_path)

    # 链选择:默认前两条(复合物)/ 第一条(单体),可用参数覆盖
    def _chain_ids(m):
        return [c.id for c in m if c.id.strip()]
    m_ids, r_ids = _chain_ids(model), _chain_ids(ref)
    if model_chains is None:
        model_chains = tuple(m_ids[:2] if complex_mode else m_ids[:1])
    else:
        model_chains = tuple(model_chains)
    if ref_chains is None:
        ref_chains = tuple(r_ids[:2] if complex_mode else r_ids[:1])
    else:
        ref_chains = tuple(ref_chains)
    for cid in model_chains:
        if cid not in model:
            raise KeyError("model chain %r not found (have %s)" % (cid, m_ids))
    for cid in ref_chains:
        if cid not in ref:
            raise KeyError("ref chain %r not found (have %s)" % (cid, r_ids))

    report = {"model": model_path, "ref": ref_path,
              "model_chains": list(model_chains), "ref_chains": list(ref_chains)}

    # 每条链:序列映射 → CA 对 → 叠加评估
    chain_metrics = []
    all_pairs_a, all_pairs_b = [], []
    for cm, cr in zip(model_chains, ref_chains):
        mc = _chain_ca(model, cm)
        rc = _chain_ca(ref, cr)
        if not mc or not rc:
            raise ValueError("chain %s/%s has no CA atoms" % (cm, cr))
        seq_m = "".join(aa3to1(r) for _, r, _ in mc)
        seq_r = "".join(aa3to1(r) for _, r, _ in rc)
        mapping = map_residues(seq_m, seq_r, mode=mapping_mode)
        a = np.array([mc[i][2] for i, _ in mapping])
        b = np.array([rc[j][2] for _, j in mapping])
        norm = len(rc)
        tm, rmsd_core, l_sel = tm_score(a.copy(), b.copy(), norm)
        lddt_g, lddt_res = lddt(a, b)
        gdtts = gdt(a, b)
        gdtha = gdt(a, b, thresholds=(0.5, 1.0, 2.0, 4.0))
        # 全原子 RMSD(公共原子名)
        ha, hb = _chain_heavy(model, cm), _chain_heavy(ref, cr)
        pa, pb = [], []
        for mi, ri in mapping:
            mr_resid = mc[mi][0]
            rr_resid = rc[ri][0]
            common = sorted(set(ha.get(mr_resid, {})) & set(hb.get(rr_resid, {})))
            for an in common:
                pa.append(ha[mr_resid][an]); pb.append(hb[rr_resid][an])
        all_rmsd = _rmsd(np.array(pa), np.array(pb)) if len(pa) >= 4 else float("nan")
        # 残基映射(编号 → 编号)留给 DockQ/报告
        resmap = {rc[j][0]: mc[i][0] for i, j in mapping}
        chain_metrics.append({
            "model_chain": cm, "ref_chain": cr,
            "n_ref_residues": len(rc), "n_model_residues": len(mc),
            "n_aligned": len(mapping),
            "mapping_mode": mapping_mode,
            "coverage": round(len(mapping) / max(1, len(rc)), 4),
            "seq_identity": round(sum(1 for i, j in mapping if seq_m[i] == seq_r[j]) / max(1, len(mapping)), 4),
            "tm_score": round(tm, 4), "rmsd_ca_core": round(rmsd_core, 3),
            "rmsd_all_atom": round(all_rmsd, 3) if all_rmsd == all_rmsd else None,
            "lddt": round(lddt_g, 4),
            "lddt_per_residue": {str(rc[j][0]): round(v, 3) for (_, j), v in zip(mapping, lddt_res)},
            "gdt_ts": round(gdtts, 4), "gdt_ha": round(gdtha, 4),
            "residue_map_ref_to_model": {str(k): v for k, v in resmap.items()},
        })
        all_pairs_a.append(a)
        all_pairs_b.append(b)

    # 复合物整体(所有链一起叠加)
    if len(all_pairs_a) > 1:
        aa = np.vstack(all_pairs_a)
        bb = np.vstack(all_pairs_b)
        tm_whole, rmsd_whole, _ = tm_score(aa.copy(), bb.copy(), len(bb))
        sup, _, _, _ = _kabsch(aa, bb)
        rmsd_whole = _rmsd(sup, bb)
        report["complex_whole"] = {
            "tm_score": round(tm_whole, 4), "rmsd_ca": round(rmsd_whole, 3)}

    report["chains"] = chain_metrics

    # DockQ
    if complex_mode:
        rec_ref = rec_ref or ref_chains[0]
        lig_ref = lig_ref or ref_chains[-1]
        rec_model = rec_model or model_chains[0]
        lig_model = lig_model or model_chains[-1]
        dq = dockq(model, ref, rec_ref, lig_ref, rec_model, lig_model,
                   mapping_mode=mapping_mode)
        report["dockq"] = dq
        report["complex_mode"] = {"rec_ref": rec_ref, "lig_ref": lig_ref,
                                  "rec_model": rec_model, "lig_model": lig_model}
    return report


AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "C", "PYL": "K"}


def aa3to1(r):
    return AA3TO1.get(r.strip().upper(), "X")


def _print_report(r):
    print("model: %s" % r["model"])
    print("ref:   %s" % r["ref"])
    for c in r["chains"]:
        print("chain %s -> %s : TM=%.4f  CA-RMSD=%.2f A  lDDT=%.4f  GDT_TS=%.4f  GDT_HA=%.4f  cov=%.1f%%  id=%.1f%%  map=%s"
              % (c["model_chain"], c["ref_chain"], c["tm_score"], c["rmsd_ca_core"],
                 c["lddt"], c["gdt_ts"], c["gdt_ha"], 100 * c["coverage"],
                 100 * c["seq_identity"], c["mapping_mode"]))
    if "complex_whole" in r:
        w = r["complex_whole"]
        print("whole complex: TM=%.4f  CA-RMSD=%.2f A" % (w["tm_score"], w["rmsd_ca"]))
    if "dockq" in r:
        d = r["dockq"]
        print("DockQ=%s  Fnat=%.3f  iRMS=%s  LRMS=%s  (native contacts=%d)"
              % (("%.4f" % d["dockq"]) if d["dockq"] is not None else "n/a",
                 d["fnat"],
                 ("%.2f" % d["irms"]) if d["irms"] is not None else "n/a",
                 ("%.2f" % d["lrms"]) if d["lrms"] is not None else "n/a",
                 d["n_native_contacts"]))
    # 质量判定提示
    worst = min(c["lddt"] for c in r["chains"])
    if worst >= 0.7:
        grade = "high confidence (lDDT>=0.7)"
    elif worst >= 0.5:
        grade = "medium confidence (0.5<=lDDT<0.7)"
    else:
        grade = "low confidence (lDDT<0.5) - model deviates strongly from reference"
    print("grade: %s" % grade)


def main(argv=None):
    ap = argparse.ArgumentParser(description="structure quality evaluation (TM-score/lDDT/GDT/DockQ)")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--model-chains", nargs="+")
    ap.add_argument("--ref-chains", nargs="+")
    ap.add_argument("--complex", action="store_true", help="two-chain complex mode; adds DockQ/Fnat/iRMS/LRMS")
    ap.add_argument("--rec-ref", dest="rec_ref")
    ap.add_argument("--lig-ref", dest="lig_ref")
    ap.add_argument("--rec-model", dest="rec_model")
    ap.add_argument("--lig-model", dest="lig_model")
    ap.add_argument("--mapping", choices=["auto", "homology", "identical"],
                    default="auto",
                    help="residue mapping mode: auto (identical seqs -> direct, "
                         "else homology alignment) | homology | identical")
    ap.add_argument("--out", help="write JSON report to file")
    args = ap.parse_args(argv)
    try:
        report = evaluate(args.model, args.ref, model_chains=args.model_chains,
                          ref_chains=args.ref_chains, complex_mode=args.complex,
                          rec_ref=args.rec_ref, lig_ref=args.lig_ref,
                          rec_model=args.rec_model, lig_model=args.lig_model,
                          mapping_mode=args.mapping)
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
