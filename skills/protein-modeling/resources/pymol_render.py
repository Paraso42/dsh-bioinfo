#!/usr/bin/env python3
r"""pymol_render.py — PyMOL 无头科研级渲染(PDB → PNG)

运行环境(D:\bioai\venv:PyMOL 3.1.0 open-source,cgohlke cp313 wheel):
  & 'D:\bioai\venv\Scripts\python.exe' pymol_render.py complex.pdb --out render.png
  & 'D:\bioai\venv\Scripts\python.exe' pymol_render.py complex.pdb --out pub.png --style publication `
      --chains A+B --hetatm sticks --surface B --width 2400 --height 1800
  & 'D:\bioai\venv\Scripts\python.exe' pymol_render.py model.pdb --style rainbow --width 1600 --height 1200

输出约定:
  PNG(默认 300 dpi 光栅 + ray-trace);控制台打印 out/atoms 摘要。
  --style cartoon     默认:二级结构着色(螺旋红/折叠黄/环绿)
  --style publication 按链着色 + 配体按元素着色(可叠加 --surface)
  --style rainbow     B-factor 谱着色
  --style surface     溶剂可及面
  --style line        线框
  --chains A+B        只显示指定链(逗号或 + 分隔)
  --hetatm sticks|lines|off  非聚合物(配体/离子/水)显示方式,默认 sticks(按元素着色)
  --surface SEL       额外叠加半透明表面(如 A+B 或 chain B)
  --width/--height    光栅尺寸(默认 1600x1200;ray 分辨率)
  --bg white|black|transparent
  --no-ray            跳过 ray-trace(快速预览)
"""
import argparse
import os
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description="PyMOL headless publication renderer")
    ap.add_argument("pdb", help="input PDB/mmCIF path")
    ap.add_argument("--out", help="output PNG (default: <pdb stem>.png)")
    ap.add_argument("--style", default="cartoon",
                    choices=["cartoon", "publication", "rainbow", "surface", "line"])
    ap.add_argument("--chains", help="chain selection, e.g. A+B or A,B")
    ap.add_argument("--hetatm", default="sticks", choices=["sticks", "lines", "off"])
    ap.add_argument("--surface", help="selection to render as semi-transparent surface")
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=1200)
    ap.add_argument("--bg", default="white", choices=["white", "black", "transparent"])
    ap.add_argument("--no-ray", action="store_true")
    args = ap.parse_args(argv)

    if not os.path.exists(args.pdb):
        print("ERROR: input not found: %s" % args.pdb, file=sys.stderr)
        return 2
    out = args.out or (os.path.splitext(args.pdb)[0] + ".png")

    import pymol
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    sel = "all"
    if args.chains:
        chains = [c.strip() for c in args.chains.replace("+", ",").split(",") if c.strip()]
        sel = " or ".join("chain %s" % c for c in chains)

    cmd.load(args.pdb, "obj")
    cmd.hide("everything", "all")

    if args.style == "surface":
        cmd.show("surface", sel)
        cmd.color("gray80", sel)
    elif args.style == "line":
        cmd.show("lines", sel)
    elif args.style == "rainbow":
        cmd.show("cartoon", sel)
        cmd.spectrum("b", "blue_white_red", sel)
    elif args.style == "publication":
        cmd.show("cartoon", sel)
        cmd.util.cbc(sel)                       # 按链着色
    else:                                       # cartoon
        cmd.show("cartoon", sel)
        cmd.util.cbss(sel, "red", "yellow", "green")   # 螺旋/折叠/环

    het = "(%s) and hetatm" % sel
    if args.hetatm in ("sticks", "lines") and cmd.count_atoms(het) > 0:
        cmd.show(args.hetatm, het)
        cmd.util.cnc(het)          # PyMOL 3.x 移除 byelement;cnc = 碳绿 + 杂原子元素色
        if args.hetatm == "sticks":
            cmd.set("stick_radius", 0.18, het)
    else:
        cmd.hide("everything", het)

    if args.surface:
        cmd.create("surf_obj", args.surface)
        cmd.show("surface", "surf_obj")
        cmd.color("gray70", "surf_obj")
        cmd.set("transparency", 0.4, "surf_obj")

    cmd.set("ray_opaque_background", 1 if args.bg != "transparent" else 0)
    cmd.bg_color(args.bg if args.bg != "transparent" else "white")
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.orient(sel)

    if args.no_ray:
        cmd.set("ray_trace_mode", 0)
        cmd.draw(args.width, args.height)
    else:
        cmd.ray(args.width, args.height)
    cmd.png(out, dpi=300)

    n_atoms = cmd.count_atoms(sel)
    print("rendered %s (%d atoms, style=%s) -> %s" % (args.pdb, n_atoms, args.style, out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
