// protein-tools.js — protein structure prediction & interaction tools for
// DeepSeek Harness (DSH) agent presets.
// Copyright (c) 2026 dsh-bioinfo contributors. SPDX-License-Identifier: MIT.
//
// Registers 7 model tools:
//   esmfold_predict / pp_interact / vina_dock / af2_predict / struct_eval /
//   vscreen_run / md_run
// Every tool executes a configurable backend script through ctx.shell — this
// file is the frontend; the per-tool backend contract (argv, output JSON,
// dependencies, licenses) is documented in BACKENDS.md. When a backend is
// missing, the shell error is surfaced to the model instead of crashing.
// The plugin only consumes host services (shell/fs/tools) and publishes none,
// so the preset row needs no isolate realm.
//
// Tool registration is dual-path:
//   - Sandboxed load (a `harness` global exists): ctx.tools.register accepts
//     only tools returned by harness.defineTool(...) (DYNAMIC_TOOL marker).
//   - Preset file-plugin load (raw ctx, no `harness`): ctx.tools.register
//     stores the definition verbatim and the model layer projects
//     definition.parameters straight to the API. parameters must therefore be
//     the compiled root-object form { type:'object', properties, required }
//     (equivalent to dsh-tools parameterSchemaSpecToJsonSchema output); a flat
//     per-property map makes the API reject the tool with
//     "Invalid schema ... got 'type: null'".
// defineToolDef() produces the same compiled shape on both paths.

module.exports = {
  name: 'protein-tools',
  inject: ['tools', 'shell'],
  apply(ctx) {
    const sandboxPolicy = ctx.get('sandboxPolicy')
    const shellEnv = ctx.get('shellEnv')
    const fs = ctx.get('fs')

    // ── configuration ──────────────────────────────────────────────────────
    // Every machine-specific location can be overridden via environment
    // variables; the defaults reproduce the original D:\bioai deployment
    // layout, so an unconfigured machine keeps working unchanged.
    const env = (typeof process !== 'undefined' && process.env) || {}
    const cfg = {
      python: env.BIO_TOOLS_PYTHON || 'C:\\Program Files\\Python313\\python.exe',
      venvPy: env.BIO_TOOLS_VENV_PY || 'D:\\bioai\\venv\\Scripts\\python.exe',
      resDir: env.BIO_TOOLS_RES_DIR || 'C:\\deepseek-harness\\.dsh\\.agent-presets\\bioinfo\\skills\\protein-modeling\\resources',
      resPqDir: env.BIO_TOOLS_RES_PQ_DIR || 'C:\\deepseek-harness\\.dsh\\.agent-presets\\bioinfo\\skills\\protein-quality\\resources',
      resCiDir: env.BIO_TOOLS_RES_CI_DIR || 'C:\\deepseek-harness\\.dsh\\.agent-presets\\bioinfo\\skills\\chem-informatics\\resources',
      jobsDir: env.BIO_TOOLS_JOBS_DIR || 'D:\\bioai\\jobs',
      biopython: env.BIO_TOOLS_BIOPYTHON || 'D:\\biopython',
    }

    const sq = (s) => String(s).replace(/'/g, "''")

    // Compile the author-facing flat parameter map into a root-object JSON
    // Schema (equivalent to dsh-tools' compiled output): per-property
    // required:true is lifted into the root `required` array.
    const compileParameters = (flat) => {
      const properties = {}
      const required = []
      for (const key of Object.keys(flat)) {
        const copy = { ...flat[key] }
        if (copy.required === true) {
          required.push(key)
          delete copy.required
        }
        properties[key] = copy
      }
      return required.length > 0
        ? { type: 'object', properties, required }
        : { type: 'object', properties }
    }

    const defineToolDef = (def) => {
      if (typeof harness !== 'undefined' && harness !== null
        && typeof harness.defineTool === 'function') {
        return harness.defineTool(def)
      }
      return {
        name: def.name,
        description: def.description,
        parameters: compileParameters(def.parameters),
        output: def.output,
        execute: def.execute,
      }
    }

    const runShell = async (command, exec, timeoutMs) => {
      const request = { command, timeoutMs }
      if (sandboxPolicy !== undefined) {
        request.sandboxPolicy = sandboxPolicy.resolve(
          exec !== undefined && exec.agent !== undefined ? { session: exec.agent.session } : {}
        )
      }
      if (shellEnv !== undefined && exec !== undefined) {
        try { request.dshEnv = shellEnv.collect(exec) } catch (e) { /* optional */ }
      }
      const spec = ctx.shell.resolve(request)
      return ctx.shell.run(spec)
    }

    const outOf = (result) => ({
      exitCode: result.exitCode,
      stdout: result.stdout ? result.stdout.text : '',
      stderr: result.stderr ? result.stderr.text : '',
    })

    const renderText = (_args, value) => {
      let body = value.stdout || ''
      if (value.stderr && value.stderr.length > 0) {
        if (body.length > 0 && !body.endsWith('\n')) body += '\n'
        body += '[stderr]\n' + value.stderr
      }
      if (value.report !== undefined) {
        if (body.length > 0 && !body.endsWith('\n')) body += '\n'
        body += '[report]\n' + JSON.stringify(value.report, null, 2)
      }
      if (body.length === 0) body = '(no output)'
      if (value.exitCode !== 0 && value.exitCode !== null && value.exitCode !== undefined) {
        body += '\n[exit code: ' + value.exitCode + ']'
      }
      return [{ type: 'text', text: body }]
    }

    const readJson = async (path) => {
      try {
        const target = await fs.resolve(path)
        const raw = await fs.readText(target)
        return JSON.parse(raw)
      } catch (e) {
        return null
      }
    }

    ctx.effect(() => {
      const disposers = []

      // ── 1. esmfold_predict: sequence → PDB via the free ESM Atlas API ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'esmfold_predict',
        description: 'Fold a protein sequence into a PDB structure via the free ESM Metagenomic Atlas API (no key). Fast (minutes) cloud channel for rough models; browser-UA + retry built in. For research-grade local prediction use af2_predict instead. Returns stdout/stderr and the written PDB path.',
        parameters: {
          sequence: { type: 'string', required: true, description: 'Protein sequence (letters only).' },
          out: { type: 'string', description: 'Output PDB path. Default: <jobs dir>\\esmfold\\esm_<timestamp>.pdb (BIO_TOOLS_JOBS_DIR).' },
          retries: { type: 'number', description: 'Retry count for transient network failures. Default 4.' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 300000.' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          const out = typeof args.out === 'string' && args.out.length > 0
            ? args.out
            : cfg.jobsDir + '\\esmfold\\esm_' + Date.now() + '.pdb'
          const retries = Number.isFinite(args.retries) ? ' --retries ' + Math.max(1, Math.floor(args.retries)) : ''
          const command = [
            "$env:PYTHONPATH='" + cfg.biopython + "'",
            "& '" + cfg.python + "' '" + cfg.resDir + '\\esmfold_api.py\' \'' + sq(args.sequence) + '\' --out \'' + sq(out) + '\'' + retries,
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 300000)
          const value = outOf(result)
          value.pdbPath = result.exitCode === 0 ? out : undefined
          return value
        },
      })))

      // ── 2. pp_interact: protein-protein interface analysis (Bio.PDB) ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'pp_interact',
        description: 'Analyze the interaction interface of a protein-protein complex PDB: atom contacts (NeighborSearch, 5 A default), interface residues, and buried surface area (BSA via ShrakeRupley delta-SASA). Runs Biopython 1.87 locally; returns the full JSON report.',
        parameters: {
          complex: { type: 'string', required: true, description: 'Path to the complex PDB file.' },
          chains: { type: 'array', items: { type: 'string' }, description: 'Two chain ids, e.g. ["A","B"]. Defaults to ["A","B"].' },
          cutoff: { type: 'number', description: 'Contact cutoff in Angstrom. Default 5.0.' },
          out: { type: 'string', description: 'Optional JSON report path; default <jobs dir>\\pp_interact\\interface_<timestamp>.json.' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 120000.' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          const jsonPath = typeof args.out === 'string' && args.out.length > 0
            ? args.out
            : cfg.jobsDir + '\\pp_interact\\interface_' + Date.now() + '.json'
          const chains = Array.isArray(args.chains) && args.chains.length === 2
            ? args.chains.map(sq).join("' '")
            : "A' 'B"
          const cutoff = Number.isFinite(args.cutoff) ? ' --cutoff ' + args.cutoff : ''
          const command = [
            "$env:PYTHONPATH='" + cfg.biopython + "'",
            "& '" + cfg.python + "' '" + cfg.resDir + "\\pp_interact.py' --complex '" + sq(args.complex) + "' --chains '" + chains + "'" + cutoff + " --out '" + sq(jsonPath) + "'",
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 120000)
          const value = outOf(result)
          const report = await readJson(jsonPath)
          value.report = report !== null ? report : undefined
          value.jsonPath = jsonPath
          return value
        },
      })))

      // ── 3. vina_dock: AutoDock Vina docking (protein–small molecule) ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'vina_dock',
        description: 'Dock a small-molecule ligand into a protein receptor with AutoDock Vina: meeko ligand preparation (SMILES or .sdf/.mol/.mol2), docking via the local vina binary, pose scoring/ranking and JSON report. Requires the Vina backend toolchain (RDKit + meeko in the venv interpreter, local vina binary) — see BACKENDS.md.',
        parameters: {
          receptorPdbqt: { type: 'string', description: 'Path to the prepared receptor PDBQT.' },
          receptor: { type: 'string', description: 'Receptor PDB path; auto-converted to rigid PDBQT (pdb_to_pdbqt.py). Alternative to receptorPdbqt.' },
          smiles: { type: 'string', description: 'Ligand SMILES string.' },
          ligand: { type: 'string', description: 'Ligand file (.sdf/.mol/.mol2).' },
          ligandPdbqt: { type: 'string', description: 'Pre-made ligand PDBQT (skip meeko prep).' },
          center: { type: 'array', items: { type: 'number' }, required: true, description: 'Search box center [x, y, z] in Angstrom.' },
          size: { type: 'array', items: { type: 'number' }, required: true, description: 'Search box size [x, y, z] in Angstrom.' },
          exhaustiveness: { type: 'number', description: 'Vina exhaustiveness. Default 16.' },
          outdir: { type: 'string', description: 'Output directory. Default <jobs dir>\\dock_<timestamp>.' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 3600000.' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          if (!Array.isArray(args.center) || args.center.length !== 3) throw new Error('center must be [x, y, z]')
          if (!Array.isArray(args.size) || args.size.length !== 3) throw new Error('size must be [x, y, z]')
          const outdir = typeof args.outdir === 'string' && args.outdir.length > 0
            ? args.outdir
            : cfg.jobsDir + '\\dock_' + Date.now()
          const jsonPath = outdir + '\\report.json'
          let recArg = ''
          if (typeof args.receptorPdbqt === 'string' && args.receptorPdbqt.length > 0) {
            recArg = ' --receptor-pdbqt \'' + sq(args.receptorPdbqt) + '\''
          } else if (typeof args.receptor === 'string' && args.receptor.length > 0) {
            recArg = ' --receptor \'' + sq(args.receptor) + '\''
          } else {
            throw new Error('provide receptorPdbqt or receptor (PDB path)')
          }
          let src = ''
          if (typeof args.smiles === 'string' && args.smiles.length > 0) src = ' --smiles \'' + sq(args.smiles) + '\''
          else if (typeof args.ligand === 'string' && args.ligand.length > 0) src = ' --ligand \'' + sq(args.ligand) + '\''
          else if (typeof args.ligandPdbqt === 'string' && args.ligandPdbqt.length > 0) src = ' --ligand-pdbqt \'' + sq(args.ligandPdbqt) + '\''
          else throw new Error('provide one of smiles / ligand / ligandPdbqt')
          const exh = Number.isFinite(args.exhaustiveness) ? ' --exhaustiveness ' + Math.floor(args.exhaustiveness) : ''
          const command = [
            "$env:PYTHONPATH='" + cfg.biopython + "'",
            "& '" + cfg.venvPy + "' '" + cfg.resDir + "\\vina_dock.py'" + recArg + src +
            ' --center ' + args.center.join(' ') + ' --size ' + args.size.join(' ') + exh +
            ' --outdir \'' + sq(outdir) + '\' --name dock --out \'' + sq(jsonPath) + '\'',
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 3600000)
          const value = outOf(result)
          const report = await readJson(jsonPath)
          value.report = report !== null ? report : undefined
          value.outdir = outdir
          return value
        },
      })))

      // ── 4. af2_predict: LocalColabFold (AF2 / AF2-Multimer) prediction ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'af2_predict',
        description: 'Predict protein structure (monomer or complex) with LocalColabFold (AlphaFold2 / AF2-Multimer), typically on a GPU inside WSL2. For a complex, join chains with ":" in one fasta record. Results land in outdir (unrelaxed_rank_001_*.pdb + PAE/pLDDT plots). Requires the colabfold backend — see BACKENDS.md.',
        parameters: {
          fasta: { type: 'string', required: true, description: 'Path to the input fasta file (multi-chain complex: join with ":").' },
          modelType: { type: 'string', description: 'alphafold2_multimer_v3 (complex, default) | alphafold2_ptm (monomer) | alphafold2 | alphafold2_multimer_v1/v2 | auto. Default alphafold2_multimer_v3.' },
          numModels: { type: 'number', description: 'Number of models. Default 1 (8GB VRAM discipline).' },
          numRecycle: { type: 'number', description: 'Recycle iterations. Default 3.' },
          msaMode: { type: 'string', description: 'MSA mode; default "mmseqs2_uniref_env" (MMseqs2 server, needs internet). Use "single_sequence" for offline/fast runs (lower accuracy).' },
          outdir: { type: 'string', description: 'Output directory. Default <jobs dir>\\af2_<timestamp>.' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 3600000 (1h); run in background for big jobs.' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          const outdir = typeof args.outdir === 'string' && args.outdir.length > 0
            ? args.outdir
            : cfg.jobsDir + '\\af2_' + Date.now()
          const mt = typeof args.modelType === 'string' && args.modelType.length > 0 ? args.modelType : 'alphafold2_multimer_v3'
          const nm = Number.isFinite(args.numModels) ? ' -NumModels ' + Math.max(1, Math.floor(args.numModels)) : ''
          const nr = Number.isFinite(args.numRecycle) ? ' -NumRecycle ' + Math.max(0, Math.floor(args.numRecycle)) : ''
          const ms = typeof args.msaMode === 'string' && args.msaMode.length > 0 ? ' -MsaMode \'' + sq(args.msaMode) + '\'' : ''
          const command = [
            "& '" + cfg.resDir + "\\run_colabfold.ps1' -Fasta '" + sq(args.fasta) + "' -OutDir '" + sq(outdir) + "' -ModelType " + mt + nm + nr + ms,
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 3600000)
          const value = outOf(result)
          value.outdir = outdir
          return value
        },
      })))

      // ── 5. struct_eval: quality vs reference (TM-score/lDDT/GDT/DockQ) ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'struct_eval',
        description: 'Evaluate a predicted structure against a reference (crystal/native): TM-score (validated against official TMalign), CA/all-atom RMSD, lDDT, GDT-TS/GDT-HA, and for complexes DockQ (Fnat/iRMS/LRMS). Residue numbering is mapped automatically by sequence alignment. Returns the full JSON report plus a confidence grade.',
        parameters: {
          model: { type: 'string', required: true, description: 'Path to the predicted/model PDB.' },
          ref: { type: 'string', required: true, description: 'Path to the reference/native PDB.' },
          complex: { type: 'boolean', description: 'Two-chain complex mode; adds DockQ/Fnat/iRMS/LRMS.' },
          modelChains: { type: 'array', items: { type: 'string' }, description: 'Model chain ids, e.g. ["A","B"]. Default: first chain(s).' },
          refChains: { type: 'array', items: { type: 'string' }, description: 'Reference chain ids, e.g. ["A","D"]. Default: first chain(s).' },
          recRef: { type: 'string', description: 'Complex mode: reference receptor chain (default first).' },
          ligRef: { type: 'string', description: 'Complex mode: reference ligand chain (default last).' },
          recModel: { type: 'string', description: 'Complex mode: model receptor chain (default first).' },
          ligModel: { type: 'string', description: 'Complex mode: model ligand chain (default last).' },
          out: { type: 'string', description: 'JSON report path. Default <jobs dir>\\struct_eval\\eval_<timestamp>.json' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 900000 (large proteins need ~1 min).' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          const jsonPath = typeof args.out === 'string' && args.out.length > 0
            ? args.out
            : cfg.jobsDir + '\\struct_eval\\eval_' + Date.now() + '.json'
          const mc = Array.isArray(args.modelChains) && args.modelChains.length > 0
            ? " --model-chains '" + args.modelChains.map(sq).join("' '") + "'" : ''
          const rc = Array.isArray(args.refChains) && args.refChains.length > 0
            ? " --ref-chains '" + args.refChains.map(sq).join("' '") + "'" : ''
          const dq = typeof args.recRef === 'string' || typeof args.ligRef === 'string'
            || typeof args.recModel === 'string' || typeof args.ligModel === 'string'
            ? (typeof args.recRef === 'string' ? " --rec-ref '" + sq(args.recRef) + "'" : '')
              + (typeof args.ligRef === 'string' ? " --lig-ref '" + sq(args.ligRef) + "'" : '')
              + (typeof args.recModel === 'string' ? " --rec-model '" + sq(args.recModel) + "'" : '')
              + (typeof args.ligModel === 'string' ? " --lig-model '" + sq(args.ligModel) + "'" : '')
            : ''
          const command = [
            "& '" + cfg.venvPy + "' '" + cfg.resPqDir + "\\struct_eval.py' --model '" + sq(args.model)
              + "' --ref '" + sq(args.ref) + "'" + (args.complex ? ' --complex' : '') + mc + rc + dq
              + " --out '" + sq(jsonPath) + "'",
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 900000)
          const value = outOf(result)
          const report = await readJson(jsonPath)
          value.report = report !== null ? report : undefined
          value.jsonPath = jsonPath
          return value
        },
      })))

      // ── 6. vscreen_run: batch virtual screening (library → ranked) ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'vscreen_run',
        description: 'Batch virtual screening with AutoDock Vina: receptor PDB (auto rigid PDBQT, co-ligand/ions excluded), ligand library CSV (one SMILES per row, meeko preparation), docking, incremental results.csv with resume support, ranked summary and top poses as PDB. Positive control validated on 3PTB (benzamidine ranked #1).',
        parameters: {
          receptor: { type: 'string', required: true, description: 'Receptor PDB path (HETATM water/ligand auto-excluded).' },
          ligands: { type: 'string', required: true, description: 'Library CSV path (column with "smiles", plus optional id/name).' },
          refLigand: { type: 'string', description: 'PDB of a co-crystallized ligand; search box auto-derived (preferred over explicit center/size).' },
          center: { type: 'array', items: { type: 'number' }, description: 'Explicit box center [x,y,z] (alternative to refLigand).' },
          size: { type: 'array', items: { type: 'number' }, description: 'Explicit box size [x,y,z] (alternative to refLigand).' },
          excludeRes: { type: 'string', description: 'Comma list of residue names excluded from receptor, e.g. "HOH,WAT,BEN,CA,SO4".' },
          exhaustiveness: { type: 'number', description: 'Vina exhaustiveness. Default 8 (screening); 16 for re-docking hits.' },
          top: { type: 'number', description: 'Top N poses exported as PDB. Default 5.' },
          outdir: { type: 'string', description: 'Output directory (resumable). Default <jobs dir>\\vscreen_<timestamp>.' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 3600000 (1h).' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          const outdir = typeof args.outdir === 'string' && args.outdir.length > 0
            ? args.outdir
            : cfg.jobsDir + '\\vscreen_' + Date.now()
          let box = ''
          if (typeof args.refLigand === 'string' && args.refLigand.length > 0) {
            box = " --ref-ligand '" + sq(args.refLigand) + "'"
          } else if (Array.isArray(args.center) && args.center.length === 3
            && Array.isArray(args.size) && args.size.length === 3) {
            box = ' --center ' + args.center.join(' ') + ' --size ' + args.size.join(' ')
          } else {
            throw new Error('provide refLigand or center+size')
          }
          const excl = typeof args.excludeRes === 'string' && args.excludeRes.length > 0
            ? " --exclude-res '" + sq(args.excludeRes) + "'" : ''
          const exh = Number.isFinite(args.exhaustiveness) ? ' --exhaustiveness ' + Math.floor(args.exhaustiveness) : ''
          const top = Number.isFinite(args.top) ? ' --top ' + Math.floor(args.top) : ''
          const jsonPath = outdir + '\\report.json'
          const command = [
            "& '" + cfg.venvPy + "' '" + cfg.resCiDir + "\\virtual_screen.py' --receptor '" + sq(args.receptor)
              + "' --ligands '" + sq(args.ligands) + "'" + box + excl + exh + top
              + " --outdir '" + sq(outdir) + "' --out '" + sq(jsonPath) + "'",
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 3600000)
          const value = outOf(result)
          const report = await readJson(jsonPath)
          value.report = report !== null ? report : undefined
          value.outdir = outdir
          return value
        },
      })))

      // ── 7. md_run: OpenMM MM-GBSA binding free energy / explicit MD ──
      disposers.push(ctx.tools.register(defineToolDef({
        name: 'md_run',
        description: 'OpenMM binding free energy and MD: mode "gb" = MM-GBSA (OBC2 implicit solvent; dG_bind with internal/nonbonded decomposition, minutes) for protein-protein complexes; mode "md" = explicit-solvent protocol (TIP3P box, 0.15 M NaCl, heating -> equilibration -> production, DCD trajectory, RMSD/RMSF analysis with plots). PDB sanitization (MSE, waters, incomplete residues, OXT capping) built in.',
        parameters: {
          mode: { type: 'string', required: true, description: '"gb" (MM-GBSA) or "md" (explicit solvent).' },
          complex: { type: 'string', required: true, description: 'Complex PDB path (protein chains).' },
          recChains: { type: 'array', items: { type: 'string' }, description: 'Receptor chain ids (gb mode). Default ["A"].' },
          ligChains: { type: 'array', items: { type: 'string' }, description: 'Ligand chain ids (gb mode). Default ["B"].' },
          steps: { type: 'number', description: 'md mode: production steps (2 fs each). Default 100000 (200 ps); real projects >= 25M (50 ns).' },
          platform: { type: 'string', description: 'OpenMM platform. Default CPU; OpenCL for speed (large systems).' },
          out: { type: 'string', description: 'gb mode: JSON report path. Default <jobs dir>\\md\\mmgbsa_<ts>.json' },
          outdir: { type: 'string', description: 'md mode: output directory. Default <jobs dir>\\md_<ts>.' },
          timeoutMs: { type: 'number', description: 'Timeout in milliseconds. Default 3600000 (1h); long MD needs background execution.' },
        },
        output: {
          schema: { type: 'object', additionalProperties: true },
          render: renderText,
        },
        async execute(args, exec) {
          if (args.mode !== 'gb' && args.mode !== 'md') throw new Error('mode must be "gb" or "md"')
          if (args.mode === 'gb') {
            const jsonPath = typeof args.out === 'string' && args.out.length > 0
              ? args.out
              : cfg.jobsDir + '\\md\\mmgbsa_' + Date.now() + '.json'
            const rc = Array.isArray(args.recChains) && args.recChains.length > 0
              ? " --rec-chains '" + args.recChains.map(sq).join("' '") + "'" : ''
            const lc = Array.isArray(args.ligChains) && args.ligChains.length > 0
              ? " --lig-chains '" + args.ligChains.map(sq).join("' '") + "'" : ''
            const pf = typeof args.platform === 'string' && args.platform.length > 0
              ? ' --platform ' + args.platform : ''
            const command = [
              "& '" + cfg.venvPy + "' '" + cfg.resDir + "\\md_mmgbsa.py' --mode gb --complex '" + sq(args.complex)
                + "'" + rc + lc + pf + " --out '" + sq(jsonPath) + "'",
              'exit $LASTEXITCODE',
            ].join('\n')
            const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 3600000)
            const value = outOf(result)
            const report = await readJson(jsonPath)
            value.report = report !== null ? report : undefined
            value.jsonPath = jsonPath
            return value
          }
          const outdir = typeof args.outdir === 'string' && args.outdir.length > 0
            ? args.outdir
            : cfg.jobsDir + '\\md_' + Date.now()
          const steps = Number.isFinite(args.steps) ? ' --steps ' + Math.floor(args.steps) : ''
          const pf = typeof args.platform === 'string' && args.platform.length > 0
            ? ' --platform ' + args.platform : ''
          const command = [
            "& '" + cfg.venvPy + "' '" + cfg.resDir + "\\md_mmgbsa.py' --mode md --complex '" + sq(args.complex)
              + "'" + steps + pf + " --outdir '" + sq(outdir) + "'",
            'exit $LASTEXITCODE',
          ].join('\n')
          const result = await runShell(command, exec, args.timeoutMs !== undefined ? args.timeoutMs : 3600000)
          const value = outOf(result)
          const report = await readJson(outdir + '\\md_report.json')
          value.report = report !== null ? report : undefined
          value.outdir = outdir
          return value
        },
      })))

      return () => disposers.forEach((d) => d())
    })
  },
}
