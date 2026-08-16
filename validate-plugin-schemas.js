// validate-plugin-schemas.js — offline validator for plugins/protein-tools.js.
// Mocks a raw (unsandboxed) preset-file ctx (no `harness` global) and checks
// that every registered tool definition compiles to a root-object JSON Schema
// — the same shape dsh-tools' parameterSchemaSpecToJsonSchema produces. This
// is the regression guard for the upstream "Invalid schema ... got
// 'type: null'" API rejection.
// Usage: node validate-plugin-schemas.js

const plugin = require('./plugins/protein-tools.js')

const registered = []
const ctx = {
  get: () => undefined,
  tools: {
    register: (def) => {
      registered.push(def)
      return () => {}
    },
  },
  shell: {},
  effect: (fn) => {
    const stop = fn()
    return () => { if (typeof stop === 'function') stop() }
  },
}

plugin.apply(ctx)

const EXPECTED = [
  'esmfold_predict', 'pp_interact', 'vina_dock', 'af2_predict',
  'struct_eval', 'vscreen_run', 'md_run',
]

for (const name of EXPECTED) {
  if (!registered.some((d) => d.name === name)) {
    console.error('FAIL missing tool: ' + name)
    process.exit(1)
  }
}

for (const def of registered) {
  const p = def.parameters
  if (!p || typeof p !== 'object' || p.type !== 'object') {
    console.error('FAIL ' + def.name + ': parameters.type must be "object", got ' + (p && p.type))
    process.exit(1)
  }
  if (!p.properties || typeof p.properties !== 'object') {
    console.error('FAIL ' + def.name + ': missing parameters.properties')
    process.exit(1)
  }
  for (const key of Object.keys(p.properties)) {
    if (p.properties[key] && p.properties[key].required !== undefined) {
      console.error('FAIL ' + def.name + '.' + key + ': per-property `required` must be lifted to the root array')
      process.exit(1)
    }
  }
  if (typeof def.execute !== 'function') {
    console.error('FAIL ' + def.name + ': missing execute')
    process.exit(1)
  }
}

const af2 = registered.find((d) => d.name === 'af2_predict')
if (!Array.isArray(af2.parameters.required) || af2.parameters.required.indexOf('fasta') === -1) {
  console.error("FAIL af2_predict: required must include 'fasta'")
  process.exit(1)
}

const schema = af2.parameters
console.log('sample af2_predict wire schema:')
console.log(JSON.stringify({ type: schema.type, required: schema.required, keys: Object.keys(schema.properties) }))
console.log('ALL SCHEMAS OK (' + registered.length + ' tools)')
