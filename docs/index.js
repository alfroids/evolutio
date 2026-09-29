importScripts("https://cdn.jsdelivr.net/pyodide/v0.29.3/full/pyodide.js");

function sendPatch(patch, buffers, msg_id) {
  self.postMessage({
    type: 'patch',
    patch: patch,
    buffers: buffers
  })
}

async function startApplication() {
  console.log("Loading pyodide...");
  self.postMessage({type: 'status', msg: 'Loading pyodide'})
  self.pyodide = await loadPyodide();
  self.pyodide.globals.set("sendPatch", sendPatch);
  console.log("Loaded pyodide!");
  const data_archives = ['index.resources.zip'];
  for (const archive of data_archives) {
    let zipResponse = await fetch(archive);
    let zipBinary = await zipResponse.arrayBuffer();
    self.postMessage({type: 'status', msg: `Unpacking ${archive}`})
    self.pyodide.unpackArchive(zipBinary, "zip");
  }
  await self.pyodide.loadPackage("micropip");
  self.postMessage({type: 'status', msg: `Installing environment`})
  try {
    await self.pyodide.runPythonAsync(`
      import micropip
      await micropip.install(['https://cdn.holoviz.org/panel/wheels/bokeh-3.9.2-py3-none-any.whl', 'https://cdn.holoviz.org/panel/1.9.4/dist/wheels/panel-1.9.4-py3-none-any.whl', 'pyodide-http', 'numpy', 'plotly', 'pandas']);
    `);
  } catch(e) {
    console.log(e)
    self.postMessage({
      type: 'status',
      msg: `Error while installing packages`
    });
  }
  console.log("Environment loaded!");
  self.postMessage({type: 'status', msg: 'Executing code'})
  try {
    const [docs_json, render_items, root_ids] = await self.pyodide.runPythonAsync(`\nimport asyncio\n\nfrom panel.io.pyodide import init_doc, write_doc\n\ninit_doc()\n\n"""\nInteractive web app for the evolutionary models in evolution_models.\n\nEach tab is built by its own build_tab() function in tabs/, which creates\nall of its widgets, panes, and callbacks as local variables and returns\n(title, layout). This file only wires those tabs together - it holds no\nwidget or pane instances of its own, so there is no way for two tabs to\naccidentally end up sharing the same underlying Bokeh model (the bug\nbehind a tab appearing during the loading skeleton and then vanishing\nonce rendering completes).\n"""\n\nimport sys\n\n# Plain relative paths, not something built from __file__: under panel\n# convert's pyodide-worker target, this script is exec'd from a string\n# rather than run as a real file, so __file__ is never defined there. The\n# --resources bundle unpacks preserving this project's layout relative to\n# the working directory, so relative paths work in both that environment\n# and a normal \`panel serve\` run from the project root.\nsys.path.insert(0, "src")\nsys.path.insert(0, ".")\n\nimport panel as pn\n\nfrom tabs import (\n\tone_locus_drift,\n\tone_locus_drift_with_bottleneck,\n\tone_locus_drift_with_selection,\n\tone_locus_selection,\n\ttwo_loci_selection,\n)\n\npn.extension("plotly")\n\napp = pn.Tabs(\n\tone_locus_selection.build_tab(),\n\ttwo_loci_selection.build_tab(),\n\tone_locus_drift.build_tab(),\n\tone_locus_drift_with_bottleneck.build_tab(),\n\tone_locus_drift_with_selection.build_tab(),\n)\napp.servable(title="Evolutio")\n\n\nawait write_doc()`)
    self.postMessage({
      type: 'render',
      docs_json: docs_json,
      render_items: render_items,
      root_ids: root_ids
    })
  } catch(e) {
    const traceback = `${e}`
    const tblines = traceback.split('\n')
    self.postMessage({
      type: 'status',
      msg: tblines[tblines.length-2]
    });
    throw e
  }
}

self.onmessage = async (event) => {
  const msg = event.data
  if (msg.type === 'rendered') {
    self.pyodide.runPythonAsync(`
    from panel.io.state import state
    from panel.io.pyodide import _link_docs_worker

    _link_docs_worker(state.curdoc, sendPatch, setter='js')
    `)
  } else if (msg.type === 'patch') {
    self.pyodide.globals.set('patch', msg.patch)
    self.pyodide.runPythonAsync(`
    from panel.io.pyodide import _convert_json_patch
    state.curdoc.apply_json_patch(_convert_json_patch(patch), setter='js')
    `)
    self.postMessage({type: 'idle'})
  } else if (msg.type === 'location') {
    self.pyodide.globals.set('location', msg.location)
    self.pyodide.runPythonAsync(`
    import json
    from panel.io.state import state
    from panel.util import edit_readonly
    if state.location:
        loc_data = json.loads(location)
        with edit_readonly(state.location):
            state.location.param.update({
                k: v for k, v in loc_data.items() if k in state.location.param
            })
    `)
  }
}

startApplication()