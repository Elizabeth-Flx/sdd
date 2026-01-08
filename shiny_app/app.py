from shiny import App, ui, render, reactive
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------
# 1) Preloaded / Demo dataset
# ---------------------------------------------------------------------
DEMO_PATH = "data/adata.h5ad"

try:
    demo_adata = sc.read_h5ad(DEMO_PATH)
except Exception as e:
    print(f"Could not load demo dataset from {DEMO_PATH}: {e}")
    demo_adata = None


def read_visium_h5ad(path: str):
    return sc.read_h5ad(path)


# ---------------------------------------------------------------------
# 2) UI
# ---------------------------------------------------------------------
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h4("1. Data input"),
        ui.input_file(
            "file",
            "Upload a spatial dataset (.h5ad)",
            multiple=False,
            accept=[".h5ad"],
        ),
        ui.br(),
        ui.output_text_verbatim("info"),
        ui.h4("2. Tools"),
        ui.input_checkbox_group(
            "tools",
            "Select tools to run:",
            choices={
                "bayesspace": "BayesSpace",
                "bass": "BASS",
                "graphst": "GraphST",
                "scanit": "SCAN-IT",
                "stlearn": "stLearn",
            },
        ),
        ui.input_action_button("run_tools", "Run selected tools"),
        ui.br(),
        ui.output_text_verbatim("tools_status"),
    ),
    ui.layout_column_wrap(
        ui.card(
            ui.card_header("Visualization"),
            ui.output_ui("color_by_ui"),
            ui.br(),
            ui.output_text_verbatim("plot_info"),
            ui.output_plot("visium_plot"),
        ),
        ui.card(
            ui.card_header("SCAN-IT result"),
            ui.output_text_verbatim("scanit_plot_info"),
            ui.output_plot("scanit_plot"),
        ),
        ui.card(
            ui.card_header("GraphST result"),
            ui.output_text_verbatim("graphst_plot_info"),
            ui.output_plot("graphst_plot"),
        ),
        ui.card(
            ui.card_header("stLearn result"),
            ui.output_text_verbatim("stlearn_plot_info"),
            ui.output_plot("stlearn_plot"),
        ),

        width=1 / 2,
    ),
    title="Spatial Domain Detection",
)


# ---------------------------------------------------------------------
# 3) Server
# ---------------------------------------------------------------------
def server(input, output, session):
    # -----------------------
    # Tool execution handling
    # -----------------------
    tools_result = reactive.Value({})  

    def run_tool(tool_key: str, adata_path: str) -> dict:
        project_root = Path(__file__).resolve().parents[1] 
        out_dir = (project_root / "results")
        out_dir.mkdir(parents=True, exist_ok=True)

        outputs = {
            "scanit": str((out_dir / "scanit_out.h5ad").resolve()),
            "graphst": str((out_dir / "graphst_out.h5ad").resolve()),
            "stlearn": str((out_dir / "stlearn_out.h5ad").resolve()),
        }

        cmds = {
            "scanit": [
                "micromamba", "run", "-n", "t_scanit",
                "python", "-m", "tools.t_scanit",
                "--input", adata_path,
                "--output", outputs["scanit"],
            ],
            "graphst": [
                "micromamba", "run", "-n", "t_graphst",
                "python", "-m", "tools.t_graphst",
                "--input", adata_path,
                "--output", outputs["graphst"]
            ],
            "stlearn": [
                "micromamba", "run", "-n", "t_stlearn",
                "python", "-m", "tools.t_stlearn",
                "--input", adata_path,
                "--output", outputs["stlearn"]
            ],

        }

        cmd = cmds[tool_key]

        try:
            p = subprocess.run(
                cmd,
                cwd=project_root,
                capture_output=True,
                text=True,
                check=True,
            )
            return {
                "ok": True,
                "cmd": cmd,
                "returncode": p.returncode,
                "stdout": (p.stdout or "")[-4000:],
                "stderr": (p.stderr or "")[-4000:],
                "output": outputs.get(tool_key),
            }

        except subprocess.CalledProcessError as e:
            return {
                "ok": False,
                "cmd": cmd,
                "returncode": e.returncode,
                "stdout": (e.stdout or "")[-4000:],
                "stderr": (e.stderr or "")[-4000:],
                "output": outputs.get(tool_key),
            }

        except FileNotFoundError as e:
            return {
                "ok": False,
                "cmd": cmd,
                "error": f"FileNotFoundError: {e}",
                "output": outputs.get(tool_key),
            }

        except Exception as e:
            return {
                "ok": False,
                "cmd": cmd,
                "error": repr(e),
                "output": outputs.get(tool_key),
            }

    @reactive.Effect
    @reactive.event(input.run_tools)
    def _run_selected_tools():
        a = adata()
        if a is None:
            tools_result.set({"error": {"ok": False, "error": "No dataset loaded."}})
            return

        selected = input.tools() or []
        if not selected:
            tools_result.set({})
            return

        file_list = input.file()
        adata_path = file_list[0]["datapath"] if file_list else DEMO_PATH

        results = {}
        for tool_key in selected:
            if tool_key not in ("scanit", "graphst", "stlearn"):
                results[tool_key] = {"ok": False, "error": "Tool not wired up yet."}
                continue
            results[tool_key] = run_tool(tool_key, adata_path)

        tools_result.set(results)

    @output
    @render.text
    def tools_status():
        res = tools_result.get()
        if not res:
            return "No tools executed yet."

        blocks = []
        for tool_name, r in res.items():
            lines = [f"=== {tool_name.upper()} ==="]
            lines.append("[OK]" if r.get("ok") else "[FAIL]")

            if r.get("cmd"):
                lines.append("\nCMD:")
                lines.append(" ".join(map(str, r["cmd"])))

            if r.get("returncode") is not None:
                lines.append(f"\nRETURN CODE: {r['returncode']}")

            if r.get("error"):
                lines.append("\nERROR:")
                lines.append(str(r["error"]))

            if r.get("stdout"):
                lines.append("\nSTDOUT:")
                lines.append(str(r["stdout"]))

            if r.get("stderr"):
                lines.append("\nSTDERR:")
                lines.append(str(r["stderr"]))

            if r.get("output"):
                lines.append("\nOUTPUT FILE:")
                lines.append(str(r["output"]))

            blocks.append("\n".join(lines))

        return "\n\n".join(blocks)

    # -----------------------
    # Data
    # -----------------------
    @reactive.Calc
    def adata():
        file_list = input.file()
        if not file_list:
            return demo_adata
        file_info = file_list[0]
        return read_visium_h5ad(file_info["datapath"])

    @output
    @render.text
    def info():
        a = adata()
        if a is None:
            return (
                "No dataset available.\n"
                "Check that the demo dataset can be loaded, or upload a .h5ad file."
            )

        file_list = input.file()
        source = "Uploaded dataset" if file_list else f"Demo dataset ({DEMO_PATH})"

        return (
            f"{source}\n"
            f"Cells/Spots: {a.n_obs}\n"
            f"Genes: {a.n_vars}\n"
            f"obs columns: {list(a.obs.columns)}"
        )

    # -----------------------
    # Main visualization (current adata)
    # -----------------------
    @output
    @render.ui
    def color_by_ui():
        a = adata()
        if a is None:
            return ui.div("No dataset available.")
        cols = list(a.obs.columns)
        if not cols:
            return ui.div("No obs columns found.")
        default_col = "leiden" if "leiden" in cols else cols[0]

        return ui.input_select(
            "color_by",
            "Color by obs column:",
            {c: c for c in cols},
            selected=default_col,
        )

    @output
    @render.text
    def plot_info():
        a = adata()
        if a is None:
            return "No data loaded."
        if "spatial" not in a.obsm.keys():
            return "No 'spatial' coordinates found in adata.obsm."

        color_by = None
        try:
            color_by = input.color_by()
        except Exception:
            pass

        text = "Visium Spatial Plot (points = spots)"
        if color_by and color_by in a.obs.columns:
            text += f"\nColoring by obs['{color_by}']."
        else:
            text += "\nNo valid obs column selected for coloring."
        return text

    @output
    @render.plot
    def visium_plot():
        a = adata()
        if a is None or "spatial" not in a.obsm.keys():
            return

        coords = a.obsm["spatial"]
        x = coords[:, 0]
        y = coords[:, 1]

        color_by = None
        try:
            color_by = input.color_by()
        except Exception:
            pass

        fig, ax = plt.subplots(figsize=(5, 5))

        if color_by and color_by in a.obs.columns:
            s = a.obs[color_by]

            if s.dtype.name == "category" or s.dtype == object:
                categories = s.astype("category")
                codes = categories.cat.codes.to_numpy()
                sc_plt = ax.scatter(x, y, c=codes, s=10, alpha=0.85)

                handles = []
                max_code = codes.max() if codes.size else 0
                for cat_value in categories.cat.categories:
                    idxs = np.where(categories.to_numpy() == cat_value)[0]
                    if len(idxs) == 0:
                        continue
                    idx = idxs[0]
                    norm_value = codes[idx] / max_code if max_code > 0 else 0
                    color = sc_plt.cmap(norm_value)
                    handles.append(
                        plt.Line2D([], [], marker="o", color=color, linestyle="",
                                   label=str(cat_value), markersize=6)
                    )

                ax.legend(handles=handles, title=color_by, bbox_to_anchor=(1.05, 1), loc="upper left")
            else:
                sc_plt = ax.scatter(x, y, c=s.to_numpy(), s=10, alpha=0.85)
                cbar = plt.colorbar(sc_plt, ax=ax)
                cbar.set_label(color_by)
        else:
            ax.scatter(x, y, s=10, alpha=0.85)

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Visium Spatial Plot")
        ax.invert_yaxis()
        ax.set_aspect("equal", "box")
        return fig

    # -----------------------
    # SCAN-IT result plotting
    # -----------------------
    @reactive.Calc
    def scanit_result_adata():
        res = tools_result.get()
        r = res.get("scanit") if isinstance(res, dict) else None
        out_path = r.get("output") if (r and r.get("ok")) else None
        if not out_path:
            return None

        p = Path(out_path)
        if not p.exists():
            return None

        return sc.read_h5ad(str(p))

    @output
    @render.text
    def scanit_plot_info():
        a = scanit_result_adata()
        if a is None:
            return "No SCAN-IT result yet. Run SCAN-IT first."
        if "spatial" not in a.obsm.keys():
            return "SCAN-IT output has no spatial coordinates in adata.obsm['spatial']."
        if "scanit_louvain" not in a.obs.columns:
            return "SCAN-IT output has no obs['scanit_louvain'] column."
        return "Showing SCAN-IT clusters (scanit_louvain)."

    @output
    @render.plot
    def scanit_plot():
        a = scanit_result_adata()
        if a is None:
            return
        if "spatial" not in a.obsm.keys() or "scanit_louvain" not in a.obs.columns:
            return

        coords = a.obsm["spatial"]
        x = coords[:, 0]
        y = coords[:, 1]

        clusters = a.obs["scanit_louvain"].astype("category")
        codes = clusters.cat.codes.to_numpy()

        fig, ax = plt.subplots(figsize=(5, 5))
        sc_plt = ax.scatter(x, y, c=codes, s=12, alpha=0.9)

        handles = []
        max_code = codes.max() if codes.size else 0
        for cat_value in clusters.cat.categories:
            idxs = np.where(clusters.to_numpy() == cat_value)[0]
            if len(idxs) == 0:
                continue
            idx = idxs[0]
            norm_value = codes[idx] / max_code if max_code > 0 else 0
            color = sc_plt.cmap(norm_value)
            handles.append(
                plt.Line2D([], [], marker="o", color=color, linestyle="",
                           label=str(cat_value), markersize=6)
            )

        ax.legend(handles=handles, title="scanit_louvain", bbox_to_anchor=(1.05, 1), loc="upper left")

        ax.set_title("scanit_louvain")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.invert_yaxis()
        ax.set_aspect("equal", "box")
        return fig

    # -----------------------
    # GraphST result plotting
    # -----------------------
    @reactive.Calc
    def graphst_result_adata():
        res = tools_result.get()
        r = res.get("graphst") if isinstance(res, dict) else None
        out_path = r.get("output") if (r and r.get("ok")) else None

        if not out_path:
            return None

        p = Path(out_path)
        if not p.exists():
            return None

        return sc.read_h5ad(str(p))

    @output
    @render.text
    def graphst_plot_info():
        a = graphst_result_adata()
        if a is None:
            return "No GraphST result yet. Run GraphST first."

        if "spatial" not in a.obsm.keys():
            return "GraphST output has no spatial coordinates in adata.obsm['spatial']."

        if "graphst_cluster" not in a.obs.columns:
            return "GraphST output has no obs['graphst_cluster'] column."

        return "Showing GraphST clusters (graphst_cluster)."

    @output
    @render.plot
    def graphst_plot():
        a = graphst_result_adata()
        if a is None:
            return

        if "spatial" not in a.obsm.keys() or "graphst_cluster" not in a.obs.columns:
            return

        coords = a.obsm["spatial"]
        x = coords[:, 0]
        y = coords[:, 1]

        clusters = a.obs["graphst_cluster"].astype("category")
        codes = clusters.cat.codes.to_numpy()

        fig, ax = plt.subplots(figsize=(5, 5))
        sc_plt = ax.scatter(x, y, c=codes, s=12, alpha=0.9)

        handles = []
        max_code = codes.max() if codes.size else 0
        for cat_value in clusters.cat.categories:
            idxs = np.where(clusters.to_numpy() == cat_value)[0]
            if len(idxs) == 0:
                continue
            idx = idxs[0]
            norm_value = codes[idx] / max_code if max_code > 0 else 0
            color = sc_plt.cmap(norm_value)
            handles.append(
                plt.Line2D([], [], marker="o", color=color,
                        linestyle="", label=str(cat_value), markersize=6)
            )

        ax.legend(
            handles=handles,
            title="graphst_cluster",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
        )

        ax.set_title("graphst_cluster")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.invert_yaxis()
        ax.set_aspect("equal", "box")
        return fig

    # -----------------------
    # stLearn result plotting
    # -----------------------
    @reactive.Calc
    def stlearn_result_adata():
        res = tools_result.get()
        r = res.get("stlearn") if isinstance(res, dict) else None
        out_path = r.get("output") if (r and r.get("ok")) else None

        if not out_path:
            return None

        p = Path(out_path)
        if not p.exists():
            return None

        return sc.read_h5ad(str(p))

    @output
    @render.text
    def stlearn_plot_info():
        a = stlearn_result_adata()
        if a is None:
            return "No stLearn result yet. Run stLearn first."

        if "spatial" not in a.obsm.keys():
            return "stLearn output has no spatial coordinates in adata.obsm['spatial']."

        if "stlearn_louvain" not in a.obs.columns:
            return "stLearn output has no obs['stlearn_louvain'] column."

        return "Showing stLearn clusters (stlearn_louvain)."

    @output
    @render.plot
    def stlearn_plot():
        a = stlearn_result_adata()
        if a is None:
            return

        if "spatial" not in a.obsm.keys() or "stlearn_louvain" not in a.obs.columns:
            return

        coords = a.obsm["spatial"]
        x = coords[:, 0]
        y = coords[:, 1]

        clusters = a.obs["stlearn_louvain"].astype("category")
        codes = clusters.cat.codes.to_numpy()

        fig, ax = plt.subplots(figsize=(5, 5))
        sc_plt = ax.scatter(x, y, c=codes, s=12, alpha=0.9)

        handles = []
        max_code = codes.max() if codes.size else 0
        for cat_value in clusters.cat.categories:
            idxs = np.where(clusters.to_numpy() == cat_value)[0]
            if len(idxs) == 0:
                continue
            idx = idxs[0]
            norm_value = codes[idx] / max_code if max_code > 0 else 0
            color = sc_plt.cmap(norm_value)
            handles.append(
                plt.Line2D([], [], marker="o", color=color,
                        linestyle="", label=str(cat_value), markersize=6)
            )

        ax.legend(
            handles=handles,
            title="stlearn_louvain",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
        )

        ax.set_title("stlearn_louvain")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.invert_yaxis()
        ax.set_aspect("equal", "box")
        return fig





app = App(app_ui, server)
