from shiny import App, ui, render, reactive
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
MPL_BASE = mpl.rcParams.copy()
import numpy as np
import subprocess
from pathlib import Path
import tempfile
import shutil
from uuid import uuid4
import os

import sys

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))


from evaluation_metrics.utils import load_csvs_to_arrays
from evaluation_metrics.plots import plot_means_lineplot, plot_sankey_diagram


# ---------------------------------------------------------------------
# 1) Preloaded / Demo dataset
# ---------------------------------------------------------------------
DEMO_PATH = "data/adata.h5ad"

try:
    project_root = Path(__file__).resolve().parents[1]
    demo_adata = sc.read_h5ad(str((project_root / DEMO_PATH).resolve()))

except Exception as e:
    print(f"Could not load demo dataset from {DEMO_PATH}: {e}")
    demo_adata = None


def read_visium_h5ad(path: str):
    return sc.read_h5ad(path)

def run_metrics(anndata_dir: Path, csv_dir: Path, metrics_dir: Path) -> dict:
    project_root = Path(__file__).resolve().parents[1]

    metrics_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "micromamba", "run", "-n", "sdd_metrics",  
        "python", "-m", "evaluation_metrics.metrics_main",   
        "--input_anndata", str(anndata_dir),  
        "--input_csv", str(csv_dir),            
        "--output", str(metrics_dir),            
    ]

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
            "output_dir": str(metrics_dir),
        }
    except subprocess.CalledProcessError as e:
        return {
            "ok": False,
            "cmd": cmd,
            "returncode": e.returncode,
            "stdout": (e.stdout or "")[-4000:],
            "stderr": (e.stderr or "")[-4000:],
            "output_dir": str(metrics_dir),
        }
    except Exception as e:
        return {"ok": False, "cmd": cmd, "error": repr(e), "output_dir": str(metrics_dir)}

def read_bass_n_clusters(csv_dir: Path) -> int | None:
    p = csv_dir / "bass_n_clusters.txt"
    if not p.exists():
        return None
    try:
        txt = p.read_text().strip()
        digits = "".join([c for c in txt if c.isdigit()])
        return int(digits) if digits else None
    except Exception:
        return None




# ---------------------------------------------------------------------
# 2) UI
# ---------------------------------------------------------------------
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h4("1. Data input"),
        ui.input_file(
            "file",
            "Upload a spatial dataset (.h5ad)",
            multiple=True,
            accept=[".h5ad"],
        ),
        ui.input_file(
            "file_atac",
            "Upload an ATAC dataset (.h5ad)",
            multiple=True,
            accept=[".h5ad"],
        ),
        ui.output_ui("dataset_picker_ui"),
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
        ui.input_checkbox("run_metrics", "Calculate metrics after tools run", value=False),
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
        ui.card(
            ui.card_header("BayesSpace result"),
            ui.output_text_verbatim("bayesspace_plot_info"),
            ui.output_plot("bayesspace_plot"),
        ),
        ui.card(
            ui.card_header("BASS result"),
            ui.output_text_verbatim("bass_plot_info"),
            ui.output_plot("bass_plot"),
        ),
        ui.card(
            ui.card_header("Metrics plots"),
            ui.output_ui("gt_col_ui"),
            ui.input_select(
                "metrics_plot_kind",
                "Plot type:",
                choices={
                    "sankey": "Sankey (methods flow)",
                    "means_all": "Mean metrics (all)",
                    "means_acc": "Mean metrics (accuracy; needs ground truth metrics in CSV)",
                    "means_spatial": "Mean metrics (spatial)",
                    "means_sil": "Mean metrics (silhouette)",
                },
                selected="means_all",
            ),
            ui.output_ui("metrics_panel_ui"),
        ),



        width=1 / 2,
    ),
    title="Spatial Domain Detection",
)


# ---------------------------------------------------------------------
# 3) Server
# ---------------------------------------------------------------------
def server(input, output, session):
    # ----------------------- Session workdir ------------------------
    base_tmp = Path(tempfile.gettempdir()) / "spatial-shiny"
    base_tmp.mkdir(parents=True, exist_ok=True)

    session_dir = base_tmp / f"session-{uuid4().hex}"
    session_dir.mkdir(parents=True, exist_ok=True)

    # Clean up session directory on session end
    def clean_session_dir():
        shutil.rmtree(session_dir, ignore_errors=True)

    session.on_ended(clean_session_dir)


    # -----------------------
    # Tool execution handling
    # -----------------------

    tools_result = reactive.Value({})

    uploads_state = reactive.Value([])

    atac_uploads_state = reactive.Value([])

    @reactive.Effect
    def _track_uploads():
        file_list = input.file() or []
        if not file_list:
            return

        current = uploads_state.get()
        seen_srcs = {u["src_datapath"] for u in current}

        new_entries = []
        for info in file_list:
            src_datapath = info["datapath"]
            if src_datapath in seen_srcs:
                continue

            new_entries.append({
                "id": uuid4().hex,
                "name": info["name"],
                "src_datapath": src_datapath,
                "path": src_datapath,  
            })

        if new_entries:
            uploads_state.set(current + new_entries)
    
    @reactive.Effect
    def _track_atac_uploads():
        file_list = input.file_atac() or []
        if not file_list:
            return

        current = atac_uploads_state.get()
        seen_srcs = {u["src_datapath"] for u in current}

        new_entries = []
        for info in file_list:
            src_datapath = info["datapath"]
            if src_datapath in seen_srcs:
                continue

            new_entries.append({
                "id": uuid4().hex,
                "name": info["name"],
                "src_datapath": src_datapath,
                "path": src_datapath,
            })

        if new_entries:
            atac_uploads_state.set(current + new_entries)



    @output
    @render.ui
    def dataset_picker_ui():
        ups = uploads_state.get()
        if not ups:
            return ui.div("No uploads yet (using demo).")

        choices = {u["id"]: f"{u['name']} ({u['id'][:6]})" for u in ups}
        selected = ups[-1]["id"]   
        return ui.input_select("active_upload", "Active dataset:", choices, selected=selected)

    @reactive.Calc
    def active_adata_path():
        ups = uploads_state.get()
        if not ups:
            return DEMO_PATH

        try:
            active_id = input.active_upload()
        except Exception:
            active_id = ups[-1]["id"]

        for u in ups:
            if u["id"] == active_id:
                return u["path"]

        return ups[-1]["path"]

    @reactive.Calc
    def adata():
        res = tools_result.get()
        combined = res.get("_combined_h5ad") if isinstance(res, dict) else None

        if combined and Path(combined).exists():
            return sc.read_h5ad(combined)

        path = active_adata_path()
        if path == DEMO_PATH:
            return demo_adata
        return read_visium_h5ad(path)


  

    def run_tool(tool_key: str, adata_path: str, csv_dir: Path, r_input_dir: Path | None, n_clusters: int | None) -> dict:
        project_root = Path(__file__).resolve().parents[1]

        if tool_key in ("bayesspace", "bass") and r_input_dir is None:
            return {"ok": False, "error": "Missing r_input_dir (adata_to_csv not run)."}

        outputs = {
            "scanit": str((csv_dir / "scanit.csv").resolve()),
            "graphst": str((csv_dir / "graphst.csv").resolve()),
            "stlearn": str((csv_dir / "stlearn.csv").resolve()),
            "bayesspace": str((csv_dir / "bayesspace.csv").resolve()),
            "bass": str((csv_dir / "bass.csv").resolve()),
        }

        cmds = {
            "bass": [
                "micromamba", "run", "-n", "t_bayesspace_bass",
                "Rscript", "tools/t_bass.R",
                str(r_input_dir), outputs["bass"], "True"
            ],
            "scanit": [
                "micromamba", "run", "-n", "t_scanit",
                "python", "-m", "tools.t_scanit",
                "--input", adata_path,
                "--output", outputs["scanit"]
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
            "bayesspace": [
                "micromamba", "run", "-n", "t_bayesspace_bass",
                "Rscript", "tools/t_bayesspace.R",
                str(r_input_dir), outputs["bayesspace"], "True"
            ],
        }

        cmd = cmds[tool_key]

        if n_clusters is not None and tool_key in ("scanit", "graphst", "stlearn"):
            cmd += ["--k_clusters", str(n_clusters)]

        if tool_key == "bayesspace" and n_clusters is not None:
            cmd += [str(n_clusters)]


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
        n_clusters = None

        if a is None:
            tools_result.set({"error": {"ok": False, "error": "No dataset loaded."}})
            return

        selected = input.tools() or []

        if not selected:
            tools_result.set({})
            return

        run_id = f"run-{uuid4().hex}"
        run_dir = session_dir / run_id

        results = {"_run_dir": str(run_dir)}

        upload_visium_dir = run_dir / "upload_visium"
        upload_visium_dir.mkdir(parents=True, exist_ok=True)

        upload_atac_dir = run_dir / "upload_atac"
        upload_atac_dir.mkdir(parents=True, exist_ok=True)
        
        csv_dir = run_dir / "csv"
        csv_dir.mkdir(parents=True, exist_ok=True)

        metrics_dir = run_dir / "metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)

        project_root = Path(__file__).resolve().parents[1]

        src_path = active_adata_path() or DEMO_PATH
        src = Path(src_path)

        if src_path == DEMO_PATH:
            src = (project_root / DEMO_PATH).resolve()
        else:
            src = src.resolve()

        dst = upload_visium_dir / src.name
        shutil.copy2(src, dst)
        adata_path = str(dst)

        # --- copy ATAC upload (optional) into upload_atac_dir ---
        atac_path = None
        atac_ups = atac_uploads_state.get()

        if atac_ups:
            atac_src_path = atac_ups[-1]["path"]
            atac_src = Path(atac_src_path).resolve()

            atac_dst = upload_atac_dir / atac_src.name
            shutil.copy2(atac_src, atac_dst)
            atac_path = str(atac_dst)

        # helpful for debugging / downstream use
        results["_upload_visium_dir"] = str(upload_visium_dir)
        results["_upload_atac_dir"] = str(upload_atac_dir)
        results["_atac_h5ad"] = atac_path



        needs_r_input = any(t in ("bayesspace", "bass") for t in selected)
        r_input_dir = None

        if needs_r_input:
            r_input_dir = run_dir / "r_input"
            r_input_dir.mkdir(parents=True, exist_ok=True)

            adata_to_csv_cmd = [
                "micromamba", "run", "-n", "t_stlearn",  
                "python", "-m", "tools.adata_to_csv",
                "--input", adata_path,
                "--output_dir", str(r_input_dir),
            ]

            p_csv = subprocess.run(
                adata_to_csv_cmd,
                cwd=project_root,
                capture_output=True,
                text=True,
            )

            results["_adata_to_csv"] = {
                "ok": p_csv.returncode == 0,
                "returncode": p_csv.returncode,
                "stdout": (p_csv.stdout or "")[-4000:],
                "stderr": (p_csv.stderr or "")[-4000:],
                "input_dir": str(r_input_dir),
            }

            if p_csv.returncode != 0:
                tools_result.set(results)
                return

        # 1) Run Bass first to get n_clusters if Bass is selected
        if "bass" in selected:
            results["bass"] = run_tool("bass", adata_path, csv_dir, r_input_dir, n_clusters=None)

            n_clusters = read_bass_n_clusters(csv_dir)
            results["_bass_n_clusters"] = {
                "ok": n_clusters is not None,
                "value": n_clusters,
                "path": str((csv_dir / "bass_n_clusters.txt").resolve()),
            }

        # 2) Run other selected tools with n_clusters if available
        for tool_key in selected:
            if tool_key == "bass":
                continue

            if tool_key not in ("scanit", "graphst", "stlearn", "bayesspace"):
                results[tool_key] = {"ok": False, "error": "Tool not wired up yet."}
                continue

            results[tool_key] = run_tool(tool_key, adata_path, csv_dir, r_input_dir, n_clusters=n_clusters)


        combined_path = run_dir / "combined.h5ad"

        combine_cmd = [
            "micromamba", "run", "-n", "t_stlearn",
            "python", "-m", "tools.combine",
            "--input", adata_path,
            "--csv_dir", str(csv_dir),
            "--output", str(combined_path),
        ]

        p = subprocess.run(
            combine_cmd,
            cwd=project_root,
            capture_output=True,
            text=True,
        )

        combined_ok = (p.returncode == 0) and combined_path.exists()

        results["_combine"] = {
            "ok": combined_ok,
            "returncode": p.returncode,
            "stdout": (p.stdout or "")[-4000:],
            "stderr": (p.stderr or "")[-4000:],
        }

        if combined_ok:
            results["_combined_h5ad"] = str(combined_path)
        else:
            results["_combined_h5ad"] = None

        # -----------------------
        # Metrics (optional)
        # -----------------------
        do_metrics = False
        try:
            do_metrics = bool(input.run_metrics())
        except Exception:
            do_metrics = False

        if do_metrics:
            results["_metrics"] = run_metrics(
                anndata_dir=upload_visium_dir,
                csv_dir=csv_dir,
                metrics_dir=metrics_dir,
            )
        else:
            results["_metrics"] = {"ok": False, "skipped": True, "reason": "run_metrics unchecked"}


        tools_result.set(results)


    @output
    @render.text
    def tools_status():
        res = tools_result.get()
        if not res:
            return "No tools executed yet."

        run_dir = res.get("_run_dir") if isinstance(res, dict) else None

        blocks = []
        if run_dir:
            blocks.append(f"RUN DIR:\n{run_dir}\n")

        for tool_name, r in res.items():
            if tool_name in ("_run_dir", "_combined_h5ad"):
                continue

            if not isinstance(r, dict):
                blocks.append(f"=== {tool_name.upper()} ===\n{r}")
                continue

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

            if r.get("output_dir"):
                lines.append("\nOUTPUT DIR:")
                lines.append(str(r["output_dir"]))


            blocks.append("\n".join(lines))

        return "\n\n".join(blocks)


    @reactive.Calc
    def current_csv_dir():
        res = tools_result.get()
        run_dir = res.get("_run_dir") if isinstance(res, dict) else None
        if not run_dir:
            return None
        p = Path(run_dir) / "csv"
        return p if p.exists() else None

    @reactive.Calc
    def dict_all_methods():
        cdir = current_csv_dir()
        if not cdir:
            return {}

        files = [f for f in os.listdir(cdir) if f.endswith(".csv")]
        if not files:
            return {}

        dict_all_methods, groundtruth, ref_barcodes = load_csvs_to_arrays(str(cdir), files)
        return dict_all_methods


    @reactive.Calc
    def df_all_metrics_mean():
        res = tools_result.get()
        run_dir = res.get("_run_dir") if isinstance(res, dict) else None
        if not run_dir:
            return None
        p = Path(run_dir) / "metrics" / "all_metrics_summary.csv"
        if not p.exists():
            return None
        return pd.read_csv(p, index_col=0)

    # -----------------------
    # Ground truth (optional)
    # -----------------------
    @output
    @render.ui
    def gt_col_ui():
        a = adata()
        if a is None:
            return ui.div("No data.")
        cols = list(a.obs.columns)
        choices = {"__none__": "(none)"} | {c: c for c in cols}
        return ui.input_select("gt_col", "Ground truth column (optional):", choices, selected="__none__")

    def get_groundtruth_or_none(a):
        try:
            gt = input.gt_col()
        except Exception:
            return None
        if (not gt) or gt == "__none__":
            return None
        if a is None or gt not in a.obs.columns:
            return None
        return a.obs[gt].to_numpy()

    # -----------------------
    # Metrics plots (matplotlib / plotnine)
    # -----------------------
    @output
    @render.plot
    def metrics_matplotlib_plot():
        kind = input.metrics_plot_kind()

        if kind == "sankey":
            return

        with mpl.rc_context(MPL_BASE):
            df = df_all_metrics_mean()
            if df is None or df.empty:
                fig, ax = plt.subplots(figsize=(5, 5))
                ax.text(
                    0.5, 0.5,
                    "No metrics summary found yet.\nRun tools with 'Calculate metrics' enabled.",
                    ha="center", va="center", transform=ax.transAxes
                )
                ax.axis("off")
                return fig

            df_acc = df[df.columns.intersection(["Method", "ARI", "Purity", "AMI"])]
            df_spa = df[df.columns.intersection(["Method", "LISI", "PAS"])]
            df_sil = df[df.columns.intersection(["Method", "Silhouette_exp", "Silhouette_spat"])]

            if kind == "means_acc":
                p = plot_means_lineplot(df_acc, title="Mean Metric Scores (Accuracy)")
            elif kind == "means_spatial":
                p = plot_means_lineplot(df_spa, title="Mean Metric Scores (Spatial)")
            elif kind == "means_sil":
                p = plot_means_lineplot(df_sil, title="Mean Metric Scores (Silhouette)")
            else:
                p = plot_means_lineplot(df, title="Mean Metric Scores (All)")

            return p.draw()

    # -----------------------
    # Sankey (Plotly) as inline HTML
    # -----------------------
    @output
    @render.ui
    def metrics_plotly_ui():
        if input.metrics_plot_kind() != "sankey":
            return ui.div()

        methods = dict_all_methods()
        if not methods:
            return ui.div("No method CSVs found yet (run tools first).")

        a = adata()
        gt = get_groundtruth_or_none(a)  # None allowed

        fig = plot_sankey_diagram(
            methods,
            ground_truth=gt,
            title="Sankey" + (" vs Ground Truth" if gt is not None else ""),
        )

        html = fig.to_html(include_plotlyjs="inline", full_html=False, config={"responsive": True})
        return ui.HTML(html)

    @output
    @render.ui
    def metrics_panel_ui():
        kind = input.metrics_plot_kind()

        if kind == "sankey":
            return ui.card(
                ui.card_header("Sankey"),
                ui.output_ui("metrics_plotly_ui"),
            )

        return ui.card(
            ui.card_header("Metrics"),
            ui.output_plot("metrics_matplotlib_plot"),
        )




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

    def plot_spatial_clusters(a, col: str, title: str):
        mpl.rcParams.update(MPL_BASE)
        if a is None or "spatial" not in a.obsm:
            return None

        coords = a.obsm["spatial"]
        x, y = coords[:, 0], coords[:, 1]

        fig, ax = plt.subplots(figsize=(5, 5))

        if col in a.obs.columns:
            s = a.obs[col]
            cats = s.astype("category")
            codes = cats.cat.codes.to_numpy()
            sc_plt = ax.scatter(x, y, c=codes, s=12, alpha=0.9)

            handles = []
            max_code = codes.max() if codes.size else 0
            for cat_value in cats.cat.categories:
                idxs = np.where(cats.to_numpy() == cat_value)[0]
                if len(idxs) == 0:
                    continue
                idx = idxs[0]
                norm_value = codes[idx] / max_code if max_code > 0 else 0
                color = sc_plt.cmap(norm_value)
                handles.append(
                    plt.Line2D([], [], marker="o", color=color, linestyle="",
                            label=str(cat_value), markersize=6)
                )

            ax.legend(handles=handles, title=col, bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            ax.text(0.5, 0.5, f"Missing obs['{col}']",
                    ha="center", va="center", transform=ax.transAxes)

        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.invert_yaxis()
        ax.set_aspect("equal", "box")
        return fig


    # -----------------------
    # Main visualization 
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
        mpl.rcParams.update(MPL_BASE)
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
    @output
    @render.text
    def scanit_plot_info():
        a = adata()
        if a is None:
            return "No data loaded."
        if "spatial" not in a.obsm:
            return "No spatial coordinates."
        if "scanit" not in a.obs.columns:
            return "No SCAN-IT result yet (obs['scanit'] missing). Run SCAN-IT first."
        return "Showing SCAN-IT clusters (obs['scanit'])."

    @output
    @render.plot
    def scanit_plot():
        a = adata()
        fig = plot_spatial_clusters(a, "scanit", "SCAN-IT clusters")
        if fig is not None:
            return fig

    # -----------------------
    # GraphST result plotting
    # -----------------------
    @output
    @render.text
    def graphst_plot_info():
        a = adata()
        if a is None:
            return "No data loaded."
        if "spatial" not in a.obsm:
            return "No spatial coordinates."
        if "graphst" not in a.obs.columns:
            return "No GraphST result yet (obs['graphst'] missing). Run GraphST first."
        return "Showing GraphST clusters (obs['graphst'])."

    @output
    @render.plot
    def graphst_plot():
        a = adata()
        fig = plot_spatial_clusters(a, "graphst", "GraphST clusters")
        if fig is not None:
            return fig

    # -----------------------
    # stLearn result plotting
    # -----------------------
    @output
    @render.text
    def stlearn_plot_info():
        a = adata()
        if a is None:
            return "No data loaded."
        if "spatial" not in a.obsm:
            return "No spatial coordinates."
        if "stlearn" not in a.obs.columns:
            return "No stLearn result yet (obs['stlearn'] missing). Run stLearn first."
        return "Showing stLearn clusters (obs['stlearn'])."

    @output
    @render.plot
    def stlearn_plot():
        a = adata()
        fig = plot_spatial_clusters(a, "stlearn", "stLearn clusters")
        if fig is not None:
            return fig

    # -----------------------
    # BayesSpace result plotting
    # -----------------------
    @output
    @render.text
    def bayesspace_plot_info():
        a = adata()
        if a is None:
            return "No data loaded."
        if "spatial" not in a.obsm:
            return "No spatial coordinates."
        if "bayesspace" not in a.obs.columns:
            return "No BayesSpace result yet (obs['bayesspace'] missing). Run BayesSpace first."
        return "Showing BayesSpace clusters (obs['bayesspace'])."

    @output
    @render.plot
    def bayesspace_plot():
        a = adata()
        fig = plot_spatial_clusters(a, "bayesspace", "BayesSpace clusters")
        if fig is not None:
            return fig

    # -----------------------
    # Bass result plotting
    # -----------------------

    @output
    @render.text
    def bass_plot_info():
        a = adata()
        if a is None:
            return "No data loaded."
        if "spatial" not in a.obsm:
            return "No spatial coordinates."
        if "bass" not in a.obs.columns:
            return "No BASS result yet (obs['bass'] missing). Run BASS first."
        return "Showing BASS clusters (obs['bass'])."

    @output
    @render.plot
    def bass_plot():
        a = adata()
        fig = plot_spatial_clusters(a, "bass", "BASS clusters")
        if fig is not None:
            return fig


app = App(app_ui, server)
