from shiny import App, ui, render, reactive
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------
# 1) Preloaded / Demo dataset
# ---------------------------------------------------------------------
# Pfad zu deinem Demo-.h5ad (anpassen!)
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
    ),

    ui.layout_column_wrap(
        ui.card(
            ui.card_header("Overview of obs"),
            ui.output_table("obs_preview"),
        ),

        ui.card(
            ui.card_header("2. Visualization"),
            ui.output_ui("color_by_ui"),
            ui.br(),
            ui.output_text_verbatim("plot_info"),
            ui.output_plot("visium_plot"),
        ),

        width=1/2,
    ),

    title="Spatial Domain Detection",
)


# ---------------------------------------------------------------------
# 3) Server
# ---------------------------------------------------------------------
def server(input, output, session):

    @reactive.Calc
    def adata():
        """Gibt entweder das Demo-Dataset oder das hochgeladene zurück."""
        file_list = input.file()

        # Kein Upload -> Demo-Dataset
        if not file_list:
            return demo_adata

        # Upload vorhanden -> hochgeladenes Dataset benutzen
        file_info = file_list[0]
        path = file_info["datapath"]
        return read_visium_h5ad(path)

    @output
    @render.text
    def info():
        a = adata()
        if a is None:
            return (
                "No dataset available.\n"
                "Check that the demo dataset can be loaded, "
                "or upload a .h5ad file."
            )

        # check, ob gerade ein File hochgeladen ist
        file_list = input.file()
        if file_list:
            source = "Uploaded dataset"
        else:
            source = f"Demo dataset ({DEMO_PATH})"

        n_cells = a.n_obs
        n_genes = a.n_vars
        return (
            f"{source}\n"
            f"Cells/Spots: {n_cells}\n"
            f"Genes: {n_genes}\n"
            f"obs columns: {list(a.obs.columns)}"
        )

    @output
    @render.table
    def obs_preview():
        a = adata()
        if a is None:
            return None
        obs_df: pd.DataFrame = a.obs.copy()
        obs_df = obs_df.reset_index().rename(columns={"index": "cell_id"})
        return obs_df.head(5)

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

        try:
            color_by = input.color_by()
        except Exception:
            color_by = None

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
        if a is None:
            return

        if "spatial" not in a.obsm.keys():
            return

        coords = a.obsm["spatial"]
        x = coords[:, 0]
        y = coords[:, 1]

        try:
            color_by = input.color_by()
        except Exception:
            color_by = None

        fig, ax = plt.subplots(figsize=(5, 5))

        if color_by and color_by in a.obs.columns:
            s = a.obs[color_by]

            # categorical
            if s.dtype.name == "category" or s.dtype == object:
                categories = s.astype("category")
                codes = categories.cat.codes.to_numpy()
                sc_plt = ax.scatter(x, y, c=codes, s=10, alpha=0.8)

                handles = []
                for cat_value in categories.cat.categories:
                    idx = (categories == cat_value).to_numpy().nonzero()[0][0]
                    norm_value = codes[idx] / codes.max() if codes.max() > 0 else 0
                    color = sc_plt.cmap(norm_value)
                    handles.append(
                        plt.Line2D(
                            [], [], marker="o", color=color,
                            linestyle="", label=str(cat_value), markersize=6
                        )
                    )
                ax.legend(
                    handles=handles,
                    title=color_by,
                    bbox_to_anchor=(1.05, 1),
                    loc="upper left",
                )
            else:
                sc_plt = ax.scatter(x, y, c=s.to_numpy(), s=10, alpha=0.8)
                cbar = plt.colorbar(sc_plt, ax=ax)
                cbar.set_label(color_by)
        else:
            ax.scatter(x, y, s=10, alpha=0.8)

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Visium Spatial Plot")
        ax.invert_yaxis()
        ax.set_aspect("equal", "box")
        return fig


app = App(app_ui, server)
