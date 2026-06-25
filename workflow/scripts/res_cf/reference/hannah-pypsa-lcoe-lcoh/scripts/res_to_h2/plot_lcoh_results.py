import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


RESULTS_DIR = Path("results")
PLOTS_DIR = Path("results/plots/lcoh_plots")


def plot_scenario(file_path):
    df = pd.read_csv(file_path)

    # keep only feasible rows with valid LCOH
    df = df[df["storage_feasible"] == True].copy()
    df = df[df["lcoh_with_storage_eur_per_kg"].notna()].copy()

    if df.empty:
        print(f"Skipped (no feasible rows): {file_path.name}")
        return

    x = df["res_to_el_ratio"]
    lcoh = df["lcoh_with_storage_eur_per_kg"]
    storage_duration = df["storage_capacity_required_mwh"] / df["electrolyser_mw"]

    # find optimum
    best_idx = lcoh.idxmin()
    best_ratio = df.loc[best_idx, "res_to_el_ratio"]
    best_lcoh = df.loc[best_idx, "lcoh_with_storage_eur_per_kg"]
    best_storage_duration = storage_duration.loc[best_idx]

    fig, ax1 = plt.subplots(figsize=(10, 5))

    # left axis: LCOH
    ax1.plot(x, lcoh, label="LCOH (€/kg)", color="blue")
    ax1.set_xlabel("RES / Electrolyser ratio")
    ax1.set_ylabel("LCOH [€/kg]")
    ax1.grid(True)

    # mark optimum on LCOH curve
    ax1.scatter(best_ratio, best_lcoh, zorder=5)

    # right axis: storage duration
    ax2 = ax1.twinx()
    ax2.plot(x, storage_duration, linestyle="--", color="red", label="Storage duration (h)")
    ax2.set_ylabel("Storage duration [h]")

    # combined legend
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper left")

    # fixed info box
    ax2.scatter(best_ratio, best_storage_duration, zorder=5)
    ax1.text(
        0.98,
        0.98,
        f"Optimal point\nratio = {best_ratio:.1f}\nLCOH = {best_lcoh:.2f} €/kg\n{best_storage_duration:.1f} h storage",
        transform=ax1.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9),
    )

    plt.title(file_path.stem)
    fig.tight_layout()

    save_path = PLOTS_DIR / f"{file_path.stem}.png"
    plt.savefig(save_path, dpi=150)
    plt.close()

    print(f"Saved: {save_path}")


def main_lineplots():
    for file in RESULTS_DIR.glob("*_stage3b_sweep.csv"):
        plot_scenario(file)


def extract_best_point(file_path):
    df = pd.read_csv(file_path)

    # keep only feasible rows
    df = df[df["storage_feasible"] == True].copy()
    df = df[df["lcoh_with_storage_eur_per_kg"].notna()].copy()

    if df.empty:
        return None

    best_idx = df["lcoh_with_storage_eur_per_kg"].idxmin()
    row = df.loc[best_idx]

    
    # use precomputed with-storage LCOH component contributions
    res_lcoh = row["lcoh_res_component_with_storage_eur_per_kg"]
    el_lcoh = row["lcoh_electrolyser_component_with_storage_eur_per_kg"]
    storage_lcoh = row["lcoh_storage_component_with_storage_eur_per_kg"]

    parts = file_path.stem.replace("_stage3b_sweep", "").split("_")

    # remove year if present
    if parts[-1] == "2023":
        parts = parts[:-1]

    country = parts[0]
    variant = parts[-2] + "_" + parts[-1] if parts[-2] == "bestsite" else parts[-1]

    if variant == "avg":
        tech = "_".join(parts[1:-1])
    else:
        tech = "_".join(parts[1:-2])

    return {
        "scenario": file_path.stem.replace("_stage3b_sweep", ""),
        "country": country,
        "tech": tech,
        "variant": variant,
        "res": res_lcoh,
        "electrolyser": el_lcoh,
        "storage": storage_lcoh,
        "total": res_lcoh + el_lcoh + storage_lcoh,
    }


def main_stacked():
    records = []

    for file in RESULTS_DIR.glob("*_stage3b_sweep.csv"):
        rec = extract_best_point(file)
        if rec:
            records.append(rec)

    df = pd.DataFrame(records)

    if df.empty:
        print("No valid scenarios found.")
        return

    variant_order = {"avg": 0, "bestsite_p95": 1}
    df["variant_order"] = df["variant"].map(variant_order)

    df = df.sort_values(["country", "tech", "variant_order"]).copy()

    df["label"] = (
        df["country"] + "_" +
        df["tech"] + "_" +
        df["variant"]
    )

    x = range(len(df))

    fig, ax = plt.subplots(figsize=(11, 6))

    # stacked bars
    ax.bar(x, df["res"], label="RES", color="green")
    ax.bar(x, df["electrolyser"], bottom=df["res"], label="Electrolyser", color="blue")
    ax.bar(
        x,
        df["storage"],
        bottom=df["res"] + df["electrolyser"],
        label="Storage",
        color="orange",
    )

    # give extra headroom for total labels
    ax.set_ylim(0, df["total"].max() * 1.22)

    # annotate totals on top of bars
    offset = df["total"].max() * 0.01
    for i, total in enumerate(df["total"]):
        ax.text(
            i,
            total + offset,
            f"{total:.1f}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )

    # labels
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["label"], rotation=45, ha="right")
    ax.set_ylabel("LCOH [€/kg]")
    ax.set_title("LCOH breakdown by scenario")
    ax.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    save_path = PLOTS_DIR / "lcoh_stacked_comparison.png"
    plt.savefig(save_path, dpi=150)
    plt.close()

    print(f"Saved: {save_path}")


def main():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    main_lineplots()
    main_stacked()


if __name__ == "__main__":
    main()