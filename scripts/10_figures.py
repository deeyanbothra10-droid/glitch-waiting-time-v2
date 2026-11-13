"""Figures -> figures/*.pdf (vector PDF, greyscale, units on the axes), from results/ and data/processed/.

Allowed before freezing (no waiting times): F1 exclusion flow (logs/exclusions.csv, results/method_counts.json)
and F5 validation by scenario (results/validation.json, synthetic). Only when the guard passes
(glitchlib.guards): F2 ln T against ln R per pulsar with slope-1 lines (M1 intercepts), F3 Kaplan-Meier
curves for the four pulsars with the most intervals, F4 model comparison (AIC and LOO-CV for M0-M4), F6
implied beta per pulsar with its standard error.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import models  # noqa: E402
from glitchlib.analysis import kaplan_meier, primary_data  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, unblinding_allowed  # noqa: E402

FIG = ROOT / "figures"
plt.rcParams.update({"font.size": 9, "axes.prop_cycle": plt.cycler(color=["0.0", "0.45", "0.7"])})


def f1_exclusion_flow():
    log = pd.read_csv(ROOT / "logs" / "exclusions.csv")
    counts = json.loads((ROOT / "results" / "method_counts.json").read_text())["included_pulsars_totals"]
    steps = [f"Catalogue\n{int(log.loc[0, 'rows'])} glitches, {int(log.loc[0, 'pulsars'])} pulsars",
             f"Non-pulsar objects removed\n{int(log.loc[1, 'rows'])} glitches, {int(log.loc[1, 'pulsars'])} pulsars",
             f"Usable (dff > 0, df1 > 0)\n{int(log.loc[2, 'usable'])} glitches",
             f"Admissible Δν̇ method\n{int(log.loc[4, 'admissible'])} glitches",
             f"Included pulsars (D13)\n{counts['pulsars']} pulsars, {counts['admissible_completed']} completed "
             f"+ {counts['censored_at_obs_end'] + counts['censored_at_gap']} censored intervals"]
    fig, ax = plt.subplots(figsize=(4.2, 5.2))
    ax.axis("off")
    for k, text in enumerate(steps):
        y = 1 - k / (len(steps) - 0.6)
        ax.text(0.5, y, text, ha="center", va="center", bbox={"boxstyle": "round", "fc": "white", "ec": "black"})
        if k:
            ax.annotate("", xy=(0.5, y + 0.07), xytext=(0.5, y + 1 / (len(steps) - 0.6) - 0.07),
                        arrowprops={"arrowstyle": "->", "color": "black"})
    fig.savefig(FIG / "F1_exclusion_flow.pdf")
    plt.close(fig)


def f5_validation():
    record = json.loads((ROOT / "results" / "validation.json").read_text())
    sigmas = list(record["per_sigma"])
    names = list(record["per_sigma"][sigmas[0]]["scenarios"])
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6), sharey=True)
    markers = ["o", "s", "^", "D", "v", "<", ">"][:len(sigmas)]    # one per sigma of the validation grid (D14 + D4k)
    if len(markers) < len(sigmas):
        raise ValueError(f"F5 has markers for {len(markers)} sigma values, not {len(sigmas)}.")
    for k, sigma in enumerate(sigmas):
        s = record["per_sigma"][sigma]["scenarios"]
        y = np.arange(len(names))
        axes[0].plot([s[n].get("T1_rejection_rate", np.nan) for n in names], y, markers[k], mfc="white",
                     color="black", label=f"σ = {sigma}")
        axes[1].plot([s[n].get("b_mean", np.nan) for n in names], y, markers[k], mfc="white", color="black")
    axes[0].axvline(record["per_test_alpha"], ls=":", color="0.4")
    axes[0].set(yticks=np.arange(len(names)), yticklabels=names, xlabel="T1 rejection rate (fraction of replicates)")
    axes[1].axvline(0, ls=":", color="0.4")
    axes[1].set(xlabel="mean slope b̂ (dimensionless)")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "F5_validation.pdf")
    plt.close(fig)


def guarded_figures(cfg):
    data = primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL))
    data["ln_t"] = np.log(data["t"])
    fit = models.fit("M1", data, cfg["D6.primary_errors_in_lnR"])
    pulsars = data.groupby("pulsar").size().sort_values(ascending=False).index
    # F2
    n = len(pulsars)
    cols = 5
    fig, axes = plt.subplots(int(np.ceil(n / cols)), cols, figsize=(10, 2.1 * np.ceil(n / cols)), squeeze=False)
    for ax, p in zip(axes.flat, pulsars):
        d = data[data["pulsar"] == p]
        done = d["event"].astype(bool)
        ax.plot(d.loc[done, "ln_R"], d.loc[done, "ln_t"], "o", color="black", ms=3)
        ax.plot(d.loc[~done, "ln_R"], d.loc[~done, "ln_t"], "o", mfc="white", color="black", ms=3)
        x = np.linspace(d["ln_R"].min(), d["ln_R"].max(), 2)
        ax.plot(x, x + fit.estimates[f"c[{p}]"], "-", color="0.45")
        ax.set_title(p, fontsize=8)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.supxlabel("ln R (R = (Δν/ν)/(Δν̇/ν̇), dimensionless)")
    fig.supylabel("ln T (T in days; open: censored)")
    fig.tight_layout()
    fig.savefig(FIG / "F2_lnT_vs_lnR.pdf")
    plt.close(fig)
    # F3
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    for k, p in enumerate(pulsars[:4]):
        d = data[data["pulsar"] == p]
        times, surv = kaplan_meier(d["t"], d["event"])
        ax.step(np.r_[0, times], np.r_[1, surv], where="post", color=["0", "0.3", "0.55", "0.75"][k],
                ls=["-", "--", "-.", ":"][k], label=p)
    ax.set(xlabel="waiting time T (days)", ylabel="fraction still waiting", xscale="log")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "F3_kaplan_meier.pdf")
    plt.close(fig)
    # F4
    prim = json.loads((ROOT / "results" / "primary.json").read_text())
    sec = json.loads((ROOT / "results" / "secondary.json").read_text())
    ms = ["M0", "M2", "M3", "M4"]
    fig, axes = plt.subplots(1, 2, figsize=(7, 3))
    axes[0].bar(ms, [prim["aic"][m] - prim["aic"]["M0"] for m in ms], color="0.6", edgecolor="black")
    axes[0].set(ylabel="AIC − AIC(M0)")
    loo = ["M0", "M1", "M2", "M3", "M4"]
    axes[1].bar(loo, [sec["loo"][m] for m in loo], color="0.85", edgecolor="black")
    axes[1].set(ylabel="mean LOO log score (per interval)")
    fig.tight_layout()
    fig.savefig(FIG / "F4_model_comparison.pdf")
    plt.close(fig)
    # F6
    beta = pd.DataFrame(sec["implied_beta"]).sort_values("beta")
    fig, ax = plt.subplots(figsize=(4.5, 0.25 * len(beta) + 1))
    ax.errorbar(beta["beta"], np.arange(len(beta)), xerr=beta["se_beta"], fmt="o", color="black", mfc="white")
    ax.axvline(0, ls=":", color="0.4")
    ax.set(yticks=np.arange(len(beta)), yticklabels=beta["pulsar"], xlabel="implied β (dimensionless)")
    fig.tight_layout()
    fig.savefig(FIG / "F6_implied_beta.pdf")
    plt.close(fig)


def main():
    cfg = load_config()
    FIG.mkdir(exist_ok=True)
    f1_exclusion_flow()
    print("Wrote figures/F1_exclusion_flow.pdf")
    if (ROOT / "results" / "validation.json").exists():
        f5_validation()
        print("Wrote figures/F5_validation.pdf")
    else:
        print("Skipped F5: results/validation.json does not exist (run scripts/06_validate.py).")
    if unblinding_allowed():
        guarded_figures(cfg)
        print("Wrote figures/F2, F3, F4, F6")
    else:
        print("Skipped F2, F3, F4, F6: they use real waiting times and the guard does not pass (not frozen).")


if __name__ == "__main__":
    main()
