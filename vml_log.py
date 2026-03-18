#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import argparse
import numpy as np
from matplotlib import pyplot as plt
import my_pyplot as mpt

WDIR = Path.cwd()

def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect and plot ML_LOGFILE diagnostics."
    )
    parser.add_argument(
        "--sys",
        "-s",
        dest="sys_list",
        nargs="+",
        default=["4000"],
        help="List of system directory names. Default: 4000",
    )
    parser.add_argument(
        "--idx",
        "-i",
        dest="i_list",
        nargs="+",
        default=[""],
        help='List of subdirectory names. Default: "" (i.e. use system directory itself)',
    )
    parser.add_argument(
        "--force-window",
        type=int,
        default=81,
        help="Moving-average window for BEEF and BEFPS plots. Default: 81",
    )
    parser.add_argument(
        "--sf-window",
        type=int,
        default=81,
        help="Moving-average window for SFF/SPFFPS plots. Default: 81",
    )
    parser.add_argument(
        "--force-skip",
        dest="force_skip_steps",
        type=int,
        default=20,
        help="Skip steps smaller than this value in BEEF and BEFPS plots. Default: 20",
    )
    parser.add_argument(
        "--err-skip",
        dest="err_skip_steps",
        type=int,
        default=20,
        help="Skip steps smaller than this value in ERR plots. Default: 20",
    )
    parser.add_argument(
        "--sf-skip",
        dest="sf_skip_steps",
        type=int,
        default=20,
        help="Skip steps smaller than this value in SFF/SPFFPS plots. Default: 20",
    )
    return parser.parse_args()

# -----------------------------
# Configurations to process
# -----------------------------
args = parse_args()
sys_list = args.sys_list
i_list = args.i_list

err_dir = WDIR / "_results/err"
beef_dir = WDIR / "_results/beef"
befps_dir = WDIR / "_results/befps"
sff_dir = WDIR / "_results/sff"
spffps_dir = WDIR / "_results/spffps"

for d in [err_dir, beef_dir, befps_dir, sff_dir, spffps_dir]:
    d.mkdir(exist_ok=True, parents=True)

force_window = args.force_window
sf_window = args.sf_window

force_skip_steps = args.force_skip_steps
err_skip_steps = args.err_skip_steps
sf_skip_steps = args.sf_skip_steps

pat_old_n = re.compile(r"^ML_LOGFILE\.old\.(\d+)$")


def make_param_name(sys_name, i_name):
    return sys_name if i_name == "" else f"{sys_name}-{i_name}"


def get_logfiles(workdir: Path):
    old_n = []
    old = None
    cur = None

    for f in workdir.iterdir():
        if not f.is_file():
            continue
        m = pat_old_n.match(f.name)
        if m:
            old_n.append((int(m.group(1)), f))
        elif f.name == "ML_LOGFILE.old":
            old = f
        elif f.name == "ML_LOGFILE":
            cur = f

    old_n.sort(reverse=True)  # ML_LOGFILE.old.N -> ... -> ML_LOGFILE.old.1
    files = [f for _, f in old_n]
    if old is not None:
        files.append(old)
    if cur is not None:
        files.append(cur)
    return files


def grep_lines(logfile: Path, key: str):
    """Return comment lines starting with '# key' and data lines starting with 'key'."""
    out = []
    key_head = f"# {key}"
    with logfile.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.lstrip()
            if s.startswith(key_head) or s.startswith(key):
                out.append(line.rstrip("\n"))
    return out


def split_lines(lines, key):
    header, data = [], []
    for line in lines:
        s = line.lstrip()
        if s.startswith("#"):
            header.append(line)
        elif s.startswith(key):
            data.append(line)
    return header, data


def shift_data_lines(data_lines, key, offset):
    out = []
    last_nstep = None

    for line in data_lines:
        parts = line.split()
        if parts[0] != key:
            continue

        raw_nstep = int(parts[1])

        # Skip restart/initial frame with nstep = 0 for all outputs
        if raw_nstep == 0:
            continue

        nstep = raw_nstep + offset
        parts[1] = str(nstep)
        out.append(" ".join(parts))
        last_nstep = nstep

    return out, last_nstep


def moving_average(y, window=81):
    if window < 1:
        raise ValueError("window must be >= 1")
    if window == 1:
        return y.copy()
    if window % 2 == 0:
        window += 1

    pad = window // 2
    ypad = np.pad(y, pad, mode="edge")
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(ypad, kernel, mode="valid")


def append_shifted_block(out_lines, logfile, key, offset, first_block):
    lines = grep_lines(logfile, key)
    if not lines:
        return first_block, None

    header, data = split_lines(lines, key)
    shifted_data, last_nstep = shift_data_lines(data, key, offset)

    out_lines.append(f"# ===== source: {logfile.name} =====")

    if first_block:
        out_lines.extend(header)
        first_block = False

    out_lines.extend(shifted_data)
    return first_block, last_nstep


def write_text_file(path: Path, lines):
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def collect_logs():
    for sys_name in sys_list:
        for i_name in i_list:
            workdir = WDIR / sys_name / i_name
            param_name = make_param_name(sys_name, i_name)

            outputs = {
                "BEEF": beef_dir / f"{param_name}.dat",
                "BEFPS": befps_dir / f"{param_name}.dat",
                "ERR": err_dir / f"{param_name}.dat",
                "SFF": sff_dir / f"{param_name}.dat",
                "SPFFPS": spffps_dir / f"{param_name}.dat",
            }

            if not workdir.is_dir():
                print(f"Skip missing directory: {workdir}")
                continue

            logfiles = get_logfiles(workdir)
            if not logfiles:
                for path in outputs.values():
                    write_text_file(path, [])
                print(f"No ML_LOGFILE found in: {workdir}")
                continue

            collected = {key: [] for key in outputs}
            first_flags = {key: True for key in outputs}
            prev_last = 0

            for logfile in logfiles:
                offset = prev_last

                first_flags["BEEF"], last_beef = append_shifted_block(
                    collected["BEEF"], logfile, "BEEF", offset, first_flags["BEEF"]
                )
                if last_beef is not None:
                    prev_last = last_beef

                for key in ["BEFPS", "ERR", "SFF", "SPFFPS"]:
                    first_flags[key], _ = append_shifted_block(
                        collected[key], logfile, key, offset, first_flags[key]
                    )

            for key, path in outputs.items():
                write_text_file(path, collected[key])

            print(
                f"Done: {sys_name}/{i_name}. "
                f"Log files processed: {len(logfiles)}. "
                f"BEEF lines: {len(collected['BEEF'])}. "
                f"BEFPS lines: {len(collected['BEFPS'])}. "
                f"ERR lines: {len(collected['ERR'])}. "
                f"SFF lines: {len(collected['SFF'])}. "
                f"SPFFPS lines: {len(collected['SPFFPS'])}."
            )


def get_params():
    return [make_param_name(sys_name, i_name) for sys_name in sys_list for i_name in i_list]


def apply_skip(x, *ys, skip_steps=0):
    x = np.atleast_1d(x)
    ys = [np.atleast_1d(y) for y in ys]

    if skip_steps > 0:
        mask = x >= skip_steps
        x = x[mask]
        ys = [y[mask] for y in ys]

    return (x, *ys)


def plot_raw_and_smooth(ax, x, y, window, label=None):
    raw = ax.plot(x, y, ls="--", alpha=0.2)[0]
    line = ax.plot(
        x,
        moving_average(y, window=window),
        ls="-",
        c=raw.get_color(),
        label=label,
    )[0]
    return line


def read_triplet_species_file(filepath: Path, key: str, value1_name: str, value2_name: str):
    data_out = {}

    with filepath.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            parts = s.split()
            if parts[0] != key:
                continue

            x = float(parts[1])

            for i in range(2, len(parts), 3):
                el = parts[i]
                v1 = float(parts[i + 1])
                v2 = float(parts[i + 2])

                if el not in data_out:
                    data_out[el] = {"x": [], value1_name: [], value2_name: []}

                data_out[el]["x"].append(x)
                data_out[el][value1_name].append(v1)
                data_out[el][value2_name].append(v2)

    for el in data_out:
        for k in data_out[el]:
            data_out[el][k] = np.asarray(data_out[el][k])

    return data_out


def read_sff(filepath: Path):
    data = np.loadtxt(filepath, comments="#", dtype=str)
    data = np.atleast_2d(data)

    mask = data[:, 0] == "SFF"
    data = data[mask]

    if data.size == 0:
        return {
            "x": np.array([]),
            "sfmax": np.array([]),
            "sfmin": np.array([]),
            "sfmean": np.array([]),
            "sfvar": np.array([]),
            "threshold": np.array([]),
        }

    return {
        "x": data[:, 1].astype(float),
        "sfmax": data[:, 2].astype(float),
        "sfmin": data[:, 3].astype(float),
        "sfmean": data[:, 4].astype(float),
        "sfvar": data[:, 5].astype(float),
        "threshold": data[:, 6].astype(float),
    }


def read_spffps(filepath: Path):
    return read_triplet_species_file(
        filepath, "SPFFPS", "sfmax", "sfmean"
    )


def read_befps(filepath: Path):
    return read_triplet_species_file(
        filepath, "BEFPS", "bee_max_force", "bee_ave_force"
    )


def init_species_panel_figure(n_panels):
    nrows, ncols = 2, 3
    max_panels = nrows * ncols

    if n_panels > max_panels:
        print(
            f"Warning: number of panels ({n_panels}) exceeds 2x3 layout. "
            f"Only the first {max_panels} panels will be shown."
        )
        n_panels = max_panels

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16.0, 10.0),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = axes.ravel()
    fig.subplots_adjust(wspace=0.03, hspace=0.03)
    return fig, axes, n_panels


def finalize_species_panel_figure(fig, axes, n_panels, legend_handles, legend_labels, save_name):
    for i in range(n_panels):
        mpt.minor(axes[i])

    for i in range(n_panels, len(axes)):
        axes[i].set_visible(False)

    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.95),
        ncol=min(len(legend_labels), 8),
        frameon=False,
    )

    mpt.savepdf(save_name)
    plt.show()


def plot_beef(params, smooth_window=81, skip_steps=0):
    fig, axes = mpt.init_plot_triple(wspace=0.25)

    for param in params:
        f = beef_dir / f"{param}.dat"
        if not f.exists() or f.stat().st_size == 0:
            continue

        try:
            x, e, fmax, thr = np.loadtxt(
                f,
                usecols=(1, 2, 3, 5),
                unpack=True,
                comments="#",
            )
        except Exception as exc:
            print(f"Failed to read {f}: {exc}")
            continue

        x, e, fmax, thr = apply_skip(x, e, fmax, thr, skip_steps=skip_steps)

        if x.size == 0:
            continue

        axes[0].plot(x, e * 1e3, label=param)
        plot_raw_and_smooth(axes[1], x, fmax, smooth_window, label=param)
        axes[2].plot(x, thr, label=param)

    axes[0].set_ylabel("Energy (meV/atom)")
    axes[1].set_ylabel("Max Force (eV/Å)")
    axes[2].set_ylabel("Threshold (eV/Å)")
    axes[0].legend(fancybox=False, edgecolor="black", ncol=2)

    for ax in axes:
        mpt.minor(ax)
        ax.set_xlabel("Time (fs)")

    mpt.savepdf("BEEF")
    plt.show()


def plot_err(params, skip_steps=0):
    fig, axes = mpt.init_plot_triple(wspace=0.25)

    for param in params:
        f = err_dir / f"{param}.dat"
        if not f.exists() or f.stat().st_size == 0:
            continue

        try:
            x, e, force, stress = np.loadtxt(
                f,
                usecols=(1, 2, 3, 4),
                unpack=True,
                comments="#",
            )
        except Exception as exc:
            print(f"Failed to read {f}: {exc}")
            continue

        x, e, force, stress = apply_skip(
            x, e, force, stress, skip_steps=skip_steps
        )

        if x.size == 0:
            continue

        axes[0].plot(x, e * 1e3, marker="s", label=param)
        axes[1].plot(x, force, marker="s", label=param)
        axes[2].plot(x, stress / 10.0, marker="s", label=param)

    axes[0].set_ylabel("Energy RMSE (meV/atom)")
    axes[1].set_ylabel("Force RMSE (eV/Å)")
    axes[2].set_ylabel("Stress RMSE (GPa)")
    axes[0].legend(fancybox=False, edgecolor="black", ncol=2)

    for ax in axes:
        mpt.minor(ax)
        ax.set_xlabel("Time (fs)")

    mpt.savepdf("RMSE")
    plt.show()


def plot_befps(params, smooth_window=81, skip_steps=0):
    all_elements = set()

    for param in params:
        f_befps = befps_dir / f"{param}.dat"
        if not f_befps.exists() or f_befps.stat().st_size == 0:
            continue
        befps = read_befps(f_befps)
        all_elements.update(befps.keys())

    all_elements = sorted(all_elements)
    n_panels = 1 + len(all_elements)

    fig, axes, n_panels = init_species_panel_figure(n_panels)
    all_elements = all_elements[: n_panels - 1]

    legend_handles = []
    legend_labels = []

    for param in params:
        f_beef = beef_dir / f"{param}.dat"
        f_befps = befps_dir / f"{param}.dat"

        if (not f_beef.exists() or f_beef.stat().st_size == 0) and \
           (not f_befps.exists() or f_befps.stat().st_size == 0):
            continue

        try:
            if f_beef.exists() and f_beef.stat().st_size > 0:
                x_beef, _, fmax_beef, _ = np.loadtxt(
                    f_beef,
                    usecols=(1, 2, 3, 5),
                    unpack=True,
                    comments="#",
                )
                x_beef, fmax_beef = apply_skip(
                    x_beef, fmax_beef, skip_steps=skip_steps
                )
            else:
                x_beef = np.array([])
                fmax_beef = np.array([])

            befps = (
                read_befps(f_befps)
                if f_befps.exists() and f_befps.stat().st_size > 0
                else {}
            )
        except Exception as exc:
            print(f"Failed to read {param}: {exc}")
            continue

        if x_beef.size > 0:
            line = plot_raw_and_smooth(
                axes[0], x_beef, fmax_beef, smooth_window, label=param
            )
            if param not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(param)

        for iel, el in enumerate(all_elements):
            ax = axes[iel + 1]
            if el not in befps:
                continue

            x = befps[el]["x"]
            y = befps[el]["bee_max_force"]
            x, y = apply_skip(x, y, skip_steps=skip_steps)

            if x.size == 0:
                continue

            plot_raw_and_smooth(ax, x, y, smooth_window, label=param)

    axes[0].set_ylabel("Max Force (eV/Å)")
    axes[3].set_ylabel("Max Force (eV/Å)")
    axes[0].text(
        0.03, 0.97, "BEEF",
        transform=axes[0].transAxes,
        ha="left", va="top", fontsize=20, fontweight="bold"
    )

    for iel, el in enumerate(all_elements):
        axes[iel + 1].text(
            0.03, 0.97, el,
            transform=axes[iel + 1].transAxes,
            ha="left", va="top", fontsize=20, fontweight="bold"
        )

    for ax in axes[3:]:
        ax.set_xlabel("Time (fs)")

    finalize_species_panel_figure(
        fig, axes, n_panels, legend_handles, legend_labels, "BEFPS"
    )


def plot_sff(params, smooth_window=81, skip_steps=0):
    all_elements = set()

    for param in params:
        f_spffps = spffps_dir / f"{param}.dat"
        if not f_spffps.exists() or f_spffps.stat().st_size == 0:
            continue
        spffps = read_spffps(f_spffps)
        all_elements.update(spffps.keys())

    all_elements = sorted(all_elements)
    n_panels = 1 + len(all_elements)

    fig, axes, n_panels = init_species_panel_figure(n_panels)
    all_elements = all_elements[: n_panels - 1]

    legend_handles = []
    legend_labels = []

    for param in params:
        f_sff = sff_dir / f"{param}.dat"
        f_spffps = spffps_dir / f"{param}.dat"

        if (not f_sff.exists() or f_sff.stat().st_size == 0) and \
           (not f_spffps.exists() or f_spffps.stat().st_size == 0):
            continue

        try:
            sff = (
                read_sff(f_sff)
                if f_sff.exists() and f_sff.stat().st_size > 0
                else {"x": np.array([]), "sfmax": np.array([])}
            )
            spffps = (
                read_spffps(f_spffps)
                if f_spffps.exists() and f_spffps.stat().st_size > 0
                else {}
            )
        except Exception as exc:
            print(f"Failed to read {param}: {exc}")
            continue

        if sff["x"].size > 0:
            x, y = apply_skip(sff["x"], sff["sfmax"], skip_steps=skip_steps)
            if x.size > 0:
                line = plot_raw_and_smooth(
                    axes[0], x, y, smooth_window, label=param
                )
                if param not in legend_labels:
                    legend_handles.append(line)
                    legend_labels.append(param)

        for iel, el in enumerate(all_elements):
            ax = axes[iel + 1]
            if el not in spffps:
                continue

            x = spffps[el]["x"]
            y = spffps[el]["sfmax"]
            x, y = apply_skip(x, y, skip_steps=skip_steps)

            if x.size == 0:
                continue

            plot_raw_and_smooth(ax, x, y, smooth_window, label=param)

    axes[0].set_ylabel("SFF max")
    axes[3].set_ylabel("SFF max")
    axes[0].text(
        0.03, 0.97, "SFF",
        transform=axes[0].transAxes,
        ha="left", va="top", fontsize=20, fontweight="bold"
    )

    for iel, el in enumerate(all_elements):
        axes[iel + 1].text(
            0.03, 0.97, el,
            transform=axes[iel + 1].transAxes,
            ha="left", va="top", fontsize=20, fontweight="bold"
        )

    for ax in axes[3:]:
        ax.set_xlabel("Time (fs)")

    finalize_species_panel_figure(
        fig, axes, n_panels, legend_handles, legend_labels, "SPFFPS"
    )


def main():
    collect_logs()
    params = get_params()
    plot_beef(params, smooth_window=force_window, skip_steps=force_skip_steps)
    plot_err(params, skip_steps=err_skip_steps)
    plot_befps(params, smooth_window=force_window, skip_steps=force_skip_steps)
    plot_sff(params, smooth_window=sf_window, skip_steps=sf_skip_steps)


if __name__ == "__main__":
    main()
