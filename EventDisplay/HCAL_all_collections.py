import os
import argparse
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

###############################################################################
# Constants
###############################################################################

# PF rechit energy thresholds by depth [GeV] (applied per rechit, HB or HE)
EMIN_HB_BY_DEPTH = {1: 0.6, 2: 0.4, 3: 0.5, 4: 0.5}
EMIN_HE_BY_DEPTH = {1: 0.2, 2: 0.3, 3: 0.3, 4: 0.3, 5: 0.3, 6: 0.3, 7: 0.3}

# Special values that carry no timing measurement -> drawn in gray
TDC_INVALID_HB = [3]
TDC_INVALID_HE = [62, 63]
TIME_INVALID = -999

# Colorbar ranges for TDC codes once the special codes are removed
TDC_MAX_HB = 2
TDC_MAX_HE = 61

###############################################################################
# HCAL Event Display
###############################################################################

class HCALEventDisplay:
    def __init__(self, data, start_event=0):
        self.data = data
        self.event = start_event
        self.cluster_index = 0

        # Display configuration
        self.rechit_type    = self.data.get("rechit_type", "pfrh")
        self.time_mode      = self.data.get("time_mode", "tdc")
        self.pf_energy_cuts = self.data.get("pf_energy_cuts", True)
        self.cluster_emin   = self.data.get("cluster_emin", None)
        # TDC logic only applies to HBHE rechits read in TDC mode
        self.use_tdc = (self.rechit_type == "hbhe") and (self.time_mode == "tdc")

        self.hcal_eta    = self.data["hcal_eta"]
        self.hcal_phi    = self.data["hcal_phi"]
        self.hcal_energy = self.data["hcal_energy"]
        self.hcal_depth  = self.data["hcal_depth"]

        self.rh_eta      = self.data["rh_eta"]
        self.rh_phi      = self.data["rh_phi"]
        self.rh_energy   = self.data["rh_energy"]
        self.rh_time     = self.data["rh_time"]
        self.rh_depth    = self.data["rh_depth"]
        self.rh_clusterIndex = self.data["rh_clusterIndex"]
        self.rh_isHE     = self.data["rh_isHE"]

        ############################################################################
        # Set up figure (panels are built in build_axes, 4 or 7 depths)
        ############################################################################
        self.fig = plt.figure(figsize=(20, 10))

        self.time_norm = colors.Normalize(vmin=0, vmax=TDC_MAX_HE)
        self.time_sm = cm.ScalarMappable(norm=self.time_norm, cmap="plasma")
        self.time_sm.set_array([])

        self.energy_norm = LogNorm(vmin=1e-3, vmax=1.0)
        self.energy_sm = cm.ScalarMappable(norm=self.energy_norm, cmap="viridis")
        self.energy_sm.set_array([])

        self.ax = []
        self.n_depth = None
        self.cbar_energy = None
        self.cbar_time = None
        self.gray_legend = None
        self.cluster_summary = []

        # Navigation buttons
        self.ax_prevEvent = self.fig.add_axes([0.05, 0.01, 0.12, 0.05])
        self.ax_nextEvent = self.fig.add_axes([0.18, 0.01, 0.12, 0.05])
        self.ax_prevClus  = self.fig.add_axes([0.70, 0.01, 0.12, 0.05])
        self.ax_nextClus  = self.fig.add_axes([0.83, 0.01, 0.12, 0.05])

        self.b_prevEvent = Button(self.ax_prevEvent, "Prev Event")
        self.b_nextEvent = Button(self.ax_nextEvent, "Next Event")
        self.b_prevClus  = Button(self.ax_prevClus, "Prev Cluster")
        self.b_nextClus  = Button(self.ax_nextClus, "Next Cluster")

        self.b_prevEvent.on_clicked(self.prev_event)
        self.b_nextEvent.on_clicked(self.next_event)
        self.b_prevClus.on_clicked(self.prev_cluster)
        self.b_nextClus.on_clicked(self.next_cluster)

        self.build_axes(4)
        self.fig.canvas.mpl_connect("resize_event", self.on_resize)
        self.draw_display()
        plt.show()

    ############################################################################
    # Panel layout: 2 x n_depth (top row energy, bottom row time)
    ############################################################################

    def build_axes(self, n_depth):
        # Remove the previous layout (panels + colorbars)
        for cb in (self.cbar_energy, self.cbar_time):
            if cb is not None:
                cb.remove()
        for axx in self.ax:
            axx.remove()

        # Axes are placed by position_axes(); all panels share the same zoom
        # window, so only the first column carries phi tick labels
        self.ax = [self.fig.add_axes([0, 0, 1, 1]) for _ in range(2 * n_depth)]
        self.n_depth = n_depth

        for d in range(n_depth):
            for axx in (self.ax[d], self.ax[d + n_depth]):
                axx.set_title(f"Depth {d + 1}", fontsize=10)
                axx.set_xlabel(r"$\eta$")
                if d == 0:
                    axx.set_ylabel(r"$\phi$")
                else:
                    axx.tick_params(labelleft=False)
                axx.set_aspect("equal", adjustable="box")
                axx.xaxis.set_major_locator(MaxNLocator(3))
                axx.yaxis.set_major_locator(MaxNLocator(3))
                axx.tick_params(labelsize=8)

        self.cbar_energy = self.fig.colorbar(self.energy_sm, cax=self.fig.add_axes([0, 0, 1, 1]))
        self.cbar_energy.set_label("Energy [GeV] (log)")
        self.cbar_time = self.fig.colorbar(self.time_sm, cax=self.fig.add_axes([0, 0, 1, 1]))

        self.position_axes()

    # Spacings in inches, so the layout holds when the window is resized
    LAYOUT_IN = dict(
        summary=3.0,   # cluster list on the left
        ylabel=0.55,   # phi tick labels + label of the first column
        right=1.4,     # colorbar + its labels
        top=0.95,      # suptitle + panel titles
        row_gap=0.85,  # x label of the energy row + titles of the time row
        col_gap=0.35,
        margin=0.25,   # left margin when the cluster list sits above the panels
        list_min=1.8,  # minimum height kept for the cluster list above the panels
        line=0.18,     # cluster list line spacing
        cbar_pad=0.25,
        cbar_w=0.18,
    )

    def position_axes(self):
        """Square panels, as large as the window allows, with fixed gaps."""
        W, H = self.fig.get_size_inches()
        L = self.LAYOUT_IN
        n = self.n_depth

        x_right = W - L["right"]
        y_bottom = 0.07 * H + 0.65          # above the buttons + x labels
        y_top = H - L["top"]

        # 4 depths: cluster list beside the panels.
        # 7 depths: cluster list above the panels, so they can use the full width.
        self.list_above = n == 7
        if self.list_above:
            x_left = L["margin"] + L["ylabel"]
            y_top -= L["list_min"]
        else:
            x_left = L["summary"] + L["ylabel"]

        side = min((x_right - x_left - (n - 1) * L["col_gap"]) / n,
                   (y_top - y_bottom - L["row_gap"]) / 2)
        side = max(side, 0.3)

        block_w = n * side + (n - 1) * L["col_gap"]
        block_h = 2 * side + L["row_gap"]
        x0 = x_left + max(0.0, (x_right - x_left - block_w) / 2)
        # List above: panels sit in the lower part (40% of the spare height
        # below them) so the list still gets most of the remaining height
        slack = max(0.0, y_top - y_bottom - block_h)
        y0 = y_bottom + (0.4 * slack if self.list_above else slack / 2)

        for row in range(2):
            y = y0 + (1 - row) * (side + L["row_gap"])
            for col in range(n):
                x = x0 + col * (side + L["col_gap"])
                self.ax[row * n + col].set_position([x / W, y / H, side / W, side / H])
            cbar = self.cbar_energy if row == 0 else self.cbar_time
            x_cbar = x0 + block_w + L["cbar_pad"]
            cbar.ax.set_position([x_cbar / W, y / H, L["cbar_w"] / W, side / H])

        # Lowest point the cluster list may reach (inches)
        self.summary_ymin = (y0 + block_h + 0.40) if self.list_above else 0.09 * H

    def on_resize(self, _event):
        self.position_axes()
        self.draw_summary(self.valid_clusters())
        self.fig.canvas.draw_idle()

    def clear_axes(self):
        for axx in self.ax:
            for coll in axx.collections[:]:
                coll.remove()
            for line in axx.lines[:]:
                line.remove()
            for txt in axx.texts[:]:
                txt.remove()

    ############################################################################
    # Navigation
    ############################################################################

    def valid_clusters(self):
        # Clusters must have E > 0; --cluster-emin adds E > emin on top
        energies = ak.to_numpy(self.hcal_energy[self.event])
        mask = energies > 0.0
        if self.cluster_emin is not None:
            mask &= energies > self.cluster_emin
        return np.where(mask)[0]

    def next_event(self, _event):
        self.event = (self.event + 1) % len(self.hcal_eta)
        self.cluster_index = 0
        self.draw_display()

    def prev_event(self, _event):
        self.event = (self.event - 1) % len(self.hcal_eta)
        self.cluster_index = 0
        self.draw_display()

    def step_cluster(self, step):
        # Step through clusters that pass the cluster energy cut only
        valid = self.valid_clusters()
        if len(valid) == 0:
            print("No clusters passing the cluster energy cut in this event.")
            return
        pos = np.searchsorted(valid, self.cluster_index)
        if pos < len(valid) and valid[pos] == self.cluster_index:
            pos = (pos + step) % len(valid)
        else:
            pos = pos % len(valid) if step > 0 else (pos - 1) % len(valid)
        self.cluster_index = int(valid[pos])
        self.draw_display()

    def next_cluster(self, _event):
        self.step_cluster(+1)

    def prev_cluster(self, _event):
        self.step_cluster(-1)

    ############################################################################
    # Labels
    ############################################################################

    def collection_label(self):
        if self.rechit_type == "pfrh":
            return "PFRecHits"
        return "HBHE rechits, " + ("TDC" if self.use_tdc else "time")

    def cluster_cut_label(self):
        if self.cluster_emin is None:
            return "E > 0 GeV"
        return f"E > {self.cluster_emin:g} GeV"

    ############################################################################
    # Draw Display
    ############################################################################

    def draw_summary(self, valid):
        for t in self.cluster_summary:
            t.remove()
        self.cluster_summary = []

        # Show a window of clusters around the current one if the list is long
        H = self.fig.get_size_inches()[1]
        dy, y0 = self.LAYOUT_IN["line"] / H, 0.98
        y_min = self.summary_ymin / H
        max_lines = max(int((y0 - y_min) / dy), 3)
        rows = list(valid)
        above = below = 0
        if len(rows) > max_lines:
            n_show = max_lines - 2
            cur = int(np.searchsorted(valid, self.cluster_index))
            start = min(max(cur - n_show // 2, 0), len(rows) - n_show)
            above, below = start, len(rows) - (start + n_show)
            rows = rows[start:start + n_show]

        lines = []
        if above:
            lines.append((f"  ... {above} more above", "normal"))
        for i in rows:
            lines.append((
                f"{i:2d}: "
                f"E={float(self.hcal_energy[self.event][i]):7.2f}  "
                f"$\\eta$={float(self.hcal_eta[self.event][i]):+5.2f}  "
                f"$\\phi$={float(self.hcal_phi[self.event][i]):+5.2f}",
                "bold" if i == self.cluster_index else "normal"))
        if below:
            lines.append((f"  ... {below} more below", "normal"))

        for counter, (text, weight) in enumerate(lines):
            t = self.fig.text(0.01, y0 - counter * dy, text, ha="left", va="top",
                              fontsize=9, family="monospace", fontweight=weight)
            self.cluster_summary.append(t)

    def draw_empty(self, message):
        self.clear_axes()
        for axx in self.ax:
            axx.text(0.5, 0.5, "No clusters", ha="center", va="center",
                     fontsize=10, transform=axx.transAxes)
        self.fig.suptitle(f"Event {self.event}: {message} — {self.collection_label()}", fontsize=16)
        self.fig.canvas.draw_idle()

    def draw_display(self):
        ev = self.event
        nClusters = len(self.hcal_eta[ev])
        valid = self.valid_clusters()

        if len(valid) > 0 and self.cluster_index not in valid:
            nxt = valid[valid >= self.cluster_index]
            self.cluster_index = int(nxt[0]) if len(nxt) else int(valid[0])

        self.draw_summary(valid)

        if self.gray_legend is not None:
            self.gray_legend.remove()
            self.gray_legend = None

        if nClusters == 0:
            self.draw_empty("No HCAL clusters")
            return
        if len(valid) == 0:
            self.draw_empty(f"No clusters with {self.cluster_cut_label()}")
            return

        idx = self.cluster_index

        ###########################################################################
        # Cluster center
        ###########################################################################

        c_eta = float(self.hcal_eta[ev][idx])
        c_phi = float(self.hcal_phi[ev][idx])
        c_depth = float(self.hcal_depth[ev][idx])

        ###########################################################################
        # Select rechits matched to this cluster
        ###########################################################################

        mask_cluster = self.rh_clusterIndex[ev] == idx

        eta    = ak.to_numpy(self.rh_eta[ev][mask_cluster])
        phi    = ak.to_numpy(self.rh_phi[ev][mask_cluster])
        depth  = ak.to_numpy(self.rh_depth[ev][mask_cluster])
        energy = ak.to_numpy(self.rh_energy[ev][mask_cluster])
        time   = ak.to_numpy(self.rh_time[ev][mask_cluster])
        is_HE  = ak.to_numpy(self.rh_isHE[ev][mask_cluster])

        # Layout: 7 depths if the cluster has any HE rechit (before cuts), else 4
        cluster_is_HE = bool(np.any(is_HE))
        n_depth = 7 if cluster_is_HE else 4
        if n_depth != self.n_depth:
            self.build_axes(n_depth)
            self.draw_summary(valid)
        else:
            self.clear_axes()

        ###########################################################################
        # Rechit cuts (HB/HE decided per rechit)
        ###########################################################################

        keep = np.isfinite(energy) & (energy > 0) & np.isfinite(time)

        if self.use_tdc:
            # Negative TDC codes are dropped
            keep &= time >= 0

        if self.rechit_type == "hbhe" and self.pf_energy_cuts:
            for d, emin in EMIN_HB_BY_DEPTH.items():
                keep &= ~((~is_HE) & (depth == d) & (energy < emin))
            for d, emin in EMIN_HE_BY_DEPTH.items():
                keep &= ~((is_HE) & (depth == d) & (energy < emin))

        ###########################################################################
        # Geometry cut: deltaEta / deltaPhi < 0.4
        ###########################################################################

        deta = eta - c_eta
        dphi = (phi - c_phi + np.pi) % (2 * np.pi) - np.pi
        keep &= (np.abs(deta) < 0.4) & (np.abs(dphi) < 0.4)
        # Draw phi relative to the cluster so hits across phi = +-pi stay in view
        phi = c_phi + dphi

        eta, phi, depth = eta[keep], phi[keep], depth[keep]
        energy, time, is_HE = energy[keep], time[keep], is_HE[keep]

        ###########################################################################
        # Rechits without a timing measurement -> gray, off the colorbar
        ###########################################################################

        if self.use_tdc:
            mask_bad_time = ((~is_HE) & np.isin(time, TDC_INVALID_HB)) | \
                            ((is_HE) & np.isin(time, TDC_INVALID_HE))
            gray_label = "no TDC measurement (HB = 3, HE = 62/63)"
        else:
            mask_bad_time = time == TIME_INVALID
            gray_label = f"time = {TIME_INVALID}"
        mask_good_time = ~mask_bad_time

        ###########################################################################
        # Energy color scale (log)
        ###########################################################################

        if len(energy) > 0:
            vminE, vmaxE = np.min(energy), np.max(energy)
        else:
            vminE, vmaxE = 1e-3, 1.0
        if vminE == vmaxE:
            vmaxE = vminE * 1.01

        self.energy_norm.vmin = vminE
        self.energy_norm.vmax = vmaxE

        ###########################################################################
        # Time color scale
        ###########################################################################

        if self.use_tdc:
            # HE scale if any displayed rechit is in HE (pure HE or mixed HB+HE),
            # HB scale only for pure HB clusters
            use_he_scale = bool(np.any(is_HE)) if len(is_HE) > 0 else cluster_is_HE
            vminT = 0
            vmaxT = TDC_MAX_HE if use_he_scale else TDC_MAX_HB
            time_label = f"Time [TDC codes, {'HE' if use_he_scale else 'HB'} scale]"
        else:
            good_times = time[mask_good_time]
            if len(good_times) > 0:
                vminT, vmaxT = np.min(good_times), np.max(good_times)
                if vminT == vmaxT:
                    vminT, vmaxT = vminT - 0.5, vmaxT + 0.5
            else:
                vminT, vmaxT = -1.0, 1.0
            time_label = f"Time [ns], excluding {TIME_INVALID}"

        self.time_norm.vmin = vminT
        self.time_norm.vmax = vmaxT

        ###########################################################################
        # Draw scatter plots and cluster outlines
        ###########################################################################

        # Zoom: square window around the cluster, at least +-0.27 and large
        # enough that every selected rechit (up to +-0.4) is fully in view
        half = 0.27
        if len(eta) > 0:
            half = max(half, np.max(np.abs(eta - c_eta)), np.max(np.abs(phi - c_phi)))
        half += 0.06
        # Red box keeps the original box/zoom ratio (0.25 / 0.27)
        box = half * 0.25 / 0.27

        for d in range(1, n_depth + 1):
            axE = self.ax[d - 1]
            axT = self.ax[d - 1 + n_depth]

            hit = depth == d
            hit_good_time = hit & mask_good_time
            hit_bad_time  = hit & mask_bad_time

            if np.sum(hit) > 0:
                # Energy plot: all rechits
                axE.scatter(eta[hit], phi[hit], s=80, c=energy[hit],
                            cmap="viridis", norm=self.energy_norm)

                # Time plot: valid times on the plasma color scale
                if np.sum(hit_good_time) > 0:
                    axT.scatter(eta[hit_good_time], phi[hit_good_time], s=80,
                                c=time[hit_good_time], cmap="plasma", norm=self.time_norm)

                # Time plot: no timing measurement -> gray
                if np.sum(hit_bad_time) > 0:
                    axT.scatter(eta[hit_bad_time], phi[hit_bad_time], s=80,
                                color="gray", edgecolors="black", linewidths=0.5)
            else:
                for axx in (axE, axT):
                    axx.text(0.5, 0.5, "No hits", ha="center", va="center",
                             fontsize=9, transform=axx.transAxes)

            # Cluster outline at the closest depth(s)
            if abs(c_depth - d) < 1:
                for axx in (axE, axT):
                    axx.plot(
                        [c_eta - box, c_eta + box, c_eta + box, c_eta - box, c_eta - box],
                        [c_phi - box, c_phi - box, c_phi + box, c_phi + box, c_phi - box],
                        "r--",
                    )

        for axx in self.ax:
            axx.set_xlim(c_eta - half, c_eta + half)
            axx.set_ylim(c_phi - half, c_phi + half)

        ###########################################################################
        # Colorbars and gray-point legend
        ###########################################################################

        self.cbar_energy.update_normal(self.energy_sm)
        self.cbar_time.update_normal(self.time_sm)
        self.cbar_time.set_label(time_label)
        if self.use_tdc:
            # TDC codes are integers
            self.cbar_time.ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        if np.any(mask_bad_time):
            handle = Line2D([], [], marker="o", linestyle="", markersize=8,
                            markerfacecolor="gray", markeredgecolor="black", label=gray_label)
            self.gray_legend = self.fig.legend(handles=[handle], loc="lower center",
                                               bbox_to_anchor=(0.5, 0.01), fontsize=10)

        self.fig.suptitle(
            f"Event {ev} — Cluster {idx} ({'HE' if cluster_is_HE else 'HB'}), "
            f"cluster depth {c_depth:.2f} — {self.collection_label()}",
            fontsize=16,
        )
        self.fig.canvas.draw_idle()


###############################################################################
# Helpers
###############################################################################

def parse_cluster_emin(value):
    if value.strip().upper() == "NONE":
        return None
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a number in GeV or NONE, got '{value}'")


def classify_HE(ieta, depth, eta):
    """Per-rechit HE flag. |ieta| 16 is shared: HB owns depths 1-3, HE owns depth 4."""
    if ieta is not None:
        ieta = np.abs(ieta)
        return (ieta > 16) | ((ieta == 16) & (depth >= 4))
    # Fallback when no ieta branch is stored: same boundary expressed in |eta|
    aeta = np.abs(eta)
    return (aeta > 1.392) | ((aeta > 1.305) & (depth >= 4))


def check_split(tree, prefix, is_HE):
    """Cross-check the HB/HE split against the separate hb_*/he_* branches, if present."""
    hb, he = f"hb_{prefix}_eta", f"he_{prefix}_eta"
    if hb not in tree.keys() or he not in tree.keys():
        print("  HB/HE split: no hb_*/he_* branches to cross-check against")
        return
    split = tree.arrays([hb, he], library="ak")
    ok = ak.all(ak.sum(~is_HE, axis=1) == ak.num(split[hb])) and \
         ak.all(ak.sum(is_HE, axis=1) == ak.num(split[he]))
    if ok:
        print(f"  HB/HE split: matches hb_{prefix}_* / he_{prefix}_* in every event")
    else:
        print(f"  WARNING: HB/HE split does not match hb_{prefix}_* / he_{prefix}_* counts")


###############################################################################
# Main
###############################################################################

def main():
    ###############################################################################
    # Command line arguments
    ###############################################################################

    parser = argparse.ArgumentParser(description="HCAL Event Display")

    parser.add_argument(
        "filename",
        type=str,
        help="ROOT file from the ntupler (relative to --input-dir unless it is an absolute path)",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        default=os.getcwd(),
        help="Directory containing the ROOT file (default: the current working directory)",
    )

    parser.add_argument(
        "--tree-dir",
        type=str,
        default="pfObjectsNtupler",
        help="Directory inside the ROOT file containing pfTree",
    )

    parser.add_argument(
        "--rechit-type",
        choices=["hbhe", "pfrh"],
        default="pfrh",
        help="Which rechit collection to display: 'pfrh' for PFRecHits or 'hbhe' for HBHE rechits",
    )

    parser.add_argument(
        "--hbhe-time",
        choices=["tdc", "time"],
        default=None,
        help="HBHE rechits only: show raw TDC codes ('tdc', default) or reconstructed time in ns ('time')",
    )

    parser.add_argument(
        "--pf-energy-cuts",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="HBHE rechits only: apply the depth-dependent PF energy thresholds (default: on)",
    )

    parser.add_argument(
        "--cluster-emin",
        type=parse_cluster_emin,
        default=None,
        metavar="GEV|NONE",
        help="Only show clusters with E > this value in GeV (default: NONE, i.e. every cluster with E > 0)",
    )

    parser.add_argument(
        "--start-event",
        type=int,
        default=5,
        help="Event index to start from",
    )

    args = parser.parse_args()

    ###############################################################################
    # Resolve options that only apply to one collection
    ###############################################################################

    if args.rechit_type == "pfrh":
        if args.pf_energy_cuts is not None:
            print("NOTE: currently using pfrh collections, PF cuts on rechits already applied "
                  "(--pf-energy-cuts / --no-pf-energy-cuts is ignored)")
        if args.hbhe_time is not None:
            print("NOTE: --hbhe-time only applies to --rechit-type hbhe; "
                  "PFRecHits always use hbhe_pfrh_time [ns]")
        time_mode = "time"
        pf_energy_cuts = False
    else:
        time_mode = args.hbhe_time or "tdc"
        pf_energy_cuts = True if args.pf_energy_cuts is None else args.pf_energy_cuts

    ###############################################################################
    # Load ROOT file
    ###############################################################################

    path = os.path.join(args.input_dir, args.filename)
    file = uproot.open(path)
    directory = file[args.tree_dir]
    tree = directory["pfTree"]

    ###############################################################################
    # HCAL clusters
    ###############################################################################

    hcal_eta    = tree["hcal_eta"].array(library="ak")
    hcal_phi    = tree["hcal_phi"].array(library="ak")
    hcal_energy = tree["hcal_energy"].array(library="ak")
    hcal_depth  = tree["hcal_depth"].array(library="ak")

    ###############################################################################
    # Choose rechit branches
    ###############################################################################

    if args.rechit_type == "hbhe":
        prefix                 = "rechit"
        rh_eta_branch          = "hbhe_rechit_eta"
        rh_phi_branch          = "hbhe_rechit_phi"
        rh_energy_branch       = "hbhe_rechit_energy"
        rh_time_branch         = "hbhe_rechit_tdc" if time_mode == "tdc" else "hbhe_rechit_time"
        rh_depth_branch        = "hbhe_rechit_depth"
        rh_clusterIndex_branch = "hbheRechit_clusterIdx"
        rh_ieta_branch         = "hbhe_rechit_ietaAbs"
    else:
        prefix                 = "pfrh"
        rh_eta_branch          = "hbhe_pfrh_eta"
        rh_phi_branch          = "hbhe_pfrh_phi"
        rh_energy_branch       = "hbhe_pfrh_energy"
        rh_time_branch         = "hbhe_pfrh_time"
        rh_depth_branch        = "hbhe_pfrh_depth"
        rh_clusterIndex_branch = "hbhe_pfrh_clusterIdx"
        rh_ieta_branch         = "hbhe_pfrh_ieta"

    print("File:", path)
    print("Using rechit type:", args.rechit_type)
    if args.rechit_type == "hbhe":
        print("  time          :", time_mode)
        print("  PF energy cuts:", "on" if pf_energy_cuts else "off")
    print("  cluster cut   :", "E > 0 GeV" if args.cluster_emin is None else f"E > {args.cluster_emin:g} GeV")
    print("Rechit branches:")
    print("  eta          :", rh_eta_branch)
    print("  phi          :", rh_phi_branch)
    print("  energy       :", rh_energy_branch)
    print("  time         :", rh_time_branch)
    print("  depth        :", rh_depth_branch)
    print("  cluster index:", rh_clusterIndex_branch)

    ###############################################################################
    # Load rechits
    ###############################################################################

    rh_eta    = tree[rh_eta_branch].array(library="ak")
    rh_phi    = tree[rh_phi_branch].array(library="ak")
    rh_energy = tree[rh_energy_branch].array(library="ak")
    rh_time   = tree[rh_time_branch].array(library="ak")
    rh_depth  = tree[rh_depth_branch].array(library="ak")
    rh_clusterIndex = tree[rh_clusterIndex_branch].array(library="ak")

    # Per-rechit HB/HE flag
    if rh_ieta_branch in tree.keys():
        rh_ieta = tree[rh_ieta_branch].array(library="ak")
    else:
        print(f"  WARNING: {rh_ieta_branch} not found, deciding HB/HE from |eta| and depth")
        rh_ieta = None
    rh_isHE = classify_HE(rh_ieta, rh_depth, rh_eta)
    check_split(tree, prefix, rh_isHE)

    data = dict(
        hcal_eta=hcal_eta,
        hcal_phi=hcal_phi,
        hcal_energy=hcal_energy,
        hcal_depth=hcal_depth,
        rh_eta=rh_eta,
        rh_phi=rh_phi,
        rh_energy=rh_energy,
        rh_time=rh_time,
        rh_depth=rh_depth,
        rh_clusterIndex=rh_clusterIndex,
        rh_isHE=rh_isHE,
        rechit_type=args.rechit_type,
        time_mode=time_mode,
        pf_energy_cuts=pf_energy_cuts,
        cluster_emin=args.cluster_emin,
    )

    ###############################################################################
    # Start Event Display
    ###############################################################################

    start_event = args.start_event % len(hcal_eta)
    viewer = HCALEventDisplay(data, start_event=start_event)


if __name__ == "__main__":
    main()