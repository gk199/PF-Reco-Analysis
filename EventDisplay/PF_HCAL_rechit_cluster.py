import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import argparse
###############################################################################
# HCAL Event Display
###############################################################################

class HCALEventDisplay:
    def __init__(self, data, start_event=0):
        self.data = data
        self.event = start_event
        self.cluster_index = 0
        event = self.event
        
        self.hcal_eta    = self.data['hcal_eta']
        self.hcal_phi    = self.data['hcal_phi']
        self.hcal_energy = self.data['hcal_energy']
        self.hcal_time   = self.data['hcal_time']
        self.hcal_depth  = self.data['hcal_depth']

        self.rh_eta      = self.data['rh_eta']
        self.rh_phi      = self.data['rh_phi']
        self.rh_energy   = self.data['rh_energy']
        self.rh_time     = self.data['rh_time']  
        #self.nClusters = len(self.hcal_eta[event])
        self.rh_depth    = self.data['rh_depth']
        self.rh_clusterIndex = self.data['rh_clusterIndex']

        ############################################################################
        # Set up figure and axes
        ############################################################################
        # Figure with 2x4 subplots (top row is energy, lower is timing)
        self.fig = plt.figure(figsize=(18, 10))
        self.ax = [self.fig.add_subplot(2, 4, i+1) for i in range(8)]
        plt.subplots_adjust(left=0.25, right=0.95)
        plt.subplots_adjust(wspace=0.2, hspace=0.2)

        # Initialize axes: titles, labels, aspect, limits
        for d in range(4):
            axE = self.ax[d]
            axT = self.ax[d + 4]
            for axx in [axE, axT]:
                axx.set_title(f"Depth {d+1}")
                axx.set_xlabel("$\eta$")
                axx.set_ylabel("$\phi$")
                axx.set_aspect("equal", adjustable="box")
                axx.set_xlim(-1.5, 1.5)
                axx.set_ylim(-1.5, 1.5)

        # Navigation buttons
        self.ax_prevEvent = self.fig.add_axes([0.05, 0.01, 0.12, 0.05])
        self.ax_nextEvent = self.fig.add_axes([0.18, 0.01, 0.12, 0.05])
        self.ax_prevClus  = self.fig.add_axes([0.70, 0.01, 0.12, 0.05])
        self.ax_nextClus  = self.fig.add_axes([0.83, 0.01, 0.12, 0.05])

        self.b_prevEvent = Button(self.ax_prevEvent, 'Prev Event')
        self.b_nextEvent = Button(self.ax_nextEvent, 'Next Event')
        self.b_prevClus  = Button(self.ax_prevClus, 'Prev Cluster')
        self.b_nextClus  = Button(self.ax_nextClus, 'Next Cluster')

        self.b_prevEvent.on_clicked(self.prev_event)
        self.b_nextEvent.on_clicked(self.next_event)
        self.b_prevClus.on_clicked(self.prev_cluster)
        self.b_nextClus.on_clicked(self.next_cluster)

        # Colorbars (initially None)
        self.cbar_energy = None
        self.cbar_time = None

        # Draw first display
        self.draw_display()
        plt.show()

    ############################################################################
    # Navigation
    ############################################################################

    def next_event(self,_event):
        self.event = (self.event + 1) % len(self.hcal_eta)
        self.cluster_index = 0
        self.draw_display()

    def prev_event(self,_event):
        self.event = (self.event - 1) % len(self.hcal_eta)
        self.cluster_index = 0
        self.draw_display()

    def next_cluster(self,_event):
        nClusters = len(self.hcal_eta[self.event])
        if nClusters == 0:
            print("No clusters in this event.")
            return
        self.cluster_index = (self.cluster_index + 1) % nClusters
        self.draw_display()

    def prev_cluster(self,_event):
        nClusters = len(self.hcal_eta[self.event])
        if nClusters == 0:
            print("No clusters in this event.")
            return
        self.cluster_index = (self.cluster_index - 1) % nClusters
        self.draw_display()

    ############################################################################
    # Draw Display
    ############################################################################

    def draw_display(self):
        nClusters = len(self.hcal_eta[self.event])

        # Add cluster summary text (top-right of figure)
        # Remove old annotation if it exists
        if hasattr(self, "cluster_summary"): 
            for t in self.cluster_summary: t.remove()
        self.cluster_summary = []

        # only plot clusters over 2 GeV 
        counter = 0
        energies = ak.to_numpy(self.hcal_energy[self.event])
        valid_clusters = np.where(energies > 2.0)[0]
        for i in valid_clusters:
            summary_lines = (
                f"{i:2d}: "
                f"E={float(self.hcal_energy[self.event][i]):5.2f}  "
                f"$\eta$={float(self.hcal_eta[self.event][i]):+5.2f}  "
                f"$\phi$={float(self.hcal_phi[self.event][i]):+5.2f}  "
                f"t={float(self.hcal_time[self.event][i]):+6.2f} ns"
            )

            # make current cluster bolded
            if i == self.cluster_index: weight = "bold"
            else: weight = "normal"

            dy = 0.018 # offset for each line
            # Add new annotation (figure coordinates)
            t = self.fig.text(
                0.03, 0.98 - counter * dy,        # x,y in figure coordinates, dy to move down a line each time
                summary_lines,
                ha="left", va="top",
                fontsize=10,
                family="monospace",
                fontweight=weight
            )
            self.cluster_summary.append(t)
            counter += 1


        # Remove old scatter plots and outlines
        for axx in self.ax:
            for coll in axx.collections[:]:
                coll.remove()
            for line in axx.lines[:]:
                line.remove()

        if nClusters == 0:
            # Empty event
            for axx in self.ax:
                axx.text(0.5, 0.5, "No clusters in this event",
                         ha='center', va='center', fontsize=14)
            self.fig.suptitle(f"Event {self.event}: No HCAL Clusters", fontsize=16)
            self.fig.canvas.draw_idle()
            return

        idx = self.cluster_index
        # ---------------------------------------------------------
        # Skip clusters with energy < 2 GeV
        # ---------------------------------------------------------
        energies = ak.to_numpy(self.hcal_energy[self.event])
        valid = np.where(energies > 2.0)[0]

        if len(valid) == 0:
            # No clusters above threshold in this event
            self.fig.suptitle(f"Event {self.event}: No clusters with E > 1 GeV", fontsize=16)
            self.fig.canvas.draw_idle()
            return

        # Snap cluster_index to nearest valid cluster
        if idx not in valid:
            # find closest valid index after current
            next_valid = valid[valid >= idx]
            if len(next_valid) == 0:
                idx = valid[0]
            else:
                idx = next_valid[0]
            self.cluster_index = idx

        if idx >= nClusters:
            idx = 0
            self.cluster_index = 0

        # Cluster center
        c_eta = float(self.hcal_eta[self.event][idx])
        c_phi = float(self.hcal_phi[self.event][idx])
        c_depth = float(self.hcal_depth[self.event][idx])

        # Select HBHE rechits for this cluster
        mask_cluster = (self.rh_clusterIndex[self.event] == idx)
        eta   = ak.to_numpy(self.rh_eta[self.event][mask_cluster])
        phi   = ak.to_numpy(self.rh_phi[self.event][mask_cluster])
        depth = ak.to_numpy(self.rh_depth[self.event][mask_cluster])
        energy = ak.to_numpy(self.rh_energy[self.event][mask_cluster])
        time   = ak.to_numpy(self.rh_time[self.event][mask_cluster])

        # Geometry cut: deltaEta/Phi < 0.4
        # should already be done based on rechits that match to the cluster but in case:
        deta = np.abs(eta - c_eta)
        dphi = np.abs((phi - c_phi + np.pi) % (2*np.pi) - np.pi)
        mask_geo = (deta < 0.4) & (dphi < 0.4)

        eta   = eta[mask_geo]
        phi   = phi[mask_geo]
        depth = depth[mask_geo]
        energy = energy[mask_geo]
        time   = time[mask_geo]

        # Mask rechits with invalid sentinel time (-999)
        valid_time = time > -100

        # Determine color scale
        vminE, vmaxE = (np.min(energy), np.max(energy)) if len(energy) > 0 else (0,1)
        valid_t_vals = time[valid_time]
        if len(valid_t_vals) > 0:
            vminT, vmaxT = np.min(valid_t_vals), np.max(valid_t_vals)
            if vmaxT - vminT < 1.0:  # degenerate range → center on cluster time ±0.5 ns
                mid = (vminT + vmaxT) / 2
                vminT, vmaxT = mid - 0.5, mid + 0.5
        else:
            vminT, vmaxT = -1.0, 1.0

        # Draw scatter and cluster outlines
        scE_sample = None
        scT_sample = None
        for d in [1,2,3,4]:
            axE = self.ax[d-1]
            axT = self.ax[d-1+4]

            hit = (depth == d)
            hit_invalid = hit & ~valid_time
            hit_t = hit & valid_time

            # Draw hits
            if np.sum(hit) > 0:
                sc = axE.scatter(eta[hit], phi[hit], s=80, c=energy[hit],
                                 cmap="viridis", vmin=vminE, vmax=vmaxE)
                scE_sample = sc
            else:
                axE.text(0.5, 0.5, "No hits", ha='center', va='center')

            if np.sum(hit) == 0:
                axT.text(0.5, 0.5, "No hits", ha='center', va='center')
            else:
                if np.sum(hit_invalid) > 0:
                    axT.scatter(eta[hit_invalid], phi[hit_invalid], s=80,
                                c="grey", marker='x', zorder=2)
                if np.sum(hit_t) > 0:
                    sc = axT.scatter(eta[hit_t], phi[hit_t], s=80, c=time[hit_t],
                                     cmap="plasma", vmin=vminT, vmax=vmaxT, zorder=3)
                    scT_sample = sc

            # Cluster outline
            if c_depth == d or (c_depth < d and c_depth > d-1) or (c_depth > d and c_depth < d+1): # draw cluster at the closest depth
                for axx in (axE, axT):
                    axx.plot(
                        [c_eta-0.25, c_eta+0.25, c_eta+0.25, c_eta-0.25, c_eta-0.25],
                        [c_phi-0.25, c_phi-0.25, c_phi+0.25, c_phi+0.25, c_phi-0.25],
                        'r--'
                    )

        for axx in self.ax:
            # Zoom to cluster +-0.27 in eta and phi, slightly bigger than dEta dPhi used
            axx.set_xlim(c_eta - 0.27, c_eta + 0.27)
            axx.set_ylim(c_phi - 0.27, c_phi + 0.27)

        # Add/update colorbars
        if scE_sample is None:
            scE_sample = self.ax[0].scatter([], [], c=[], cmap="viridis", vmin=vminE, vmax=vmaxE)
        if scT_sample is None:
            scT_sample = self.ax[4].scatter([], [], c=[], cmap="plasma", vmin=vminT, vmax=vmaxT)

        if self.cbar_energy is None:
            self.cbar_energy = self.fig.colorbar(scE_sample, ax=self.ax[0:4], location="right", fraction=0.02)
            self.cbar_energy.set_label("Energy [GeV]")
        else:
            self.cbar_energy.mappable.set_clim(vminE, vmaxE)

        if self.cbar_time is None:
            self.cbar_time = self.fig.colorbar(scT_sample, ax=self.ax[4:8], location="right", fraction=0.02)
            self.cbar_time.set_label("Time [ns]")
        elif len(valid_t_vals) > 0:
            self.cbar_time.mappable.set_clim(vminT, vmaxT)

        self.fig.suptitle(f"Event {self.event} — Cluster {idx}, cluster depth {c_depth}", fontsize=16)
        self.fig.canvas.draw_idle()


def main():

    ###############################################################################
    # Creates a parser for command line arguments
    ###############################################################################
    parser = argparse.ArgumentParser(description="HCAL Event Display")
    
    # Add a required positional argument
    parser.add_argument("filename", type=str, nargs="?", help="Path to the ROOT file from the ntupler",default="pfObjectsNtuple_tdc.root")
    args = parser.parse_args()

    ###############################################################################
    # Load ROOT file
    ###############################################################################

    # file = uproot.open("../../Downloads/pfObjectsNtuple.root") # path to your root file from the ntupler
    # directory = file["pfObjectsNtupler"] 
    file = uproot.open(args.filename) # path to your root file from the ntupler
    # directory = file["pfObjectsNtuplertdc"] 
    directory = file["pfObjectsNtupler"] 
    tree = directory["pfTree"]

    # HCAL clusters
    hcal_eta    = tree["hcal_eta"].array(library="ak")
    hcal_phi    = tree["hcal_phi"].array(library="ak")
    hcal_energy = tree["hcal_energy"].array(library="ak")
    hcal_time   = tree["hcal_time"].array(library="ak")
    hcal_depth  = tree["hcal_depth"].array(library="ak")

    # HBHE rechits
    rh_eta    = tree["hbhe_rechit_eta"].array(library="ak")
    rh_phi    = tree["hbhe_rechit_phi"].array(library="ak")
    rh_energy = tree["hbhe_rechit_energy"].array(library="ak")
    rh_time   = tree["hbhe_rechit_time"].array(library="ak")
    # rh_time   = tree["hbhe_rechit_tdc"].array(library="ak")
    rh_depth  = tree["hbhe_rechit_depth"].array(library="ak")
    # rh_clusterIndex = tree["hbheRechit_clusterIdx"].array(library="ak")
    rh_clusterIndex = tree["hbhe_rechit_clusterIndex"].array(library="ak")

    data = dict(
        hcal_eta=hcal_eta,
        hcal_phi=hcal_phi,
        hcal_energy=hcal_energy,
        hcal_time=hcal_time,
        hcal_depth=hcal_depth,
        rh_eta=rh_eta,
        rh_phi=rh_phi,
        rh_energy=rh_energy,
        rh_time=rh_time,
        rh_depth=rh_depth,
        rh_clusterIndex=rh_clusterIndex)
    ###############################################################################
    # Start Event Display
    ###############################################################################
    # Example: start at event 5
    viewer = HCALEventDisplay(data, start_event=20)

if __name__ == "__main__":
    main()
