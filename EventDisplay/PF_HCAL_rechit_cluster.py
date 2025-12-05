import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

###############################################################################
# Load ROOT file
###############################################################################

file = uproot.open("../../Downloads/pfObjectsNtuple.root") # path to your root file from the ntupler
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
rh_depth  = tree["hbhe_rechit_depth"].array(library="ak")
rh_clusterIndex = tree["hbhe_rechit_clusterIndex"].array(library="ak")

###############################################################################
# HCAL Event Display
###############################################################################

class HCALEventDisplay:
    def __init__(self, start_event=0):
        self.event = start_event
        self.cluster_index = 0

        # Figure with 2x4 subplots (top row is energy, lower is timing)
        self.fig = plt.figure(figsize=(14, 10))
        self.ax = [self.fig.add_subplot(2, 4, i+1) for i in range(8)]
        plt.subplots_adjust(wspace=0.25, hspace=0.25)

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

    def next_event(self, event):
        self.event = (self.event + 1) % len(hcal_eta)
        self.cluster_index = 0
        self.draw_display()

    def prev_event(self, event):
        self.event = (self.event - 1) % len(hcal_eta)
        self.cluster_index = 0
        self.draw_display()

    def next_cluster(self, event):
        nClusters = len(hcal_eta[self.event])
        if nClusters == 0:
            print("No clusters in this event.")
            return
        self.cluster_index = (self.cluster_index + 1) % nClusters
        self.draw_display()

    def prev_cluster(self, event):
        nClusters = len(hcal_eta[self.event])
        if nClusters == 0:
            print("No clusters in this event.")
            return
        self.cluster_index = (self.cluster_index - 1) % nClusters
        self.draw_display()

    ############################################################################
    # Draw Display
    ############################################################################

    def draw_display(self):
        event = self.event
        nClusters = len(hcal_eta[event])

        # Add cluster summary text (top-right of figure)
        # Remove old annotation if it exists
        if hasattr(self, "cluster_summary"): 
            for t in self.cluster_summary: t.remove()
        self.cluster_summary = []

        for i in range(nClusters):
            summary_lines = (
                f"{i:2d}: "
                f"E={float(hcal_energy[event][i]):5.2f}  "
                f"$\eta$={float(hcal_eta[event][i]):+5.2f}  "
                f"$\phi$={float(hcal_phi[event][i]):+5.2f}"
            )

            # make current cluster bolded
            if i == self.cluster_index: weight = "bold"
            else: weight = "normal"

            dy = 0.018 # offset for each line
            # Add new annotation (figure coordinates)
            t = self.fig.text(
                0.03, 0.98 - i * dy,        # x,y in figure coordinates, dy to move down a line each time
                summary_lines,
                ha="left", va="top",
                fontsize=10,
                family="monospace",
                fontweight=weight
            )
            self.cluster_summary.append(t)


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
            self.fig.suptitle(f"Event {event}: No HCAL Clusters", fontsize=16)
            self.fig.canvas.draw_idle()
            return

        idx = self.cluster_index
        if idx >= nClusters:
            idx = 0
            self.cluster_index = 0

        # Cluster center
        c_eta = float(hcal_eta[event][idx])
        c_phi = float(hcal_phi[event][idx])
        c_depth = float(hcal_depth[event][idx])

        # Select HBHE rechits for this cluster
        mask_cluster = (rh_clusterIndex[event] == idx)
        eta   = ak.to_numpy(rh_eta[event][mask_cluster])
        phi   = ak.to_numpy(rh_phi[event][mask_cluster])
        depth = ak.to_numpy(rh_depth[event][mask_cluster])
        energy = ak.to_numpy(rh_energy[event][mask_cluster])
        time   = ak.to_numpy(rh_time[event][mask_cluster])

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

        # Determine color scale
        vminE, vmaxE = (np.min(energy), np.max(energy)) if len(energy) > 0 else (0,1)
        vminT, vmaxT = (np.min(time), np.max(time)) if len(time) > 0 else (0,1)

        # Draw scatter and cluster outlines
        for d in [1,2,3,4]:
            axE = self.ax[d-1]
            axT = self.ax[d-1+4]

            hit = (depth == d)

            # Draw hits
            if np.sum(hit) > 0:
                axE.scatter(eta[hit], phi[hit], s=80, c=energy[hit],
                            cmap="viridis", vmin=vminE, vmax=vmaxE)
                axT.scatter(eta[hit], phi[hit], s=80, c=time[hit],
                            cmap="plasma", vmin=vminT, vmax=vmaxT)
            else:
                axE.text(0.5, 0.5, "No hits", ha='center', va='center')
                axT.text(0.5, 0.5, "No hits", ha='center', va='center')

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
        scE_sample = axE.collections[0] if np.sum(hit)>0 else axE.scatter([], [], cmap="viridis")
        scT_sample = axT.collections[0] if np.sum(hit)>0 else axT.scatter([], [], cmap="plasma")

        if self.cbar_energy is None:
            self.cbar_energy = self.fig.colorbar(scE_sample, ax=self.ax[0:4], location="right", fraction=0.02)
            self.cbar_energy.set_label("Energy [GeV]")
        else:
            self.cbar_energy.mappable.set_array(energy)
            self.cbar_energy.mappable.set_clim(vminE, vmaxE)

        if self.cbar_time is None:
            self.cbar_time = self.fig.colorbar(scT_sample, ax=self.ax[4:8], location="right", fraction=0.02)
            self.cbar_time.set_label("Time [ns]")
        else:
            self.cbar_time.mappable.set_array(time)
            self.cbar_time.mappable.set_clim(vminT, vmaxT)

        self.fig.suptitle(f"Event {event} — Cluster {idx}, cluster depth {c_depth}", fontsize=16)
        self.fig.canvas.draw_idle()


###############################################################################
# Run
###############################################################################
# Example: start at event 5
viewer = HCALEventDisplay(start_event=5)