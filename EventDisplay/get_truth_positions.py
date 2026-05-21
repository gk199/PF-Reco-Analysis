import sys
import math
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Events, Handle

gen_file   = "/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src/DiPionGun_DR0.1_DT5.0_GEN-SIM-RAW.root"
target_event = 20  # 0-indexed

events = Events(gen_file)
hepmc_handle  = Handle("edm::HepMCProduct")
genp_handle   = Handle("vector<reco::GenParticle>")

for i, event in enumerate(events):
    if i != target_event:
        continue

    # --- GenParticles (eta, phi, pT) ---
    event.getByLabel("genParticles", genp_handle)
    genps = genp_handle.product()
    print(f"Event {i} — GenParticles (pions, pdgId = ±211):")
    print(f"  {'pdgId':>6}  {'pT':>8}  {'eta':>8}  {'phi':>8}  {'vtx_x':>8}  {'vtx_y':>8}  {'vtx_z':>8}")
    for p in genps:
        if abs(p.pdgId()) != 211:
            continue
        vx = p.vx(); vy = p.vy(); vz = p.vz()
        print(f"  {p.pdgId():>6}  {p.pt():8.3f}  {p.eta():8.4f}  {p.phi():8.4f}  {vx:8.2f}  {vy:8.2f}  {vz:8.2f}")

    # --- HepMC vertex times (generatorSmeared = after beam spot smearing) ---
    event.getByLabel("generatorSmeared", hepmc_handle)
    hepmc_evt = hepmc_handle.product().GetEvent()

    c_mm_per_ns = 299.792458  # speed of light in mm/ns
    print(f"\nHepMC particles with production vertex time (generatorSmeared):")
    print(f"  {'barcode':>7}  {'pdgId':>6}  {'pT [GeV]':>9}  {'eta':>8}  {'phi':>8}  {'vtx_t [ns]':>10}  {'t_prop [ns]':>11}")
    for bc in range(1, hepmc_evt.particles_size() + 1):
        p = hepmc_evt.barcode_to_particle(bc)
        if p is None:
            continue
        prod_vtx = p.production_vertex()
        if prod_vtx is None:
            t_ns = t_prop = float('nan')
        else:
            pos = prod_vtx.position()  # mm
            t_ns  = pos.t() / c_mm_per_ns
            t_prop = math.sqrt(pos.x()**2 + pos.y()**2 + pos.z()**2) / c_mm_per_ns
        px, py = p.momentum().px(), p.momentum().py()
        pt = math.sqrt(px**2 + py**2)  # already in GeV
        print(f"  {bc:>7}  {p.pdg_id():>6}  {pt:9.3f}  {p.momentum().eta():8.4f}  {p.momentum().phi():8.4f}  {t_ns:10.3f}  {t_prop:11.3f}")
    break
