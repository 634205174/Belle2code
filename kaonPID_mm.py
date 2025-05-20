#!/usr/bin/env python3

import ROOT
import os

# --- Configuration ---
generic_dirs = {
    'qqbar': ["JpsiKs_dd_0508", "JpsiKs_uu_0508", "JpsiKs_ss_0508"],
    'mixed': ["JpsiKs_mixed_0508"],
    'tau':   ["JpsiKs_tau_0508"]
}
signal_dir = "signal_kp_mm"
root_filename = "signal_fourchannel_merged.root"

# Function to find subdirectories matching 'sub*'
def find_subdirs(base_dir):
    try:
        return [d for d in os.listdir(base_dir)
                if os.path.isdir(os.path.join(base_dir, d)) and d.startswith('sub')]
    except FileNotFoundError:
        return []

# Initialize TChains for each component
chains = {label: ROOT.TChain("cut_kaons_all") for label in list(generic_dirs.keys()) + ['signal']}

# Add files for generic components
def add_generic(label, base_list):
    for base in base_list:
        subs = find_subdirs(base)
        if not subs:
            print(f"Warning: no subdirs found in {base}")
        for sub in subs:
            path = os.path.join(base, sub, root_filename)
            if os.path.exists(path):
                chains[label].Add(path)
            else:
                print(f"Warning: {label} file not found: {path}")

for label, bases in generic_dirs.items():
    add_generic(label, bases)

# Add files for signal component
subs = find_subdirs(signal_dir)
if not subs:
    print(f"Warning: no subdirs found in {signal_dir}")
for sub in subs:
    path = os.path.join(signal_dir, sub, root_filename)
    if os.path.exists(path):
        chains['signal'].Add(path)
    else:
        print(f"Warning: signal file not found: {path}")

# Debug: print entry counts
for label, chain in chains.items():
    print(f"Entries in {label}: {chain.GetEntries()}")

# Histogram settings
bin_count = 100
x_min, x_max = 0.0, 1.0

# Create and fill histograms
hists = {}
colors = {'qqbar': ROOT.kRed, 'mixed': ROOT.kGreen+2, 'tau': ROOT.kMagenta, 'signal': ROOT.kBlue}
for label, chain in chains.items():
    hist_name = f"h_{label}"
    title = f"{'kaon mumu channel PID' if label=='qqbar' else label.capitalize()}" + ";kaonIDNN;Entries"
    h = ROOT.TH1F(hist_name, title, bin_count, x_min, x_max)
    cut = 'isSignal==0' if label in ['qqbar', 'mixed', 'tau'] else 'isSignal==1'
    chain.Draw(f"kaonIDNN_my >> {hist_name}", cut, 'goff')
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(colors[label])
    h.SetLineWidth(2)
    h.SetStats(0)  # disable stats box
    hists[label] = h

# Style canvas
draw_order = ['qqbar', 'mixed', 'tau', 'signal']
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas('c', 'kaon PID ', 900, 650)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetBottomMargin(0.12)
c.SetGrid()

# Draw histograms
drawn = False
for label in draw_order:
    opt = 'HIST' if not drawn else 'HIST SAME'
    hists[label].Draw(opt)
    drawn = True

# Legend
leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
for label in draw_order:
    display = 'QQbar' if label=='qqbar' else label.capitalize()
    leg.AddEntry(hists[label], display, 'l')
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw()

# Save and keep canvas open
output = 'kaonIDNN_components.png'
c.SaveAs(output)
print(f"Saved plot as {output}")
ROOT.gApplication.Run()

