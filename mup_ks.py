#!/usr/bin/env python3
import ROOT
import os
# --- Configuration ---
generic_dirs = {
    'qqbar':   ["JpsiKs_dd_0417", "JpsiKs_uu_0417", "JpsiKs_ss_0417"],
    'mixed':   ["JpsiKs_mixed_0417"],
    'tau':     ["JpsiKs_tau_0417"],
    'charged': ["JpsiKs_charged_0417"]  
}
signal_dir   = "signal_ks_mm"
root_filename = "signal_fourchannel_merged.root"
# Function to find subdirectories matching 'sub*'
def find_subdirs(base_dir):
    try:
        return [d for d in os.listdir(base_dir)
                if os.path.isdir(os.path.join(base_dir, d)) and d.startswith('sub')]
    except FileNotFoundError:
        return []
# Initialize TChains for each component
chains = { label: ROOT.TChain("cut_muons_all")
           for label in list(generic_dirs.keys()) + ['signal'] }
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
x_min, x_max = 0.0, 4.0

# 定义 goodtrack 条件
goodtrack = 'dr < 0.5 && abs(dz) < 2 && thetaInCDCAcceptance && muonIDNN_my > 0.8'

# Create and fill histograms
hists = {}
colors = {
    'qqbar':   ROOT.kRed,
    'mixed':   ROOT.kGreen+2,
    'tau':     ROOT.kMagenta,
    'charged': ROOT.kCyan,
    'signal':  ROOT.kBlue
}
for label, chain in chains.items():
    hist_name = f"h_{label}"
    title = f"{'Muon momentum at Jpsi Ks channel' if label=='qqbar' else label.capitalize()}" + ";Muon momentum;Entries"
    h = ROOT.TH1F(hist_name, title, bin_count, x_min, x_max)
    
    # 根据不同样本设置 cut
    if label == 'signal':
        sample_cut = 'isSignal==1'
    elif label in ['mixed', 'charged']:
        sample_cut = 'isSignal!=1'
    else:
        sample_cut = 'isSignal!=1'
    
    # 将样本 cut 与 goodtrack 条件结合
    combined_cut = f"({sample_cut}) && ({goodtrack})"
    
    chain.Draw(f"p >> {hist_name}", combined_cut, 'goff')
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(colors[label])
    h.SetLineWidth(2)
    h.SetStats(0)  # disable stats box
    hists[label] = h

# Style canvas
draw_order = ['qqbar', 'mixed', 'charged', 'tau', 'signal']
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas('c', 'muon momentum at Jpsi Ks channel', 900, 650)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetBottomMargin(0.12)
c.SetGrid()

# 找出所有直方图中的最大值
y_max = 0
for label in draw_order:
    hist_max = hists[label].GetMaximum()
    if hist_max > y_max:
        y_max = hist_max

# 添加一些边距（增加 20%）
y_max *= 1.2

# 绘制直方图并设置 Y 轴范围
drawn = False
for label in draw_order:
    if not drawn:
        # 对第一个直方图设置 Y 轴范围
        hists[label].GetYaxis().SetRangeUser(0, y_max)
        hists[label].Draw("HIST")
        drawn = True
    else:
        hists[label].Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
for label in draw_order:
    display = label
    leg.AddEntry(hists[label], display, 'l')
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw()

# 保存并保持画布打开
output = 'ep_components.png'
c.SaveAs(output)
print(f"Saved plot as {output}")
ROOT.gApplication.Run()
