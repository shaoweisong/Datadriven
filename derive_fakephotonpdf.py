#/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
import ROOT
from array import array
import awkward
import sys
sys.path.append("/eos/user/s/shsong/pkgs_condor/parquet_to_root-0.3.0/")
from parquet_to_root import parquet_to_root
import pandas
#run mk_pp_gjet.py firstly to get the Data_2018.root, GJet.root
GJets1=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/GJet/GJet.root")
# GJets1=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/UL16dd/GJet16.root")
GJets1_tree=GJets1.Get("GJet")
cuts="(Diphoton_mass<115 || Diphoton_mass>135)*(minID_genPartFlav!=1)"
# cuts="(Diphoton_mass<115 || Diphoton_mass>135)*(minID_genPartFlav!=1)"
h_minphotonID=ROOT.TH1F("h_minphotonID_gjet","h_minphotonID_gjet",19,-0.9,1)
GJets1_tree.Project("h_minphotonID_gjet","Diphoton_minID_modified","weight_central*"+cuts)
photonIDPDF_fake=ROOT.TF1("photonIDPDF_fake","pol10",-0.9,1.)

h_minphotonID.Fit(photonIDPDF_fake,"R")
fit_result = h_minphotonID.Fit(photonIDPDF_fake, "RS")  
chi2 = photonIDPDF_fake.GetChisquare()
ndf = photonIDPDF_fake.GetNDF()
p_value = ROOT.TMath.Prob(chi2, ndf)
print("chisquare:" ,chi2)
print("freedom degree:", ndf)
print("p_value:", p_value)
# 构造并打印多项式表达式
poly_terms = []
for i in range(11):
    coeff = photonIDPDF_fake.GetParameter(i)
    if abs(coeff) < 1e-10:  # 忽略极小值
        continue
    term = f"{coeff:.5g}"
    if i == 1:
        term += "*x"
    elif i > 1:
        term += f"*x^{i}"
    poly_terms.append(term)

polynomial_expr = " + ".join(poly_terms)

display_terms = []
for i in range(11):  # 只取前6项以避免太长
    coeff = photonIDPDF_fake.GetParameter(i)
    if abs(coeff) < 1e-10:
        continue
    coeff_str = f"{coeff:.2g}"
    if i == 0:
        term = f"{coeff_str}"
    elif i == 1:
        term = f"{coeff_str}x"
    else:
        term = f"{coeff_str}x^{{{i}}}"  # 用 LaTeX 格式加大括号
    display_terms.append(term)

latex_poly = " + ".join(display_terms)
latex_poly = "f(x) = " + latex_poly

# 创建画布并绘图
c1 = ROOT.TCanvas("c1", "c1", 600, 800)
h_minphotonID.Draw("E1")
photonIDPDF_fake.Draw("SAME")  # 叠加拟合曲线

# 在图上写上多项式公式
latex = ROOT.TLatex()
latex.SetNDC()  # 使用归一化坐标 (0~1)
latex.SetTextSize(0.03)
latex.SetTextAlign(13)  # 左上角对齐
latex.DrawLatex(0.15, 0.85, latex_poly)


# c1=ROOT.TCanvas("c1","c1",600,800)
# h_minphotonID.Draw("E1")
c1.SaveAs("/eos/user/s/shsong/www/datadriven/fakephoton_pdf16postnew.png")
