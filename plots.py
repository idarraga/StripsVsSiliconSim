from ROOT import *

def getDistribution(fn, toDraw, histoName, histoParameters, cut):
    T = fn.Get("T")
    toDraw += ">>"
    toDraw += histoName+histoParameters
    T.Draw(toDraw, cut)
    return gROOT.FindObject(histoName)
    

fStrips = TFile("./datasave/output_Strips.root", "READ")
fSi = TFile("./datasave/output_Si.root", "READ")

print fStrips
print fSi

# Get distributions
hStrips = getDistribution(fStrips,"Sum$(scatteredAngle)*180/TMath::Pi()", "h_scatteredAngle_Strips", "(300,0,3)", "")
hSi = getDistribution(fSi,"Sum$(scatteredAngle)*180/TMath::Pi()", "h_scatteredAngle_Si", "(300,0,3)", "")

c1 = TCanvas("SiVsStrips","SiVsStrips")
c1.cd()
c1.SetLogy()

hSi.SetTitle("")
hSi.Draw()
hSi.SetLineColor(ROOT.kBlue)
hSi.SetLineWidth(2)

hSi.GetXaxis().SetTitle("Total deviation angle (theta) for primaries (degrees)")
hSi.GetYaxis().SetTitle("entries")

hStrips.SetTitle("")
hStrips.Draw("same")
hStrips.SetLineColor(ROOT.kRed)
hStrips.SetLineWidth(2)

leg = TLegend(0.6, 0.6, 0.9, 0.9)
leg.SetFillColor(ROOT.kWhite)
#leg.SetBorder(1)
leg.AddEntry(hSi, "Silicon setup", "L")
leg.AddEntry(hStrips, "Strips setup", "L")
leg.Draw()

iTotal_hSi = hSi.Integral(  )
iTotal_hStrips = hStrips.Integral( )

i_hSi = hSi.Integral( hSi.FindBin(1), hSi.GetNbinsX() )
i_hStrips = hStrips.Integral( hStrips.FindBin(1), hStrips.GetNbinsX() )

i_hSi = i_hSi/iTotal_hSi
i_hStrips = i_hStrips/iTotal_hStrips

print( "Si = %f , Strips = %f "%(i_hSi, i_hStrips) )


leg = TLegend(0.6, 0.6, 0.9, 0.9)
leg.SetFillColor(ROOT.kWhite)
#leg.SetBorder(1)
leg.AddEntry(hSi, "Si | frac = " + str ( TString.Format("%.2f "%(i_hSi*100.0)).Data() ) + '%', "L")
leg.AddEntry(hStrips, "Strips | frac = " + str( TString.Format("%.2f "%(i_hStrips *100.0)).Data() ) + '%', "L")
leg.Draw()

lat1 = TLatex()
lat1.DrawLatex(1.1, 2000, "Si/Strips = " + str( TString.Format("%.2f"%(i_hSi/i_hStrips)).Data() ) )

line1 = TLine(1,0, 1, 1E4)
line1.SetLineStyle(2)
line1.Draw()

c1.Update()

##########################################################
# Energy
h_edep_Strips = getDistribution(fStrips,"Sum$(scatteredAngle)*180/TMath::Pi():Sum$(edep)", "h_edep_Strips", "(200,0,1500,100,0,3)", "")
h_edep_Si = getDistribution(fSi,"Sum$(scatteredAngle)*180/TMath::Pi():Sum$(edep)", "h_edep_Si", "(200,0,1500,100,0,3)", "")

c2 = TCanvas("SiVsStrips_Edep","SiVsStrips_Edep")
c2.Divide(1,2)
c2.cd(1)
#c2.SetLogy()

h_edep_Si.SetTitle("Si")
h_edep_Si.Draw("colz")
h_edep_Si.GetXaxis().SetTitle("edep (keV)")
h_edep_Si.GetYaxis().SetTitle("angle")

c2.cd(2)

h_edep_Strips.SetTitle("Strips")
h_edep_Strips.Draw("colz")
h_edep_Strips.GetXaxis().SetTitle("edep (keV)")
h_edep_Strips.GetYaxis().SetTitle("angle")



