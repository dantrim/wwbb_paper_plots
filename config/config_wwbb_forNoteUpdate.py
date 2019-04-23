import sys

from dantrimania.python.analysis.utility.plotting.plot1d import plot1d
import dantrimania.python.analysis.utility.samples.sample as sample
import dantrimania.python.analysis.utility.samples.region as region

##############################################################################
# additional variables
##############################################################################

##################################################################################
# sample definition
##################################################################################
filelist_dir = "/data/uclhc/uci/user/dantrim/n0307val/susynt-read/filelists/"
filelist_dir_data = "/data/uclhc/uci/user/dantrim/n0307val/susynt-read/filelists/"
h5_dir_mc = "/data/uclhc/uci/user/dantrim/ntuples/n0307/a_feb10/mc/h5/"
#h5_dir_mc = "/data/uclhc/uci/user/dantrim/ntuples/n0307/b_feb18/mc/h5/"
h5_dir_data = "/data/uclhc/uci/user/dantrim/ntuples/n0307/a_feb10/data/h5/"

lumi_factor = 140.48

tags = ["mc16a", "mc16d", "mc16e"]

which_data = '1516'
sf_vals_tt = { '1516' : 1.00, '151617' : 1.00 }
sf_vals_wt = { '1516' : 1.00, '151617' : 1.00 }
sf_vals_tt = { '1516' : 1.00, '151617' : 1.00 }
sf_vals_wt = { '1516' : 1.00, '151617' : 1.00 }

# backgrounds
#top = sample.Sample("topDS", "Top (DS)")
#top.scalefactor = lumi_factor * 0.93
#top.fillstyle = 0
#top.linestyle = '-'
#top.color = "#0057FF"
#top.load(filelist_dir + "topDS_mc16a", h5_dir_mc, tags = tags)
#loaded_samples.append(top)

top_sf =  1.0
z_sf = 1.0

top = sample.Sample("Top", "Top")
top.scalefactor = lumi_factor * top_sf
top.fillstyle = 0
top.linestyle = "-"
top.color = "#057390"
top.load(filelist_dir + "top_mc16a", h5_dir_mc, tags = tags)
loaded_samples.append(top)

zhf = sample.Sample("ZjetsHF", "$Z/\\gamma*$+jets HF")
zhf.scalefactor = lumi_factor * z_sf
zhf.fillstyle = 0
zhf.linestyle = "-"
zhf.color = "#fc8f1e"
zhf.load(filelist_dir + "zjets_and_dy_sherpa_mc16a", h5_dir_mc, tags = tags)
loaded_samples.append(zhf)

zlf = sample.Sample("ZjetsLF", "$Z/\\gamma*$+jets LF")
zlf.scalefactor = lumi_factor
zlf.fillstyle = 0
zlf.linestyle = "-"
zlf.color = "#fc8f1a"
zlf.load(filelist_dir + "zjets_and_dy_sherpa_mc16a", h5_dir_mc, tags = tags)
loaded_samples.append(zlf)

higgs = sample.Sample("SingleHiggs", "Higgs")
higgs.scalefactor = lumi_factor
higgs.fillstyle = 0
higgs.linestyle = '-'
higgs.color = '#d93b3b'
higgs.load(filelist_dir + "higgs_mc16a", h5_dir_mc, tags = tags) #['mc16a', 'mc16d', 'mc16e'])
loaded_samples.append(higgs)

dib = sample.Sample("Diboson", "$VV$")
dib.scalefactor = lumi_factor
dib.fillstyle = 0
dib.linestyle = '-'
dib.color = "#f6f7eb"
dib.load(filelist_dir + "diboson_sherpa_mc16a", h5_dir_mc, tags = tags) #["mc16a", "mc16d", "mc16e"])
loaded_samples.append(dib)

ttv = sample.Sample("TTbarV", "$t\\bar{t} + V$")
ttv.scalefactor = lumi_factor
ttv.fillstyle = 0
ttv.linestyle = '-'
ttv.color = "#353531"
ttv.load(filelist_dir + "ttV_mc16a", h5_dir_mc, tags = tags) #["mc16a", "mc16d", "mc16e"])
loaded_samples.append(ttv)

hh = sample.Sample("hhWWbb", "SM $hh$ (arbitrary $\\sigma$)")
hh.is_signal = True
hh.scalefactor = lumi_factor * 5000#* 350
hh.fillstyle = 0
hh.linestyle = '--'
hh.color = 'r'
hh.load(filelist_dir + "hh_wwbb_mc16a", h5_dir_mc, tags = tags)

#wjets = sample.Sample("WjetsFull", "$W + jets$")
#wjets.scalefactor = lumi_factor
#wjets.fillstyle = 0
#wjets.linestyle = '-'
#wjets.color = "#726e97"
#wjets.load(filelist_dir + "wjets_sherpa_mc16a", h5_dir_mc, tags = tags)
#loaded_samples.append(wjets)

## data
data = sample.Sample("data15161718", "Data")
data.is_data = True
data.scalefactor = 1.0
data.fillstyle = 0
data.linestyle = '-'
data.color = 'k'
data.load(filelist_dir_data + "n0307_data15161718", h5_dir_data)
loaded_samples.append(data)

#############################################################
# region definitions
#############################################################
trigger = "( ( year == 2015 && trig_tight_2015 == 1 ) || ( year == 2016 && trig_tight_2016 == 1 ) || ( year == 2017 && trig_tight_2017rand == 1 ) || ( year == 2018 && trig_tight_2018 == 1 ))" 

r = region.Region("srIncNoDhh", "srIncNoDhh")
r.tcut = "mll>20 && mll<60 && nBJets>=2 && mbb>110 && mbb<140"
loaded_regions.append(r)
