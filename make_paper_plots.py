#!/usr/bin/env python

import sys, os
import argparse

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import dantrimania.python.analysis.utility.samples.sample as sample
import dantrimania.python.analysis.utility.samples.sample_utils as sample_utils
import dantrimania.python.analysis.utility.samples.region_utils as region_utils
import dantrimania.python.analysis.utility.utils.utils as utils
import dantrimania.python.analysis.utility.samples.sample_cacher as sample_cacher
import dantrimania.python.analysis.utility.plotting.m_py.errorbars as errorbars
import dantrimania.python.analysis.utility.utils.plib_utils as plib
from dantrimania.python.analysis.utility.plotting.histogram1d import histogram1d
from dantrimania.python.analysis.utility.plotting.histogram_stack import histogram_stack
from dantrimania.python.analysis.utility.plotting.ratio_canvas import ratio_canvas


plt = plib.import_pyplot()
import matplotlib
from matplotlib.ticker import ScalarFormatter 
import numpy as np
import json

# helvetica
#from matplotlib import rc
#plt.rcParams['ps.useafm'] = True
#plt.rcParams['font.sans-serif'] = 'Helvetica'
#plt.rcParams['pdf.fonttype'] = 42
#print plt.rcParams
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

class PlotDescription :

    def __init__(self, descriptor = "", abs_val = False) :
        self.descriptor = descriptor

        # variable stuff
        self.var_to_plot = ""
        self.is_abs = abs_val

        # region struff
        self.region_name = ""
        self.is_cr = False

        self.load()

    def load(self) :

        fields = self.descriptor.strip().split(":")
        if len(fields) > 2 :
            print "PlotDescription::load   Input descriptor (=%s) is not in expected format" % self.descriptor
            sys.exit()
        input_variable = fields[0]
        region = fields[1]

        if "abs(" in input_variable :
            self.is_abs = True
        input_variable = input_variable.replace("abs(","").replace(")","")
        self.var_to_plot = input_variable

        self.region_name = region
        if "cr" or "vr" in region :
            self.is_cr = True

    def __str__(self) :
        return "PlotDescription    variable: %s (abs? %s), region: %s" % (self.var_to_plot, self.is_abs, self.region_name)

def get_variables_from_tcut(tcut) :

    operators = ["==", ">=", "<=", ">", "<", "!=", "*", "-"]
    logics = ["&&", "||", ")", "(", "abs"]
    vars_only = tcut
    for op in operators :
        vars_only = vars_only.replace(op, " ")
    for log in logics :
        vars_only = vars_only.replace(log, " ")
    vars_only = vars_only.split()
    out = []
    for v in vars_only :
        if v not in out and not v.isdigit() :
            try :
                flv = float(v)
            except :
                out.append(v)
    #vars_only = [v for v in vars_only if not v.isdigit()]

    return out

def get_required_variables(region_to_plot, plot_description) :

    variables = []
    variables.append(plot_description.var_to_plot)
    
    tcut = region_to_plot.tcut
    selection_variables = get_variables_from_tcut(tcut)
    for sv in selection_variables :
        if sv not in variables :
            variables.append(sv)
    
    # we always need the eventweight, which will not show up in the tcut
    variables.append("eventweight")
    variables.append("eventweightNoPRW")
    variables.append("eventweightbtag")
    variables.append("eventweightbtagNoPRW")
    variables.append("eventweightBtagJvt")
    variables.append("eventweight_multi")
    variables.append("eventweightNoPRW_multi")
    variables.append("eventweightbtag_multi")
    variables.append("eventweightBtagJvt_multi")
    variables.append("pupw_period")
    variables.append("mt2_bb")
    variables.append("NN_d_hh")
    variables.append("NN_d_top")
    variables.append("NN_d_zsf")
    variables.append("NN_d_ztt")
    variables.append("isEE")
    variables.append("isMM")
    variables.append("isSF")
    variables.append("isDF")
    variables.append("year")
    variables.append("trig_tight_2015")
    variables.append("trig_tight_2016")
    variables.append("trig_tight_2017rand")
    variables.append("trig_tight_2018")
    variables.append('mll')
    variables.append('mbb')
    variables.append("isBB")
    variables.append("isCC")
    variables.append("isBC")
    variables.append("isBL")
    variables.append("isCL")
    variables.append("isLL")
    variables.append("dRll")
    variables.append("HT2Ratio")
    variables.append("dphi_ll")
    #variables.append("dphi_bb")
    variables.append("met")
    variables.append("mt2_bb")
    #variables.append("pTll")
    #variables.append("dphi_met_ll")
    #variables.append("met_pTll")
    variables.append("l0_pt")
    variables.append("l1_pt")
    variables.append("bj0_pt")
    variables.append("bj1_pt")
    return variables

def get_vrz_mll_cut(arr) :

    flav = arr["isSF"]
    mll = arr["mll"]
    idx_low = (mll>71.2) & (mll<81.2)
    idx_high = (mll>101.2) & (mll<115) 
    idx_flav = flav == 1
    idx_pass_vrz = (idx_low | idx_high) & idx_flav
    return idx_pass_vrz

def get_trigger_idx(arr) :

    idx_15 = (arr['year'] == 2015) & (arr['trig_tight_2015'] == 1)
    idx_16 = (arr['year'] == 2016) & (arr['trig_tight_2016'] == 1)
    idx_17 = (arr['year'] == 2017) & (arr['trig_tight_2017rand'] == 1)
    idx_18 = (arr['year'] == 2018) & (arr['trig_tight_2018'] == 1)
    idx = (idx_15 | idx_16 | idx_17 | idx_18)
    return idx

def legend_and_labels_type(region_name = "", var_name = "") :

    key = "{}_{}".format(region_name, var_name)
    ll = {}
    pass

def bounds_dict() :

    v = {}
    sr_sf_cut_bins = np.arange(5.45, 12+0.45, 0.45)
    v["NN_d_hh"] = {
                    "srIncNoDhh" : [1, -11, 11],
                    "srSFNoDhh" : [1, -11, 12],
                    "srDFNoDhh" : [1, -12, 10],
                        "srPreSel" : [1, -11, 11],
                    "srIncNoDhhClose" : [1, 0, 12],
                    "srSFNoDhhClose" : [0.45, 5.45, 12],
                    "srDFNoDhhClose" : [0.5, 4, 9],
                    #"srSFNoDhhCloseCut" : [5.45, 7, 8, 9, 10, 11, 12], #, 6.5, 7, 7.5, 8, 9, 10, 11, 12], #, 5.45, 12], #sr_sf_cut_bins,
                    #"srSFNoDhhCloseCut" : [5.45, 5.90, 6.35, 7, 8, 10, 12], #, 8.15, 8.60, 9.05],
                    "srSFNoDhhCloseCut" : [1, 5.45, 11.45],
                    "srDFNoDhhCloseCut" : [1, 5.55, 9.55],
                    "crTopNoDhh" : [1, -12, 12],
                    "crTop" : [1, 4.5, 11.5],
                    "vrTop" : [1, 4.5, 11.5],
                    "vrZ" : [1, 0, 6],
                    "crZNoDhh" : [1, -12, 12],
                    "crZ" : [1, 0, 8],
    }
    v["mbb"] = {
        "srIncNoDhhNoMbb" : [10, 0, 300]
        ,"srIncNoDhh" : [5, 110, 140]
        ,"crTop" : [50, 0, 500]
        ,"vrTop" : [5, 100, 140]
        ,"crZ" : [5, 100, 140]
        ,"vrZ" : [5, 100, 140]
        ,"srIncNoMbbDhh" : [20, 20, 300]
    }
    v["HT2Ratio"] = {
        "srIncNoDhh" : [0.1, 0, 1]
        ,"crTop" : [0.1, 0.2, 1.0]
        ,"vrTop" : [0.1, 0.2, 1.0]
        ,"crZ"   : [0.1, 0.5, 1 ]
        ,"vrZ"   : [0.1, 0.5, 1.0 ]
        ,"srIncNoMbbDhh" : [0.1, 0.6, 1.0]
    }
    v["dRll"] = {
        "srIncNoDhh" : [0.2, 0, 3.6]
        ,"crTop" : [0.2, 0.2, 2.2]
        ,"vrTop" : [0.2, 0, 2.4]
        ,"crZ"   : [0.4, 0, 2.8]
        ,"vrZ"   : [0.4, 0, 3.2]
    }
    v["mt2_bb"] = {
        "srIncNoDhh" : [30, 0, 300]
        ,"crTop" : [40, 0, 360]
        ,"vrTop" : [40, 0, 360]
        ,"crZ"   : [20, 0, 160]
        ,"vrZ"   : [20, 0, 200]
    }
    v["dphi_ll"] = {
        "srIncNoDhh" : [0.2, 0, 3.2]
        ,"crTop" : [0.2, 0, 1.8]
        ,"vrTop" : [0.2, 0, 1.8]
        ,"crZ" : [0.2, 0, 1.8]
        ,"vrZ" : [0.2, 0, 1.8]
    }

    v["mll"] = {
        "srIncNoMllDhh" : [5, 20, 80]
    }
    #v["NN_d_hh"] = { "srIncNoDhh" : [-11, -5, -4, -3, 11] }
    return v

def region_nice_names_dict() :

    dhh_text = r"\textit{d}$_{\mbox{\small{\textit{HH}}}}$"
    mbb_text = r"\textit{m}$_{\mbox{\small{\textit{bb}}}}$"
    mll_text = r"\textit{m}$_{\textit{\normalsize{ll}}}$"
    n = {}
    n["srIncNoDhhNoMbb"] = "$m_{\\ell\\ell} \\in [20,60]$ GeV"
    n["srIncNoDhh"] = "SR, SF+DF and no %s cut" % dhh_text
    n["srSFNoDhh"] = "SR-SF, no %s cut" % dhh_text
    n["srDFNoDhh"] = "SR-DF, no %s cut" % dhh_text
    n["srPreSel"] = "Pre-selection"
    n["srIncNoDhhClose"] = "SR, SF+DF and no %s cut" % dhh_text
    n["srSFNoDhhClose"] = "SR-SF, no %s cut" % dhh_text
    n["srDFNoDhhClose"] = "SR-DF, no %s cut" % dhh_text
    n["srSFNoDhhCloseCut"] = "SR-SF"
    n["srDFNoDhhCloseCut"] = "SR-DF"
    n["crTopNoDhh"] = "CR-Top, no %s cut" % dhh_text
    n["crTop"] = "CR-Top"
    #n["vrTop"] = "VR-Top"
    n["vrTop"] = "VR-1"
    #n["vrZ"] = "VR-Z+HF"
    n["vrZ"] = "VR-2"
    n["crZNoDhh"] = "CR-Z+HF, no %s cut" % dhh_text
    n["crZ"] = "CR-Z+HF"
    n["srIncNoMbbDhh"] = "SR, SF+DF and no %s cut" % mbb_text
    n["srIncNoMllDhh"] = "SR, SF+DF and no %s cut" % mll_text
    return n

def nice_names_dict() :

    n = {}
    n["NN_d_hh"] = "$d_{HH}$"
    n["mbb"] = "$m_{bb}$ [GeV]"
    n["HT2Ratio"] = "$H_{T2}^{R}$"
    n["dRll"] = "$\\Delta R _{\\ell \\ell}$"
    n["mt2_bb"] = "$m_{T2}^{bb}$ [GeV]"
    n["dphi_ll"] = "$|\\Delta \\phi_{\\ell \\ell}|$ [rad.]"
    n["mll"] = "$m_{\\ell\\ell}$ [GeV]"

    n["NN_d_hh"] = r"\textit{d}$_{\mbox{\small{\textit{HH}}}}$"
    n["mbb"] = r"\textit{m}$_{\mbox{\small{\textit{bb}}}}$ [GeV]"
    n["HT2Ratio"] = r"\textit{H}$_{\mbox{\textit{\normalsize{T2}}}}^{\mbox{\textit{\normalsize{R}}}}$"
    n["dRll"] = r"$\Delta$\textit{R}$_{\textit{\normalsize{ll}}}$"
    n["mt2_bb"] = r"\textit{m}$_{\mbox{\textit{\normalsize{T2}}}}^{\mbox{\textit{\normalsize{bb}}}}$ [GeV]"
    n["dphi_ll"] = r"$|\Delta\phi_{\mbox{\textit{\normalsize{ll}}}}|$ [rad.]"
    n["mll"] = r"\textit{m}$_{\textit{\normalsize{ll}}}$ [GeV]"
    return n

def add_labels(pad, region_name = "", var_name = "") :


    x_atlas = 0.04
    y_atlas = 0.96
    x_type_offset = 0.22
    y_type = 0.96

    x_lumi = 0.04
    y_lumi = 0.88

    x_region = 0.042
    y_region = 0.80

    if (region_name == "srSFNoDhhCloseCut" or region_name == "srDFNoDhhCloseCut") and var_name == "NN_d_hh" :
        x_atlas = 0.29
        y_atlas = 0.96

        y_type = 0.96

        x_lumi = 0.29
        y_lumi = 0.87

        x_region = 0.29
        y_region = 0.79

    if (region_name == "crZ" or region_name == "crTop" or region_name == "vrTop" or region_name == "vrZ") and (var_name == "NN_d_hh" or var_name == "dRll" or var_name == "dphi_ll" or var_name == "mt2_bb"):
        x_atlas = 0.33
        y_atlas = 0.96
        y_type = 0.96

        x_lumi = 0.33
        y_lumi = 0.88

        x_region = 0.33
        y_region = 0.8

    if (region_name == "crZ" or region_name == "vrZ") and (var_name == "dRll" or var_name == "dphi_ll") :
        x_atlas = 0.04
        y_atlas = 0.97
        x_type_offset = 0.24
        y_type = 0.97

        x_lumi = 0.04
        y_lumi = 0.88

        x_region = 0.042
        y_region = 0.80

    if (region_name == "crTop") and (var_name == "mt2_bb") :
        x_atlas = 0.04
        y_atlas = 0.97
        x_type_offset = 0.24
        y_type = 0.97

        x_lumi = 0.04
        y_lumi = 0.88

        x_region = 0.042
        y_region = 0.80

    # ATLAS label
    size = 24
    #text = 'ATLAS'
    text = '\\textbf{\\textit{ATLAS}}'
    opts = dict(transform = pad.transAxes)
    opts.update( dict(va = 'top', ha = 'left') )
    #pad.text(x_atlas, y_atlas, text, size = size, style = 'italic', weight = 'bold', **opts)
    pad.text(x_atlas, y_atlas, text, size = size, **opts) #style = 'italic', weight = 'bold', **opts)

    #what_kind = 'Preliminary'
    what_kind = ''
    pad.text(x_atlas + x_type_offset, y_type, what_kind, size = size, **opts)

    lumi = '139'#.5'
    pad.text(x_lumi, y_lumi, '$\\sqrt{s} =$ 13 TeV, %s fb$^{-1}$' % lumi, size = 0.75 * size, **opts)

    if "no" in region_nice_names_dict()[region_name] and "cut" in region_nice_names_dict()[region_name] :
        size = 0.9 * size
        

    # region
    if (region_name == "srIncNoDhh" or region_name == "srIncNoMbbDhh" or region_name == "srSFNoDhh" or region_name == "srDFNoDhh") and (var_name == "NN_d_hh" or var_name == "mbb") :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        #region_text = "\t%s" % region_nice_names_dict()[region_name]
        region_text = "%s" % region_nice_names_dict()[region_name]
        pad.text( 1.0 * x_region, 0.91 * y_region, region_text, size = 0.75 * size, **opts)
    elif (region_name == "srIncNoMllDhh") and (var_name == "mll") :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        #region_text = "\t%s" % region_nice_names_dict()[region_name]
        region_text = "   %s" % region_nice_names_dict()[region_name]
        pad.text( 1.0 * x_region, 0.91 * y_region, region_text, size = 0.75 * size, **opts)
    elif (region_name == "srIncNoDhh" or region_name == "srIncNoMbbDhh") and var_name == "HT2Ratio" :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        region_text = "   %s" % region_nice_names_dict()[region_name]
        pad.text( 1.0 * x_region, 0.91 * y_region, region_text, size = 0.75 * size, **opts)
    elif region_name == "srIncNoDhh" and var_name == "dRll" :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        region_text = "   %s" % region_nice_names_dict()[region_name]
        pad.text( 1.0 * x_region, 0.91 * y_region, region_text, size = 0.75 * size, **opts)
    elif region_name == "srIncNoDhh" and var_name == "mt2_bb" :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        region_text = "   %s" % region_nice_names_dict()[region_name]
        pad.text( 1.0 * x_region, 0.91 * y_region, region_text, size = 0.75 * size, **opts)
    elif region_name == "srIncNoDhh" and var_name == "dphi_ll" :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        region_text = "  %s" % region_nice_names_dict()[region_name]
        pad.text(x_region, 0.91 * y_region, region_text, size = 0.75 * size, **opts)
    else :
        region_text = "Selection: %s" % region_nice_names_dict()[region_name]
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)

def make_legend(ordered_labels, n_cols, var_name, region_name, pad) :

    handles, labels = pad.get_legend_handles_labels()
    new_handles = []
    new_labels = []
    for l in ordered_labels :
        for il, label in enumerate(labels) :
            #if label == l :
            if l in label :
                new_labels.append(l.replace("SIG",""))
                new_handles.append(handles[il])

    leg_x, leg_y = 0.45, 0.75
    legend_fontsize = 15

    if n_cols == 1 :
        leg_x = 0.58
        leg_y = 0.58

    if n_cols == 1 and (region_name == "srIncNoDhh" or region_name == "srIncNoMbbDhh" or region_name == "srSFNoDhh" or region_name == "srDFNoDhh") and (var_name == "dRll" or var_name == "HT2Ratio"
                        or var_name == "NN_d_hh" or var_name == "mbb" or var_name == "mt2_bb"
                        or var_name == "dphi_ll"
    ) :
        leg_x = 0.58
        #leg_y = 0.71
        leg_y = 0.71

    if n_cols == 1 and (region_name == "srIncNoMllDhh") and (var_name == "mll") :
        leg_x = 0.58
        leg_y = 0.71

    # sr
    if n_cols == 2 and (region_name == "srSFNoDhhCloseCut" or region_name == "srDFNoDhhCloseCut") and var_name == "NN_d_hh" :
        leg_x = 0.27
        leg_y = 0.57
        legend_fontsize = 14

    # crZ
    if n_cols == 1 and (region_name == "crZ") and var_name == "NN_d_hh" :
        leg_x = 0.48
        leg_y = 0.3
        legend_fontsize = 16
    if n_cols == 2 and (region_name == "crZ" or region_name == "crTop" or region_name == "vrTop" or region_name == "vrZ") and var_name == "NN_d_hh" :
        leg_x = 0.31
        #leg_y = 0.5
        leg_y = 0.57
        legend_fontsize = 14

    #crTop
    if n_cols == 2 and (region_name == "crTop" or region_name == "vrTop" or region_name == "crZ") and (var_name == "HT2Ratio") :
        leg_x = 0.02
        leg_y = 0.57

    if n_cols == 2 and (region_name == "crTop" or region_name == "vrTop") and (var_name == "dRll" or var_name == "dphi_ll" or var_name == "mt2_bb") :
        leg_x = 0.31
        #leg_y = 0.5
        leg_y = 0.6
        legend_fontsize = 14

    if n_cols == 2 and (region_name == "crZ" or region_name == "vrZ") and (var_name == "mt2_bb") :
        leg_x = 0.31
        #leg_y = 0.5
        leg_y = 0.6
        legend_fontsize = 14

    if (region_name == "crTop") and (var_name == "mt2_bb") :
        n_cols = 1 
        leg_x = 0.58
        leg_y = 0.71

    if (region_name == "crZ") and (var_name == "dRll" or var_name == "dphi_ll") :
        n_cols = 1 
        leg_x = 0.58
        leg_y = 0.71

    if n_cols == 2 and ("vrZ") and (var_name == "HT2Ratio") :
        leg_x = 0.02
        leg_y = 0.57

    if (region_name == "vrZ") and (var_name == "dRll" or var_name == "dphi_ll") :
        n_cols = 1 
        leg_x = 0.58
        leg_y = 0.71
        


    if len(ordered_labels) > 10 :
        leg_y = 0.95 * leg_y
    pad.legend(new_handles,
                new_labels,
                loc = (leg_x, leg_y),
                frameon = False,
                ncol = n_cols,
                fontsize = legend_fontsize,
                numpoints = 1,
                labelspacing = 0.2,
                columnspacing = 0.4)

    return leg_x, leg_y

def colors_dict() :

    palette = 5
    colors = {}
    colors[0] = { "other" : "#F1FAEE",
                    "zjetshf" : "#A8DADC",
                    "top" : "#457B9D",
                    "singlehiggs" : "#E63946",
                    "ttbarv" : "#1D3557" }
    colors[1] = { "top" : "#3E6990",
                    "zjetshf" : "#AABD8C",
                    "ttbarv" : "#381D2A",
                    "singlehiggs" : "#F39B6D",
                    "other" : "#E9E3B4" }
    colors[2] = { "top" : "#5BC0EB",
                    "zjetshf" : "#FDE74C",
                    "ttbarv" : "#9BC53D",
                    "singlehiggs" : "#C3423F",
                    "other" : "#211A1E" }
    #colors[3] = { "top" : "#3F88C5",
    colors[3] = { "top" : "#016FB9",
                    "zjetshf" : "#FF9505",
                    #"zjetshf" : "#44BBA4",
                    "ttbarv" : "#393E41",
                    "singlehiggs" : "#E94F37",
    }
                    #"other" : "#44BBA4" }
                    #"other" : "#F6F7EB" }
    colors[4] = {
        "top" : "#2C7BB6",
        "zjetshf" : "#D7191C",
        "other" : "#953AC8"
        #"other" : "#ABCBFB"
        #"other" : "#FDAE61"
    }

    colors[5] = {
        "top" : "#009FFD",
        #"top" : "#3F88C5",
        "zjetshf" : "#E94F37",
        "other" : "#F7B32B"
        #"other" : "#44BBA4"
    }
                
    return colors[palette]

def get_legend_order(var_name, region_name) :

    standard_order = {}
    standard_order[1] = ["Data", "Top", r'\textit{Z}/\textbf{$\gamma^*$}+jets HF', "Other"]
    standard_order[2] = ["Data", "Top", r'\textit{Z}/\textbf{$\gamma^*$}+jets HF', "$t\\bar{t} + V$", "Total SM", "Higgs", "Other"]

    n_cols = {}
    n_cols["srIncNoDhhNoMbb"] = 1
    n_cols["srIncNoDhh"] = 1
    n_cols["srSFNoDhh"] = 1
    n_cols["srDFNoDhh"] = 1
    n_cols["srIncNoDhhClose"] = 1
    n_cols["srSFNoDhhClose"] = 1
    n_cols["srDFNoDhhClose"] = 1
    n_cols["srSFNoDhhCloseCut"] = 2
    n_cols["srDFNoDhhCloseCut"] = 2
    n_cols["crTopNoDhh"] = 1
    n_cols["crTop"] = 2
    n_cols["vrTop"] = 2
    n_cols["vrZ"] = 2
    n_cols["crZNoDhh"] = 1
    n_cols["crZ"] = 2
    n_cols["srIncNoMbbDhh"] = 1
    n_cols["srIncNoMllDhh"] = 1

    order_dict = {}
    order_dict["NN_d_hh"] = {
        "srIncNoDhh" : standard_order[n_cols[region_name]],
        "srSFNoDhh" : standard_order[n_cols[region_name]],
        "srDFNoDhh" : standard_order[n_cols[region_name]],
        "srIncNoDhhClose" : standard_order[n_cols[region_name]],
        "srSFNoDhhClose" : standard_order[n_cols[region_name]],
        "srDFNoDhhClose" : standard_order[n_cols[region_name]],
        "srSFNoDhhCloseCut" : standard_order[n_cols[region_name]],
        "srDFNoDhhCloseCut" : standard_order[n_cols[region_name]],
        "crTopNoDhh" : standard_order[n_cols[region_name]],
        "crTop" : standard_order[n_cols[region_name]],
        "vrTop" : standard_order[n_cols[region_name]],
        "vrZ" : standard_order[n_cols[region_name]],
        "crZNoDhh" : standard_order[n_cols[region_name]],
        "crZ" : standard_order[n_cols[region_name]],
    }
    order_dict["mbb"] = {
        "srIncNoDhhNoMbb" : standard_order[n_cols[region_name]]
        ,"srIncNoDhh" : standard_order[n_cols[region_name]]
        ,"crTop" : standard_order[n_cols[region_name]]
        ,"vrTop" : standard_order[n_cols[region_name]]
        ,"crZ" : standard_order[n_cols[region_name]]
        ,"vrZ" : standard_order[n_cols[region_name]]
        ,"srIncNoMbbDhh" : standard_order[n_cols[region_name]]
    }
    order_dict["HT2Ratio"] = {
        "srIncNoDhh" : standard_order[n_cols[region_name]]
        ,"crTop" : standard_order[n_cols[region_name]]
        ,"vrTop" : standard_order[n_cols[region_name]]
        ,"crZ" : standard_order[n_cols[region_name]]
        ,"vrZ" : standard_order[n_cols[region_name]]
        ,"srIncNoMbbDhh" : standard_order[n_cols[region_name]]
    }
    order_dict["dRll"] = {
        "srIncNoDhh" : standard_order[n_cols[region_name]]
        ,"crTop" : standard_order[n_cols[region_name]]
        ,"vrTop" : standard_order[n_cols[region_name]]
        ,"crZ" : standard_order[n_cols[region_name]]
        ,"vrZ" : standard_order[n_cols[region_name]]
    }
    order_dict["mt2_bb"] = {
        "srIncNoDhh" : standard_order[n_cols[region_name]]
        ,"crTop" : standard_order[n_cols[region_name]]
        ,"vrTop" : standard_order[n_cols[region_name]]
        ,"crZ" : standard_order[n_cols[region_name]]
        ,"vrZ" : standard_order[n_cols[region_name]]
    }
    order_dict["dphi_ll"] = {
        "srIncNoDhh" : standard_order[n_cols[region_name]]
        ,"crTop" : standard_order[n_cols[region_name]]
        ,"vrTop" : standard_order[n_cols[region_name]]
        ,"crZ" : standard_order[n_cols[region_name]]
        ,"vrZ" : standard_order[n_cols[region_name]]
    }
    order_dict["mll"] = {
        "srIncNoMllDhh" : standard_order[n_cols[region_name]]
    }

    return order_dict[var_name][region_name], n_cols[region_name]

class OOMFormatter(ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

def make_paper_plot(region, backgrounds, signals, data, plot_description, args) :

    histograms_bkg = {}
    labels_bkg = {}
    colors_bkg = {}
    top_bkg = ["ttbar", "wt"]
    other_bkg = ["zjetslf", "diboson", "fakest3", "singlehiggs", "ttbarv"]
    labels_bkg["other"] = "Other"
    colors_bkg["other"] = "g"
    labels_bkg["top"] = "Top"
    colors_bkg["top"] = "#057390"


    x_bounds = bounds_dict()[plot_description.var_to_plot][region.name]
    is_variable_width = len(x_bounds) > 3
    bin_width = x_bounds[0]
    hatch_bin_width = 1.0 * bin_width
    x_lo = x_bounds[1]
    x_hi = x_bounds[2]
    if is_variable_width :
        x_lo = x_bounds[0]
        x_hi = x_bounds[-1]
        hatch_bin_width = 1.0 * np.array(x_bounds)[1:] - 1.0 * np.array(x_bounds)[0:-1]

    ##
    ## Fill background histograms
    ##

    h_other = histogram1d("histo_other", binning = x_bounds)
    h_top = histogram1d("histo_top", binning = x_bounds)
    for ibkg, bkg in enumerate(backgrounds) :

        # initialize the histogram
        h = histogram1d("histo_%s" % (bkg.name), binning = x_bounds)
        if bkg.name.lower() not in other_bkg :
            labels_bkg[bkg.name] = bkg.displayname
            colors_bkg[bkg.name.lower()] = bkg.color

        # start filling
        chain = bkg.chain()
        weight_str = "eventweightNoPRW_multi"
        for ic, c in enumerate(chain) :

            # project out flavor component of Z+jets
            if "zjets" in bkg.name.lower() :
                is_bb = c["isBB"]
                is_cc = c["isCC"]
                is_bc = c["isBC"]
                is_bl = c["isBL"]
                is_cl = c["isCL"]
                is_ll = c["isLL"]
                if "hf" in bkg.name.lower() :
                    idx_hf = (is_bb == 1) | (is_cc == 1) | (is_bc == 1)
                    c = c[idx_hf]
                elif "lf" in bkg.name.lower() :
                    idx_lf = (is_bl == 1) | (is_cl == 1) | (is_ll == 1)
                    c = c[idx_lf]

            trigger_idx = get_trigger_idx(c)
            c = c[trigger_idx]
            if plot_description.region_name == "vrZ" :
                vrz_idx = get_vrz_mll_cut(c)
                c = c[vrz_idx]
            weights = c[weight_str]

            if plot_description.region_name == "crTopNoDhh" and bkg.name.lower() == "zjetslf" :
                #print "zlf weights: %s" % weights
                idx_good_weight =  weights > 1e-2
                c = c[idx_good_weight]
                weights = weights[idx_good_weight]

            if plot_description.is_abs :
                plot_data = np.absolute(c[plot_description.var_to_plot])
            else :
                plot_data = c[plot_description.var_to_plot]
            lumis = bkg.scalefactor * np.ones(len(plot_data))

            #if bkg.name.lower() == "ttbar" and ('_multi' in weight_str and 'NoPRW' not in weight_str) :
            #    print "******** WARNING APPLYING PRW SUMW CORRECTION FACTORS TO SAMPLES %s ********" % bkg.name
            #    idx_a = (c['year'] == 2015) | (c['year'] == 2016)
            #    idx_d = (c['year'] == 2017)
            #    idx_e = (c['year'] == 2018)
            #    lumis[idx_a] *= 0.93
            #    lumis[idx_d] *= 0.61
            #    lumis[idx_e] *= 1.23

            if args.apply_sr_sf  :
                dhh_cut = 1.0
                if plot_description.region_name == "srSFNoDhhCloseCut" :
                    dhh_cut = 5.45
                elif plot_description.region_name == "srDFNoDhhCloseCut" :
                    dhh_cut = 5.55
                elif "crz" in plot_description.region_name.lower() :
                    dhh_cut = 0
                elif plot_description.region_name.lower() == "crtop" :
                    dhh_cut = 4.5
                elif (plot_description.region_name == "srIncNoDhh" or plot_description.region_name == "srIncNoMbbDhh") :
                    dhh_cut = 5.45
                elif (plot_description.region_name == "srIncNoMllDhh") :
                    dhh_cut = 5.45
                elif (plot_description.region_name == "srSFNoDhh") :
                    dhh_cut = 5.45
                elif (plot_description.region_name == "srDFNoDhh") :
                    dhh_cut = 5.55
                elif plot_description.region_name == "srIncNoDhhNoMbb" :
                    dhh_cut = 5.45
                if bkg.name.lower() == "ttbar" or bkg.name.lower() == "wt" or "top" in bkg.name.lower() :
                    scales = np.ones(len(plot_data))
                    dhh = c['NN_d_hh']
                    sr_cut_idx = dhh > dhh_cut
                    lumis[sr_cut_idx] *= 0.79
                    print " *** SCALING TTBAR IN SR ONLY ***"
                    
                if "zjets" in bkg.name.lower() and "hf" in bkg.name.lower() :
                    scales = np.ones(len(plot_data))
                    dhh = c['NN_d_hh']
                    sr_cut_idx = dhh > dhh_cut
                    lumis[sr_cut_idx] *= 1.34
                    print " *** SCALING ZJETS IN SR ONLY ***"

            if "fakes" in bkg.name.lower() :
                rname = plot_description.region_name
                fake_scaling = 1.0
                if "srdf" in rname.lower() :
                    fake_scaling = 1.1798
                elif "srsf" in rname.lower() :
                    fake_scaling = 1.1300
                elif "srinc" in rname.lower() :
                    fake_scaling = 1.1513
                elif "crt" in rname.lower() :
                    #fake_scaling = 1.476 * 2.47
                    #print "SCALING FAKES FUCK"
                    fake_scaling = 1.476
                elif "crz" in rname.lower() :
                    fake_scaling = 0.8299

                scales = np.ones(len(plot_data))
                lumis *= fake_scaling

            weights = lumis * weights
            if bkg.name.lower() not in other_bkg and bkg.name.lower() not in top_bkg :
                h.fill(plot_data, weights)
            elif bkg.name.lower() in other_bkg :
                if "fakes" in bkg.name.lower() :
                    print "ADDING FAKES FUCK HISTO"
                h_other.fill(plot_data, weights)
            elif bkg.name.lower() in top_bkg :
                h_top.fill(plot_data, weights)
                h.fill(plot_data, weights)

        if bkg.name.lower() not in other_bkg :
            histograms_bkg[bkg.name] = h

    # adjust the colors
    for name in colors_dict() :
        colors_bkg[name] = colors_dict()[name]

    histograms_bkg["other"] = h_other
    histograms_bkg["top"] = h_top

    add_overflow = False
    if add_overflow :
        for bkg in histograms_bkg :
            histograms_bkg[bkg].add_overflow()

    stack = histogram_stack("bkg_stack_%s" % plot_description.var_to_plot, binning = x_bounds)
    ordered_bkg_labels = []
    ordered_bkg_colors = []
    stack_order = ["top", "zjetshf", "ttbarv", "higgs", "other"]
#    for bkg in backgrounds :
#            if bkg.name.lower() in other_bkg : continue
#            if bkg.name.lower() in top_bkg : continue
#            stack.add(histograms_bkg[bkg.name])
    #stack.add(histograms_bkg["other"])
    for so in stack_order[::-1] :
        for bkg_name in histograms_bkg :
            if so in bkg_name.lower() :
                stack.add(histograms_bkg[bkg_name])
                break
#    stack.add(histograms_bkg["other"])
#    stack.add(histograms_bkg["top"])
    #stack.sort(reverse = True)
    for bkg_name in stack.order :
        name = bkg_name.replace("histo_", "")
        name = name.split("_")[0]
        ordered_bkg_labels.append(labels_bkg[name])
        ordered_bkg_colors.append(colors_bkg[name.lower()])

    ##
    ## Fill data histogram
    ##

    h_data = histogram1d("histo_data_%s" % plot_description.var_to_plot, binning = x_bounds)
    chain = data.chain()
    for idc, dc in enumerate(chain) :
        trigger_idx = get_trigger_idx(dc)
        dc = dc[trigger_idx]
        if plot_description.region_name == "vrZ" :
            vrz_idx = get_vrz_mll_cut(dc)
            dc = dc[vrz_idx]
        if plot_description.is_abs :
            plot_data = np.absolute(dc[plot_description.var_to_plot])
        else :
            plot_data = dc[plot_description.var_to_plot]
        h_data.fill(plot_data)
    if add_overflow :
        h_data.add_overflow()

    ##
    ## start drawing
    ##

    canvas = ratio_canvas("ratio_canvas_%s" % plot_description.var_to_plot, logy = args.logy)
    x_label = nice_names_dict()[plot_description.var_to_plot]
    bw_str = "%.1f" % bin_width
    if bw_str.strip().split(".")[-1] == "0" :
        bw_str = bw_str.strip().split(".")[0]
    y_label = "Events / %s" % bw_str
    if "GeV" in x_label :
        y_label += " GeV"
    elif "rad" in x_label.lower() :
        y_label += " rad."
    else :
        y_label = "Events / %s" % bw_str
    #y_label = "Events / Bin"

    #if plot_description.var_to_plot == "dRll" and plot_description.region_name == "srIncNoDhh" :
    #    y_label = "Events / 0.2"
    #elif plot_description.var_to_plot == "HT2Ratio" and plot_description.region_name == "srIncNoDhh" :
    #    y_label = "Events / 0.1"
    #elif plot_description.var_to_plot == "mbb" and plot_description.region_name = "srIncNoMbbDhh" :
    #    y_label = "Events / 20 GeV"
    #elif plot_description.var_to_plot == "HT2Ratio" and plot_description.region_name = "srIncNoMbbDhh" :
    #    y_label = "Events / 0.1"

    canvas.labels = [x_label, y_label]
    canvas.x_bounds = [x_lo, x_hi]
    canvas.build()

    upper_pad = canvas.upper_pad
    lower_pad = canvas.lower_pad
    
    if not is_variable_width :
        ticks = list(np.arange(x_lo, x_hi + bin_width, bin_width))
        ticks_plot = ticks
        if len(ticks) > 10 :
            ticks_plot = ticks[::2]
        upper_pad.set_xticks(ticks_plot)
        lower_pad.set_xticks(ticks_plot)
    else :
        ticks = x_bounds

    

    histo_total = stack.total_histo
    maxy = histo_total.maximum()
    miny = 0.0
    if args.logy :
        miny = 1e-1
    multiplier = 1.75
    if len(signals) :
        multiplier = 1.8
    if args.logy :
        multiplier = 1e4
        if "srsf" in plot_description.region_name.lower() or "srdf" in plot_description.region_name.lower() :
            multiplier = 1e2
    maxy = multiplier * maxy

    if plot_description.region_name == "srSFNoDhhCloseCut" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 10
    elif plot_description.region_name == "srDFNoDhhCloseCut" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 10
    elif plot_description.region_name == "crZ" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 400
    elif plot_description.region_name == "crTop" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 70
    elif plot_description.region_name == "vrTop" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 85
    elif plot_description.region_name == "vrZ" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 95
    elif plot_description.region_name == "srIncNoDhh" and plot_description.var_to_plot == "HT2Ratio" and args.logy :
        maxy = 1e8
    elif plot_description.region_name == "srIncNoDhh" and plot_description.var_to_plot == "dRll" and args.logy :
        maxy = 1e7
    elif plot_description.region_name == "srDFNoDhh" and plot_description.var_to_plot == "NN_d_hh" and args.logy :
        maxy = 5e6
    elif plot_description.region_name == "srSFNoDhh" and plot_description.var_to_plot == "NN_d_hh" and args.logy :
        maxy = 1e7
    elif plot_description.region_name == "srIncNoMbbDhh" and plot_description.var_to_plot == "mbb" and not args.logy :
        maxy = 50
    elif plot_description.region_name == "srIncNoMbbDhh" and plot_description.var_to_plot == "HT2Ratio" and not args.logy :
        maxy = 260
    elif plot_description.region_name == "srIncNoMllDhh" and plot_description.var_to_plot == "mll" and not args.logy :
        maxy = 18
    elif plot_description.region_name == "crTop" and plot_description.var_to_plot == "HT2Ratio" and not args.logy :
        maxy = 70
    elif plot_description.region_name == "crTop" and plot_description.var_to_plot == "dRll" and not args.logy :
        maxy = 42
    elif plot_description.region_name == "crTop" and plot_description.var_to_plot == "mt2_bb" and not args.logy :
        maxy = 40
    elif plot_description.region_name == "crTop" and plot_description.var_to_plot == "dphi_ll" and not args.logy :
        maxy = 42
    elif plot_description.region_name == "vrTop" and plot_description.var_to_plot == "HT2Ratio" and not args.logy :
        maxy = 87
    elif plot_description.region_name == "vrTop" and plot_description.var_to_plot == "dRll" and not args.logy :
        maxy = 70
    elif plot_description.region_name == "vrTop" and plot_description.var_to_plot == "mt2_bb" and not args.logy :
        maxy = 65
    elif plot_description.region_name == "vrTop" and plot_description.var_to_plot == "dphi_ll" and not args.logy :
        maxy = 78
    elif plot_description.region_name == "crZ" and plot_description.var_to_plot == "HT2Ratio" and not args.logy :
        maxy = 570
    elif plot_description.region_name == "crZ" and plot_description.var_to_plot == "dRll" and not args.logy :
        maxy = 670
    elif plot_description.region_name == "crZ" and plot_description.var_to_plot == "mt2_bb" and not args.logy :
        maxy = 350
    elif plot_description.region_name == "crZ" and plot_description.var_to_plot == "dphi_ll" and not args.logy :
        maxy = 285
    elif plot_description.region_name == "vrZ" and plot_description.var_to_plot == "HT2Ratio" and not args.logy :
        maxy = 100
    elif plot_description.region_name == "vrZ" and plot_description.var_to_plot == "mt2_bb" and not args.logy :
        maxy = 46
    elif plot_description.region_name == "vrZ" and plot_description.var_to_plot == "dphi_ll" and not args.logy :
        maxy = 48



    # stat error
    sm_x_error = np.zeros(len(histo_total.y_error()))
    sm_y_error = histo_total.y_error()

    if not args.add_syst :
        stat_error_band = errorbars.error_hatches(histo_total.bins[:-1], histo_total.histogram, \
                sm_x_error, sm_y_error, hatch_bin_width)

    # total SM line
    sm_line = histo_total.bounding_line()

    # print yields
    stack.print_counts()
    histos = []
    weights = []
    for name in stack.order :
        for h in stack.histograms :
            if name != h.name.replace("hist_", "") : continue
            histos.append(h.data)
            weights.append(h.weights)
    bin_vals = ticks
    if is_variable_width :
        bin_vals = x_bounds
    upper_pad.hist( histos,
                    weights = weights,
                    bins = bin_vals,
                    color = ordered_bkg_colors,
                    label = ordered_bkg_labels,
                    stacked = True,
                    histtype = "stepfilled",
                    lw = 0.5,
                    edgecolor = "k",
                    alpha = 1.0
    )

    # print yields
    sm_total_yield = histo_total.integral()
    print 30 * "*"
    print histo_total.count_str(name = "Total SM")

    # draw error band
    if not args.add_syst :
        upper_pad.add_collection(stat_error_band)

    # draw total SM
    #upper_pad.plot(sm_line[0], sm_line[1], ls = "-", color = "k", label = "Total SM", lw = 2)

    # draw data
    if not is_variable_width :
        data_x = np.array(h_data.bin_centers())
    else :
        x_left_edges = 1.0 * np.array(x_bounds[:-1])
        x_left_edges += 0.5 * hatch_bin_width
        data_x = x_left_edges
    data_y = h_data.histogram
    data_y[data_y == 0.] = -5
    upper_pad.plot(data_x, data_y, "ko", label = "Data", zorder = 1e6, markersize = 8)

    # poisson errors on data
    data_err_low, data_err_high = errorbars.poisson_interval(data_y)
    data_err_low = data_y - data_err_low
    data_err_high = data_err_high - data_y
    data_err = [data_err_low, data_err_high]
    upper_pad.errorbar(data_x, data_y, yerr = data_err, fmt = "none", color = "k", zorder = 1e6)

    print h_data.count_str(name = "Data")
    print " > Data / SM : %5.2f" % ( h_data.integral() / sm_total_yield )

    ##
    ## signal histograms
    ##
    #signal_color = "#29fd2f" # color used for first submission
    #signal_color = "#00ffff" # bright blue
    #signal_color = "#6fff00" # neon green
    #signal_color = "r"
    #signal_color = "#ffff00" # bright yellow
    #signal_color = "#33cc33" # green
    #signal_color = "k"
    signal_color = "#ff00ff"#
    if len(signals) > 0 and "sr" in plot_description.region_name.lower() :
        for signal in signals :
            h = histogram1d("signal_histo_%s" % signal.name, binning = x_bounds)
            chain = signal.chain()
            for isc, sc in enumerate(chain) :
                trigger_idx = get_trigger_idx(sc)
                sc = sc[trigger_idx]

                if plot_description.region_name == "vrZ" :
                    vrz_idx = get_vrz_mll_cut(sc)
                    sc = sc[vrz_idx]

                weights = sc[weight_str]

                if plot_description.is_abs :
                    plot_data = np.absolute(sc[plot_description.var_to_plot])
                else :
                    plot_data = sc[plot_description.var_to_plot]
                lumis = signal.scalefactor * np.ones(len(plot_data))

                weights = lumis * weights
                #if "inc" in plot_description.region_name.lower() :
                #    weights *= 50
                h.fill(plot_data, weights)

            # signal color
            upper_pad.hist( h.data, weights = h.weights,
                        bins = bin_vals,
                        #color = 'k', #signal.color,
                        #color = "#72ff02",
                        #color = "#29fd2f", # color used for first submission
                        color = signal_color,
                        label = "$hh$",
                        #ls = ":",
                        ls = "-",
                        stacked = False,
                        histtype = "step",
                        lw = 2,
                        zorder = 1e5
            )
            print 35 * "- "
            print "hh counts"
            print h.count_str()


    ##
    ## ratio
    ##
    ratio_y_min = 0.5
    ratio_y_max = 1.5
    lower_pad.set_yticks([0.5, 0.75, 1.0, 1.25, 1.5])
    lower_pad.set_ylim([ratio_y_min, ratio_y_max])

    upper_pad.set_ylim(miny, maxy)
    upper_pad.set_xlim(x_lo, x_hi)
    lower_pad.set_xlim(x_lo, x_hi)

    # helvetic x and y axis labels
    y_tick_loc = upper_pad.get_yticks()
    y_tick_loc = [x for x in y_tick_loc if x >= miny]
    y_tick_loc = [x for x in y_tick_loc if x <= maxy]
    y_tick_labs = ["%s" % x for x in upper_pad.get_yticks() ]
    y_tick_labs = [x.replace(".0","").replace(".00","") for x in y_tick_labs]
    upper_pad.set_yticks(y_tick_loc)
    upper_pad.set_yticklabels(y_tick_labs)
    if args.logy :
        #print "FOOBS %s" % ["{:.0e}".format(x) for x in upper_pad.get_yticks()]
        #upper_pad.yaxis.set_major_formatter(ScalarFormatter())
        #upper_pad.yaxis.set_major_formatter(OOMFormatter(9, "%1.1f"))
        #upper_pad.ticklabel_format(axis = "y", style = "scientific", scilimits=(2,2))
        y_tick_labs = ["{:.2e}".format(x) for x in upper_pad.get_yticks()]
        y_tick_labs = [r"10$^{\mbox{%s}}$" % str(int(x.split("e")[-1].replace("+",""))) for x in y_tick_labs]
        upper_pad.set_yticklabels(y_tick_labs)
        #y_tick_labs = ["{:.2e}".format(x) for x in upper_pad.get_yticks() ]

    

    r_tick_loc = lower_pad.get_yticks()
    r_tick_labs = ["%s" % x for x in lower_pad.get_yticks() ]
    lower_pad.set_yticks(r_tick_loc)
    lower_pad.set_yticklabels(r_tick_labs)

    #x_tick_loc = lower_pad.get_xticks()
    bw, xlo, xhi = bounds_dict()[plot_description.var_to_plot][plot_description.region_name]
    #print "FUCK bw xlo xhi = %s %s %s" % (bw, xlo, xhi)
    x_tick_loc = [round(x, 2) for x in np.arange(x_lo, xhi+bw, bw) if x <= xhi]
    x_tick_loc_orig = [round(x, 2) for x in np.arange(x_lo, xhi+bw, bw) if x <= xhi]
    x_tick_labs = ["%.2f" % x for x in x_tick_loc] # if x < (xhi + 1e-5)]
    if (plot_description.var_to_plot == "NN_d_hh" or plot_description.var_to_plot == "mt2_bb" or plot_description.var_to_plot == "mbb" or plot_description.var_to_plot == "mll") and ("srInc" in plot_description.region_name or "srDF" in plot_description.region_name or "srSF" in plot_description.region_name) and ("Cut" not in plot_description.region_name) : 
        x_tick_labs = ["%d" % int(x) for x in x_tick_loc] # if x < (xhi + 1e-5)]
    x_tick_labs_orig = ["%.2f" % x for x in x_tick_loc] # if x < (xhi + 1e-5)]
    if len(x_tick_labs) > 8 :
        x_tick_loc = x_tick_loc[::2]
        x_tick_labs = x_tick_labs[::2]
    #print "FUCK3 %s" % x_tick_labs_orig[-1]
    #if x_tick_labs_orig[-1] not in x_tick_labs :
    #    x_tick_loc.append(x_tick_loc_orig[-1])
    #    x_tick_labs.append(x_tick_labs_orig[-1])
    #print "FUCK loc %s" % x_tick_loc
    #print "FUCK lab %s" % x_tick_labs
    lower_pad.set_xlim([xlo,xhi])
    upper_pad.set_xlim([xlo,xhi])
    
    lower_pad.set_xticks(x_tick_loc)
    lower_pad.set_xticklabels(x_tick_labs)

    x_tick_labs = ["%s" % x for x in lower_pad.get_xticks()]
    x_tick_labs = [x.replace(".00","").replace(".50",".5").replace(".20",".2") for x in x_tick_labs]
    all_zeros = True
    found_split = False
    for x in x_tick_labs :
        x_split = x.split(".")
        if len(x_split) == 2 :
            found_split = True
            if float(x_split[1]) != 0.0 :
                all_zeros = False
    if found_split and all_zeros :
        x_tick_labs = [x.split(".")[0] for x in x_tick_labs]
                
            
    print "x_tick_labs = %s" % x_tick_labs
    x_tick_loc = lower_pad.get_xticks()
    lower_pad.set_xticks(x_tick_loc)
    lower_pad.set_xticklabels(x_tick_labs)
    upper_pad.set_xticks(x_tick_loc)

    pred_y = histo_total.histogram
    ratio_y = h_data.divide(histo_total)
    ratio_y[ ratio_y == 0. ] = -1
    ratio_x = np.array(h_data.bin_centers())

    ratio_data_err_low = -1 * np.ones(len(ratio_y))
    ratio_data_err_high = -1 * np.ones(len(ratio_y))
    for idata, d in enumerate(ratio_y) :
        prediction = pred_y[idata]
        if ratio_y[idata] == 0.0 or ratio_y[idata] < 0 :
            ratio_data_err_low[idata] = 0
            ratio_data_err_high[idata] = 0
        else :
            ratio_data_err_low[idata] = data_err_low[idata] / prediction
            ratio_data_err_high[idata] = data_err_high[idata] / prediction
    lower_pad.plot(ratio_x, ratio_y, "ko", zorder = 1000, markersize = 8)
    yerr = [ratio_data_err_low, ratio_data_err_high]
    lower_pad.errorbar(ratio_x, ratio_y, yerr = yerr, fmt = "none", color = "k")

    ##
    ## draw arrows for points outside of the y-axis range on the ratio
    ##

    # arrow props

    # get indices for points outside y-axix
    idx_ratio_up = ((ratio_y - ratio_data_err_low) > ratio_y_max) & (ratio_y > 0)
    idx_ratio_dn = ((ratio_y + ratio_data_err_high) < ratio_y_min) & (ratio_y > 0)
    x_vals_arrow_up = np.extract(idx_ratio_up, ratio_x)
    x_vals_arrow_dn = np.extract(idx_ratio_dn, ratio_x)


    arrow_transform = lower_pad.transAxes
    for direction in [-1.0, 1.0] :
        x_arrow = { -1.0 : x_vals_arrow_dn,
                    1.0 : x_vals_arrow_up }[direction]
        y_ratio_val = { -1.0 : ratio_y_min,
                        1.0 : ratio_y_max }[direction]
        arrow_length = 0.2
        length = direction * arrow_length
        head_width = 0.012 * (abs(x_hi - x_lo))
        head_length = 0.045 * (abs(ratio_y_max - ratio_y_min))
        y_arrow_pos = y_ratio_val + (-1.0 * direction) * 1.38 * arrow_length
        for ix, x_arrow_pos in enumerate(x_arrow) :
            arrow_opts = {"head_width" : head_width, "head_length" : head_length, "edgecolor" : "r", "facecolor" : "r", "linewidth" : 2.5, "zorder" : 1e9 }
            lower_pad.arrow(x_arrow_pos, y_arrow_pos, 0, length, **arrow_opts)

    ##
    ## SM error on ratio
    ##
    sm_ratio_err = []
    for ism, sm in enumerate(pred_y) :
        sm_y_error_ratio = sm_y_error[ism]
        relative_error = 0.0
        if sm != 0 :
            relative_error = float(sm_y_error_ratio) / float(sm)
        if ratio_y[ism] == 0 :
            relative_error = 0
        sm_ratio_err.append(relative_error)
    sm_x_error_ratio = [bin_width for a in ratio_x]
    if is_variable_width :
        sm_x_error_ratio = []
        for ix, a in enumerate(ratio_x) :
            sm_x_error_ratio.append( hatch_bin_width[ix] )

    # syst
    if args.add_syst and is_variable_width :
        print "ERROR Systematics and variable bin widths are not yet implmeneted"
        sys.exit()

    if args.add_syst :
        syst_config_dir = "/data/uclhc/uci/user/dantrim/n0307val/dantrimania/python/analysis/wwbb/error_band/"
        syst_files = {}

        if plot_description.var_to_plot == "NN_d_hh" :
            #syst_files = { "crTop" : "%s/process_uncerts_NN_d_hh_crTop_paper.json" % syst_config_dir,
            #                "crZ" : "%s/process_uncerts_NN_d_hh_crZ_paper.json" % syst_config_dir,
            #                "vrTop" : "%s/process_uncerts_NN_d_hh_top_vr_paper.json" % syst_config_dir,
            #                "vrZ" : "%s/process_uncerts_NN_d_hh_z_vr_paper.json" % syst_config_dir,
            #                "srIncNoDhh" : "%s/process_uncerts_NN_d_hh_srIncNoDhh_paper.json" % syst_config_dir,
            #                "srSFNoDhhCloseCut" : "%s/process_uncerts_NN_d_hh_srSFNoDhhCloseCut_paper.json" % syst_config_dir,
            #                "srDFNoDhhCloseCut" : "%s/process_uncerts_NN_d_hh_srDFNoDhhCloseCut_paper.json" % syst_config_dir

            #}
            syst_files = {
                "crTop" : "process_uncerts_NN_d_hh_top_cr_paper_may29.json"
                ,"vrTop" : "process_uncerts_NN_d_hh_top_vr_paper_may29.json"
                ,"crZ" : "process_uncerts_NN_d_hh_z_cr_paper_may29.json"
                ,"vrZ" : "process_uncerts_NN_d_hh_z_vr_paper_may29.json"
                ,"srIncNoDhh" : "process_uncerts_NN_d_hh_srIncNoDhh_paper_may29.json"
                ,"srSFNoDhh" : "process_uncerts_NN_d_hh_srSFNoDhh_paper_may29.json"
                ,"srDFNoDhh" : "process_uncerts_NN_d_hh_srDFNoDhh_paper_may29.json"
                ,"srSFNoDhhCloseCut" : "process_uncerts_NN_d_hh_srSFNoDhhCloseCut_paper_may29.json"
                ,"srDFNoDhhCloseCut" : "process_uncerts_NN_d_hh_srDFNoDhhCloseCut_paper_may29.json"
            }

        if plot_description.var_to_plot == "HT2Ratio" :
            #syst_files = {
            #    "srIncNoDhh" : "%s/process_uncerts_HT2Ratio_srIncNoDhh_paper.json" % syst_config_dir
            #}
            syst_files = {
                "crTop" :   "process_uncerts_HT2Ratio_top_cr_paper_may29.json"
                ,"vrTop" :  "process_uncerts_HT2Ratio_top_vr_paper_may29.json"
                #,"crZ" :    "process_uncerts_HT2Ratio_z_cr_paper_may29.json"
                ,"crZ" :    "process_uncerts_HT2Ratio_z_cr_paper_test.json"
                ,"vrZ" :    "process_uncerts_HT2Ratio_z_vr_paper_may29.json"
                ,"srIncNoDhh" : "process_uncerts_HT2Ratio_srIncNoDhh_paper_may29.json"
                ,"srIncNoMbbDhh" : "process_uncerts_HT2Ratio_srIncNoMbbDhh_paper_may29.json"
            }

        if plot_description.var_to_plot == "dRll" :
            #syst_files = {
            #    "srIncNoDhh" : "%s/process_uncerts_dRll_srIncNoDhh_paper.json" % syst_config_dir
            #}
            syst_files = {
                "crTop" :       "process_uncerts_dRll_top_cr_paper_may29.json"
                ,"vrTop" :      "process_uncerts_dRll_top_vr_paper_may29.json"
                ,"crZ" :        "process_uncerts_dRll_z_cr_paper_may29.json"
                ,"vrZ" :        "process_uncerts_dRll_z_vr_paper_may29.json"
                ,"srIncNoDhh" : "process_uncerts_dRll_srIncNoDhh_paper_may29.json"
                ,
            }

        if plot_description.var_to_plot == "mt2_bb" :
            syst_files = {
                "crTop" :       "process_uncerts_mt2_bb_top_cr_paper_may29.json"
                ,"vrTop" :      "process_uncerts_mt2_bb_top_vr_paper_may29.json"
                ,"crZ" :        "process_uncerts_mt2_bb_z_cr_paper_may29.json"
                ,"vrZ" :        "process_uncerts_mt2_bb_z_vr_paper_may29.json"
                ,"srIncNoDhh" : "process_uncerts_mt2_bb_srIncNoDhh_paper_may29.json"
                ,
            }

        if plot_description.var_to_plot == "dphi_ll" :
            syst_files = {
                "crTop" :       "process_uncerts_dphi_ll_top_cr_paper_may29.json"
                ,"vrTop" :      "process_uncerts_dphi_ll_top_vr_paper_may29.json"
                ,"crZ" :        "process_uncerts_dphi_ll_z_cr_paper_may29.json"
                ,"vrZ" :        "process_uncerts_dphi_ll_z_vr_paper_may29.json"
                ,"srIncNoDhh" : "process_uncerts_dphi_ll_srIncNoDhh_paper_may29.json"
                ,
            }
        if plot_description.var_to_plot == "mbb" :
            syst_files = {
                "srIncNoMbbDhh" : "process_uncerts_mbb_srIncNoMbbDhh_paper_may29.json"
                ,
            }
        if plot_description.var_to_plot == "mll" :
            syst_files = {
                "srIncNoMllDhh" : "process_uncerts_mll_srIncNoMllDhh_paper_may29.json"
            }

        syst_proc_names = { "TTbar" : "ttbar", "Wt" : "wt", "ZjetsHF" : "z" }
        process_fractions = {}
        for process in syst_proc_names :
            process_fractions[process] = np.array(histograms_bkg[process].histogram / histo_total.histogram)
            neg_frac = process_fractions[process] < 0
            process_fractions[process][neg_frac] = 0

        if plot_description.region_name in syst_files :
            json_file = "%s/%s" % (syst_config_dir, syst_files[plot_description.region_name])
            with open(json_file, "r") as input_file :
                syst_data = json.load(input_file)
            process_errors = {}
            if syst_data["variable"] == plot_description.var_to_plot :
                processes = syst_data["processes"]
                for process in syst_proc_names :
                    for process_err_group in processes :
                        if process_err_group["name"] == syst_proc_names[process] :
                            #if process_err_group["name"] == "ttbar" :
                            #    print("SKIPPING TTBAR SYSTEMATICS IN ERROR BAR")
                            #    continue
                            err = process_err_group["errors"]
                            if process_err_group["name"] == "z" :
                                print("Z err: %s" % err)
                            if plot_description.var_to_plot == "dphi_ll" and plot_description.region_name == "srIncNoDhh" :
                                print("WARNING Only using scale factor  uncertainties!")
                                err = np.zeros(len(histo_total.histogram))

                            # add Top and Z+HF mu uncerts
                            if process_err_group["name"] == "ttbar" or process_err_group["name"] == "wt" :
                                mu_top_uncert = 0.098 * np.ones(len(err))
                                err = np.sqrt(np.power(err,2) + np.power(mu_top_uncert,2))
                            elif process_err_group["name"] == "z" :
                                mu_z_uncert = 0.066 * np.ones(len(err))
                                err = np.sqrt(np.power(err,2) + np.power(mu_z_uncert,2))
                            err_times_frac = err * process_fractions[process]
                            process_errors[process] = err_times_frac
            total_err = np.zeros(len(histo_total.histogram))
            for proc_name, proc_err in process_errors.iteritems() :
                total_err += np.power(proc_err,2)
            total_err = np.sqrt(total_err)

            # this is the fraction of the already-present error that the sytematic error is
            syst_err_factor = total_err / sm_ratio_err
            sm_ratio_err = np.sqrt( np.power(sm_ratio_err,2) + np.power(total_err,2) )

            ##
            ## update the error band in the upper-pad
            ##

            # update the already-calculated stat-only errors
            sm_y_error *= syst_err_factor
            hist_error_band = errorbars.error_hatches(histo_total.bins[:-1], histo_total.histogram, \
                sm_x_error, sm_y_error, hatch_bin_width)
            upper_pad.add_collection(hist_error_band)

    ratio_error_band = errorbars.error_hatches(
        [xv - 0.5 * bin_width for xv in ratio_x],
        np.ones(len(ratio_y)),
        sm_x_error_ratio,
        sm_ratio_err,
        hatch_bin_width
    )
    lower_pad.add_collection(ratio_error_band)

    # line at unity
    lower_pad.plot([x_lo, x_hi], [1.0, 1.0], 'r--', lw = 1, alpha = 0.5)
    
    #upper_pad.get_yaxis().set_label_coords(-0.16, 1.0)
    #lower_pad.get_yaxis().set_label_coords(-0.16, 0.5)
    upper_pad.get_yaxis().set_label_coords(-0.11, 1.0)
    lower_pad.get_yaxis().set_label_coords(-0.11, 0.5)

    ##
    ## legend
    ##
    legend_order, n_cols = get_legend_order(plot_description.var_to_plot, plot_description.region_name)
    leg_x, leg_y = make_legend(legend_order, n_cols, plot_description.var_to_plot, plot_description.region_name, upper_pad)

    ##
    ## signal legend
    ##
    if len(signals) > 0 and "sr" in plot_description.region_name :
            # add signal to legend (here we assume that there is only one signal sample provided)

            # label
            signal_label = "\\textit{HH} ($\\times$20)" 
            #signal_label = "$hh \\rightarrow b\\bar{b}\\ell\\nu\\ell\\nu$" # ($\\times 20$)"
            y_text = leg_y - 0.04
            #x_text = leg_x + 0.15
            x_text = leg_x + 0.12
            upper_pad.text(x_text, y_text, signal_label,
                        transform = upper_pad.transAxes,
                        fontsize = 14,
                        ha = "left"
            )
            #upper_pad.text(x_text + 0.21, y_text + 0.0025, "($\\times 20$)",
            #        transform = upper_pad.transAxes,
            #        fontsize = 12, ha = "left"
            #)

            # line
            y_line = 1.03 * y_text
            #x_line_0 = leg_x + 0.007 
            #x_line_1 = leg_x + 0.08
            x_line_0 = leg_x + 0.0155
            x_line_1 = leg_x + 0.09

            if n_cols == 2 :
                y_line = 1.02 * y_text
                x_line_0 = leg_x + 0.016
                x_line_1 = leg_x + 0.083
            upper_pad.plot([x_line_0, x_line_1], [y_line, y_line],
                    #":",
                    "-",
                    lw = 2,
                    #color = "k", #signals[0].color,
                    color = signal_color,
                    transform = upper_pad.transAxes
            )

    ##
    ## labels
    ##
    add_labels(upper_pad, plot_description.region_name, plot_description.var_to_plot)

    ##
    ## save
    ##
    output_dir = args.outdir
    utils.mkdir_p(output_dir)
    if not output_dir.endswith("/") :
        output_dir += "/"
    if args.suffix != "" :
        args.suffix = "_" + args.suffix
    save_name = output_dir + "%s_%s%s.pdf" % (plot_description.region_name, plot_description.var_to_plot, args.suffix)
    print " >>> Saving: %s" % os.path.abspath(save_name)
    canvas.fig.savefig(save_name, bbox_inches = "tight", dpi = 200)
    

def main() :

    parser = argparse.ArgumentParser( description = "Make some bespoke plots" )
    parser.add_argument("-c", "--config", required = True,
        help = "Provide the plotting configuration file"
    )
    parser.add_argument("--plot", required = True,
        help = "Provide the plot description for the plot to make: <region>_<variable>"
    )
    parser.add_argument("--abs", action = "store_true", default = False,
        help = "Take the absolute value of the variable"
    )
    parser.add_argument("--suffix", default = "",
        help = "Provide a suffix to append to any outputs"
    )
    parser.add_argument("--cache-dir", default = "./sample_cache",
        help = "Directory to place/look for the cached-samples"
    )
    parser.add_argument("--logy", default = False, action = "store_true",
        help = "Make plot have log y-scale"
    )
    parser.add_argument("--outdir", default = "./",
        help = "Provide an output directory to dump plots"
    )
    parser.add_argument("--apply-sr-sf", default = False, action = "store_true",
        help = "Apply the norm factors for Z+HF and Top in d_hh tail regions"
    )
    parser.add_argument("--add-syst", default = False, action = "store_true",
        help = "Add systematic uncertainty bands"
    )
    args = parser.parse_args()

    if not utils.file_exists(args.config) :
        sys.exit()

    plot_desc = PlotDescription(args.plot, args.abs)
    if plot_desc.var_to_plot not in bounds_dict() :
        print "ERROR Requested variable (=%s) not configured in bounds dict" % plot_desc.var_to_plot
        sys.exit()
    if plot_desc.var_to_plot not in nice_names_dict() :
        print "ERROR Requested variable (=%s) not in configured names dict" % plot_desc.var_to_plot
        sys.exit()
    print plot_desc

    global loaded_samples, loaded_regions, selected_region, loaded_plots
    selected_region = plot_desc.region_name
    loaded_samples = []
    loaded_regions = []
    loaded_plots = []
    execfile(args.config, globals(), locals())

    if len(loaded_samples) == 0 :
        print "ERROR No loaded samples found in provided configuration (=%s)" % args.config
        sys.exit()

    if not region_utils.regions_unique(loaded_regions) :
        print "ERROR Loaded regions are not unique, here are the counts:"
        for rname, count in region_utils.region_counts(loaded_regions).iteritems() :
            print " > %s : %d" % (rname, count)
        sys.exit()

    backgrounds, signals, data = sample_utils.categorize_samples(loaded_samples)

    region_to_plot = None
    for r in loaded_regions :
        if r.name == selected_region :
            region_to_plot = r
            break

    # cache
    cacher = sample_cacher.SampleCacher(args.cache_dir)
    cacher.samples = loaded_samples
    cacher.region = region_to_plot
    required_variables = get_required_variables(region_to_plot = region_to_plot, plot_description = plot_desc)
    cacher.fields = required_variables
    print cacher
    cacher.cache()

    make_paper_plot(region_to_plot, backgrounds, signals, data, plot_desc, args)

if __name__ == "__main__" :
    main()
