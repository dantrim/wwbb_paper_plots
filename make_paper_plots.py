#!/usr/bin/env python

import sys, os
import argparse

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
import numpy as np
import json

class PlotDescription :

    def __init__(self, descriptor = "") :
        self.descriptor = descriptor

        # variable stuff
        self.var_to_plot = ""
        self.is_abs = False

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
    return variables

def get_trigger_idx(arr) :

    idx_15 = (arr['year'] == 2015) & (arr['trig_tight_2015'] == 1)
    idx_16 = (arr['year'] == 2016) & (arr['trig_tight_2016'] == 1)
    idx_17 = (arr['year'] == 2017) & (arr['trig_tight_2017rand'] == 1)
    idx_18 = (arr['year'] == 2018) & (arr['trig_tight_2018'] == 1)
    idx = (idx_15 | idx_16 | idx_17 | idx_18)
    return idx

def bounds_dict() :

    v = {}
    sr_sf_cut_bins = np.arange(5.45, 12+0.45, 0.45)
    v["NN_d_hh"] = { "srIncNoDhh" : [1, -11, 11],
                        "srPreSel" : [1, -11, 11],
                    "srIncNoDhhClose" : [1, 0, 12],
                    "srSFNoDhhClose" : [0.45, 5.45, 12],
                    "srDFNoDhhClose" : [0.5, 4, 9],
                    #"srSFNoDhhCloseCut" : [5.45, 7, 8, 9, 10, 11, 12], #, 6.5, 7, 7.5, 8, 9, 10, 11, 12], #, 5.45, 12], #sr_sf_cut_bins,
                    #"srSFNoDhhCloseCut" : [5.45, 5.90, 6.35, 7, 8, 10, 12], #, 8.15, 8.60, 9.05],
                    "srSFNoDhhCloseCut" : [1, 5.5, 11.5],
                    "srDFNoDhhCloseCut" : [1, 5.5, 9.5],
                    "crTopNoDhh" : [1, -12, 12],
                    "crTop" : [1, 4.5, 11.5],
                    "crZNoDhh" : [1, -12, 12],
                    "crZ" : [1, 0, 8],
    }
    #v["NN_d_hh"] = { "srIncNoDhh" : [-11, -5, -4, -3, 11] }
    return v

def region_nice_names_dict() :

    n = {}
    n["srIncNoDhh"] = "SR, $\\ell \\ell$-inc., no $d_{hh}$"
    n["srPreSel"] = "Pre-selection"
    n["srIncNoDhhClose"] = "SR, $\\ell \\ell$-inc., no $d_{hh}$"
    n["srSFNoDhhClose"] = "SR-SF, no $d_{hh}$"
    n["srDFNoDhhClose"] = "SR-DF, no $d_{hh}$"
    n["srSFNoDhhCloseCut"] = "SR-SF"
    n["srDFNoDhhCloseCut"] = "SR-DF"
    n["crTopNoDhh"] = "CR-Top, no $d_{hh}$"
    n["crTop"] = "CR-Top"
    n["crZNoDhh"] = "CR-Z+HF, no $d_{hh}$"
    n["crZ"] = "CR-Z+HF"
    return n

def nice_names_dict() :

    n = {}
    n["NN_d_hh"] = "$d_{hh}$"
    return n

def add_labels(pad, region_name = "", var_name = "") :


    x_atlas = 0.04
    y_atlas = 0.97
    x_type_offset = 0.24
    y_type = 0.97

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

    if (region_name == "crZ" or region_name == "crTop") and var_name == "NN_d_hh" :
        x_atlas = 0.33
        y_atlas = 0.96
        y_type = 0.96

        x_lumi = 0.33
        y_lumi = 0.88

        x_region = 0.33
        y_region = 0.8

    # ATLAS label
    size = 24
    text = 'ATLAS'
    opts = dict(transform = pad.transAxes)
    opts.update( dict(va = 'top', ha = 'left') )
    pad.text(x_atlas, y_atlas, text, size = size, style = 'italic', weight = 'bold', **opts)

    what_kind = 'Internal'
    pad.text(x_atlas + x_type_offset, y_type, what_kind, size = size, **opts)

    lumi = '140'#.5'
    pad.text(x_lumi, y_lumi, '$\\sqrt{s} = 13$ TeV, %s fb$^{-1}$' % lumi, size = 0.75 * size, **opts)

    # region
    if region_name == "srIncNoDhh" and var_name == "NN_d_hh" :
        region_text = "Selection:"
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)
        region_text = "\t%s" % region_nice_names_dict()[region_name]
        pad.text(x_region, 0.92 * y_region, region_text, size = 0.75 * size, **opts)
    else :
        region_text = "Selection: %s" % region_nice_names_dict()[region_name]
        pad.text(x_region, y_region, region_text, size = 0.75 * size, **opts)

def make_legend(ordered_labels, n_cols, var_name, region_name, pad) :

    handles, labels = pad.get_legend_handles_labels()
    new_handles = []
    new_labels = []
    for l in ordered_labels :
        for il, label in enumerate(labels) :
            if label == l :
                new_labels.append(l.replace("SIG",""))
                new_handles.append(handles[il])

    leg_x, leg_y = 0.45, 0.75
    legend_fontsize = 15

    if n_cols == 1 :
        leg_x = 0.58
        leg_y = 0.58

    # sr
    if n_cols == 2 and (region_name == "srSFNoDhhCloseCut" or region_name == "srDFNoDhhCloseCut") and var_name == "NN_d_hh" :
        leg_x = 0.27
        leg_y = 0.5
        legend_fontsize = 14

    # crZ
    if n_cols == 1 and (region_name == "crZ") and var_name == "NN_d_hh" :
        leg_x = 0.48
        leg_y = 0.3
        legend_fontsize = 16
    if n_cols == 2 and (region_name == "crZ" or region_name == "crTop") and var_name == "NN_d_hh" :
        leg_x = 0.31
        leg_y = 0.5
        legend_fontsize = 14


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

    palette = 3
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
                    "other" : "#44BBA4" }
                    #"other" : "#F6F7EB" }
                
    return colors[palette]

def get_legend_order(var_name, region_name) :

    standard_order = {}
    standard_order[1] = ["Data", "Total SM", "Top", "$Z/\\gamma*$+jets HF", "$t\\bar{t} + V$", "Higgs", "Other"]
    standard_order[2] = ["Data", "Top", "$Z/\\gamma*$+jets HF", "$t\\bar{t} + V$", "Total SM", "Higgs", "Other"]

    n_cols = {}
    n_cols["srIncNoDhh"] = 1
    n_cols["srIncNoDhhClose"] = 1
    n_cols["srSFNoDhhClose"] = 1
    n_cols["srDFNoDhhClose"] = 1
    n_cols["srSFNoDhhCloseCut"] = 2
    n_cols["srDFNoDhhCloseCut"] = 2
    n_cols["crTopNoDhh"] = 1
    n_cols["crTop"] = 2
    n_cols["crZNoDhh"] = 1
    n_cols["crZ"] = 2

    order_dict = {}
    order_dict["NN_d_hh"] = {
        "srIncNoDhh" : standard_order[n_cols[region_name]],
        "srIncNoDhhClose" : standard_order[n_cols[region_name]],
        "srSFNoDhhClose" : standard_order[n_cols[region_name]],
        "srDFNoDhhClose" : standard_order[n_cols[region_name]],
        "srSFNoDhhCloseCut" : standard_order[n_cols[region_name]],
        "srDFNoDhhCloseCut" : standard_order[n_cols[region_name]],
        "crTopNoDhh" : standard_order[n_cols[region_name]],
        "crTop" : standard_order[n_cols[region_name]],
        "crZNoDhh" : standard_order[n_cols[region_name]],
        "crZ" : standard_order[n_cols[region_name]],
    }

    return order_dict[var_name][region_name], n_cols[region_name]

def make_paper_plot(region, backgrounds, signals, data, plot_description, args) :

    histograms_bkg = {}
    labels_bkg = {}
    colors_bkg = {}
    top_bkg = ["ttbar", "wt"]
    other_bkg = ["zjetslf", "diboson"]
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
                elif plot_description.region_name == "srIncNoDhh" :
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

            weights = lumis * weights
            if bkg.name.lower() not in other_bkg and bkg.name.lower() not in top_bkg :
                h.fill(plot_data, weights)
            elif bkg.name.lower() in other_bkg :
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
    y_label = "Events / Bin"
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
    print "maxy = %.2f" % maxy

    if plot_description.region_name == "srSFNoDhhCloseCut" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 12
    elif plot_description.region_name == "srDFNoDhhCloseCut" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 10
    elif plot_description.region_name == "crZ" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 400
    elif plot_description.region_name == "crTop" and plot_description.var_to_plot == "NN_d_hh" and not args.logy :
        maxy = 70

    upper_pad.set_ylim(miny, maxy)

    # stat error
    sm_x_error = np.zeros(len(histo_total.y_error()))
    sm_y_error = histo_total.y_error()
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
                    lw = 1,
                    edgecolor = "none",
                    alpha = 1.0
    )

    # print yields
    sm_total_yield = histo_total.integral()
    print 30 * "*"
    print histo_total.count_str(name = "Total SM")

    # draw error band
    upper_pad.add_collection(stat_error_band)

    # draw total SM
    upper_pad.plot(sm_line[0], sm_line[1], ls = "-", color = "k", label = "Total SM", lw = 2)

    # draw data
    if not is_variable_width :
        data_x = np.array(h_data.bin_centers())
    else :
        x_left_edges = 1.0 * np.array(x_bounds[:-1])
        x_left_edges += 0.5 * hatch_bin_width
        data_x = x_left_edges
    data_y = h_data.histogram
    data_y[data_y == 0.] = -5
    upper_pad.plot(data_x, data_y, "ko", label = "Data", zorder = 1e6)

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
    if len(signals) > 0 and "sr" in plot_description.region_name.lower() :
        for signal in signals :
            h = histogram1d("signal_histo_%s" % signal.name, binning = x_bounds)
            chain = signal.chain()
            for isc, sc in enumerate(chain) :
                trigger_idx = get_trigger_idx(sc)
                sc = sc[trigger_idx]
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

            upper_pad.hist( h.data, weights = h.weights,
                        bins = bin_vals,
                        color = signal.color,
                        label = "$hh$",
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
    #lower_pad.set_yticks([0.5, 0.75, 1.0, 1.25, 1.5])
    #lower_pad.set_ylim([0.5, 1.5])

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
    lower_pad.plot(ratio_x, ratio_y, "ko", zorder = 1000)
    yerr = [ratio_data_err_low, ratio_data_err_high]
    lower_pad.errorbar(ratio_x, ratio_y, yerr = yerr, fmt = "none", color = "k")

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
            syst_files = { "crTop" : "%s/process_uncerts_NN_d_hh_crTop_paper.json" % syst_config_dir,
                            "crZ" : "%s/process_uncerts_NN_d_hh_crZ_paper.json" % syst_config_dir,
                            "srIncNoDhh" : "%s/process_uncerts_NN_d_hh_srIncNoDhh_paper.json" % syst_config_dir,
                            "srSFNoDhhCloseCut" : "%s/process_uncerts_NN_d_hh_srSFNoDhhCloseCut_paper.json" % syst_config_dir,
                            "srDFNoDhhCloseCut" : "%s/process_uncerts_NN_d_hh_srDFNoDhhCloseCut_paper.json" % syst_config_dir
            }

        syst_proc_names = { "TTbar" : "ttbar", "Wt" : "wt", "ZjetsHF" : "z" }
        process_fractions = {}
        for process in syst_proc_names :
            process_fractions[process] = np.array(histograms_bkg[process].histogram / histo_total.histogram)
            neg_frac = process_fractions[process] < 0
            process_fractions[process][neg_frac] = 0
        if plot_description.region_name in syst_files :
            json_file = syst_files[plot_description.region_name]
            with open(json_file, "r") as input_file :
                syst_data = json.load(input_file)
            process_errors = {}
            if syst_data["variable"] == plot_description.var_to_plot :
                processes = syst_data["processes"]
                for process in syst_proc_names :
                    for process_err_group in processes :
                        if process_err_group["name"] == syst_proc_names[process] :
                            err = process_err_group["errors"]

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

            sm_ratio_err = np.sqrt( np.power(sm_ratio_err,2) + np.power(total_err,2) )

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
            signal_label = "SM $hh \\rightarrow bb\\ell\\nu\\ell\\nu$"
            y_text = leg_y - 0.04
            #x_text = leg_x + 0.15
            x_text = leg_x + 0.12
            upper_pad.text(x_text, y_text, signal_label,
                        transform = upper_pad.transAxes,
                        fontsize = 14,
                        ha = "left"
            )

            # line
            y_line = 1.02 * y_text
            #x_line_0 = leg_x + 0.007 
            #x_line_1 = leg_x + 0.08
            x_line_0 = leg_x + 0.0155
            x_line_1 = leg_x + 0.09

            if n_cols == 2 :
                y_line = 1.02 * y_text
                x_line_0 = leg_x + 0.016
                x_line_1 = leg_x + 0.083
            upper_pad.plot([x_line_0, x_line_1], [y_line, y_line],
                    "-",
                    lw = 2,
                    color = signals[0].color,
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

    plot_desc = PlotDescription(args.plot)
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
