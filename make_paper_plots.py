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
    v["NN_d_hh"] = { "srIncNoDhh" : [1, -11, 11] }
    return v

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
    args = parser.parse_args()

    if not utils.file_exists(args.config) :
        sys.exit()

    plot_desc = PlotDescription(args.plot)
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

if __name__ == "__main__" :
    main()
