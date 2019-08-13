#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# ampayne_tools/plot_bokeh                                                    #
#                                                                             #
# Copyright (c) 2019- Mark Payne (ampayne@mit.edu)                            #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Contains functions for quickly creating plots using the Bokeh package https://bokeh.pydata.org/en/latest/
"""

from __future__ import absolute_import, print_function, division

import bokeh.plotting as bp
import bokeh.palettes as bpp
import bokeh.models as bm
from bokeh.io import export_png


def plot2d(xs, ys, xlabel=None, ylabel=None, log_x=False, log_y=False, legend=None, path=None, line_width=5,
           plot_width=800, plot_height=500):
    """
    Generate a 2D line plot of one or more x/y series on the same plot

    Args:
        xs (list): [np.ndarray, ...]
        ys (list): [np.ndarray, ...]
        xlabel (str): x-axis label
        ylabel (str): y-axis label
        log_x (bool): If `True` the x axis will be log-scale
        log_y (bool): If `True` the y axis will be log-scale
        legend (list): [label_1, label_2, ...]
        path (str): If specified the generate plot will be saved to disk at this path
        line_width (float): Line thickness for every series
        plot_width (int): number of pixels
        plot_height (int): number of pixels
    """
    color_palette = bpp.Paired[12]

    x_axis_type = "log" if log_x else "linear"
    y_axis_type = "log" if log_y else "linear"

    fig = bp.figure(x_axis_label=xlabel, y_axis_label=ylabel, x_axis_type=x_axis_type, y_axis_type=y_axis_type,
                    plot_width=plot_width, plot_height=plot_height)

    legend_items = []

    for i, label in enumerate(legend):
        legend_items.append([label, [fig.line(xs[i], ys[i], line_width=line_width, color=color_palette[i])]])

    legend = bm.Legend(items=legend_items)
    fig.add_layout(legend, 'right')

    if path:
        fig.toolbar.logo = None
        fig.toolbar_location = None
        export_png(fig, path)

    else:
        bp.show(fig)


if __name__ == '__main__':
    pass
