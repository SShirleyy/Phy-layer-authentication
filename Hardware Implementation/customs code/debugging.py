#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 <+YOU OR YOUR COMPANY+>.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is  USE_PYSIDE,distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

import numpy as np
from numpy import matlib as mat
from gnuradio import gr
import pmt
from pyqtgraph.Qt import QtGui, QtCore, USE_PYQT5
import pyqtgraph as pg
from pyqtgraph.ptime import time
import pyqtgraph.opengl as gl

#import sugar

class debugging(gr.basic_block):
    """
    docstring for block debugging
    """

    symbols_in = pmt.intern('symbols_in')
    meta_in = pmt.intern('meta_in')
    channel_in = pmt.intern('channel')
    PMT_OUT = pmt.intern('out')

    snr_file = open('/tmp/snr','w')
    offset_file = open('/tmp/offset','w')
    chanfile = open('/tmp/channel','w')

    def __init__(self, write_bool, feature):
        gr.basic_block.__init__(self,
            name="debugging",
            in_sig=[],
            out_sig=[])
        self.feature = feature

        self.message_port_register_in(self.symbols_in)
        self.set_msg_handler(self.symbols_in, self.handle_symbols)

        self.message_port_register_in(self.channel_in)
        self.set_msg_handler(self.channel_in, self.handle_channel)

        self.message_port_register_in(self.meta_in)
        self.set_msg_handler(self.meta_in, self.handle_meta)

        self.message_port_register_out(self.PMT_OUT)


        #self.mac_frame_count = 0

        self.channel = np.array([])
        self.channel_fifo = []

        self.snr_upper = 12.8
        self.snr_lower = 9
        self.ofs_upper = -15000
        self.ofs_lower = -30000
        self.chan_limit = 0.6

        with open('/media/karl/HYPERX/m_alice.txt') as f:
           self.m_alice = f.read().splitlines()
           self.m_alice = map(complex,self.m_alice)

        #self.symbols = np.array([1-1j, 1+1j, -1-1j, -1+1j])

        self.snr_curve_list = []
        self.snr_curve_datas_list = []
        self.snr_mean_curve_list = []
        self.snr_mean_datas_list = []

        self.ofs_curve_list = []
        self.ofs_curve_datas_list = []

        self.chan_data = []

        self.mac_list = []

        self.app = QtGui.QApplication([])
        self.app.aboutToQuit.connect(self.stop)
        self.win = pg.GraphicsWindow(title="eq3440 plotting")
        self.win.resize(1000,1000)
        #self.win.setWindowTitle('eq3440 plotting')

        self.snr = self.win.addPlot(title="SNR plot")
        self.snr_legend = self.snr.addLegend()
        self.snr_curve = self.snr.plot()
        self.snr.addItem(pg.InfiniteLine(self.snr_upper,0))
        self.snr.addItem(pg.InfiniteLine(self.snr_lower,0))

        self.vb = self.win.addViewBox()
        self.vb.setMaximumWidth(300)
        self.snr_legend.setParentItem(self.vb)
        self.snr_legend.anchor((0,0),(0,0))

        self.win.nextRow()

        self.ofs = self.win.addPlot(title="offset plot")
        self.ofs_curve = self.ofs.plot()
        self.ofs.addItem(pg.InfiniteLine(self.ofs_upper,0))
        self.ofs.addItem(pg.InfiniteLine(self.ofs_lower,0))

        self.win.nextRow()

        #self.qi = self.win.addPlot(title="constellation plot")
        #self.qi.setRange(xRange=[-2, 2], yRange=[-2, 2])
        #self.qi_curve = self.qi.plot(x=self.symbols.real,y=self.symbols.imag,pen=None,symbol='o',symbolSize=6)

        self.chan = self.win.addPlot(title="channel plot")
        self.chan_curve = self.chan.plot()
        self.chan.addItem(pg.InfiniteLine(self.chan_limit,0))

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_plot)
        self.timer.start(100)



    def __exit__(self, exc_type, exc_value, traceback):
        self.snr_file.close()
        self.offset_file.close()

    def stop(self):
        self.snr_file.close()
        self.offset_file.close()
        self.chanfile.close()

    def handle_symbols(self, sym_pdu):
        self.symbols = pmt.to_python(pmt.cdr(sym_pdu))

    def update_plot(self):
        for index in xrange(len(self.snr_curve_list)):
            self.snr_curve_list[index].setData(self.snr_curve_datas_list[index])
            self.ofs_curve_list[index].setData(self.ofs_curve_datas_list[index])
            self.snr_mean_curve_list[index].setData(self.snr_mean_datas_list[index])

        self.chan_curve.setData(self.chan_data)


    def handle_meta(self,meta_pdu):
        meta_car = pmt.to_python(pmt.car(meta_pdu))
        mac_frame = pmt.to_python(pmt.cdr(meta_pdu))

        dest_mac = ':'.join(map(hex,mac_frame[4:10].tolist()))
        src_mac = ':'.join(map(hex,mac_frame[10:16].tolist()))
        single_snr = meta_car['snr']
        single_ofs = meta_car['freqofs']

        #self.mac_frame_count = self.mac_frame_count + 1

        if src_mac not in self.mac_list and not src_mac=='' :
            self.mac_list.append(src_mac)

            self.snr_curve_list.append(pg.PlotCurveItem(name=src_mac, pen=(1,len(self.mac_list))))
            self.snr_curve_datas_list.append([single_snr])
            self.snr_curve_list[-1].setData(self.snr_curve_datas_list[-1])
            self.snr_mean_curve_list.append(pg.PlotCurveItem(pen=(1,len(self.mac_list))))
            self.snr_mean_datas_list.append([single_snr])

            self.ofs_curve_list.append(pg.PlotCurveItem(name=src_mac, pen=(1,len(self.mac_list))))
            self.ofs_curve_datas_list.append([single_ofs])
            self.ofs_curve_list[-1].setData(self.snr_curve_datas_list[-1])

            self.snr.addItem(self.snr_curve_list[-1])
            self.snr.addItem(self.snr_mean_curve_list[-1])
            self.ofs.addItem(self.ofs_curve_list[-1])

        elif src_mac in self.mac_list:
            index = self.mac_list.index(src_mac)
            self.snr_curve_datas_list[index] = np.append(self.snr_curve_datas_list[index],single_snr)
            self.ofs_curve_datas_list[index] = np.append(self.ofs_curve_datas_list[index],single_ofs)
            #self.snr_curve_list[index].setData(self.snr_curve_datas_list[index])
            #self.ofs_curve_list[index].setData(self.ofs_curve_datas_list[index])
            #self.snr_mean_curve_list[index].setData(self.snr_mean_datas_list[index])

            #snr_single_mean = ((self.mac_frame_count-1) * self.snr_mean_datas_list[index][-1] / self.mac_frame_count) + single_snr/self.mac_frame_count
            snr_mean = np.mean(self.snr_curve_datas_list[index])
            self.snr_mean_datas_list[index] = np.append(self.snr_mean_datas_list[index],snr_mean)

            if len(self.snr_curve_datas_list[index]) > 200:
                self.snr_curve_datas_list[index] = np.delete(self.snr_curve_datas_list[index],0)
                self.snr_mean_datas_list[index] = np.delete(self.snr_mean_datas_list[index],0)
                self.ofs_curve_datas_list[index] = np.delete(self.ofs_curve_datas_list[index],0)

        if (self.feature != 'none'):
            self.decide(single_snr,single_ofs)
        self.snr_file.write(str(meta_car['snr']) + '\n')
        self.offset_file.write(str(meta_car['freqofs']) + '\n')


    def handle_channel(self, chan_pdu):
        self.channel = pmt.to_python(pmt.cdr(chan_pdu))
        self.channel = self.channel[6:58]
        for sym in self.channel:
            self.chanfile.write('('+str(sym.real)+','+str(sym.imag)+') ')
        self.chanfile.write('\n')
        m_alice = self.m_alice
        #phi = np.angle(self.channel.conj().T*m_alice)
        #delta = self.channel - np.dot(m_alice, np.exp(1j*phi))
        delta = self.channel - m_alice

        L = np.abs(np.sum(delta))
        self.chan_data.append(L)
        #self.chan_curve.setData(self.chan_data)


        if len(self.chan_data) > 200:
            self.chan_data.pop(0)

        if L < self.chan_limit:
            self.channel_fifo.append(True)
        else:
            self.channel_fifo.append(False)


        if len(self.channel_fifo) > 1:
            self.channel_fifo.pop(0)

        #     print '\n'+str(L)+'\n'


    def decide(self,single_snr,single_ofs):
        snr_alice = False
        ofs_alice = False
        chan_alice = False

        snr_weight = 0
        ofs_weight = 0
        chan_weight = 0

        if (self.feature == 'all'):
            snr_weight = 1
            ofs_weight = 1
            chan_weight = 1
        elif (self.feature == 'snr'):
            snr_weight = 1
        elif (self.feature == 'offset'):
            ofs_weight = 1
        elif (self.feature == 'channel'):
            chan_weight = 1


        if single_snr < self.snr_upper and single_snr > self.snr_lower:
            snr_alice = True
        if single_ofs < self.ofs_upper and single_ofs > self.ofs_lower:
            ofs_alice = True
        if sum(self.channel_fifo)/len(self.channel_fifo) > 0.5:
            chan_alice = True

        if float(snr_alice*snr_weight + ofs_alice*ofs_weight + chan_alice*chan_weight)/float(snr_weight+ofs_weight+chan_weight) > 0.5:
            out_dict = pmt.from_double(1)
            #print '\n'+'Alice: '+'snr:'+str(snr_alice)+' ofs:'+str(ofs_alice)+' chan:'+str(chan_alice)+'\n'
            #print 'weights: '+str(snr_weight)+str(ofs_weight)+str(chan_weight)
        else:
            out_dict = pmt.from_double(0)
            #print '\n'+'eve: '+'snr:'+str(snr_alice)+' ofs:'+str(ofs_alice)+' chan:'+str(chan_alice)+'\n'
            #print 'weights: '+str(snr_weight)+str(ofs_weight)+str(chan_weight)

        self.message_port_pub(self.PMT_OUT, out_dict)


    def set_decision(self,feature):
        self.feature = feature
