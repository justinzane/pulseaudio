#!/usr/bin/python
# encoding: utf-8
'''
@author:    Justin Chudgar <justin@justinzane.com>
@updated:    2013-03-14 22:41:01
@license:
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import matplotlib
matplotlib.use("qt4agg", warn=True)
import matplotlib.pyplot as pp
from matplotlib.font_manager import FontProperties
pp.ioff()
import numpy as np
import re
import tempfile
import time
import sys
import spectrum
import struct
import subprocess as sp

Fs = 44100
MAX_SAMPLES = Fs * 30
CHANS = 6
FREQS = [249]# 39, 79, 119, 249]

COLORS = ["#7f0000", "#5f1f00", "#5f0f0f", "#7f001f", "#007f00", "#00007f"]

SRC_FILE = "./stereo_20-20000Hz.wav"
REC_FILE = "./output/test-%d.pcm"
PLOT_FILE = "./output/test_plot.png"
CONF_SRC_FILE = "./lfe-lp.pa"

# PLAY_CMD = "../src/pacat -d alsa_output.pci-0000_00_1b.0.analog-stereo.lfe_lp " + SRC_FILE
# REC_CMD = "../src/pacat -d alsa_output.pci-0000_00_1b.0.analog-stereo.monitor -r > " + REC_FILE

PLAY_CMD = "../src/pacat -d null.lfe_lp " + SRC_FILE
REC_CMD = "../src/pacat --rate=%d --format=s16le --channels=%d " % (Fs, CHANS)
REC_CMD += "--channel-map=\"front-left,front-right,rear-left,rear-right,front-center,lfe\" "
REC_CMD += "-d null.monitor -r > " + REC_FILE

DAEMON_CMD = "../src/pulseaudio -n -F %s"
KILL_CMD_1 = "../src/pulseaudio -k; /usr/bin/pulseaudio -k"
KILL_CMD_3 = "killall lt-pacat"

FREQ_RE = re.compile(r'lpfreq=[0-9.]+')

###############################################################################


def run_test(f):
    # setup temp config file
    dummy_handle, tf_name = tempfile.mkstemp()
    tf = open(tf_name, 'w+b')
    lines = open(CONF_SRC_FILE, 'r').readlines()
    for line in lines:
        if FREQ_RE.search(line):
            line = FREQ_RE.sub("lpfreq=%d.0" % (f), line)
        tf.writelines(line)
    tf.close()

    # kill existing daemons
    kill1_p = sp.Popen(args=KILL_CMD_1, shell=True,
                       bufsize= -1, stdout=sys.stdout, stderr=sys.stderr)
    kill1_p.wait()

    # run daemon
    daemon_p = sp.Popen(args=DAEMON_CMD % (tf_name,), shell=True,
                        bufsize= -1, stdout=sys.stdout, stderr=sys.stderr)
    time.sleep(2)

    # record
    rec_p = sp.Popen(args=REC_CMD % (f), shell=True,
                     bufsize= -1, stdout=None, stderr=None)
    time.sleep(0.001)
    # play
    play_p = sp.Popen(args=PLAY_CMD, shell=True,
                      bufsize= -1, stdout=None, stderr=None)
    play_p.wait()

    rec_p.kill()
    time.sleep(1)
    daemon_p.kill()
    time.sleep(1)
    # kill stray pacats
    kill3_p = sp.Popen(args=KILL_CMD_3, shell=True,
                       bufsize= -1, stdout=sys.stdout, stderr=sys.stderr)
    kill3_p.wait()
    time.sleep(1)
    print("Finished run with lpfreq=%d" % (f))

    # kill existing daemons
    kill1_p = sp.Popen(args=KILL_CMD_1, shell=True,
                       bufsize= -1, stdout=sys.stdout, stderr=sys.stderr)
    kill1_p.wait()


def analyze_test(f, ax):
    test_data = np.zeros([MAX_SAMPLES, CHANS], dtype=np.int16)

    # read data from test file
    tf = open(REC_FILE % f, 'r+b', 1536)
    counter = 0
    started = False
    frame = tf.read(CHANS * 2)
    while not (frame == '' or
               counter >= MAX_SAMPLES):
        samples = list(struct.unpack("<hhhhhh", frame))
        if (started):
            test_data[counter] = np.array(samples, dtype=np.int16)
            counter += 1
        elif not (max(samples) == 0 and min(samples) == 0):
            started = True
        frame = tf.read(CHANS * 2)
    print "Finished reading data."

    # calculate psd's
    spectra = []
    for c in range(CHANS):
#         spectra.append(spectrum.pdaniell(data=test_data.transpose()[c],
#                                        P=44,
#                                        sampling=Fs,
#                                        window='hann',
#                                        NFFT='nextpow2',
#                                        scale_by_freq=True,
#                                        detrend='none'))
#         spectra.append(spectrum.pmusic(data=test_data.transpose()[c],
#                                        sampling=Fs,
#                                        IP=15, NFFT='nextpow2'))
        nfft = pow(2, int(np.log2(test_data.transpose()[c].size)) + 1)
        spectra.append(spectrum.Periodogram(data=test_data.transpose()[c],
                                            sampling=Fs, window='kaiser',
                                            NFFT=nfft, scale_by_freq=True,
                                            detrend='none'))
        spectra[c].run()
        spectra[c].plot(norm=True,
                        axes=ax,
                        label="chan %d" % (c),
                        linestyle='-',
                        linewidth=1,
                        color=COLORS[c],
                        alpha=1.0)
        print "finished chan %d" % (c)

    font1 = FontProperties(family=["Anka/Coder Condensed", "monospace"],
                           style="normal",
                           variant="normal",
                           weight="normal",
                           size="small")

    ax.axvline(x=f, color=[0.5, 0.5, 0.5, 0.5])
    ax.xaxis.label.fontproperties = font1
    ax.yaxis.label.fontproperties = font1
    ax.tick_params(axis='both', which='both', direction='in', length=3,
                   pad=2, labelsize='x-small')
    ax.grid(True, which='both', axis='both', color=[0.5, 0.5, 0.5, 0.5])
    ax.set_xscale("log")
    ax.set_xlim(20, 20000)
    ax.set_ylim(-96, 0)
    leg = ax.legend(loc="lower center", ncol=CHANS / 2, prop=font1)
    leg.get_frame().set_alpha(0.25)

#    fig.set_dpi(600)
#    fig.set_size_inches(10, 8)
#    fig.savefig(filename=PLOT_FILE)

###############################################################################

if __name__ == "__main__":
    for f in FREQS:
        run_test(f)
    fig, ax = pp.subplots(nrows=len(FREQS), ncols=1,
                         sharex=True, sharey=True, squeeze=False)
    for i in range(len(FREQS)):
        analyze_test(FREQS[i], ax[i][0])

    pp.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
    pp.show()

