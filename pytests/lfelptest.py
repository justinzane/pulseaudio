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
FREQS = [250]
POLES = [1, 2, 3, 4, 5, 6, 7, 8]

SRC_FILE = "./stereo_20-20000Hz.wav"
REC_FILE = "./output/test-%d-%d.pcm"
PLOT_FILE = "./output/test_plot.png"
CONF_SRC_FILE = "./lfe-lp.pa"

PLAY_CMD = "../src/pacat -d alsa_output.pci-0000_00_1b.0.analog-stereo.lfe_lp " + SRC_FILE
REC_CMD = "../src/pacat -d alsa_output.pci-0000_00_1b.0.analog-stereo.monitor -r > " + REC_FILE

DAEMON_CMD = "../src/pulseaudio -n -F %s"
KILL_CMD_1 = "../src/pulseaudio -k; /usr/bin/pulseaudio -k"
KILL_CMD_3 = "killall lt-pacat"

FREQ_RE = re.compile(r'lpfreq=[0-9.]+')
POLE_RE = re.compile(r'lppoles=[0-9]+')


def run_test():
    for f in FREQS:
        for p in POLES:
            # setup temp config file
            dummy_handle, tf_name = tempfile.mkstemp()
            tf = open(tf_name, 'w+b')
            lines = open(CONF_SRC_FILE, 'r').readlines()
            for line in lines:
                if FREQ_RE.search(line):
                    line = FREQ_RE.sub("lpfreq=%d.0" % (f), line)
                    line = POLE_RE.sub("lppoles=%d" % (p), line)
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
            rec_p = sp.Popen(args=REC_CMD % (f, p), shell=True,
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
            print("Finished run with lpfreq=%d and lppoles=%d" % (f, p))

    # kill existing daemons
    kill1_p = sp.Popen(args=KILL_CMD_1, shell=True,
                       bufsize= -1, stdout=sys.stdout, stderr=sys.stderr)
    kill1_p.wait()


def analyze_test():
    tests_left = np.zeros([len(FREQS), len(POLES), MAX_SAMPLES],
                          dtype=np.int16)
#    tests_right = np.zeros([len(FREQS), len(POLES), MAX_SAMPLES], dtype=np.int16)

    font1 = FontProperties(family=["Anka/Coder Condensed", "monospace"],
                           style="normal",
                           variant="normal",
                           weight="normal",
                           size="x-small")

    fig, axar = pp.subplots(nrows=1, # len(POLES),
                            ncols=len(FREQS),
                            sharex=True,
                            sharey=True,
                            squeeze=False)

    for f in range(len(FREQS)):
        for p in range(len(POLES)):
            tf = open(REC_FILE % (FREQS[f], POLES[p]), 'r+b')
            counter = 0
            started = False
            frame = tf.read(4)
            while not (frame == '' or counter >= MAX_SAMPLES):
                l, r = struct.unpack("<hh", frame)
                if (started):
                    tests_left[f][p][counter] = l
#                    tests_right[f][p][counter] = r
                    counter += 1
                elif not (l == 0 and r == 0):
                    started = True
                frame = tf.read(4)

            mean_l = np.sqrt(np.mean(np.square(tests_left[f][p])))
            mean_db_l = 20.0 * np.log10(mean_l / pow(2, 15))
            print ("f=%d p=%d mean=%f" % (FREQS[f], POLES[p], mean_db_l))

            pp.sca(axar[0][f].axes)
            pp.hold(True)
#            if f == 0:
#            spectra_r = spectrum.Periodogram(data=tests_right[f][p],
#                                             sampling=Fs,
#                                             window='hann',
#                                             NFFT=pow(2, 15),
#                                             scale_by_freq=False,
#                                             detrend='mean')
#            spectra_r.run()
#            pp.sca(axar[0][f].axes)
#            spectra_r.plot(norm=True,
#                           axes=axar[0][f].axes,
#                           label="unfiltered",
#                           color=[1.0, 0.0, 0.0, 0.25])

            spectra_l = spectrum.pdaniell(data=tests_left[f][p], P=32,
                                             sampling=Fs,
                                             window='hann',
                                             NFFT='nextpow2',
                                             scale_by_freq=False,
                                             detrend='mean')
            spectra_l.run()
            spec_max = spectra_l.frequencies()[np.argmax(spectra_l.psd)]
            spectra_l.plot(norm=True,
                           axes=axar[0][f].axes,
                           label="%d poles" % (POLES[p]),
                           color=[(1.0 - (float(p) / len(POLES))) * ((p + 1) % 2),
                                  (1.0 - (float(p) / len(POLES))) * (p % 2),
                                  (float(p) / len(POLES)),
                                  0.85])
            axar[0][f].axes.axvline(x=spec_max,
                                    label="max=%04.1fHz" % (spec_max),
                                    color=[(1.0 - (float(p) / len(POLES))) * ((p + 1) % 2),
                                  (1.0 - (float(p) / len(POLES))) * (p % 2),
                                  (float(p) / len(POLES)),
                                  0.85])
            print "finished run %d %d" % (FREQS[f], POLES[p])

        axar[0][f].axes.axvline(x=FREQS[f],
                                color=[1.0, 0.0, 0.0, 0.75])
        axar[0][f].axes.legend(loc="best",
                               prop=font1,
                               frameon=False,
                               title="Cutoff %d" % FREQS[f])
        axar[0][f].axes.xaxis.label.fontproperties = font1
        axar[0][f].axes.tick_params(axis='both',
                                    which='both',
                                    direction='in',
                                    length=3,
                                    pad=2,
                                    labelsize='x-small')
        axar[0][f].axes.grid(True,
                             which='both',
                             axis='both',
                             color=[0.5, 0.5, 0.5, 0.25])
        axar[0][f].axes.set_xscale("log")
        axar[0][f].axes.set_xlim(20, 20000)
        axar[0][f].axes.set_ylim(-80, 0)

    pp.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
    pp.show()

#    fig.set_dpi(600)
#    fig.set_size_inches(10, 8)
#    fig.savefig(filename=PLOT_FILE)

###############################################################################

if __name__ == "__main__":
#    run_test()
    analyze_test()
