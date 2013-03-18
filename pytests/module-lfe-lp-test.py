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

import subprocess as sp
import re, math, tempfile, time, sys
import numpy as np
import numpy.fft as fft
import struct
import matplotlib.pyplot as pp
pp.ioff()

Fs = 44100
MAX_SAMPLES = int(math.pow(2, 18))
test_freqs = [100]
test_poles = range(1, 6)

fig_file = "/home/justin/tmp/test-figure.png"
rec_file = "/home/justin/tmp/test-%d-%d.pcm"
conf_file = "/home/justin/src/pulseaudio/src/lfe-lp.pa"

play_cmd = "/home/justin/src/pulseaudio/src/pacat"
play_cmd += " -d alsa_output.pci-0000_00_1b.0.analog-stereo.lfe_lp"
play_cmd += " /home/justin/shared-music/Misc/1000-20-stereo.wav"

daemon_cmd = "/home/justin/src/pulseaudio/src/pulseaudio -n -F %s"
kill1_cmd = "/home/justin/src/pulseaudio/src/pulseaudio -k"
kill2_cmd = "/usr/bin/pulseaudio -k"

rec_cmd = "/home/justin/src/pulseaudio/src/pacat"
rec_cmd += " -d alsa_output.pci-0000_00_1b.0.analog-stereo.monitor -r > "
rec_cmd += rec_file

freq_re = re.compile(r'lpfreq=[0-9.]+')
pole_re = re.compile(r'lppoles=[0-9]+')

def fnp2rgb(f, p):
    r = int(255.0 * float(f) / float(len(test_freqs)))
    if (p % 2 == 0):
        g = int(255.0 * float(p) / float(len(test_poles)))
        b = 0
    else:
        b = int(255.0 * float(p) / float(len(test_poles)))
        g = 0
    return("#%x%x%x" % (r, g, b))

def run_test():
    for f in test_freqs:
        for p in test_poles:
            # setup temp config file
            tf_handle, tf_name = tempfile.mkstemp()
            tf = open(tf_name, 'w+b')
            lines = open(conf_file, 'r').readlines()
            for line in lines:
                if freq_re.search(line):
                    line = freq_re.sub("lpfreq=%d.0" % (f), line)
                    line = pole_re.sub("lppoles=%d" % (p), line)
                tf.writelines(line)
            tf.close()
            # kill existing daemons
            kill1_p = sp.Popen(args=kill1_cmd,
                               shell=True,
                               bufsize= -1,
                               stdout=sys.stdout,
                               stderr=sys.stderr)
            kill1_p.wait()
            kill2_p = sp.Popen(args=kill2_cmd,
                               shell=True,
                               bufsize= -1,
                               stdout=sys.stdout,
                               stderr=sys.stderr)
            kill2_p.wait()
            # run daemon
            daemon_p = sp.Popen(args=daemon_cmd % (tf_name,),
                                shell=True,
                                bufsize= -1,
                                stdout=sys.stdout,
                                stderr=sys.stderr)
            time.sleep(5)
            # record
            rec_p = sp.Popen(args=rec_cmd % (f, p),
                             shell=True,
                             bufsize= -1,
                             stdout=sys.stdout,
                             stderr=sys.stderr)
            # play
            play_p = sp.Popen(args=play_cmd,
                              shell=True,
                              bufsize= -1,
                              stdout=sys.stdout,
                              stderr=sys.stderr)
            play_p.wait()
            rec_p.kill()
            daemon_p.kill()
            print("Finished run with lpfreq=%d and lppoles=%d" % (f, p))

def analyze_test():
    ffts = np.empty([len(test_freqs), len(test_poles), MAX_SAMPLES / 2], dtype=np.float)
    for f in range(len(test_freqs)):
        for p in range(len(test_poles)):
            left = []
            tf = open(rec_file % (test_freqs[f], test_poles[p]), 'r+b')
            frame = tf.read(4)
            while not frame == '':
                l, r = struct.unpack("<hh", frame)
                left.append(np.int16(l))
                frame = tf.read(4)
            left = np.array(left, dtype=np.int16)[0:MAX_SAMPLES]

            start = time.clock()
            lft = fft.fft(left) / MAX_SAMPLES
            print("fft took %f" % ((time.clock() - start),))
            ffts[f][p] = np.real(lft[range(MAX_SAMPLES / 2)])

    pp.ioff()
    fig = pp.figure()
    frq = np.arange(MAX_SAMPLES, dtype=np.int16) / (float(MAX_SAMPLES) / float(Fs))
    frq = frq[range(MAX_SAMPLES / 2)]
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    for f in range(len(test_freqs)):
        for p in range(len(test_poles)):
            ax.plot(frq, ffts[f][p], alpha=0.5,
                    label="%d-%d" % (test_freqs[f], test_poles[p]))
    ax.legend()
    pp.show()

if __name__ == '__main__':
    run_test()
    analyze_test()
