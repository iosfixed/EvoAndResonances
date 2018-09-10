# -*- coding: utf-8 -*-
import os
import re

import numpy as np
import time
from numpy import sqrt, pi, dot, sin, cos, arctan2, cross, arcsin, arctanh
import pandas as pd
# import matplotlib.pyplot as plt

from celestial_main import mu_solar, ephem
from celestial_resonances import resonances, critical, resonances_labels, total as total_resonances
from libcelestial import Sat, readeph
from linecache import getline
from utils import sat_counter, runcmd, sat_line

rad = pi / 180

completed = []


def evo(infile, selected_sat, plot=True):
    infile = 'EPH/%s' % infile
    N_sats = sat_counter(infile)
    wc_out = runcmd(['wc', '-l', '%s' % infile])
    total_lines = int(wc_out.split(' ')[0])
    total_records = total_lines // (N_sats * 2 + 1)

    sat = Sat()
    sun_sat = Sat(mu=mu_solar)
    moon_sat = Sat()

    evolution = np.ndarray(shape=(total_records, 7), dtype='float64')

    moon_res = np.ndarray(shape=(total_records, total_resonances + 1), dtype='float64')
    moon_crit = np.ndarray(shape=(total_records, total_resonances + 1), dtype='float64')
    sun_res = np.ndarray(shape=(total_records, total_resonances + 1), dtype='float64')
    sun_crit = np.ndarray(shape=(total_records, total_resonances + 1), dtype='float64')

    for rec in range(total_records):
        dateline = sat_line(N_sats, 0, rec) - 1
        q_line = sat_line(N_sats, selected_sat, rec)
        v_line = sat_line(N_sats, selected_sat, rec) + 1

        date_str = getline(infile, dateline)
        q_str = getline(infile, q_line)
        v_str = getline(infile, v_line)

        date_match = readeph(date_str[0], 'date')
        JD = date_match[0] + date_match[1] / 86400
        year = date_match[2] + date_match[3] / 12 + date_match[4] / 365

        if q_str == '' or v_str == '':
            # TODO: Emergency save
            pass

        q_match, v_match = readeph(q_str, 'q'), readeph(v_str, 'v')

        MEGNO, mMEGNO = q_match[4], v_match[3]
        q, v = q_match[1:4], v_match[0:3]

        moonqv = ephem(JD, 'moon')
        earthqv = ephem(JD, 'earth')
        sunqv = -1 * earthqv

        sat.set(q + v, 'qv')
        sun_sat.set(sunqv, 'qv')
        moon_sat.set(moonqv, 'qv')

        evolution[rec, :] = (JD, year, sat.a, sat.e, sat.i, MEGNO, mMEGNO)
        moon_res[rec, 1:] = resonances(sat, 'moon')
        moon_crit[rec, 1:] = critical(sat, moon_sat)
        sun_res[rec, 1:] = resonances(sat, 'sun')
        sun_crit[rec, 1:] = critical(sat, sun_sat)

    moon_res[:, 0] = moon_crit[:, 0] = sun_res[:, 0] = sun_crit[:, 0] = evolution[:, 1]  # Time column

    a, i, mMEGNO = round(evolution[0, 2]), round(evolution[0, 4] / rad), evolution[-1, 6]
    e = evolution[-1, 3]
    if e < 1:
        sat_dir = 'EVO/%d - %d (%.1f)/' % (a, i, mMEGNO)
    else:
        sat_dir = 'EVO/%d - %d (%.1f) H/' % (a, i, mMEGNO)

    os.makedirs(sat_dir, exist_ok=True)
    os.makedirs(sat_dir + 'MOON/', exist_ok=True)
    os.makedirs(sat_dir + 'SUN/', exist_ok=True)

    np.savetxt(sat_dir + 'evo.csv', evolution, delimiter=',', header='JD,year,a,e,i,MEGNO,mMEGNO')
    np.savetxt(sat_dir + 'moon_res.csv', moon_res, delimiter=',')
    np.savetxt(sat_dir + 'moon_crit.csv', moon_crit, delimiter=',')
    np.savetxt(sat_dir + 'sun_res.csv', sun_res, delimiter=',')
    np.savetxt(sat_dir + 'sun_crit.csv', sun_crit, delimiter=',')

    if plot:
        labels = resonances_labels('Moon')
        for expression in range(1, moon_res.shape[1]):
            f = plt.figure()
            plt.subplot(211)
            plt.scatter(moon_crit[:, 0], moon_crit[:, expression], s=1)
            plt.ylabel('Крит. аргумент, рад')
            plt.title(labels[expression - 1])
            plt.xticks([], [])

            plt.subplot(212)
            plt.plot(moon_res[:, 0], moon_res[:, expression])
            if np.sign(min(moon_res[:, expression]) * max(moon_res[:, expression])) == -1:
                plt.plot([moon_res[0, 0], moon_res[0, -1]], [0, 0], c='r')

            plt.ylabel('Рез. соотн., рад/с')
            plt.xlabel('Время, годы')
            plt.xlim((moon_res[0, 0], moon_res[-1, 0]))
            # f.subplots_adjust(hspace=0)
            plt.savefig(sat_dir + 'MOON/%d.png' % expression)
            plt.close()

        labels = resonances_labels('Sun')
        for i in range(1, sun_res.shape[1]):
            f = plt.figure()
            plt.subplot(211)
            plt.scatter(sun_crit[:, 0], sun_crit[:, i], s=1)
            plt.ylabel('Крит. аргумент, рад')
            plt.title(labels[i - 1])
            plt.xticks([], [])

            plt.subplot(212)
            plt.plot(sun_res[:, 0], sun_res[:, i])
            if np.sign(min(sun_res[:, expression]) * max(sun_res[:, expression])) == -1:
                plt.plot([sun_res[0, 0], sun_res[0, -1]], [0, 0], c='r')

            plt.ylabel('Рез. соотн., рад/с')
            plt.xlabel('Время, годы')
            plt.xlim((moon_res[0, 0], moon_res[-1, 0]))

            # f.subplots_adjust(hspace=0)
            plt.savefig(sat_dir + 'SUN/%d.png' % i)
            plt.close()


satlist = pd.read_csv('satlist.csv')

for a in np.arange(220000, 320001, 10000):
    for i in [5, 10, 15, 30, 45, 60, 75, 90]:
        tic = time.time()
        print('a = %d km, i = %d° =========================' % (a, i))
        print('Start: ', time.strftime('%Y-%m-%d %H:%M:%S'))

        target = satlist[(satlist['a'] == a) & (satlist['i'] == i)]
        evo(target['file'].values[0], target['sat_N'].values[0])

        print('End: ', time.strftime('%Y-%m-%d %H:%M:%S'))
        print('Time elapsed: %d seconds ===================\n' % (time.time() - tic))
