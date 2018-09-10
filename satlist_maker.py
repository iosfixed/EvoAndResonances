from tqdm import tqdm
from numpy import sqrt, pi, dot, sin, cos, arctan2, cross, arcsin, arctanh
import pandas as pd
from libcelestial import Sat, readeph
from pathlib import Path
from utils import sat_counter, runcmd, tail
from multiprocessing import Pool

rad = pi / 180
EPH_DIR = Path('/home/dmitry/DATA/GPS_inclinated/EPH/')


def parse_file(infile):
    sat = Sat()
    N_sats = sat_counter(infile)

    with open(infile, 'r') as file:
        Houtput = [file.readline() for i in range(2 * N_sats + 1)]

    Toutput = tail(infile, 2*N_sats)

    sat_dicts = []
    for sat_N in range(N_sats):
        q0_str = Houtput[1 + 2 * sat_N]
        v0_str = Houtput[2 + 2 * sat_N]

        q0, v0 = readeph(q0_str, 'q'), readeph(v0_str, 'v')
        q0, v0 = q0[1:4], v0[0:3]
        qv = q0 + v0

        sat.set(qv, 'qv')
        init_a, init_i = round(sat.a), round(sat.i / rad)

        q1_str = Toutput[2 * sat_N]
        v1_str = Toutput[2 * sat_N + 1]

        q1, v1 = readeph(q1_str, 'q'), readeph(v1_str, 'v')
        if v1 == None or q1 == None:
            break
        qv1 = q1[1:4] + v1[0:3]

        sat.set(qv1, 'qv')
        e, mMEGNO = sat.e, v1[3]

        sat_dicts.append({'file': infile.name, 'sat_N': sat_N, 'a': init_a, 'i': init_i, 'e': e, 'mMEGNO': mMEGNO})

    return sat_dicts


files = sorted(EPH_DIR.glob('EPH_*.OUT'))
satlist = pd.DataFrame(columns=['file', 'sat_N', 'a', 'i', 'mMEGNO'])
sat_dict_list = []
with Pool(4) as p:
    with tqdm(total=len(files), desc='Parse files...') as progress_bar:
        for sat_dict in p.imap_unordered(parse_file, files):
            sat_dict_list += sat_dict
            progress_bar.update()

for sat_dict in tqdm(sat_dict_list, desc='Appending to dataframe'):
    satlist = satlist.append(sat_dict, ignore_index=True)

satlist.to_csv(EPH_DIR.parent / 'satlist.csv', index=False)
