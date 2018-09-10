from subprocess import Popen, PIPE
from libcelestial import readeph
import os

def tail(infile, lines):
    with open(infile, 'rb') as f:
        new_line_counter = 0
        f.seek(-2, os.SEEK_END)
        while new_line_counter < lines:
            f.seek(-2, os.SEEK_CUR)
            if f.read(1) == b'\n':
                new_line_counter += 1

        return [f.readline().decode() for i in range(lines)]

# Определяет количество спутников в файле. Возвращает целое число.
def sat_counter(filename='EPH_000.OUT'):
    infile = open(filename, 'r')
    infile.readline()  # Вхолостую считываем строку с датой

    N = None
    while True:
        qstr = infile.readline()
        infile.readline()  # Вхолостую считываем строку со скоростями

        match = readeph(qstr, 'q')
        if match:
            N = match[0]
        else:
            break

    infile.close()
    return int(N)


def sat_line(N_sats, n_sat, rec):
    return (2 * N_sats + 1) * rec + 2 * (n_sat + 1)


def stdout_to_str(bstring):
    bstring = bstring.split(b'\r\n')
    stringlist = [bs.decode() for bs in bstring]
    if len(stringlist) == 1:
        return stringlist[0]
    return stringlist[:-1]  # Обрезаем пустую строку в коце.


def runcmd(cmd):
    proc = Popen(cmd, stdout=PIPE)
    b_out = proc.stdout.read()
    return stdout_to_str(b_out)
