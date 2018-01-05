# -*- coding: utf-8 -*-
import os
import re
from subprocess import Popen, PIPE
import numpy as np
import time
from numpy import sqrt, pi, dot, sin, cos, arctan2, cross, arcsin, arctanh
import pandas as pd
# import matplotlib.pyplot as plt

from celestial_main import mu_solar, ephem
from celestial_resonances import resonances, critical, resonances_labels, total as total_resonances

rad = pi / 180

completed = []


class Sat:
	t0 = 2451545.0

	def __init__(self, init=None, type='kep', t0=2451545.0, mu='geo'):
		self.t0 = t0
		if (mu == 'geo'):
			self.mu = 398600.4415  # км/с
		else:
			self.mu = 132712440018  # км/с

		if (type == 'kep'):
			if (init == None):
				init = [42164, 0, 0, 0, 0, 0]  # Нехай по умолчанию будет геостационарный спутник
			self.a, self.e, self.i, self.W, self.w, self.M0 = init
			self.cartesian()

		if (type == 'qv'):
			self.x, self.y, self.z, self.vx, self.vy, self.vz = init
			self.findOrbit()

	# Процедуры вычисления:
	def findOrbit(self):
		q = [self.x, self.y, self.z]
		v = [self.vx, self.vy, self.vz]

		r = sqrt(dot(q, q))
		V = sqrt(dot(v, v))
		s = dot(q, v)
		h = cross(q, v)

		k = sqrt(self.mu)

		# Большая полуось
		a = 1 / abs(2. / r - V ** 2 / k ** 2)

		# Эксцентрисистет:
		e = sqrt((s / k) ** 2 / a + (1 - r / a) ** 2)

		# Гиперболический случай:
		hyperbolic = True if e > 1 else False  # Флаг гиперболического случая.

		if hyperbolic:
			dx = s / (k * sqrt(a))
			dy = (r / a + 1)
			H0 = arctanh(dx / dy)

		else:
			# Средняя аномалия:
			dy = s / (e * k * sqrt(a))
			dx = (a - r) / (a * e)
			E0 = arctan2(dy, dx)
			M0 = E0 - e * sin(E0)

		# Долгота восходящего узла:
		W = arctan2(h[0], -h[1])

		# Наклонение
		i = arctan2(sqrt(h[0] ** 2 + h[1] ** 2), h[2])

		# Аргумент перицентра:


		p = a * (e ** 2 - 1) if hyperbolic else a * (1 - e ** 2)

		dy = sqrt(p) * s
		dx = k * (p - r)
		vv = arctan2(dy, dx)

		if (sin(i) != 0):
			dy = self.z / sin(i)
			dx = self.x * cos(W) + self.y * sin(W)
			uu = arctan2(dy, dx)
		else:
			uu = 0

		w = uu - vv

		while (w < 0):
			w += 2 * pi
		if hyperbolic:
			self.a, self.e, self.i, self.W, self.w, self.H0 = a, e, i, W, w, H0
			return [a, e, i, W, w, H0]
		else:
			self.a, self.e, self.i, self.W, self.w, self.M0 = a, e, i, W, w, M0
			return [a, e, i, W, w, M0]

	def cartesian(self, t=t0, dt=None):
		a, e, i, W, w, M0 = self.get('kep')
		mu = self.mu
		t0 = self.t0
		# Поворотные матрицы:
		A = [[cos(W), -sin(W), 0],
			 [sin(W), cos(W), 0],
			 [0, 0, 1]]

		A = np.matrix(A)

		B = [[1, 0, 0], [0, cos(i), -sin(i)], [0, sin(i), cos(i)]]
		B = np.matrix(B)

		C = [[cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]]
		C = np.matrix(C)

		R = A * B * C  # Конечная поворотная матрица

		n = sqrt(mu * a ** (-3))

		if (dt):
			M = M0 + n * dt
		else:
			M = M0 + n * (t - t0) * 86400

		E = M
		# Численное решение уравнения Кеплера
		while (abs(E - e * sin(E) - M) > 1e-12):
			E = E - (E - e * sin(E) - M) / (1 - e * cos(E))

		q = np.matrix([a * (cos(E) - e), a * sqrt(1 - e ** 2) * sin(E), 0])
		dq = np.matrix([-a * n * sin(E) / (1 - e * cos(E)), a * sqrt(1 - e ** 2) * cos(E) * n / (1 - e * cos(E)), 0])

		q = R * q.T  # T - преобразование строки к столбцу, суть транспозиция
		dq = R * dq.T

		self.x, self.y, self.z = q.A1  # A1 - преобразование к человеческому типу
		self.vx, self.vy, self.vz = dq.A1
		return [self.x, self.y, self.z, self.vx, self.vy, self.vz]

	# Процедуры установки новых параметров
	def set(self, settings, type='qv'):
		if (type == 'qv'):
			self.x, self.y, self.z, self.vx, self.vy, self.vz = settings
			self.findOrbit()
		if (type == 'kep'):
			self.a, self.e, self.i, self.W, self.w, self.M0 = settings
			self.cartesian()

	# Процедуры возвращающие параметры:
	def get(self, type='qv'):
		if (type == 'qv'):
			return [self.x, self.y, self.z, self.vx, self.vy, self.vz]

		if (type == 'kep'):
			return [self.a, self.e, self.i, self.W, self.w, self.M0]


# Чиитаем выходной файл численной модели движения ИСЗ
def readeph(str, type, mode='MEGNO'):
	if mode == 'MEGNO':
		# ------------#-----x--------------------y--------------------z--------------------MEGNO-----------------
		qrx = r'\s+(\d+)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s*'

		# ----------vx-------------------vy-------------------vz------------------ mMEGNO-------------
		vrx = r'\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s*'

		# --------------JD0---------dt---------------yr------mnth-----dy------HR------MIN------SEC--------
		daterx = r'\s?(\d+\.\d?)\s+(\d+\.\d+)\s+\(\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)H\s+(\d+)M\s+(\d+\.\d+)S\)\s?'

	else:
		# ------------#-----x--------------------y--------------------z--------------------MEGNO-----------------
		qrx = r'\s*(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s*'
		# ----------vx-------------------vy-------------------vz------------------ mMEGNO-------------
		vrx = r'\s*(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s+(-?\d+\.\d+E?-?\d*)\s*'
		# --------------JD0---------dt---------------yr------mnth-----dy------HR------MIN------SEC--------
		daterx = r'\s?(\d+\.\d?)\s+(\d+.\d+)\s+TT \(\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)H\s+(\d+)M\s+(\d+\.\d+)S UTC\)\s?'

	if (type == 'date'):
		match = re.match(daterx, str)
	elif (type == 'q'):
		match = re.match(qrx, str)
	elif (type == 'v'):
		match = re.match(vrx, str)
	else:
		print("Couldn't recognize string: ", str)
		return None

	if (match):  # Проверяем, есть ли совпадение:
		out = []  # Создаём пустой массив чтобы записать в него float и вернуть
		for s in match.groups():
			out.append(float(s))  # Преобразуем из str в float
		return out
	else:
		return None


# Определяет количество спутников в файле. Возвращает целое число.
def SatCount(filename='EPH_000.OUT'):
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


# Возвращает номер строки с координатами спутника. N_sats - всего спутников в файл, n_sat - номер спутника, rec - номер записи
def satLine(N_sats, n_sat, rec):
	return (2 * N_sats + 1) * rec + 2 * (n_sat + 1)

def stdoutToStr(bstring):
	bstring = bstring.split(b'\r\n')
	stringlist = [bs.decode() for bs in bstring]
	if len(stringlist) == 1:
		return stringlist[0]
	return stringlist[:-1]  # Обрезаем пустую строку в коце.

def runcmd(cmd):
	proc = Popen(cmd, stdout=PIPE)
	b_out = proc.stdout.read()
	return stdoutToStr(b_out)

def evo(infile, selected_sat, plot=True):
	infile = 'EPH/%s' % infile
	N_sats = SatCount(infile)
	wc_out = runcmd(['pwd'])
	wc_out = runcmd(['wc', '-l','%s' % infile])
	total_lines = int(wc_out.split(' ')[0])
	total_records = total_lines // (N_sats*2 + 1)

	sat = Sat()
	sun_sat = Sat(mu=mu_solar)
	moon_sat = Sat()

	evolution = np.ndarray(shape=(total_records, 7), dtype='float64')

	moon_res = np.ndarray(shape=(total_records, total_resonances+1), dtype='float64')
	moon_crit = np.ndarray(shape=(total_records, total_resonances+1), dtype='float64')
	sun_res = np.ndarray(shape=(total_records, total_resonances+1), dtype='float64')
	sun_crit = np.ndarray(shape=(total_records, total_resonances+1), dtype='float64')

	for rec in range(total_records):
		dateline = satLine(N_sats, 0, rec) - 1
		date_str = runcmd(['sed', '-n', '%dp' % dateline, infile])
		line = satLine(N_sats, selected_sat, rec)
		qv_str = runcmd(['sed', '-n', '%d,%dp' % (line, line + 1), infile])

		date_match = readeph(date_str[0], 'date')
		JD = date_match[0] + date_match[1] / 86400
		year = date_match[2] + date_match[3] / 12 + date_match[4] / 365

		if qv_str == []:
			# TODO: Emergency save
			pass

		q_match, v_match = readeph(qv_str[0], 'q'), readeph(qv_str[1], 'v')

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

	moon_res[:, 0] = moon_crit[:, 0] = sun_res[:, 0] = sun_crit[:, 0] = evolution[:, 1] # Time column

	a, i, mMEGNO = round(evolution[0, 2]), round(evolution[0, 4]/rad), evolution[-1, 6]
	e = evolution[-1, 3]
	if e < 1:
		sat_dir = 'EVO/%d - %d (%.1f)/'% (a, i, mMEGNO)
	else:
		sat_dir = 'EVO/%d - %d (%.1f) H/' % (a, i, mMEGNO)

	os.makedirs(sat_dir, exist_ok=True)
	os.makedirs(sat_dir + 'MOON/', exist_ok=True)
	os.makedirs(sat_dir + 'SUN/', exist_ok=True)

	np.savetxt(sat_dir + 'evo.csv' , evolution, delimiter=',', header='JD,year,a,e,i,MEGNO,mMEGNO')
	np.savetxt(sat_dir + 'moon_res.csv' , moon_res, delimiter=',')
	np.savetxt(sat_dir + 'moon_crit.csv' , moon_crit, delimiter=',')
	np.savetxt(sat_dir + 'sun_res.csv' , sun_res, delimiter=',')
	np.savetxt(sat_dir + 'sun_crit.csv' , sun_crit, delimiter=',')

	if plot:
		labels = resonances_labels('Moon')
		for expression in range(1, moon_res.shape[1]):
			f = plt.figure()
			plt.subplot(211)
			plt.scatter(moon_crit[:, 0], moon_crit[:, expression], s=1)
			plt.ylabel('Крит. аргумент, рад')
			plt.title(labels[expression-1])
			plt.xticks([], [])

			plt.subplot(212)
			plt.plot(moon_res[:, 0], moon_res[:, expression])
			if np.sign(min(moon_res[:, expression]) * max(moon_res[:, expression])) == -1:
				plt.plot([moon_res[0, 0], moon_res[0, -1]], [0, 0], c='r')

			plt.ylabel('Рез. соотн., рад/с')
			plt.xlabel('Время, годы')
			plt.xlim((moon_res[0,0], moon_res[-1, 0]))
			# f.subplots_adjust(hspace=0)
			plt.savefig(sat_dir + 'MOON/%d.png' % expression)
			plt.close()

		labels = resonances_labels('Sun')
		for i in range(1, sun_res.shape[1]):
			f = plt.figure()
			plt.subplot(211)
			plt.scatter(sun_crit[:, 0], sun_crit[:, i], s=1)
			plt.ylabel('Крит. аргумент, рад')
			plt.title(labels[i-1])
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