import os
import re
from subprocess import Popen, PIPE
import numpy as np
from numpy import sqrt, pi, dot, sin, cos, arctan2, cross, arcsin, arctanh
import pandas as pd

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
rad = pi / 180
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

def stdoutToStr(bstring):
	bstring = bstring.split(b'\r\n')
	stringlist = [bs.decode() for bs in bstring]
	return stringlist[:-1]  # Обрезаем пустую строку в коце.

def runcmd(cmd):
	proc = Popen(cmd, stdout=PIPE)
	return stdoutToStr(proc.stdout.read())
sat = Sat()
files = os.chdir('EPH')
files = os.listdir()
files = list(filter(lambda fname: re.match(r'EPH.*OUT', fname), files))

files.sort()

satlist = pd.DataFrame(columns=['file', 'sat_N', 'a', 'i', 'mMEGNO'])
for infile in files:
	N_sats = SatCount(infile)

	Houtput = runcmd(['head', '-n %d' % (2 * N_sats + 1), infile])
	Toutput = runcmd(['tail', '-n %d' % (2 * N_sats), infile])

	date_match = readeph(Houtput[0], 'date')
	JD = date_match[0]

	for sat_N in range(N_sats):
		q0_str = Houtput[1 + 2 * sat_N]
		v0_str = Houtput[2 + 2 * sat_N]

		q0, v0 = readeph(q0_str, 'q'), readeph(v0_str, 'v')
		q0, v0 = q0[1:4], v0[0:3]
		qv = q0 + v0

		sat.set(qv, 'qv')
		init_a, init_i = round(sat.a), round(sat.i/rad)

		q1_str = Toutput[2 * sat_N]
		v1_str = Toutput[2 * sat_N + 1]

		q1, v1 = readeph(q1_str, 'q'), readeph(v1_str, 'v')
		if v1 == None or q1 == None:
			break
		qv1 = q1[1:4] + v1[0:3]

		sat.set(qv1, 'qv')
		e, mMEGNO = sat.e, v1[3]
		satlist = satlist.append({'file': infile,
		                'sat_N': sat_N,
		                'a': init_a,
		                'i': init_i,
		                'e': e,
		                'mMEGNO': mMEGNO}, ignore_index=True)

satlist.to_csv('../satlist.csv')