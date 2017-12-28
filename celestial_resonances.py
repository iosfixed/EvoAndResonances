# coding=utf-8
from numpy import pi, sin, cos, sqrt
from celestial_main import Sat, siderial
rad = pi/180
# Константы для вычисления резонансных аргументов:
dM_L = (1739527262.8478*rad)/(3600.0*36525*86400)
dM_S = (129602768.13*rad)/(3600.0*36525*86400)

dw_L = 2.25e-8
dw_S = 9.2e-14

dW_L = (6962890.5431*rad)/(3600.0*36525*86400)
dW_S = 3.32e-10

inc_S = 23.45*rad
inc_L = 23.45*rad

Moon2Earth = 1/81.3005690699  # Отношение масс Луна/Земля
Sun2Earth = 332946.048166  # Солнце/Земля

a_L = 384748
a_S = 149597868

J20 = 1.0826359e-3
r0 = 6363672.6e-3 # Средний экваториальный? радиус

# Критические аргументы для апсидально-нодальных резонансов
def critical(sat1, body, type='rad'):
    dW = sat1.W - body.W
    wSat = sat1.w
    wBody = body.w

    out = []

    out.append(dW + wSat - wBody)
    out.append(dW - wSat + wBody)
    out.append(dW + wSat + wBody) # 3

    out.append(dW + 2*wSat - 2*wBody)
    out.append(dW - 2 * wSat + 2 * wBody)
    out.append(dW + 2 * wSat + 2 * wBody) # 6

    out.append(dW + wSat)
    out.append(dW + 2*wSat) # 8

    out.append(wSat - wBody)
    out.append(wSat + wBody)
    out.append(wSat) # 11

    #  Причёсываем выхлоп:
    for i, r in enumerate(out):
        while (r < 0):
            r += 2*pi
        while (r > 2*pi):
            r -= 2*pi
        if (type == 'deg'): # Криво, но по умолчанию радианы
            r *= 180/pi
        out[i] = r

    return out

def resonances(sat, body='sun'):
    a, e, i = sat.a, sat.e, sat.i   # Нечего тысячу раз дёргать параметры объекта
    n = sqrt(sat.mu * a**(-3))      # Среднее движение спутника

    # Все величиниы без указания тела относятся к текущему спутнику:
    WJ2 = -3/2 * J20 * n * r0**2 / a**2 * cos(i) * (1 - e**2)**-2
    WL = -3/16 * n * Moon2Earth * (a / a_L)**3 * (2 + 3*e**2)/sqrt(1 - e**2) * (2 - 3*sin(inc_L)**2) * cos(i)
    WS = -3/16 * n * Sun2Earth * (a / a_S)**3 * (2 + 3*e**2)/sqrt(1 - e**2) * (2 - 3*sin(inc_S)**2) * cos(i)

    WSat = WJ2 + WL + WS

    wJ2 = 3/4 * J20 * n * (r0/a)**2 * (5*cos(i)**2 - 1)*(1 - e**2)**(-2)
    wL = 3/16 * n * Moon2Earth * (a/a_L)**3 * (4 - 5*sin(i)**2 + e**2)/sqrt(1 - e**2) * (2 - 3*sin(inc_L)**2)
    wS = 3/16 * n * Sun2Earth * (a/a_S)**3 * (4 - 5*sin(i)**2 + e**2)/sqrt(1 - e**2) * (2 - 3*sin(inc_S)**2)

    wSat = wJ2 + wL + wS

    
    if body == 'sun':
        dW = WSat - WS  # D for delta
        wBody = wS
    else:
        dW = WSat - WL  # D for delta
        wBody = wL


    out = []
    out.append(dW + wSat - wBody)
    out.append(dW - wSat + wBody)
    out.append(dW + wSat + wBody)  # 3

    out.append(dW + 2 * wSat - 2 * wBody)
    out.append(dW - 2 * wSat + 2 * wBody)
    out.append(dW + 2 * wSat + 2 * wBody)  # 6

    out.append(dW + wSat)
    out.append(dW + 2 * wSat)  # 8

    out.append(wSat - wBody)
    out.append(wSat + wBody)
    out.append(wSat)  # 11


    return out

def resonances_labels(body):
	out = []
	out.append('dW + wSat - w%s' % body)
	out.append('dW - wSat + %s' % body)
	out.append('dW + wSat + %s' % body)  # 3

	out.append('dW + 2 * wSat - 2 * %s' % body)
	out.append('dW - 2 * wSat + 2 * %s' % body)
	out.append('dW + 2 * wSat + 2 * %s' % body)  # 6

	out.append('dW + wSat')
	out.append('dW + 2 * wSat')  # 8

	out.append('wSat - %s' % body)
	out.append('wSat + %s' % body)
	out.append('wSat')  # 11

	return out