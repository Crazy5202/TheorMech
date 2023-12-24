import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sympy as sp
from scipy.integrate import odeint

def odesys(y, t, r, a, l, m1, m2, c, g): # Функция создания системы диффуров
    dy = np.zeros(4)
    dy[0] = y[2]
    dy[1] = y[3]
    a11 = m1*r**2+2*m2*a**2
    a12 = 2*m2*a*l*np.cos(y[0]+y[1])
    a21 = m2*l*a*np.cos(y[0]+y[1])
    a22 = m2*l**2
    b1 = 2*m2*a*l*np.sin(y[0]+y[1])*y[3]**2+2*m2*g*a*np.sin(y[0])-2*c*(y[0]+y[1])
    b2 = m2*l*a*np.sin(y[0]+y[1])*y[2]**2-m2*g*l*np.sin(y[1])-c*(y[0]+y[1])
    dy[2] = (b1*a22-b2*a12)/(a11*a22-a12*a21)
    dy[3] = (b2*a11-b1*a21)/(a11*a22-a12*a21)
    return dy

### ИЗМЕНЯЕМЫЕ ПАРАМЕТРЫ СИСТЕМЫ

R = 10 # радиус диска
A = 2 # расстояние между шарниром и центром диска
L = 10 # длина стержня, на котором шарнирно прикреплён груз
M1 = 10 # масса диска
M2 = 10 # масса груза
C = 200 # жёсткость спиральной пружины
G = 9.81 # ускорение свободного падения

K = 1 # величина для расчёта R (график)

### НАЧАЛЬНЫЕ ЗНАЧЕНИЯ

PHI0 = math.pi
PSI0 = 0
DPHI0 = 5
DPSI0 = 6
Y0 = [PHI0, PSI0, DPHI0, DPSI0]

### СТАТИЧЕСКАЯ ЧАСТЬ

X_C = R+A+L # координаты центра диска
Y_C = R+A+L
RM = R/20 # радиус маленького круга в центре диска

ang = np.linspace(0, 2*math.pi, 80) # углы для отрисовки кругов
X_Disk = X_C + R*np.cos(ang) # координаты диска
Y_Disk = Y_C + R*np.sin(ang)
X_Sm = X_C + RM*np.cos(ang) # координаты маленького круга в центре диска
Y_Sm = Y_C + RM*np.sin(ang)

X_Side_1 = [X_C+RM*np.cos(math.pi*5/4), X_C+RM*3*np.cos(math.pi*5/4)] # боковые линии (центр)
Y_Side_1 = [X_C+RM*np.sin(math.pi*5/4), Y_C+RM*3*np.sin(math.pi*5/4)]
X_Side_2 = [X_C+RM*np.cos(math.pi/-4), X_C+RM*3*np.cos(math.pi/-4)]
Y_Side_2 = [X_C+RM*np.sin(math.pi/-4), Y_C+RM*3*np.sin(math.pi/-4)]

X_Bottom = [X_Side_1[1]-R/40, X_Side_2[1]+R/40] # линия-закреп центра
Y_Bottom = [Y_Side_1[1], Y_Side_2[1]]

X_Lines_1 = np.linspace(float(X_Bottom[0])+R/50, float(X_Bottom[1])-R/50, 5) # полоски на линии-закрепа центра диска
X_Lines_2 = X_Lines_1 + R/20*np.cos(math.pi*9/8)
Y_Lines_1 = np.full(5, Y_Bottom[0])
Y_Lines_2 = Y_Lines_1 + R/20*np.sin(math.pi*9/8)

### ДИНАМИЧЕСКАЯ ЧАСТЬ

Steps = 1000
t_fin = 20
t = np.linspace(0, t_fin, Steps) # время
X_Sh = np.zeros_like(t) # координаты шарнира
Y_Sh = np.zeros_like(t)
X_Gr = np.zeros_like(t) # координаты груза
Y_Gr = np.zeros_like(t)

Sol = odeint(odesys, Y0, t, (R, A, L, M1, M2, C, G)) # решение диффура

phi = Sol[:,0] # угол между вертикальной осью и радиус-вектором к шарниру
psi = Sol[:,1] # угол между вертикальной осью и стержнем
dphi = Sol[:,2] # угловые скорости
dpsi = Sol[:,3]
ddphi = [odesys(y, t, R, A, L, M1, M2, C, G)[2] for y,t in zip(Sol,t)] # угловые ускорения
ddpsi = [odesys(y, t, R, A, L, M1, M2, C, G)[3] for y,t in zip(Sol,t)]

Rp = -K*dphi # сила сопротивления

for i in np.arange(len(t)): # просчёт основных величин
    X_Sh[i] = X_C + A*np.cos(phi[i]+math.pi/2)
    Y_Sh[i] = Y_C + A*np.sin(phi[i]+math.pi/2)
    X_Gr[i] = X_Sh[i] + L*np.cos(-psi[i]-math.pi/2)
    Y_Gr[i] = Y_Sh[i] + L*np.sin(-psi[i]-math.pi/2)

### ПЕРЕХОД К ОТРИСОВКЕ
    
fig = plt.figure() # задаём пространство для отрисовки
ax = fig.add_subplot(1, 1, 1)
ax.axis('equal')
ax.set(xlim = [0, X_C+R+A+L], ylim = [0, Y_C+R+A+L])
ax.set(xlabel="x", ylabel="y")

### СТАТИЧЕСКАЯ ОТРИСОВКА

ax.plot(X_C, Y_C, marker = 'o', markersize=2, color = 'blue') # отрисовка центра диска
ax.plot(X_Disk, Y_Disk, color = 'blue') # отрисовка диска
ax.plot(X_Sm, Y_Sm, color = 'blue') # отрисовка кружка вокруг центра диска
ax.plot(X_Side_1, Y_Side_1, color = 'blue') # отрисовка боковых линий от центра диска
ax.plot(X_Side_2, Y_Side_2, color = 'blue')
ax.plot(X_Bottom, Y_Bottom, color = 'blue') # отрисовка линии-закрепа центра диска
for i in np.arange(len(X_Lines_1)): # отрисовка штрихов на линии-закрепе центра диска
    ax.plot([X_Lines_1[i], X_Lines_2[i]], [Y_Lines_1[i], Y_Lines_2[i]], color = 'darkblue')

### ДИНАМИЧЕСКАЯ ОТРИСОВКА
    
LEN = R/6 # длина линии-закрепа пружинки
WIDE = R/8 # ширина линии-закрепа пружинки
X_DSHT = R/20*np.cos(math.pi/4) # сдвиги штрихов по координатам
Y_DSHT = R/20*np.cos(math.pi/4)
R1 = R/8 # радиусы спиральной пружины
R2 = R/64

thetta = np.linspace(0, 3/2*math.pi+psi[0], 100) # угол проворота спиральной пружины
X_SpiralSpr = (R1 + thetta*(R2-R1)/thetta[-1])*np.cos(thetta) # координаты точек спиральной пружины
Y_SpiralSpr = -(R1 + thetta*(R2-R1)/thetta[-1])*np.sin(thetta)

spr, = ax.plot(X_SpiralSpr+X_Sh[0], Y_SpiralSpr+Y_Sh[0], color = 'green') # отрисовка спиральной пружины
pl1, = ax.plot([X_Sh[0]+R1-WIDE/2, X_Sh[0]+R1-WIDE/2+X_DSHT], [Y_Sh[0]+LEN, Y_Sh[0]+LEN+Y_DSHT], color = 'darkgreen') # штрихи на линии-закрепе спиральки
pl2, = ax.plot([X_Sh[0]+R1, X_Sh[0]+R1+X_DSHT], [Y_Sh[0]+LEN, Y_Sh[0]+LEN+Y_DSHT], color = 'darkgreen')
pl3, = ax.plot([X_Sh[0]+R1+WIDE/2, X_Sh[0]+R1+WIDE/2+X_DSHT], [Y_Sh[0]+LEN, Y_Sh[0]+LEN+Y_DSHT], color = 'darkgreen')
hl, = ax.plot([X_Sh[0]+R1-WIDE/2-R/10, X_Sh[0]+R1+WIDE/2+R/10], [Y_Sh[0]+LEN, Y_Sh[0]+LEN], color = 'green') # отрисовка вертикальной линии от спиральки
upl, = ax.plot([X_Sh[0]+R1, X_Sh[0]+R1], [Y_Sh[0], Y_Sh[0]+LEN], color = 'green') # отрисовка линии-закрепа спирали
sh, = ax.plot(X_Sh[0], Y_Sh[0], marker='o', markersize = 5, color = 'orange') # отрисовка шарнира
st, = ax.plot([X_Sh[0], X_Gr[0]], [Y_Sh[0], Y_Gr[0]], color = 'orange') # отрисовка стержня
gr, = ax.plot(X_Gr[0], Y_Gr[0], marker = 'o', markersize = 20, color = 'orange') # отрисовка грузика

def anima(i): # функция анимации
    thetta = np.linspace(0, 3/2*math.pi+psi[i], 100)
    X_SpiralSpr = (R1 + thetta*(R2-R1)/thetta[-1])*np.cos(thetta)
    Y_SpiralSpr = -(R1 + thetta*(R2-R1)/thetta[-1])*np.sin(thetta)
    spr.set_data(X_SpiralSpr+X_Sh[i], Y_SpiralSpr+Y_Sh[i])
    pl1.set_data([X_Sh[i]+R1-WIDE/2, X_Sh[i]+R1-WIDE/2+X_DSHT], [Y_Sh[i]+LEN, Y_Sh[i]+LEN+Y_DSHT])
    pl2.set_data([X_Sh[i]+R1, X_Sh[i]+R1+X_DSHT], [Y_Sh[i]+LEN, Y_Sh[i]+LEN+Y_DSHT])
    pl3.set_data([X_Sh[i]+R1+WIDE/2, X_Sh[i]+R1+WIDE/2+X_DSHT], [Y_Sh[i]+LEN, Y_Sh[i]+LEN+Y_DSHT])
    hl.set_data([X_Sh[i]+R1-WIDE/2, X_Sh[i]+R1+WIDE/2], [Y_Sh[i]+LEN, Y_Sh[i]+LEN])
    upl.set_data([X_Sh[i]+R1, X_Sh[i]+R1], [Y_Sh[i], Y_Sh[i]+LEN])
    sh.set_data(X_Sh[i], Y_Sh[i])
    st.set_data([X_Sh[i], X_Gr[i]], [Y_Sh[i], Y_Gr[i]])
    gr.set_data(X_Gr[i], Y_Gr[i])
    return spr, pl1, pl2, pl3, hl, upl, sh, st, gr

anim = FuncAnimation(fig, anima, frames=Steps, interval=40, repeat=False) # создаём разовую анимацию
fig.suptitle('Kutsenko LW3', fontsize=14) # добавляем название
anim.save("./Animation.mp4", writer="ffmpeg") # сохраняем анимацию

### ГРАФИКИ ЗАВИСИМОСТЕЙ ВЕЛИЧИН ОТ ВРЕМЕНИ

pls = plt.figure()
p1 = pls.add_subplot(3, 1, 1) # строим графики величин
p1.set(xlim = [0,t_fin])
p1.plot(t, phi)
p1.grid()
plt.title('Phi(t)')

p2 = pls.add_subplot(3, 1, 2)
p2.plot(t, psi)
p2.grid()
p2.set(xlim = [0,t_fin])
plt.title('Psi(t)')

p3 = pls.add_subplot(3, 1, 3)
p3.plot(t, Rp)
p3.grid()
p3.set(xlim = [0,t_fin])
plt.title('R(t)')

plt.tight_layout() # чтобы не накладывались названия
plt.savefig('Plots.png') # сохраняем графики