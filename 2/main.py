import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

### ИЗМЕНЯЕМЫЕ ДАННЫЕ

X_C = 5 # координаты центра диска
Y_C = 5
R = 4 # радиус диска
A = 2 # расстояние между шарниром и центром диска
L = 2 # длина стержня, на котором шарнирно прикреплён груз

### СТАТИЧЕСКАЯ ЧАСТЬ

ang = np.linspace(0, 2*math.pi, 80) # углы для отрисовки кругов
X_Disk = X_C + R*np.cos(ang) # координаты диска
Y_Disk = Y_C + R*np.sin(ang)
X_Sm = X_C + 0.2*np.cos(ang) # координаты маленького круга в центре диска
Y_Sm = Y_C + 0.2*np.sin(ang)

X_Side_1 = [X_C+0.2*np.cos(math.pi*5/4), X_C+0.5*np.cos(math.pi*5/4)] # боковые линии (центр)
Y_Side_1 = [X_C+0.2*np.sin(math.pi*5/4), Y_C+0.5*np.sin(math.pi*5/4)]
X_Side_2 = [X_C+0.2*np.cos(math.pi/-4), X_C+0.5*np.cos(math.pi/-4)]
Y_Side_2 = [X_C+0.2*np.sin(math.pi/-4), Y_C+0.5*np.sin(math.pi/-4)]

X_Bottom = [X_Side_1[1]-0.1, X_Side_2[1]+0.1] # линия-закреп центра
Y_Bottom = [Y_Side_1[1], Y_Side_2[1]]

X_Lines_1 = np.linspace(float(X_Bottom[0])+0.05, float(X_Bottom[1])-0.05, 5) # полоски на линии-закрепа центра диска
X_Lines_2 = X_Lines_1 + 0.3*np.cos(math.pi*9/8)
Y_Lines_1 = np.full(5, Y_Bottom[0])
Y_Lines_2 = Y_Lines_1 + 0.3*np.sin(math.pi*9/8)

### ДИНАМИЧЕСКАЯ ЧАСТЬ

Steps = 250
t_fin = 5
t = np.linspace(0, t_fin, Steps) # время
phi = np.zeros_like(t) # угол между вертикальной осью и радиус-вектором к шарниру
psi = np.zeros_like(t) # угол между вертикальной осью и стержнем
X_Sh = np.zeros_like(t) # координаты шарнира
Y_Sh = np.zeros_like(t)
X_Gr = np.zeros_like(t) # координаты груза
Y_Gr = np.zeros_like(t)

for i in np.arange(len(t)): # просчёт основных величин
    phi[i] = 2*np.sin(1.7*t[i]) + 5*np.cos(1.2*t[i])
    psi[i] = 2*np.sin(1.7*t[i]) + 5*np.cos(1.2*t[i])
    X_Sh[i] = X_C + A*np.cos(phi[i]+math.pi/2)
    Y_Sh[i] = Y_C + A*np.sin(phi[i]+math.pi/2)
    X_Gr[i] = X_Sh[i] + L*np.cos(-psi[i]-math.pi/2)
    Y_Gr[i] = Y_Sh[i] + L*np.sin(-psi[i]-math.pi/2)

### ПЕРЕХОД К ОТРИСОВКЕ
    
fig = plt.figure() # задаём пространство для отрисовки
ax = fig.add_subplot(1, 1, 1)
ax.axis('equal')
ax.set(xlim=[0, 10], ylim=[0, 10])

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
    
Len = 0.4 # длина линии-закрепа пружинки
Wide = 0.2 # ширина линии-закрепа пружинки
X_Dsht = 0.1*np.cos(math.pi/4) # сдвиги штрихов по координатам
Y_Dsht = 0.1*np.cos(math.pi/4)
R1 = 0.3 # радиусы спиральной пружины
R2 = 0.1

thetta = np.linspace(0, max(3/2*math.pi,3/2*math.pi+psi[0]), 100) # угол проворота спиральной пружины
X_SpiralSpr = (R1 + thetta*(R2-R1)/thetta[-1])*np.cos(thetta) # координаты точек спиральной пружины
Y_SpiralSpr = -(R1 + thetta*(R2-R1)/thetta[-1])*np.sin(thetta)

spr, = ax.plot(X_SpiralSpr+X_Sh[0], Y_SpiralSpr+Y_Sh[0], color = 'green') # отрисовка спиральной пружины
pl1, = ax.plot([X_Sh[0]+R1-Wide/2, X_Sh[0]+R1-Wide/2+X_Dsht], [Y_Sh[0]+Len, Y_Sh[0]+Len+Y_Dsht], color = 'darkgreen') # штрихи на линии-закрепе спиральки
pl2, = ax.plot([X_Sh[0]+R1, X_Sh[0]+R1+X_Dsht], [Y_Sh[0]+Len, Y_Sh[0]+Len+Y_Dsht], color = 'darkgreen')
pl3, = ax.plot([X_Sh[0]+R1+Wide/2, X_Sh[0]+R1+Wide/2+X_Dsht], [Y_Sh[0]+Len, Y_Sh[0]+Len+Y_Dsht], color = 'darkgreen')
hl, = ax.plot([X_Sh[0]+R1-Wide/2, X_Sh[0]+R1-Wide/2], [Y_Sh[0]+Len, Y_Sh[0]+Len]) # отрисовка вертикальной линии от спиральки
upl, = ax.plot([X_Sh[0]+R1, X_Sh[0]+R1], [Y_Sh[0], Y_Sh[0]+Len], color = 'green') # отрисовка линии-закрепа спирали
sh, = ax.plot(X_Sh[0], Y_Sh[0], marker='o', markersize = 5, color = 'orange') # отрисовка шарнира
st, = ax.plot([X_Sh[0], X_Gr[0]], [Y_Sh[0], Y_Gr[0]], color = 'orange') # отрисовка стержня
gr, = ax.plot(X_Gr[0], Y_Gr[0], marker = 'o', markersize = 20, color = 'orange') # отрисовка грузика

def anima(i):
    thetta = np.linspace(0, max(3/2*math.pi,3/2*math.pi+psi[i]), 100)
    X_SpiralSpr = (R1 + thetta*(R2-R1)/thetta[-1])*np.cos(thetta)
    Y_SpiralSpr = -(R1 + thetta*(R2-R1)/thetta[-1])*np.sin(thetta)
    spr.set_data(X_SpiralSpr+X_Sh[i], Y_SpiralSpr+Y_Sh[i])
    pl1.set_data([X_Sh[i]+R1-Wide/2, X_Sh[i]+R1-Wide/2+X_Dsht], [Y_Sh[i]+Len, Y_Sh[i]+Len+Y_Dsht])
    pl2.set_data([X_Sh[i]+R1, X_Sh[i]+R1+X_Dsht], [Y_Sh[i]+Len, Y_Sh[i]+Len+Y_Dsht])
    pl3.set_data([X_Sh[i]+R1+Wide/2, X_Sh[i]+R1+Wide/2+X_Dsht], [Y_Sh[i]+Len, Y_Sh[i]+Len+Y_Dsht])
    hl.set_data([X_Sh[i]+R1-Wide/2, X_Sh[i]+R1-Wide/2], [Y_Sh[i]+Len, Y_Sh[i]+Len])
    upl.set_data([X_Sh[i]+R1, X_Sh[i]+R1], [Y_Sh[i], Y_Sh[i]+Len])
    sh.set_data(X_Sh[i], Y_Sh[i])
    st.set_data([X_Sh[i], X_Gr[i]], [Y_Sh[i], Y_Gr[i]])
    gr.set_data(X_Gr[i], Y_Gr[i])
    return spr, pl1, pl2, pl3, hl, upl, sh, st, gr

anim = FuncAnimation(fig, anima, frames=Steps, interval=50, repeat=False) # создаём разовую анимацию
fig.suptitle('Kutsenko LW2', fontsize=14) # добавляем название
anim.save("./Animation.mp4", writer="ffmpeg") # сохраняем анимацию
plt.close()
plt.show()