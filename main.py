from math import *
from numpy import *
import matplotlib.pyplot as plt


def main():
    # test_cd()
    # test_ci()
    test_jacobiana()
    # test_tra()


def test_cd():
    # casos de la a2
    m_a = cin_dir(0, 0, 0, 0)
    p_a = [round(m_a[0][3], 2), round(m_a[1][3], 2), round(m_a[2][3], 2)]
    o_a = [[round(m_a[0][0], 2), round(m_a[0][1], 2), round(m_a[0][2], 2)],
           [round(m_a[1][0], 2), round(m_a[1][1], 2), round(m_a[1][2], 2)],
           [round(m_a[2][0], 2), round(m_a[2][1], 2), round(m_a[2][2], 2)]]
    print("a) q(0, 0, 0, 0):")
    print("     Posición: ")
    print(array(p_a))
    print("     Orientación: ")
    print(array(o_a))
    m_b = cin_dir(0, pi/2, 0, 0)
    p_b = [round(m_b[0][3], 2), round(m_b[1][3], 2), round(m_b[2][3], 2)]
    o_b = [[round(m_b[0][0], 2), round(m_b[0][1], 2), round(m_b[0][2], 2)],
           [round(m_b[1][0], 2), round(m_b[1][1], 2), round(m_b[1][2], 2)],
           [round(m_b[2][0], 2), round(m_b[2][1], 2), round(m_b[2][2], 2)]]
    print("b) q(0, pi/2, 0, 0):")
    print("     Posición: ")
    print(array(p_b))
    print("     Orientación: ")
    print(array(o_b))
    m_c = cin_dir(-pi/2, pi/2, 0, 50)
    p_c = [round(m_c[0][3], 2), round(m_c[1][3], 2), round(m_c[2][3], 2)]
    o_c = [[round(m_c[0][0], 2), round(m_c[0][1], 2), round(m_c[0][2], 2)],
           [round(m_c[1][0], 2), round(m_c[1][1], 2), round(m_c[1][2], 2)],
           [round(m_c[2][0], 2), round(m_c[2][1], 2), round(m_c[2][2], 2)]]
    print("c) q(-pi/2, pi/2, 0, 50):")
    print("     Posición: ")
    print(array(p_c))
    print("     Orientación: ")
    print(array(o_c))
    m_d = cin_dir(pi, 0, pi/2, 25)
    p_d = [round(m_d[0][3], 2), round(m_d[1][3], 2), round(m_d[2][3], 2)]
    o_d = [[round(m_d[0][0], 2), round(m_d[0][1], 2), round(m_d[0][2], 2)],
           [round(m_d[1][0], 2), round(m_d[1][1], 2), round(m_d[1][2], 2)],
           [round(m_d[2][0], 2), round(m_d[2][1], 2), round(m_d[2][2], 2)]]
    print("a) q(pi, 0, pi/2, 25):")
    print("     Posición: ")
    print(array(p_d))
    print("     Orientación: ")
    print(array(o_d))


def cin_dir(th1, th2, th4, d3):
    #   Tabla DyH:
    #   i   alfa    a       theta   d
    #   1   0       105     th1     120
    #   2   0       75      th2     0
    #   3   0       0       0       -70-d3
    #   4   0       0       th4     0
    m01 = mat_dh(0, 105, th1, 120)
    m12 = mat_dh(0, 75, th2, 0)
    m23 = mat_dh(0, 0, 0, -70-d3)
    m34 = mat_dh(0, 0, th4, 0)
    m04 = dot(dot(dot(m01, m12), m23), m34)     # M01*M12*M23*M34
    return m04


def mat_dh(alfa, a, theta, d):
    mdh = [[cos(theta), -sin(theta) * cos(alfa), sin(theta) * sin(alfa), a * cos(theta)],
           [sin(theta), cos(theta) * cos(alfa), -cos(theta) * sin(alfa), a * sin(theta)],
           [0, sin(alfa), cos(alfa), d],
           [0, 0, 0, 1]]
    return mdh


def cin_inv(x, y, z, orientacion):
    brazo = 105
    antebrazo = 75
    d3 = [50 - z, 50 - z]
    r = sqrt(x**2 + y**2)
    epsilon = atan2(y, x)
    aux1 = brazo**2 + r**2 - antebrazo**2
    aux2 = 2 * r * brazo
    fi = acos(aux1 / aux2)
    th1 = [epsilon - fi, epsilon + fi]
    th_aux = asin((r * sin(fi)) / antebrazo)
    if r < sqrt(brazo ** 2 + antebrazo ** 2):
        th_aux = pi - th_aux
    th2 = [th_aux, -th_aux]
    th4 = [orientacion - th1[0] - th2[0], orientacion - th1[1] - th2[1]]
    return [th1, th2, th4, d3]


def test_ci():
    # casos de la a2
    conf_a = cin_inv(0, 0, 0, 0)
    conf_b = cin_inv(0, 0, 0, 0)
    conf_c = cin_inv(0, 0, 0, 0)
    conf_d = cin_inv(0, 0, 0, 0)
    # mas casos


def jacobiana(th1, th2, th4, d3):
    j = [[-75*sin(th1 + th2) - 105*sin(th1), -75*sin(th1 + th2), 0, 0],
         [75*cos(th1 + th2) + 105*cos(th1), 75*cos(th1 + th2), 0, 0],
         [0, 0, -1, 0],
         [1, 1, 0, 1]]
    return j


def test_jacobiana():
    # casos de la a2
    mja = jacobiana(0, -pi/2, pi/2, 25)
    for i in range(len(mja)):
        for j in range(len(mja[i])):
            mja[i][j] = round(mja[i][j], 2)
    mjb = jacobiana(pi/4, -pi/4, 0, 0)
    for i in range(len(mja)):
        for j in range(len(mjb[i])):
            mjb[i][j] = round(mjb[i][j], 2)

    print("a) q(0, -pi/2, pi/2, 25):")
    print(array(mja))
    print("b) q(pi/4, -pi/4, 0, 0):")
    print(array(mjb))

    # mas casos
    mjp1 = jacobiana(-pi/2, pi/2, -pi/2, 0)
    for i in range(len(mja)):
        for j in range(len(mjp1[i])):
            mjp1[i][j] = round(mjp1[i][j], 2)
    mjp2 = jacobiana(pi/4, -pi/2, 0, 0)
    for i in range(len(mja)):
        for j in range(len(mjp2[i])):
            mjp2[i][j] = round(mjp2[i][j], 2)
    mjp3 = jacobiana(0, pi/2, pi/4, 30)
    for i in range(len(mjp3)):
        for j in range(len(mjp3[i])):
            mjp3[i][j] = round(mjp3[i][j], 2)

    print("Prueba1) q(-pi/2, pi/2, -pi/2, 0):")
    print(array(mjp1))
    print("Prueba2) q(pi/4, -pi/2, 0, 0):")
    print(array(mjp2))
    print("Prueba3) q(0, pi/2, pi/4, 30):")
    print(array(mjp3))


def trayectoria(x0, v_x, a_x, y0, v_y, a_y, z0, v_z, a_z):
    puntos = [[x0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [y0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [z0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
    configuraciones = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
    for t in range(1, 21):
        puntos[0][t] = x0 + v_x * t + a_x * t**2
        puntos[1][t] = y0 + v_y * t + a_y * t**2
        puntos[2][t] = z0 + v_z * t + a_z * t**2
    print(array(puntos))
    for t in range(0, 21):
        # calcular la orientacion del eje x: angulo formado por dx y dy
        dx = v_x + a_x * t
        dy = v_y + a_y * t
        angulo = atan2(dy, dx)
        config = cin_inv(puntos[0][t], puntos[1][t], puntos[2][t], angulo)
        print("config")
        print(array(config))
        """if t == 0:
            if puntos[1][t] > 0:
                configuraciones[0][t] = config[0][0]
                configuraciones[1][t] = config[1][0]
                configuraciones[2][t] = config[2][0]
                configuraciones[3][t] = config[3][0]
            else:
                configuraciones[0][t] = config[0][1]
                configuraciones[1][t] = config[1][1]
                configuraciones[2][t] = config[2][1]
                configuraciones[3][t] = config[3][1]
        elif abs(config[0][0] - configuraciones[0][t-1]) < abs(config[0][1] - configuraciones[0][t-1]):
            configuraciones[0][t] = config[0][0]
            configuraciones[1][t] = config[1][0]
            configuraciones[2][t] = config[2][0]
            configuraciones[3][t] = config[3][0]
        else:"""
        configuraciones[0][t] = config[0][1]
        configuraciones[1][t] = config[1][1]
        configuraciones[2][t] = config[2][1]
        configuraciones[3][t] = config[3][1]
    return configuraciones


def test_tra():
    config = trayectoria(-75, 8, 0, 105, 0, 0, 50, 0, 0)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(0, 21):
        print('theta1: {}, theta2: {}, theta4: {}, d3: {}'.format(round(config[0][i], 2), round(config[1][i], 2), round(config[2][i], 2), round(config[3][i], 2)))
        matriz = cin_dir(config[0][i], config[1][i], config[2][i], config[3][i])
        pos = [round(matriz[0][3], 2), round(matriz[1][3], 2), round(matriz[2][3], 2)]
        print(array(pos))
        ax.scatter(pos[0], pos[1], pos[2], c='r', marker='o')

    config = trayectoria(-75, 8, 0, 105, 10, -0.5, 0, 10, -0.5)
    for i in range(0, 21):
        print('theta1: {}, theta2: {}, theta4: {}, d3: {}'.format(round(config[0][i], 2), round(config[1][i], 2), round(config[2][i], 2), round(config[3][i], 2)))
        matriz = cin_dir(config[0][i], config[1][i], config[2][i], config[3][i])
        pos = [round(matriz[0][3], 2), round(matriz[1][3], 2), round(matriz[2][3], 2)]
        print(array(pos))
        ax.scatter(pos[0], pos[1], pos[2], c='b', marker='^')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
