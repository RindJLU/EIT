# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# define some constant:
i = 1j
n = pow(10, 19)
e = 1.6*pow(10, -19)
h = 1.05*pow(10, -34)
gama = pow(10, 5)
int_ab = 1.0*pow(10, -30)
eps0 = 8.854*pow(10, -12)
Gama = pow(10, 8)
delta = pow(10, 6)*80
OmegaC = pow(10, 8)

x = np.arange(-3, 3, 0.0001)

# eps_r =(i*n*pow(int_ab, 2)*gama)/(eps0*h*(1 + i*x)) #
# eps_r_dl = i/(1 + i*x) #

def find_max(vector, a): # find the max of a matrix, as well as get the index in the vector
    if a == '1': # only one peak value--two-energy level
        v_max = np.max(vector) # v_max is the biggest number of a vector.
        delta = 0.000000001
        for n in range(len(vector)):
            if abs(vector[n] - v_max) <= delta:
                return(v_max, n)
                break
    elif a == '2': # two peak values--three-energy level
        leng = len(vector)
        half_leng = int(len(vector)/2)
        vector1 = vector[0:half_leng] # split the vector equally to find the max of each part.
        vector2 = vector[half_leng:leng]
        (v1_max, n1) = find_max(vector1, '1') # recursion
        (v2_max, n2) = find_max(vector2, '1')
        n2 = n2 + half_leng
        return(v1_max, n1, v2_max, n2)

# Main =================================================================================================================
print('=================================================================')
print('Please select a mode: 1: EIT, 2: Parameters change')
b = input()
if b == '1':
    print('Please select a mode: 1(Two-energy level); 2(Three-energy level)')
    a = input()
    print('Start------------------------------------------------------------')
    if a == '1':
        print('二能级系统')
        chi = Gama * (i * n * pow(int_ab, 2) / (eps0 * h)) / ((i * x + 1 / 2))

    elif a == '2':
        print('三能级系统')
        chi = Gama * (i * n * pow(int_ab, 2) / (eps0 * h)) * (i * (x * Gama - delta + gama / 2) / (
        (i * x + 1 / 2) * (i * (x * Gama - delta) + gama / 2) + OmegaC * OmegaC / Gama))

    else:
        print('Input error!')

    norm_factor = max(np.real(chi))/0.7 # normization factor
    chi = chi/norm_factor
    print(chi)
    if a == '1':
        (im_chi_max, n) = find_max(np.imag(chi), a)
        print('相对介电常数虚部最大值（相对）：' + str(im_chi_max) + ', 最大值出现在： 相对失谐 = ' + str(x[n]))
    elif a == '2':
        (im_chi_max1, n1, im_chi_max2, n2) = find_max(np.imag(chi), a)
        print('相对介电常数虚部最大值1（相对）：' + str(im_chi_max1) + ', 最大值出现在： 相对失谐 = ' + str(x[n1]))
        print('相对介电常数虚部最大值2（相对）：' + str(im_chi_max2) + ', 最大值出现在： 相对失谐 = ' + str(x[n2]))

    plt.plot(x, np.imag(chi), color='red', ls='-.')
    plt.plot(x, 1 + np.real(chi), color='blue', ls='-')
    plt.legend('虚实')

    plt.grid(True)
    plt.title('失谐和相对介电常数的关系')
    plt.xlabel('相对失谐')
    plt.ylabel('相对介电常数实部和虚部')
    plt.show()
if b == '2':
    print('select parameters: 1:gama_h, 2:delta')
    p = input()
    if p == '1':
        para_change = np.arange(0.1, 50, 0.2)
        gama = gama*para_change
        x1_list = []
        x2_list = []
        for g in gama:
            chi = Gama * (i * n * pow(int_ab, 2) / (eps0 * h)) * (i * (x * Gama - delta + g / 2) / (
                (i * x + 1 / 2) * (i * (x * Gama - delta) + g / 2) + OmegaC * OmegaC / Gama))
            (im_chi_max1, n1, im_chi_max2, n2) = find_max(np.imag(chi), '2')
            x1_list.append(x[n1])
            x2_list.append(x[n2])
        grad1 = (x1_list[len(x1_list) - 1] - x1_list[0]) / (gama[len(gama) - 1] - gama[0])
        grad2 = (x2_list[len(x2_list) - 1] - x2_list[0]) / (gama[len(gama) - 1] - gama[0])

        grad1 = int(grad1*pow(10, 12))/pow(10, 12)
        grad2 = int(grad2 * pow(10, 12)) / pow(10, 12)
        print(grad1, grad2)
        grad_total = grad2 - grad1
        plt.grid(True)

        plt.subplot(121)
        plt.title('左极值位置随gama_h的关系（斜率：' + str(grad1) + '）')
        plt.scatter(gama, x1_list)
        plt.text(10, .025, r'$mu=100, sigma=15$')
        plt.xlabel('gama_h')
        plt.ylabel('极值位置')

        plt.subplot(122)
        plt.scatter(gama, x2_list)
        plt.title('右极值位置随gama_h的关系（斜率：' + str(grad2) + '）')
        plt.xlabel('gama_h')
        plt.ylabel('极值位置')
        plt.show()
    elif p == '2': # delta 的影响
        delta_list = np.arange(-50, 50, 1)*pow(10, 6)
        x1_list = []
        x2_list = []
        chi_tuple = ()
        for d in delta_list:
            chi = Gama * (i * n * pow(int_ab, 2) / (eps0 * h)) * (i * (x * Gama - d + gama / 2) / (
                (i * x + 1 / 2) * (i * (x * Gama - d) + gama / 2) + OmegaC * OmegaC / Gama))
            chi_tuple = chi_tuple + (chi, )
            (im_chi_max1, n1, im_chi_max2, n2) = find_max(np.imag(chi), '2')
            x1_list.append(x[n1])
            x2_list.append(x[n2])
            x_t = np.array(x2_list) - np.array(x1_list)

        grad1 = (x1_list[len(x1_list) - 1] - x1_list[0]) / (delta_list[len(delta_list) - 1] - delta_list[0])
        grad2 = (x2_list[len(x2_list) - 1] - x2_list[0]) / (delta_list[len(delta_list) - 1] - delta_list[0])
        grad1 = int(grad1 * pow(10, 12)) / pow(10, 12)
        grad2 = int(grad2 * pow(10, 12)) / pow(10, 12)
        print(grad1, grad2)
        grad_total = grad2 - grad1


        plt.grid(True)

        plt.subplot(121)
        plt.title('左极值位置随delta变化（斜率：' + str(grad1) + '）')
        plt.scatter(delta_list, x1_list)
        plt.text(10, .025, r'$mu=100, sigma=15$')
        plt.xlabel('delta')
        plt.ylabel('极值位置')

        plt.subplot(122)
        plt.scatter(delta_list, x2_list)
        plt.title('右极值位置随delta变化（斜率：' + str(grad2) + '）')
        plt.xlabel('delta')
        plt.ylabel('极值位置')
        plt.show()

        # 缀饰态
        plt.figure()
        plt.plot(delta_list, x_t)

        # y = (delta_list+0.5*pow(10, 6))/(2*Gama)
        # plt.plot(delta_list, 2*(1+0.5*pow(y, 2)), 'd')

        plt.title('两种方法求三能级系统两峰值位置差')
        plt.xlabel('delta')
        plt.ylabel('两峰值对应探测场频率差值')

        plt.show()

        # PLOT THREE DISTINGERISHED PICTURES:
        # plt.figure()
        # plt.plot(x, np.imag(chi_tuple[0]), ls='-')
        # plt.plot(x, np.imag(chi_tuple[int(len(chi_tuple)/2)]), ls='-.')
        # plt.plot(x, np.imag(chi_tuple[len(chi_tuple) - 1]), ls=':')
        #
        # plt.legend(('delta = -0.1*omega(ab)', 'delta = 0', 'delta = 0.1*omega(ab)'))
        #
        # plt.grid(True)
        # plt.title('不同delta相对介电常数虚部图像')
        # plt.xlabel('相对失谐')
        # plt.ylabel('相对介电常数虚部')
        # plt.show()
print('End--------------------------------------------------------------')
print('=================================================================')
# ======================================================================================================================

# gama_v = gama*np.arange(0.5, 1.5, 0.3)
# for m in range(len(gama_v)):
#     chi = Gama * (i * n * pow(int_ab, 2) / (eps0 * h)) * (i * (x * Gama - delta + gama_v[m] / 2) / (
#     (i * x + 1 / 2) * (i * (x * Gama - delta) + gama_v[m] / 2) + OmegaC * OmegaC / Gama))
#     norm_factor = max(np.real(chi)) / 0.7
#     plt.subplot(2, 2, m+1)
#     plt.plot(x, np.imag(chi) / norm_factor, ls='-.')
#     plt.xlabel('相对失谐')
#     plt.ylabel('相对介电常数虚部')
# plt.show()