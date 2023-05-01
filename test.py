import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def model(x, a, b, c):
    x1, x2 = x
    return a * np.exp(-b * x1) + c * x2

# 生成模拟数据
np.random.seed(42)
x1_data = np.linspace(0, 4, 50)
x2_data = np.linspace(0, 4, 50)
X1, X2 = np.meshgrid(x1_data, x2_data)
Y = 3.0 * np.exp(-1.0 * X1) + 2.0 * X2 + np.random.normal(size=X1.shape, scale=0.2)

# 将数据转换为1D数组
x1_data_1d = X1.ravel()
x2_data_1d = X2.ravel()
y_data_1d = Y.ravel()

popt, pcov = curve_fit(model, (x1_data_1d, x2_data_1d), y_data_1d)
a, b, c = popt
print("拟合参数: a =", a, "b =", b, "c =", c)

fig = plt.figure()
ax = fig.gca(projection="3d")

# 绘制原始数据
ax.scatter(X1, X2, Y, c='blue', label='原始数据')

# 绘制拟合曲面
Y_fit = model((X1, X2), a, b, c)
ax.plot_surface(X1, X2, Y_fit, alpha=0.5, color='red', label='拟合曲面')

ax.set_xlabel('X1')
ax.set_ylabel('X2')
ax.set_zlabel('Y')
#plt.legend()
plt.savefig('test.pdf')
