import numpy as np
import matplotlib.pyplot as plt

# # prepare some data
# x = np.linspace(-1, 1, num=100)
# y = np.arccos(x)

# figure, axes = plt.subplots(1, 2, constrained_layout=False)

# axes[0].plot(x, y)
# axes[0].set_title("Title for plot 01", fontsize=16)
# axes[0].set_xlabel("X", color="red")
# axes[0].set_xlim([-2, 2])
# ticks = np.linspace(-2, 2, num=10)
# axes[0].set_xticks(ticks)
# axes[0].set_xticklabels([f"{value:.2f}" for value in ticks])

# axes[0].axhline(y=1, xmin=-2, xmax=2, color = 'red', linewidth = 0.5, linestyle='--')

x = np.linspace(-2*np.pi, 2*np.pi, num=100, endpoint=True)
y = np.sin(x)
plt.plot(x, y, color='red', linewidth=1.5)
plt.
plt.show()
