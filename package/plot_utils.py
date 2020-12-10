import numpy as np
import matplotlib.pyplot as plt


def plot_with_axhline(x: np.array, y: np.array, axis=None, yline=0, xmin=0, xmax=0, legends=None, title="Figure", color='red') -> None:

    if axis is None:
        plt.title(title)
        plt.plot(x, y, linewidth=1, color=color)
        # plt.axhline(y=yline, xmin=xmin, xmax=xmax, color='red', linewidth=0.5, linestyle='--')
        if legends is not None:
            plt.legend(legends)
    else:
        axis.plot(x, y, linewidth=1, color='blue')
        # axis.axhline(y=yline, xmin=xmin, xmax=xmax, color='red', linewidth=0.5, linestyle='--')
        if legends is not None:
            plt.legend(legends)


def show_plots():
    plt.show()
