import random
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from package.utility import TypeLatLong, TypeOE, getDiffAngle, normalize
from package.orbital import Simulation, OrbitCalculate
from package.orbital_elements import calc_oe_from_sv
from package.plot_utils import plot_with_axhline, show_plots
"""
TISAT 1
1 36799U 10035E   20215.46271709  .00000102  00000-0  16380-4 0  9994
2 36799  98.1977  27.5125 0012384 172.8977 187.2418 14.91391657 545784
"""

if __name__ == "__main__":
    # OE = TypeOE(eccentricity=0.0012384, inclination=98.1977, right_ascension=27.5125, argument_of_perigee=172.8977, mean_anomaly=187.2418, mean_motion=14.91391657)
    OE = TypeOE(semimajor_axis=6876.6644, eccentricity=0.0020122, inclination=1.6971285,
                right_ascension=2.7791209, argument_of_perigee=1.6541096, mean_anomaly=0.2653564)
    time = datetime(year=2020, month=4, day=15, hour=18, minute=22, second=4)
    observer = TypeLatLong(0, 0)

    # r2_expected = np.array([-6392.0838, 2491.2094, 412.8915])
    # v2_expected = np.array([0.7804, 0.7196, 7.5057])
    incl = 1.6966172451134904
    right_ascension = 2.7772465497672525
    argument_of_perigee = 3.007028209392104

    r2_expected = np.array([-6384.8773242, 2497.81798463, 463.89455914])
    v2_expected = np.array([0.83321341, 0.6980882, 7.50396268])

    # Preparation for calculate from different observers
    linspace = np.linspace(0, 10, num=161, endpoint=True)

    r2s = np.zeros(len(linspace) * 3)
    r2s.resize((len(r2s) // 3, 3))

    v2s = np.zeros(len(linspace) * 3)
    v2s.resize((len(v2s) // 3, 3))

    # diff_r = np.zeros(len(linspace))
    # diff_v = np.zeros(len(linspace))
    diff_incls = np.zeros(len(linspace))
    diff_right_ascensions = np.zeros(len(linspace))
    diff_argument_of_perigees = np.zeros(len(linspace))
    for i in range(0, len(linspace)):
        sim = Simulation(OE, observer, time=time)
        orbit_calculate = OrbitCalculate(sim)
        alias = np.array([random.random(), random.random(), random.random()])
        alias = normalize(alias)
        orbit_calculate.preCalculate([3900, 4000, 4100], observer_alias=linspace[i] * alias)
        r2s[i], v2s[i] = orbit_calculate.calculate([3900, 4000, 4100], shouldPrecalculate=False)
        # diff_r[i] = getDiffAngle(r2s[i], r2_expected, False)
        # diff_v[i] = getDiffAngle(v2s[i], v2_expected)
        # r2s[i] -= r2_expected
        # v2s[i] -= v2_expected
        temp_OE = calc_oe_from_sv(r2s[i], v2s[i])
        diff_incls[i] = temp_OE.inclination - incl
        diff_right_ascensions[i] = temp_OE.right_ascension - right_ascension
        diff_argument_of_perigees[i] = temp_OE.argument_of_perigee - argument_of_perigee

    # print(diff_r[10])
    # print(diff_v[10])
    # print(r2s[10])
    # plt.plot(linspace, diff_incls)
    # plt.show()
    title = "Sai số random theo cả ba trục"
    plt.xlim(left=0,right=linspace[-1])
    # plt.ylim(bottom=0, top=diff_argument_of_perigees[-1])
    # plt.ylim(bottom=0, top=diff_incls[-1] * 2)
    plot_with_axhline(linspace, diff_incls, color='red')
    # plot_with_axhline(linspace, diff_right_ascensions, color='green')
    plot_with_axhline(linspace, diff_right_ascensions, color='green', title=title, legends=['Góc nghiêng', 'Góc điểm lên'])
    # plot_with_axhline(linspace, diff_argument_of_perigees, title="Difference when change X axis", color='blue', legends=['\u0394 inclination', '\u0394 right_ascension', '\u0394 argument_of_perigee'])
    # plot_with_axhline(linspace, diff_argument_of_perigees, title=title, color='blue', legends=['góc cận điểm'])
    plt.xlabel("Đơn vị: km")
    show_plots()
    # print(linspace[1])
    # print(linspace[0])
    # OE = calc_oe_from_sv(r2s[0], v2s[0])
    # print(r2s[0])
    # print(v2s[0])
    # print(OE)
