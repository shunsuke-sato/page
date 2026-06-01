# /// script
# dependencies = [
#   "pandas",
#   "matplotlib",
#   "numpy"
# ]
# ///

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    filename = "td.general/total_current"

    df = pd.read_csv(
        filename,
        comment="#",
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 2, 3, 4],
        names=["Iter", "t", "jx", "jy", "jz"],
        engine="python",
    )

    print("Data preview:")
    print(df.head())
    print("\nData shape:", df.shape)

    cell_volume = 263.7445

    tt = df["t"].to_numpy()
    jxt = df["jx"].to_numpy() / cell_volume
    jyt = df["jy"].to_numpy() / cell_volume
    jzt = df["jz"].to_numpy() / cell_volume

    if len(tt) < 2:
        raise ValueError("Need at least two time points to compute dt.")

    dt = tt[1] - tt[0]

    kick_k = 0.005
    kick_dir = np.array([1.0, 0.0, 0.0])

    ev = 27.2114
    nomega = 400
    omega_ini = 0.1 / ev
    omega_fin = 20.0 / ev
    gamma = 0.3 / ev

    omega = np.linspace(omega_ini, omega_fin, nomega)

    sigma = np.zeros(nomega, dtype=complex)

    for iw in range(nomega):
        jxw = 0.0 + 0.0j
        jyw = 0.0 + 0.0j
        jzw = 0.0 + 0.0j

        for it in range(len(tt)):
            phase = np.exp((1j * omega[iw] - gamma) * tt[it])
            jxw += phase * jxt[it]
            jyw += phase * jyt[it]
            jzw += phase * jzt[it]

        sigma[iw] = -dt * (
            kick_dir[0] * jxw
            + kick_dir[1] * jyw
            + kick_dir[2] * jzw
        ) / kick_k

    epsilon_w = 1.0 + 4.0j * np.pi * sigma / omega

    print("\nFirst 10 omega values:")
    print(omega[:10])

    print("\nFirst 10 sigma values:")
    print(sigma[:10])

    plt.figure()
    plt.plot(omega * ev, epsilon_w.real, label="Re(epsilon)")
    plt.plot(omega * ev, epsilon_w.imag, label="Im(epsilon)")
    plt.xlim(0.0, 20.0)
    plt.ylim(-20, 20.0)
    plt.xlabel("Omega (eV)")
    plt.ylabel("epsilon")
    plt.legend()
    plt.tight_layout()
    plt.savefig("dielectric_function.png")


if __name__ == "__main__":
    main()    
