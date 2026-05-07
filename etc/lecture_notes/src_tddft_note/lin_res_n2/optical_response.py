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
    # Input file
    filename = "td.general/multipoles"

    # Read data (ignore lines starting with '#')
    df = pd.read_csv(
        filename,
        comment="#",
        sep=r"\s+",
        header=None,
        engine="python"
    )

    # Assign column names
    df.columns = [
        "Iter",
        "t",
        "Electronic_charge",
        "x",
        "y",
        "z",
    ]

    # Show basic info
    print("Data preview:")
    print(df.head())
    print("\nData shape:", df.shape)

    # Extract columns as NumPy arrays
    tt = df["t"].to_numpy()
    xt = df["x"].to_numpy()
    yt = df["y"].to_numpy()
    zt = df["z"].to_numpy()

    dt = tt[1]-tt[0]

    kick_k = 0.01
    kick_dir = np.zeros(3)
    kick_dir[0] = 0.0
    kick_dir[1] = 0.0
    kick_dir[2] = 1.0

    ev = 27.2114
    nomega = 400
    omega_ini = 0.1 / ev
    omega_fin = 20.0 / ev
    gamma = 0.5 / ev  # damping factor; adjust as needed

    omega = np.linspace(omega_ini, omega_fin, nomega)
    domega = omega[1] - omega[0]

    alpha = np.zeros(nomega, dtype=complex)

    for iw in range(nomega):
        xw = 0.0 + 0.0j
        yw = 0.0 + 0.0j
        zw = 0.0 + 0.0j

        for it in range(len(tt)):
            phase = np.exp((1j * omega[iw] - gamma) * tt[it])
            xw += phase * xt[it]
            yw += phase * yt[it]
            zw += phase * zt[it]

        alpha[iw] = -dt*(
            kick_dir[0] * xw +
            kick_dir[1] * yw +
            kick_dir[2] * zw
        ) / kick_k

    # Example output
    print("\nFirst 10 omega values:")
    print(omega[:10])

    print("\nFirst 10 alpha values:")
    print(alpha[:10])

    # Optional plot
    plt.figure()
    plt.plot(omega * ev, alpha.real, label="Re(alpha)")
    plt.plot(omega * ev, alpha.imag, label="Im(alpha)")
    plt.xlim(0.0,20.0)
    plt.xlabel("Omega (eV)")
    plt.ylabel("Alpha")
    plt.legend()
    plt.tight_layout()
#    plt.show()
    plt.savefig("polarizability.png")


if __name__ == "__main__":
    main()
    
