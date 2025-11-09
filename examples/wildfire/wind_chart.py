import numpy as np
import matplotlib.pyplot as plt

def compute_wind(y, t):
    speed = 10.0
    if y < 30.0:
        return 0.0, -speed
    angle = (np.pi / 2.0) * (1.0 - np.exp(-t / 1500.0))
    bx = -speed * np.cos(angle)
    by = -speed * np.sin(angle)
    return bx, by

# --- dane ---
times = np.linspace(0, 6000, 100)
y = 50.0
U, V = [], []
for t in times:
    bx, by = compute_wind(y, t)
    U.append(bx)
    V.append(by)
U, V = np.array(U), np.array(V)

# --- przygotowanie wykresu ---
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-11, 1)
ax.set_ylim(-11, 1)
ax.set_xlabel("Składowa pozioma (Zachód–Wschód)")
ax.set_ylabel("Składowa pionowa (Południe–Północ)")
ax.set_title("Ewolucja kierunku wiatru w czasie")

# trajektoria końcówki wektora
ax.plot(U, V, color="lightgray", linewidth=1.5, label="ścieżka końca wektora")

# strzałki co pewien krok czasu
step = 10
for i in range(0, len(times), step):
    ax.arrow(0, 0, U[i], V[i],
             head_width=0.5, head_length=0.7,
             fc="royalblue", ec="royalblue", alpha=0.7)

# punkt początkowy i końcowy
ax.scatter(U[0], V[0], color="green", s=50, label="start (t=0)")
ax.scatter(U[-1], V[-1], color="red", s=50, label="koniec (t=6000)")
ax.legend()

# siatka i proporcje
ax.grid(True, linestyle="--", alpha=0.6)
ax.set_aspect("equal", adjustable="box")

plt.tight_layout()
plt.savefig("wind_chart.png", dpi=200)
print("✅ Zapisano: wind_chart.png")
