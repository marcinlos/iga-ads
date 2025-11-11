import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os

# Parametry
size = 20
timesteps = [0, 25, 50, 75, 100]   # momenty w czasie
os.makedirs("wind_maps", exist_ok=True)

def wind_field(x, y, t):
    """Model pola wiatru zmienny w czasie."""
    # bazowy kierunek SW -> NE
    bx_base = 30 + 0.5 * t
    by_base = 60 + 0.7 * t

    # lokalne zaburzenia
    noise_x = 10 * np.sin(0.1 * x + 0.05 * t) * np.cos(0.1 * y)
    noise_y = 10 * np.cos(0.1 * y + 0.03 * t) * np.sin(0.1 * x)

    # wczesna faza: burzowe porywy
    if t < 30:
        bx = bx_base + noise_x * 1.5
        by = by_base + noise_y * 1.5
    # środkowa: stabilizacja
    elif t < 70:
        bx = bx_base + noise_x * 0.8
        by = by_base + noise_y * 0.6
    # późna: silny, ustalony wiatr
    else:
        bx = bx_base + noise_x * 0.3
        by = by_base + noise_y * 0.3
    return bx, by

# Siatka
x = np.linspace(0, 100, size)
y = np.linspace(0, 100, size)
X, Y = np.meshgrid(x, y)

for t in timesteps:
    Bx, By = wind_field(X, Y, t)
    magnitude = np.sqrt(Bx**2 + By**2)

    # Skala szarości do BMP
    norm = ((magnitude - magnitude.min()) / (magnitude.max() - magnitude.min()) * 255).astype(np.uint8)
    Image.fromarray(norm, mode='L').save(f"wind_maps/wildfire_wind_t{t:03}.bmp")

    # Wizualizacja
    plt.figure(figsize=(6, 5))
    plt.quiver(X, Y, Bx, By, magnitude, cmap="coolwarm", scale=2000)
    plt.colorbar(label="Wind magnitude")
    plt.title(f"Wind Field at moment t={t}\nmean(Bx)={Bx.mean():.2f}, mean(By)={By.mean():.2f}, mean(|W|)={magnitude.mean():.2f}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.tight_layout()
    plt.savefig(f"wind_maps/wind_field_t{t:03}.png", dpi=200)
    plt.close()

print("✅ Wygenerowano mapy wiatru i wizualizacje w folderze 'wind_maps/'")
