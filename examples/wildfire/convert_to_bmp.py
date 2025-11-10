# note: just use `python3 convert_to_bmp.py`. default input: map.png, default output: wildfire_fuel.bmp saved both to current dir and 2 levels above

from PIL import Image
import numpy as np
import os
from PIL import ImageEnhance

input_path = "map.png"
output_name = "wildfire_fuel.bmp"

output_local = os.path.join(".", output_name)
output_upper = os.path.join("..", "..", output_name)

img = Image.open(input_path).convert("RGB")
r, g, b = [np.array(c, dtype=np.float32) / 255.0 for c in img.split()]

green_index = g / (0.3 + 0.7 * (r + b) / 2)
green_index = np.clip(green_index, 0, 1)

enhanced = green_index ** 0.5

enhanced = np.power(enhanced, 2.2)

enhanced = (enhanced - enhanced.min()) / (enhanced.max() - enhanced.min() + 1e-8)
enhanced = np.clip((enhanced - 0.2) * 2.5, 0, 1)
gamma = 1.5  
enhanced = np.power(enhanced, gamma)

gray = (enhanced * 255).astype(np.uint8)
gray_img = Image.fromarray(gray)

enhancer = ImageEnhance.Contrast(gray_img)
gray_img = enhancer.enhance(3.0)

gray_resized = gray_img.resize((100, 100), Image.BILINEAR)

gray_resized.save(output_local)
gray_resized.save(output_upper)

print(f"Fuel map saved:\n - {output_local}\n - {output_upper}")
