import random
import os
import platform

# Generate two random numbers up to 1e4
num1 = random.randint(0, 10000)
num2 = random.randint(0, 10000)

# Add them
sum1 = num1 + num2

# Generate a third random number and add it to the sum
num3 = random.randint(0, 10000)
total_sum = sum1 + num3

# Print the progress
print(f"First number: {num1}")
print(f"Second number: {num2}")
print(f"Sum of first two numbers: {sum1}")
print(f"Third number: {num3}")
print(f"Total sum: {total_sum}")

# Get system information
system_info = {
    "OS": platform.system(),
    "OS Version": platform.version(),
    "Machine": platform.machine(),
    "Processor": platform.processor(),
    "CPU Count": os.cpu_count(),
    "Memory": os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024. ** 3)
}

# Print system information
print("\nSystem Information:")
for key, value in system_info.items():
    print(f"{key}: {value}")

