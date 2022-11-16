import numpy as np

z = np.array([8, 2, 5, 1, 6])
print(z.dtype)

maxper = 6

ticks = list(set(z))
labels = []
for i in ticks:
    if (i < maxper):
        name = f"{i}T"
        labels.append(name)
    elif (i == maxper):
        name = f"MP"
        labels.append(name)

print(ticks)
print(labels)
