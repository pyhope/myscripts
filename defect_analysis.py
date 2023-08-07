
from ase.io import read
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

frame = read("selected.dump", index=0, format="lammps-dump-text")

chemical_symbols = []

for atom in frame:
    if atom.symbol == 'H':
        chemical_symbols.append('Mg')
    elif atom.symbol == 'He':
        chemical_symbols.append('Si')
    elif atom.symbol == 'Li':
        chemical_symbols.append('O')

frame.set_chemical_symbols(chemical_symbols)

x_lst = []
y_lst = []
z_lst = []
for atom in frame.positions:
    x_lst.append(atom[0])
    y_lst.append(atom[1])
    z_lst.append(atom[2])

df = pd.DataFrame({"Index": [i for i in range(1, len(frame) + 1)], "Type": frame.get_chemical_symbols(), "X": x_lst, "Y": y_lst, "Z": z_lst})

df = df.loc[(15.51785 < df.Z) & (df.Z < 95.11695) & ((df.Type == "Mg") | (df.Type == "Si"))]

df.sort_values(by=["Z"], inplace=True, ignore_index=True)
df.to_csv("sorted.csv")

z_np = df.Z.to_numpy()
dz_np_sorted = np.sort(z_np[1:] - z_np[:-1])
ddz_np_sorted = dz_np_sorted[1:] - dz_np_sorted[:-1]

ddz_max_index = ddz_np_sorted.argmax()
dz_thres = (dz_np_sorted[ddz_max_index + 1] + dz_np_sorted[ddz_max_index] * 1.414) / 2.414

layers = []

df_iter = df.iterrows()
current_layer = [next(df_iter)[0]]
for index, row in df_iter:
    if row.Z - df.loc[current_layer[-1], "Z"] > dz_thres:
        layers.append(current_layer)
        current_layer = [index]
    else:
        current_layer.append(index)
if current_layer:
    layers.append(current_layer)

atom_num_per_layer = 40

# quantities to calculate
si_to_mg = 0
mg_to_si = 0
vac_in_si = 0
vac_in_mg = 0
total_si_perfect = 0
total_mg_perfect = 0

print("Total number of layers:", len(layers))
print()
for layer in layers:
    print("Number of atoms in layer:", len(layer))
    si_num = 0
    mg_num = 0
    for id in layer:
        if df.loc[id, "Type"] == "Si":
            si_num += 1
        elif df.loc[id, "Type"] == "Mg":
            mg_num += 1
    if si_num > mg_num:
        print("Si layer")
        total_si_perfect += atom_num_per_layer
        vac_in_si += atom_num_per_layer - (si_num + mg_num)
        si_to_mg += mg_num
    else:
        print("Mg layer")
        total_mg_perfect += atom_num_per_layer
        vac_in_mg += atom_num_per_layer - (si_num + mg_num)
        mg_to_si += si_num
    print("Number of Si atoms:", si_num)
    print("Number of Mg atoms:", mg_num)
    print()

print("Total number of Si atoms in perfect struct:", total_si_perfect)
print("Total number of Mg atoms in perfect struct:", total_mg_perfect)
print("Si to Mg:", si_to_mg)
print("Mg to Si:", mg_to_si)
print("Vacancy in Si:", vac_in_si)
print("Vacancy in Mg:", vac_in_mg)




