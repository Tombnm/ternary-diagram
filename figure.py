#!/usr/bin/env python3

import ternary
from pymatgen import MPRester, Composition
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.entries.computed_entries import ComputedEntry
#from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.analysis.phase_diagram import *
import matplotlib.pyplot as plt
from ternary_diagram import TernaryDiagram
import pandas as panda

MAPI_KEY = None  # You must change this to your Materials API key! (or set MAPI_KEY env variable)
system = ["H", "Y", "Ba"]  # system we want to get PD for

# Create phase diagram!
processed_entries = []

struct_sta,struct_unsta,struct_metasta = [],[],[]
element_ter,tmp = [],[]
#element_fig = []
element_cor = {}
element_tercor = {}

with open("int_aga.dat", "r") as fin:
    lines = fin.readlines()

    tmp.append(lines[0].split()[1])
    tmp.append(lines[0].split()[2])
    tmp.append(lines[0].split()[3])
    #print(tmp)
    element_ter.append(tmp[0].split('_')[1])
    element_ter.append(tmp[1].split('_')[1])
    element_ter.append(tmp[2].split('_')[1])
    tmp.clear()
    #print("element_ter = ",element_ter)

    for i in range(1, len(lines)):
        com = lines[i].split()[0]

        element_cor[com] = []
        element_cor[com].append(int(lines[i].split()[1]))
        element_cor[com].append(int(lines[i].split()[2]))
        element_cor[com].append(int(lines[i].split()[3]))

        ene_per = float(lines[i].split()[4])
        MY_COMPOSITION = Composition(com)  # put the formula of your prediction
        MY_ENERGY_PER_ATOM = ene_per  # note: this must be energy normalized PER ATOM in eV!
        MY_PARAMETERS = {"potcar_symbols": ['pbe H', 'pbe Y', 'pbe Ba']}

        my_entry = ComputedEntry(MY_COMPOSITION, MY_ENERGY_PER_ATOM * MY_COMPOSITION.num_atoms, parameters=MY_PARAMETERS)
        #my_entry = ComputedEntry(MY_COMPOSITION, MY_ENERGY_PER_ATOM * MY_COMPOSITION.num_atoms)
        processed_entries.append(my_entry)
#print("element_cor = ",element_cor)
#print(len(processed_entries))
pd = PhaseDiagram(processed_entries)
# Plot!
plotter = PDPlotter(pd, show_unstable= True) # you can also try show_unstable=True
#plotter.show()
#plotter.write_image("{}-aga.png".format('-'.join(system)), "png")  # save figure

print('Stable Entries\n--------')
for e in pd.stable_entries:
    struct_sta.append(e.composition.reduced_formula)
    print(e.composition.reduced_formula)
#print("struct_sta = ",struct_sta)

print('\nUnstable Entries\n--------')
for e in pd.unstable_entries:
    decomp, e_above_hull = pd.get_decomp_and_e_above_hull(e)
    pretty_decomp = [("{}:{}".format(k.composition.reduced_formula, k.entry_id), round(v, 2)) for k, v in decomp.items()]
    if float("%.3f" % e_above_hull) <= 0.05:
        struct_metasta.append(e.composition.reduced_formula)
    elif float("%.3f" % e_above_hull) > 0.05:
        struct_unsta.append(e.composition.reduced_formula)
    print(e.composition.reduced_formula, "%.3f" % e_above_hull)
#print("struct_unsta = ",struct_unsta)
#print("struct_unstaener = ",struct_unstaener)

#-------------------------------------------



## Boundary and Gridlines
scale = 1
figure, tax = ternary.figure(scale=scale)

# Draw Boundary and Gridlines
tax.boundary(linewidth=2)
#tax.gridlines(color="black", multiple=1)
tax.gridlines(color="blue", multiple=0.2, linewidth=0.5)

# Set Axis labels and Title
fontsize = 20
tax.set_title("", fontsize=fontsize)
tax.left_axis_label(element_ter[0]+" fraction", fontsize=fontsize)
tax.right_axis_label(element_ter[1]+" fraction", fontsize=fontsize)
tax.bottom_axis_label(element_ter[2]+" fraction", fontsize=fontsize)
tax.right_corner_label(element_ter[0],fontsize=fontsize)
tax.left_corner_label(element_ter[1],fontsize=fontsize)
tax.top_corner_label(element_ter[2],fontsize=fontsize)


#matplotlib operation
#print("***************************")
#print(dir(tax))
ax = tax.get_axes()
ax.axis('off')

tax.ticks(axis='lbr', multiple=0.2, linewidth=1, offset=0.015,tick_formats="%.1f")
tax.clear_matplotlib_ticks()

#points
#tax.scatter([[0.1,0.1,0.8],[0.125,0.125,0.75]], marker='D', color='green', label="Green Diamonds")
#stable struct points
sum = 0
tmpsta = []
for i in range(0,len(struct_sta)):
    sum = sum + element_cor[struct_sta[i]][0] + element_cor[struct_sta[i]][1] + element_cor[struct_sta[i]][2]
    s = panda.Series(element_cor[struct_sta[i]])
    tmpsta.append(s/sum)
    sum = 0
tax.scatter(tmpsta,s=50,c='g',marker='s',label="stable")

#unstable struct points
tmpunsta = []
for i in range(0,len(struct_unsta)):
    sum = sum + element_cor[struct_unsta[i]][0] + element_cor[struct_unsta[i]][1] + element_cor[struct_unsta[i]][2]
    s = panda.Series(element_cor[struct_unsta[i]])
    tmpunsta.append(s/sum)
    sum = 0
tax.scatter(tmpunsta,s=30,c='r',label="unstable")

#metastable struct points
tmpmetasta = []
for i in range(0,len(struct_metasta)):
    sum = sum + element_cor[struct_metasta[i]][0] + element_cor[struct_metasta[i]][1] + element_cor[struct_metasta[i]][2]
    s = panda.Series(element_cor[struct_metasta[i]])
    tmpmetasta.append(s/sum)
    sum = 0
#print("tmpmetasta",tmpmetasta)
tax.scatter(tmpmetasta,s=30,c='b',label="metastable")

#legend()
tax.legend()


# Set ticks
tax.ticks(axis='lbr', linewidth=1)

# Remove default Matplotlib Axes
tax.clear_matplotlib_ticks()

ternary.plt.show()
