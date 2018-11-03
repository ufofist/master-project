
import glob

from IPython.display import Image
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import pandas as pd

import openmc
get_ipython().run_line_magic('matplotlib', 'inline')


#uranium oxide
uo2 = openmc.Material(1, "uo2")
uo2.add_nuclide('U235', 0.06,'wo')
uo2.add_nuclide('U238', 0.803,'wo')
uo2.add_nuclide('O16', 0.12,'wo')
uo2.add_nuclide('Pu238', 0.001781,'wo')
uo2.add_nuclide('Pu239', 0.08357,'wo')
uo2.add_nuclide('Pu240', 0.033428,'wo')
uo2.add_nuclide('Pu241', 0.011371,'wo')
uo2.add_nuclide('Pu242', 0.00685,'wo')
uo2.set_density('g/cm3', 11.0)
#print(uo2)

#stainless steel for clad
stainless = openmc.Material(2, "stainless")
stainless.add_element('Fe',0.74,'wo')
stainless.add_element('Cr',0.18,'wo')
stainless.add_element('Ni',0.08,'wo')
stainless.set_density('g/cm3',7.92)
#print(stainless)

#gap gas hellium
He = openmc.Material(3,"He")
He.add_nuclide('He4',1.0)
He.set_density('g/cm3',0.145)
#print(He)

#water
water = openmc.Material(4,"water")
water.add_nuclide('H1',2.0)
water.add_nuclide('O16',1.0)
water.set_density('g/cm3',1.0)
#H2o alpha and beta cs
water.add_s_alpha_beta('c_H_in_H2O')
#print(water)

#Na coolant
na = openmc.Material(5,"sodium")
na.add_element('Na',0.228,'wo')
na.add_element('K',0.778,'wo')
na.set_density('g/cm3',0.88508)    #0.968
#print(na)

#lead-Bi coolant
lead = openmc.Material(6,"leed")
lead.add_element('Pb',0.445,'wo')
lead.add_nuclide('Bi209',0.555,'wo')
lead.set_density('g/cm3',10.204)
#print(lead)

#control
boron = openmc.Material(7,'boron')
boron.add_nuclide('B10',0.3915,'wo')
boron.add_nuclide('B11',0.3915,'wo')
boron.add_nuclide('C0',0.217,'wo')
boron.set_density('g/cm3',2.52)
#print(boron)

#uo2_2
uo2_2 = openmc.Material(8, "uo2")
uo2_2.add_nuclide('U235', 0.03,'wo')
uo2_2.add_nuclide('U238', 0.833,'wo')
uo2_2.add_nuclide('O16', 0.12,'wo')
uo2_2.add_nuclide('Pu238', 0.001781,'wo')
uo2_2.add_nuclide('Pu239', 0.08357,'wo')
uo2_2.add_nuclide('Pu240', 0.033428,'wo')
uo2_2.add_nuclide('Pu241', 0.011371,'wo')
uo2_2.add_nuclide('Pu242', 0.00685,'wo')
uo2_2.set_density('g/cm3', 11.0)
#print(uo2_2)

#uranium
uranium = openmc.Material(9,'uranium')
uranium.add_nuclide('U238',1)
uranium.set_density('g/cm3',19.05)
#print(u238)

#Be
be=openmc.Material(10,'be')
be.add_element('Be',1,'ao')
be.set_density('g/cm3',1.85)
#print(be)

#BeO
beo= openmc.Material(11,'beo')
beo.add_element('Be',1,'ao')
beo.add_element('O',1,'ao')
beo.set_density('g/cm3',3.025)
#print(beo)

#graphite
gc= openmc.Material(12,'gc')
gc.add_element('C',1,'ao')
gc.set_density('g/cm3',2.25)
#print(gc)

#316 s.s 
ss= openmc.Material(13,'wo')
ss.add_element('C',0.0008,'wo')
ss.add_element('Si',0.01,'wo')
ss.add_element('Mn',0.02,'wo')
ss.add_element('Cr',0.185,'wo')
ss.add_element('Ni',0.14,'wo')
ss.add_element('Mo',0.03,'wo')
ss.add_element('Fe',0.6142,'wo')
ss.set_density('g/cm3',7.92)
#print(ss)

#uo2_3
uo2_3 = openmc.Material(14, "uo2")
uo2_3.add_nuclide('U235', 0.01,'wo')
uo2_3.add_nuclide('U238', 0.853,'wo')
uo2_3.add_nuclide('O16', 0.12,'wo')
uo2_3.add_nuclide('Pu238', 0.001781,'wo')
uo2_3.add_nuclide('Pu239', 0.08357,'wo')
uo2_3.add_nuclide('Pu240', 0.033428,'wo')
uo2_3.add_nuclide('Pu241', 0.011371,'wo')
uo2_3.add_nuclide('Pu242', 0.00685,'wo')
uo2_3.set_density('g/cm3', 11.0)
#print(uo2_2)

#uo2_4
uo2_4 = openmc.Material(15, "uo2")
uo2_4.add_nuclide('U235', 0.007,'wo')
uo2_4.add_nuclide('Pu238', 0.001781,'wo')
uo2_4.add_nuclide('Pu239', 0.08357,'wo')
uo2_4.add_nuclide('Pu240', 0.033428,'wo')
uo2_4.add_nuclide('Pu241', 0.011371,'wo')
uo2_4.add_nuclide('Pu242', 0.00685,'wo')
uo2_4.add_nuclide('U238', 0.856,'wo')
uo2_4.add_nuclide('O16', 0.12,'wo')
uo2_4.set_density('g/cm3', 11.0)
#print(uo2_2)

#k
k= openmc.Material(16)
k.add_element('K',1.0)
k.set_density('g/cm3',0.86)

#na and k 1
nak1=openmc.Material(17)
nak1.add_element('Na',0.44)
nak1.add_element('K',0.56)
nak1.set_density('g/cm3',0.9084)

#na and k 2
nak2=openmc.Material(18)
nak2.add_element('Na',0.228)
nak2.add_element('K',0.772)
nak2.set_density('g/cm3',0.88508)

#create materials file
mats = openmc.Materials([uo2,stainless,He,water,na,lead,boron,uo2_2,uo2_3,uo2_4,uranium,be,beo,gc,ss,k,nak1,nak2])
mats.export_to_xml()
# !cat materials.xml

#define geometry
#plane
fuel_or = openmc.ZCylinder(R=0.54)
clad_ir = openmc.ZCylinder(R=0.56)
clad_or = openmc.ZCylinder(R=0.60)
z1=openmc.ZPlane(z0=0)
z2=openmc.ZPlane(z0=60)
z3=openmc.ZPlane(z0=150)
z4=openmc.ZPlane(z0=165)
z5=openmc.ZPlane(z0=-0.4)
z6=openmc.ZPlane(z0=165.4)
z7=openmc.ZPlane(z0=82.5)
z8=openmc.ZPlane(z0=123.95)
z9=openmc.ZPlane(z0=41.05)
top = openmc.ZPlane(z0=200, boundary_type='vacuum')
bottom = openmc.ZPlane(z0=-200, boundary_type='vacuum')
box = openmc.get_rectangular_prism(width=400, height=400,boundary_type='vacuum')
type(box)

#region
fuel_region= -fuel_or& +z2& -z3
gap_region= +fuel_or& -clad_ir& +z2& -z3
clad_region = +clad_ir& -clad_or& +z1& -z4|-clad_or& -z1&+z5|-clad_or&-z6&+z4
gas_region = -clad_ir& +z1& -z2|-clad_ir& -z4& +z3
pin_region = -clad_or& +z5& -z6
coolant_region = ~pin_region & -top& +bottom &box
control_region = -clad_or& -z6& +z5
steel_region = -z6& +z5
controlvoid_region= ~control_region

#cell
fuel= openmc.Cell(1,'fuel')
fuel.fill= uo2
fuel.region=fuel_region

fuel_2= openmc.Cell(15,'fuel2')
fuel_2.fill= uo2_2
fuel_2.region=fuel_region

fuel_3= openmc.Cell(16,'fuel3')
fuel_3.fill= uo2_3
fuel_3.region=fuel_region

fuel_4= openmc.Cell(17,'fuel4')
fuel_4.fill= uo2_4
fuel_4.region=fuel_region

gap=openmc.Cell(2, 'He gap')
gap.fill= He
gap.region=gap_region

gap_2=openmc.Cell(18, 'He gap2')
gap_2.fill= He
gap_2.region=gap_region

gap_3=openmc.Cell(19, 'He gap3')
gap_3.fill= He
gap_3.region=gap_region

gap_4=openmc.Cell(20, 'He gap4')
gap_4.fill= He
gap_4.region=gap_region

clad=openmc.Cell(3, 'clad')
clad.fill= stainless
clad.region=clad_region

clad_2=openmc.Cell(21, 'clad2')
clad_2.fill= stainless
clad_2.region=clad_region

clad_3=openmc.Cell(22, 'clad3')
clad_3.fill= stainless
clad_3.region=clad_region

clad_4=openmc.Cell(23, 'clad4')
clad_4.fill= stainless
clad_4.region=clad_region

gas=openmc.Cell(4, 'gas plenum')
gas.fill= He
gas.region=gas_region

gas_2=openmc.Cell(24, 'gas plenum2')
gas_2.fill= He
gas_2.region=gas_region

gas_3=openmc.Cell(25, 'gas plenum3')
gas_3.fill= He
gas_3.region=gas_region

gas_4=openmc.Cell(26, 'gas plenum4')
gas_4.fill= He
gas_4.region=gas_region

coolant= openmc.Cell(5, 'coolant')
coolant.fill= na
coolant.region = coolant_region

coolant_2= openmc.Cell(27, 'coolant2')
coolant_2.fill= na
coolant_2.region = coolant_region

coolant_3= openmc.Cell(28, 'coolant3')
coolant_3.fill= na
coolant_3.region = coolant_region

coolant_4= openmc.Cell(29, 'coolant4')
coolant_4.fill= na
coolant_4.region = coolant_region

coolant2 = openmc.Cell(6, 'coolant2')
coolant2.fill= na

steel = openmc.Cell(7,'steel')
steel.fill = stainless

steel_2 = openmc.Cell(100,'steel2')
steel_2.fill = stainless

steel_3 = openmc.Cell(31,'steel3')
steel_3.fill = stainless

steel_4 = openmc.Cell(32,'steel4')
steel_4.fill = stainless

control = openmc.Cell(8,'control')
control.fill = na      #boron
control.region = control_region

control_2 = openmc.Cell(51,'control2')
control_2.fill = na      #boron
control_2.region = control_region

control_3 = openmc.Cell(52,'control3')
control_3.fill = na      #boron
control_3.region = control_region

control_4 = openmc.Cell(53,'control4')
control_4.fill = na      #boron
control_4.region = control_region

controlvoid = openmc.Cell(9,'controlvoid')
controlvoid.fill= na
controlvoid.region=controlvoid_region

controlvoid_2 = openmc.Cell(36,'controlvoid2')
controlvoid_2.fill= na
controlvoid_2.region=controlvoid_region

controlvoid_3 = openmc.Cell(37,'controlvoid3')
controlvoid_3.fill= na
controlvoid_3.region=controlvoid_region

controlvoid_4 = openmc.Cell(38,'controlvoid4')
controlvoid_4.fill= na
controlvoid_4.region=controlvoid_region

breedfuel = openmc.Cell(10,'breedfuel')
breedfuel.fill= uranium
breedfuel.region=fuel_region

breedcoolant = openmc.Cell(11,'breedcoolant')
breedcoolant.fill= na
breedcoolant.region=coolant_region

breedgas = openmc.Cell(12,'breedgas')
breedgas.fill= He
breedgas.region=gas_region

breedgap = openmc.Cell(13,'breedgap')
breedgap.fill= He
breedgap.region=gap_region

breedclad = openmc.Cell(14,'breedclad')
breedclad.fill= stainless
breedclad.region=clad_region

#assembly
#assemebly cell
f=openmc.Universe(name='Fuel Pin', cells=[fuel,gap,clad,gas,coolant])
all_coolant= openmc.Universe(cells=[coolant2])
s=openmc.Universe(cells=[steel])
c=openmc.Universe(name='control rod',cells=[control,controlvoid])

#lattice
assembly = openmc.HexLattice()
assembly.center = (0, 0)
assembly.pitch = [1.36]
assembly.outer=all_coolant
assembly.universes = [[f]*48,
                      [f]*42,
                      [f]*36,
                      [f]*30,
                      [f]*24,
                      [f]*18,
                      [s]*12,
                      [c]*6,
                      [c]]
                      
# Create boundary planes to surround the geometry
asseouter_sphere= openmc.Sphere(R=300,boundary_type='vacuum')
asseouter_surface = openmc.ZCylinder(R=100, boundary_type='vacuum')
min_z = openmc.ZPlane(z0=-1.4, boundary_type='vacuum')
max_z = openmc.ZPlane(z0=+208.6, boundary_type='vacuum')
max_z2= openmc.ZPlane(z0=+165.4, boundary_type='vacuum')
min_z2= openmc.ZPlane(z0=-0.5,boundary_type='vacuum')

# cell for assembly
asse = openmc.Cell(77,name='asse', fill=assembly)
c1=openmc.Cell()
c2=openmc.Cell()
# Add boundary planes
c1.fill=all_coolant
c1.region= +max_z2 & -asseouter_sphere
c2.fill=all_coolant
c2.region= -min_z2 & -asseouter_sphere
asse.region= -asseouter_sphere& +min_z2& -max_z2
# Create assembly Universe
root = openmc.Universe(name='root assembly')
root.add_cell(asse)
root.add_cell(c1)
root.add_cell(c2)
rotated_cell = openmc.Cell(fill=assembly)
rotated_universe = openmc.Universe(cells=[rotated_cell])
asse.fill= rotated_universe
asse.rotation = (0.,0.,30.)
#root.plot(width=(20,20),origin=(0,0,0),pixels=(500,500),color_by='material')#,basis= 'yz')

f2=openmc.Universe(name='Fuel Pin2', cells=[fuel_2,gap_2,clad_2,gas_2,coolant_2])
s2=openmc.Universe(cells=[steel_2])
c2=openmc.Universe(name='control rod2',cells=[control_2,controlvoid_2])
assembly_2 = openmc.HexLattice()
assembly_2.center = (0, 0)
assembly_2.pitch = [1.36]
assembly_2.outer=all_coolant
assembly_2.universes = [[f2]*48,
                        [f2]*42,
                        [f2]*36,
                        [f2]*30,
                        [f2]*24,
                        [f2]*18,
                        [s2]*12,
                        [c2]*6,
                        [c2]]
asse2 = openmc.Cell(name='asse2', fill=assembly_2)
c7=openmc.Cell()
c8=openmc.Cell()
c7.fill=all_coolant
c7.region= +max_z2 & -asseouter_sphere
c8.fill=all_coolant
c8.region= -min_z2 & -asseouter_sphere
asse2.region= -asseouter_sphere& +min_z2& -max_z2

root_2 = openmc.Universe(name='root assembly2')
root_2.add_cell(asse2)
root_2.add_cell(c7)
root_2.add_cell(c8)
rotated_cell2 = openmc.Cell(fill=assembly_2)
rotated_universe2 = openmc.Universe(cells=[rotated_cell2])
asse2.fill= rotated_universe2
asse2.rotation = (0.,0.,30.)
#root_2.plot(width=(20,20),origin=(0,0,0),pixels=(500,500),color_by='material')#,basis= 'yz')

f3=openmc.Universe(name='Fuel Pin3', cells=[fuel_3,gap_3,clad_3,gas_3,coolant_3])
s3=openmc.Universe(cells=[steel_3])
c3=openmc.Universe(name='control rod3',cells=[control_3,controlvoid_3])
assembly_3 = openmc.HexLattice()
assembly_3.center = (0, 0)
assembly_3.pitch = [1.36]
assembly_3.outer=all_coolant
assembly_3.universes = [[f3]*48,
                        [f3]*42,
                        [f3]*36,
                        [f3]*30,
                        [f3]*24,
                        [f3]*18,
                        [s3]*12,
                        [c3]*6,
                        [c3]]
asse3 = openmc.Cell(name='asse3', fill=assembly_3)
c9=openmc.Cell()
c10=openmc.Cell()
c9.fill=all_coolant
c9.region= +max_z2 & -asseouter_sphere
c10.fill=all_coolant
c10.region= -min_z2 & -asseouter_sphere
asse3.region= -asseouter_sphere& +min_z2& -max_z2

root_3 = openmc.Universe(name='root assembly3')
root_3.add_cell(asse3)
root_3.add_cell(c9)
root_3.add_cell(c10)
rotated_cell3 = openmc.Cell(fill=assembly_3)
rotated_universe3 = openmc.Universe(cells=[rotated_cell3])
asse3.fill= rotated_universe3
asse3.rotation = (0.,0.,30.)
#root_3.plot(width=(20,20),origin=(0,0,0),pixels=(500,500),color_by='material')#,basis= 'yz')

f4=openmc.Universe(name='Fuel Pin4', cells=[fuel_4,gap_4,clad_4,gas_4,coolant_4])
s4=openmc.Universe(cells=[steel_4])
c4=openmc.Universe(name='control rod4',cells=[control_4,controlvoid_4])
assembly_4 = openmc.HexLattice()
assembly_4.center = (0, 0)
assembly_4.pitch = [1.36]
assembly_4.outer=all_coolant
assembly_4.universes = [[f4]*48,
                        [f4]*42,
                        [f4]*36,
                        [f4]*30,
                        [f4]*24,
                        [f4]*18,
                        [s4]*12,
                        [c4]*6,
                        [c4]]
asse4 = openmc.Cell(name='asse4', fill=assembly_4)
c11=openmc.Cell()
c12=openmc.Cell()
c11.fill=all_coolant
c11.region= +max_z2 & -asseouter_sphere
c12.fill=all_coolant
c12.region= -min_z2 & -asseouter_sphere
asse4.region= -asseouter_sphere& +min_z2& -max_z2

root_4 = openmc.Universe(name='root assembly4')
root_4.add_cell(asse4)
root_4.add_cell(c11)
root_4.add_cell(c12)
rotated_cell4 = openmc.Cell(fill=assembly_4)
rotated_universe4 = openmc.Universe(cells=[rotated_cell4])
asse4.fill= rotated_universe4
asse4.rotation = (0.,0.,30.)
#root_4.plot(width=(20,20),origin=(0,0,0),pixels=(500,500),color_by='material')#,basis= 'yz')

#breed 
r=openmc.Universe(name='bredd Pin', cells=[breedfuel,breedgap,breedclad,breedgas,breedcoolant])
#r.plot(width=(1.3, 1.3),origin=(0,0,105),pixels=(500,500))

breed = openmc.HexLattice()
breed.center = (0, 0)
breed.pitch = [1.36]
breed.outer= all_coolant
breed.universes = [[r]*48,
                   [r]*42,
                   [r]*36,
                   [r]*30,
                   [r]*24,
                   [r]*18,
                   [r]*12,
                   [r]*6,
                   [r]]
# cell for breed
bree = openmc.Cell(name='bree', fill=breed)

# Add boundary planes
bree.region= -asseouter_sphere& +min_z2& -max_z2
# Create assembly Universe
root_breed = openmc.Universe(name='root bree')
root_breed.add_cell(bree)
bree_cell = openmc.Cell(fill=breed)
bree_universe = openmc.Universe(cells=[bree_cell])
bree.fill= bree_universe
bree.rotation = (0.,0.,30.)
#root_breed.plot(width=(10,10),origin=(0,0,0),pixels=(500,500),color_by='material')#,basis= 'yz')

a1 = openmc.Universe(name='assembly unit1')
ba = openmc.Universe(name='breed unit')
a2 = openmc.Universe(name='assembly unit2')
a3 = openmc.Universe(name='assembly unit3')
a4 = openmc.Universe(name='assembly unit4')
a1.add_cell(asse)
ba.add_cell(bree)
a2.add_cell(asse2)
a3.add_cell(asse3)
a4.add_cell(asse4)

#lattice
core = openmc.HexLattice()
core.center = (0, 0)
core.pitch = [20.4]
core.outer= all_coolant
core.universes = [[ba]*30,
                  [a4]*24,
                  [a3]*18,
                  [a2]*12,
                  [a1]*6,
                  [a1]]
# Create boundary planes to surround the geometry
coreout_sphere=openmc.Sphere(R=800,boundary_type='vacuum')
reflector_outer= openmc.ZCylinder(R=135,boundary_type='vacuum')
reflector_inner= openmc.ZCylinder(R=115,boundary_type='vacuum')
min_zz = openmc.ZPlane(z0=-117.55, boundary_type='vacuum')   #-117.55
max_zz = openmc.ZPlane(z0=+282.45, boundary_type='vacuum')   #+282.45
min_z3 = openmc.ZPlane(z0=-67.55, boundary_type='vacuum')    #-67.55
max_z3 = openmc.ZPlane(z0=+232.45, boundary_type='vacuum')   #+232.45
#max_z2= 165.4
#min_z2= -0.5

# Create core Cell
corecell = openmc.Cell(name='core cell', fill=core)
c3=openmc.Cell()
c4=openmc.Cell()
c5=openmc.Cell()
c6=openmc.Cell()
c7=openmc.Cell()
c8=openmc.Cell()
# Add boundary planes
c3.fill=be
c3.region= +max_z2 & -reflector_inner & -max_zz
c4.fill=be
c4.region= -min_z2 & -reflector_inner & +min_zz
c5.region= +reflector_inner& -reflector_outer & -max_zz & +min_zz 
c5.fill=be
c6.region= -coreout_sphere& +reflector_outer& +max_zz& -min_zz
c6.fill=None
c7.region= +max_z2 & -reflector_inner & -max_z3
c7.fill=na
c8.region= -min_z2 & -reflector_inner & +min_z3
c8.fill=na
corecell.region= -reflector_inner& +min_z2& -max_z2
# Create root Universe
root_universe = openmc.Universe(name='root core')
root_universe.add_cell(corecell)
root_universe.add_cell(c3)
root_universe.add_cell(c4)
root_universe.add_cell(c5)
root_universe.add_cell(c6)
root_universe.add_cell(c7)
root_universe.add_cell(c8)
#root_universe.plot(width=(250,250),origin=(0,0,82),pixels=(500,500),color_by='material')

wholecore = openmc.Geometry(root_universe)
wholecore.export_to_xml()
#!cat geometry.xml

p=openmc.Plot()
p.filename = 'coreplot'
p.width=(220,220)
p.pixels=(3000,3000)
p.origin=(0,0,105)
p.color_by = 'material'
p.colors = {uo2:'yellow',na:'blue',stainless:'black',He:'grey'}

plots= openmc.Plots([p])
plots.export_to_xml()

openmc.plot_geometry()

get_ipython().system('convert coreplot.ppm core1.png')

tallies= openmc.Tallies()
# Create mesh which will be used for tally
mesh_core= openmc.Mesh()
mesh_core.dimension = [21,1,1]
mesh_core.lower_left=(-100,-100,0)
mesh_core.upper_right=(100,100,160)
# Create filter for tally
mesh_filtercore = openmc.MeshFilter(mesh_core)
core_filter=openmc.UniverseFilter([a1,a2,a3,a4,ba])
energy_filter = openmc.EnergyFilter([0.0,1.0,1.0e5,2.0e7])
# Create mesh tally to score flux and fission rate
tally_core = openmc.Tally(name='flux')
tally_core.filters = [mesh_filtercore]
tally_core.scores = ['fission-q-prompt']

tallies.append(tally_core)
tallies.export_to_xml()

coretest = openmc.Settings()
coretest.batches =40                   
coretest.inactive =10                 
coretest.particles = 10000

core_dist = openmc.stats.Box((-10,-10,60),(10,10,70),only_fissionable=True)
coretest.source = openmc.source.Source(space=core_dist)

coretest.export_to_xml()
# !cat coretest.xml

openmc.run()

