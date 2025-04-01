# -*- coding: utf-8 -*-

#USER DEFINED#
bond_length_cutoff=0.05 #.05
bond_angle_cutoff=1   #5

###############

folder = '.' #Must set to your operating folder

import shutil
from shutil import copyfile
#from google.colab import files
import re
import bisect
import os
import glob
import math
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from ase.io import read, write

def calculate_bond_distances(coordinates):
    # Calculate pairwise distances
    distances = pdist(coordinates)
    
    # Convert to square matrix
    distance_matrix = squareform(distances)
    
    return distance_matrix


with open('sys.txt', 'r') as file:
    system = file.read().replace('\n', '')




K=np.loadtxt('k.txt')
#outpath= os.path.join(wd, 'OUT')
#wd= os.getcwd()
#target_out_file='OUT/'+system+'.out'
#shutil.copy(target_out_file, os.path.join(wd, system+'.out'))

#wd=os.getcwd()
#os.chdir(os.path.join(wd,'xyz'))
#for item in os.listdir():
#    os.remove(item)
#os.chdir(wd)


filestest = glob.glob("*")
listOfFiles = filestest
 
listOfOut_temp=[0]*len(listOfFiles)
listOfXYZ_temp=[0]*len(listOfFiles) #creates a blank list to get filled with the filenames of files that worked

for file in listOfFiles:
    if '.out' in listOfFiles[listOfFiles.index(file)]:
        listOfOut_temp[listOfFiles.index(file)]=file   
    if '.xyz' in listOfFiles[listOfFiles.index(file)]:
        listOfXYZ_temp[listOfFiles.index(file)]=file   
listOfOut=list(filter(lambda a: a != 0, listOfOut_temp))  #removes the file name of the csv output from the files that get searched for energies https://www.geeksforgeeks.org/lambda-filter-python-examples/
complete_temp=[0]*len(listOfFiles) #creates a blank list to get filled with the filenames of files that worked
XYZList = listOfXYZ_temp
XYZList=list(filter(lambda a: a != 0, XYZList))
#print('XYZ: ', XYZList)

for entry in listOfOut: 
    file=open(entry) #opens each file in the folder
    if 'int' not in entry:
      for line in file: #reads through the file lines
         if line.startswith('        *  Thank you'): #checks if file is complete
             complete_temp[listOfOut.index(entry)]=entry #replaces a 0 with the file name if it worked
    else:
      complete_temp[listOfOut.index(entry)]=entry
file.close()
complete=complete_temp
complete=list(filter(lambda a: a != 0, complete_temp))  #removes the file name of the csv output from the files that get searched for energies
failed = [i for i in complete + listOfOut if i not in complete or i not in listOfOut] #finds unique values between the files that have thank you and those that dont. From https://www.geeksforgeeks.org/python-difference-two-lists/
#print('files that didnt work')
#print(failed)
print("irc: ", complete)
temp_name = [0]*len(complete)
for entry2 in complete: 
  temp_name[complete.index(entry2)] = entry2[:-4]
  #mkdir(temp_name[complete.index(entry2)])  #use Selin's hint to check if it exists first and warn you about overwriting or add a number to the end

def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect.bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

startmatch = 'Reaction path following'
endmatch = 'Nuclear Repulsion Energy'
coordmatch = 'Standard Nuclear Orientation (Angstroms)'

for entry2 in complete: 
  filename = entry2
  #subfolder =  subfolder = temp_name[complete.index(entry2)]+'/'
  coordlines = []
  ircsteplines = []
  f = open(filename, "r")
  all_lines=f.readlines()

  with open(filename) as ircfile:
      for num, line in enumerate(ircfile, 0):
          if startmatch in line:
              ircsteplines.append(num)
          if coordmatch in line:
              coordlines.append(num)

  actual_coordlines = []
  for ircstep in ircsteplines:
      coordindex = find_lt(coordlines,ircstep)
      actual_coordlines.append(coordindex)

  with open(filename) as ircfile, open('irc_temp.xyz', 'a') as outfile: 
      copy = False
      geomlines = []
      for num, line in enumerate(ircfile, 0):
          if num in actual_coordlines:
              copy= True
          elif endmatch in line: #startswith(endmatch): 
              copy= False
          elif copy:
              if '----------------------------------------------------------------' not in line and 'Standard Nuclear Orientation' not in line: 
                  geomlines.append(line)

          if num in ircsteplines: 
              outfile.write(all_lines[num+1])
              for geomline in geomlines:
                  outfile.write(geomline)
              geomlines=[]


for entry2 in complete: 
    #subfolder = temp_name[complete.index(entry2)]+'/'
    txt = open('irc_temp.xyz').read()
    #print(txt)
    data_list = txt.split('\n') # create lines from the after text
    #print(after_lines)

    # Initialize empty lists and matrices to store data
    step_list = []
    step_atoms = []
    step_matrices = []
    atoms=[]
    atoms_list=[]
    s_list=[]
    a=1
    # Parse the input list
    for line in data_list:
        if line.startswith(' Step'):
            if int(line.split()[1])  in step_list:
                a=-1
            step_list.append(a*int(line.split()[1]))
            s_list.append(a*float(line.split()[-3]))
            if step_matrices:
                step_matrices = np.array(step_matrices)
                step_atoms.append(step_matrices)
                step_matrices = []
                atoms = np.array(atoms)
                atoms_list.append(atoms)
                atoms=[]
        elif line.startswith('    I') or line.strip() == '':
            # Skip the header line and empty lines
            continue
        else:
            atom_info = line.split()
            atoms.append([(atom_info[1])])
            step_matrices.append([float(atom_info[2]), float(atom_info[3]), float(atom_info[4].replace('\\', ''))])
    

    step_matrices = np.array(step_matrices)
    step_atoms.append(step_matrices)
    atoms = np.array(atoms)
    atoms_list.append(atoms)
    
    df = pd.DataFrame({
        'Step': step_list,
        'Atoms':atoms_list,
        'Coordinates': step_atoms,
        's': s_list
    })
wd=os.getcwd()
shutil.copy('irc_temp.xyz', os.path.join(wd,'xyz'))
os.remove('irc_temp.xyz')

#Calculate distances between the R and P for all atoms

for i in range(len(df.Atoms[0])):
    atom1=df['Coordinates'][df['Step'].idxmin()][i]
    atom2=df['Coordinates'][df['Step'].idxmax()][i]
    #print(calculate_3d_distance(atom1,atom2))

IPythonConsole.ipython_3d = True
#if os.path.exists('xyz'):
#    shutil.rmtree('xyz')    

#create xyz files for the end points of the irc and the TS
from R_P_TS_xyz import process_irc_output_files
process_irc_output_files(len(np.array(df.Coordinates)[0]))
files=glob.glob('xyz/*')
for file in files:
    struc=read(file, index=-1)
    print(len(struc))
    os.remove(file)
    write(file, struc)
    print(len(read(file, index=':')))


xyz_files=glob.glob("xyz/*")

raw_mol = Chem.MolFromXYZFile(xyz_files[0])

from rdkit.Chem import rdDetermineBonds
conn_mol = Chem.Mol(raw_mol)
rdDetermineBonds.DetermineConnectivity(conn_mol)#,charge=-2)

def distance(point1, point2):
    return np.linalg.norm(np.array(point1) - np.array(point2))

bond_atom_indices=[]

for i in xyz_files:
    raw_mol = Chem.MolFromXYZFile(i)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    
    for bond in conn_mol.GetBonds():
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        atom1_pos = conn_mol.GetConformer().GetAtomPosition(atom1_idx)
        atom2_pos = conn_mol.GetConformer().GetAtomPosition(atom2_idx)
        bond_length = distance(atom1_pos, atom2_pos)
        #print(f"Bond between atom {atom1_idx} and atom {atom2_idx} - Length: {bond_length:.5f} Å")
        if [atom1_idx,atom2_idx] not in bond_atom_indices:
            bond_atom_indices.append([atom1_idx,atom2_idx])

distances=[]
bond_length_change=[]
for pair in bond_atom_indices:
    dist=[]
    for i in range(len(df)):
        pairwise_distances=calculate_bond_distances(df['Coordinates'][i])
        dist.append(pairwise_distances[pair[0]][pair[1]])
    #if (max(dist)-min(dist))>bond_length_cutoff: 
        #plt.scatter(s_list,dist,label=str(pair[0]+1)+'-'+str(pair[1]+1))
    distances.append(dist)
    bond_length_change.append(max(dist)-min(dist))
'''    
plt.xlabel('s')
plt.ylabel('Bond Length (A)')
plt.title(complete[0].split('_')[0]+' cutoff: '+str(bond_length_cutoff))
plt.legend()

plt.show()
'''

newdf = pd.DataFrame({
    'Atom_indices': bond_atom_indices,
    'Bond_lengths':distances,
    'Bond_length_changes':bond_length_change})

sorted_newdf = newdf.sort_values(by='Bond_length_changes', ascending=False)
#Fit splines to the bond lengths along the reaction path
from scipy.interpolate import CubicSpline

def sort_and_remove_duplicates(x, y):
    # Zip x and y together
    zipped = list(zip(x, y))
    
    # Sort the zipped list based on the values of x
    zipped.sort(key=lambda pair: pair[0])
    
    # Initialize lists to store sorted x and y without duplicates
    sorted_x = []
    sorted_y = []
    
    # Iterate through the sorted zipped list to remove duplicates
    prev_x = None
    for pair in zipped:
        current_x, current_y = pair
        if current_x != prev_x:
            sorted_x.append(current_x)
            sorted_y.append(current_y)
            prev_x = current_x
    
    return sorted_x, sorted_y




splines=[]
for dist in sorted_newdf.Bond_lengths:

    sorted_s, dist = sort_and_remove_duplicates(s_list, dist)
    
    cs = CubicSpline(sorted_s, dist)
    splines.append(CubicSpline(sorted_s, dist))
    # Plot data points
  #  plt.scatter(sorted_s, dist, label='Data')
    
    # Plot spline curve
    s = np.linspace(min(sorted_s), max(sorted_s), 100)
 #   plt.plot(s, cs(s), label='Cubic Spline',color='k')

    
# Add labels and legend
#plt.xlabel('s')
#plt.ylabel('Bond length')
#plt.title('Cubic Spline Fit to X-Y Data')
#plt.legend()

# Show plot
#plt.grid(True)
#plt.show()


def find_crossing_points(spline1, spline2, s_list, num_points=1000, tolerance=1):
    x_values = np.linspace(min(s_list), max(s_list), num_points)
    crossing_points = []
    #plt.plot(x_values,spline1(x_values))
    #plt.plot(x_values,spline2(x_values))
    for x in x_values:
        y1 = spline1(x)
        y2 = spline2(x)
        if abs(y1 - y2)/abs(max(y1,y2))*100 < tolerance:
            crossing_points.append((x, y1))
    closest_x_points = []
    for crossing_point in crossing_points:
        x_crossing = crossing_point[0]
        closest_x_point = min(s_list, key=lambda x: abs(x - x_crossing))
        if closest_x_point!=0:
            closest_x_points.append(closest_x_point)
    index_list = [closest_x_points.index(point) for point in closest_x_points]
    new_index_list=np.unique(index_list)
    #print(new_index_list)
    consecutive_count =0
    my_list=[]
    for ind in new_index_list:
        if consecutive_count ==0:
            my_list.append(ind)
            consecutive_count += 1
        elif ind - prev_index in np.arange(1,50,1):
            consecutive_count += 1
        else: 
            consecutive_count = 0

        prev_index = ind
    my_list2=[]
    for j in my_list:
        if closest_x_points[j] not in [min(s_list),max(s_list),0]:
            my_list2.append(j)
    
    return [closest_x_points[i] for i in my_list2]


def most_separated(numbers, K):
    sorted_numbers = sorted(numbers)
    n = len(sorted_numbers)
    
    # Calculate the segment size
    segment_size = n // (K + 1)
    
    # Initialize list to store selected points
    selected_points = [sorted_numbers[0], sorted_numbers[-1]]  # Select the lowest and highest points
    
    # Select additional points closest to segment boundaries
    for i in range(1, K + 1):
        segment_index = i * segment_size
        selected_points.append(sorted_numbers[segment_index])
    
    return sorted(selected_points)

def unique_elements(arr):
    unique_list = []
    for item in arr:
        if item not in unique_list:
            unique_list.append(item)
    return unique_list


def neglect_close_points(all_crossing_points, new_points,s, tolerance=5):
    updated_points = all_crossing_points.copy()
    for new_point in new_points:
        close = False
        for existing_point in all_crossing_points:
            if abs(np.argmin(abs(s -existing_point)) - np.argmin(abs(s -new_point))) < tolerance:
                close = True
                break
        if not close:
            updated_points.append(new_point)
    return updated_points


from itertools import combinations


index_combinations=[[0,1]]
comb = 3
while comb < K + 5:
    indices = np.arange(comb)
    new_combinations = np.array(list(combinations(indices, 2)))
    index_combinations = np.vstack([index_combinations, new_combinations])
    _, unique_indices = np.unique(index_combinations, axis=0, return_index=True)
    index_combinations = index_combinations[np.sort(unique_indices)]
    comb += 1
    
all_crossing_points = []
k=0
while len(all_crossing_points)<K and k<len(index_combinations):
    i=index_combinations[k][0]
    j=index_combinations[k][1]
    #print(i,j)
    new_points = find_crossing_points(splines[i], splines[j], s_list, tolerance=1)
    all_crossing_points = neglect_close_points(all_crossing_points, new_points,s, tolerance=2.5)

    k += 1

print((all_crossing_points))
#print('all_crossing_points',all_crossing_points)
#print('all_crossing_points',most_separated(all_crossing_points,K),sorted(s_list))

for i in range(len(all_crossing_points)):
    y=all_crossing_points[i]
    #plt.plot((y,y),(1,3),color='k')

#plt.title(complete)


def find_element_indices(subvector, vector):
    indices = []
    for element in subvector:
        try:
            indices.append(vector.index(element)+1) #plus one for MATLAB
        except ValueError:
            pass  # If element not found in vector, ignore

    return indices
print(sorted(find_element_indices(all_crossing_points,sorted(s_list))))
np.savetxt('bd_s.txt', all_crossing_points)






