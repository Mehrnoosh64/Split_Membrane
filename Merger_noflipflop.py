#################################################
# Merging the split lipids to the standard ones #
#################################################

#Import necessary libraries
import sys,os,shutil,re,math,scipy
import numpy as np
import pandas as pd
import math
from scipy.optimize import linear_sum_assignment
scipy.optimize.linear_sum_assignment
import warnings
warnings.filterwarnings('ignore')

#Read the initial .gro file
gro_file = input("Please enter the name of your .gro file (Example.gro):")
dg = pd.read_csv(gro_file, skiprows=1)
dg.drop(dg.tail(1).index,inplace=True) # drop the last row
dg.columns.values[0] = 'A'

#Find and fix if there is any bug in the .gro file
dg['A']=dg['A'].str.slice(stop=44)
for line in dg['A']:
    if line[-1] == ' ':
        bug_index = dg.index.get_loc(dg.index[dg['A'] == line][0])
        line = ' ' + line
        edited_line = line[-1] + line[1:-1] + line[0]
        dg.iloc[bug_index] = edited_line
dg['A']=dg['A'].str.slice(stop=44)

#Define size of each group based on the .gro file format:
#residue number(5 positions, integer)
#residue name(5 characters)
#atom name(5 characters)
#atom number(5 positions, integer)
#position(in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
sizes = [5, 5, 5, 5, 8, 8, 8]
#Generate the pattern string and compile it
pattern = re.compile(''.join([ f'(?P<Column{idx}>.{{{n}}})' for idx, n in enumerate(sizes, start=1)]))
#Generate the dataframe
df = dg.A.str.extract(pattern)
df.columns = ['residue number','residue name','atom name','atom number','X','Y','Z']
df['residue number'] = df['residue number'].values.astype(int)
df['residue name'] = df['residue name'].values.astype(str)
df['atom name'] = df['atom name'].values.astype(str) 
df['atom number'] = df['atom number'].values.astype(int)
df['X'] = df['X'].values.astype(float)
df['Y'] = df['Y'].values.astype(float)
df['Z'] = df['Z'].values.astype(float)
 
#Read the box size
with open(gro_file) as file:
    for line in file:
        pass
    footer = line

def Convert(string):
    xyz = list(string.split("  "))
    return xyz

box_x = float(Convert(footer)[1])
box_y = float(Convert(footer)[2])
box_z = float(Convert(footer)[3])       

#Divide dataframe based on molecule types
print("Processing...")

#Functions
def get_type1(dataframe):
    df_lenght=0
    for index in range (len(dataframe)):
        if not ( ("PO16" in dataframe['residue name'][index])\
        or ("SM" in dataframe['residue name'][index])\
        or ("SOL" in dataframe['residue name'][index])\
        or ("CL" in dataframe['residue name'][index])\
        or ("NA" in dataframe['residue name'][index])\
        or ("CHL1" in dataframe['residue name'][index]) ):
        	df_lenght+=1
    return df_lenght


def get_type2(dataframe):
    df_lenght=0
    for index in range (len(dataframe)):
        if (("PO16" in dataframe['residue name'][index])\
        or ("SM" in dataframe['residue name'][index])):
        	df_lenght+=1
    return df_lenght


def get_type3(dataframe):
    df_lenght=0
    for index in range (len(dataframe)):
        if (("SOL" in dataframe['residue name'][index])\
        or ("CL" in dataframe['residue name'][index])\
        or ("NA" in dataframe['residue name'][index])\
        or ("CHL1" in dataframe['residue name'][index])):
        	df_lenght+=1
    return df_lenght 


#Variables
df1_lenght = get_type1(df)
df2_lenght = get_type2(df)
df3_lenght = get_type3(df)

df1 = df[0:df1_lenght].copy()
df2 = df[0:df2_lenght].copy()
df3 = df[0:df3_lenght].copy()


#Start
from datetime import timedelta
from datetime import datetime as dtt

def TimestampMillisec64():
    return int((dtt.utcnow() - dtt(1970, 1, 1)).total_seconds() * 1000) 

ii=0
start = TimestampMillisec64()
df1_ii=0
df2_ii=0
df3_ii=0

for index in range (len(df)):
    if ("PO16" not in df.iloc[index,1])\
    and ("SM" not in df.iloc[index,1])\
	and ("SOL" not in df.iloc[index,1])\
	and ("CL" not in df.iloc[index,1])\
	and ("NA" not in df.iloc[index,1])\
	and ("CHL1" not in df.iloc[index,1]) :        
        df1.iloc[df1_ii] = df.iloc[index]
        df1_ii += 1   
        df1.reset_index(drop=True, inplace=True)         
                         
    if (("PO16" in df['residue name'][index]) or ("SM" in df['residue name'][index])):
        df2.iloc[df2_ii] = df.iloc[index] 
        df2_ii += 1
        df2.reset_index(drop=True, inplace=True)


    if (("SOL" in df['residue name'][index]) or ("CL" in df['residue name'][index]) or ("NA" in df['residue name'][index]) or ("CHL1" in df['residue name'][index])):
        df3.iloc[df3_ii] = df.iloc[index]
        df3_ii += 1
        df3.reset_index(drop=True, inplace=True)    

            
    ii+=1
    if ii % 10000 == 0:
        elapsed_time = (TimestampMillisec64() - start)/1000
        steps_made = ii
        steps_to_finish = len(df) - steps_made
        time_to_finish = (steps_to_finish / steps_made) * elapsed_time
        print ("Time to finish:", str(timedelta(seconds=time_to_finish)).split(".")[0])
#End


#Calculate minimum distance between C1 and C2 atoms
print("Merging the pairs!")
print("Processing...")

#Return [x,y,z], index_in_df
def get_coordinate(df, string):
    array = []
    index_array = []
    for i in range (len(df)):
        if df.loc[i]['atom name'] == string:
            array.append([df.loc[i]['X'], df.loc[i]['Y'], df.loc[i]['Z']])
            index_array.append(i)
    return array, index_array

#Calculate distance between C1 and C2 considering periodic box
def get_length(df, pos1, pos2):
    x_i = pos1[0]
    x_j = pos2[0]
    y_i = pos1[1]
    y_j = pos2[1]    
    z_i = pos1[2]
    z_j = pos2[2]
    
    #Convert to unit box
    x_i /= box_x
    x_j /= box_x
    y_i /= box_y
    y_j /= box_y
    z_i /= box_z
    z_j /= box_z
    
    #Minimal image distance, we omited the priodic condition: "if(x_len > box_x/2)"
    x_len = abs( x_i - x_j ) - math.floor( abs( x_i - x_j ) )
    y_len = abs( y_i - y_j ) - math.floor( abs( y_i - y_j ) )
    z_len = abs( z_i - z_j ) - math.floor( abs( z_i - z_j ) )
    
    length = math.sqrt( x_len**2 + y_len**2 + z_len**2 )
    
    return length
        
def get_distances(df):
    global p
    global p_index
    global q
    global q_index
    global minimal_assignment
    p, p_index = get_coordinate(df, "   DH")
    q, q_index = get_coordinate(df, "  DDH")
    distances = np.random.rand(len(p),len(q)) 
    
    for j in range (len(p)):
        for k in range (len(q)):
            distances[j][k] = get_length(df, p[j], q[k])
            # introduce bias for minimalization
            distances[j][k] *= distances[j][k]
    minimal_assignment = linear_sum_assignment(distances)
    return minimal_assignment
            
def get_pair(df):
    for i in range (len(p)):
        head_index = p_index[minimal_assignment[0][i]]
        tail_index = q_index[minimal_assignment[1][i]]
        print (minimal_assignment[0][i],":",head_index, "-->", minimal_assignment[1][i], ":", tail_index)
        
#Rearrange heads and tails
tmp_df = pd.DataFrame(columns = ['A'])
tmp3 = tmp_df.A.str.extract(pattern)
tmp3.columns = ['residue number','residue name','atom name','atom number','X','Y','Z']
tmp3.columns = tmp3.columns.str.replace('^ +', '')
tmp3['residue number'] = tmp3['residue number'].values.astype(float)
tmp3['residue name'] = tmp3['residue name'].values.astype(str)
tmp3['atom name'] = tmp3['atom name'].values.astype(str) 
tmp3['atom number'] = tmp3['atom number'].values.astype(float)
tmp3['X'] = tmp3['X'].values.astype(float)
tmp3['Y'] = tmp3['Y'].values.astype(float)
tmp3['Z'] = tmp3['Z'].values.astype(float) 
tmp1_ul = tmp3.copy()
tmp1_ll = tmp3.copy()
tmp2_ul = tmp3.copy()
tmp2_ll = tmp3.copy()

#Arrange POPC, POPS, POPE, and POPA
#Upper Leaflet
list_heads = []
list_tails = []

for i in range(len(df1)):
    if (" DH" in df1['atom name'][i] and df1['Z'][i] > box_z/2):
        if ( "PC" in df1['residue name'][i] ):
            list_heads.append(df1[i-27: i+1])
        if ( "PS" in df1['residue name'][i] ):
            list_heads.append(df1[i-20: i+1])
        if ( "PE" in df1['residue name'][i] ):
            list_heads.append(df1[i-18: i+1])                      
        if ( "PA" in df1['residue name'][i] ):
            list_heads.append(df1[i-9: i+1])     
    if ("DDH" in df1['atom name'][i] and df1['Z'][i] > box_z/2):
        if ( "PO " in df1['residue name'][i] ):
            list_tails.append(df1[i-1: i+107])                      
                       
get_distances(df1)        
get_pair(df1)                

for head,tail in zip(range(int(len(minimal_assignment[0])/2)), range(int(len(minimal_assignment[1])/2))):
    tmp1_ul = tmp1_ul._append( list_heads[head] )
    tmp1_ul = tmp1_ul._append( list_tails[tail] )
    
#Fix atom numbers
tmp1_ul.reset_index(drop=True, inplace=True) 

#Arrange POPC, POPS, POPE, and POPA
#Lower Leaflet
list_heads = []
list_tails = []

for i in range(len(df1)):
    if (" DH" in df1['atom name'][i] and df1['Z'][i] <= box_z/2):
        if ( "PC" in df1['residue name'][i] ):
            list_heads.append(df1[i-27: i+1])
        if ( "PS" in df1['residue name'][i] ):
            list_heads.append(df1[i-20: i+1])
        if ( "PE" in df1['residue name'][i] ):
            list_heads.append(df1[i-18: i+1])                      
        if ( "PA" in df1['residue name'][i] ):
            list_heads.append(df1[i-9: i+1])     
    if ("DDH" in df1['atom name'][i] and df1['Z'][i] <= box_z/2):
        if ( "PO " in df1['residue name'][i] ):
            list_tails.append(df1[i-1: i+107])                      
                       
get_distances(df1)        
get_pair(df1)                

for head,tail in zip(range(int(len(minimal_assignment[0])/2)), range(int(len(minimal_assignment[1])/2))):
    tmp1_ll = tmp1_ll._append( list_heads[head] )
    tmp1_ll = tmp1_ll._append( list_tails[tail] )
    
#Fix atom numbers
tmp1_ll.reset_index(drop=True, inplace=True) 

#Arrange SM16
#Upper Leaflet
list_heads = []
list_tails = []

for i in range(len(df2)):
    if (" N" in df2['atom name'][i] and df2['Z'][i] > box_z/2):
        if ( "SM" in df2['residue name'][i] ):
            list_heads.append(df2[i: i+28])
    if ("DDH" in df2['atom name'][i] and df2['Z'][i] > box_z/2):           
        if ( "PO16" in df2['residue name'][i] ):
            list_tails.append(df2[i-1: i+100])

get_distances(df2)        
get_pair(df2)         
            
for head,tail in zip(range(int(len(minimal_assignment[0])/2)), range(int(len(minimal_assignment[1])/2))):
    tmp2_ul = tmp2_ul._append( list_heads[head] )
    tmp2_ul = tmp2_ul._append( list_tails[tail] )   

#Fix atom numbers
tmp2_ul.reset_index(drop=True, inplace=True)

#Arrange SM16
#Lower Leaflet
list_heads = []
list_tails = []

for i in range(len(df2)):
    if (" N" in df2['atom name'][i] and df2['Z'][i] <= box_z/2):
        if ( "SM" in df2['residue name'][i] ):
            list_heads.append(df2[i: i+28])
    if ("DDH" in df2['atom name'][i] and df2['Z'][i] <= box_z/2):           
        if ( "PO16" in df2['residue name'][i] ):
            list_tails.append(df2[i-1: i+100])

get_distances(df2)        
get_pair(df2)         
            
for head,tail in zip(range(int(len(minimal_assignment[0])/2)), range(int(len(minimal_assignment[1])/2))):
    tmp2_ll = tmp2_ll._append( list_heads[head] )
    tmp2_ll = tmp2_ll._append( list_tails[tail] )   
    
#Fix atom numbers
tmp2_ll.reset_index(drop=True, inplace=True)

#Arrange other molecules in the box (Cholesterol, Water, CL, and NA)
list_cholesterol = []
list_water = []
list_ions = [] 
            
for i in range(len(df3)):
    if ("CHL1" in df3['residue name'][i]):            
        list_cholesterol.append(df3[i:i+1]) 
    if "SOL" in df3['residue name'][i]:            
        list_water.append(df3[i:i+1])
    if ("NA" in df3['residue name'][i] or "CL" in df3['residue name'][i]):            
        list_ions.append(df3[i:i+1])

tmp3 = tmp3.append( list_cholesterol )
tmp3 = tmp3.append( list_water )
tmp3 = tmp3.append( list_ions )

#Fix atom numbers
tmp3.reset_index(drop=True, inplace=True)

#Make the new dataframe
dataframes = [tmp1_ul, tmp1_ll, tmp2_ul, tmp2_ll, tmp3]
tmp = pd.concat(dataframes)
tmp.reset_index(drop=True, inplace=True)

for i in range(len(tmp)):
    if ('    N' in tmp.loc[i,('atom name')] and 'PC   ' in tmp.loc[i,('residue name')]):
        tmp['residue name'][i:i+136] = 'POPC '
        tmp['residue number'][i:i+136] = tmp['residue number'][i]
    if ('    N' in tmp.loc[i,('atom name')] and 'PS   ' in tmp.loc[i,('residue name')]):
        tmp['residue name'][i:i+129] = 'POPS '
        tmp['residue number'][i:i+129] = tmp['residue number'][i]  
    if ('    N' in tmp.loc[i,('atom name')] and 'PE   ' in tmp.loc[i,('residue name')]):
        tmp['residue name'][i:i+127] = 'POPE '
        tmp['residue number'][i:i+127] = tmp['residue number'][i] 
    if ('    P' in tmp.loc[i,('atom name')] and 'PA   ' in tmp.loc[i,('residue name')]):
        tmp['residue name'][i:i+118] = 'POPA '
        tmp['residue number'][i:i+118] = tmp['residue number'][i]
    if ('    N' in tmp.loc[i,('atom name')] and 'SM   ' in tmp.loc[i,('residue name')]):
        tmp['residue name'][i:i+129] = 'SM16 '
        tmp['residue number'][i:i+129] = tmp['residue number'][i]    
        
#Delete virtual sites
for i in range(len(tmp)):
    if ("DH" in tmp['atom name'][i]):
        tmp = tmp.drop(i)   

        
#Order based on Residue names, use the same order for the topol.top file!        
s1=tmp[(tmp['residue name'] == 'CHL1 ')]
s2=tmp[(tmp['residue name'] == 'POPE ')]
s3=tmp[(tmp['residue name'] == 'POPC ')]
s4=tmp[(tmp['residue name'] == 'POPS ')]
s5=tmp[(tmp['residue name'] == 'POPA ')]
s6=tmp[(tmp['residue name'] == 'SM16 ')]
s7=tmp[(tmp['residue name'] == 'SOL  ')]
s8=tmp[(tmp['residue name'] == 'NA   ')]
s9=tmp[(tmp['residue name'] == 'CL   ')]
dataframes = [s1,s2,s3,s4,s5,s6,s7,s8,s9]
ordered = pd.concat(dataframes)
ordered.reset_index(drop=True, inplace=True)

#New Topol.top
print("Please edit the topol.top file:")
num = 0
count = []
for resname,x in zip(['CHL1 ', 'POPC ', 'POPE ', 'POPS ', 'POPA ', 'SM16 ', 'SOL  ', 'NA   ', 'CL   '], [74, 134, 125, 127, 116, 127, 3, 1, 1]):
    if resname not in ordered.values:
        print(resname, "does not exists in system")
    else :
        count.append(int(ordered['residue name'].value_counts()[resname]/x))
        print(resname, count[num])
        num += 1
        
#Save the merged.gro file
def main():
    fmt = '%5d','%5s','%5s','%5d','%8.3f','%8.3f','%8.3f'
    header = '{0:^13s}\n{1:^6d}'.format('merged bilayer', len(ordered))
    with open(gro_file) as file:
        for line in file:
            pass
        footer = line
    np.savetxt('merged.gro', ordered.values, fmt=fmt, delimiter="",header=header,footer=footer,comments='')
        
if __name__ == "__main__":
    main() 
    
print("Merging is done! -->> merged.gro")
