#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLParser import getionimage

# read in the .imzML file and get an image at a desired m/z
f = ImzMLParser('./20220304_ctrl_neg.imzML')
img = getionimage(f, 419.2562, tol=0.1, z=1,reduce_func=sum) # img stored as 2D numpy array
# plt.imshow(img, cmap='viridis', interpolation ='nearest')
px.imshow(img)

import cv2
blur = cv2.blur(img, (1,3))
blur = cv2.blur(blur, (3,1))
blur=diag_filter(blur)
px.imshow(blur)

import copy

def diag_filter(img):
    temp_img=copy.deepcopy(img)
    for i in range(1,img.shape[0]-1):
        temp_img[i][i]=(img[i-1][i-1]+img[i][i]+img[i+1][i+1])/3
    
    temp_img2=copy.deepcopy(temp_img)
    for m in range(temp_img.shape[0]-2,0,-1):
        temp_img2[m][temp_img.shape[0]-1-m]=(temp_img[m+1][temp_img.shape[0]-2-m]+temp_img[m][temp_img.shape[0]-1-m]+temp_img[m-1][temp_img.shape[0]-m])/3
    return temp_img2

test=np.array([[1,3,8,10,5],[1,3,8,10,5],[1,3,8,10,5],[1,3,8,10,5],[1,3,8,10,5]])
diag_filter(test)


for num in range(6, 2, -1) :
    print(num)

# dict - {index: x,y}
coord_dict = {}
for i in range(len(f.coordinates)):
  coord_dict[i] = (f.coordinates[i][0], f.coordinates[i][1])


# feature m/z
mz_list = [255.2289,281.2476,283.2657,303.2279,305.2434,329.2473,331.2618,417.2382,419.2562,478.2967,559.3104]

def mean_filter(img):
  # apply mean filter with 3*3 kernal to smooth image, return as 2D numpy array
  img_f = img
  for i in range(2,img.shape[0]-2):
    for j in range(2,img.shape[1]-2):
      block = img[i-1:i+1+1,j-1:j+1+1]
      m = np.mean(block,dtype=np.float32)
      img_f[i][j] = int(m)
  return(img_f)

img_f = mean_filter(img)
px.imshow(img_f)

def plot_by_mz(f,mz_list):
  # Only apply to <= 12 m/z values since will plot at most 2*6 images
  img_f_list = list()
  for mz in mz_list:
    img = getionimage(f, mz, tol=0.1, z=1,reduce_func=sum)
    img_f = mean_filter(img)
    img_f_list.append(img_f)
  fig = make_subplots(rows=2, cols=-(len(img_f_list)//-2),subplot_titles=mz_list)
  for i in range(len(img_f_list)):
    # fig.add_trace(px.imshow(img_f_list[i]).data[0], i%2+1, -((i+1)//-2))
    if i+1 <= -(len(img_f_list)//-2):
      fig.add_trace(px.imshow(img_f_list[i]).data[0],1,i+1)
    else:
      fig.add_trace(px.imshow(img_f_list[i]).data[0],2,(i+(len(img_f_list)//-2)+1))
  fig.show()
  return()

plot_by_mz(f,mz_list) # this method can be used to examine if a filter is suitable for some m/z


# In[ ]:


def find_cell (img,i):
  # find cell areas - pixel with intensity >= given i
  # intensities of non-cell pixels are set to 0, return as 2D numpy array
  img_cell = img
  img_cell[img_cell<i] = 0
  return(img_cell)


# In[ ]:


mz_list_img = list()
for mz in mz_list:
  img = getionimage(f, mz, tol=0.1, z=1,reduce_func=sum)
  img_f = mean_filter(img)
  img_cell = find_cell(img_f,1000)
  mz_list_img.append(img_cell)


# In[ ]:


mz_list_img = list()
for mz in mz_list:
  img = getionimage(f, mz, tol=0.1, z=1,reduce_func=sum)
  img_f = blur
  img_cell = find_cell(img_f,1000)
  mz_list_img.append(img_cell)


# In[ ]:


def union (l):
  # union (x,y) of cells identified from different m/z values
  # set cell=1 non-cell=0
  dim = l[0].shape
  union_cell = np.zeros((dim[0],dim[1]))
  non_zeros = list()
  for n in l:
    non_zeros.append(np.argwhere(n>0))
  ones = list()
  for i in non_zeros:
    for (x,y) in i:
      if((x,y) not in ones):
        ones.append((x,y))
  for x, y in ones:
    union_cell[x][y]=1
  return(union_cell)


# In[ ]:


uc = union(mz_list_img)
uc


# In[ ]:


cc = dfs(uc)
len(cc)


# In[ ]:


img_cell = find_cell(img_f,1400)
plt.imshow(img_cell, cmap='viridis', interpolation ='nearest')


# In[ ]:


def dfs(img_cell):
  # depth first search used to find (x,y) of each cell
  # return a list of (x,y) coordinates for each cell
  # first set cell=1 and non-cell=0
  grid = img_cell
  grid[grid != 0]=1
  seen = set()
  l_all = list()
  for r0, row in enumerate(grid):
    for c0, val in enumerate(row):
      l = list()
      if val and (r0, c0) not in seen:
        stack = [(r0, c0)]
        seen.add((r0, c0))
        l.append((r0,c0))
        while stack:
          r, c = stack.pop()
          for nr, nc in ((r-1, c), (r+1, c), (r, c-1), (r, c+1)):
            if (0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and grid[nr][nc] and (nr, nc) not in seen):
              stack.append((nr, nc))
              seen.add((nr, nc))
              l.append((nr,nc))
      l_all.append(l)
  return list(filter(None, l_all))


# In[ ]:


coord_cell = dfs(img_cell)


# In[ ]:


import pandas as pd
df = pd.DataFrame(columns=['x','y','cell_num'])


# In[ ]:


i = 1 # cell number
for l in coord_cell:
    for x,y in l:
        df.loc[len(df.index)] = [y, x, i] 
    i = i+1


# In[ ]:


df['t'] = df['cell_num'].astype('category')
px.scatter(df, x='x',y='y',color='t')


# In[ ]:


x = list()
y = list()
for i,j in coord_cell:
    x.append(i)
    y.appen


# In[ ]:


spec = ImzMLParser.getspectrum(f,117)
spec2 = ImzMLParser.getspectrum(f,118)
spec[0]


# In[ ]:


spec2[0]


# In[ ]:


# get mass spectra using coord_dict
for coords in coord_cell:
  for (x,y) in coords:
    index = coord_dict[(x,y)]
    spectrum = ImzMLParser.getspectrum(f,index)


# In[ ]:


# TODO: based on coords, extract mass spectra
def mass_spectra(coord_cell):
  for coords in coord_cell:
    


# In[ ]:


from google.colab import drive
drive.mount('/content/drive')
import numpy as np

f = open("/content/drive/My Drive/cell_filter.csv")
cell_filter = np.genfromtxt(f, delimiter=",")

cell_filter


# In[ ]:


import csv
def write_to_csv(l_all):
  # additionally, flip x and y 
  with open('/content/drive/My Drive/cell_pixel_coords_py.csv', 'w') as f:
    write = csv.writer(f)
    i = 1 # cell number
    for l in l_all:
      for x,y in l:
        row = [y,x,i]
        write.writerow(row)
      i = i+1
  return


# In[ ]:


write_to_csv(cc)


# In[ ]:




