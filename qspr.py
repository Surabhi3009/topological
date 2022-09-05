

#!pip install kora -q
#import kora.install.rdkit

import numpy as np
import pandas as pd
from rdkit import Chem
import networkx as nx
import matplotlib.pyplot as plt
import warnings
#import streamlit as st
warnings.filterwarnings('ignore')
def mol_to_nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                  atomic_num=atom.GetAtomicNum(),
                  is_aromatic=atom.GetIsAromatic(),
                  atom_symbol=atom.GetSymbol())
        
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                  bond.GetEndAtomIdx(),
                  bond_type=bond.GetBondType())
        
       
    return G

from urllib.request import urlopen
from urllib.parse import quote

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'

import requests


CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"


def smiles_to_iupac(smiles):
    rep = "iupac_name"
    url = CACTUS.format(smiles, rep)
    response = requests.get(url)
    response.raise_for_status()
    return response.text

def topological_ind(compound):
  # define the smiles string and covert it into a molecule sturcture ------------
  #def smiles(compound):
      #c_smiles = CIRconvert(compound)
  c_mol = Chem.MolFromSmiles(compound)

  # conver rdkit object to networkx object --------------------------------------
  c_nx = mol_to_nx(c_mol)

  c_atom = nx.get_node_attributes(c_nx, 'atom_symbol')

  color_map = {'C': 'cyan',
              'O': 'orange',
              'N': 'magenta',
              'H': 'blue'}  

  c_colors = []
  for idx in c_nx.nodes():
      if (c_nx.nodes[idx]['atom_symbol'] in color_map):
          c_colors.append(color_map[c_nx.nodes[idx]['atom_symbol']])
      else:
          c_colors.append('gray')
      
  nx.draw(c_nx,
          labels=c_atom,
          with_labels = True,
          node_color=c_colors,
          node_size=800)

  plt.show()

  # print out the adjacency matrix ---------------------------------------------- 
  matrix = nx.to_numpy_matrix(c_nx)
  #print(matrix)
  print(nx.wiener_index(c_nx) )
 
  # print Wiener index -----------------------------------------------------------
  G=c_nx
  n=G.number_of_nodes()
  #print(n)
  
  
#degree based
  #first zagreb index
  index=[]
  index.append(compound)
  degree=np.sum(matrix, axis = 1)
  sq_deg=np.square(degree)
  #sq_deg
  M1=np.sum(sq_deg)
  #print("first zagreb index:", M1)
  index.append(M1)
  #print(degree)
  
  #randic---------------------------------------------------
  degree_mat=[]
  for i in degree:
    #for j in i:
    a=int(i)
    degree_mat.append(a)
  #print(degree_mat)
  
  randic=[] 
  rr=[] 
  m2=[]
  abc=[]
  azi=[]
  g=[]
  s=[]
  
  for i in range(0,n):
    for j in range(0,n):
      matrix_upper=np.triu(matrix, 1)
      s_adj=0
      if matrix_upper[i,j]!=0:
        p_adj=degree_mat[i]*degree_mat[j]
        s_adj=degree_mat[i]+degree_mat[j]
        ab=np.sqrt((s_adj-2)/p_adj)
        if s_adj<=2:
            az="Augemented Zagreb index not applicable for this compound"
        else:
            az=((p_adj)/(s_adj-2))**3
        ag=(2*np.sqrt(p_adj)/s_adj)
        sa=np.sqrt(1/s_adj)
        m2.append(p_adj)
        
        r=np.sqrt(1/p_adj)
        reci=np.sqrt(p_adj)
        randic.append(r)
        
        rr.append(reci)
        
        abc.append(ab)
        
        azi.append(az)
        
        g.append(ag)
        
        s.append(sa)
        
  #print("randic index:",np.sum(randic))   
  randic1=np.sum(randic)
  index.append(randic1) 
#recipocal randic--------------------------------------------------------------
  #print("recipocal randic:",np.sum(rr))
  ranrec=np.sum(rr)
  index.append(ranrec)
# second zagreb index----------------------------------------------------------
  #print("second zagreb index: ",np.sum(m2))
  randic1=np.sum(m2)
  index.append(randic1)

#forgetten index---------------------------------------------------------------
  f = []
  for i in degree_mat:
    f.append(pow(i, 3))
  
# printing result
  #print("forgettern index:",np.sum(f))
  forg=np.sum(f)
  index.append(forg)
#ABC index------------------------------------------------------------------------
  #print("ABC index:",np.sum(abc))  
  abc1=np.sum(abc)
  index.append(abc1)
#AZI index------------------------------------------------------------------------
  #print("augmented zagreb:",np.sum(azi))
  if s_adj<=2:
    azi1=" not applicable for this compound"
  else:
    azi1=np.sum(azi)
  index.append(azi1)
#”Geometric - Arithmetic index--------------------------------------------------
  #print("Geometric - Arithmetic index:",np.sum(g))
  g1=np.sum(g)
  index.append(g1)
#sum connectivity index ---------------------------------------------------------
  #print("sum connectivity index:",np.sum(s))
  s1=np.sum(s)
  index.append(s1)
  #print(index)
  #global df
  #df = pd.DataFrame(columns=["coumpund",'wiener index','hyper wiener',"first zagreb",'randic','recipocal','second zagreb','forgetten','ABC','augemented','ga','sum'])
  #df = df.append(pd.DataFrame([index], columns=df.columns), ignore_index=True)
  #print(df)
  #df.to_csv('file1.csv')
# define the function for coverting rdkit object to networkx object -----------
  print(index)
  return plt,index



#c=input("enter the isomelies of compound")
#c=smiles_to_iupac(c)
#c=CIRconvert(c)
#topological_ind(c)
