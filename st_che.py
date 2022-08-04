import qspr
import streamlit as st



def display(f):
    l=["First zagreb index ",'Randic index ','Recipocal Randic index ',
       'Second zagreb index ',
       'Forgetten Topological index ','Atom -bond connectivity index',
       'Augemented Zagreb index ','Geometric - Arithmetic index ','Sum connectivity index ']
    
    st.subheader("Degree based Topological Indices")
    for i in range(0,len(l)):
      st.write(l[i],f[1][i+1])
    
try:    
    #if option=="cmpd name":
        c= st.text_input("Enter cmpd name")
        #st.write("Hello", your_name)
        #c=smiles_to_iupac(c)
        #st.success("Success")

        c=qspr.CIRconvert(c)
        f=qspr.topological_ind(c)
        st.subheader("Molecular graph")
        st.pyplot(fig=f[0])
        display(f)
        
    """else:
        c= st.text_input("Enter SMILES")
        #st.write("Hello", your_name)
        #c=smiles_to_iupac(c)
        #c=CIRconvert(c)
        #st.success("Success")

        r=qspr.topological_ind(c)
        st.subheader("Molecular graph")
        st.pyplot(fig=r[0])
        dispaly(r)"""
except AttributeError:
    pass
        
        
