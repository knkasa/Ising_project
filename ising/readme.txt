The file 'ising' contains various c++ code for this project.  The file 
name represents as follows   


 
ising.cpp  =  general ising model for a cubic lattice which a user can 
              specify a temperature to obtain magnetization.  The transition 
              temperature is around T= 2.2 



cube.cpp  =  same cubic lattice as ising.cpp but the code will produce a 
             plot 'magnetization vs. temperature' in cube.dat 



iron.cpp  =  this is the main code for the project 'Iron-Aluminum aloys'
             which has a body-center-cubic crystal structure.  The code 
             will produce a plot 'magnetization vs. temperature' in ironX.OUT 
             where X correspondes to the concentration of Aluminum 
             (e.g. iron2.OUT takes q=0.2)  (note that the temperature is not 
             scaled yet)  

external.cpp  =  this code is same as the  cubic lattice, but includes
                 external magnetic fields in the hamiltonial. It will
                 produce M vs. T plot in external.dat 



anti.cpp  =  this code will demonstrate anti-ferromagnets using cubic lattice.
             a user can specify a tempeture to visualize the spin configuratio
        
model.cpp  =  this code uses the modified ising model (Blume Capel Model) 
              using cubic lattice.  Please see ising.pdf for details.
      
model2.cpp  =  this code uses another modified ising model, spin cross over 
               nanochanin, using cubic lattice.  Please see ising.pdf for 
               details.  
