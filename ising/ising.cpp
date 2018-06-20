#include <iostream>
#include <cmath>
#include "nr3.h"
#include "ran.h"
#include <vector>
#include <fstream>
using namespace std;

int main ()
{ Ran  ran(67); //initializing random number, pick any integer of you choice
  

  double T;  // specifying temperature  (T_critical ~ 2.2)
  cout << "inpute temperature : " << endl;
  cin >> T;


    int dim=10;     // specifying the dimension of the array
    int    m=1000000;  // specigying # of iteration for random sampling
    double b=1;      //  interaction parameter E = b*sum(Si*Sj)  
    int  a=0;        // number of Aluminum
   
    double  x ; 
    int arr[dim][dim][dim];
    for ( int i=0; i<dim ; i=i+1 ){    //creating cystal as an array 
       for (int j=0; j<dim ; j=j+1){
          for (int c=0; c<dim ; c=c+1){
            x  =  2*ran.doub() -1  ;
               if ( x>0)              //  specifying spin up or down
                  arr[i][j][c] = 1;
               else    arr[i][j][c] = -1;    }}}

   int s=0;
   while ( s < a ){             //  placing Al randomly
    int l = ran.doub()*dim;
    int k = ran.doub()*dim;
    int q = ran.doub()*dim;
       arr[l][k][q] = 0;
       s=s+1; }
 
    for ( int i=0; i<dim; i=i+1){   //  viewing the spin configuration 
      for ( int j=0; j<dim;  j=j+1){
          cout << arr[i][j][0] << "  " ;  }
          cout << endl;  }

   int n=0;             
   while ( n < m ) {           // m  number of iteration for random sampling
    int  i = ran.doub()*dim ;   // picking some random site
    int  j = ran.doub()*dim ;
    int  k = ran.doub()*dim ; 
      if ( arr[i][j][k] == 0) goto loop; // if the site is Al, skip calculation
    double  u, top, right, left, bottom, above, below; 
      if ( i==0 )                        // applying boundary condition
         top = arr[i][j][k]*arr[dim-1][j][k] ;
        else top = arr[i][j][k]*arr[i-1][j][k] ;
      if ( j==0 )
         left = arr[i][j][k]*arr[i][dim-1][k] ;
        else  left = arr[i][j][k]*arr[i][j-1][k] ;
      if ( i==dim-1 )
           bottom = arr[i][j][k]*arr[0][j][k] ;
        else  bottom = arr[i][j][k]*arr[i+1][j][k] ;
      if ( j==dim-1 )
           right = arr[i][j][k]*arr[i][0][k] ;
        else right = arr[i][j][k]*arr[i][j+1][k] ;
      if ( k==0 )
           below = arr[i][j][k]*arr[i][j][dim-1] ;
        else  below = arr[i][j][k]*arr[i][j][k-1] ;
      if ( k==dim-1 )
           above = arr[i][j][k]*arr[i][j][0] ;
        else above = arr[i][j][k]*arr[i][j][k+1] ;      
    // calculating total energy at a site i,j,k
    u = top + right + bottom + left + above + below; 
       if ( u <= 0 ) 
           // if u<=0 flipping reduce energy
           arr[i][j][k] = -arr[i][j][k] ; 
       else if ( ran.doub() < exp(-u/T) ) // boltzmann factor decides flipping
           arr[i][j][k] = -arr[i][j][k] ;
       else  arr[i][j][k] = arr[i][j][k] ;  // otherwise no flipping
     loop:    n=n+1 ; }

    cout << endl;          
   
    for ( int i=0; i<dim; i=i+1)  {   //  viewing the spin configuration
       for ( int j=0; j<dim;  j=j+1)  {
       cout << arr[i][j][0] << "  " ;  }
        cout << endl;  }

    int sum=0;                    // calculating magnetization 
    for ( int i=0; i<dim; i++) {
       for ( int j=0; j<dim; j++ ) {
          for ( int k=0; k<dim; k++ ) {
           sum = sum + arr[i][j][k] ;  }}}
     double  mg  = abs ( ( (double)sum ) / ( dim*dim*dim -a ) )  ;     

    double sum2 = 0;                 //  calculating average energy
    for ( int i=0; i<dim-1; i++) {
       for ( int j=0; j<dim-1; j++) {
          for ( int k=0; k<dim-1; k++) {
           sum2 = sum2 + arr[i][j][k]*arr[i+1][j][k] + arr[i][j][k]*arr[i][j+1][k] +arr[i][j][k+1] ; }}}
   // caluculatinge energy at 3 surface
    for ( int i=0; i<dim-1 ; i++) {
       for ( int j=0; j<dim-1; j++) {
        sum2 = sum2 + arr[i][j][dim-1]*arr[i+1][j][dim-1] + arr[i][j][dim-1]*arr[i][j+1][dim-1] ; }}
    for ( int i=0; i<dim-1 ; i++) {
       for ( int j=0; j<dim-1; j++) {
        sum2 = sum2 + arr[dim-1][i][j]*arr[dim-1][i+1][j] + arr[dim-1][i][j]*arr[dim-1][i][j+1] ; }}
    for ( int i=0; i<dim-1; i++) {
       for ( int j=0; j<dim-1; j++) {
        sum2 = sum2 + arr[i][dim-1][j]*arr[i+1][dim-1][j] + arr[i][dim-1][j]*arr[i][dim-1][j+1] ; }}
   // calculating energy at 3 edges
    for ( int i=0; i<dim-1 ; i++ ) {
        sum2 = sum2 + arr[i][dim-1][dim-1] ;
        sum2 = sum2 + arr[dim-1][i][dim-1] ;
        sum2 = sum2 + arr[dim-1][dim-1][i] ; }
   //  getting avg energy
     double E = b * ((double)sum2) / ( ( dim*dim*dim - a ) );

     cout << "magnetization = " << mg << endl;
     cout << "avg energy = " << E << endl;  
     cout << "temperature = " << T << endl;
     cout << "# of Al = " << a << endl;

return 0;

} 


