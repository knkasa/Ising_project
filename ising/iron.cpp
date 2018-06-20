#include <iostream>
#include <cmath>
#include "nr3.h"
#include "ran.h"
#include <vector>
#include <fstream>
using namespace std;

// Ising model of Body-center cube diluted with non-magnetic Aluminum

int main ()
{ Ran  ran(21);

 ofstream myfile;
 myfile.open("iron.OUT");
 for ( double T=0.0; T<5; T=T+0.1 ) {

    int dim=21;     // specifying the dimension of the array (should be odd #)
 //   double  T;  // specifing temperature T
   // cout << "enter temperature :  ";
   // cin >> T;    
    int    m=20000000;  // specigying # of iteration for random sampling
    double b=1;      //  interaction parameter E = b*sum(Si*Sj)  
    int  a=0;        // number of Aluminum

    int arr[dim][dim][dim] ;
    int num;        // number of Fe
     num = pow(double(dim-1)/2,3) + pow(double(dim+1)/2,3) -a; 

    double  x ;     
    for ( int i=0; i<dim ; i=i+1 ){    //creating cystal as an array 
       for (int j=0; j<dim ; j=j+1){
          for (int c=0; c<dim ; c=c+1){
            x  =  2*ran.doub() -1  ;
               if ( x>0)              //  specifying spin up or down
                  arr[i][j][c] = 1;
               else    arr[i][j][c] = -1;    }}}

    for ( int i=1; i<dim-1 ; i=i+2){  //constructing BCC structure
       for ( int j=0; j<dim ; j=j+1){  // working on odd layers 
          for ( int k=0; k<dim ; k=k+2){  
            arr[i][j][k] = 0 ;  }}}    // puting 0 for empty site
    for ( int i=0; i<dim; i=i+1){
       for ( int j=1; j<dim-1; j=j+2){
          for ( int k=0; k<dim; k=k+2){
            arr[i][j][k] = 0 ;  }}}
    for ( int i=0; i<dim; i=i+2){        // working on even layers
       for ( int j=0; j<dim; j=j+1){
          for ( int k=1; k<dim-1; k=k+2){
            arr[i][j][k]=0;  }}}
    for ( int i=0; i<dim; i=i+1){
       for ( int j=0; j<dim; j=j+2 ){
          for ( int k=1; k<dim-1; k=k+2 ){
            arr[i][j][k]=0;   }}}  
          
    for ( int i=0; i<dim; i=i+1){   //  viewing the spin configuration 
      for ( int j=0; j<dim;  j=j+1){
          cout << arr[i][j][0] << "  " ;  }
          cout << endl;  }  

   int s=0;
   while ( s < a ){             //  placing Al randomly a= # of Al
    loop2:
    int l = ran.doub()*dim;
    int k = ran.doub()*dim;
    int q = ran.doub()*dim;
   if ( arr[l][k][q] == 0 )    goto loop2;
     else arr[l][k][q] = 0 ;
         s=s+1;   }

   int n=0;
   while ( n < m ) {           // m  number of iteration for random sampling
    int  i = ran.doub()*dim ;   // picking some random site
    int  j = ran.doub()*dim ;
    int  k = ran.doub()*dim ; 
      if ( arr[i][j][k] == 0) goto loop; // if the site is Al or empty, skip calculation
    //  BCC structure has 8 nearest neighbors
    double  u, first, second, third, fourth, fifth, sixth, seventh, eighth; 
      if ( (i==0 || i==dim-1) && (j==0 || j==dim-1) && (k==0 || k==dim-1)  )  
     {  first = arr[i][j][k]*arr[dim-2][dim-2][dim-2]  ;
       second = arr[i][j][k]*arr[1][1][dim-2] ;
       third = arr[i][j][k]*arr[1][dim-2][dim-2] ;
       fourth = arr[i][j][k]*arr[dim-2][1][dim-1] ;
       fifth =  arr[i][j][k]*arr[dim-2][dim-2][1]  ;
       sixth = arr[i][j][k]*arr[1][1][1] ;
       seventh = arr[i][j][k]*arr[1][dim-2][1] ;
       eighth = arr[i][j][k]*arr[dim-2][1][1] ;  }
       if  ( (i==0 || i==dim-1) && (j==0 || j==dim-1) ) ;
     {  first = arr[i][j][k]*arr[dim-1][dim-2][k-1] ;
      second = arr[i][j][k]*arr[1][dim-2][k-1] ;
      third = arr[i][j][k]*arr[dim-2][1][k-1] ;
      fourth = arr[i][j][k]*arr[1][1][k-1] ;
      fifth = arr[i][j][k]*arr[dim-2][dim-2][k+1] ;
      sixth =  arr[i][j][k]*arr[1][dim-2][k+1] ;
      seventh =  arr[i][j][k]*arr[dim-2][1][k+1]  ;
      eighth =  arr[i][j][k]*arr[1][1][k+1]  ;  }
      if (  (i==0 || i==(dim-1)) &&  (k==0 || k==(dim-1)) )
    {   first = arr[i][j][k]*arr[dim-2][j-1][dim-2] ;
      second = arr[i][j][k]*arr[1][j-1][dim-2] ;
      third = arr[i][j][k]*arr[dim-2][j+1][dim-2] ;
      fourth = arr[i][j][k]*arr[1][j+1][dim-2] ;
      fifth = arr[i][j][k]*arr[dim-2][j-1][1] ;
      sixth =  arr[i][j][k]*arr[dim-2][j+1][1] ;
      seventh =  arr[i][j][k]*arr[1][j-1][1]  ;
      eighth =  arr[i][j][k]*arr[1][j+1][1]  ;  }
     if (  (j==0 || j==dim-1) &&  (k==0 || k==dim-1) )
     {  first = arr[i][j][k]*arr[i-1][dim-2][dim-2] ;
      second = arr[i][j][k]*arr[i+1][dim-2][dim-2] ;
      third = arr[i][j][k]*arr[i-1][1][dim-2] ;
      fourth = arr[i][j][k]*arr[i+1][1][dim-2] ;
      fifth = arr[i][j][k]*arr[i-1][dim-2][1] ;
      sixth =  arr[i][j][k]*arr[i-1][1][1] ;
      seventh =  arr[i][j][k]*arr[i+1][dim-2][1]  ;
      eighth =  arr[i][j][k]*arr[i+1][1][1] ;  }
     if (  (i==0 || i==dim-1) )
    {  first = arr[i][j][k]*arr[1][j-1][k-1] ;
      second = arr[i][j][k]*arr[1][j-1][k-1] ;
      third = arr[i][j][k]*arr[dim-2][j+1][k-1] ;
      fourth = arr[i][j][k]*arr[dim-2][j+1][k-1] ;
      fifth = arr[i][j][k]*arr[1][j-1][k+1] ;
      sixth =  arr[i][j][k]*arr[1][j+1][k+1] ;
      seventh =  arr[i][j][k]*arr[dim-2][j-1][k+1]  ;
      eighth =  arr[i][j][k]*arr[dim-2][j+1][k+1] ;  }
    if ( (j==0  || j==dim-1) )
    {  first = arr[i][j][k]*arr[i-1][1][k-1] ;
      second = arr[i][j][k]*arr[i+1][1][k-1] ;
      third = arr[i][j][k]*arr[i-1][dim-2][k-1] ;
      fourth = arr[i][j][k]*arr[i+1][dim-2][k-1] ;
      fifth = arr[i][j][k]*arr[i-1][1][k+1] ;
      sixth =  arr[i][j][k]*arr[i-1][1][k+1] ;
      seventh =  arr[i][j][k]*arr[i+1][dim-2][k+1]  ;
      eighth =  arr[i][j][k]*arr[i+1][dim-2][k+1] ;  }
    if ( (k==0 ||  k==dim-1) )
    {  first = arr[i][j][k]*arr[i-1][j-1][1] ;
      second = arr[i][j][k]*arr[i+1][j-1][1] ;
      third = arr[i][j][k]*arr[i-1][j+1][dim-2] ;
      fourth = arr[i][j][k]*arr[i+1][j+1][dim-2] ;
      fifth = arr[i][j][k]*arr[i-1][j-1][1] ;
      sixth =  arr[i][j][k]*arr[i-1][j+1][1] ;
      seventh =  arr[i][j][k]*arr[i+1][j-1][dim-2]  ;
      eighth =  arr[i][j][k]*arr[i+1][j+1][dim-2] ;  }
    else
   {   first = arr[i][j][k]*arr[i-1][j-1][k-1] ;
      second = arr[i][j][k]*arr[i+1][j-1][k-1] ;
      third = arr[i][j][k]*arr[i-1][j+1][k-1] ;
      fourth = arr[i][j][k]*arr[i+1][j+1][k-1] ;
      fifth = arr[i][j][k]*arr[i-1][j-1][k+1] ;
      sixth =  arr[i][j][k]*arr[i-1][j+1][k+1] ;
      seventh =  arr[i][j][k]*arr[i+1][j-1][k+1]  ;
      eighth =  arr[i][j][k]*arr[i+1][j+1][k+1] ;    }
         
    // calculating total energy at a site i,j,k
    u = first + second + third + fourth + fifth + sixth + seventh + eighth ; 
       if ( u <= 0 ) 
           // if u<=0 flipping reduce energy
           arr[i][j][k] = -arr[i][j][k] ; 
       else if ( ran.doub() < exp(-b*u/T) ) // boltzmann factor decides flipping
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
    double  mg  = abs ( ((double)sum) / (num) )  ; 
  
// calculating E
    double sum2 = 0;
    for ( int i=0; i<=dim-2; i++ ) {
       for ( int j=0; j<=dim-2; j++ ) {
          for ( int k=0; k<=dim-2; k++) {
            sum2 = sum2 +  arr[i][j][k]*arr[i+1][j+1][k+1] ;   }}}
   for ( int i=dim-1; i>=1; i=i-1) {
      for ( int j=dim-1; j>=1; j=j-1 ) {
         for ( int k=0; k<=dim-2; k=k+1) {
            sum2 = sum2 + arr[i][j][k]*arr[i-1][j-1][k+1];   }}}
   for ( int i=dim-1; i>=1; i=i-1 ) {
      for ( int j=0; j<=dim-2; j=j+1 ) {
         for ( int k=0; k<=dim-2; k=k+1 ) {
            sum2 = sum2 + arr[i][j][k]*arr[i+1][j-1][k+1] ; }}}
   for ( int i=0; i<=dim-2; i=i+1 ) {
      for ( int j=dim-1; j>=1; j=j-1 ) {
         for ( int k=0; k<=dim-2; k=k+1 ) {
            sum2 = sum2 + arr[i][j][k]*arr[i-1][j+1][k+1] ; }}}
   double E;
    E = ( ((double)sum2) / (num) ) ;


   cout << " numbero of Fe : " << num << endl;
   cout << "magnetization = " << mg << endl;      
   cout << "temperature = " << T << endl;
   cout << "# of Al = " << a << endl;
   cout << "ave E = " << E << endl;
 myfile << T << " " <<  mg << endl; }

myfile.close();
return 0;

} 


