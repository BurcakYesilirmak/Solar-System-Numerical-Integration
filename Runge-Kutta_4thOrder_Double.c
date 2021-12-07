#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void energies(double r[3],double v[3], double *E, double *H)  {
    double v2 =v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    double r2 =r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    double norm_r=sqrt(r2);
    // jupiter mass : 9.547919e-4 wrt Sun
    
    double ek=0.5*v2;           // kinetic energy 
    double epot = -1/norm_r;    // potential energy
    *E= ek + epot;              // mechanical energy

    double i = r[1] * v[2] - r[2] * v[1] ;  // i = R1*V2-R2V1 
    double j = r[2] * v[0] - r[0] * v[2] ;  // j = R2*V0-R0V2 
    double k = r[0] * v[1] - r[1] * v[0] ;  // k = R0*V1-R1V0
    *H = sqrt(i*i+j*j+k*k);                 // angular momentum                              
    } 

void semimajor_eccentricity(double r[3], double v[3], double *a, double *e) {
   double k2= 2.9619474286664206E-04;
   double norm_r,norm_v, r1, v1, d1;
   double x, y2, z;
   double var1, var2, var3;
   
   r1= r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
   v1= v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; 
   d1= r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

   norm_r = sqrt(r1);
   norm_v = sqrt(v1);
   var1=  v1/k2 - 1.0/norm_r;
   var2= d1/k2;
   x  = r[0]*var1 - v[0]*var2; 
   y2 = r[1]*var1 - v[1]*var2;
   z  = r[2]*var1 - v[2]*var2;
   var3= x*x + y2*y2+ z*z;

   *e = sqrt(var3);                           // eccentricity vector norm
   *a = 1/((2.0/norm_r) - v1/k2);             // semi-major axis length   
   }    

// equation of motion's right-hand side evaluation as function:
void f(double r[3], double v[3], double *kr_i, double *kv_i)  {
    double k2= 2.9619474286664206E-04;
    double dot = r[0]*r[0]+ r[1]*r[1]+ r[2]*r[2];
    double r3=dot*sqrt(dot); 

        for(int i=0; i<3; i++)                {
         kv_i[i] = -r[i]*k2/r3; 
         kr_i[i] = v[i];                      }               }    


 // 4th Order Runge Kutta Method :
void RK4(double r[3], double v[3], double dt) {
    double k_r1[3], k_v1[3], k_r2[3], k_v2[3];
    double k_r3[3], k_v3[3], k_r4[3], k_v4[3];
    double v_temp1[3], v_temp2[3], v_temp3[3];
    double r_temp1[3], r_temp2[3], r_temp3[3];
    double kr[3], kv[3];

    f(r, v, k_r1, k_v1);
        for(int i=0; i<3; i++)               {
            v_temp1[i]= v[i]+k_v1[i]*0.5*dt;
            r_temp1[i]= r[i]+k_r1[i]*0.5*dt; 
                                             }

    f(r_temp1, v_temp1, k_r2, k_v2);
        for(int i=0; i<3; i++)               {
            v_temp2[i]= v[i]+k_v2[i]*0.5*dt;
            r_temp2[i]= r[i]+k_r2[i]*0.5*dt;
                                             }

    f(r_temp2, v_temp2, k_r3, k_v3);
        for(int i=0; i<3; i++)               {
            v_temp3[i]= v[i]+k_v3[i]*dt;
            r_temp3[i]= r[i]+k_r3[i]*dt;    
                                             }

    f(r_temp3, v_temp3, k_r4, k_v4);

        for(int i=0; i<3; i++)               {
           kv[i] = (k_v1[i] +2*(k_v2[i]+k_v3[i])+ k_v4[i])/6.0;
           kr[i] = (k_r1[i] +2*(k_r2[i]+k_r3[i])+ k_r4[i])/6.0;
           v[i] += kv[i]*dt;
           r[i] += kr[i]*dt;                 }
    }

int main() {
    double t, dt, T, a, e, energy, angmomentum; 
    //  r[3]={X,Y,Z}  V[3]={Vx,Vy,Vz}  
    // DATA: NASA HORIZON https://ssd.jpl.nasa.gov/horizons/app.html#/
    // JUPITER INITIAL VECTORS 2018-Jun-05 
    double r[]= {-3.460167504309613E+00, -4.149454064629457E+00, 9.465721330038770E-02};
    double v[]= {5.709741990408655E-03, -4.481465873394258E-03, -1.091471606521913E-04};

    // time parameters setting:
    T=  2*M_PI*1000;
    dt= 1e-2;
    
    // error measure/comparison with long double 
    long double  A=  5.20273584355200957629L;       // semi-major axis lenght
    long double EC=  4.88056797545033895064e-02L;   // eccentricity 
    long double  E= -1.85032942351832700980e-01L;   // mechanical energy of system
    long double  H=  3.92090843714996407606e-02L;   // angular momentum of system
    
   semimajor_eccentricity(r, v, &a, &e);
   energies(r, v, &energy, &angmomentum);
   // printf("%.20lf %.20lf %.20lf %.20lf\n", a , e, energy, angmomentum);
   // exit(1);
            
    FILE *kk, *pp;
    kk=fopen("RK4jupiter_double_pos_vel.txt","w");
    pp=fopen("RK4jupiter_double_a_e_energy.txt","w");
    t=0.0 ;   
    int Count=0.0;
    while ( t <= T )  {          // ( t <= tmax/dt )
       RK4(r, v, dt);     
       semimajor_eccentricity(r, v, &a, &e); 
       energies(r, v, &energy, &angmomentum);   
       fprintf(kk,"%d %.20lf %.20lf %.20lf %.20lf %.20lf %.20lf \n", Count, r[0], r[1], r[2], v[0], v[1], v[2]);
       fprintf(pp,"%d %.15Le %.15Le %.15Le %.15Le\n", Count, (A-a)/A, (EC-e)/EC, (E-energy)/E, (H-angmomentum)/H);
       t +=dt;       
       Count++; }
    fclose(kk);
    fclose(kk);
    return 0;  }
