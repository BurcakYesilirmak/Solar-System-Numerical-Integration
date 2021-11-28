#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define C  134217729 

  typedef struct {
   double  xh1[3], xl1[3],   xh2[3],  xl2[3], xh3[3], xl3[3], xh4[3], xl4[3];
   double  xh5[3], xl5[3],   xh6[3],  xl6[3], xh7[3], xl7[3], xh8[3], xl8[3];

   double xh11[3], xl11[3], xh22[3], xl22[3], xh33[3], xl33[3], xh44[3], xl44[3];
   double xh55[3], xl55[3], xh66[3], xl66[3], xh77[3], xl77[3], xh88[3], xl88[3];
   double  xh9[3],  xl9[3], xh10[3], xl10[3], sch1[3], scl1[3], sch2[3], scl2[3];  
  } value;


  // Functions of Accurate Arithmetic Operations:
 void VeltkampSplit(double x, double *xh, double *xl)   {        
    double p = C * x;
    double q = x - p;
    *xh = p + q; 
    *xl = x - *xh;                                      }


 void F2Sum(double a, double b,double *xh, double *xl)  {
    double z;
    *xh= a+b;
    z = *xh-a;
    *xl= b-z;                                           }


 void TwoSum(double a, double b,double *xh, double *xl) {
    double a1 ,b1 ,del_a ,del_b;
    *xh= a+b;
    a1= *xh-b;
    b1= *xh-a1;
    del_a = a-a1;
    del_b = b-b1;
    *xl = del_a+del_b;                                   }

  //Pichat and Neumaier’s summation algorithm
 void Summation2(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3 ){  
    double t,s,u,v,z;
    TwoSum(xh1, xh2, &t, &s);
    TwoSum(xl1, xl2, &u, &v);
    z= v+s+u; 
    //F2Sum(t, z, xh3, xl3);
    *xh3=t;
    *xl3=z;                                                                               }   

  // Dekker's & Polynomial Multiplication
 void Mult(double xh1, double xl1, double xh2, double xl2,double *t, double *z) {
    double xh_1, xl_1, xh_2, xl_2, xh3, xl3;
    VeltkampSplit(xh1, &xh_1, &xl_1); 
    VeltkampSplit(xh2, &xh_2, &xl_2); 
    xh3 =  xh1*xh2;
    xl3 =  xh_1*xh_2 - xh3;
    xl3 +=  xh_1*xl_2;
    xl3 +=  xl_1*xh_2;
    xl3 +=  xl_1*xl_2+ xl1*xh2 + xh1*xl2+ xl1*xl2;
    F2Sum(xh3, xl3, t, z);
    return;                                                                     }

 void Division(double xh1, double xl1, double xh2, double xl2,double *xh3, double *xl3){
    double u1, u0, t, z;
    t =xh1/xh2;
    Mult(t, 0, xh2, 0, &u0, &u1); 
    z = xh1- u0;
    z -= u1;
    z += xl1;
    z -= t * xl2;
    z /= xh2;  
    //*xh3= t;
    //*xl3= z;
    F2Sum(t, z, xh3, xl3);
    return; } 
 
 // the square of the vector's itself
  void square_product(double xh[3], double xl[3],double *xh3, double *xl3) {        
      double u0, u1, v0, v1, t0, t1 ,t2, t3;
      Mult(xh[0], xl[0], xh[0], xl[0], &u0, &u1); // xh^2 & xl^2
      Mult(xh[1], xl[1], xh[1], xl[1], &v0, &v1); // yh^2 & yl^2
      Mult(xh[2], xl[2], xh[2], xl[2], &t0, &t1); // zh^2 & zl^2
      Summation2(u0, u1, v0, v1, &t2, &t3);
      Summation2(t0, t1, t2, t3, xh3, xl3);                               }

  void SquareRoot(double R_h[3], double R_l[3], double *xh3, double *xl3)  {  //as pairs
      double y0, y1;
      double p0, p1, t4, t5;

     square_product(R_h, R_l, &t4, &t5);
     y0 = sqrt(t4);                          //starting point to the square root 
     Mult(y0, 0, y0, 0, &p0, &p1); 
     y1 = t4 - p0;
     y1 -= p1;
     y1 += t5 ;
     if(y0<0) { 
     printf("error \n"); }
     y1 = (y1 / y0)*0.5;                   // y0 and y1 are pairs of the square root operation
     *xh3= y0;                             // Mult(y0, y1, t4, t5, xh3, xl3);    // k^2*sqrt(k)    
     *xl3= y1;                                                            } 

   void scalar_product(double xh[3], double xl[3],double yh[3], double yl[3], double *xh3, double *xl3) { 
      double u0, u1, v0, v1, t0, t1 ,t2, t3;
      Mult(xh[0], xl[0], yh[0], yl[0], &u0, &u1); // X^2
      Mult(xh[1], xl[1], yh[1], yl[1], &v0, &v1); // Y^2
      Mult(xh[2], xl[2], yh[2], yl[2], &t0, &t1); // Z^2
      Summation2(u0, u1, v0, v1, &t2, &t3);
      Summation2(t0, t1, t2, t3, xh3, xl3);                                                             }

   void semimajor_eccentricity(double rh[3], double rl[3], double vh[3], double vl[3], double *a, double *e) {

   double sc_h, sc_l, sc_h2, sc_l2, norm_rh, norm_rl, norm_vh, norm_vl;
   double div1, div2, div3, div4, div5, div6;
   double mult_h1, mult_l1, mult_h2, mult_l2, mult_h3, mult_l3;
   double mult_h4, mult_l4, mult_h5, mult_l5, mult_h6, mult_l6;
   double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8;
   double div_h, div_l, div_h1, div_l1, sum_h2, sum_l2;
   double mu_h=2.9619474286664206E-04, mu_l=-0.0000000000000000497815544e-4; // Gravitational Parameters
   
   //Mult(k_h, k_l, k_h, k_l, &kh, &kl); 
   scalar_product(vh, vl, vh, vl, &sc_h, &sc_l);        // v.v = v1
   scalar_product(rh, rl, vh, vl, &sc_h2, &sc_l2);      // r.v = d1
   SquareRoot(rh, rl, &norm_rh, &norm_rl);              // norm_r
   SquareRoot(vh, vl, &norm_vh, &norm_vl);              // norm_v

   Division(sc_h, sc_l, mu_h, mu_l, &div3, &div4);      //  v1/k^2 
   Division(1.0, 0, norm_rh, norm_rl, &div1, &div2);    //  1/norm_r

   Summation2(div3, div4, -div1, -div2, &sum1, &sum2);  //  var1= (v1/k^2)-1/norm_r
   Division(sc_h2, sc_l2, mu_h, mu_l, &div5, &div6);    //  var2=  d1/k^2      

   Mult(rh[0], rl[0], sum1, sum2, &mult_h1, &mult_l1);  //  r_vec[]*var1
   Mult(rh[1], rl[1], sum1, sum2, &mult_h2, &mult_l2);
   Mult(rh[2], rl[2], sum1, sum2, &mult_h3, &mult_l3);

   Mult(vh[0], vl[0], div5, div6, &mult_h4, &mult_l4);  //  v_vec[]*var2
   Mult(vh[1], vl[1], div5, div6, &mult_h5, &mult_l5);
   Mult(vh[2], vl[2], div5, div6, &mult_h6, &mult_l6);

   Summation2(mult_h1, mult_l1, -mult_h4, -mult_l4, &sum3, &sum4); // x
   Summation2(mult_h2, mult_l2, -mult_h5, -mult_l5, &sum5, &sum6); // y
   Summation2(mult_h3, mult_l3, -mult_h6, -mult_l6, &sum7, &sum8); // z 
   double sum_h[]= {sum3, sum5, sum7}, sq_h;
   double sum_l[]= {sum4, sum6, sum8}, sq_l;
   SquareRoot(sum_h, sum_l, &sq_h, &sq_l);              //  { ecc_h, ecc_l } => sqrt(x^2+y^2+z^2)

    // semi-major axis = 1 / ((2.0/norm_r) - v1/k2);
    Division(2.0, 0.0, norm_rh, norm_rl, &div_h, &div_l); 
    Summation2(div_h, div_l, -div3, -div4, &sum_h2, &sum_l2); 
    Division(1.0, 0.0, sum_h2, sum_l2, &div_h1, &div_l1); 
    a[0]= div_h1;
    a[1]= div_l1;
    e[0]= sq_h;
    e[1]= sq_l;
   }
   
   void energies(double rh[3], double rl[3], double vh[3], double vl[3], double *energy, double *angmomentum){
      double norm_rh, norm_rl, sc_h, sc_l;
      double ek_h, ek_l, pot_h, pot_l, mech_h, mech_l;
      double mulh, mull, mulh1, mull1, mulh2, mull2;
      double mulh3, mull3, mulh4, mull4, mulh5, mull5; 
      double sumh, suml, sumh1, suml1, sumh2, suml2;
      
      // kinetic and potential energies:
      scalar_product(vh, vl, vh, vl, &sc_h, &sc_l); // v.v
      SquareRoot(rh, rl, &norm_rh, &norm_rl); //norm_r
      Mult(sc_h, sc_l, 0.5, 0.0, &ek_h, &ek_l);
      Division(-1.0, 0.0, norm_rh, norm_rl, &pot_h, &pot_l);
      Summation2( ek_h, ek_l, pot_h, pot_l, &mech_h, &mech_l);
      
      // angular momentum: 
                         // i = R1*V2-R2V1 
      Mult(rh[1], rl[1], vh[2], vl[2], &mulh, &mull);
      Mult(rh[2], rl[2], vh[1], vl[1], &mulh1, &mull1);
      Summation2(mulh, mull, -mulh1, -mull1, &sumh, &suml);

                         // j = R2*V0-R0V2 
      Mult(rh[2], rl[2], vh[0], vl[0], &mulh2, &mull2);
      Mult(rh[0], rl[0], vh[2], vl[2], &mulh3, &mull3);
      Summation2(mulh2, mull2, -mulh3, -mull3, &sumh1, &suml1);

                          // k = R0*V1-R1V0 
      Mult(rh[0], rl[0], vh[1], vl[1],  &mulh4, &mull4);
      Mult(rh[1], rl[1], vh[0], vl[0],  &mulh5, &mull5);
      Summation2(mulh4, mull4, -mulh5, -mull5, &sumh2, &suml2);
   
       // ang_h, ang_l => sqrt(i^2 + j^2 + k^2) :
      double sumh3[]= {sumh, sumh1, sumh2}, sq_h;
      double suml3[]= {suml, suml1, suml2}, sq_l;
          SquareRoot(sumh3, suml3, &sq_h, &sq_l); 
     
      //results:
      energy[0]= mech_h;
      energy[1]= mech_l;
      angmomentum[0]= sq_h;
      angmomentum[1]=sq_l;
    } 


   void f(double r[3], double rl[3], double v[3], double vl[3], double t, double *krh, double *krl, double *kvh, double *kvl)  {
      double sc_h1, sc_l1, sc_h2, sc_l2, sc_h3, sc_l3, x, y;
      value  inp[3];
      double mu_h=2.9619474286664206E-04, mu_l=-0.0000000000000000497815544e-4; // gravitational parameters
      
      square_product(r, rl, &sc_h1, &sc_l1);             // k = y[0]*y[0]+y[1]*y[1]+y[2]*y[2] 
      SquareRoot(r, rl, &sc_h2, &sc_l2);                 // sqrt(k) 
      Mult(sc_h1, sc_l1, sc_h2, sc_l2, &sc_h3, &sc_l3);  // k*sqrt(k)
      Mult(mu_h, mu_l, -1.0, 0.0, &x, &y);
    
    for(int i=0; i<3; i++)   {
        Mult(x, y, r[i], rl[i], &inp[i].sch1[i], &inp[i].scl1[i]); 
        Division(inp[i].sch1[i], inp[i].scl1[i], sc_h3, sc_l3, &kvh[i], &kvl[i]); // -mu*R/(x^2+y^2+z^2)^3/2  
         krh[i]=  v[i];
         krl[i]= vl[i];   }  }

   void RK4(double r[3], double rl[3], double v[3], double vl[3], double dt) { 
       double  d1, d2;
       double kvh1[3], kvh2[3], kvh3[3], kvh4[3];
       double kvl1[3], kvl2[3], kvl3[3], kvl4[3];
       double krh1[3], krh2[3], krh3[3], krh4[3];
       double krl1[3], krl2[3], krl3[3], krl4[3];
       double  rh2[3],  rl2[3],  vh2[3],  vl2[3];
        double rh3[3],  rl3[3],  vh3[3],  vl3[3];
       double  rh4[3],  rl4[3],  vh4[3],  vl4[3]; 
        double  Rh[3],   Rl[3],   Vh[3],   Vl[3];
       double kvfh[3], kvfl[3], krfh[3], krfl[3]; 

      value  inp[3];

      f(r, rl, v, vl, dt, krh1, krl1, kvh1, kvl1);                                  // k1_v[]= { a1[0], a1[1], a1[2]                                                            

        for(int i=0; i<3; i++)           {
        Mult(kvh1[i], kvl1[i], 0.5, 0.0, &inp[i].xh1[i], &inp[i].xl1[i]);
         Mult(inp[i].xh1[i], inp[i].xl1[i], dt, 0.0, &inp[i].xh2[i], &inp[i].xl2[i]);
        Summation2(v[i], vl[i], inp[i].xh2[i], inp[i].xl2[i], &vh2[i], &vl2[i]);   // k2_ri

        Mult(krh1[i], krl1[i], 0.5, 0.0, &inp[i].xh3[i], &inp[i].xl3[i]);           
         Mult(inp[i].xh3[i], inp[i].xl3[i],  dt, 0.0, &inp[i].xh4[i], &inp[i].xl4[i]);             
        Summation2(r[i], rl[i], inp[i].xh4[i], inp[i].xl4[i], &rh2[i], &rl2[i]);    // r[i]+k2_ri*0.5*dt;
                                        }

      f(rh2, rl2, vh2, vl2, dt, krh2, krl2, kvh2, kvl2);                           //k2_v[]= { a2[0], a2[1], a2[2] }
                                                                                                                                 
      for(int i=0; i<3; i++)             {
         Mult(kvh2[i], kvl2[i], 0.5, 0.0, &inp[i].xh5[i], &inp[i].xl5[i]);
            Mult(inp[i].xh5[i], inp[i].xl5[i], dt, 0.0, &inp[i].xh6[i], &inp[i].xl6[i]);
        Summation2(v[i], vl[i], inp[i].xh6[i], inp[i].xl6[i], &vh3[i], &vl3[i]);   // k3_ri

         Mult(krh2[i], krl2[i], 0.5, 0.0, &inp[i].xh7[i], &inp[i].xl7[i]);           
            Mult(inp[i].xh7[i], inp[i].xl7[i], dt, 0.0, &inp[i].xh8[i], &inp[i].xl8[i]);
        Summation2(r[i], rl[i], inp[i].xh8[i], inp[i].xl8[i], &rh3[i], &rl3[i]);   // r[i]+k2_ri*0.5*dt;
                                        }

      f(rh3, rl3, vh3, vl3, dt, krh3, krl3, kvh3, kvl3);   
                                                                                   // k3_v[]= { a3[0], a3[1], a3[2] }
        for(int i=0; i<3; i++)           {
           Mult(kvh3[i], kvl3[i], dt, 0.0, &inp[i].xh9[i], &inp[i].xl9[i]);
        Summation2(v[i], vl[i],  inp[i].xh9[i],  inp[i].xl9[i], &vh4[i], &vl4[i]); // k4_ri
          Mult( krh3[i], krl3[i], dt, 0.0, &inp[i].xh10[i], &inp[i].xl10[i]);          
        Summation2(r[i], rl[i], inp[i].xh10[i], inp[i].xl10[i], &rh4[i], &rl4[i]); // r[i]+k3_ri*0.5*dt; 
                                         }
  
      f(rh4, rl4, vh4, vl4, dt, krh4, krl4, kvh4, kvl4);                           // k4_v[]= { a3[0], a3[1], a3[2] }
                                                        
        for(int i=0; i<3; i++)                                {
        //  v[i] = v[i] + dt*(k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])/6.0;   
        Summation2(kvh2[i], kvl2[i], kvh3[i], kvl3[i], &inp[i].xh11[i], &inp[i].xl11[i]);
        Mult(2.0, 0.0, inp[i].xh11[i], inp[i].xl11[i], &inp[i].xh22[i], &inp[i].xl22[i]);
        Summation2(kvh1[i], kvl1[i], kvh4[i], kvl4[i], &inp[i].xh33[i], &inp[i].xl33[i]);
        Summation2(inp[i].xh22[i], inp[i].xl22[i], inp[i].xh33[i], inp[i].xl33[i], &inp[i].xh44[i], &inp[i].xl44[i]);
        Division(inp[i].xh44[i], inp[i].xl44[i], 6.0, 0.0, &kvfh[i], &kvfl[i]);
        Mult(kvfh[i], kvfl[i], dt, 0.0, &Vh[i], &Vl[i]);
        Summation2(v[i], vl[i], Vh[i], Vl[i], &v[i], &vl[i]);

        // r[i] = r[i] + dt*(k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i])/6.0; 
        Summation2(krh2[i], krl2[i],  krh3[i], krl3[i], &inp[i].xh55[i], &inp[i].xl55[i]);
        Mult(2.0, 0.0, inp[i].xh55[i],  inp[i].xl55[i], &inp[i].xh66[i], &inp[i].xl66[i]);
        Summation2(krh1[i], krl1[i],  krh4[i], krl4[i], &inp[i].xh77[i], &inp[i].xl77[i]);
        Summation2(inp[i].xh66[i], inp[i].xl66[i], inp[i].xh77[i], inp[i].xl77[i], &inp[i].xh88[i], &inp[i].xl88[i]);
        Division(inp[i].xh88[i],  inp[i].xl88[i], 6.0, 0.0, &krfh[i], &krfl[i]);
        Mult(krfh[i], krfl[i], dt, 0.0, &Rh[i], &Rl[i]);
        Summation2(r[i],rl[i], Rh[i], Rl[i], &r[i], &rl[i]);  } }

    int main() {
      double t, dt, T, EKP[2], ANGM[2], a[2], e[2]; 

      //  Rh[3]={Xh,Yh,Zh} Rl[3]={Xl,Yl,Zl}  Vh[3]={Vxh,Vyh,Vzh} Vl[3]={Vxl,Vyl,Vzl} 
      // DATA: NASA HORIZON https://ssd.jpl.nasa.gov/horizons/app.html#/
      // JUPITER INITIAL VECTORS 2018-Jun-05 

       double  r[]= {-3.460167504309613E+00, -4.149454064629457E+00, 9.465721330038770E-02}; 
       double rl[]= {-0.0000000000000000486257802,  -0.0000000000000001269077053, -0.0000000000000004757454799e-2}; 
       double  v[]= {5.709741990408655E-03, -4.481465873394258E-03, -1.091471606521913E-04};
       double vl[]= {-0.0000000000000002100430652e-3, 0.0000000000000004062370784e-3, 0.0000000000000000152006560e-4};
    
    //  energies(r, rl, v, vl, EKP, ANGM);
    //  semimajor_eccentricity(r, rl, v, vl, a, e);  
    //  printf("%.20le %.20le %.20le %.20le %.20le %.20le\n", a[0] , a[1], e[0], e[1], energy[0], energy[1]);

      T=  2*M_PI*1000;
      dt= 1e-2;
      // error measure with long double 
      long double  A=  5.20273584355200957629L;       // semi-major axis lenght
      long double EC=  4.88056797545033895064e-02L;   // eccentricity 
      long double  E= -1.85032942351832700980e-01L;   // energy of system
      long double  H=  3.92090843714996407606e-02L;   // angular momentum of system

        FILE *kk, *pp;
        kk=fopen("RK4jupiter_accurate_pos.txt","w");
        pp=fopen("RK4jupiter_accurate_a_e_energy.txt","w");
        t=0.0 ;   
        int Count=0.0;
        while ( t <= T )  {    
        RK4(r, rl, v, vl, dt);                        
        energies(r, rl, v, vl, EKP, ANGM);              // energies
        semimajor_eccentricity(r, rl, v, vl, a, e);     // orbital elements
        long double energy= (EKP[0]-E+EKP[1])/E;
        long double angmomentum = (ANGM[0]-H+ANGM[1])/H;
        fprintf(kk,"%d %.20lf %.20lf %.20lf %.20lf %.20lf %.20lf \n", Count, r[0], rl[0], r[1], rl[1], r[2], rl[2]);
        fprintf(pp,"%d %.15Le %.15Le %.15Le %.15Le\n", Count, (a[0]-A+a[1])/A, (e[0]-EC+e[1])/EC, energy, angmomentum);
                                                                    
        t +=dt;
        Count++;
        }
        fclose(kk);
        fclose(pp);
        return (0);  }
