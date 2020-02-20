#include <math.h>
#include "nrutil.c"
#include "nrutil.h"

#define ITMAX 200
#define EPS 1.0e-10


//Here ITMAX is the maximum allowed number of iterations, while EPS is a small number to
//rectify the special case of converging to exactly zero function value.
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);
void frprmn(float p[], int n, float ftol, int *iter, float *fret,
            float (*func)(float [], atom *atm_Helix, int int_NumerOfAtoms), 
            void (*dfunc)(float [], float [], atom *atm_Helix, int
                    int_NumerOfAtoms), atom *atm_Helix, int int_NumerOfAtoms, void (*NormalizeZetas)(float p[]))
//Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a
//function func, using its gradient as calculated by a routine dfunc. The convergence tolerance
//on the function value is input as ftol. Returned quantities are p (the location of the minimum),
//iter (the number of iterations that were performed), and fret (the minimum value of the
//function). The routine linmin is called to perform line minimizations.
{
  void dlinmin(float p[], float xi[], int n, float *fret,
              float (*func)(float [], atom *atm_Helix, int int_NumerOfAtoms), 
            void  (*dfunc)(float [], float [], atom *atm_Helix, int
                    int_NumerOfAtoms), atom *atm_Helix, int int_NumerOfAtoms);
  
  
  
  int j,its;
  float gg,gam,fp,dgg;
  float *g,*h,*xi;
  g=vector(1,n);
  h=vector(1,n);
  xi=vector(1,n);
  
  NormalizeZetas(p);
  fp=(*func)(p, atm_Helix, int_NumerOfAtoms); //Initializations.
  
  (*dfunc)(p,xi,atm_Helix, int_NumerOfAtoms);
  //printf ("\n first eval: zetas  %f %f %f %f% f %fgives %f, \n with grad %f %f %f %f% f %f\n\n", p[1],p[2], p[3], p[4], p[5], p[6],fp, xi[1],xi[2], xi[3], xi[4],xi[5],xi[6]);
  for (j=1;j<=n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=ITMAX;its++) { //Loop over iterations.
    NormalizeZetas(p);
    *iter=its;
    //printf("\n %d, ",its);
    //printf ("\n dlinmin: zetas  %f %f %f %f% f %f gives %f, \n with grad %f %f %f %f% f %f\n\n", p[1],p[2], p[3], p[4], p[5], p[6],fp, xi[1],xi[2], xi[3], xi[4],xi[5],xi[6]);
    
    dlinmin(p, xi, n, fret, func, dfunc, atm_Helix, int_NumerOfAtoms); //Next statement is the normal return:
    
    //printf ("\n dlinmin: zetas  %f %f %f %f% f %f gives %f, \n with grad %f %f %f %f% f %f\n\n", p[1],p[2], p[3], p[4], p[5], p[6],fp, xi[1],xi[2], xi[3], xi[4],xi[5],xi[6]);
    
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      //printf("2.0*fabs(*fret (%f)-fp (%f)) <= ftol*(fabs(*fret)+fabs(fp)+EPS(%f))", *fret, fp,EPS); 
      FREEALL
      return;
    }
    fp= *fret;
    (*dfunc)(p,xi,atm_Helix, int_NumerOfAtoms);
    
    //printf ("\n dfunc: zetas  %f %f %f %f% f %f gives %f, \n with grad %f %f %f %f% f %f\n\n", p[1],p[2], p[3], p[4], p[5], p[6],fp, xi[1],xi[2], xi[3], xi[4],xi[5],xi[6]);
    
    dgg=gg=0.0;
    for (j=1;j<=n;j++) {
      gg += g[j]*g[j];
      // dgg += xi[j]*xi[j];  //This statement for Fletcher-Reeves.
      dgg += (xi[j]+g[j])*xi[j]; //This statement for Polak-Ribiere.
    }
    if (gg == 0.0) { //Unlikely. If gradient is exactly zero then
    //FREEALL we are already done.
      return;
    }
    gam=dgg/gg;
    for (j=1;j<=n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }//10.3 One-Dimensional Search with First Derivatives 407
  }
  printf("**********************\n Too many iterations in frprmn \n ");
}



/* dlinmin *********************************************************************/

/*Note on Line Minimization Using Derivatives
Kindly reread the last part of §10.5. We here want to do the same thing, but
using derivative information in performing the line minimization.
The modified version of linmin, called dlinmin, and its required companion
routine df1dim follow:*/

#define TOL 2.0e-4 //Tolerance passed to dbrent.
int ncom; //Global variables communicate with df1dim.
float *pcom,*xicom,(*nrfunc)(float [], atom *atm_Helix, int int_NumerOfAtoms);
void (*nrdfun)(float [], float [], atom *atm_Helix, int int_NumerOfAtoms);
void dlinmin(float p[], float xi[], int n, float *fret,
              float (*func)(float [], atom *atm_Helix, int int_NumerOfAtoms), 
            void  (*dfunc)(float [], float [], atom *atm_Helix, int
                    int_NumerOfAtoms), atom *atm_Helix, int int_NumerOfAtoms)
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p. This is actually all accomplished by calling the
routines mnbrak and dbrent.*/
{
  float dbrent(float ax, float bx, float cx,
        float (*f)(float, atom *atm_Helix, int int_NumerOfAtoms), float (*df)(float, atom *atm_Helix, int int_NumerOfAtoms), float tol, float *xmin, atom *atm_Helix, int int_NumerOfAtoms);
  float f1dim(float x, atom *atm_Helix, int int_NumerOfAtoms);
  float df1dim(float x, atom *atm_Helix, int int_NumerOfAtoms);
  void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
              float *fc, float (*func)(float, atom *atm_Helix, int int_NumerOfAtoms), atom *atm_Helix, int int_NumerOfAtoms);
  
  int j;
  float xx,xmin,fx,fb,fa,bx,ax;
  ncom=n; // Define the global variables.
  
  pcom = vector(1,n);
  xicom = vector(1,n);
  nrfunc=func;
  nrdfun=dfunc;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  
  ax=0.0; //Initial guess for brackets.
  xx=1.0;
  //printf("%f,%f,%f,%f\n", p[0],p[1],xi[0],xi[1]);
  
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim,atm_Helix, int_NumerOfAtoms);
  
  //printf("\n in dlinmin: ax = %f, xx = %f,bx = %f, fa=%f, fx=%f,fb=%f \n",ax,xx,bx,fa,fx,fb);
  
  *fret=dbrent(ax,xx,bx,*f1dim,df1dim,TOL,&xmin,atm_Helix, int_NumerOfAtoms);
  
  //printf("\n in dlinmin: ax = %f, xx = %f,bx = %f, fret=%f, xmin=%f (post dbrent)\n",ax,xx,bx,*fret,xmin);
  
  for (j=1;j<=n;j++) { //Construct the vector results to return.
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  //printf("%f,%f,%f,%f\n", p[0],p[1],xi[0],xi[1]);
  // now need to renormalize
  
  free_vector(xicom,1,n);
  free_vector(pcom,1,n);
}

/* df1dim *********************************************************************/

extern int ncom; //Defined in linmin.
extern float *pcom,*xicom,(*nrfunc)(float [], atom *atm_Helix, int int_NumerOfAtoms);
float f1dim(float x, atom *atm_Helix, int int_NumerOfAtoms)
//Must accompany linmin.
{
  int j;
  float f,*xt;
  xt=vector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt,atm_Helix, int_NumerOfAtoms);
  free_vector(xt,1,ncom);
  return f;
}

//#include "nrutil.h"
extern int ncom; //Defined in dlinmin.
extern float *pcom,*xicom,(*nrfunc)(float [], atom *atm_Helix, int int_NumerOfAtoms);
extern void (*nrdfun)(float [], float [], atom *atm_Helix, int int_NumerOfAtoms);
float df1dim(float x, atom *atm_Helix, int int_NumerOfAtoms)
{
  int j;
  float df1=0.0;
  float *xt,*df;
  xt=vector(1,ncom);
  df=vector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  (*nrdfun)(xt,df,atm_Helix, int_NumerOfAtoms);
  for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
  free_vector(df,1,ncom);
  free_vector(xt,1,ncom);
  return df1;
}

/* dbrent. Calls: func and dfunc *********************************************/ 
/* all taken from numerical recipies in C */
#include <math.h>
//#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
float dbrent(float ax, float bx, float cx, float (*f)(float, atom *atm_Helix, int int_NumerOfAtoms),
      float (*df)(float, atom *atm_Helix, int int_NumerOfAtoms), float tol, float *xmin, atom *atm_Helix, int int_NumerOfAtoms)
/*Given a function f and its derivative function df, and given a bracketing triplet of abscissas ax,
bx, cx [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)],
this routine isolates the minimum to a fractional precision of about tol using a modification of
Brent’s method that uses derivatives. The abscissa of the minimum is returned as xmin, and

the minimum function value is returned as dbrent, the returned function value.
10.3 One-Dimensional Search with First Derivatives 407
*/
{
  int iter,ok1,ok2; //Will be used as flags for whether profloat
  float a,b,d,d1,d2,du,dv,dw,dx,e=0.0; //posed steps are acceptable or not.
  float fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  //Comments following will point out only differences from the routine brent. Read that
  //routine first.
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x,atm_Helix, int_NumerOfAtoms);
  dw=dv=dx=(*df)(x,atm_Helix, int_NumerOfAtoms); //All our housekeeping chores are doubled
  //by the necessity of moving
  //derivative values around as well
  //as function values.
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      d1=2.0*(b-a); //Initialize these d’s to an out-of-bracket
      d2=d1; //value.
      if (dw != dx) d1=(w-x)*dx/(dx-dw); //Secant method with one point.
      if (dv != dx) d2=(v-x)*dx/(dx-dv); //And the other.
      //Which of these two estimates of d shall we take? We will insist that they be within
      //the bracket, and on the side pointed to by the derivative at x:
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e; //Movement on the step before last.
      e=d;
      if (ok1 || ok2) { /*Take only an acceptable d, and if
        both are acceptable, then take
        the smallest one.*/
        if (ok1 && ok2)
          d=(fabs(d1) < fabs(d2) ? d1 : d2);
        else if (ok1)
          d=d1;
        else
          d=d2;
        if (fabs(d) <= fabs(0.5*olde)) {
          u=x+d;
          if (u-a < tol2 || b-u < tol2)
            d=SIGN(tol1,xm-x);
        } else {  //Bisect, not golden section.
          d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
          //Decide which segment by the sign of the derivative.
        }
      } else {
        d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= tol1) {
      u=x+d;
      fu=(*f)(u,atm_Helix, int_NumerOfAtoms);
    } else {
      u=x+SIGN(tol1,d);
      fu=(*f)(u,atm_Helix, int_NumerOfAtoms);
      if (fu > fx) { //If the minimum step in the downhill
        //direction takes us uphill, then we are done.
        *xmin=x;
        return fx;
      }
    } 
    du=(*df)(u,atm_Helix, int_NumerOfAtoms); //Now all the housekeeping, sigh.
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw)
      MOV3(w,fw,dw, x,fx,dx)
      MOV3(x,fx,dx, u,fu,du)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        MOV3(v,fv,dv, w,fw,dw)
        MOV3(w,fw,dw, u,fu,du)
      } else if (fu < fv || v == x || v == w) {
        MOV3(v,fv,dv, u,fu,du)
      }
    }
  }
  printf("Too many iterations in routine dbrent");
  return 0.0; //Never get here.
}




/* mnbrak *********************************************************************/

#include <math.h>
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
/*Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa, fb, and fc.*/
/*Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the
maximum magnification allowed for a parabolic-fit step.*/

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,float (*func)(float, atom *atm_Helix, int int_NumerOfAtoms), atom *atm_Helix, int int_NumerOfAtoms)
{
    float ulim,u,r,q,fu,dum;
    *fa=(*func)(*ax,atm_Helix, int_NumerOfAtoms);
    *fb=(*func)(*bx,atm_Helix, int_NumerOfAtoms);
    if (*fb > *fa) {
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx,atm_Helix, int_NumerOfAtoms);
    while (*fb > *fc) {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
    /*Switch roles of a and b so that we can go
    downhill in the direction from a to b.
    First guess for c.
    Keep returning here until we bracket.
    Compute u by parabolic extrapolation from
    a, b, c. TINY is used to prevent any pos-
    sible division by zero.
    10.1 Golden Section Search in One Dimension*/
   
        (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        //We won’t go farther than this. Test various possibilities:
        if ((*bx-u)*(u-*cx) > 0.0) {
            //Parabolic u is between b and c: try it.
            fu=(*func)(u,atm_Helix, int_NumerOfAtoms);
            if (fu < *fc) {
            //Got a minimum between b and c.
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            } else if (fu > *fb) {
            //Got a minimum between between a and u.
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx);
            //Parabolic fit was no use. Use default magnification. 
            fu=(*func)(u,atm_Helix, int_NumerOfAtoms);

        } else if ((*cx-u)*(u-ulim) > 0.0) {
        //Parabolic fit is between c and its allowed limit.
            fu=(*func)(u,atm_Helix, int_NumerOfAtoms);

            if (fu < *fc) {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                SHFT(*fb,*fc,fu,(*func)(u,atm_Helix, int_NumerOfAtoms))
            }
        } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
        //Limit parabolic u to maximum allowed value. 
            u=ulim;
            fu=(*func)(u,atm_Helix, int_NumerOfAtoms);
        } else {
        //Reject parabolic u, use default magnification.
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u,atm_Helix, int_NumerOfAtoms);
        }
        SHFT(*ax,*bx,*cx,u)
        //Eliminate oldest point and continue.
        SHFT(*fa,*fb,*fc,fu)
    }
}  
