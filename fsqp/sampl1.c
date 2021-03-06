/***********************************************************/
/*    Sample 1 from CFSQP distribution			   */
/*    Problem 32 from Hock/Schittkowski 1981               */
/***********************************************************/

#include "cfsqp.h"

void obj32();
void cntr32();
void grob32();
void grcn32();

int main() {
   int nparam,nf,nineq,neq,mode,iprint,miter,neqn,nineqn,
       ncsrl,ncsrn,nfsr,mesh_pts[1],inform;
   double bigbnd,eps,epsneq,udelta;
   double *x,*bl,*bu,*f,*g,*lambda;
   void *cd;

   mode=100;
   iprint=1;
   miter=500;  
   bigbnd=1.e10;
   eps=1.e-8;
   epsneq=0.e0;
   udelta=0.e0;
   nparam=3;
   nf=1;
   neqn=0;
   nineqn=1;
   nineq=1;
   neq=1;
   ncsrl=ncsrn=nfsr=mesh_pts[0]=0;
   bl=(double *)calloc(nparam,sizeof(double));
   bu=(double *)calloc(nparam,sizeof(double));
   x=(double *)calloc(nparam,sizeof(double));
   f=(double *)calloc(nf,sizeof(double));
   g=(double *)calloc(nineq+neq,sizeof(double));
   lambda=(double *)calloc(nineq+neq+nf+nparam,sizeof(double));
   
   bl[0]=bl[1]=bl[2]=0.e0;
   bu[0]=bu[1]=bu[2]=bigbnd;

   x[0]=0.1e0;
   x[1]=0.7e0;
   x[2]=0.2e0;

   cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
         mode,iprint,miter,&inform,bigbnd,eps,epsneq,udelta,bl,bu,x,
         f,g,lambda,obj32,cntr32,grob32,grcn32,cd);

   free(bl);
   free(bu);
   free(x);
   free(f);
   free(g);
   free(lambda);
   return 0;
}
   
void obj32(nparam,j,x,fj,cd)
int nparam,j;
double *x,*fj;
void *cd;
{
   *fj=pow((x[0]+3.e0*x[1]+x[2]),2.e0)+4.e0*pow((x[0]-x[1]),2.e0);
   return;
}

void grob32(nparam,j,x,gradfj,dummy,cd)
int nparam,j;
double *x,*gradfj;
void (* dummy)();
void *cd;
{
   double fa,fb;

   fa=2.e0*(x[0]+3.e0*x[1]+x[2]);
   fb=8.e0*(x[0]-x[1]);
   gradfj[0]=fa+fb;
   gradfj[1]=fa*3.e0-fb;
   gradfj[2]=fa;
   return;
}

void cntr32(nparam,j,x,gj,cd)
int nparam,j;
double *x,*gj;
void *cd;
{
   switch (j) {
      case 1:
	 *gj=pow(x[0],3.e0)-6.e0*x[1]-4.e0*x[2]+3.e0;
	 break;
      case 2:
	 *gj=1.e0-x[0]-x[1]-x[2];
	 break;
   }
   return;
}

void grcn32(nparam,j,x,gradgj,dummy,cd)
int nparam,j;
double *x,*gradgj;
void (* dummy)();
void *cd;
{
   switch (j) {
      case 1:
	 gradgj[0]=3.e0*x[0]*x[0];
	 gradgj[1]=-6.e0;
	 gradgj[2]=-4.e0;
	 break;
      case 2:
	 gradgj[0]=gradgj[1]=gradgj[2]=-1.e0;
	 break;
   }
   return;
}   
