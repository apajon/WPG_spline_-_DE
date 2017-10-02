#include "mex.h"
#include "cfsqp.h"


//Variables globales utilisées par les sous expressions
char obj_name[20],constr_name[20],gradob_name[20],gradcn_name[20];
mxArray *plhs_obj[1], *prhs_obj[3];
mxArray *plhs_constr[1], *prhs_constr[3];
mxArray *plhs_gradob[1], *prhs_gradob[3];
mxArray *plhs_gradcn[1], *prhs_gradcn[3];
double *save_f, *save_g;
double **save_gradf, **save_gradg;
double *ptr;
int M, N, length_f, length_g, stop_f_neg;
int nineqn, nineq, neq, neqn;
int rank_calculated_fi, rank_calculated_gi, rank_calculated_gradfi, rank_calculated_gradgi, call_matlab;
  
// Fonction objectif
void obj(int nparam, int j, double *x, double *fj, void *cd){

int i;

mexPrintf("Appel objective function %d\n",j);

if (x_is_new==TRUE) {
   x_is_new=FALSE;
   call_matlab=1;
   rank_calculated_fi=0;
   rank_calculated_gi=0;
   rank_calculated_gradfi=0;
   rank_calculated_gradgi=0;
} else {
   if (rank_calculated_fi>=j) {
      fj[0]=save_f[j-1];
      call_matlab=0;
   } else {
      call_matlab=1;
   }
}

if (call_matlab) {
   mxGetPr(prhs_obj[0])[0] = j;
   ptr = mxGetPr(prhs_obj[1]);
   for (i=0; i<nparam; i++)
   {
      ptr[i] = x[i];
   }
   prhs_obj[2] = (mxArray *)cd;
   mexCallMATLAB(1,plhs_obj,3,prhs_obj,obj_name);
   
   N=mxGetN(plhs_obj[0]);
   M=mxGetM(plhs_obj[0]);
   
   if (M==1) {
      M=N;
      N=1;
   } else if (N!=1) {
      mexErrMsgTxt("f is not a vector");
   }
   
   ptr = mxGetPr(plhs_obj[0]);

   if (M>length_f) { //if the function give more than expected
      M=length_f;
   }
   
   if (M>=j) {
      for (i=rank_calculated_fi; i<M; i++)
         save_f[i]=ptr[i];
      rank_calculated_fi=M;
      fj[0] = save_f[j-1];
   } else if (M==1) {
      fj[0] = ptr[0];
   } else { //M==0 or 1<M<j
      mexErrMsgTxt("f is empty or wrong dimension");
   }
   
   //test the stop if the criteria is negative (used for the optimization of multi-level constraints)
   if ((stop_f_neg) & (rank_calculated_fi==length_f)) {
      i=0;
      while ((save_f[i]<=0) & (i<length_f-1)) {
         i=i+1;
      }
      if ((save_f[i]<=0) & (i==length_f-1)) {
         nstop=0;
      }
   }

}

}

// Fonction des contraintes
void constr(int nparam, int j, double *x, double *gj, void *cd){

int i,jn;

mexPrintf("Appel constraints function %d\n",j);

if (j<=nineqn) { //if j is a nonlinear constraint
   jn=j+nineq-nineqn;
} else if (j<=nineq) { //if j is for a linear constraint
   jn=j-nineqn;
} else if (j<=nineq+neqn) { //if j for a nonlinear equality constraint
   jn=j+neq-neqn;
} else if (j<=nineq+neq) { //if j is for q lineqr equality constraint
   jn=j-neqn;
}

if (x_is_new==TRUE) {
   x_is_new=FALSE;
   call_matlab=1;
   rank_calculated_fi=0;
   rank_calculated_gi=0;
   rank_calculated_gradfi=0;
   rank_calculated_gradgi=0;
} else {
   if (rank_calculated_gi>=jn) {
      gj[0]=save_g[jn-1];
      call_matlab=0;
   } else {
      call_matlab=1;
   }
}

if (call_matlab) {
   mxGetPr(prhs_constr[0])[0] = jn;
   ptr = mxGetPr(prhs_constr[1]);
   for (i=0; i<nparam; i++)
   {
      ptr[i] = x[i];
   }
   prhs_constr[2] = (mxArray *)cd;
   
   mexCallMATLAB(1,plhs_constr,3,prhs_constr,constr_name);
   
   N=mxGetN(plhs_constr[0]);
   M=mxGetM(plhs_constr[0]);
   
   if (M==1) {
      M=N;
      N=1;
   } else if (N!=1) {
      mexErrMsgTxt("g is not a vector");
   }
   
   ptr = mxGetPr(plhs_constr[0]);
   
   if (M>length_g) { //if the function give more than expected
      M=length_g;
   }
   
   if (M>=jn) {
      for (i=rank_calculated_gi; i<M; i++)
         save_g[i]=ptr[i];
      rank_calculated_gi=M;
      gj[0] = save_g[jn-1];
   } else if (M==1) {
      gj[0] = ptr[0];
   } else { //M==0 or 1<M<jn
      mexErrMsgTxt("g is empty or wrong dimension");
   }
}

}


// Fonction du gradient de l'objectif
void gradob(int nparam, int j, double *x, double *gradfj, void (* dummy)(), void *cd){
	
int i,k,nb_ligne;

mexPrintf("Appel objective gradient function %d\n",j);

if (x_is_new==TRUE) {
   x_is_new=FALSE;
   call_matlab=1;
   rank_calculated_fi=0;
   rank_calculated_gi=0;
   rank_calculated_gradfi=0;
   rank_calculated_gradgi=0;
} else {
   if (rank_calculated_gradfi>=j) {
      for (k=0; k<nparam; k++) {
         gradfj[k]=save_gradf[j-1][k];
      }
      call_matlab=0;
   } else {
      call_matlab=1;
   }
}

if (call_matlab) {
   mxGetPr(prhs_gradob[0])[0] = j;
   ptr = mxGetPr(prhs_gradob[1]);
   for (i=0; i<nparam; i++)
   {
      ptr[i] = x[i];
   }
   prhs_gradob[2] = cd;
   
   mexCallMATLAB(1,plhs_gradob,3,prhs_gradob,gradob_name);
   
   N=mxGetN(plhs_gradob[0]);
   M=mxGetM(plhs_gradob[0]);
   nb_ligne=M;
   
   if (N!=nparam) {
      mexErrMsgTxt("gradf is not a matrix with nparam columns");
   }
   
   ptr = mxGetPr(plhs_gradob[0]);

   if (M>length_f) { //if the function give more than expected
      M=length_f;
   }
   
   if (M>=j) {
      for (i=rank_calculated_gradfi; i<M; i++) {
         for (k=0; k<nparam; k++) {
            save_gradf[i][k]=ptr[i+k*nb_ligne];
         }
      }
      rank_calculated_gradfi=M;
      for (k=0; k<nparam; k++) {
         gradfj[k] = save_gradf[j-1][k];
      }
   } else if (M==1) {
      for (k=0; k<nparam; k++) {
         gradfj[k] = ptr[k];
      }
   } else { //M==0 or 1<M<j
      mexErrMsgTxt("gradf is empty or wrong dimension");
   }
}

}


// Fonction du gradient des contraintes
void gradcn(int nparam, int j, double *x, double *gradgj, void (* dummy)(), void *cd){
	
int i,k,jn,nb_ligne;

mexPrintf("Appel constraints gradient function %d\n",j);

if (j<=nineqn) { //if j is a nonlinear constraint
   jn=j+nineq-nineqn;
} else if (j<=nineq) { //if j is for a linear constraint
   jn=j-nineqn;
} else if (j<=nineq+neqn) { //if j for a nonlinear equality constraint
   jn=j+neq-neqn;
} else if (j<=nineq+neq) { //if j is for q lineqr equality constraint
   jn=j-neqn;
}

if (x_is_new==TRUE) {
   x_is_new=FALSE;
   call_matlab=1;
   rank_calculated_fi=0;
   rank_calculated_gi=0;
   rank_calculated_gradfi=0;
   rank_calculated_gradgi=0;
} else {
   if (rank_calculated_gradgi>=jn) {
      for (k=0; k<nparam; k++) {
         gradgj[k]=save_gradg[jn-1][k];
      }
      call_matlab=0;
   } else {
      call_matlab=1;
   }
}

if (call_matlab) {
   mxGetPr(prhs_gradcn[0])[0] = jn;
   ptr = mxGetPr(prhs_gradcn[1]);
   for (i=0; i<nparam; i++)
   {
      ptr[i] = x[i];
   }
   prhs_gradcn[2] = cd;
   
   mexCallMATLAB(1,plhs_gradcn,3,prhs_gradcn,gradcn_name);
   
   N=mxGetN(plhs_gradcn[0]);
   M=mxGetM(plhs_gradcn[0]);
   nb_ligne=M;
   
   if (N!=nparam) {
      mexErrMsgTxt("gradg is not a matrix with nparam columns");
   }
   
   ptr = mxGetPr(plhs_gradcn[0]);

   if (M>length_g) { //if the function give more than expected
      M=length_g;
   }
   
   if (M>=jn) {
      for (i=rank_calculated_gradgi; i<M; i++) {
         for (k=0; k<nparam; k++) {
            save_gradg[i][k]=ptr[i+k*nb_ligne];
         }
      }
      rank_calculated_gradgi=M;
      for (k=0; k<nparam; k++) {
         gradgj[k] = save_gradg[jn-1][k];
      }
   } else if (M==1) {
      for (k=0; k<nparam; k++) {
         gradgj[k] = ptr[k];
      }
   } else { //M==0 or 1<M<jn
      mexErrMsgTxt("gradg is empty or wrong dimension");
   }
}

}



void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
// Variables de la fonction fsqp
int     nparam,nf,nfsr,ncsrl,ncsrn,mode,
        iprint,miter,*mesh_pts,inform[1];
double  bigbnd,eps,epseqn,udelta;
double  *bl,*bu,*x,*f,*g,*lambda;
void    *cd;
// Variables propres à mexfsqp
  int i, Np, length_lambda;
  int fd_gradob, fd_gradcn;
       
  /* Examine input (right-hand-side) arguments. */
  if (nrhs < 25) {
    mexErrMsgTxt("mexfsqp requires 25 input arguments.");
  }
  
  /* Examine output (left-hand-side) arguments. */
  if (nlhs != 5) {
    mexErrMsgTxt("mexfsqp requires 5 output arguments.");
  }
  
  //conversion des parametres d'entree
  ptr = mxGetPr(prhs[0]);
  nparam = ptr[0];
  
  //mexPrintf("nparam=%d\n",nparam);
  
  ptr = mxGetPr(prhs[1]);
  nf = ptr[0];
  
  ////mexPrintf("nf=%d\n",nf);
  
  ptr = mxGetPr(prhs[2]);
  nfsr = ptr[0];
  
  ////mexPrintf("nfsr=%d\n",nfsr);
  
  ptr = mxGetPr(prhs[3]);
  nineqn = ptr[0];
  
  //mexPrintf("nineqn=%d\n",nineqn);
  
  ptr = mxGetPr(prhs[4]);
  nineq = ptr[0];
  
  //mexPrintf("nineq=%d\n",nineq);
  
  ptr = mxGetPr(prhs[5]);
  neqn = ptr[0];
  
  //mexPrintf("neqn=%d\n",neqn);
  
  ptr = mxGetPr(prhs[6]);
  neq = ptr[0];
  
  //mexPrintf("neq=%d\n",neq);
  
  ptr = mxGetPr(prhs[7]);
  ncsrl = ptr[0];
  
  //mexPrintf("ncsrl=%d\n",ncsrl);
  
  ptr = mxGetPr(prhs[8]);
  ncsrn = ptr[0];
  
  //mexPrintf("ncsrn=%d\n",ncsrn);
  
  N=mxGetN(prhs[9]);
  M=mxGetM(prhs[9]);
  Np=nfsr+ncsrl+ncsrn;
  //mexPrintf("N=%d \nM=%d \nNp=%d\n",N,M,Np);
  if ( ( (N!=Np) & !((Np==0) & (N==1)) ) | (M!=1) )
     mexErrMsgTxt("mesh_pts is not of size max{1,nfsr+ncsrl+ncsrn}");
  if (Np==0) Np=1;
  if ((mesh_pts = (int *) malloc(Np*sizeof(int))) == NULL)
     mexErrMsgTxt("Error on malloc for mesh_pts\n");
  ptr = mxGetPr(prhs[9]);
  for (i=0; i<N; i++)
  {
      mesh_pts[i] = ptr[i];
      //mexPrintf("mesh_pts[%d]=%d\n",i,mesh_pts[i]);
  }
  
  ptr = mxGetPr(prhs[10]);
  mode = ptr[0];
  
  //mexPrintf("mode=%d\n",mode);
  
  ptr = mxGetPr(prhs[11]);
  iprint = ptr[0];
  
  //mexPrintf("iprint=%d\n",iprint);
  
  ptr = mxGetPr(prhs[12]);
  miter = ptr[0];
  
  //mexPrintf("miter=%d\n",miter);
  
  ptr = mxGetPr(prhs[13]);
  bigbnd = ptr[0];
  
  //mexPrintf("bigbnd=%22.14e\n",bigbnd);
  
  ptr = mxGetPr(prhs[14]);
  eps = ptr[0];
  
  //mexPrintf("eps=%22.14e\n",eps);
  
  ptr = mxGetPr(prhs[15]);
  epseqn = ptr[0];
  
  //mexPrintf("epseqn=%22.14e\n",epseqn);
  
  ptr = mxGetPr(prhs[16]);
  udelta = ptr[0];
  
  //mexPrintf("udelta=%22.14e\n",udelta);
  
  N=mxGetN(prhs[17]);
  M=mxGetM(prhs[17]);
  if ((N!=nparam) | (M!=1))
     mexErrMsgTxt("bl is not of size nparam");
  if ((bl = (double *) malloc(N*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for bl\n");
  ptr = mxGetPr(prhs[17]);
  for (i=0; i<N; i++)
  {
      bl[i] = ptr[i];
      //mexPrintf("bl[%d]=%22.14e\n",i,bl[i]);
  }
  
  N=mxGetN(prhs[18]);
  M=mxGetM(prhs[18]);
  if ((N!=nparam) | (M!=1))
     mexErrMsgTxt("bu is not of size nparam");
  if ((bu = (double *) malloc(N*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for bl\n");
  ptr = mxGetPr(prhs[18]);
  for (i=0; i<N; i++)
  {
      bu[i] = ptr[i];
      //mexPrintf("bu[%d]=%22.14e\n",i,bu[i]);
  }
  
  N=mxGetN(prhs[19]);
  M=mxGetM(prhs[19]);
  if ((N!=nparam) | (M!=1))
     mexErrMsgTxt("x is not of size nparam");
  if ((x = (double *) malloc(N*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for bl\n");
  ptr = mxGetPr(prhs[19]);
  for (i=0; i<N; i++)
  {
      x[i] = ptr[i];
      //mexPrintf("x[%d]=%22.14e\n",i,x[i]);
  }

  N=mxGetN(prhs[20]);
  M=mxGetM(prhs[20]);
  //mexPrintf("N=%d \nM=%d\n",N,M);
  if ((N>19) | (M!=1))
     mexErrMsgTxt("Name of objective function must be less than 19 caracters");
  mxGetString(prhs[20],obj_name,N+1);
  i=0;
  while ((obj_name[i]!=NULL) & (i<20)) {
     //mexPrintf("%c",obj_name[i]);
     i=i+1;
  }
  //mexPrintf("\n");
    
  N=mxGetN(prhs[21]);
  M=mxGetM(prhs[21]);
  if ((N>19) | (M!=1))
     mexErrMsgTxt("Name of constraints function must be less than 19 caracters");
  mxGetString(prhs[21],constr_name,N+1);
  i=0;
  while ((constr_name[i]!=NULL) & (i<20)) {
     //mexPrintf("%c",constr_name[i]);
     i=i+1;
  }
  //mexPrintf("\n");
    
  N=mxGetN(prhs[22]);
  M=mxGetM(prhs[22]);
  if ((N>19) | (M!=1))
     mexErrMsgTxt("Name of objective gradient function must be less than 19 caracters");
  mxGetString(prhs[22],gradob_name,N+1);
  i=0;
  while ((gradob_name[i]!=NULL) & (i<20)) {
     //mexPrintf("%c",gradob_name[i]);
     i=i+1;
  }
  //mexPrintf("\n");
    
  N=mxGetN(prhs[23]);
  M=mxGetM(prhs[23]);
  if ((N>19) | (M!=1))
     mexErrMsgTxt("Name of constraints gradient function must be less than 19 caracters");
  mxGetString(prhs[23],gradcn_name,N+1);
  i=0;
  while ((gradcn_name[i]!=NULL) & (i<20)) {
     //mexPrintf("%c",gradcn_name[i]);
     i=i+1;
  }
  //mexPrintf("\n");
  
  cd = prhs[24];
  
  //ptr = mxGetPr(prhs[25]);
  //stop_f_neg = ptr[0];


  //allocate memmory for f, g and lambda
  length_f=nf-nfsr;
  for (i=0; i<nfsr; i++)
  {
      length_f += mesh_pts[i];
  }
  if (length_f==0) length_f=1;
  if ((f = (double *) malloc(length_f*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for f\n");

  length_g=nineq+neq-ncsrl-ncsrn;
  for (i=nfsr; i<ncsrl+ncsrn; i++)
  {
      length_g += mesh_pts[i];
  }
  if (length_g==0) length_g=1;
  if ((g = (double *) malloc(length_g*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for g\n");

  length_lambda=nparam+length_f+length_g;
  if ((lambda = (double *) malloc(length_lambda*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for lambda\n");
  
  //allocate memmory for plhs and prhs of obj, constr, gradob, gradcn
  prhs_obj[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  prhs_obj[1]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  plhs_obj[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  
  prhs_constr[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  prhs_constr[1]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  plhs_constr[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  
  prhs_gradob[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  prhs_gradob[1]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  plhs_gradob[0]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  
  prhs_gradcn[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  prhs_gradcn[1]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  plhs_gradcn[0]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  
  
  //allocate memmory for the save of data at the same x
  if ((save_f=(double *) malloc(length_f*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for save_f\n");
  if ((save_g=(double *) malloc(length_g*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for save_g\n");
  if ((save_gradf=(double **) malloc(length_f*sizeof(double *))) == NULL)
     mexErrMsgTxt("Error on malloc for save_gradf\n");
  if ((save_gradf[0]=(double *) malloc(nparam*length_f*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for save_gradf\n");  
  for (i=1; i<length_f; i++)
     save_gradf[i]=save_gradf[i-1]+nparam;
  if ((save_gradg=(double **) malloc(length_g*sizeof(double *))) == NULL)
     mexErrMsgTxt("Error on malloc for save_gradg\n");
  if ((save_gradg[0]=(double *) malloc(nparam*length_g*sizeof(double))) == NULL)
     mexErrMsgTxt("Error on malloc for save_gradg\n");  
  for (i=1; i<length_g; i++)
     save_gradg[i]=save_gradg[i-1]+nparam;
  
  
  
  // Call to fsqp
  // -----------------------------------------------
  
  // initialization
  x_is_new=TRUE;
  rank_calculated_fi=0;
  rank_calculated_gi=0;
  rank_calculated_gradfi=0;
  rank_calculated_gradgi=0;


  // si les fonctions de gradient ne sont pas fournies, calcul par difference finie
  fd_gradob=(gradob_name[0]==NULL);
  fd_gradcn=(gradcn_name[0]==NULL);
  
  if (!fd_gradob & !fd_gradcn)
    cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
      mode,iprint,miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,
      f,g,lambda,obj,constr,gradob,gradcn,cd);
  if (fd_gradob & !fd_gradcn)
    cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
      mode,iprint,miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,
      f,g,lambda,obj,constr,grobfd,gradcn,cd);
  if (!fd_gradob & fd_gradcn)
    cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
      mode,iprint,miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,
      f,g,lambda,obj,constr,gradob,grcnfd,cd);
  if (fd_gradob & fd_gradcn)
    cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
      mode,iprint,miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,
      f,g,lambda,obj,constr,grobfd,grcnfd,cd);
      
      
  // Conversion des parametres de sortie
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  mxGetPr(plhs[0])[0]=(double)inform[0];
  
  plhs[1]=mxCreateDoubleMatrix(1,nparam,mxREAL);
  ptr = mxGetPr(plhs[1]);
  for (i=0; i<nparam; i++)
  {
      ptr[i] = x[i];
  }
  
  plhs[2]=mxCreateDoubleMatrix(1,length_f,mxREAL);
  ptr = mxGetPr(plhs[2]);
  for (i=0; i<length_f; i++)
  {
      ptr[i] = f[i];
  }
  
  plhs[3]=mxCreateDoubleMatrix(1,length_g,mxREAL);
  ptr = mxGetPr(plhs[3]);
  for (i=0; i<length_g; i++)
  {
      ptr[i] = g[i];
  }
  
  plhs[4]=mxCreateDoubleMatrix(1,length_lambda,mxREAL);
  ptr = mxGetPr(plhs[4]);
  for (i=0; i<length_lambda; i++)
  {
      ptr[i] = lambda[i];
  }
  
  //free memory
  free(mesh_pts);
  free(bl);
  free(bu);
  free(x);
  free(f);
  free(g);
  free(lambda);
  free(save_f);
  free(save_g);
  free(save_gradf[0]);
  free(save_gradg[0]);
  free(save_gradf);
  free(save_gradg);
  mxDestroyArray(prhs_obj[0]);
  mxDestroyArray(prhs_obj[1]);
  mxDestroyArray(plhs_obj[0]);
  mxDestroyArray(prhs_constr[0]);
  mxDestroyArray(prhs_constr[1]);
  mxDestroyArray(plhs_constr[0]);
  mxDestroyArray(prhs_gradob[0]);
  mxDestroyArray(prhs_gradob[1]);
  mxDestroyArray(plhs_gradob[0]);
  mxDestroyArray(prhs_gradcn[0]);
  mxDestroyArray(prhs_gradcn[1]);
  mxDestroyArray(plhs_gradcn[0]);
}
