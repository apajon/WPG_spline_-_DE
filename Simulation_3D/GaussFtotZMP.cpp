/*==========================================================
 * GaussFtotZMP.cpp - Gauss-Seidel contacts to obtain sole
 * orientation and displacement for a desired ZMP position
 * and force
 *==========================================================
*/

#include <iostream>
#include <math.h>
#include "mex.h"
#include "Gausstot.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <fstream>

using namespace std;
using namespace Eigen;

extern void _main();

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

static void Gauss(int contAngle,double fric,int m,Vector3d& displ_in,Vector3d& angleact_in,Vector6d& FtotZMPdes_in, const VectorXd& F3_in,const MatrixXd& Pfree_in, const MatrixXd& PabsOld_in, const MatrixXd& Cs_in,Vector3d& displini_in,double *displ_out,double *angleact_out,double *FtotZMPOut,double *Fc_mc3_out,double *PSurf3,double *ind_cont_out,double *ind_slip_out,double *Kcart_out)
{
  GaussSeidel *GS = new GaussSeidel;

  GS->setStaticDataSim(Cs_in,FtotZMPdes_in,fric,m);

  GS->setStaticDataStep(displ_in,angleact_in,F3_in,Pfree_in,PabsOld_in,displini_in);

  GS->run(contAngle);

  GS->set_data_output(displ_out,angleact_out,FtotZMPOut,Fc_mc3_out,PSurf3,ind_cont_out,ind_slip_out,Kcart_out); // Set data members to incoming

  delete(GS);
  flush(cout);
  return;
}

void mexFunction(int num_output, mxArray *output[], int num_input, const mxArray *input[])
{
  int contAngle;
  double fric;
  int m;                       
  double *displ_in;  
  double *displini_in;
  double *angleact_in;               
  double *FtotZMPdes_in;          
  double *F3_in;                    
  double *Pfree_in;                    
  double *PabsOld_in;                    
  double *Cs_in;
  double *displ_out;             
  double *angleact_out;             
  double *FtotZMP_out;          
  double *Fc_mc3_out;           
  double *PSurf3;           
  double *ind_cont_out;
  double *ind_slip_out;
  double *Kcart_out;
  double *B_omega_out;
  double *G_out;
  double *A_out;
  double *B_out;

  if (num_input != 11) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargin",
            "Gauss requires 10 input arguments.");
  } else if (num_output >= 13) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargout",
            "Gauss requires 13 output argument.");
  }

  contAngle = (int)mxGetScalar(input[0]);
  fric = mxGetScalar(input[1]);

  m = (int)mxGetScalar(input[2]);

  displ_in = mxGetPr(input[3]);
  Vector3d displ_in1 = Map<Vector3d>(displ_in);

  angleact_in = mxGetPr(input[4]);
  Vector3d angleact_in1 = Map<Vector3d>(angleact_in);

  FtotZMPdes_in = mxGetPr(input[5]);
  Vector6d FtotZMPdes_in1 = Map<Vector6d>(FtotZMPdes_in);

  size_t rows = mxGetM(input[6]);
  size_t cols = mxGetN(input[6]);
  F3_in = mxGetPr(input[6]);

  VectorXd F3_in1 = Map<VectorXd>(F3_in, rows);

  rows = mxGetM(input[7]);
  cols = mxGetN(input[7]);
  Pfree_in = mxGetPr(input[7]);
  MatrixXd Pfree_in1 = Map<MatrixXd>(Pfree_in, rows, cols);

  rows = mxGetM(input[8]);
  cols = mxGetN(input[8]);
  PabsOld_in = mxGetPr(input[8]);
  MatrixXd PabsOld_in1 = Map<MatrixXd>(PabsOld_in, rows, cols);

  rows = mxGetM(input[9]);
  cols = mxGetN(input[9]);
  Cs_in = mxGetPr(input[9]);
  MatrixXd Cs_in1 = Map<MatrixXd>(Cs_in, rows, cols);
  
  displini_in = mxGetPr(input[10]);
  Vector3d displini_in1 = Map<Vector3d>(displini_in);  

  output[0] = mxCreateDoubleMatrix(3,1,mxREAL);  //
  displ_out = mxGetPr(output[0]);

  output[1] = mxCreateDoubleMatrix(3,1,mxREAL);  //
  angleact_out = mxGetPr(output[1]);

  output[2] = mxCreateDoubleMatrix(6,1,mxREAL);  
  FtotZMP_out = mxGetPr(output[2]); 

  output[3] = mxCreateDoubleMatrix(3*m,1,mxREAL);  
  Fc_mc3_out = mxGetPr(output[3]);

  output[4] = mxCreateDoubleMatrix(3*m,1,mxREAL);
  PSurf3 = mxGetPr(output[4]);
  
  output[5] = mxCreateDoubleMatrix(m,1,mxREAL);
  ind_cont_out = mxGetPr(output[5]);

  output[6] = mxCreateDoubleMatrix(m,1,mxREAL);
  ind_slip_out = mxGetPr(output[6]);
  
  output[7] = mxCreateDoubleMatrix(6,6,mxREAL);
  Kcart_out = mxGetPr(output[7]);
  
//   output[8] = mxCreateDoubleMatrix(3*m*6,1,mxREAL);
//   B_omega_out = mxGetPr(output[8]);
// 
//   output[9] = mxCreateDoubleMatrix(3*m*6,1,mxREAL);
//   G_out = mxGetPr(output[9]);  
// 
//   output[10] = mxCreateDoubleMatrix(3*m*3*m,1,mxREAL);
//   A_out = mxGetPr(output[10]);  
// 
//   output[11] = mxCreateDoubleMatrix(3*m*3*m,1,mxREAL);
//   B_out = mxGetPr(output[11]); 
  
  Gauss(contAngle,fric,m,displ_in1,angleact_in1,FtotZMPdes_in1,F3_in1,Pfree_in1,PabsOld_in1,Cs_in1,displini_in1,displ_out,angleact_out,FtotZMP_out,Fc_mc3_out,PSurf3,ind_cont_out,ind_slip_out,Kcart_out);

  return;
}
