/*
 * GaussSeidel.h
 */

#include <Eigen/Dense>
#include <Eigen/StdVector>
// #ifndef _GAUSS_SEIDEL_H_
// #define _GAUSS_SEIDEL_H_
//#include "mex.h"

using namespace std;
using namespace Eigen;

class GaussSeidel {
    public:
	  typedef std::vector<Matrix3d,aligned_allocator<Matrix3d> > MatrixArray;
	  typedef Matrix<double,5,1> Vector5d;
	  typedef Matrix<double,6,1> Vector6d;
	  typedef Matrix<double,5,5> Matrix5d;

      GaussSeidel();
      ~GaussSeidel();

      void set_data_output(double *displ_out,double *angleact_out,double *FtotZMPOut,double *Fc_mc3_out,double *PSurf3,double *ind_cont_out,double *ind_slip_out,double *Kcart_out,double *B_omega_out,double *G_out,double *A_out,double *B_out);
      void setStaticDataSim(const MatrixXd& Cs_in,const Vector6d& FtotZMPdes_in,double friction,int m);
      void setStaticDataStep(const Vector3d& displ_in,const Vector3d& angleact_in,const VectorXd& F3_in,const MatrixXd& Pfree_in,const MatrixXd& PabsOld_in,const Vector3d& displini_in);
	  void run(int contAngle);
	  void gradFtotZ(int contAngle);
      void StiffCart_Omega();
      void NewtonSliding(int index);

    private:
    	//Variables
		int m;
		int mslip;
		int mstick;
		int mc;
		int contAngle;
        double psiini;
		const static double epsiGauss; // Desired precision
        const static double epsiZMP;
		const static double epsiSig; // Signorini's law tolerance
		const static double epsiCou; // Coulomb contact precision
		const static int D; //Dimension
		double fric;
		double theta;
		double phi;
		double psi;
        double critFZMP;
		MatrixXd J;
		MatrixXd J2;
		Matrix3d RPfreet_hat_mc;
		Vector6d FtotZMP;
        Vector6d FtotZMPdes;
        Vector6d FtotZMPerr;
		Vector3d displ;
        Vector3d displ_ini;
        Vector3d angleact;
		MatrixXd displ_m;
		MatrixXd displ_ini_m;
        MatrixXd Q;
		Vector6d displangle;
        Vector6d d_displangle;
		MatrixXd Cs;
		VectorXd F3;
		VectorXd Fold3;
        Matrix3d R;
        Matrix3d Rini;
		Matrix3d Rtheta;
		Matrix3d Rphi;
		Matrix3d Rpsi;
		Matrix3d Rpsiini;
		MatrixXd PabsOld;
		MatrixXd RpabsOld;
		MatrixXd Pfree;
		MatrixXd RPfree;
		MatrixXd RiniPfree;
        VectorXd RPfree3;
        VectorXd delta_sliding3;
        VectorXd P3;
		MatrixArray ci;
		MatrixArray di;
		MatrixArray W;
		MatrixArray inv_W;
		MatrixXd RtF;
		VectorXd RtF3;
		MatrixXd deltaIni;
		VectorXd qt;
		VectorXd qn;
		Vector3d deltaTest;
		Matrix3d A;
		Matrix2d A22;
		double ad;
		double ad12;
		double lambda_min;
		int i3;
		double norm_f_diff;
		double norm_f_new;
		double lambda_max;
		Matrix3d B;
		std::vector<int> ind_slip;
		std::vector<int> ind_c_slip;
		std::vector<int> ind_stick;
		std::vector<int> ind_c_stick;
		std::vector<int> ind_cont;
		std::vector<int> ind_c_cont;
		Vector3d Fcfric;
		Vector3d deltaSliding;
		Matrix2d Bigpi;
		Matrix2d alpha;
		Matrix2d dsMat;
		Vector3d phi1;
		Vector3d d2;
		Vector3d d1;
		Vector3d phi2;
		Vector2d ds;
		Matrix3d G21;
		Matrix3d G22;
		double diffFric;
		int kn;
		double fricTest;
		Vector3d Ftot;
		Vector3d tmpRCRF;
		Vector3d Mo;
		Vector2d Z;
		Vector3d ZMPdes;
		Vector3d Mzmp;

		MatrixXd Pfree_mc;
		MatrixXd W_mc;
		VectorXd Pfree_mc3;
		MatrixXd RPfree_mc;
		VectorXd RPfree_mc3;
		VectorXd Fc_mc3;
		VectorXd RtFc_mc3;
		VectorXd P_mc3;
        VectorXd delta_mc3;

		MatrixXd A_grad;
        MatrixXd A_grad1;        
		MatrixXd B_grad;
        MatrixXd B_grad1;
        MatrixXd B_omega;
        MatrixXd B_omega1;
		MatrixXd D_grad;
        MatrixXd D_grad1;
		MatrixXd Fc_mc_hat;
		MatrixXd tmpRCRF_mc_hat;
		VectorXd tmpRCRF_mc;
        MatrixXd A_1B;
        MatrixXd A_1B1;
        MatrixXd Xi;
        MatrixXd G;
        MatrixXd G1;
        MatrixXd Kcart;
		Vector3d btheta;
		Vector3d bphi;
		Vector3d bpsi;
		Vector3d PfreeFt1;
		Vector3d PfreeFt2;
		Vector3d PfreeFn;
		Matrix3d dRdtheta;
		Matrix3d dRdphi;
		Matrix3d dRdpsi;
		Matrix3d dRinidtheta;
		Matrix3d dRinidphi;
};
