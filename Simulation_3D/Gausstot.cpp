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

const double GaussSeidel::epsiGauss = 1e-12;
const double GaussSeidel::epsiZMP = 1e-12;
const double GaussSeidel::epsiSig = 1e-12;
const double GaussSeidel::epsiCou = 1e-12;
const int GaussSeidel::D = 3;
using namespace std;
using namespace Eigen;

extern void _main();

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

/****************************/

GaussSeidel::GaussSeidel(){}

GaussSeidel::~GaussSeidel(){}

void GaussSeidel::setStaticDataSim(const MatrixXd& Cs_in,const Vector6d& FtotZMPdes_in,double friction,int m1){

    Cs = Cs_in;
    FtotZMPdes = FtotZMPdes_in;
    fric = friction;
    m = m1; // number of contact surface nodes
	FtotZMP.Zero();
	//init
    J.resize(2*D,2*D);
    J2.resize(2*D,2*D);
	displ_m.resize(D,m);
	displ_ini_m.resize(D,m);
    F3.resize(D*m);
    Fold3.resize(D*m);
	RPfree.resize(D,m);
	RiniPfree.resize(D,m);
	RpabsOld.resize(D,m);
	W.resize(m);
	inv_W.resize(m);
    RtF.resize(D,m);
    RtF3.resize(D*m);
	P3.resize(D*m);
	deltaIni.resize(D,m);
	qt.resize(m);
	qn.resize(m);
	ci.resize(m);
	di.resize(m);
	ind_slip.reserve(m);
	ind_c_slip.reserve(m);
	ind_stick.reserve(m);
	ind_c_stick.reserve(m);
	ind_cont.reserve(m);
	ind_c_cont.reserve(m);
	Xi.resize(2*D,2*D);
    Kcart.resize(2*D,2*D);
}

void GaussSeidel::setStaticDataStep(const Vector3d& displ_in,const Vector3d& angleact_in,const VectorXd& F3_in,const MatrixXd& Pfree_in,const MatrixXd& PabsOld_in,const Vector3d& displini_in){
	displ = displ_in;
    angleact = angleact_in;
    F3 = F3_in;
    PabsOld = PabsOld_in;
    Pfree = Pfree_in;
    displ_first = displini_in;
}

void GaussSeidel::run(int contAngle){
    displ_ini = displ;
    angleact_ini = angleact;
	FtotZMP = VectorXd::Zero(6);
	theta = angleact(0);
	phi = angleact(1);
	psi = angleact(2);
	displangle << displ(0),displ(1),displ(2),theta,phi,psi;
    ZMPdes << FtotZMPdes(3),FtotZMPdes(4),0.;
	//mexPrintf("FtotZMPdes = %g %g %g %g %g %g\n", FtotZMPdes(0), FtotZMPdes(1), FtotZMPdes(2), FtotZMPdes(3), FtotZMPdes(4), FtotZMPdes(5));
	int ite = 0;
    int ite_tot = 0;
    int cont_just = 0;
    int flag_Kcart = 0;
	double sum_gauss = 1;
	double norm_delta_displ;
	double norm_delta_angle;
	double max_norm_ddispl = 0.001;/// mm
	double max_norm_dangle = 1.*3.14159/180.;/// Â°
	//// Compute Rotation
	Rtheta << 1, 0, 0,
		    0, cos(theta), -sin(theta),
			0, sin(theta), cos(theta);
	Rphi << cos(phi), 0, sin(phi),
		    0, 1, 0,
			-sin(phi), 0, cos(phi);
	Rpsi << cos(psi), -sin(psi), 0,
		    sin(psi), cos(psi), 0,
			0, 0, 1;
	R = Rtheta * Rphi * Rpsi;
	RPfree = R * Pfree;
	//RpabsOld = R * PabsOld;
	for (int i = 0; i < m; ++i)
	{
		i3 = 3 * i;
		A = R * Cs.block<3, 3>(i3, i3) * R.transpose();
		A22 = A.block<2,2>(0,0);
		W[i] = A;
		inv_W[i] = A.inverse();
		//Find the eigen values
		ad = A22(0,0)/2 + A22(1,1)/2;
		ad12 = sqrt(A22(0,0)*A22(0,0) - 2*A22(0,0)*A22(1,1) + A22(1,1)*A22(1,1) + 4*A22(0,1)*A22(1,0))/2;
		lambda_min = ad - ad12;
		lambda_max = ad + ad12;
		qt(i) = lambda_min/(lambda_max*lambda_max);
		qn(i) = 1/A(2,2);
	}
	if (contAngle==1){
        ///Later define displ_first and psiini as inputs, when optimizing foot position and orientation
		double psiini = 0.;
        Rpsiini << cos(psiini), -sin(psiini), 0,
                sin(psiini), cos(psiini), 0,
                0, 0, 1;
        Rini = Rtheta * Rphi * Rpsiini;
        RiniPfree = Rini * Pfree;
		displ_ini_m = displ_first.replicate(1,m);
		displ_m = displ.replicate(1,m);
		deltaIni.row(0) = displ_m.row(0) + RPfree.row(0) - displ_ini_m.row(0) - RiniPfree.row(0);
		deltaIni.row(1) = displ_m.row(1) + RPfree.row(1) - displ_ini_m.row(1) - RiniPfree.row(1);
		deltaIni.row(2) = displ_m.row(2) + RPfree.row(2);
	}
	else{
		displ_m = displ.replicate(1,m);
		deltaIni.row(0) = RPfree.row(0) - PabsOld.row(0) + displ_m.row(0);
		deltaIni.row(1) = RPfree.row(1) - PabsOld.row(1) + displ_m.row(1);
		deltaIni.row(2) = RPfree.row(2) + displ_m.row(2);
	}
	critFZMP = (FtotZMPdes-FtotZMP).norm();
    while((sum_gauss > epsiGauss) || (critFZMP > epsiZMP))
	{
		int cont_stick = 0;
		int cont_slip = 0;
		int cont_c = 0;
		Fold3 = F3;
		Map<MatrixXd> F = Map<MatrixXd>(F3.data(),D,m);

		RtF = R.transpose() * F;
		RtF3 = Map<VectorXd>(RtF.data(),D*m);
		ind_slip.clear();
		ind_c_slip.clear();
		ind_stick.clear();
		ind_c_stick.clear();
		ind_cont.clear();
		ind_c_cont.clear();
		ite = ite + 1;
		for (int i = 0; i < m; ++i)
		{

			deltaTest = deltaIni.col(i);
			int i3 = 3 * i;

			deltaTest = deltaTest  + R*(Cs.middleCols<3>(i3).transpose()*RtF3 - Cs.block<3, 3>(i3, i3)*RtF3.segment<3>(i3));
			if (deltaTest(2)<epsiSig){
				F3.segment<3>(i3) = -inv_W[i] * deltaTest;
				if (F3.segment<2>(i3).norm() > (fric * F3(i3+2))){
					Fcfric = Vector3d::Zero();
					if (fric>0.2){
						Fcfric = F3.segment<3>(i3);
						deltaNew = Vector3d::Zero();
					}
					else{
						Fcfric(2) = -qn(i)*deltaTest(2); //deltaTest/W(3,3)
						deltaNew = W[i] * Fcfric + deltaTest;
					}

					//Newton Contact for sliding contacts
					diffFric = 1;
					kn = 0;
					G21 = Matrix3d::Zero();
					G22 = Matrix3d::Zero();
					fricTest;
					while(diffFric > epsiCou){
						phi1 = deltaNew-deltaTest-W[i]*Fcfric;
						ds = (Fcfric.head<2>()-qt(i)*deltaNew.head<2>());
						phi2.head<2>() = Fcfric.head<2>()-(fric*Fcfric(2)*ds/ds.norm());
						phi2(2) = qn(i)*deltaNew(2);
						dsMat << pow(ds(1),2), -ds(0)*ds(1),
								-ds(0)*ds(1), pow(ds(0),2);
						Bigpi = pow(1/ds.norm(),3)*dsMat;
						alpha = fric*Fcfric(2)*Bigpi;
						G21.block<2, 2>(0, 0) = qt(i)*alpha;
						G21(2,2) = qn(i);
						G22.block<2, 2>(0, 0) = MatrixXd::Identity(2,2)-alpha;
						G22.block<2, 1>(0, 2) = -fric*ds/ds.norm();
						d2 = (G21*W[i]+G22).inverse() * (phi2-G21*phi1);
						d1 = phi1 + W[i]*d2;
						deltaNew = deltaNew - d1;
						Fcfric = Fcfric - d2;
						fricTest = Fcfric.head<2>().norm()/Fcfric(2);
						diffFric = fabs(fricTest-fric);
						kn = kn+1;
                        if (kn>500){diffFric = 1e-12;}
					}
                    //Vector2d acc;
                    //acc = deltaNew.head<2>().norm()*Fcfric.head<2>() + fric * Fcfric(2) * deltaNew.head<2>();
                    //mexPrintf("FtotZMP = %g\n", acc(1));
					// Find the gradient
					ds = (Fcfric.head<2>()-qt(i)*deltaNew.head<2>());
					Bigpi = pow(1/ds.norm(),3)*dsMat;
					alpha = fric*Fcfric(2)*Bigpi;
					Matrix3d Ci0 = Matrix3d::Zero();
					Ci0.block<2, 2>(0, 0) = MatrixXd::Identity(2,2)-alpha;
					Ci0.block<2, 1>(0, 2) = -fric*ds/ds.norm();
					ci[cont_slip] = Ci0;
					Matrix3d Di0 = Matrix3d::Zero();
					Di0.block<2, 2>(0, 0) = qt(i)*alpha;
					Di0(2, 2) = qn(i);
					di[cont_slip] = Di0;

					F3.segment<3>(i3) = Fcfric;
					ind_slip.push_back(i);
					//ind_slip.push_back(cont_slip);
					ind_c_slip.push_back(cont_c);
					cont_slip = cont_slip + 1;
				}
				else{
					ind_stick.push_back(i);
					ind_c_stick.push_back(cont_c);
					cont_stick = cont_stick + 1;
				}
				ind_cont.push_back(i);
				ind_c_cont.push_back(cont_c);
				cont_c = cont_c + 1;
			}
			else{
				F3.segment<3>(i3) = Vector3d::Zero();
			}

			RtF3.segment<3>(i3) = R.transpose() * F3.segment<3>(i3);
		}

		sum_gauss = 0;
		for (int i = 0; i < m; ++i){
			norm_f_diff = (F3.segment<3>(D*i)-Fold3.segment<3>(D*i)).norm();
			norm_f_new = F3.segment<3>(D*i).norm();
			if (norm_f_new!=0){
				sum_gauss = sum_gauss + norm_f_diff/norm_f_new;
			}
		}

		mslip = (int)ind_slip.size();
		mstick = (int)ind_stick.size();
		mc = (int)ind_cont.size();
		Pfree_mc.resize(D,mc);
		Pfree_mc3.resize(D*mc);
		Fc_mc3.resize(D*mc);
		P_mc3.resize(D*mc);
		RtFc_mc3.resize(D*mc);
		Ftot = Vector3d::Zero();
//		W_mc = MatrixXd::Zero(D*mc,D*mc);
		for (int i = 0; i < mc; ++i){
			Pfree_mc3.segment<3>(D*ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Pfree_mc.col(ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Fc_mc3.segment<3>(D*ind_c_cont[i]) = F3.segment<3>(D*ind_cont[i]);
			RtFc_mc3.segment<3>(D*ind_c_cont[i]) = RtF3.segment<3>(D*ind_cont[i]);
			Ftot = Ftot + Fc_mc3.segment<3>(D*ind_c_cont[i]);
//			for (int j = 0; j < mc; ++j){
//                /// ************************* Not optimized, W_mc only necessary when computing J, that is when changing pos and rot of sole ***************************
//                /// ************************* Could be replaced by first computation of Fc in sole frame, then deformation in sole frame, then changed in absolute frame ***************************
//                /// ************************* and Could be recomputed only when R is changed ***************************
//				W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
//			}
		}
		RPfree_mc = R * Pfree_mc;
		RPfree_mc3 = Map<VectorXd>(RPfree_mc.data(),D*mc);
		for (int i = 0; i < mc; ++i){
			deltaInt = Vector3d::Zero();
			for (int j = 0; j < mc; ++j){
                /// On recalcule tous les dÃ©placements avec toutes les nouvelles forces. Ce calcul est partiellement fait au dessus avec une partie des anciennes forces
                /// Est-il vraiment nÃ©cessaire de le refaire ?
                deltaInt += Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * RtFc_mc3.segment<3>(D*ind_c_cont[j]);
			}
			P_mc3.segment<3>(D*ind_c_cont[i]) = displ + RPfree_mc3.segment<3>(D*ind_c_cont[i]) + R*deltaInt;
		}
		Mo = Vector3d::Zero();
		for (int i = 0; i < mc; ++i){
			Mo = Mo + P_mc3.segment<3>(D*i).cross(Fc_mc3.segment<3>(D*i));
		}
		Z << -Mo(1)/Ftot(2), Mo(0)/Ftot(2);
        Mzmp = Mo - ZMPdes.cross(Ftot);
        FtotZMP << Ftot(0), Ftot(1), Ftot(2), Z(0), Z(1), Mzmp(2);
        FtotZMPerr = FtotZMPdes-FtotZMP;
        critFZMP = FtotZMPerr.norm();
        if ((sum_gauss > epsiGauss) || (critFZMP > epsiZMP)){
            if (mstick>0 & mc>0){
                if ((sum_gauss<((0.0001*FtotZMPerr).norm())) || (sum_gauss<epsiGauss)){
                    ite_tot = ite_tot + 1;
                    gradFtotZ(contAngle);
                    delta_displangle = J.householderQr().solve(FtotZMPerr);
                    //delta_displangle = J.colPivHouseholderQr().solve(FtotZMPerr);
                    //displangle = displangle + J.inverse() *(FtotZMPerr);
                    norm_delta_displ = delta_displangle.segment<3>(0).norm();
                    norm_delta_angle = delta_displangle.segment<3>(3).norm();
                    if(norm_delta_displ != norm_delta_displ){ //norm_d_displ != norm_d_displ if nan value
                        mexPrintf("error \n");
                        critFZMP = 1e-25;
                        sum_gauss = 1e-25;
                    }
                    else{
                        if (norm_delta_displ>max_norm_ddispl) {
                            if (norm_delta_displ*max_norm_dangle>max_norm_ddispl*norm_delta_angle) {
                                delta_displangle *= max_norm_ddispl/norm_delta_displ;
                                //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_displ/max_norm_ddispl);
                            } else {
                                delta_displangle *= max_norm_dangle/norm_delta_angle;
                                //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_angle/max_norm_dangle);
                            }
                        } else if (norm_delta_angle>max_norm_dangle) {
                                delta_displangle *= max_norm_dangle/norm_delta_angle;
                                //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_angle/max_norm_dangle);
                        }
                        displangle = displangle + delta_displangle;
                        displ = displangle.segment<3>(0);
                        angleact = displangle.segment<3>(3);
                        theta = angleact(0);
                        phi = angleact(1);
                        psi = angleact(2);
                        //// Compute Rotation
                        Rtheta << 1, 0, 0,
                                0, cos(theta), -sin(theta),
                                0, sin(theta), cos(theta);
                        Rphi << cos(phi), 0, sin(phi),
                                0, 1, 0,
                                -sin(phi), 0, cos(phi);
                        Rpsi << cos(psi), -sin(psi), 0,
                                sin(psi), cos(psi), 0,
                                0, 0, 1;
                        R = Rtheta * Rphi * Rpsi;
                        RPfree = R * Pfree;
                        //RpabsOld = R * PabsOld;
                        for (int i = 0; i < m; ++i)
                        {
                            i3 = 3 * i;
                            A = R * Cs.block<3, 3>(i3, i3) * R.transpose();
                            A22 = A.block<2,2>(0,0);
                            W[i] = A;
                            inv_W[i] = A.inverse();
                            //Find the eigen values
                            ad = A22(0,0)/2 + A22(1,1)/2;
                            ad12 = sqrt(A22(0,0)*A22(0,0) - 2*A22(0,0)*A22(1,1) + A22(1,1)*A22(1,1) + 4*A22(0,1)*A22(1,0))/2;
                            lambda_min = ad - ad12;
                            lambda_max = ad + ad12;
                            qt(i) = lambda_min/(lambda_max*lambda_max);
                            qn(i) = 1/A(2,2);
                        }
                        if (contAngle==1){
                            ///Later define displ_first and psiini as inputs, when optimizing foot position and orientation
                            double psiini = 0.;
                            Rpsiini << cos(psiini), -sin(psiini), 0,
                                    sin(psiini), cos(psiini), 0,
                                    0, 0, 1;
                            Rini = Rtheta * Rphi * Rpsiini;
                            RiniPfree = Rini * Pfree;
                            displ_ini_m = displ_first.replicate(1,m);
                            displ_m = displ.replicate(1,m);
                            //deltaIni.row(0) = RPfree.row(0) - RpabsOld.row(0) + (displ_ini_m.row(0) + displ_m.row(0));
                            //deltaIni.row(1) = RPfree.row(1) - RpabsOld.row(1) + (displ_ini_m.row(1) + displ_m.row(1));
                            deltaIni.row(0) = displ_m.row(0) + RPfree.row(0) - displ_ini_m.row(0) - RiniPfree.row(0);
                            deltaIni.row(1) = displ_m.row(1) + RPfree.row(1) - displ_ini_m.row(1) - RiniPfree.row(1);
                            deltaIni.row(2) = displ_m.row(2) + RPfree.row(2);
                        }
                        else{
                            displ_m = displ.replicate(1,m);
                            deltaIni.row(0) = RPfree.row(0) - PabsOld.row(0) + displ_m.row(0);
                            deltaIni.row(1) = RPfree.row(1) - PabsOld.row(1) + displ_m.row(1);
                            deltaIni.row(2) = RPfree.row(2) + displ_m.row(2);
                        }
                    }
                } else {
                    //mexPrintf("FtotZMP    = %g %g %g %g %g %g\n", FtotZMP(0), FtotZMP(1), FtotZMP(2), FtotZMP(3), FtotZMP(4), FtotZMP(5));
                }

                if (ite_tot>10000 || ite>1000000){
                     mexPrintf("critFZMP = %g \n",critFZMP);
                     mexPrintf("sum_gauss = %g \n",sum_gauss);
                     critFZMP = 1e-25;
                     sum_gauss = 1e-25;
                 }
            }else{
                if (cont_just>50){
                    critFZMP = 1e-25;
                    sum_gauss = 1e-25;
                    flag_Kcart = 1;
                }
                else{
                    mexPrintf("Just sliding contact \n");
                    critFZMP = 1e-25;
                    sum_gauss = 1e-25;
                    flag_Kcart = 1;
                    /*cont_just = cont_just + 1;

                    FtotZMP = VectorXd::Zero(6);
                    F3= VectorXd::Zero(D*m);
                    displangle << displ(0),displ(1),displ(2),theta,phi,psi;
                    ZMPdes << FtotZMPdes(3),FtotZMPdes(4),0.;
                    Vector3d displ_ad;
                    displ_ad << 0,0,(cont_just*0.001);
                    displ = displ_ini - displ_ad;
                    Vector3d angleact_ad;
                    angleact_ad << 0,0,0;
                    angleact = angleact_ini - angleact_ad;
                    theta = angleact(0);
                    phi = angleact(1);
                    psi = angleact(2);
                    // Compute Rotation
                    Rtheta << 1, 0, 0,
                              0, cos(theta), -sin(theta),
                              0, sin(theta), cos(theta);
                    Rphi << cos(phi), 0, sin(phi),
                            0, 1, 0,
                            -sin(phi), 0, cos(phi);
                    Rpsi << cos(psi), -sin(psi), 0,
                            sin(psi), cos(psi), 0,
                            0, 0, 1;
                    R = Rtheta * Rphi * Rpsi;
                    RPfree = R * Pfree;
                    mexPrintf("Ciao \n");
                    //RpabsOld = R * PabsOld;
                    for (int i = 0; i < m; ++i)
                    {
                        i3 = 3 * i;
                        A = R * Cs.block<3, 3>(i3, i3) * R.transpose();
                        A22 = A.block<2,2>(0,0);
                        W[i] = A;
                        inv_W[i] = A.inverse();
                        //Find the eigen values
                        ad = A22(0,0)/2 + A22(1,1)/2;
                        ad12 = sqrt(A22(0,0)*A22(0,0) - 2*A22(0,0)*A22(1,1) + A22(1,1)*A22(1,1) + 4*A22(0,1)*A22(1,0))/2;
                        lambda_min = ad - ad12;
                        lambda_max = ad + ad12;
                        qt(i) = lambda_min/(lambda_max*lambda_max);
                        qn(i) = 1/A(2,2);
                    }
                    if (contAngle==1){
                        ///Later define displ_first and psiini as inputs, when optimizing foot position and orientation
                        double psiini = 0.;
                        Rpsiini << cos(psiini), -sin(psiini), 0,
                                   sin(psiini), cos(psiini), 0,
                                   0, 0, 1;
                        Rini = Rtheta * Rphi * Rpsiini;
                        RiniPfree = Rini * Pfree;
                        displ_ini_m = displ_first.replicate(1,m);
                        displ_m = displ.replicate(1,m);
                        //deltaIni.row(0) = RPfree.row(0) - RpabsOld.row(0) + (displ_ini_m.row(0) + displ_m.row(0));
                        //deltaIni.row(1) = RPfree.row(1) - RpabsOld.row(1) + (displ_ini_m.row(1) + displ_m.row(1));
                        deltaIni.row(0) = displ_m.row(0) + RPfree.row(0) - displ_ini_m.row(0) - RiniPfree.row(0);
                        deltaIni.row(1) = displ_m.row(1) + RPfree.row(1) - displ_ini_m.row(1) - RiniPfree.row(1);
                        deltaIni.row(2) = displ_m.row(2) + RPfree.row(2);
                    }
                    else{
                        displ_m = displ.replicate(1,m);
                        deltaIni.row(0) = RPfree.row(0) - PabsOld.row(0) + displ_m.row(0);
                        deltaIni.row(1) = RPfree.row(1) - PabsOld.row(1) + displ_m.row(1);
                        deltaIni.row(2) = RPfree.row(2) + displ_m.row(2);
                    }*/
                }
            }
        } else {
        }

	}
    if(flag_Kcart==0){
        StiffCart_Omega();
    }else{
        Kcart << 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0;
    }

    //Find position of the surface points
    RPfree = R * Pfree;
    RPfree3 = Map<VectorXd>(RPfree.data(),D*m);
    for (int i = 0; i < m; ++i){
        deltaInt = Vector3d::Zero();
        for (int j = 0; j < m; ++j){
            /// On recalcule tous les déplacements avec toutes les nouvelles forces. Ce calcul est partiellement fait au dessus avec une partie des anciennes forces
            /// Est-il vraiment nécessaire de le refaire ?
            deltaInt += Cs.block<3, 3>(D*i,D*j) * RtF3.segment<3>(D*j);
        }
        P3.segment<3>(D*i) = displ + RPfree3.segment<3>(D*i) + R*deltaInt;
    }
    //gradFtotZ(contAngle);
}

void GaussSeidel::gradFtotZ(int contAngle)
{

	A_grad = MatrixXd::Zero(D*mc+D*mslip,D*mc+D*mslip);
	B_grad = MatrixXd::Zero(D*mc+D*mslip,2*D);
	D_grad = MatrixXd::Zero(2*D,D*mc+D*mslip);
	int cont_slip1 = 0;
	/* Compute derivative of rotation */
	dRdtheta << 0, 0, 0,
				 0, -sin(theta), -cos(theta),
				 0, cos(theta), -sin(theta);
    if (contAngle==1){
        dRinidtheta = dRdtheta*Rphi*Rpsiini;
    }
	dRdtheta =	dRdtheta*Rphi*Rpsi;
	dRdphi << -sin(phi), 0, cos(phi),
				0, 0, 0,
				-cos(phi), 0, -sin(phi);
    if (contAngle==1){
        dRinidphi = Rtheta*dRdphi*Rpsiini;
    }
	dRdphi = Rtheta*dRdphi*Rpsi;
	dRdpsi << -sin(psi), -cos(psi), 0,
			cos(psi), -sin(psi), 0,
			0, 0, 0;
	dRdpsi = Rtheta*Rphi*dRdpsi;
    if (contAngle==1){
        PfreeFt1 << 0, 0, 0;
        PfreeFt2 << 0, 0, 0;
        PfreeFn << 0, 0, 0;
    }
	// find A_grad
    W_mc = MatrixXd::Zero(D*mc,D*mc);
    for (int i = 0; i < mc; ++i){
        for (int j = 0; j < mc; ++j){
            W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
        }
    }
	A_grad.block(0,0,D*mc,D*mc) = -W_mc;

	for (int i = 0; i < mc; ++i){
		B_grad.block<3, 3>(D*i,0) = Matrix3d::Identity();
		D_grad.block<3, 3>(0,D*i) = Matrix3d::Identity();
		D_grad.block<2, 1>(D,D*i+2) = (1/Ftot(2)) * (P_mc3.segment<2>(D*i) - Z.segment<2>(0));
        D_grad(D+2,D*i) = -P_mc3(D*i+1) + ZMPdes(1);//Z(1);
        D_grad(D+2,D*i+1) = P_mc3(D*i) - ZMPdes(0);//Z(0);
		if (contAngle==1){
            btheta = dRdtheta*Pfree_mc3.segment<3>(D*i);
            bphi = dRdphi*Pfree_mc3.segment<3>(D*i);
            bpsi = dRdpsi*Pfree_mc3.segment<3>(D*i);
            btheta(0) -= dRinidtheta.row(0)*Pfree_mc3.segment<3>(D*i);
            btheta(1) -= dRinidtheta.row(1)*Pfree_mc3.segment<3>(D*i);
            bphi(0) -= dRinidphi.row(0)*Pfree_mc3.segment<3>(D*i);
            bphi(1) -= dRinidphi.row(1)*Pfree_mc3.segment<3>(D*i);
            PfreeFt1 = PfreeFt1 + Pfree_mc3.segment<3>(D*i)*Fc_mc3(D*i);
            PfreeFt2 = PfreeFt2 + Pfree_mc3.segment<3>(D*i)*Fc_mc3(D*i+1);
            PfreeFn = PfreeFn + Pfree_mc3.segment<3>(D*i)*Fc_mc3(D*i+2);
		} else {
            btheta = dRdtheta*Pfree_mc3.segment<3>(D*i);
            bphi = dRdphi*Pfree_mc3.segment<3>(D*i);
            bpsi = dRdpsi*Pfree_mc3.segment<3>(D*i);
		}
		for (int j = 0; j < mc; ++j){
			btheta = btheta + (dRdtheta * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose() + R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * dRdtheta.transpose()) * Fc_mc3.segment<3>(D*j);
			bphi = bphi + (dRdphi * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose() + R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * dRdphi.transpose()) * Fc_mc3.segment<3>(D*j);
			bpsi = bpsi + (dRdpsi * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose() + R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * dRdpsi.transpose()) * Fc_mc3.segment<3>(D*j);
		}
		B_grad.block<3, 1>(D*i,D) = btheta;
		B_grad.block<3, 1>(D*i,D+1) = bphi;
		B_grad.block<3, 1>(D*i,D+2) = bpsi;
		if (cont_slip1 < ind_c_slip.size()){
			if (i == ind_c_slip[cont_slip1]){
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*i) = ci[cont_slip1];
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*mc+D*cont_slip1) = di[cont_slip1];
				A_grad.block<3, 3>(D*i,D*mc+D*cont_slip1) = MatrixXd::Identity(D,D);
				D_grad(D,D*mc+D*cont_slip1) = (1/Ftot(2)) * Fc_mc3(D*i+2);
				D_grad(D+1,D*mc+D*cont_slip1+1) = (1/Ftot(2)) * Fc_mc3(D*i+2);
				D_grad(D+2,D*mc+D*cont_slip1) = Fc_mc3(D*i+1);
				D_grad(D+2,D*mc+D*cont_slip1+1) = -Fc_mc3(D*i);
				cont_slip1 = cont_slip1 + 1;
			}
		}
	}
//    A_1B = A_grad.inverse() * B_grad;
    A_1B = A_grad.householderQr().solve(B_grad);
	J = D_grad * A_1B;
    if (contAngle==1){
        J2 = MatrixXd::Zero(2*D,2*D);
        J2(D,D) = dRinidtheta.row(0).dot(PfreeFn)/Ftot(2);
        J2(D,D+1) = dRinidphi.row(0).dot(PfreeFn)/Ftot(2);
        J2(D+1,D) = dRinidtheta.row(1).dot(PfreeFn)/Ftot(2);
        J2(D+1,D+1) = dRinidphi.row(1).dot(PfreeFn)/Ftot(2);
        //J2(D+2,D) = dRinidtheta.row(0).dot(PfreeFt2) - dRdtheta.row(1).dot(PfreeFt1);
        //J2(D+2,D+1) = dRinidphi.row(0).dot(PfreeFt2) - dRdphi.row(1).dot(PfreeFt1);
        J2(D+2,D) = dRinidtheta.row(0).dot(PfreeFt2) - dRinidtheta.row(1).dot(PfreeFt1);
        J2(D+2,D+1) = dRinidphi.row(0).dot(PfreeFt2) - dRinidphi.row(1).dot(PfreeFt1);        
        J = J + J2;
    }
}

/*void GaussSeidel::StiffCart(int contAngle)
{
	int cont_slip1 = 0;
	EZ = MatrixXd::Zero(D+3,D*mc+D*mslip);
	for (int i = 0; i < mc; ++i){
		EZ.block<3, 3>(0,D*i) = Matrix3d::Identity();
		EZ(D,D*i+2) = (P_mc3(D*i+1)-Z(1));
		EZ(D+1,D*i+2) = -(P_mc3(D*i)-Z(0));
		EZ(D+2,D*i) = (-P_mc3(D*i+1)+Z(1));
		EZ(D+2,D*i+1) = (P_mc3(D*i)-Z(0));
		if (cont_slip1 < ind_c_slip.size()){
			if (i == ind_c_slip[cont_slip1]){
				EZ(D,D*mc+D*cont_slip1+1) = Fc_mc3(D*i+2);
				EZ(D+1,D*mc+D*cont_slip1) = Fc_mc3(D*i+2);
				EZ(D+2,D*mc+D*cont_slip1) = Fc_mc3(D*i+1);
				EZ(D+2,D*mc+D*cont_slip1+1) = -Fc_mc3(D*i);
				cont_slip1 = cont_slip1 + 1;
			}
		}

	}
	//MatrixXd Xi;
	Xi << 1, 0, 0,      0,      0,  -Z(1),
		0, 1, 0,        0,      0,  Z(0),
		0, 0, 1,        Z(1), -Z(0), 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 1;
	Kcart = EZ * A_1B;// * Xi;
    if (contAngle==1){
        J2 = MatrixXd::Zero(2*D,2*D);
        J2(D,D) = dRdtheta.row(1).dot(PfreeFn);
        J2(D,D+1) = dRdphi.row(1).dot(PfreeFn);
        J2(D,D+2) = dRdpsi.row(1).dot(PfreeFn);
        J2(D+1,D) = -dRdtheta.row(0).dot(PfreeFn);
        J2(D+1,D+1) = -dRdphi.row(0).dot(PfreeFn);
        J2(D+1,D+2) = -dRdpsi.row(0).dot(PfreeFn);
        J2(D+2,D) = dRdtheta.row(0).dot(PfreeFt2) - dRdtheta.row(1).dot(PfreeFt1);
        J2(D+2,D+1) = dRdphi.row(0).dot(PfreeFt2) - dRdphi.row(1).dot(PfreeFt1);
        J2(D+2,D+2) = dRdpsi.row(0).dot(PfreeFt2) - dRdpsi.row(1).dot(PfreeFt1);
        Kcart += J2;
    }
}*/

void GaussSeidel::StiffCart_Omega()
{
	G = MatrixXd::Zero(2*D,D*mc+D*mslip);
	A_grad = MatrixXd::Zero(D*mc+D*mslip,D*mc+D*mslip);
	B_omega = MatrixXd::Zero(D*mc+D*mslip,2*D);
	Fc_mc_hat.resize(D*mc,D);
	delta_mc_hat.resize(D*mc,D);
	delta_mc.resize(D*mc);
	int cont_slip1 = 0;
	// find A_grad
    W_mc = MatrixXd::Zero(D*mc,D*mc);
    for (int i = 0; i < mc; ++i){
        for (int j = 0; j < mc; ++j){
            W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
        }
    }
	A_grad.block(0,0,D*mc,D*mc) = -W_mc;

    delta_mc = W_mc * Fc_mc3;
    for (int i = 0; i < mc; ++i){
//        cross(Fc_mc3.segment<3>(D*i), Fc_mc_hat.block<3,3>(D*i,0));
//        cross(delta_mc.segment<3>(D*i), delta_mc_hat.block<3,3>(D*i,0));
        Fc_mc_hat.block<3,3>(D*i,0) << 0, -Fc_mc3(D*i+2), Fc_mc3(D*i+1),
                    Fc_mc3(D*i+2), 0, -Fc_mc3(D*i),
                    -Fc_mc3(D*i+1), Fc_mc3(D*i), 0;
        delta_mc_hat.block<3,3>(D*i,0) << 0, -delta_mc(D*i+2), delta_mc(D*i+1),
                    delta_mc(D*i+2), 0, -delta_mc(D*i),
                    -delta_mc(D*i+1), delta_mc(D*i), 0;
    }

    B_omega.block(0,D,D*mc,D) = W_mc * Fc_mc_hat;
    B_omega.block(0,D,D*mc,D) -= delta_mc_hat;

	for (int i = 0; i < mc; ++i){
		B_omega.block<3, 3>(D*i,0) = Matrix3d::Identity();
		G.block<3, 3>(0,D*i) = Matrix3d::Identity();
		G(D,D*i+2) = (P_mc3(D*i+1)-Z(1));
		G(D+1,D*i+2) = -(P_mc3(D*i)-Z(0));
		G(D+2,D*i) = (-P_mc3(D*i+1)+Z(1));
		G(D+2,D*i+1) = (P_mc3(D*i)-Z(0));
        /// Finally do not consider special problem of step 1 for cartesian stiffness
        /// Once contact state obtained, compute considering nodes stuck as other steps
        B_omega(D*i,D+1) += RPfree_mc3(D*i+2);
        B_omega(D*i,D+2) -= RPfree_mc3(D*i+1);
        B_omega(D*i+1,D) -= RPfree_mc3(D*i+2);
        B_omega(D*i+1,D+2) += RPfree_mc3(D*i);
        B_omega(D*i+2,D) += RPfree_mc3(D*i+1);
        B_omega(D*i+2,D+1) -= RPfree_mc3(D*i);
		if (cont_slip1 < ind_c_slip.size()){
			if (i == ind_c_slip[cont_slip1]){
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*i) = ci[cont_slip1];
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*mc+D*cont_slip1) = di[cont_slip1];
				A_grad.block<3, 3>(D*i,D*mc+D*cont_slip1) = MatrixXd::Identity(D,D);
				G(D,D*mc+D*cont_slip1+1) = Fc_mc3(D*i+2);
				G(D+1,D*mc+D*cont_slip1) = -Fc_mc3(D*i+2);
				G(D+2,D*mc+D*cont_slip1) = Fc_mc3(D*i+1);
				G(D+2,D*mc+D*cont_slip1+1) = -Fc_mc3(D*i);
				cont_slip1 = cont_slip1 + 1;
			}
		}

	}
	//MatrixXd Xi;
	//mexPrintf("Compute Xi \n");
	Xi << 1, 0, 0,      0,      displ(2),  -displ(1)+ZMPdes(1),
		0, 1, 0,        -displ(2),      0,   displ(0)-ZMPdes(0),
		0, 0, 1,        displ(1)-ZMPdes(1), -displ(0)+ZMPdes(0), 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 1;
//    A_1B = A_grad.inverse() * B_omega;
	//mexPrintf("Compute A^-1*B \n");
    A_1B = A_grad.householderQr().solve(B_omega);
	//mexPrintf("Compute Kcart \n");
    /// Finally do not consider special problem of step 1 for cartesian stiffness
    /// Once contact state obtained, compute considering nodes stuck as other steps
	Kcart = (G * A_1B) * Xi;
}


void GaussSeidel::set_data_output(double *displ_out,double *angleact_out,double *FtotZMP_out,double *Fc_mc3_out,double *PSurf3,double *ind_cont_out,double *ind_slip_out,double *Kcart_out)
{
    Map<Vector6d>(FtotZMP_out,FtotZMP.rows()) = FtotZMP;
    Map<Vector3d>(displ_out,displ.rows()) = displ;
    Map<Vector3d>(angleact_out,angleact.rows()) = angleact;
    Map<VectorXd>(Fc_mc3_out,Fc_mc3.rows()) = Fc_mc3;
    Map<VectorXd>(PSurf3,P3.rows()) = P3;
    for (mwSize i=0; i<mc; ++i){
        ind_cont_out[i] = ind_cont[i];
    }
    for (mwSize i=0; i<mslip; ++i){
        ind_slip_out[i] = ind_slip[i];
    }    
    Map<MatrixXd>(Kcart_out,Kcart.rows(),Kcart.cols()) = Kcart;
//     for (mwSize i=0; i<B_omega.rows(); ++i){
// 		for (mwSize j=0; j<B_omega.cols(); ++j){
// 			B_omega_out[(j*(B_omega.rows()))+i] = B_omega(i,j);
// 		}
//     }
//     for (mwSize i=0; i<G.rows(); ++i){
// 		for (mwSize j=0; j<G.cols(); ++j){
// 			G_out[(j*(G.rows()))+i] = G(i,j);
// 		}
//     }
//     
//     for (mwSize i=0; i<A_grad.rows(); ++i){
// 		for (mwSize j=0; j<A_grad.cols(); ++j){
// 			A_out[(j*(A_grad.rows()))+i] = A_grad(i,j);
// 		}
//     }
//     for (mwSize i=0; i<B_grad.rows(); ++i){
// 		for (mwSize j=0; j<B_grad.cols(); ++j){
// 			B_out[(j*(B_grad.rows()))+i] = B_grad(i,j);
// 		}
//     }        
    //Map<MatrixXd>(Gtmp,B_omega.rows(),B_omega.cols()) = B_omega;
    //B_omega_out.resize(D*mc+D*mslip,2*D);
    //Map<MatrixXd>(B_omega_out,B_omega.rows(),B_omega.cols()) = B_omega;
    /*for (int i = 0; i < (D*mc+D*mslip); ++i){
        for (int j = 0; j < 6; ++j){
            B_omega_out[i,j] = Gtmp[i][j] ;
        }
    }*/
    //Map<MatrixXd>(B_omega_out,B_omega.rows(),B_omega.cols()) = B_omega;
    //Map<MatrixXd>(B_omega_out,B_omega.rows(),B_omega.cols()) = B_omega;
}

/*********************/
