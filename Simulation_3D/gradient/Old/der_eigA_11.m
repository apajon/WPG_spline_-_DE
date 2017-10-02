function der_eigA_11 = der_eigA_11(A)
A1_1=A(1,1);
A1_2=A(1,2);
A1_3=A(1,3);
A2_1=A(2,1);
A2_2=A(2,2);
A2_3=A(2,3);
A3_1=A(3,1);
A3_2=A(3,2);
A3_3=A(3,3);

%DER_EIGA_11
%    DER_EIGA_11 = DER_EIGA_11(A1_1,A1_2,A1_3,A2_1,A2_2,A2_3,A3_1,A3_2,A3_3)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    15-Jul-2015 23:01:40

t2 = A1_1+A2_2+A3_3;
t3 = t2.^2;
t4 = A1_1.*A2_2;
t5 = A1_1.*A3_3;
t6 = A2_2.*A3_3;
t19 = A1_2.*A2_1;
t20 = A1_3.*A3_1;
t21 = A2_3.*A3_2;
t7 = t4+t5+t6-t19-t20-t21;
t8 = t2.*t3.*(1.0./2.7e1);
t12 = A1_1.*A2_2.*A3_3.*(1.0./2.0);
t13 = A1_1.*A2_3.*A3_2.*(1.0./2.0);
t14 = A1_2.*A2_1.*A3_3.*(1.0./2.0);
t15 = A1_2.*A2_3.*A3_1.*(1.0./2.0);
t16 = A1_3.*A2_1.*A3_2.*(1.0./2.0);
t17 = A1_3.*A2_2.*A3_1.*(1.0./2.0);
t18 = t2.*t7.*(1.0./6.0);
t9 = t8+t12-t13-t14+t15+t16-t17-t18;
t23 = t3.*(1.0./9.0);
t24 = A1_1.*A2_2.*(1.0./3.0);
t25 = A1_2.*A2_1.*(1.0./3.0);
t26 = A1_1.*A3_3.*(1.0./3.0);
t27 = A1_3.*A3_1.*(1.0./3.0);
t28 = A2_2.*A3_3.*(1.0./3.0);
t29 = A2_3.*A3_2.*(1.0./3.0);
t10 = t23-t24+t25-t26+t27-t28+t29;
t11 = t10.^2;
t22 = t8+t12-t13-t14+t15+t16-t17-t18;
t30 = A2_2+A3_3;
t31 = A2_2.*(1.0./9.0);
t32 = A3_3.*(1.0./9.0);
t39 = A1_1.*(2.0./9.0);
t33 = t31+t32-t39;
t34 = A1_2.*A2_1.*(1.0./6.0);
t35 = A1_3.*A3_1.*(1.0./6.0);
t36 = t8+t12-t13-t14+t15+t16-t17-t18;
t37 = t8+t12-t13-t14+t15+t16-t17-t18;
t42 = A1_1.*A2_2.*(1.0./6.0);
t43 = A1_1.*A3_3.*(1.0./6.0);
t45 = t2.*t30.*(1.0./6.0);
t38 = (t8+t12-t13-t14+t15+t16-t17-t18).*(t23+t28-t29+t34+t35-t42-t43-t45).*2.0;
t40 = t11.*t33.*3.0;
t41 = t38+t40;
t44 = t8+t12-t13-t14+t15+t16-t17-t18;
t46 = t8+t12-t13-t14+t15+t16-t17-t18;
t47 = t8+t12-t13-t14+t15+t16-t17-t18;
t48 = t8+t12-t13-t14+t15+t16-t17-t18;
t49 = t8+t12-t13-t14+t15+t16-t17-t18;
t50 = sqrt(3.0);
t51 = t8+t12-t13-t14+t15+t16-t17-t18;
t52 = t8+t12-t13-t14+t15+t16-t17-t18;
t53 = t8+t12-t13-t14+t15+t16-t17-t18;
t54 = t8+t12-t13-t14+t15+t16-t17-t18;
t55 = t8+t12-t13-t14+t15+t16-t17-t18;
t56 = t8+t12-t13-t14+t15+t16-t17-t18;
t57 = t8+t12-t13-t14+t15+t16-t17-t18;
t58 = t8+t12-t13-t14+t15+t16-t17-t18;
t59 = t8+t12-t13-t14+t15+t16-t17-t18;
t60 = t8+t12-t13-t14+t15+t16-t17-t18;
t61 = t8+t12-t13-t14+t15+t16-t17-t18;
t62 = t8+t12-t13-t14+t15+t16-t17-t18;
t63 = t8+t12-t13-t14+t15+t16-t17-t18;
t64 = t8+t12-t13-t14+t15+t16-t17-t18;
t65 = t8+t12-t13-t14+t15+t16-t17-t18;
t66 = t8+t12-t13-t14+t15+t16-t17-t18;
der_eigA_11 = [-t33.*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t9.^2)).^(1.0./3.0)+1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t36.^2)).^(2.0./3.0).*(t23+t28-t29+t34+t35-t2.*t30.*(1.0./6.0)+t41.*1.0./sqrt(-t10.*t11+t22.^2).*(1.0./2.0)-A1_1.*A2_2.*(1.0./6.0)-A1_1.*A3_3.*(1.0./6.0)).*(1.0./3.0)-t10.*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t37.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t44.^2)).^(4.0./3.0).*(1.0./3.0)+1.0./3.0;(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t47.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t48.^2)).^(2.0./3.0).*(-1.0./6.0)+t33.*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t46.^2)).^(1.0./3.0).*(1.0./2.0)+t10.*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t53.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t54.^2)).^(4.0./3.0).*(1.0./6.0)+t50.*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t51.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t52.^2)).^(2.0./3.0).*1.666666666666667e-1i+t33.*t50.*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t49.^2)).^(1.0./3.0).*5.0e-1i+t10.*t50.*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t55.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t56.^2)).^(4.0./3.0).*1.666666666666667e-1i+1.0./3.0;(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t58.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t59.^2)).^(2.0./3.0).*(-1.0./6.0)+t33.*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t57.^2)).^(1.0./3.0).*(1.0./2.0)+t10.*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t63.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t64.^2)).^(4.0./3.0).*(1.0./6.0)-t50.*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t61.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t62.^2)).^(2.0./3.0).*1.666666666666667e-1i-t50.*(A1_1.*(-1.0./9.0)+A2_2.*(1.0./1.8e1)+A3_3.*(1.0./1.8e1)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t60.^2)).^(1.0./3.0).*1i-t50.*(t3.*(1.0./1.8e1)+t34+t35-t42-t43-A2_2.*A3_3.*(1.0./6.0)+A2_3.*A3_2.*(1.0./6.0)).*(t23+t28-t29+t34+t35-t42-t43-t45+t41.*1.0./sqrt(-t10.*t11+t65.^2).*(1.0./2.0)).*1.0./(t8+t12-t13-t14+t15+t16-t17-t18+sqrt(-t10.*t11+t66.^2)).^(4.0./3.0).*3.333333333333333e-1i+1.0./3.0];