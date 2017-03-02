function CSymbols = CFunc(q2,q3,q4,q5,q1D,q2D,q3D,q4D,q5D)
%CFUNC
%    CSYMBOLS = CFUNC(Q2,Q3,Q4,Q5,Q1D,Q2D,Q3D,Q4D,Q5D)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    02-Mar-2017 10:10:13

t2 = q2D.^2;
t3 = q3D.^2;
t4 = q2+q3+q4;
t5 = sin(t4);
t6 = q4D.^2;
t7 = q3+q4+q5;
t8 = sin(t7);
t9 = q5D.^2;
t10 = q2+q3+q4+q5;
t11 = sin(t10);
t12 = q2+q3;
t13 = sin(t12);
t14 = q3+q4;
t15 = sin(t14);
t16 = q4+q5;
t17 = sin(t16);
t18 = sin(q2);
t19 = sin(q3);
t20 = sin(q4);
t21 = sin(q5);
t22 = q1D.^2;
t23 = t5.*t22.*(3.0./2.5e2);
t24 = t11.*t22.*(1.0./2.5e2);
t25 = t13.*t22.*(1.0./5.0e1);
t26 = t8.*t22.*(1.0./2.5e2);
t27 = t2.*t8.*(1.0./2.5e2);
t28 = t15.*t22.*(3.0./2.5e2);
t29 = t2.*t15.*(3.0./2.5e2);
t30 = q1D.*q2D.*t8.*(1.0./1.25e2);
t31 = q1D.*q2D.*t15.*(3.0./1.25e2);
t32 = t17.*t22.*(1.0./2.5e2);
t33 = t2.*t17.*(1.0./2.5e2);
t34 = t3.*t17.*(1.0./2.5e2);
t35 = q1D.*q2D.*t17.*(1.0./1.25e2);
t36 = q1D.*q3D.*t17.*(1.0./1.25e2);
t37 = q2D.*q3D.*t17.*(1.0./1.25e2);
CSymbols = [t2.*t5.*(-3.0./2.5e2)-t3.*t5.*(3.0./2.5e2)-t3.*t8.*(1.0./2.5e2)-t5.*t6.*(3.0./2.5e2)-t2.*t11.*(1.0./2.5e2)-t3.*t11.*(1.0./2.5e2)-t6.*t8.*(1.0./2.5e2)-t2.*t13.*(1.0./5.0e1)-t3.*t13.*(1.0./5.0e1)-t6.*t11.*(1.0./2.5e2)-t8.*t9.*(1.0./2.5e2)-t3.*t15.*(3.0./2.5e2)-t2.*t18.*(7.0./2.5e2)-t9.*t11.*(1.0./2.5e2)-t6.*t15.*(3.0./2.5e2)-t3.*t19.*(1.0./5.0e1)-t6.*t17.*(1.0./2.5e2)-t6.*t20.*(3.0./2.5e2)-t9.*t17.*(1.0./2.5e2)-t9.*t21.*(1.0./2.5e2)-q1D.*q2D.*t5.*(3.0./1.25e2)-q1D.*q3D.*t5.*(3.0./1.25e2)-q1D.*q4D.*t5.*(3.0./1.25e2)-q2D.*q3D.*t5.*(3.0./1.25e2)-q2D.*q4D.*t5.*(3.0./1.25e2)-q1D.*q3D.*t8.*(1.0./1.25e2)-q3D.*q4D.*t5.*(3.0./1.25e2)-q1D.*q4D.*t8.*(1.0./1.25e2)-q2D.*q3D.*t8.*(1.0./1.25e2)-q1D.*q2D.*t11.*(1.0./1.25e2)-q1D.*q5D.*t8.*(1.0./1.25e2)-q2D.*q4D.*t8.*(1.0./1.25e2)-q1D.*q3D.*t11.*(1.0./1.25e2)-q2D.*q5D.*t8.*(1.0./1.25e2)-q3D.*q4D.*t8.*(1.0./1.25e2)-q1D.*q2D.*t13.*(1.0./2.5e1)-q1D.*q4D.*t11.*(1.0./1.25e2)-q2D.*q3D.*t11.*(1.0./1.25e2)-q3D.*q5D.*t8.*(1.0./1.25e2)-q1D.*q3D.*t13.*(1.0./2.5e1)-q1D.*q5D.*t11.*(1.0./1.25e2)-q2D.*q4D.*t11.*(1.0./1.25e2)-q4D.*q5D.*t8.*(1.0./1.25e2)-q2D.*q3D.*t13.*(1.0./2.5e1)-q2D.*q5D.*t11.*(1.0./1.25e2)-q3D.*q4D.*t11.*(1.0./1.25e2)-q1D.*q3D.*t15.*(3.0./1.25e2)-q3D.*q5D.*t11.*(1.0./1.25e2)-q1D.*q4D.*t15.*(3.0./1.25e2)-q2D.*q3D.*t15.*(3.0./1.25e2)-q4D.*q5D.*t11.*(1.0./1.25e2)-q1D.*q2D.*t18.*(7.0./1.25e2)-q2D.*q4D.*t15.*(3.0./1.25e2)-q1D.*q4D.*t17.*(1.0./1.25e2)-q3D.*q4D.*t15.*(3.0./1.25e2)-q1D.*q3D.*t19.*(1.0./2.5e1)-q1D.*q5D.*t17.*(1.0./1.25e2)-q2D.*q4D.*t17.*(1.0./1.25e2)-q2D.*q3D.*t19.*(1.0./2.5e1)-q2D.*q5D.*t17.*(1.0./1.25e2)-q3D.*q4D.*t17.*(1.0./1.25e2)-q1D.*q4D.*t20.*(3.0./1.25e2)-q3D.*q5D.*t17.*(1.0./1.25e2)-q2D.*q4D.*t20.*(3.0./1.25e2)-q4D.*q5D.*t17.*(1.0./1.25e2)-q1D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q4D.*t20.*(3.0./1.25e2)-q2D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q5D.*t21.*(1.0./1.25e2)-q4D.*q5D.*t21.*(1.0./1.25e2);t23+t24+t25+t26+t27+t28+t29+t30+t31+t2.*t19.*(1.0./5.0e1)-t6.*t17.*(1.0./2.5e2)-t6.*t20.*(3.0./2.5e2)-t9.*t17.*(1.0./2.5e2)-t9.*t21.*(1.0./2.5e2)+t19.*t22.*(1.0./5.0e1)+q1D.*q2D.*t19.*(1.0./2.5e1)-q1D.*q4D.*t17.*(1.0./1.25e2)-q1D.*q5D.*t17.*(1.0./1.25e2)-q2D.*q4D.*t17.*(1.0./1.25e2)-q2D.*q5D.*t17.*(1.0./1.25e2)-q3D.*q4D.*t17.*(1.0./1.25e2)-q1D.*q4D.*t20.*(3.0./1.25e2)-q3D.*q5D.*t17.*(1.0./1.25e2)-q2D.*q4D.*t20.*(3.0./1.25e2)-q4D.*q5D.*t17.*(1.0./1.25e2)-q1D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q4D.*t20.*(3.0./1.25e2)-q2D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q5D.*t21.*(1.0./1.25e2)-q4D.*q5D.*t21.*(1.0./1.25e2);t23+t24+t25-t3.*t8.*(1.0./2.5e2)-t6.*t8.*(1.0./2.5e2)-t8.*t9.*(1.0./2.5e2)-t3.*t15.*(3.0./2.5e2)-t6.*t15.*(3.0./2.5e2)-t3.*t19.*(1.0./5.0e1)-t6.*t17.*(1.0./2.5e2)-t6.*t20.*(3.0./2.5e2)-t9.*t17.*(1.0./2.5e2)-t9.*t21.*(1.0./2.5e2)+t18.*t22.*(7.0./2.5e2)-q1D.*q3D.*t8.*(1.0./1.25e2)-q1D.*q4D.*t8.*(1.0./1.25e2)-q2D.*q3D.*t8.*(1.0./1.25e2)-q1D.*q5D.*t8.*(1.0./1.25e2)-q2D.*q4D.*t8.*(1.0./1.25e2)-q2D.*q5D.*t8.*(1.0./1.25e2)-q3D.*q4D.*t8.*(1.0./1.25e2)-q3D.*q5D.*t8.*(1.0./1.25e2)-q4D.*q5D.*t8.*(1.0./1.25e2)-q1D.*q3D.*t15.*(3.0./1.25e2)-q1D.*q4D.*t15.*(3.0./1.25e2)-q2D.*q3D.*t15.*(3.0./1.25e2)-q2D.*q4D.*t15.*(3.0./1.25e2)-q1D.*q4D.*t17.*(1.0./1.25e2)-q3D.*q4D.*t15.*(3.0./1.25e2)-q1D.*q3D.*t19.*(1.0./2.5e1)-q1D.*q5D.*t17.*(1.0./1.25e2)-q2D.*q4D.*t17.*(1.0./1.25e2)-q2D.*q3D.*t19.*(1.0./2.5e1)-q2D.*q5D.*t17.*(1.0./1.25e2)-q3D.*q4D.*t17.*(1.0./1.25e2)-q1D.*q4D.*t20.*(3.0./1.25e2)-q3D.*q5D.*t17.*(1.0./1.25e2)-q2D.*q4D.*t20.*(3.0./1.25e2)-q4D.*q5D.*t17.*(1.0./1.25e2)-q1D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q4D.*t20.*(3.0./1.25e2)-q2D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q5D.*t21.*(1.0./1.25e2)-q4D.*q5D.*t21.*(1.0./1.25e2);t23+t24+t26+t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t2.*t20.*(3.0./2.5e2)+t3.*t20.*(3.0./2.5e2)-t9.*t21.*(1.0./2.5e2)+t20.*t22.*(3.0./2.5e2)+q1D.*q2D.*t20.*(3.0./1.25e2)+q1D.*q3D.*t20.*(3.0./1.25e2)+q2D.*q3D.*t20.*(3.0./1.25e2)-q1D.*q5D.*t21.*(1.0./1.25e2)-q2D.*q5D.*t21.*(1.0./1.25e2)-q3D.*q5D.*t21.*(1.0./1.25e2)-q4D.*q5D.*t21.*(1.0./1.25e2);t24+t26+t27+t30+t32+t33+t34+t35+t36+t37+t2.*t21.*(1.0./2.5e2)+t3.*t21.*(1.0./2.5e2)+t6.*t21.*(1.0./2.5e2)+t21.*t22.*(1.0./2.5e2)+q1D.*q2D.*t21.*(1.0./1.25e2)+q1D.*q3D.*t21.*(1.0./1.25e2)+q1D.*q4D.*t21.*(1.0./1.25e2)+q2D.*q3D.*t21.*(1.0./1.25e2)+q2D.*q4D.*t21.*(1.0./1.25e2)+q3D.*q4D.*t21.*(1.0./1.25e2)];
