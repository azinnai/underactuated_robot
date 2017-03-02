function BSymbols = BFunc(q2,q3,q4,q5)
%BFUNC
%    BSYMBOLS = BFUNC(Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    02-Mar-2017 10:10:12

t2 = q2+q3+q4;
t3 = cos(t2);
t4 = q3+q4+q5;
t5 = cos(t4);
t6 = t5.*(1.0./1.25e2);
t7 = q2+q3+q4+q5;
t8 = cos(t7);
t9 = q2+q3;
t10 = cos(t9);
t11 = q3+q4;
t12 = cos(t11);
t13 = t12.*(3.0./1.25e2);
t14 = q4+q5;
t15 = cos(t14);
t16 = t15.*(1.0./1.25e2);
t17 = cos(q2);
t18 = cos(q3);
t19 = t18.*(1.0./2.5e1);
t20 = cos(q4);
t21 = t20.*(3.0./1.25e2);
t22 = cos(q5);
t23 = t22.*(1.0./1.25e2);
t24 = t3.*(3.0./2.5e2);
t25 = t8.*(1.0./2.5e2);
t26 = t10.*(1.0./5.0e1);
t27 = t5.*(1.0./2.5e2);
t28 = t12.*(3.0./2.5e2);
t29 = t15.*(1.0./2.5e2);
t30 = t18.*(1.0./5.0e1);
t31 = t16+t21+t23+t24+t25+t26+t27+t28+t30+9.0./5.0e1;
t32 = t20.*(3.0./2.5e2);
t33 = t22.*(1.0./2.5e2);
t34 = t17.*(7.0./2.5e2);
t35 = t6+t13+t16+t19+t21+t23+t24+t25+t26+t34+3.2e1./1.25e2;
t36 = t16+t21+t23+t27+t28+t30+9.0./5.0e1;
t37 = t23+t24+t25+t27+t28+t29+t32+1.4e1./1.25e2;
t38 = t23+t27+t28+t29+t32+1.4e1./1.25e2;
t39 = t23+t29+t32+1.4e1./1.25e2;
t40 = t25+t27+t29+t33+1.3e1./2.5e2;
t41 = t27+t29+t33+1.3e1./2.5e2;
t42 = t29+t33+1.3e1./2.5e2;
t43 = t33+1.3e1./2.5e2;
BSymbols = reshape([t3.*(3.0./1.25e2)+t6+t8.*(1.0./1.25e2)+t10.*(1.0./2.5e1)+t13+t16+t17.*(7.0./1.25e2)+t19+t21+t23+1.7e1./5.0e1,t31,t35,t37,t40,t35,t36,t6+t13+t16+t19+t21+t23+3.2e1./1.25e2,t38,t41,t31,t16+t21+t23+9.0./5.0e1,t36,t39,t42,t37,t39,t38,t23+1.4e1./1.25e2,t43,t40,t42,t41,t43,1.3e1./2.5e2],[5,5]);
