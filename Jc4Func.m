function Jc4Symbols = Jc4Func(q1,q2,q3,q4)
%JC4FUNC
%    JC4SYMBOLS = JC4FUNC(Q1,Q2,Q3,Q4)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    28-Feb-2017 10:31:52

t2 = q1+q2+q3;
t3 = sin(t2);
t4 = q1+q2+q3+q4;
t5 = sin(t4);
t6 = q1+q2;
t7 = sin(t6);
t8 = cos(t2);
t9 = t8.*(1.0./5.0);
t10 = cos(t4);
t11 = t10.*(1.0./1.0e1);
t12 = cos(t6);
t13 = t12.*(1.0./5.0);
Jc4Symbols = reshape([t3.*(-1.0./5.0)-t5.*(1.0./1.0e1)-t7.*(1.0./5.0)-sin(q1).*(1.0./5.0),t9+t11+t13+cos(q1).*(1.0./5.0),t3.*(-1.0./5.0)-t5.*(1.0./1.0e1)-t7.*(1.0./5.0),t9+t11+t13,t3.*(-1.0./5.0)-t5.*(1.0./1.0e1),t9+t11,t5.*(-1.0./1.0e1),t11,0.0,0.0],[2,5]);
