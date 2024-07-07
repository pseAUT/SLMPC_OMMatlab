model tank

parameter Real Cv=0.5;


parameter Real A1=1;
parameter Real A2=1;
parameter Real A3=1;

parameter Real h1_In=0;
parameter Real h2_In=0;
parameter Real h3_In=0;

parameter Real eps=1e-6;


output Real qo1;
output Real qo2;
output Real qo3;



input Real qi1;
input Real qi2;
input Real qi3;

Real h1;
Real h2;
Real h3;


initial equation
h1=h1_In;
h2=h2_In;
h3=h3_In;

//Flimit=Flimit_In;


equation
qo1=Cv*(h1-h2)/((abs(h1-h2)+eps)^0.5);
qo2=Cv*(h2-h3)/((abs(h2-h3)+eps)^0.5);
qo3=Cv*h3/((abs(h3)+eps)^0.5);


der(h1)=(qi1-qo1)/A1;
der(h2)=(qi2+qo1-qo2)/A2  ;
der(h3)=(qi3+qo2-qo3)/A3  ; 



end tank;
