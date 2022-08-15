data test;
x=1; output;
x=2; output;
x=3; output;
run;

data out.test;
set test;
run;

data out.test2;
set test;
run;

data out.test3;
set test;
run;

data out.test4;
set test;
run;

data out.test5;
set test;
run;

proc sql;
delete from out.test2 where x=2;
quit;

proc sql;
delete from out.test3 where x>1;
quit;

proc sql;
delete from out.test4 where x=1;
quit;

proc sql;
delete from out.test5 where x<3;
quit;
