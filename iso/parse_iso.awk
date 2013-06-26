($2=="Isochrone") {
    if ( $5 < 0.001 ) {
	z = substr($5,3,4);
    } else {
	z = substr($5,3,3);
    }

    t1 = substr($8,1,4);
    t2 = substr($8,8,2);
    lt = int(100*log(t1)/log(10) + 0.5); 
    if (lt<10) {lt = "0" lt; }
    t = t2 "." lt;
    filename = "z" z "_" t;
    print filename;
}

($3=="M_ini") {}

($1 != "#") {
  printf("%13.8f %7.3f %7.3f %7.3f %7.3f\n", $2,$8,$9,$10,$12) >filename;
}
