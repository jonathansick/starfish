BEGIN {
    stage[0]  = "TO";
    stage[1]  = "Te-m";
    stage[2]  = "Te-M";
    stage[3]  = "L-M";
    stage[4]  = "RGBb";
    stage[5]  = "L-m";
    stage[6]  = "RGBt";
    stage[7]  = "BHeb";
    stage[8]  = "EHeb";
    stage[9]  = "1TP";
    stage[10] = "AGBt";
    stage[11] = "Cb";
    nstages = 12;
}
   
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
    filename = "z" z "_" t ".eep";
    print filename;
}

($3=="M_ini") {}

($1 != "#") {
    imatch = -1;
    for (i=0; i<nstages; i++) {
	if ($11==stage[i]) {
	    imatch = i;
	    break;
	}
    }
    
    if (imatch==-1) {
	print "Warning!  Undefined EEP: " $11;
	stop;
    }
    
    printf("%13.8f %7.3f %7.3f %7.3f %7.3f %3d\n", $2,$8+$9+$7,$9+$7,$7,$7-$10,imatch+1) >filename;
}
