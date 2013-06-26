BEGIN {
  while (getline < grid >0) {
    ibox[$1,$2] = $3;
  }
}

($1==cmd) {
  pop[$2] = $4 - $3 + p0;
}

END {
  for (iy=ny; iy>0; iy--) {
    for (ix=1; ix<=nx; ix++) {
      print pop[ibox[ix,iy]] > out;
      if (pop[ibox[ix,iy]]=="") { print "null grid box: ", ix,iy; }
    }
  }
}

