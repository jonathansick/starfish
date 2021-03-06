BEGIN {
  linecount = 0;
  lo = 0;
  hi = 255;
  norm = 255/max;
  i = 0;
}

{
  if (clr==1) {
    p[0,i] = $1*norm;
    p[1,i] = $2*norm;
    p[2,i++] = $3*norm;
  } else {
    p[0,i] = $1*norm;
    p[1,i] = $1*norm;
    p[2,i++] = $1*norm;
  }
}

END {
  imax = i;

  print "P3" >>out;
  print "# CREATOR: pxl2ppm 1.0" >>out;
  printf("%d %d\n", xd, yd) >>out;
  print 255 >>out;

  for (i=0; i<imax; i++) {
    for (j=0; j<3; j++) {
      if (p[j,i] < lo) {
	p[j,i] = lo;
      } else if (p[j,i] > hi) {
	p[j,i] = hi;
      }

      printf("%3d ", p[j,i]) >>out;
    }
    linecount++;
    if (linecount >4) { printf("\n") >>out; linecount = 0; }
  }  
}
