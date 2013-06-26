BEGIN {
  out = pre ".input";
  norm = 400000;  # The number of stars to attempt per unit amplitude
  printf("%2s        output file prefix\n", pre) >>out;
  printf("      1     lockflag\n") >>out;
  printf("%7d     sfr normalization\n", norm) >>out;
  printf("   0.00     delta distance modulus\n") >>out;
  printf("   1.00     E(B-V) factor\n") >>out;
  printf("  -1.35     IMF slope\n") >>out;
  printf("   0.50     binary fraction\n") >>out;

# Make sure these numbers describe the number of isochrones
# in each of your isochrone groups (the order needs to match 
# your isochrone-lock file; if you aren't locking isochrones, 
# then you should still provide values of 1 for each isochrone).
  niso[0]  = 3;  #  log(age) = 10.00;  Z = 0.008
  niso[1]  = 4;  # 9.80
  niso[2]  = 4;  # 9.60
  niso[3]  = 4;  # 9.40
  niso[4]  = 4;  # 9.20
  niso[5]  = 4;  # 9.00
  niso[6]  = 4;  # 8.80
  niso[7]  = 4;  # 8.60
  niso[8]  = 4;  # 8.40
  niso[9]  = 4;  # 8.20
  niso[10] = 4;  # 8.00
  niso[11] = 4;  # 7.80
  niso[12] = 4;  # 7.60
  niso[13] = 4;  # 7.40
  niso[14] = 4;  # 7.20
  niso[15] = 4;  # 7.00
  niso[16] = 4;  # 6.80
  niso[17] = 3;  # 6.60
  niso[18] = 3;  # log(age) = 10.00;  Z = 0.004
  niso[19] = 4;  # 9.80
  niso[20] = 4;  # 9.60
  niso[21] = 4;  # 9.40
  niso[22] = 4;  # 9.20
  niso[23] = 4;  # 9.00
  niso[24] = 4;  # 8.80
  niso[25] = 4;  # 8.60
  niso[26] = 4;  # 8.40
  niso[27] = 4;  # 8.20
  niso[28] = 4;  # 8.00
  niso[29] = 4;  # 7.80
  niso[30] = 4;  # 7.60
  niso[31] = 4;  # 7.40
  niso[32] = 4;  # 7.20
  niso[33] = 4;  # 7.00
  niso[34] = 4;  # 6.80
  niso[35] = 3;  # 6.60
  niso[36] = 3;  # log(age) = 10.00;  Z = 0.001
  niso[37] = 4;  # 9.80
  niso[38] = 4;  # 9.60
  niso[39] = 4;  # 9.40
  niso[40] = 4;  # 9.20
  niso[41] = 4;  # 9.00
  niso[42] = 4;  # 8.80
  niso[43] = 4;  # 8.60
  niso[44] = 4;  # 8.40
  niso[45] = 4;  # 8.20
  niso[46] = 3;  # 8.00
}

(1) {
  #Modify these metallicity values and filename prefixes 
  #to match your isochrone set 
  if($1 == 0.001) { 
    filestr = "z001_"; 
  } else if($1 == 0.004) {
    filestr = "z004_";
  } else if($1 == 0.008) {
    filestr = "z008_";
  }

  # prefix a "0" on the age if it is less than 10.0
  if ($2 < 10.0) {
    filestr = filestr "0" $2;
  } else {
    filestr = filestr $2;
  }

  printf("%9.4f  %6.4f  %5.2f  %2d\n", $3/norm, $1, $2, niso[NR-1]) >>out;
}
