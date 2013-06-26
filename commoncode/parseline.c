// *** parseline.c: File I/O functions
//
//     This file is part of StarFISH
//     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
//     
// -=> StarFISH is free software; you can redistribute it 
// -=> and/or modify it under the terms of the GNU General Public 
// -=> License (GPL) as published by the Free Software Foundation;
// -=> either version 2 of the License, or (at your option) any 
// -=> later version.
//
// File I/O functions that are easy to write in C, but impossible in 
// FORTRAN.  These functions are callable from FORTRAN programs.
//
// Every "const char*" arg requires a "int" arg at the end of the list
// this isn't used in the function, but is needed to export the function to 
// fortran.
//
#include <stdio.h>
#include <string.h>
extern void getreal_( const char* string, float* f, int n ) {
  sscanf( string, "%f", f );
}

extern void getint_( const char* string, int *i, int n ) {
  sscanf( string, "%d", i );
}

extern void getword_( const char* instring, char* outstring, int *outlen, int n1, int n2 ) {
  //Truncate at first instance of a space character.
  //(replace remaining chars with ' ')
  unsigned int i;
  unsigned int j;
  int iStart = -1;

  strncpy( outstring, instring, *outlen );

  for ( i = 0; i < *outlen; i++ ) {
    if ( iStart == -1 ) {
      if ( instring[i] != ' ' ) {
	iStart = i;
	if ( instring[i] != 0 ) 
	  outstring[i-iStart] = instring[i]; //skip leading whitespace
	else 
	  break;
      }
    } else {
      if ( instring[i] != ' ' && instring != 0 ) 
	outstring[i-iStart] = instring[i]; //skip leading whitespace
      else 
	break;
    }
  }

  for ( j = i-iStart; j < *outlen; j++ ) { outstring[j] = ' '; }
}
