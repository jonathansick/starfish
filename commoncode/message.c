// ** message.c: fancy console output, callable from FORTRAN.
//               Outputs a line of text to the console 
//               *without* a newline character, so the same line 
//               can be overwritten on the console.
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
//     This code is completely modular; just link it with your 
//     compiled fortran code and call it like a fortran subroutine.
//
#include <stdio.h>

extern void message_( const char* string ) {
  printf( "%s\r", string );
}

extern void realmessage_( const char* string, float* f ) {
  printf( "%s %f\r", string, *f );
}

extern void dblemessage_( const char* string, double* f ) {
  printf( "%s %f\r", string, *f );
}
