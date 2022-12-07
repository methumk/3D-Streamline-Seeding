/*

Matrix module.

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

#include <math.h>
#include <stdio.h>
#include "tmatrix.h"

#ifndef _NON_WINDOWS_
#define M_PI 3.141915
#endif

#define  MAT_MAX  20          /*  maximum number of matrices on stack  */


static float mat_stack [MAT_MAX] [4] [4];    /* matrix stack */
static int mat_num;                          /* top of stack */

Matrix mat_form;                   /* current transformation matrix */


/******************************************************************************
Set a matrix to the identity matrix.
******************************************************************************/

void mat_ident( Matrix m)
{
  int i;

  for (i = 0; i <= 3; i++) {
    m [i] [0] = 0.0;
    m [i] [1] = 0.0;
    m [i] [2] = 0.0;
    m [i] [3] = 0.0;
    m [i] [i] = 1.0;
  }
}


/******************************************************************************
Copy a matrix.
******************************************************************************/

void mat_copy (Matrix dest, Matrix source)
{
  int i,j;

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      dest[i][j] = source[i][j];
}


/******************************************************************************
Multiplies two 4 by 4 matrices.

prod = a * b;
******************************************************************************/

void mat_mult (Matrix prod, Matrix a, Matrix b)
{
  int i,j;
  Matrix m;

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      m [i] [j] = a [0] [j] * b [i] [0] + a [1] [j] * b [i] [1] +
                  a [2] [j] * b [i] [2] + a [3] [j] * b [i] [3];

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      prod [i] [j] = m [i] [j];
}


/******************************************************************************
Transpose a matrix.
******************************************************************************/

void mat_transpose (Matrix m)
{
  int i,j;
  float t;

  for (i = 1; i <= 3; i++)
    for (j = 0; j < i; j++) {
      t = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = t;
    }
}


/******************************************************************************
Create translation matrix.
******************************************************************************/

void mat_translate (Matrix m, float x, float y, float z)
{
  mat_ident (m);
  m[3][0] = x;
  m[3][1] = y;
  m[3][2] = z;
}


/******************************************************************************
Create scaling matrix.
******************************************************************************/

void mat_scale (Matrix m, float x, float y, float z)
{
  mat_ident (m);
  m[0][0] = x;
  m[1][1] = y;
  m[2][2] = z;
}


/******************************************************************************
Creat a rotation matrix.

Entry:
  angle - angle of rotation in degrees.
  axis  - axis to rotate about: 'x', 'y' or 'z'

Exit:
  mat - rotation matrix
******************************************************************************/

void mat_rotate (Matrix mat, float angle, char axis)
{
  int a,b;
  double sine,cosine;

  mat_ident (mat);

  switch (axis) {
    case 'x' : a = 1;  b = 2;  break;
    case 'y' : a = 2;  b = 0;  break;
    case 'z' : a = 0;  b = 1;  break;
    default : printf ("bad rotation axis: '%c'\n", axis);
              return; break;
  }

  cosine = cos (M_PI / 180 * angle);
  sine   = sin (M_PI / 180 * angle);

  mat[a][a] = cosine;
  mat[a][b] = sine;
  mat[b][a] = -sine;
  mat[b][b] = cosine;
}


/******************************************************************************
Calculate matrix that rotates a given unit vector to the x-axis.
******************************************************************************/



/******************************************************************************
Print out a matrix.
******************************************************************************/

void mat_print (Matrix m)
{
  int i,j;

  printf ("\n");

  for (j = 0; j <= 3; j++) {
    for (i = 0; i <= 3; i++)
        printf ("%f  ", m[i][j]);
    printf ("\n");
  }

  printf ("\n");
}


/******************************************************************************
Transform a vector by a matrix.
******************************************************************************/

void mat_apply (Matrix m, Vector v)
{
  int j;
  Vector t;

  for (j = 0; j <= 2; j++)
    t [j] = v[0] * m[0][j] + v[1] * m[1][j] + v[2] * m[2][j] + m[3][j];

  v[0] = t[0];
  v[1] = t[1];
  v[2] = t[2];
}


/******************************************************************************
Transform a plane by a matrix.
******************************************************************************/

void mat_apply_plane (Matrix m, Plane p)
{
  int j;
  Plane t;

  for (j = 0; j <= 3; j++)
    t [j] = p[0] * m[0][j] + p[1] * m[1][j] + p[2] * m[2][j] + p[3] * m[3][j];

  p[0] = t[0];
  p[1] = t[1];
  p[2] = t[2];
  p[3] = t[3];
}


/******************************************************************************
Set 'mat_form' to the identity matrix.
******************************************************************************/

void identity()
{
  mat_ident (mat_form);
}


/******************************************************************************
Post-multiply 'mat_form' by translation.
******************************************************************************/

void translate (float x, float y, float z)
{
  Matrix n;

  mat_translate (n, x, y, z);
  mat_mult (mat_form, mat_form, n);
}


/******************************************************************************
Post-multiply 'mat_form' by scale.
******************************************************************************/

void scale (float x, float y, float z)
{
  Matrix n;

  mat_scale (n, x, y, z);
  mat_mult (mat_form, mat_form, n);
}


/******************************************************************************
Post-multiply 'mat_form' by rotation.

Entry:
  angle - angle of rotation in degrees.
  axis  - axis to rotate about: 'x', 'y' or 'z'
******************************************************************************/

void rotate (float angle, char axis)
{
  Matrix n;

  mat_rotate (n, angle, axis);
  mat_mult (mat_form, mat_form, n);
}


/******************************************************************************
Initialize the matrix stack. "mat_num" is the number of matrices on the 
stack,
and always points to the next availible matrix.
******************************************************************************/

void init_matrices()
{
  mat_num = 0;
  mat_ident (mat_form);
}


/******************************************************************************
Pushes the current matrix onto the stack.
******************************************************************************/

void push()
{
  int i,j;

  if (mat_num > MAT_MAX)
    printf ("Error in push : stack overflow. No action taken.\n");
  else {
    for (i = 0; i <= 3; i++)
      for (j = 0; j <= 3; j++)
        mat_stack [mat_num] [i] [j] = mat_form [i] [j];
    mat_num++;
  }
}


/******************************************************************************
Pops the top matrix off the stack, placing it into the current matrix.
******************************************************************************/

void pop()
{
  int i,j;

  if (mat_num == 0)
    printf ("Error in pop : stack underflow. No action taken.\n");
  else {
    mat_num--;
    for (i = 0; i <= 3; i++)
      for (j = 0; j <= 3; j++)
        mat_form [i] [j] = mat_stack [mat_num] [i] [j];
  }
}


/******************************************************************************
Return a copy of the current transformation matrix.
******************************************************************************/

void get_transformation (Matrix mat)
{
  int i,j;

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      mat [i] [j] = mat_form [i] [j];
}


/******************************************************************************
Transform a vector by the current transformation matrix.
******************************************************************************/

void vtransform (Vector v)
{
  int j;
  Vector t;

  for (j = 0; j <= 2; j++)
    t [j] = v [0] * mat_form [0] [j] + v [1] * mat_form [1] [j] +
            v [2] * mat_form [2] [j] + mat_form [3] [j];

    v [0] = t [0];
    v [1] = t [1];
    v [2] = t [2];
}

