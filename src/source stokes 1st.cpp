/* CFD I Project 2 Code 
   Version Number: 4
   By Seyed MohammadAmin Taleghani
   Sharif University of Technology - Aerospace Engineering Department */

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

double tstarmax = 1.5;	/* Maximum dimensionless time. */
double err = 1e-5;		/* defines the convergence critetion for the solution. */

#define imax 101		/* Number of the nodes in y-direction from y* = 0 to y* = 1.   */
#define jmax 12000		/* Max Number of the time steps in time from t* = 0 to t* = tsmax for Re = 20 and r = 0.25 (jmax = 1.5*(imax-1)*(imax-1)/(Re*r) +1 = 3001)	*/

// Defines arrays for the family of theta methods ( each array is indexed from 1 to imax and tmax; 0 indexes are disregarded ).
double u_anal [ imax +1 ][ jmax +1 ] , u_exp [ imax +1 ][ jmax +1 ], u_crank [ imax +1 ][ jmax +1 ], u_imp [ imax +1 ][ jmax +1 ];
double diff [ imax-2 +1 ], diff2 [ imax +1 ] , L_infty [jmax];

int iter_anal , iter_exp , iter_crank , iter_imp ;

// Defines non-dimensional space and time (t* and y*)
double t [ jmax +1 ], y [ imax +1 ];


/* Function Declarations */

double maxvec(double arr[], int n) 
{ 
    int i; 
      
    // Initialize maximum element 
    double max = arr[1]; 
  
    // Traverse array elements  
    // from second and compare 
    // every element with current max  
    for (i = 2; i <= n; i++) 
        if (arr[i] > max) 
            max = arr[i]; 
  
    return max; 
} 

double maxarr(double arr[][imax], int n , int m) 
{ 
    int i,j; 
      
    // Initialize maximum element 
    double max = arr[1][1]; 
  
    // Traverse array elements  
    // from second and compare 
    // every element with current max  
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= m; j++)
		{
        if (arr[i][j] > max) 
			{
				max = arr[i][j];
			}
		}
	}
    return max;
}

double * thomas (double A[][imax-1], double B[])
{

	double ratio ; static double X[imax-1];

	for ( int i = 2 ; i <= imax-2 ; i++  )
	{
		ratio = A[i][i-1] / A[i-1][i-1];
		A[i][i] = A[i][i] - ratio * A[i-1][i];
		B [i] = B[i] - ratio * B[i-1];
	}

	X [imax-2 ] = B[imax-2 ]/A[imax-2 ][imax-2 ];

	for ( int i = imax-2 -1 ; i >= 1 ; i--)
	{

		X [ i ] = ( B[i] - A[i][i+1]*X[i+1] )/A[i][i];

	}

	return X;
}

void run_exp_method (double	Re, double	r, string prefix)
{

	int	tmax = int(1.5*(imax-1)*(imax-1)/(Re*r) + 1);

//// Defines non-dimensional time steps and space steps (delta t* and delta y*)
	double dt = tstarmax/(tmax-1); double dy = 1.0/(imax-1);

int i,j;

for ( i = 1; i <= jmax ; i++)
{

	t [ i ] = (i-1)*dt;

}
for ( i = 1; i <= imax ; i++)
{

	y [ i ] = (i-1)*dy;

}

////// Initial boundary and initial conditions:////

////// Initial Conditions:////
for ( i = 2; i <= imax-1 ; i++)
{

		u_exp	  [ i ] [ 1 ] = 0.0;

}

////// Initial Conditions:////
for ( j = 1; j <= jmax ; j++)
{

		u_exp	  [ 1 ] [ j ] = 1.0;

		u_exp	  [ imax ] [ j ] = 0.0;

}

////// Main Loop //////

j = 1;
int iter_exp = 0;

while ( 1 )
{
	for ( i = 2; i <= imax-1 ; i++ )
	{

		u_exp [i] [j+1] = u_exp [i] [j] + r * ( u_exp [i-1] [j] -2.0 * u_exp [i] [j] + u_exp [i+1] [j] );

		diff [i-1] = abs(u_exp [i] [j+1] - u_exp [i] [j]);

	}

	iter_exp++;	

	L_infty[iter_exp] = maxvec(diff,imax-2);

	if (L_infty[iter_exp]<err | j == jmax-1)
		{

			break;
		}

	j++;

}
printf("Numerical solution using the Explicit method with Re = %1.2f and r = %1.2f has", Re , r );
if ( iter_exp == jmax-1 )
{
	printf(" diverged.\n\n");
}
else 
{
	printf("finished. Number of iterations for this method: %d\n\n", iter_exp);
}
	


FILE *dat_exp;

std::string str,Line2 ;
Line2 = "dat_exp.dat" ;
str = prefix + Line2  ;
const char * TotalName = str.c_str();

/* Tecplot Output for The Explicit Method*/
/* Header */

dat_exp=fopen(TotalName, "w");
fprintf ( dat_exp,"Variables = " );
fprintf ( dat_exp, "u<sup>*</sup>");
fprintf ( dat_exp, "," );
fprintf ( dat_exp, "y<sup>*</sup>");
fprintf ( dat_exp, "\n");
/* Data */
fprintf ( dat_exp, "Zone T = \"Explicit Re = %2.0f\" , I = %d", Re, imax );
fprintf ( dat_exp, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_exp, "%e", u_exp [i][tmax]);
		fprintf (dat_exp, "		" );
		fprintf (dat_exp, "%e", y [i]);
		fprintf (dat_exp, "\n");
	}

fclose( dat_exp );


Line2 = "residue_exp.dat" ;
str = prefix + Line2  ;
TotalName = str.c_str() ;

FILE *residue_exp;

/* Tecplot Residue Outputs for the Explicit method*/
/* Header */
residue_exp=fopen(TotalName,"w");
fprintf ( residue_exp,"Variables = " );
fprintf ( residue_exp, "\"Number of Iterations\"");
fprintf ( residue_exp, "," );
fprintf ( residue_exp, "Log<sub>10</sub>(L<sub><math>5</math></sub>)");
fprintf ( residue_exp, "\n");
fprintf ( residue_exp, "Zone T=\"r = %1.2f Re = %2.0f\" , I = %d" , r ,Re , iter_exp );
fprintf ( residue_exp, "\n" );
/* Data */
for ( i = 1; i <= iter_exp; i++ )
	{
		fprintf (residue_exp, "%d", i);
		fprintf (residue_exp, "		" );
		fprintf (residue_exp, "%e", log10(L_infty[i]) );
		fprintf (residue_exp, "\n");
	}



FILE *dat_anim;
int ratio = (tmax-1)/(1.5*(imax-1)*(imax-1)/(100*5));

Line2 = "anim_exp.dat" ;   
str = prefix + Line2  ;
TotalName = str.c_str() ;

/* Tecplot Output for animation */
/* Header */
dat_anim=fopen(TotalName,"w");
fprintf ( dat_anim,"Variables = " );
fprintf ( dat_anim, "y<sup>*</sup>");
fprintf ( dat_anim, "," );
fprintf ( dat_anim, "u<sup>*</sup>");
fprintf ( dat_anim, "\n");
/* Data */

for ( j = 1; j <= jmax; j = j + ratio*1 )
{
fprintf ( dat_anim, "Zone T =\"t<sup>*</sup> = %1.2f\" , I = %d" , t [j], imax );
fprintf ( dat_anim, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_anim, "%e", y [i]);
		fprintf (dat_anim, "		" );
		fprintf (dat_anim, "%e", u_exp [i][j]);
		fprintf (dat_anim, "\n" );
	}

fprintf ( dat_anim, "\n" );
fprintf ( dat_anim, "\n" );

}

}


void run_crank_method (double	Re, double	r, string prefix)
{

	int	tmax = int (1.5*(imax-1)*(imax-1)/(Re*r) + 1 ) ;

////// Defines non-dimensional time steps and space steps (delta t* and delta y*)
	double dt = tstarmax/(tmax-1); double dy = 1/(imax-1);

////// Initial boundary and initial conditions:////

	int i,j,k;

////// Initial Conditions:////
for ( i = 2; i <= imax-1 ; i++)
{

	u_crank	  [ i ] [ 1 ] = 0.0;

}

////// Initial Conditions:////
for ( j = 1; j <= jmax ; j++)
{

		u_crank	  [ 1 ] [ j ] = 1.0;

		u_crank	  [ imax ] [ j ] = 0.0;

}


////// Defines A, B and X matrices for A*X=B
	double A [imax-2 +1][imax-2 +1] , B [imax-2 +1];
	double	*X;


	j = 1;
	iter_crank = 0;

	while (1)
	{



	for ( i = 1; i <= imax - 2; i++)
	{
	for ( k = 1; k <= imax - 2; k++)
		{
			A [ i ] [ k ] = 0.0;
		}

	}

	/// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop)

	A [1][1] = 1.0 + r ;
	A [1][2] = -0.5*r ;
	B [1]    = 0.5*r*u_crank [1][j] + (1.0 - r)* u_crank [2][j] + 0.5*r*u_crank [3][j] + 0.5*r*u_crank [1][j+1] ;

	for ( i = 2 ; i <= imax-2 -1 ; i++)
	{

		A [i][i-1] = -0.5*r ;
		A [i][i]   =  1.0 + r ;
		A [i][i+1] = -0.5*r ;
		B [i]      =  0.5*r*u_crank [i][j] + (1.0 - r)* u_crank [i+1][j] + 0.5*r*u_crank [i+2][j] ;

	}

	A [imax-2 ][imax-2 -1 ] = -0.5*r ;
	A [imax-2 ][imax-2 ]	= 1.0 + r ;
	B [imax-2 ]				= 0.5*r*u_crank [imax-2][j] + (1.0 - r)* u_crank [imax-1][j] + 0.5*r*u_crank [imax][j] + 0.5*r*u_crank [imax][j+1] ;

		X = thomas (A,B);

		for ( i = 1 ; i <= imax-2 ; i++)
		{
			u_crank [i+1][j+1] = *(X+i);
			diff [i] = abs(u_crank [i+1] [j+1] - u_crank [i+1] [j]);
		}

		iter_crank++;

		L_infty[iter_crank] = maxvec(diff,imax-2);

	if (L_infty[iter_crank]<err)
		{

			break;

		}

	j++;

	}
	printf("Numerical solution using the Crank-Nicholson method with Re = %1.2f and r = %1.2f has finished. Number of iterations for this method: %d\n\n", Re , r, iter_crank);


	FILE *dat_crank;

	std::string str,Line2 ;
	Line2 = "dat_crank.dat" ;
	str = prefix + Line2  ;
	const char * TotalName = str.c_str();


	/* Tecplot Output for The Analytical Solution*/
	/* Header */
	dat_crank=fopen(TotalName , "w");
	fprintf ( dat_crank,"Variables = " );
	fprintf ( dat_crank, "u<sup>*</sup>");
	fprintf ( dat_crank, "," );
	fprintf ( dat_crank, "y<sup>*</sup>");
	fprintf ( dat_crank, "\n");
	/* Data */
	fprintf ( dat_crank, "Zone T = \"Crank-Nicholson Re = %2.0f\" , I = %d", Re, imax );
	fprintf ( dat_crank, "\n" );
		for ( i = 1; i <= imax; i++ )
		{
			fprintf (dat_crank, "%e", u_crank [i][tmax]);
			fprintf (dat_crank, "		" );
			fprintf (dat_crank, "%e", y [i]);
			fprintf (dat_crank, "\n");
		}

	fclose( dat_crank );


	Line2 = "residue_crank.dat" ;
	str = prefix + Line2  ;
	TotalName = str.c_str() ;

	FILE *residue_crank;

	/* Tecplot Residue Outputs for Jacubi*/
	/* Header */
	residue_crank=fopen(TotalName,"w");
	fprintf ( residue_crank,"Variables = " );
	fprintf ( residue_crank, "\"Number of Iterations\"");
	fprintf ( residue_crank, "," );
	fprintf ( residue_crank, "Log<sub>10</sub>(L<sub><math>5</math></sub>)");
	fprintf ( residue_crank, "\n");
	fprintf ( residue_crank, "Zone T=\"r = %1.2f Re = %2.0f\" , I = %d" , r ,Re , iter_crank);
	fprintf ( residue_crank, "\n" );
	/* Data */
		for ( i = 1; i <= iter_crank; i++ )
		{
			fprintf (residue_crank, "%d", i);
			fprintf (residue_crank, " " );
			fprintf (residue_crank, "%e", log10(L_infty[i]) );
			fprintf (residue_crank, "\n");
		}



	FILE *dat_anim;
	int ratio = (tmax-1)/(1.5*(imax-1)*(imax-1)/(100*5));

	Line2 = "anim_crank.dat" ;   
	str = prefix + Line2  ;
	TotalName = str.c_str() ;

	/* Tecplot Output for animation */
	/* Header */
	dat_anim=fopen(TotalName,"w");
	fprintf ( dat_anim,"Variables = " );
	fprintf ( dat_anim, "y<sup>*</sup>");
	fprintf ( dat_anim, "," );
	fprintf ( dat_anim, "u<sup>*</sup>");
	fprintf ( dat_anim, "\n");
	/* Data */

	for ( j = 1; j <= jmax; j = j + ratio*1 )
	{
		fprintf ( dat_anim, "Zone T =\"t<sup>*</sup> = %1.2f\" , I = %d" , t [j], imax );
		fprintf ( dat_anim, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_anim, "%e", y [i]);
		fprintf (dat_anim, "		" );
		fprintf (dat_anim, "%e", u_crank [i][j]);
		fprintf ( dat_anim, "\n" );
	}

	fprintf ( dat_anim, "\n" );
	fprintf ( dat_anim, "\n" );

}


}


void run_imp_method (double	Re, double	r, string prefix)
{

	int	tmax = int (1.5*(imax-1)*(imax-1)/(Re*r) + 1 ) ;

////// Defines non-dimensional time steps and space steps (delta t* and delta y*)
	double dt = tstarmax/(tmax-1); double dy = 1/(imax-1);

////// Initial boundary and initial conditions:////

	int i,j,k;

////// Initial Conditions:////
for ( i = 2; i <= imax-1 ; i++)
{

		u_imp	  [ i ] [ 1 ] = 0.0;

}

////// Initial Conditions:////
for ( j = 1; j <= jmax ; j++)
{

		u_imp	  [ 1 ] [ j ] = 1.0;

		u_imp	  [ imax ] [ j ] = 0.0;

}


//////// Defines A, B and X matrices for A*X=B;
	double A [imax-2 +1][imax-2 +1] , B [imax-2 +1];
	double	*X;


 j = 1;
 iter_imp = 0;

while ( 1 )
{

	for ( i = 1; i <= imax - 2; i++)
	{
	for ( k = 1; k <= imax - 2; k++)
		{
			A [ i ] [ k ] = 0.0;
		}
	}

//////////// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop);

	A [1][1] = 1.0 + 2.0*r ;
	A [1][2] = -r ;
	B [1]    = u_imp [2][j] + r*u_imp [1][j+1] ;

	for ( i = 2 ; i <= imax-2 -1 ; i++)
	{

		A [i][i-1] = -r ;
		A [i][i]   = 1.0 + 2.0*r ;
		A [i][i+1] = -r ;
		B [i]      = u_imp [i+1][j] ;

	}

	A [imax-2 ][imax-2 -1 ] = -r ;
	A [imax-2 ][imax-2 ] = 1.0 + 2.0*r ;
	B [imax-2 ]			 = u_imp [imax-1][j] + r*u_imp [imax][j+1] ;

		X = thomas (A,B);

		for ( i = 1 ; i <= imax-2 ; i++)
		{
			u_imp [i+1][j+1] = *(X+i);
			diff [i] = abs(u_imp [i+1] [j+1] - u_imp [i+1] [j]);
		}

	iter_imp++;

	L_infty[iter_imp] = maxvec(diff,imax-2);

	if (L_infty[iter_imp]<err)
		{

			break;

		}

	j++;

}
printf("Numerical solution using the Implicit method with Re = %1.2f and r = %1.2f has finished. Number of iterations for this method: %d\n\n", Re , r, iter_imp);


	std::string str,Line2 ;
	Line2 = "dat_imp.dat" ;
	str = prefix + Line2  ;
	const char * TotalName = str.c_str();

	FILE *dat_imp;

/* Tecplot Output for The Implicit Method*/
/* Header */
dat_imp=fopen(TotalName , "w");
fprintf ( dat_imp,"Variables = " );
fprintf ( dat_imp, "u<sup>*</sup>");
fprintf ( dat_imp, "," );
fprintf ( dat_imp, "y<sup>*</sup>");
fprintf ( dat_imp, "\n");
//* Data */
fprintf ( dat_imp, "Zone T = \"Implicit Re = %2.0f\" , I = %d", Re, imax );
fprintf ( dat_imp, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_imp, "%e", u_imp [i][tmax]);
		fprintf (dat_imp, "			" );
		fprintf (dat_imp, "%e", y [i]);
		fprintf (dat_imp, "\n");
	}

fclose( dat_imp );


Line2 = "residue_imp.dat" ;
str = prefix + Line2  ;
TotalName = str.c_str() ;

FILE *residue_imp;

/* Tecplot Residue Outputs for Implicit Method*/
/* Header */
residue_imp=fopen(TotalName,"w");
fprintf ( residue_imp,"Variables = " );
fprintf ( residue_imp, "\"Number of Iterations\"");
fprintf ( residue_imp, "," );
fprintf ( residue_imp, "Log<sub>10</sub>(L<sub><math>5</math></sub>)");
fprintf ( residue_imp, "\n");
fprintf ( residue_imp, "Zone T=\"r = %1.2f Re = %2.0f\" , I = %d" , r ,Re , iter_imp);
fprintf ( residue_imp, "\n" );
/* Data */
for ( i = 1; i <= iter_imp; i++ )
	{
		fprintf (residue_imp, "%d", i);
		fprintf (residue_imp, " " );
		fprintf (residue_imp, "%e", log10(L_infty[i]) );
		fprintf (residue_imp, "\n");
	}



FILE *dat_anim;
int ratio = (tmax-1)/(1.5*(imax-1)*(imax-1)/(100*5));

Line2 = "anim_imp.dat" ;   
str = prefix + Line2  ;
TotalName = str.c_str() ;

/* Tecplot Output for animation */
/* Header */
dat_anim=fopen(TotalName,"w");
fprintf ( dat_anim,"Variables = " );
fprintf ( dat_anim, "y<sup>*</sup>");
fprintf ( dat_anim, "," );
fprintf ( dat_anim, "u<sup>*</sup>");
fprintf ( dat_anim, "\n");
/* Data */

for ( j = 1; j <= iter_imp; j = j + ratio*1 )
{
fprintf ( dat_anim, "Zone T =\"t<sup>*</sup> = %1.2f\" , I = %d" , t [j], imax );
fprintf ( dat_anim, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_anim, "%e", y [i]);
		fprintf (dat_anim, "		" );
		fprintf (dat_anim, "%e", u_imp [i][j]);
		fprintf (dat_anim, "\n" );
	}

fprintf ( dat_anim, "\n" );
fprintf ( dat_anim, "\n" );
}

}

void run_anal_method (double Re , double r , string prefix)
{

FILE *dat_anal;

///////////////////////////////////  Analytical Solution for Reynolds = Re   \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

int tmax = int(1.5*(imax-1)*(imax-1)/(Re*r) + 1);

// Defines non-dimensional time steps and space steps (delta t* and delta y*)
double dt = tstarmax/(tmax-1); double dy = 1.0/(imax-1);


int i,j,k;
for ( i = 1; i <= jmax ; i++)
{

	t [ i ] = (i-1)*dt;

}
for ( i = 1; i <= imax ; i++)
{

	y [ i ] = (i-1)*dy;

}


//////////////////////////// Calculates Analytical Solution With Re_H = Re ////////////////////////////////////


double	sigma;

j = 1;
iter_anal = 0;


		for ( i = 1 ; i <= imax ; i++)
		{

			sigma = 0.0;

			for ( k = 1 ; k <= 100 ; k++ )
			{

				sigma = sigma + 1.0/k * sin( k*M_PI*y[i] ) * exp ( -(k*M_PI)*(k*M_PI)* t[j]/Re );

			}

			u_anal [i][j] = 1.0 - y[i] -2.0/M_PI * sigma;

		}


while (1)
{

	j++;

for ( i = 1 ; i <= imax ; i++)
{

		sigma = 0.0;

		for ( k = 1 ; k <= 100 ; k++ )
		{

			sigma = sigma + 1.0/k * sin( k*M_PI*y[i] ) * exp ( -(k*M_PI)*(k*M_PI)* t[j]/Re );

		}

		u_anal [i][j] = 1.0 - y[i] -2.0/M_PI * sigma;
		diff2 [i] = abs(u_anal [i][j] - u_anal [i][j-1]);  

}

	iter_anal ++;

	L_infty [iter_anal] = maxvec(diff2,imax);

	if ( L_infty [iter_anal]<err ) 
		{
			break;
		}

}


	std::string str,Line2 ;
	Line2 = "dat_anal.dat" ;
	str = prefix + Line2  ;
	const char * TotalName = str.c_str();

/* Tecplot Output for The Analytical Solution*/
/* Header */
dat_anal=fopen(TotalName,"w");
fprintf ( dat_anal,"Variables = " );
fprintf ( dat_anal, "u<sup>*</sup>");
fprintf ( dat_anal, "," );
fprintf ( dat_anal, "y<sup>*</sup>");
fprintf ( dat_anal, "\n");
/* Data */
fprintf ( dat_anal, "Zone T = \"Exact Solution with Re = %2.0f \" , I = %d" , Re , imax );
fprintf ( dat_anal, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_anal, "%e", u_anal [i][tmax]);
		fprintf (dat_anal, "		" );
		fprintf (dat_anal, "%e", y [i]);
		fprintf (dat_anal, "\n");
	}

fclose(dat_anal);


FILE *dat_anim;
int ratio = (tmax-1)/(1.5*(imax-1)*(imax-1)/(100*5));

Line2 = "anim_anal.dat" ;   
str = prefix + Line2  ;
TotalName = str.c_str() ;

/* Tecplot Output for animation */
/* Header */
dat_anim=fopen(TotalName,"w");
fprintf ( dat_anim,"Variables = " );
fprintf ( dat_anim, "y<sup>*</sup>");
fprintf ( dat_anim, "," );
fprintf ( dat_anim, "u<sup>*</sup>");
fprintf ( dat_anim, "\n");
/* Data */

for ( j = 1; j <= iter_imp; j = j + ratio*1 )
{
fprintf ( dat_anim, "Zone T =\"t<sup>*</sup> = %1.2f\" , I = %d" , t [j], imax );
fprintf ( dat_anim, "\n" );
	for ( i = 1; i <= imax; i++ )
	{
		fprintf (dat_anim, "%e", y [i]);
		fprintf (dat_anim, "		" );
		fprintf (dat_anim, "%e", u_anal [i][j]);
		fprintf (dat_anim, "\n" );
	}

fprintf ( dat_anim, "\n" );
fprintf ( dat_anim, "\n" );
}

}


void main ()
{

	string prefix;
	double Re,r;

	Re = 20.0;

	r = 0.25;
	prefix = "20_0.25_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);


	r = 0.5;
	prefix = "20_0.5_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);

	r = 1.0;
	prefix = "20_1.0_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);

	r = 2.0;
	prefix = "20_2.0_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);

	r = 5.0;
	prefix = "20_5.0_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);



	Re = 100.0;

	r = 0.25;
	prefix = "100_0.25_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);


	r = 0.5;
	prefix = "100_0.5_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);

	r = 1.0;
	prefix = "100_1.0_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);

	r = 2.0;
	prefix = "100_2.0_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);

	r = 5.0;
	prefix = "100_5.0_";
	run_exp_method (Re , r , prefix);
	run_crank_method (Re , r , prefix);
	run_imp_method (Re , r , prefix);
	run_anal_method (Re , r , prefix);








printf("Press Enter to terminate the program.\n");

cin.get();

}
