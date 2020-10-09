/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_stat - statistics.
 */

#include <zm/zm_stat.h>
#include <zm/zm_sf.h>

/* zPermut
 * - permutation.
 */
double zPermut(int n, int i)
{
  double result = 1.0;

  if( i > n || i < 0 ) return 0;
  for( ; i>0; i--, n-- ) result *= n;
  return result;
}

/* zFacto
 * - factorial.
 */
double zFacto(int n)
{
  return zPermut( n, n );
}

/* zCombi
 * - combination.
 */
double zCombi(int n, int i)
{
  if( n-i < i ) i = n-i;
  return zPermut( n, i ) / zFacto( i );
}

/* zCombiSeries
 * - series of combination.
 */
double *zCombiSeries(int n, size_t size, double c[])
{
  register int i, j;

  if( n < 0 || n >= size ){
    ZRUNERROR( ZM_ERR_STAT_ILLS );
    return NULL;
  }
  for( i=0; i<=n; i++ ){
    c[0] = c[i] = 1.0;
    for( j=i-1; j>0; j-- )
      c[j] += c[j-1];
  }
  return c;
}

/* basic distribution functions */

/* zNormalDistrib
 * - normal distribution.
 */
double zNormalDistrib(double x, double mu, double sigma)
{
  double z;

  z = ( x - mu ) / sigma;
  return zND / sigma * exp( -0.5 * zSqr(z) );
}

/* zNormalCumDistrib
 * - normal cumulative distribution.
 *   requires Gauss's error function (erf).
 */
double zNormalCumDistrib(double x, double mu, double sigma)
{
  return 0.5 * ( 1 + zErf( ( x - mu ) / ( sqrt(2) * sigma ) ) );
}

/* zPoissonDistrib
 * - Poisson distribution.
 */
double zPoissonDistrib(int x, double lambda)
{
  return exp( -lambda ) * pow( lambda, x ) / zFacto( x );
}

/* zBinDistrib
 * - binomial distribution.
 */
double zBinDistrib(int x, int n, double p)
{
  return zCombi( n, x ) * pow( p, x ) * pow( 1-p, n-x );
}

/* basic statistics computation */

/* zDataMax
 * - maximum component of data.
 */
double zDataMax(double *data, int num, int *im)
{
  register int i;
  double max;

  if( im ) *im = 0;
  for( max = data[0], i=1; i<num; i++ )
    if( data[i] > max ){
      max = data[i];
      if( im ) *im = i;
    }
  return max;
}

/* zDataMin
 * - minimum component of data.
 */
double zDataMin(double *data, int num, int *im)
{
  register int i;
  double min;

  if( im ) *im = 0;
  for( min = data[0], i=1; i<num; i++ )
    if( data[i] < min ){
      min = data[i];
      if( im ) *im = i;
    }
  return min;
}

/* zDataAbsMax
 * - maximum absolute component of data.
 */
double zDataAbsMax(double *data, int num, int *im)
{
  register int i;
  double val, max;

  if( im ) *im = 0;
  for( max = fabs( data[0] ), i=1; i<num; i++ )
    if( ( val = fabs( data[i] ) ) > max ){
      max = val;
      if( im ) *im = i;
    }
  return max;
}

/* zDataAbsMin
 * - minimum component of data.
 */
double zDataAbsMin(double *data, int num, int *im)
{
  register int i;
  double val, min;

  if( im ) *im = 0;
  for( min = fabs( data[0] ), i=1; i<num; i++ )
    if( ( val = fabs( data[i] ) ) < min ){
      min = val;
      if( im ) *im = i;
    }
  return min;
}

/* zDataSum
 * - sum up all data values.
 */
double zDataSum(double *data, int num)
{
  double s=0, s_prev=0, q=0, r;
  register int i;

  for( i=0; i<num; i++ ){
    s = s_prev + data[i];
    r = s - s_prev;
    q += data[i] - r;
    s_prev = s;
  }
  return s + q;
}

/* zDataAve
 * - calculate the average of data.
 */
double zDataAve(double *data, int num)
{
  return zDataSum(data,num) / num;
}

/* zDataVar
 * - calculate the variance of data.
 */
double zDataVar(double *data, int num)
{
  register int i;
  double ave, result;

  ave = zDataAve( data, num );
  for( result=0, i=0; i<num; i++ )
    result += zSqr( data[i] - ave );
  return result / num;
}

/* zDataSD
 * - calculate the standard deviation of data.
 */
double zDataSD(double *data, int num)
{
  return sqrt( zDataVar( data, num ) );
}
