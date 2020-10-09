/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS interpolation.
 */

#include <zm/zm_nurbs.h>

static int _zNURBSSeg(zNURBS *nurbs, double t);
static double _zNURBSBasis(zNURBS *nurbs, double t, int i, int r, int k);
static double _zNURBSBasisDiff(zNURBS *nurbs, double t, int i, int r, int k, int diff);
static double _zNURBSDenDiff(zNURBS *nurbs, double t, int s, int e, int diff);

/* zNURBSCreate
 * - create a NURBS interpolator.
 */
bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int dim)
{
  register int i, j;
  zSeqListCell *cp;

  if( zListNum(seq) <= dim ){
    ZRUNERROR( ZM_ERR_NURBS_INVDIM );
    return false;
  }
  nurbs->dim = dim;
  nurbs->knot = zVecAlloc( zListNum(seq)+dim+1 );

  zArrayAlloc( &nurbs->cparray, zNURBSCPCell, zListNum(seq) );
  if( !nurbs->knot || zArrayNum(&nurbs->cparray) == 0 ){
    ZALLOCERROR();
    zNURBSDestroy( nurbs );
    return false;
  }
  nurbs->seq = seq;
  /* initialize knot vector */
  for( i=0; i<dim; i++ )
    zNURBSKnot(nurbs,i) = 0;
  for( j=1; i<=zListNum(seq); i++, j++ )
    zNURBSKnot(nurbs,i) =  j;
  for( ; i<zVecSizeNC(nurbs->knot); i++ )
    zNURBSKnot(nurbs,i) = j;
  /* assign control points & initialize weight uniformly */
  i = 0;
  zListForEach( seq, cp ){
    zNURBSWeight(nurbs,i) = 1.0;
    zNURBSCP(nurbs,i) = cp->data.v;
    i++;
  }
  return true;
}

/* zNURBSDestroy
 * - destroy a NURBS interpolator.
 */
void zNURBSDestroy(zNURBS *nurbs)
{
  nurbs->seq = NULL;
  nurbs->dim = 0;
  zVecFree( nurbs->knot );
  nurbs->knot = NULL;
  zArrayFree( &nurbs->cparray );
}

/* zNURBSKnotNormalize
 * - normalize the knot vector of a NURBS interpolator.
 */
void zNURBSKnotNormalize(zNURBS *nurbs)
{
  zVecShift( nurbs->knot, -zNURBSKnot0(nurbs) );
  zVecDivDRC( nurbs->knot, zNURBSKnotE(nurbs) );
}

/* (static)
 * _zNURBSSeg
 * - find the knot segment of interpolation.
 */
int _zNURBSSeg(zNURBS *nurbs, double t)
{
  register int i, j, k;

  if( t <  zNURBSKnot0(nurbs) ) return -1;
  if( t >  zNURBSKnotE(nurbs) ) return -2;
  for( i=0, j=zVecSizeNC(nurbs->knot)-1; ; ){
    while( zNURBSKnot(nurbs,i+1) == zNURBSKnot(nurbs,i) && i < j ) i++;
    while( zNURBSKnot(nurbs,j-1) == zNURBSKnot(nurbs,j) && j > i ) j--;
    if( ( k = ( i + j ) / 2 ) == i ) break;
    if( zNURBSKnot(nurbs,k) <= t )
      i = k;
    else
      j = k;
  }
  return i;
}

/* (static)
 * _zNURBSBasis
 * - basis function.
 */
double _zNURBSBasis(zNURBS *nurbs, double t, int i, int r, int k)
{
  double t1, tr1, b=0;

  if( r == 1 )
    return i == k ? 1 : 0;
  if( i > k - r + 1 ){
    t1  = zNURBSKnot(nurbs,i);
    tr1 = zNURBSKnot(nurbs,i+r-1);
    if( tr1 != t1 )
      b += ( t - t1 ) / ( tr1 - t1 ) * _zNURBSBasis(nurbs,t,i,r-1,k);
  }
  if( i <= k ){
    t1  = zNURBSKnot(nurbs,i+1);
    tr1 = zNURBSKnot(nurbs,i+r);
    if( tr1 != t1 )
      b += ( tr1 - t ) / ( tr1 - t1 ) * _zNURBSBasis(nurbs,t,i+1,r-1,k);
  }
  return b;
}

/* zNURBSVec
 * - compute an interpolated vector by a NURBS interpolator.
 */
zVec zNURBSVec(zNURBS *nurbs, double t, zVec v)
{
  register int s, e, i;
  double b, den;

  zVecClear( v );
  den = 0;
  s = _zNURBSSeg( nurbs, t );
  if( s == -1 )
    return zVecCopy( zNURBSCP(nurbs,0), v );
  if( s == -2 )
    return zVecCopy( zNURBSCP(nurbs,zArrayNum(&nurbs->cparray)-1), v );
  e = zMin( s+1, zArrayNum(&nurbs->cparray) );
  for( i=zMax(s-nurbs->dim,0); i<e; i++ ){
    b = zNURBSWeight(nurbs,i) * _zNURBSBasis(nurbs,t,i,nurbs->dim+1,s);
    den += b;
    zVecCatNCDRC( v, b, zNURBSCP(nurbs,i) );
  }
  return zIsTiny(den) ?
    zVecCopy( zNURBSCP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* (static)
 * _zNURBSBasisDiff
 * - differential of a basis function.
 */
double _zNURBSBasisDiff(zNURBS *nurbs, double t, int i, int r, int k, int diff)
{
  double t1, dt, b=0;

  if( diff == 0 )
    return _zNURBSBasis( nurbs, t, i, r, k );
  if( diff > nurbs->dim + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVODR );
    return NAN;
  }
  if( i > k - r + 1 ){
    t1 = zNURBSKnot(nurbs,i);
    if( !zIsTiny( dt = zNURBSKnot(nurbs,i+r-1) - t1 ) )
      b += ( r - 1 ) / dt * _zNURBSBasisDiff(nurbs,t,i,r-1,k,diff-1);
  }
  if( i <= k ){
    t1  = zNURBSKnot(nurbs,i+1);
    if( !zIsTiny( dt = zNURBSKnot(nurbs,i+r) - t1 ) )
      b -= ( r - 1 ) / dt * _zNURBSBasisDiff(nurbs,t,i+1,r-1,k,diff-1);
  }
  return b;
}

/* (static)
 * _zNURBSDenDiff
 * - differential of the denominator of NURBS function.
 */
double _zNURBSDenDiff(zNURBS *nurbs, double t, int s, int e, int diff)
{
  register int i;
  double den;

  for( den=0, i=zMax(s-nurbs->dim,0); i<e; i++ )
    den += zNURBSWeight(nurbs,i) * _zNURBSBasisDiff(nurbs,t,i,nurbs->dim+1,s,diff);
  return den;
}

/* zNURBSVecDiff
 * - compute a differential of an interpolated vector by a NURBS interpolator.
 */
zVec zNURBSVecDiff(zNURBS *nurbs, double t, zVec v, int diff)
{
  register int s, e, i;
  double den, b, t_tmp;
  zVec tmp;

  if( diff == 0 )
    return zNURBSVec( nurbs, t, v );
  if( diff > nurbs->dim + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVODR );
    return NULL;
  }
  if( ( tmp = zVecAlloc( zVecSize(zNURBSCP(nurbs,0)) ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  zVecClear( v );
  den = 0;
  t_tmp = t;
  s = _zNURBSSeg( nurbs, t_tmp );
  if( s == -1 )
    t_tmp = zNURBSKnot0(nurbs);
  else if( s == -2 )
    t_tmp = zNURBSKnotE(nurbs);
  e = zMin( s+1, zArrayNum(&nurbs->cparray) );
  for( i=zMax(s-nurbs->dim,0); i<e; i++ ){
    b = zNURBSWeight(nurbs,i) * _zNURBSBasisDiff(nurbs,t_tmp,i,nurbs->dim+1,s,diff);
    zVecCatNCDRC( v, b, zNURBSCP(nurbs,i) );
  }
  for( i=1; i<diff+1; i++ )
    zVecCatNCDRC( v, -zCombi(diff,i)*_zNURBSDenDiff(nurbs,t,s,e,i), zNURBSVecDiff(nurbs,t_tmp,tmp,diff-i) );
  zVecFree( tmp );
  den = _zNURBSDenDiff( nurbs, t_tmp, s, e, 0 );
  return zIsTiny(den) ?
    zVecCopy( zNURBSCP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* for debug */

/* zNURBSCPArrayFWrite
 * - output an array of control points.
 */
void zNURBSCPArrayFWrite(FILE *fp, zNURBS *nurbs)
{
  register int i;

  for( i=0; i<zArrayNum(&nurbs->cparray); i++ ){
    fprintf( fp, "[%03d] (%g) ", i, zNURBSWeight(nurbs,i) );
    zVecFWrite( fp, zNURBSCP(nurbs,i) );
  }
}
