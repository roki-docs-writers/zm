/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS curve.
 */

#include <zm/zm_nurbs.h>

static int _zNURBSSeg(zNURBS *nurbs, double t);

static double _zNURBSBasis(zNURBS *nurbs, double t, int i, int r, int k);
static double _zNURBSBasisDiff(zNURBS *nurbs, double t, int i, int r, int k, int diff);
static double _zNURBSDenDiff(zNURBS *nurbs, double t, int s, int e, int diff);

/* zNURBSCreate
 * - create a NURBS curve.
 */
bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int dim)
{
  register int i, j;
  zSeqListCell *cp;
  bool ret = true;

  if( zListNum(seq) <= dim ){
    ZRUNERROR( ZM_ERR_NURBS_INVDIM );
    return false;
  }
  nurbs->dim = dim;
  nurbs->knot = zVecAlloc( zListNum(seq)+dim+1 );

  zArrayAlloc( &nurbs->cparray, zNURBSCPCell, zListNum(seq) );
  if( !nurbs->knot || zNURBSCPNum(nurbs) == 0 ){
    ZALLOCERROR();
    zNURBSDestroy( nurbs );
    return false;
  }
  /* set knots & assign control points & initialize weight uniformly */
  for( j=0; j<dim; j++ )
    zNURBSKnot(nurbs,j) = 0;
  i = 0;
  zListForEachRew( seq, cp ){
    zNURBSWeight(nurbs,i) = 1.0;
    if( !( zNURBSCP(nurbs,i) = zVecClone( cp->data.v ) ) )
      ret = false;
    zNURBSKnot(nurbs,j) = zNURBSKnot(nurbs,j-1) + cp->data.dt;
    j++;
    i++;
  }
  for( ; j<zVecSizeNC(nurbs->knot); j++ )
    zNURBSKnot(nurbs,j) = zNURBSKnot(nurbs,j-1);
  if( !ret )
    zNURBSDestroy( nurbs );
  return ret;
}

/* zNURBSDestroy
 * - destroy a NURBS curve.
 */
void zNURBSDestroy(zNURBS *nurbs)
{
  register int i;

  nurbs->dim = 0;
  zVecFree( nurbs->knot );
  nurbs->knot = NULL;
  for( i=0; i<zNURBSCPNum(nurbs); i++ )
    zVecFree( zNURBSCP(nurbs,i) );
  zArrayFree( &nurbs->cparray );
}

/* zNURBSKnotNormalize
 * - normalize the knot vector of a NURBS curve.
 */
void zNURBSKnotNormalize(zNURBS *nurbs)
{
  zVecShift( nurbs->knot, -zNURBSKnot0(nurbs) );
  zVecDivDRC( nurbs->knot, zNURBSKnotE(nurbs) );
}

/* (static)
 * _zNURBSSeg
 * - find the knot segment that involves the given parameter.
 */
int _zNURBSSeg(zNURBS *nurbs, double t)
{
  register int i, j, k;

  if( t < zNURBSKnot0(nurbs) ) return -1;
  if( t >= zNURBSKnotE(nurbs) ) return -2;
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
 * - basis function of NURBS.
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
 * - compute a vector on a NURBS curve.
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
    return zVecCopy( zNURBSCP(nurbs,zNURBSCPNum(nurbs)-1), v );
  e = zMin( s+1, zNURBSCPNum(nurbs) );
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
 * - derivative of the basis function of NURBS.
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
 * - derivative of the denominator of NURBS.
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
 * - compute the derivative a NURBS curve.
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
  e = zMin( s+1, zNURBSCPNum(nurbs) );
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

/* nearest neighbor on NURBS */
#define ZNURBS_NN_DIV 21
double zNURBSVecNN(zNURBS *nurbs, zVec v, zVec nn)
{
  double s1, s2, s1old, s2old, si;
  double d, dmin1, dmin2;
  zVec vs;
  register int i;

  if( !( vs = zVecAlloc( zVecSizeNC(v) ) ) )
    return zNURBSKnot0(nurbs); /* dummy */
  s1 = zNURBSKnot0(nurbs);
  s2 = zNURBSKnotE(nurbs);
  dmin1 = dmin2 = HUGE_VAL;
  do{
    s1old = s1;
    s2old = s2;
    for( i=0; i<ZNURBS_NN_DIV; i++ ){
      si = (s2old-s1old)*i/ZNURBS_NN_DIV + s1old;
      zNURBSVec( nurbs, si, vs );
      if( ( d = zVecDist( v, vs ) ) <= dmin1 ){
        dmin2 = dmin1; s2 = s1;
        dmin1 = d;     s1 = si;
      } else
      if( d <= dmin2 ){
        dmin2 = d;
        s2 = si;
      }
    }
  } while( !zIsTiny( s1 - s2 ) && !zIsTiny( dmin1 - dmin2 ) );
  zNURBSVec( nurbs, ( si = 0.5*(s1+s2) ), nn );
  zVecFree( vs );
  return si;
}

/* for debug */

/* zNURBSCPArrayFWrite
 * - output an array of control points.
 */
void zNURBSCPArrayFWrite(FILE *fp, zNURBS *nurbs)
{
  register int i;

  for( i=0; i<zNURBSCPNum(nurbs); i++ ){
    fprintf( fp, "[%03d] (%g) ", i, zNURBSWeight(nurbs,i) );
    zVecFWrite( fp, zNURBSCP(nurbs,i) );
  }
}
