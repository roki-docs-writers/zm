/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS interpolation.
 */

#ifndef __ZM_NURBS_H__
#define __ZM_NURBS_H__

#include <zm/zm_seq.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \struct zNURBSCPCell
 * \brief cell of an array of control points for NURBS.
 *
 * zNURBSCPCell is a cell of an array of control points for NURBS.
 * It contains a control point and a weight put on it.
 *//* ******************************************************* */
typedef struct{
  zVec cp;  /*!< control point */
  double w; /*!< weight */
} zNURBSCPCell;

/* ********************************************************** */
/*! \struct zNURBSCPArray
 * \brief an array of control points for NURBS.
 *
 * zNURBSCPArray is an array of control points for NURBS.
 * It is defined with a macro zArrayClass.
 * \sa zArrayClass.
 *//* ******************************************************* */
zArrayClass( zNURBSCPArray, zNURBSCPCell );

/* ********************************************************** */
/*! \struct zNURBS
 * \brief NURBS interpolator.
 *
 * zNURBS is an NURBS interpolator of a sequence.
 *//* ******************************************************* */
typedef struct{
  int dim;   /*!< \brief dimension of a curve */
  zSeq *seq; /*!< \brief a sequence to be interpolated */
  /*! \cond */
  zVec knot;             /* knot vector */
  zNURBSCPArray cparray; /* an array of control points */
  /*! \endcond */
} zNURBS;

#define zNURBSKnot(n,i)    zVecElem((n)->knot,i)
#define zNURBSWeight(n,i)  ( zArrayElem(&(n)->cparray,i)->w )
#define zNURBSCP(n,i)      ( zArrayElem(&(n)->cparray,i)->cp )

#define zNURBSKnot0(n)     zNURBSKnot(n,0)
#define zNURBSKnotE(n)     zNURBSKnot(n,zVecSizeNC((n)->knot)-1)

/*! \brief create NURBS interpolator.
 *
 * zNURBSCreate() creates a NURBS interpolator \a nurbs for a given
 * sequence \a seq.
 * \a dim is the dimension of the interpolation curve, which has to
 * be less than the size of \a seq.
 * It is initialized as a uniform Bezier spline curve with fixed
 * boundary points. The weights on each control point and the knot
 * vector can be modified later.
 * \return the false value if \a dim is larger than the size of
 *         \a seq plus one.
 * \return the false value if it fails to allocate internal workspace.
 * \return the true value if it succeeds to create an interpolator.
 */
__EXPORT bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int dim);

/*! \brief destroy NURBS interpolator.
 *
 * zNURBSDestroy() destroys a NURBS interpolator \a nurbs.
 */
__EXPORT void zNURBSDestroy(zNURBS *nurbs);

/*! \brief normalize the knot vector of a NURBS interpolator.
 *
 * zNURBsKnotNormalize() normalizes the knot vector of a NURBS
 * interpolator \a nurbs so that it starts from 0 and ends at 1
 * by using zNURBSKnot0 and zNURBSKnotE.
 */
__EXPORT void zNURBSKnotNormalize(zNURBS *nurbs);

/*! \brief zNURBSVec
 *
 * zNURBSVec() computes an interpolated vector by a NURBS interpolator
 * \a nurbs at the abscissa \a t. The computed vector will be stored
 * where \a v points.
 * \return \a v if it succeeds to interpolate with a valid \a t.
 * \return the null vector if \a t is invalid.
 */
__EXPORT zVec zNURBSVec(zNURBS *nurbs, double t, zVec v);

/*! \brief zNURBSVecDiff
 *
 * zNURBSVecDiff() computes a differential of an interpolated vector
 * by a NURBS interpolator
 */
__EXPORT zVec zNURBSVecDiff(zNURBS *nurbs, double t, zVec v, int diff);

/*! \brief nearest neighbor on NURBS
 *
 * zNURBSVecNN() finds the nearest-neighbor vector on a NURBS curve
 * defined by \a nurbs from a vector \a v. The result is put into \a nn.
 * \return
 * zNURBSVecNN() returns the parameter corresponding to the nearest-
 * neighbor vector found by this function.
 */
__EXPORT double zNURBSVecNN(zNURBS *nurbs, zVec v, zVec nn);

/* for debug */

__EXPORT void zNURBSCPArrayFWrite(FILE *fp, zNURBS *nurbs);

__END_DECLS

#endif /* __ZM_NURBS_H__ */
