#include <zm/zm_misc.h>

void printf_ieee_fp(__ieee_fp_t val)
{
  int i;

  printf( "%f (0x", val.d );
  for( i=7; i>=0; i-- )
    printf( "%02x", val.c[i] );
  printf( ")\n" );
}

int main(void)
{
  __ieee_fp_t val;

  /* HUGE_VAL */
  printf( "<HUGE_VAL test>\n" );
  val.d = 100;
  printf( " 100=" ); printf_ieee_fp( val );
  printf( " (inf?=%s)\n", zBoolExpr(isinf(val.d)) );
  val.d = HUGE_VAL + 1;
  printf( " HUGE_VAL=" ); printf_ieee_fp( val );
  printf( " (inf?=%s)\n", zBoolExpr(isinf(val.d)) );
  val.d = -val.d;
  printf( "-HUGE_VAL=" ); printf_ieee_fp( val );
  printf( " (inf?=%s)\n", zBoolExpr(isinf(val.d)) );
  val.d = sqrt(-1);
  printf( " sqrt(-1)=" ); printf_ieee_fp( val );
  printf( " (inf?=%s)\n", zBoolExpr(isinf(val.d)) );

  printf( " HUGE_VAL(def)=%f\n", HUGE_VAL );
  printf( "-HUGE_VAL(def)=%f\n", -HUGE_VAL );
  printf( "%f > %f -> %s\n", HUGE_VAL, 0.0, zBoolExpr( HUGE_VAL > 0 ) );
  printf( "%f > %f -> %s\n",-HUGE_VAL, 0.0, zBoolExpr(-HUGE_VAL > 0 ) );

  /* NAN */
  printf( "<NAN test>\n" );
  val.d = HUGE_VAL;
  printf( " HUGE_VAL=" ); printf_ieee_fp( val );
  printf( " (nan?=%s)\n", zBoolExpr(isnan(val.d)) );
  val.d = sqrt(100);
  printf( " sqrt(100)=" ); printf_ieee_fp( val );
  printf( " (nan?=%s)\n", zBoolExpr(isnan(val.d)) );
  val.d = sqrt(-100);
  printf( " sqrt(-100)=" ); printf_ieee_fp( val );
  printf( " (nan?=%s)\n", zBoolExpr(isnan(val.d)) );
  val.d = log(-100);
  printf( " log(-100)=" ); printf_ieee_fp( val );
  printf( " (nan?=%s)\n", zBoolExpr(isnan(val.d)) );
  val.d = NAN;
  printf( " NAN=" ); printf_ieee_fp( val );
  printf( " (nan?=%s)\n", zBoolExpr(isnan(val.d)) );
  val.d = -val.d;
  printf( "-NAN=" ); printf_ieee_fp( val );
  printf( " (nan?=%s)\n", zBoolExpr(isnan(val.d)) );

  printf( " NAN(def)=%f\n", NAN );
  printf( "-NAN(def)=%f\n", -NAN );
  printf( "%f > %f -> %s\n", NAN, 0.0, zBoolExpr( NAN > 0 ) );
  printf( "%f > %f -> %s\n",-NAN, 0.0, zBoolExpr(-NAN > 0 ) );

  return 0;
}
