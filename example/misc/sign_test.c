#include <zm/zm_misc.h>

void test(double a, double b)
{
  printf( "%f %f - %s\n", a, b, zBoolExpr( zIsSgnOpp(a,b) ) );
}

int main(void)
{
  test( 1.0, 1.0 );
  test( 1.0,-1.0 );
  test(-1.0, 1.0 );
  test(-1.0,-1.0 );
  test( 0.0, 0.0 );
  test( 0.0,-0.0 );
  test(-0.0, 0.0 );
  test(-0.0,-0.0 );
  return 0;
}
