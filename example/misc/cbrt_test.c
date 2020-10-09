#include <zm/zm_misc.h>

int main(int argc, char *argv[])
{
  double x;

  if( argc <= 1 ) return 1;
  x = atof( argv[1] );
  printf( "%.15lf %.15lf\n", x, zCbrt( x ) );
  return 0;
}
