#include <zm/zm_misc.h>

int main(int argc, char *argv[])
{
  double s, c;
  double angle;

  if( argc > 1 )
    angle = atof( argv[1] );
  else{
    printf( "angle[deg]> " );
    if( scanf( "%lf", &angle ) == 0 );
  }
  angle = zDeg2Rad( angle );
  zSinCos( angle, &s, &c );
  printf( "sin(%f[deg]) = %f\n", angle, s );
  printf( "cos(%f[deg]) = %f\n", angle, c );
  return 0;
}
