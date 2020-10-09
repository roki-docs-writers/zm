#include <zm/zm_misc.h>

int main(void)
{
  void daa(int n){
    if( zIsEven( n ) )
      printf( "%d is an even number.\n", n );
    if( zIsOdd( n ) )
      printf( "%d is an odd number.\n", n );
  }

  daa(-4 );
  daa(-3 );
  daa(-2 );
  daa(-1 );
  daa( 0 );
  daa( 1 );
  daa( 2 );
  daa( 3 );
  daa( 4 );
  daa( 5 );

  return 0;
}
