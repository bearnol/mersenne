/* Factoring Mersennes with trial-division method using divisors 2kp+1, and small-prime sieve improvement (courtesy DT).
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gmp.h"

int
factor_using_trialdiv (mpz_t M, mpz_t S, mpz_t T)
{
  mpz_t ZERO, ONE, TWO, HUNDREDTHOUSAND;
  mpz_t EIGHT,THREE,FIVE,SEVEN,ELEVEN,THIRTEEN,SEVENTEEN,NINETEEN,TWENTYTHREE,TWENTYNINE,THIRTYONE,THIRTYSEVEN,FORTYONE,FORTYTHREE,FORTYSEVEN,FIFTYTHREE,FIFTYNINE;
  mpz_t i, j;
  mpz_t lm8,lm3,lm5,lm7,lm11,lm13,lm17,lm19,lm23,lm29,lm31,lm37,lm41,lm43,lm47,lm53,lm59;
  mpz_t ldm8,ldm3,ldm5,ldm7,ldm11,ldm13,ldm17,ldm19,ldm23,ldm29,ldm31,ldm37,ldm41,ldm43,ldm47,ldm53,ldm59;
  long m8,m3,m5,m7,m11,m13,m17,m19,m23,m29,m31,m37,m41,m43,m47,m53,m59;
  long dm8,dm3,dm5,dm7,dm11,dm13,dm17,dm19,dm23,dm29,dm31,dm37,dm41,dm43,dm47,dm53,dm59;
  mpz_t onep, twop, progress, divisor, result;
  mpz_t mersenne;

	int flag=0;

	mpz_init_set_si (ZERO, 0);
	mpz_init_set_si (ONE, 1);
	mpz_init_set_si (TWO, 2);
	mpz_init_set_si (HUNDREDTHOUSAND, 100000);

	mpz_init_set_si (EIGHT, 8);
	mpz_init_set_si (THREE, 3);
	mpz_init_set_si (FIVE, 5);
	mpz_init_set_si (SEVEN, 7);
	mpz_init_set_si (ELEVEN, 11);
	mpz_init_set_si (THIRTEEN, 13);
	mpz_init_set_si (SEVENTEEN, 17);
	mpz_init_set_si (NINETEEN, 19);
	mpz_init_set_si (TWENTYTHREE, 23);
	mpz_init_set_si (TWENTYNINE, 29);
	mpz_init_set_si (THIRTYONE, 31);
	mpz_init_set_si (THIRTYSEVEN, 37);
	mpz_init_set_si (FORTYONE, 41);
	mpz_init_set_si (FORTYTHREE, 43);
	mpz_init_set_si (FORTYSEVEN, 47);
	mpz_init_set_si (FIFTYTHREE, 53);
	mpz_init_set_si (FIFTYNINE, 59);


	mpz_init_set_si (i, 1);
	mpz_init_set_si (j, 1);
	
	mpz_init_set_si (lm8, -1);
	mpz_init_set_si (lm3, -1);
	mpz_init_set_si (lm5, -1);
	mpz_init_set_si (lm7, -1);
	mpz_init_set_si (lm11, -1);
	mpz_init_set_si (lm13, -1);
	mpz_init_set_si (lm17, -1);
	mpz_init_set_si (lm19, -1);
	mpz_init_set_si (lm23, -1);
	mpz_init_set_si (lm29, -1);
	mpz_init_set_si (lm31, -1);
	mpz_init_set_si (lm37, -1);
	mpz_init_set_si (lm41, -1);
	mpz_init_set_si (lm43, -1);
	mpz_init_set_si (lm47, -1);
	mpz_init_set_si (lm53, -1);
	mpz_init_set_si (lm59, -1);
	
	mpz_init_set_si (ldm8, 1);
	mpz_init_set_si (ldm3, 1);
	mpz_init_set_si (ldm5, 1);
	mpz_init_set_si (ldm7, 1);
	mpz_init_set_si (ldm11, 1);
	mpz_init_set_si (ldm13, 1);
	mpz_init_set_si (ldm17, 1);
	mpz_init_set_si (ldm19, 1);
	mpz_init_set_si (ldm23, 1);
	mpz_init_set_si (ldm29, 1);
	mpz_init_set_si (ldm31, 1);
	mpz_init_set_si (ldm37, 1);
	mpz_init_set_si (ldm41, 1);
	mpz_init_set_si (ldm43, 1);
	mpz_init_set_si (ldm47, 1);
	mpz_init_set_si (ldm53, 1);
	mpz_init_set_si (ldm59, 1);
	

	mpz_init_set_si (onep, 1);
	mpz_init_set_si (twop, 1);
	mpz_init_set_si (progress, 1);
	mpz_init_set_si (divisor, 1);
	mpz_init_set_si (result, 1);

	mpz_init_set_si (mersenne, 1);
	

    mpz_set (mersenne, M);

	mpz_set (onep, mersenne);
	mpz_add (twop, onep, onep);
	
    mpz_set (j, S);
	

	mpz_set (divisor, j);
	mpz_mul (divisor, divisor, twop);
	mpz_add (divisor, divisor, ONE);
	
      // pre-compute the modulus for 2*p and for 2kp+1 for each small prime

	mpz_mod (ldm8, twop, EIGHT);
	mpz_mod (ldm3, twop, THREE);
	mpz_mod (ldm5, twop, FIVE);
	mpz_mod (ldm7, twop, SEVEN);
	mpz_mod (ldm11, twop, ELEVEN);
	mpz_mod (ldm13, twop, THIRTEEN);
	mpz_mod (ldm17, twop, SEVENTEEN);
	mpz_mod (ldm19, twop, NINETEEN);
	mpz_mod (ldm23, twop, TWENTYTHREE);
	mpz_mod (ldm29, twop, TWENTYNINE);
	mpz_mod (ldm31, twop, THIRTYONE);
	mpz_mod (ldm37, twop, THIRTYSEVEN);
	mpz_mod (ldm41, twop, FORTYONE);
	mpz_mod (ldm43, twop, FORTYTHREE);
	mpz_mod (ldm47, twop, FORTYSEVEN);
	mpz_mod (ldm53, twop, FIFTYTHREE);
	mpz_mod (ldm59, twop, FIFTYNINE);

	dm8=mpz_get_ui(ldm8);
	dm3=mpz_get_ui(ldm3);
	dm5=mpz_get_ui(ldm5);
	dm7=mpz_get_ui(ldm7);
	dm11=mpz_get_ui(ldm11);
	dm13=mpz_get_ui(ldm13);
	dm17=mpz_get_ui(ldm17);
	dm19=mpz_get_ui(ldm19);
	dm23=mpz_get_ui(ldm23);
	dm29=mpz_get_ui(ldm29);
	dm31=mpz_get_ui(ldm31);
	dm37=mpz_get_ui(ldm37);
	dm41=mpz_get_ui(ldm41);
	dm43=mpz_get_ui(ldm43);
	dm47=mpz_get_ui(ldm47);
	dm53=mpz_get_ui(ldm53);
	dm59=mpz_get_ui(ldm59);

	mpz_mul (lm8, ldm8, j);
	mpz_add (lm8, lm8, ONE);
	mpz_mod (lm8, lm8, EIGHT);
	m8=mpz_get_ui(lm8);

	mpz_mul (lm3, ldm3, j);
	mpz_add (lm3, lm3, ONE);
	mpz_mod (lm3, lm3, THREE);
	m3=mpz_get_ui(lm3);

	mpz_mul (lm5, ldm5, j);
	mpz_add (lm5, lm5, ONE);
	mpz_mod (lm5, lm5, FIVE);
	m5=mpz_get_ui(lm5);

	mpz_mul (lm7, ldm7, j);
	mpz_add (lm7, lm7, ONE);
	mpz_mod (lm7, lm7, SEVEN);
	m7=mpz_get_ui(lm7);

	mpz_mul (lm11, ldm11, j);
	mpz_add (lm11, lm11, ONE);
	mpz_mod (lm11, lm11, ELEVEN);
	m11=mpz_get_ui(lm11);

	mpz_mul (lm13, ldm13, j);
	mpz_add (lm13, lm13, ONE);
	mpz_mod (lm13, lm13, THIRTEEN);
	m13=mpz_get_ui(lm13);

	mpz_mul (lm17, ldm17, j);
	mpz_add (lm17, lm17, ONE);
	mpz_mod (lm17, lm17, SEVENTEEN);
	m17=mpz_get_ui(lm17);

	mpz_mul (lm19, ldm19, j);
	mpz_add (lm19, lm19, ONE);
	mpz_mod (lm19, lm19, NINETEEN);
	m19=mpz_get_ui(lm19);

	mpz_mul (lm23, ldm23, j);
	mpz_add (lm23, lm23, ONE);
	mpz_mod (lm23, lm23, TWENTYTHREE);
	m23=mpz_get_ui(lm23);

	mpz_mul (lm29, ldm29, j);
	mpz_add (lm29, lm29, ONE);
	mpz_mod (lm29, lm29, TWENTYNINE);
	m29=mpz_get_ui(lm29);

	mpz_mul (lm31, ldm31, j);
	mpz_add (lm31, lm31, ONE);
	mpz_mod (lm31, lm31, THIRTYONE);
	m31=mpz_get_ui(lm31);

	mpz_mul (lm37, ldm37, j);
	mpz_add (lm37, lm37, ONE);
	mpz_mod (lm37, lm37, THIRTYSEVEN);
	m37=mpz_get_ui(lm37);

	mpz_mul (lm41, ldm41, j);
	mpz_add (lm41, lm41, ONE);
	mpz_mod (lm41, lm41, FORTYONE);
	m41=mpz_get_ui(lm41);

	mpz_mul (lm43, ldm43, j);
	mpz_add (lm43, lm43, ONE);
	mpz_mod (lm43, lm43, FORTYTHREE);
	m43=mpz_get_ui(lm43);

	mpz_mul (lm47, ldm47, j);
	mpz_add (lm47, lm47, ONE);
	mpz_mod (lm47, lm47, FORTYSEVEN);
	m47=mpz_get_ui(lm47);
	
	mpz_mul (lm53, ldm53, j);
	mpz_add (lm53, lm53, ONE);
	mpz_mod (lm53, lm53, FIFTYTHREE);
	m53=mpz_get_ui(lm53);

	mpz_mul (lm59, ldm59, j);
	mpz_add (lm59, lm59, ONE);
	mpz_mod (lm59, lm59, FIFTYNINE);
	m59=mpz_get_ui(lm59);


	while (mpz_cmp (j, T) < 0)
	    {
		
		mpz_mod (progress, j, HUNDREDTHOUSAND);
		if (mpz_cmp (progress, ZERO) == 0) {
		    mpz_out_str (stdout, 10, j);
			printf("\r");
			fflush (stdout);
		}

         // check for correct mod 8 and no small prime factor
		if ((m8 == 1 || m8 == 7) &&
			m3 != 0 && m5 != 0 && m7 != 0 && m11 != 0 &&
			m13 != 0 && m17 != 0 && m19 != 0 && m23 != 0 &&
			m29 != 0 && m31 != 0 && m37 != 0 && m41 != 0 &&
			m43 != 0 && m47 != 0 && m53 != 0 && m59 != 0) {
				mpz_powm (result, TWO, onep, divisor);
				mpz_sub (result, result, ONE);
				mpz_mod (result, result, divisor);
				if (mpz_cmp (result, ZERO) == 0)
					{
					mpz_out_str (stdout, 10, j);
					printf("\t");
					mpz_out_str (stdout, 10, divisor);
					printf ("\n");
					fflush (stdout);
					}
		}
		
		
		mpz_add (divisor, divisor, twop);

		mpz_add (j, j, ONE);
		

         // update the moduli reflecting divisor=divisor+2p		 
         m8  += dm8;  if(m8  >= 8)  m8  -= 8;
         m3  += dm3;  if(m3  >= 3)  m3  -= 3;
         m5  += dm5;  if(m5  >= 5)  m5  -= 5;
         m7  += dm7;  if(m7  >= 7)  m7  -= 7;
         m11 += dm11; if(m11 >= 11) m11 -= 11;
         m13 += dm13; if(m13 >= 13) m13 -= 13;
         m17 += dm17; if(m17 >= 17) m17 -= 17;
         m19 += dm19; if(m19 >= 19) m19 -= 19;
         m23 += dm23; if(m23 >= 23) m23 -= 23;
         m29 += dm29; if(m29 >= 29) m29 -= 29;
         m31 += dm31; if(m31 >= 31) m31 -= 31;
         m37 += dm37; if(m37 >= 37) m37 -= 37;
         m41 += dm41; if(m41 >= 41) m41 -= 41;
         m43 += dm43; if(m43 >= 43) m43 -= 43;
         m47 += dm47; if(m47 >= 47) m47 -= 47;
         m53 += dm53; if(m53 >= 53) m53 -= 53;
         m59 += dm59; if(m59 >= 59) m59 -= 59;
		
	}




	
		
  mpz_clear (ZERO);
  mpz_clear (ONE);
  mpz_clear (TWO);
  mpz_clear (HUNDREDTHOUSAND);
  
  mpz_clear (EIGHT);
  mpz_clear (THREE);
  mpz_clear (FIVE);
  mpz_clear (SEVEN);
  mpz_clear (ELEVEN);
  mpz_clear (THIRTEEN);
  mpz_clear (SEVENTEEN);
  mpz_clear (NINETEEN);
  mpz_clear (TWENTYTHREE);
  mpz_clear (TWENTYNINE);
  mpz_clear (THIRTYONE);
  mpz_clear (THIRTYSEVEN);
  mpz_clear (FORTYONE);
  mpz_clear (FORTYTHREE);
  mpz_clear (FORTYSEVEN);
  mpz_clear (FIFTYTHREE);
  mpz_clear (FIFTYNINE);

  mpz_clear (i);
  mpz_clear (j);

  mpz_clear (lm8);
  mpz_clear (lm3);
  mpz_clear (lm5);
  mpz_clear (lm7);
  mpz_clear (lm11);
  mpz_clear (lm13);
  mpz_clear (lm17);
  mpz_clear (lm19);
  mpz_clear (lm23);
  mpz_clear (lm29);
  mpz_clear (lm31);
  mpz_clear (lm37);
  mpz_clear (lm41);
  mpz_clear (lm43);
  mpz_clear (lm47);
  mpz_clear (lm53);
  mpz_clear (lm59);

  mpz_clear (ldm8);
  mpz_clear (ldm3);
  mpz_clear (ldm5);
  mpz_clear (ldm7);
  mpz_clear (ldm11);
  mpz_clear (ldm13);
  mpz_clear (ldm17);
  mpz_clear (ldm19);
  mpz_clear (ldm23);
  mpz_clear (ldm29);
  mpz_clear (ldm31);
  mpz_clear (ldm37);
  mpz_clear (ldm41);
  mpz_clear (ldm43);
  mpz_clear (ldm47);
  mpz_clear (ldm53);
  mpz_clear (ldm59);
  

  mpz_clear (onep);
  mpz_clear (twop);
  mpz_clear (progress);
  mpz_clear (divisor);
  mpz_clear (result);

  mpz_clear (mersenne);

return (flag);

}

main (int argc, char *argv[])
{
    mpz_t M, S, T;
    int i;
	int flag=0;

    mpz_init (M);
    mpz_init (S);
    mpz_init (T);

    if (argc > 1)
        {
        for (i = 1; i < argc; i++)
	        {
			if (!strncmp (argv[i], "-M", 2))
	            {
	            mpz_set_str (M, argv[i] + 2, 0);
	            }
	        else if (!strncmp (argv[i], "-S", 2))
	            {
	            mpz_set_str (S, argv[i] + 2, 0);
	            }
	        else if (!strncmp (argv[i], "-T", 2))
	            {
	            mpz_set_str (T, argv[i] + 2, 0);
	            }
            }
	

	     flag = factor_using_trialdiv (M, S, T);

	     }

	mpz_clear(M);
	mpz_clear(S);
	mpz_clear(T);

    exit (flag);
}
