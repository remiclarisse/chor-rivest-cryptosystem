#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool verbose = true;

int main () {

	/* This is an implementation of Pollard's rho algorithme for
	   finding logarithms over Z/pZ. */
	/* Before using, YOU have to check that :
		p is a prime number,
		g is primitve modulo p,
		h is non zero. */

	long unsigned g;
	long unsigned h;
	long unsigned p;

	p = 101;
	g = 7;
	h = 4;

	if (verbose) {
		fprintf (stdout, "\nMÉTHODE RHO DE POLLARD\n\n"
			"Calcul du logarithme de h = %lu dans la base g = %lu modulo p = %lu\n\n", h, g, p);
	}

	long unsigned tortoise = 1;
	long unsigned hare = 1;
	unsigned int power_g_tortoise = 0;
	unsigned int power_h_tortoise = 0;
	unsigned int power_g_hare = 0;
	unsigned int power_h_hare = 0;

	/* At first, the tortoise makes one step and the hare makes two steps.
	   Then, we carry on this pattern until a collision is found, both
	   the turtoise and the hare stand on the same element in Z/pZ.
	   This is Floyd's cycle finding algorithm. */

	int step = 0;

	do {
		/* The tortoise makes one step */
		if (tortoise < (p / 3)) {
			tortoise = (tortoise * h) % p;
			power_h_tortoise = (power_h_tortoise + 1) % (p - 1);
		} else if (tortoise < (2 * p / 3)) {
			tortoise = (tortoise * tortoise) % p;
			power_g_tortoise = (2 * power_g_tortoise) % (p - 1);
			power_h_tortoise = (2 * power_h_tortoise) % (p - 1);
		} else {
			tortoise = (tortoise * g) % p;
			power_h_tortoise = (power_h_tortoise + 1) % (p - 1);
		}

		if (verbose) {
			fprintf (stdout, "Step %d :\n"
				"Tortoise : %lu = (%lu**%u)(%lu**%u)\t(mod %lu)\n",
				 step, tortoise, g, power_g_tortoise, h, power_h_tortoise, p);

		}

		hare = tortoise;
	} while (tortoise != hare);

	fprintf (stdout, "\n");
	return EXIT_SUCCESS;
}
