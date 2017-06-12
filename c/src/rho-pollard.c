#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

/* Global variables */
static bool verbose = false, quiet = false;
static unsigned long g = 0, h = 0, p = 0;

/* Options list */
static struct option long_opts[] = {
  {"verbose", no_argument, NULL, 'v'},
  {"quiet", no_argument, NULL, 'q'},
  {"help", no_argument, NULL, 'H'},
  {"base", required_argument, NULL, 'g'},
  {"modulo", required_argument, NULL, 'p'},
  {"target", required_argument, NULL, 'h'},
  {NULL, 0, NULL, 0}
};

/* Displays the usage of the software ans quits. */
static void usage(int status) {
  if (status == EXIT_SUCCESS)
    fprintf(stdout, "Usage: rho-pollard [OPTIONS]\n"
    "This is an implementation of Pollard's rho algorithm for finding logarithms over Z/pZ\n\n"
		"-p:, --modulo\t [REQUIRED] a prime number p\n"
		"-g:, --base\t [REQUIRED] a primitive element g modulo p\n"
    "-h:, --target\t [REQUIRED] a power of the element g modulo p\n"
    "-v, --verbose\t verbose output\n"
    "-q, --quiet\t quiet output\n"
    "-H, --help\t display this help\n\n"
		"Before using, YOU have to check that:\n"
		"\t-p is a prime number,\n"
		"\t-g is primitive modulo p,\n"
		"\t-h is non zero.\n\n"
    "Try for instance:\n\t./rho-pollard -p999959 -g7 -h3\n"
    "\t./rho-pollard -p99989 -g2 -h107\n"
  );
  else
    fprintf(stderr, "Try 'rho-pollard --help' for more information.\n");

  exit(status);
}

/* Computes a beast's step in the walk */
static void walk(unsigned long* beast, unsigned long* power_g, unsigned long* power_h) {
	if (*beast < (p / 3)) {
		*beast = (*beast * g) % p;
		*power_g = (*power_g + 1) % (p - 1);
	} else if (*beast < (2 * p / 3)) {
		*beast = (*beast * *beast) % p;
		*power_g = (2 * *power_g) % (p - 1);
		*power_h = (2 * *power_h) % (p - 1);
	} else {
		*beast = (*beast * h) % p;
		*power_h = (*power_h + 1) % (p - 1);
	}
}

/* Computes Euclid's algorithm and returns Bézout's identity */
static unsigned long euclid(unsigned long a, unsigned long b, long* u, long* v) {
	unsigned long r0 = a, r1 = b;
	long u0 = 1, v0 = 0, u1 = 0, v1 = 1;
	unsigned long q, swapr;
	long swapu, swapv;

	while (r1 != 0) {
		q = r0 / r1;
		swapr = r1;
		swapu = u1;
		swapv = v1;
		r1 = r0 - q * r1;
		u1 = u0 - q * u1;
		v1 = v0 - q * v1;
		r0 = swapr;
		u0 = swapu;
		v0 = swapv;
	}

	if (u != NULL && v != NULL) {
		*u = u0;
		*v = v0;
	}

	return r0;
}

static unsigned long modPow(unsigned long base, unsigned long exponent, unsigned long modulo) {
  if (exponent == 0) { return 1; }
  else if (exponent == 1) { return base % modulo; }
  else if (exponent % 2 == 0) { return modPow((base * base) % modulo, exponent / 2, modulo) % modulo; }
  else { return (base * modPow((base * base) % modulo, (exponent - 1) / 2, modulo)) % modulo; }
}

static unsigned long getLog(unsigned long power_g, unsigned long power_h, unsigned long modulo) {

	unsigned long r = 0, l = 0;
	long u = 0, v = 0;

	r = euclid(modulo, power_h, &u, &v);

	if (power_g % r != 0) {
        if (!quiet) {
          fprintf(stdout, "Bézout's indentity:\n\tgcd(%lu, %lu) = %lu = (%li)*%lu + (%li)*%lu\n\n"
                        "Solving diophantine equation:\n\tgcd(%lu, %lu) = %lu does not divide %lu\n\n"
												"FAILING TO FIND THE LOGARITHM: returns 0!!!\n\n",
                        power_h, modulo, r, v, power_h, u, modulo,
                        power_h, modulo, r, power_g);
      }
      return 0;
		}

  l = v > 0 ? ((power_g / r) * v) % (modulo / r) : ((power_g / r) * (v + (modulo / r))) % (modulo / r);

	if (!quiet) {
		fprintf(stdout, "Bézout's indentity:\n\tgcd(%lu, %lu) = %lu = (%li)*%lu + (%li)*%lu\n\n",
										power_h, modulo, r, v, power_h, u, modulo);
		fprintf(stdout, "Solving diophantine equation:\n\ta*%lu + n*%lu = %lu\n\tgcd(%lu, %lu) = %lu divides %lu, so there exist solutions\n",
										power_h, modulo, power_g, power_h, modulo, r, power_g);
    if (verbose) {
      fprintf(stdout, "\tSolutions (solving only for a):\n");
      for (unsigned int i = 0; i < r; i++) {
        fprintf(stdout, "\t\ta = %lu (mod %lu)\n", l + (i * modulo / r), modulo);
      }
      fprintf(stdout, "\tThe solution is of the form: %lu = (%lu**%lu)(w**k) (mod %lu)\n"
                    "\twhere w is a %lu-th root of unity, i.e. w = %lu**%lu, and 0 <= k <= %lu\n",
                     h, g, l, p,
                     r, g, modulo / r, r - 1);
    }
    fprintf(stdout, "\tHence, the solution is of the form:\n\t\t%lu = %lu**(%lu + k*%lu) (mod %lu), "
                    "where 0 <= k <= %lu\n",
                     h, g, l, modulo / r, p, r - 1);
    if (!verbose) {
      fprintf(stdout, "\tTrying out each k for k between 0 and %lu: k = ", r - 1);
    }
	}

  if (r == 1) {
    return l;
  } else {
    unsigned long first = modPow(g, l, p), root = modPow(g, modulo / r, p);
    if (verbose) {
      for (unsigned long i = 0; i < r; i++) {
        fprintf(stdout, "\t\tk = %lu:\t%lu**%lu = %lu (mod %lu)\n", i, g, l + (i * modulo / r), modPow(g, l + (i * modulo / r), p), p);
        if ((first * modPow(root, i, p)) % p == h) {
          l = l + (i * modulo / r);
        }
      }
      fprintf(stdout, "\n");
      return l;
    } else {
      unsigned long i = 0;
      for (i = 0; i < r; i++) {
          if ((first * modPow(root, i, p)) % p == h) {
            l = l + (i * modulo / r);
            break;
          }
        }
      if (!quiet) {
        fprintf(stdout, "%lu\n\n", i);
      }
      return l;
    }
  }
}

/* Main function */
int main(int argc, char* argv[]) {

	int optc = 0;

	/* Options and Arguments Parser */
  while ((optc = getopt_long(argc, argv, "p:g:h:vqH", long_opts, NULL)) != -1) {
    switch (optc) {
      case 'H': /* Display the help and exit */
        usage(EXIT_SUCCESS);
        break;

      case 'v': /* Set verbode mode on */
        verbose = true;
        break;

      case 'q': /* Set quiet mode on */
        quiet = true;
        break;

      case 'p': /* Set the modulo */
        p = strtoul(optarg, NULL, 10);
        break;

      case 'g': /* Set the primitve element */
        g = strtoul(optarg, NULL, 10);
        break;

      case 'h': /* Set the target */
        h = strtoul(optarg, NULL, 10);
        break;

      default:
        usage(EXIT_FAILURE);
    }
  }

	if (!p || !g || !h) {
		fprintf(stderr, "error: Missing arguments.\n");
		usage(EXIT_FAILURE);
	}
	if (verbose && quiet) {
		fprintf(stderr, "error: Can't choose both verbose and quiet!\n");
		usage(EXIT_FAILURE);
	}
	if (!quiet) {
		fprintf(stdout, "POLLARD'S RHO ALGORITHM\n\n"
				"Target: \th = %lu\nBase: \t\tg = %lu\nModulo: \tp = %lu\n\n", h, g, p);
	}

	unsigned long tortoise = 1, hare = 1;
	unsigned long power_g_tortoise = 0, power_h_tortoise = 0;
	unsigned long power_g_hare = 0, power_h_hare = 0;

	/* At first, the tortoise makes one step and the hare makes two steps.
	   Then, we carry on this pattern until a collision is found, both
	   the turtoise and the hare stand on the same element in Z/pZ.
	   This is Floyd's cycle finding algorithm. */

	unsigned int step = 0;

	if (verbose) {
		do {
			/* The tortoise makes one step */
			walk(&tortoise, &power_g_tortoise, &power_h_tortoise);

			/* The hare makes two steps */
			walk(&hare, &power_g_hare, &power_h_hare);
			walk(&hare, &power_g_hare, &power_h_hare);

			fprintf(stdout, "Step %u:\n"
					"\tTortoise: \t%lu = (%lu**%lu)(%lu**%lu) (mod %lu)\n"
					"\tHare: \t\t%lu = (%lu**%lu)(%lu**%lu) (mod %lu)\n\n",
					 step, tortoise, g, power_g_tortoise, h, power_h_tortoise, p,
				 	 hare, g, power_g_hare, h, power_h_hare, p);

			step++;
		} while (tortoise != hare);
	} else {
		do {
			/* The tortoise makes one step */
			walk(&tortoise, &power_g_tortoise, &power_h_tortoise);

			/* The hare makes two steps */
			walk(&hare, &power_g_hare, &power_h_hare);
			walk(&hare, &power_g_hare, &power_h_hare);

			step++;
		} while (tortoise != hare);
	}

	step--;

	if (!verbose && !quiet) {
		fprintf(stdout, "Last step (%u):\n"
			"\tTortoise: \t%lu = (%lu**%lu)(%lu**%lu) (mod %lu)\n"
			"\tHare: \t\t%lu = (%lu**%lu)(%lu**%lu) (mod %lu)\n\n",
			step, tortoise, g, power_g_tortoise, h, power_h_tortoise, p,
			hare, g, power_g_hare, h, power_h_hare, p);
	}

	unsigned long power_g = power_g_tortoise > power_g_hare ? power_g_tortoise - power_g_hare : (p - 1) + (power_g_tortoise - power_g_hare);
	unsigned long power_h = power_h_hare > power_h_tortoise ? power_h_hare - power_h_tortoise : (p - 1) + (power_h_hare - power_h_tortoise);

	if (!quiet) {
		fprintf(stdout, "COLLISION FOUND\n\nEquation solving: finding a where %lu = %lu**a (mod %lu)\n"
											"\t(%lu**%lu)(%lu**%lu) = (%lu**%lu)(%lu**%lu) (mod %lu)\n",
											h, g, p, g, power_g_tortoise, h, power_h_tortoise, g, power_g_hare, h, power_h_hare, p);
		if (verbose) {
			fprintf(stdout, "<=>\t%lu**(%lu-%lu) = %lu**(%lu-%lu) (mod %lu)\n"
										"<=>\t%lu**(%lu) = %lu**(%lu) (mod %lu)\n"
										"<=>\t%lu**(%lu) = %lu**(a*%lu) (mod %lu)\n",
										g, power_g_tortoise, power_g_hare, h, power_h_hare, power_h_tortoise, p,
										g, power_g, h, power_h, p,
										g, power_g, g, power_h, p);
		}
		fprintf(stdout, "<=>\t%lu = a*%lu (mod %lu)\n\n", power_g, power_h, p - 1);
	}

	unsigned long loga = getLog(power_g, power_h, p - 1);

	if (quiet) {
		fprintf(stdout, "%lu\n", loga);
	} else {
		fprintf(stdout, "Logarithm: %lu\n\t%lu = %lu**%lu (mod %lu)\n", loga, h, g, loga, p);
	}
	return EXIT_SUCCESS;
}
