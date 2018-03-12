#include "reed-solomon.h"

static void reed_solomon::encode(ZZ* a, int k, int s, ZZ* &y, int &n) {

	n = k + 2 * s;
	y = new ZZ[n];
	for(int i = 0; i < k; i++)
		y[i] = a[i];

}
