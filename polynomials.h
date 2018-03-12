#include <bits/stdc++.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>

using namespace std;
using namespace NTL;

template<class T>
T eval(const Vec<T> &v, T x) {
	T ans;
	for(int i = 0; i < v.length();i++)
		ans += v[i] * power(x, i);
	return ans;
}
