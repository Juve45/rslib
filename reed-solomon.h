#include <bits/stdc++.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "polynomials.h"

#define dbg(x) cerr<<#x": "<<x<<endl
#define dbg_p(x) cerr<<#x": "<<x.first<<","<<x.second<<"\n"
#define dbg_v(x, n) do{cerr<<#x"[]: ";for(int _=0;_<n;++_)cerr<<x[_]<<" ";cerr<<'\n';}while(0)
#define dbg_ok cerr<<"OK!\n"

using namespace std;
using namespace NTL;


template<class T>
ostream& operator<<(ostream& out, vector<T> v)
{
	out << v.size() << '\n';
	for(auto e : v) out << e << ' ';
	return out;
}

class reed_solomon {
	ZZ p;
	int p_size;
	int bk_size;
	int s;

public:
	
	reed_solomon(int sz, int _s);
	// constructor, which forces the user to initialise 
	// the size of the prime number which RS is using

	void encode(const Vec<ZZ_p> &a, Vec<ZZ_p> &b);
	// Encodes a vector of numbers modulo p using RS
	// b is the resulting vector

	Vec<ZZ_p> encode(const Vec<ZZ_p> &a);
	// Encodes a vector of numbers modulo p using RS
	// the resulting vector is given via return
	
	Vec<ZZ_p> encode(const unsigned char* str, int len);

	void decode(const Vec<ZZ_p> &a, Vec<ZZ_p> &b);
	Vec<ZZ_p> decode(const Vec<ZZ_p> &a);
	void decode(Vec<ZZ_p> &y, unsigned char* m, int &len);
	
private:
	Vec<ZZ_p> encode(const unsigned char* str, int len, int block_size);

	ZZ_p compute_fc(const Vec<ZZ_p> &z, vector<int> &a);
	ZZ_p compute_fc_naive(const Vec<ZZ_p> &z, vector<int> &a);
	ZZ_p compute_fc_0_inv(const Vec<ZZ_p> &z, vector<int> &a);

	void compute_polynomial_naive(const Vec<ZZ_p> &z, vector <int> &a, Vec<ZZ_p> &y);
	void compute_polynomial_1_inv(const Vec<ZZ_p> &z, vector <int> &a, Vec<ZZ_p> &y);
};

reed_solomon::reed_solomon(int sz, int _s) {
	p_size = sz;
	bk_size = sz/8 - 1;
	p = GenPrime_ZZ(p_size);	
	dbg(p);
	ZZ_p::init(p);
	s = _s;
}

void reed_solomon::encode(const Vec<ZZ_p> &a, Vec<ZZ_p> &b) {
	b.SetLength(a.length() + 2 * s);
	for(int i = 1; i <= a.length() + 2 * s; i++)
		b[i - 1] = eval(a, ZZ_p(i));
}

Vec<ZZ_p> reed_solomon::encode(const Vec<ZZ_p> &a) {
	Vec<ZZ_p> b;
	encode(a, b);
	return b;
}


ZZ_p reed_solomon::compute_fc(const Vec<ZZ_p> &z, vector<int> &a) {
	return compute_fc_0_inv(z, a);

}

ZZ_p reed_solomon::compute_fc_naive(const Vec<ZZ_p> &z, vector<int> &a) {
	
	ZZ_p fc = ZZ_p(0);
	for(int i = 0; i < a.size(); i++) {
		ZZ_p tmp = z[a[i] - 1];
		for(int j = 0; j < a.size(); j++) {
			if(i == j) continue;
			tmp = tmp * ZZ_p(a[j]) * inv(ZZ_p(a[j] - a[i]));
		}
		fc += tmp;
	}	
	// dbg_ok;
	return fc;
}


ZZ_p reed_solomon::compute_fc_0_inv(const Vec<ZZ_p> &z, vector<int> &a) {
	
	ZZ_p fc = ZZ_p(0);
	ZZ_p in[a.size()];
	ZZ_p pre[a.size()];
	ZZ_p suf[a.size()];

	for(int i = 0; i < a.size(); i++) {
		ZZ_p I = ZZ_p(1);

		for(int j = 0; j < a.size(); j++) {
			if(i == j) continue;
			I *= ZZ_p(a[j] - a[i]);
		}
		in[i] = I;
	}	

	pre[0] = ZZ_p(1);
	for(int i = 1; i < a.size(); i++)
		pre[i] = pre[i - 1] * in[i - 1];

	suf[a.size() - 1] = ZZ_p(1);
	for(int i = (int) a.size() - 1; i >=0; i--) {
		suf[i] = suf[i + 1] * in[i + 1];
	}

	for(int i = 0; i < a.size(); i++) {
		ZZ_p tmp = z[a[i] - 1], I = ZZ_p(1);

		for(int j = 0; j < a.size(); j++) {
			if(i == j) continue;
			tmp = tmp * ZZ_p(a[j]);
			I *= ZZ_p(a[j] - a[i]);
		}
		in[i] = I;
		fc += tmp * pre[i] * suf[i];
	}	
	// dbg_ok;
	return fc;
}


Vec<ZZ_p> mul_pol(Vec<ZZ_p> &a, ZZ_p b) {
	Vec<ZZ_p> ans;
	ans.SetLength(a.length() + 1);
	for(int i = 1; i <= a.length(); i++)
		ans[i] = a[i - 1];
	for(int i = 0; i < a.length(); i++)
		ans[i] += a[i] * b;
	return ans;
}

void reed_solomon::compute_polynomial_naive(const Vec<ZZ_p> &z, vector <int> &a, Vec<ZZ_p> &y) {

	y.SetLength(a.size());
	Vec<ZZ_p> t;

	for(auto i : a) {
		
		t.SetLength(1);	t[0] = 1;
		ZZ_p total_inv(1);

		for(auto j : a) {
			if(i == j) continue;
			t = mul_pol(t, ZZ_p(-j));
			total_inv *=  ZZ_p(i - j);
		}

		for(int ii = 0; ii < a.size(); ii++)
			y[ii] += t[ii] * z[i - 1] * inv(total_inv);
	}

}


Vec<ZZ_p> reed_solomon::encode(const unsigned char* str, int len, int block_size){

	int k = (len + block_size - 1) / block_size; //number of blocks
	Vec<ZZ_p> a;
	a.SetLength(k + 1);
	dbg(block_size);
	dbg(k);
	a[0] = ZZ_p(0);
	for(int i = 0; i + 1 < k; i++) {
		// dbg(i);
		a[i + 1] = to_ZZ_p(ZZFromBytes((const unsigned char*)(str + i * block_size), block_size));
	}
	a[k] = to_ZZ_p(ZZFromBytes((const unsigned char*)(str + (k - 1) * block_size), len - (k - 1) * block_size));
	dbg(a);
	return encode(a);
}



Vec<ZZ_p> reed_solomon::encode(const unsigned char* str, int len) {
	return encode(str, len, bk_size);
}

void reed_solomon::decode(const Vec<ZZ_p> &z, Vec<ZZ_p> &y) {

	srand(time(NULL));
	int max_nr = 100;
	vector <int> p, a;
	for(int i = 1; i <= z.length(); i++)
		p.push_back(i);
	for(int i = 0; i < max_nr; i++) {
		dbg(i);
		a.clear();
		int nr = p.size();
		for(int j = 0; j < z.length() - 2 * s; j++) {
			int ind = rand() % nr;
			nr--;
			swap(p[ind], p[nr]);
			a.push_back(p[nr]);
		}
		ZZ_p fc = compute_fc_naive(z, a);
		// cout << a << fc << '\n';
		if(fc == 0) {
			compute_polynomial_naive(z, a, y);
			// dbg(y);
			return ;
		}
	}
	// vector <int> v({1, 3, 4});
	return;
}

Vec<ZZ_p> reed_solomon::decode(const Vec<ZZ_p> &a) {
	Vec<ZZ_p> b;
	decode(a, b);
	return b;
}

void reed_solomon::decode(Vec<ZZ_p> &b, unsigned char* m, int &len) {

	Vec<ZZ_p> y = decode(b);
	len = bk_size * (y.length() - 1);
	for(int i = 1; i < y.length(); i++) {
		unsigned char buff[bk_size];
		BytesFromZZ(buff, rep(y[i]), bk_size);
		memcpy(m + bk_size * (i - 1), buff, bk_size);
	}
}
