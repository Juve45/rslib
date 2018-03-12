#include <bits/stdc++.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include "reed-solomon.h"

#define dbg(x) cerr<<#x": "<<x<<endl
#define dbg_p(x) cerr<<#x": "<<x.first<<","<<x.second<<"\n"
#define dbg_v(x, n) do{cerr<<#x"[]: ";for(int _=0;_<n;++_)cerr<<x[_]<<" ";cerr<<'\n';}while(0)
#define dbg_ok cerr<<"OK!\n"

#define st first
#define nd second

#define DMAX 
#define NMAX 
#define MMAX 

using namespace std;
using namespace NTL;

int n, k;
double l, v1, v2;
string s;


int main()
{
	ios_base::sync_with_stdio(false);
	
	reed_solomon rs(512, 2); //p = 11, s = 1
	Vec<ZZ_p> a, b;

	unsigned char msg[1000], *ss = (unsigned char*)"asta e rpimul cuvant ass =asdasdggMMtyb=asdasdggMMtb=asdasdggMMtb=asda sdggMMt ybfsadfwasta e rpimul cuvant ass =asdasdggMMtyb=asdasdggMMtb=asdasdggMMtb=asda sdggMMt ybfsadfwasta e rpimul cuvant ass =asdasdggMMtyb=asdasdggMMtb=asdasdggMMtb=asda sdggMMt ybfsadfwasta e rpimul cuvant ass =asdasdggMMtyb=asdasdggMMtb=asdasdggMMtb=asda sdggMMt ybfsadfwqersadf=asdasdggMasdasdggMMtyb=asdasdggMMtybfsadfwqersadf=asdasdggM A lalala asdasdggM A lalala Mtyb=asdasdggMMtyb=asdas Salut ce mai faci? mai faci? mai faci? mai faci? mai faci? mai faci? mai faci? mai faci? mai faci? dggMMtyb=asdasdggMMtybsadfwqerwioejrowi";
	int lg = strlen((char *)ss), len;
	
	dbg(lg);
	
	b = rs.encode(ss, lg);
	
	//stric putin
	b[2] = 1117;
	b[0] = 117;
	
	rs.decode(b, msg, len);
	
	cout << len;
	cout << "[" << msg << "]\n";
	

	// b[0] = 117;
	// Vec<ZZ_p> y = rs.decode(b);
	// for(int i = 1; i < y.length(); i++) {
	// 	unsigned char buff[3];
	// 	BytesFromZZ(buff, rep(y[i]), 2);
	// 	cout << buff;
	// }
}