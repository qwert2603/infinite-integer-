#include <iostream>
#include "infinite_integer.h"
using namespace std;


int main() {
	typedef infinite_integer::Digit D; 
	D d1 = "-349657539278576539476405850968365684535645346785646453";
	D d2 = 91345747845601;
	cout << d1 + d2 << '\n' << d1 - d2 << endl;
	cout << d1 / d2 << '\n' << d1 % d2 << endl;
	cout << d2.to_long_long() << endl;
	D d3 = d1 * d2;
	cout << d3 << endl;
	d3 /= string("10000000000000000000000000000000000000000000000");
	cout << d3 << endl;
	cout << (d1 == d2);
	
	return 0;
}
