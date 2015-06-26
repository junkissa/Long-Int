#include "cmath"
#include "cstring"
#include "string"
#include "fstream"
#include "iostream"
#include "time.h"
#include <complex>
#include <vector>


using namespace std;
class LI;
class StackInt;
class StackLI;

//void fft(vcd &a, bool invert);

class LI
{
private:
	static const int p = 1000;  // systema chislennia
	static const  int max_n = 1000;  // max dlina v zifrah
	int *digits; // masiv zifer	
	int active_n; // kilkist' zajniatyh iachejok

public:
	LI();
	LI(const LI & old);
	~LI();
	void From_Int(const long long int x);
	void From_String(const string & str);
	void From_String_bin(const string & str);
	void From_String_Double(const string & str);
	void From_p_system_to_10();
	void From_10_to_binary();
	void Print_String();
	void Print();
	void Print_Native_Inverse(); // for stack
	void Print_Native();
	void Shift(const int & m); // for * p^m
	void Shift_left(int m); // for Stack
	void Nul(); 

	LI operator + (const LI & second);
	LI operator - (const LI & second);
	LI operator - ();
	LI operator * (const LI & second);
	const LI operator * (int & second);
	const LI operator / (const int & second);
	const LI operator / (const LI & second);
	const LI operator % (const LI & second);
	const int operator % (const int & second);
	const int binmod(const int & second);
	const LI bindiv  (const int & second);
	LI operator = (const LI & old);
	bool operator > (const LI & secind);
	bool operator >= (const LI & second);
	const bool operator == (const LI & second);
	bool operator != (LI & second);
	const LI powLI (const int & n);
	const LI powLI ( LI & n);
	const LI pow_mod (LI & n, const LI & mod);

	const LI gcd(const LI & second);
	const LI Karatsuba(LI & b);
	const LI T_K(const LI & b);
	int ident_max_n(const LI & a, const LI & b);
	int ident_min_n(const LI & a, const LI & b);
	const int Jacobi (const LI & n);
	void divide(StackLI & st, int s, int t); // for Toom-Kook
	void swap (LI & second);

	
	const void S_S();
	const void R_M();
	const void L_L();
	const LI invCook(const string & str);
	LI Strassen(LI & second);
};


class StackInt
{
private:
	int * stack_int;
	static const int SIZE = 100;
	int head;
public:
	StackInt() : head(0) { stack_int = new int [SIZE]; }
	StackInt(const StackInt & old);
	~StackInt() { head = 0; if (stack_int != NULL) delete []stack_int; stack_int = NULL;}
	void push(const int & elem);
	const int get_head() { return head; }
	int pop();
};

class StackLI
{
private:
	LI * stack;
	static const int SIZE = 100;
	int  head;
public:
	StackLI();
	StackLI(const StackLI & old);
	~StackLI();
	void push(const LI & elem);
	void inverse();
	LI pop();
	void Nul();
	void change(int n, const LI & li); // n z verhu
	LI pick(int n); //const не хочe множити в Т_К (крок 6)
	// posmotret' n-toe sverhu 
	const void Print();
	const int get_head()
	{
		return head;
	}
};

/*typedef complex<double> cd;
typedef vector<cd> vcd;
void fft(vcd &a, bool invert)
{
	int n = a.size();
	if (n == 1)
		return;

	vcd a0, a1;
	a0.resize(n / 2); a1.resize(n / 2);
	for (int i = 0, j = 0; i < n; i += 2, ++j)
	{
		a0[j] = a[i];
		a1[j] = a[i + 1];
	}
	fft(a0, invert);
	fft(a1, invert);

	double ang = 2 * 3.14159265359 / n * (invert ? -1 : 1);
	complex<double> w(1), wn(cos(ang), sin(ang));

	for (int i = 0; i < n / 2; ++i)
	{
		a[i] = a0[i] + (w * a1[i]);
		a[i + n / 2] = a0[i] - (w * a1[i]);
		w = w * wn;
	}
	if (invert)
	{
		for (int i = 0; i < n; ++i)
			a[i] /= 2;
	}
}
*/