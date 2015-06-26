#include "Classes.h"





int main()
{
	LI a, b, c, d;
	LI mod, nod, po;
	string a1, b1;
	int Jac;
	cout << " Methods: \n";
	cout << "1 - Karatsuba\n";
	cout << "4 - Strassen\n";
	cout << "7 - S_S\n";
	cout << "8 - R_M\n";
	cout << "9 - L_L\n";

	cout << "Enter a and b than methods 1 and 4 will count the product \n 7 and 8 will check if a is the prime number \n 9 will check if 2^a - 1 is prime (a must be from 3 to 3000) \n try a = 2281=)\n";
	/*LI M;
	M.From_String("2281");
	M.L_L();*/

	
while (cin >> a1 && cin >> b1)
	{
		a.From_String(a1);
		b.From_String(b1);

		c = a.Karatsuba(b); 
					
		c.Print_String();
		d = a.Strassen(b);
		c.Print_String();
		a.S_S();
		a.R_M();
		a.L_L();

	}

	

	cin.get();
	cin.get();
	system("pause");

	return 0;
}

