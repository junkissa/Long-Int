#include "Classes.h"
const int MIN_KARATSUBA = 4;



StackInt :: StackInt (const StackInt & old) : head(old.head)
{
	stack_int = new int[SIZE];
	for (int i = 0; i < head; ++i) 
		stack_int[i] = old.stack_int[i];
}

void StackInt::push(const int & elem)
{
	if (head < SIZE)
		stack_int[head++] = elem;
	else {
		cout << "Stack is full\n"; exit(1);
	}
	return;
}


int StackInt::pop()
{
	if (head <= 0)
	{
		cout << "Stack is empty\n"; exit(1);
	}
	return stack_int[--head];
}

StackLI::StackLI()
{
	stack = new LI[SIZE];
	head = 0;
}

StackLI::StackLI(const StackLI & old)
{
	stack = new LI [SIZE];
	for (int i = 0; i <= old.head; ++i)
		stack[i] = old.stack[i];
	head = old.head;
}

StackLI :: ~StackLI()
{
	if (stack != NULL)
		delete[]stack;
	head = 0;
}

const void  StackLI::Print()
{
	for (int i = 1; i < head + 1; ++i)
	{
		// (stack[SIZE - i]).Print();
		cout << "i = " << i << "  ";
		(stack[head - i]).Print_Native_Inverse();
	}

}

void StackLI::push(const LI & elem)
{
	if (head < SIZE)
		stack[head++] = elem;
	else {
		cout << "Stack is full\n"; exit(1);
	}
	return;
}

LI StackLI::pop()
{
	if (head <= 0)
	{
		cout << "Stack is empty\n"; exit(1);
	}
	return stack[--head];
}

LI StackLI::pick(int n)
{
	return stack[head - 1 - n];
}

void StackLI::change(int n, const LI & li)/// Wrong work!!!! inverse li need
{

	stack[head - n - 1] = li;
	return;
}

void StackLI::Nul()
{
	for (int i = 0; i < head; ++i)
		pop();
	head = 0;

}

void StackLI::inverse()
{
	StackLI dop = (*this);
	for (int i = 0; i < head; ++i)
		dop.push((*this).pop());
	(*this) = dop;
	return;
}

LI::LI()
{
	active_n = 0;
	digits = new int [max_n];
	for (int i = 0; i < max_n; ++i)
		digits[i] = 0;
}

LI::LI(const LI & old)
{
	active_n = old.active_n;
	digits = new int[max_n];
	for (int i = 0; i < max_n; ++i)
		digits[i] = old.digits[i];
}



LI :: ~LI()

{
	if (digits != NULL)
		delete  []digits;
	digits = NULL;
	active_n = 0;
}

void LI::From_Int(const long long x)
{
	if (x == 0) 
	{
	active_n = 1;
	return;
	}
	active_n = 0;
	long long int a;
	int b;
	a = x;
	while (a != 0)
	{
		b = a % p;
		digits[active_n] = b;
		a /= p;
		++active_n;
	}
	return;
}

void LI::Print()
{
	long long int res = 0, k = 1;
	int a;
	for (int i = 0; i < active_n; ++i)
	{
		a = digits[i];
		res += a*k;
		k *= p;
	}
	cout << res << endl;
}


void LI::Nul()
{
	for (int i = 0; i < max_n; ++i)
		digits[i] = 0;
	active_n = 0;
}


void LI::Shift(const int & m)
{
	for (int i = 1; i <= active_n; ++i)
		digits[active_n + m - i] = digits[active_n - i];
	for (int i = 0; i < m; ++i)
		digits[i] = 0;
	active_n += m;
}


void LI::Shift_left(int m)
{
	for (int i = 0; i < active_n; ++i)
		digits[i] = digits[m + i];
	for (int i = 1; i <= m; ++i)
		digits[active_n - i] = 0;	
	
	active_n -= m;
}

int LI::ident_max_n(const LI & a, const LI & b)
{
	if (a.active_n >= b.active_n)
		return a.active_n;
	return b.active_n;
}

int LI::ident_min_n(const LI & a, const LI & b)
{
	if (a.active_n <= b.active_n)
		return a.active_n;
	return b.active_n;
}




const LI LI::Karatsuba(LI & b)
{
	LI res;

	if (active_n <= MIN_KARATSUBA || b.active_n <= MIN_KARATSUBA)
		return (b*(*this));
	else
	{
		
		int N = ident_min_n((*this), b);
		int m = N / 2, Nu, Nv;

		Nu = (*this).active_n - m; // менш значуща частина числа a
		Nv = b.active_n - m; // більш значуща частина числа б
		LI U, u, V, v;
		LI c1, c2, c3, c4, c5;

		// заповнення чисел U u V v
		for (int i = 0; i < active_n - m; ++i, ++U.active_n)
			U.digits[i] = digits[m + i];
		
		for (int i = 0; i < b.active_n - m; ++i, ++V.active_n)
			V.digits[i] = b.digits[m + i];
		
		for (int i = 0; i < m; ++i, ++u.active_n)
			u.digits[i] = digits[i];
		
		for (int i = 0; i < m; ++i, ++v.active_n)
			v.digits[i] = b.digits[i];
		
		c1 = u.Karatsuba(v);
		c2 = U.Karatsuba(V);
		c3 = (u + U).Karatsuba((v + V));
		
		c5 = c3 - c2; 
		c4 = c5 - c1; 

		c2.Shift(2 * m);
		c4.Shift(m);

		res = c2 + c4 + c1;

		// підрахунок формули
		res.active_n = (*this).active_n + b.active_n;
				
		while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
			--res.active_n;
		return res;
	}
}
/*
const LI LI::T_K(const LI & b)
{
	LI u = (*this), v = b, w;
	int n = ident_max_n (*this, b);
	u.active_n = n;
	v.active_n = n;
	/// dopisat' 0
	
	StackLI C, U, V, W, Z, dop_stack;
	StackInt stan;
	long int q[11], r[11];
	int K, Q, R, s, t, pTK, A;
// step 1
	C.Nul();
	U.Nul();
	V.Nul();
	W.Nul();
	K = 1;
	q[0] = p*p*p*p;
	q[1] = q[0];
	r[0] = p*p;
	r[1] = r[0];
	Q = p*p;
	R = p;

// step 2
	while (K < 10 && q[K - 1] + q[K] < n)
	{	
		++K;
		Q += R;
		if ((R + 1)*(R + 1) <= Q)
			++R;
		q[K] = pow (p, Q);
		r[K] = pow (p, R);

	}
	if (K == 10)
	{
		cout << "Your number is too long \n";
		int j; cin >> j; exit(1); system("pause");
	}
	// q[k - 1] + q[k] >= n; prodovzenja

	// ????? k = K on site
	C.push(v);
	C.push(u); // may be swap
	stan.push(0); // 0 = stop

// step 3 
	while (stan.get_head() > 0)
	{
// step  4
		while (K > 1)
		{
// step 5
			--K;
			s = q[K];
			t = r[K];
			pTK = q[K - 1] + q[K];
// step 6
			LI timeLI;
			timeLI = C.pop();
			timeLI.divide(dop_stack, s, t);
			LI u_st = dop_stack.pop(), u_st1;			
			U.push(u_st);
			for (int i = 1; i < 2 * t + 1; ++i)
			{
				u_st1.From_Int(0);
				for (int a = 0; a < t + 1; ++a)				
					u_st1 = u_st1 + dop_stack.pick(dop_stack.get_head() - 1 - a) * i;
				u_st1 = u_st1 + u_st;
				U.push(u_st1);
			}

// step 7
			dop_stack.Nul();
			timeLI = C.pop();
			timeLI.divide(dop_stack, s, t);
			LI v_st = dop_stack.pop(), v_st1;
			V.push(v_st);
			for (int i = 1; i < 2 * t + 1; ++i)
			{
				v_st1.From_Int(0);
				for (int a = 0; a < t + 1; ++a)
					v_st1 = v_st1 + dop_stack.pick(dop_stack.get_head() - 1 - a) * i;
				v_st1 = v_st1 + v_st;
				V.push(v_st1);
			}

// step 8
			for (int i = 0; i <= 2 * t; ++i)
			{
				C.push(V.pop());
				C.push(U.pop());
			}
			stan.push(2); // 2 = interpolate
			for (int i = 0; i < 2 * t; ++i)
				stan.push(1); // 1 = save
		}

// step 9
		K = 0;
		u = C.pop();
		v = C.pop();
		w = u * v;

// step 10
		A = stan.pop();
		if (A == 2)
		{
// step 11
			W.push(w);
			s = q[K];
			t = r[K];
			pTK = q[K]; // + q[K - 1], but K == 0 :'(
			Z.Nul();
			w.divide(Z, s, 2*t); // было просто т
		//	dop_stack.inverse(); // for Z

// step 12
			LI timeLI, timeLI1;
			for (int i = 1; i < p; ++i) // may be i < 2
// step 13
				for (int j = 2 * t; j >= i; --j) // nema 2t-1, prosto 2t
				{
					int j1 = j-1;
					timeLI = Z.pick(j);
					timeLI1 = Z.pick(j1); ///  -1
					timeLI = timeLI - timeLI1;
					timeLI = timeLI / i;
					Z.change(j, timeLI);
				}
				Z.Print(); // DELETE

			int i = 2 * t;
			timeLI = (Z.pick(i) - Z.pick(i - 1)) / i;
			Z.change(i, timeLI);
			Z.Print(); // DELETE

// spet 14
			for (int i = 2 * t - 1; i > 0; --i)
// step 15
			for (int j = i; j < 2 * t; ++j)
			{	
				timeLI = (Z.pick(j) - Z.pick(j + 1) * i);
				Z.change(j, timeLI);
			}
			W.pop();

// step 16
			LI res;
			res.From_Int(0);
			for (int j = 2 * t; j > 0; --j)
			{
				timeLI = Z.pick(j);
				timeLI.Shift(s);
				res = res + timeLI;
			}
			res = res + Z.pick(0);
			Z.Nul();
			++K;
		}
// step 17
		else
		{
			if (A == 0) // 0 = stop
				return w;
		}
// step 18
		if (A == 1) // 1 = save	
			{
				++K;
				W.push(w);
			}
		else
		{
			cout << "A != save\n";  	int j; cin >> j; exit(1); system("pause");
		}
	}

/*	if (stan.get_head() == 0)
	{
		cout << "Stan stack is empty\n";

		int j; cin >> j;  exit(1); system("pause");
	}
	

		






}*/

void LI::divide(StackLI & st, int s, int t) // розділяе число на т+1 частину і записує у "правильному" напрямку
{
	LI timeli;
	int r;
	int act = active_n;
	active_n = s * (t + 1);
	for (int j = 0; j < t + 1; ++j)
	{
		
		
		timeli.Nul();
		for (int i = 0; i < s; ++i)
		{
			r = digits[active_n - 1 - i - 3 * j];
			timeli.digits[i] = r;
			++timeli.active_n;
		}
		int k = 0;
		int act_time = timeli.active_n;
		for (; act_time > 1 && timeli.digits[k] == 0; ++k, --act_time){}

		timeli.Shift_left(k);       
		st.push(timeli);
	}
	active_n = act;
	//delete []timeli.digits;
	//timeli.digits = NULL;
}

void LI::Print_Native()
{
	for (int i = active_n - 1; i >= 0; --i)
		cout << digits[i] << "  ";
	cout << endl;
}

void LI::Print_Native_Inverse()
{
	for (int i = 0; i < active_n; ++i)
		cout << digits[i] << "  ";
	cout << endl;
}


LI LI :: operator = (const LI & old)
{
	(*this).digits = new int[max_n];
	(*this).active_n = old.active_n;
	for (int i = 0; i < max_n; ++i)
		(*this).digits[i] = old.digits[i];
	return (*this);
}

LI LI :: operator + (const LI & second)
{
	LI res;
	int N = ident_max_n(*this, second);
	int z, c, ost = 0;

	for (int i = 0; i < N || ost; ++i)
	{
		z = (*this).digits[i] + second.digits[i] + ost;
		c = z % p;
		ost = z / p;
		res.digits[i] = c;
		++res.active_n;
	}
	return res;
}

LI LI :: operator * (const LI & second)
{
	LI res;
	int c, z;
	res.active_n = active_n + second.active_n;
	for (int i = 0; i < second.active_n; ++i)
	for (int j = 0; j < active_n; ++j)
	{
		c = digits[j] * second.digits[i]; // добуток
		z = c % p; // остача, що сумуеться
		c /= p; // остача від множення, що переноситься
		res.digits[i + j] += z; // сумування
		z = res.digits[i + j] / p; // остача від додавання, що переноситься			
		res.digits[i + j + 1] += z + c; // перенос остачі
		res.digits[i + j] %= p; // виділення модуля р
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0) // занулення неактивних розрядів
		--res.active_n;
	return res;
}

const LI LI :: operator * (int & second) // тут коротке чиссло за замовчуванням мае туж систему числення, що і довге!!!
{
	LI res;
	int c;
	int k = 1;
	int sec = second;
	while ((sec /= p) != 0)
		++k;
	res.active_n = (*this).active_n + k;
	for (int i = 0; i < (*this).active_n; ++i)
	{
		c = (*this).digits[i] * second;
		res.digits[i] += c%p;
		res.digits[i + 1] += c / p;
		res.digits[i + 1] += res.digits[i] / p;
		res.digits[i] %= p;
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0) // занулення неактивних розрядів
		--res.active_n;
	return res;
}

const LI LI :: operator / (const int & second)
{
	LI res;
	LI zero;
	zero.From_Int(0);
	if (second == 0)
	{
		cout << "Div by 0\n";
		exit(1);
	}
	if ((*this).active_n == 1)
	{
		int r = (*this).digits[0] / second;
		res.From_Int(r);
		return res;
	}
	
	res.active_n = (*this).active_n;
	int ost = 0, dig;
	for (int i = 0; i < (*this).active_n; ++i)
	{
		dig = ost*p + (*this).digits[(*this).active_n - 1 - i];
		res.digits[(*this).active_n - 1 - i] = dig / second;
		ost = dig - res.digits[(*this).active_n - 1 - i] * second;
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;
	return res;
}


const int LI :: operator % (const int & second)
{
	if (second == 0)
	{
		cout << "Div by 0\n";
		exit(1);
	}
	LI res;
	res.active_n = (*this).active_n;
	int ost = 0, dig;
	for (int i = 0; i < (*this).active_n; ++i)
	{
		dig = ost*p + (*this).digits[(*this).active_n - 1 - i];
		res.digits[(*this).active_n - 1 - i] = dig / second;
		ost = dig - res.digits[(*this).active_n - 1 - i] * second;
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;
	
	return ost;
}

const LI LI :: operator / (const LI & second)
{
	LI res, ost;

	if ((*this).active_n < second.active_n)
	{
		res.From_Int(0);
		return res;
	}

	// якщо дільник невеликий (вміщаеться в базу) 
	if (second.active_n == 1)
	{
		res = (*this) / (second.digits[0]);
		return res;
	}

	// тимчасовий додатковий масив
	LI U = (*this), B = second;
	U.active_n = (*this).active_n + 1;
	long n = B.active_n, m = U.active_n - B.active_n;
	long uJ, vJ, i;
	long temp1, temp2, temp;
	int scale; // коеф нормализації

	int qGuess, r; // підбір частки і остача для неї
	int borrow, carry; // переноси

	// Norm
	scale = p / (B.digits[n - 1] + 1);
	if (scale > 1)
	{
		U = U * scale;
		B =  B * scale;
	}

	for (vJ = m, uJ = n + vJ; vJ >= 0; --vJ, --uJ)
	{
		qGuess = (U.digits[uJ] * p + U.digits[uJ - 1]) / B.digits[n - 1];
		r = (U.digits[uJ] * p + U.digits[uJ - 1]) % B.digits[n - 1];

		while (r < p)
		{
			temp2 = B.digits[n - 2] * qGuess;
			temp1 = r * p + U.digits[uJ - 2];

			if ((temp2 > temp1) || (qGuess == p))
			{ // умови не виконані, зменшити і порахувати остачу
				--qGuess;
				r += B.digits[n - 1];
			}
			else break;
		} // 1!
		// частку вгадана правильно або на 1 більша за таку

		carry = 0; borrow = 0;
		int *uShift = U.digits + vJ;

		// po B
		for (i = 0; i < n; ++i) // 2!
		{// отримати в темп1 цифру добутку B*qGuess
			temp1 = B.digits[i] * qGuess + carry;
			carry = temp1 / p;
			temp1 -= carry * p;

			temp2 = uShift[i] - temp1 + borrow;
			if (temp2 < 0)
			{
				uShift[i] = temp2 + p;
				borrow = -1;  // 3!
			}
			else
			{
				uShift[i] = temp2;
				borrow = 0;
			}
		}
		temp2 = uShift[i] - carry + borrow; // 4!
		if (temp2 < 0)
		{
			uShift[i] = temp2 + p;
			borrow = -1;
		}
		else
		{
			uShift[i] = temp2 ; //5!
			borrow = 0;
		}

		// is OK???
		if (borrow == 0) // Yesss
			res.digits[vJ] = qGuess;
		else // no, qGuess -= 1
		{
			res.digits[vJ] = qGuess - 1;
			// dodaty 
			carry = 0;
			for (i = 0; i < n; ++i) //6!
			{	
				temp = uShift[i] + B.digits[i] + carry;
				if (temp >= p)
				{
					uShift[i] = temp - p;
					carry = 1;
				}
				else { // 7!
					uShift[i] = temp;
					carry = 0;
				}
			}
			uShift[i] = uShift[i] + carry - p;
		}
		i = U.active_n - 1;
		while ((i > 0) && (U.digits[i] == 0)) --i;
		U.active_n = i + 1;
	}


	if (scale > 1)
	{
		B = B / scale; // ??? B = B / scale
		ost = U / scale;
	}
	else ost = U;
	res.active_n = (*this).active_n; // переделать актив
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;
	return res;
}


const LI LI :: operator % (const LI & second)
{
	LI res, ost;

	if ((*this).active_n < second.active_n)
	{
		
		return (*this);
	}

	// якщо дільник невеликий (вміщаеться в базу) 
	if (second.active_n == 1)
	{int rr;
		rr = (*this) % (second.digits[0]);
		res.From_Int(rr);
		return res;
	}

	// тимчасовий додатковий масив
	LI U = (*this), B = second;
	U.active_n = (*this).active_n + 1;
	long n = B.active_n, m = U.active_n - B.active_n;
	long uJ, vJ, i;
	long temp1, temp2, temp;
	int scale; // коеф нормализації

	int qGuess, r; // підбір частки і остача для неї
	int borrow, carry; // переноси

	// Norm
	scale = p / (B.digits[n - 1] + 1);
	if (scale > 1)
	{
		U = U * scale;
		B = B * scale;
	}

	for (vJ = m, uJ = n + vJ; vJ >= 0; --vJ, --uJ)
	{
		qGuess = (U.digits[uJ] * p + U.digits[uJ - 1]) / B.digits[n - 1];
		r = (U.digits[uJ] * p + U.digits[uJ - 1]) % B.digits[n - 1];

		while (r < p)
		{
			temp2 = B.digits[n - 2] * qGuess;
			temp1 = r * p + U.digits[uJ - 2];

			if ((temp2 > temp1) || (qGuess == p))
			{ // умови не виконані, зменшити і порахувати остачу
				--qGuess;
				r += B.digits[n - 1];
			}
			else break;
		} // 1!
		// частку вгадана правильно або на 1 більша за таку

		carry = 0; borrow = 0;
		int *uShift = U.digits + vJ;

		// po B
		for (i = 0; i < n; ++i) // 2!
		{// отримати в темп1 цифру добутку B*qGuess
			temp1 = B.digits[i] * qGuess + carry;
			carry = temp1 / p;
			temp1 -= carry * p;

			temp2 = uShift[i] - temp1 + borrow;
			if (temp2 < 0)
			{
				uShift[i] = temp2 + p;
				borrow = -1;  // 3!
			}
			else
			{
				uShift[i] = temp2;
				borrow = 0;
			}
		}
		temp2 = uShift[i] - carry + borrow; // 4!
		if (temp2 < 0)
		{
			uShift[i] = temp2 + p;
			borrow = -1;
		}
		else
		{
			uShift[i] = temp2; //5!
			borrow = 0;
		}

		// is OK???
		if (borrow == 0) // Yesss
			res.digits[vJ] = qGuess;
		else // no, qGuess -= 1
		{
			res.digits[vJ] = qGuess - 1;
			// dodaty 
			carry = 0;
			for (i = 0; i < n; ++i) //6!
			{
				temp = uShift[i] + B.digits[i] + carry;
				if (temp >= p)
				{
					uShift[i] = temp - p;
					carry = 1;
				}
				else { // 7!
					uShift[i] = temp;
					carry = 0;
				}
			}
			uShift[i] = uShift[i] + carry - p;
		}
		i = U.active_n - 1;
		while ((i > 0) && (U.digits[i] == 0)) --i;
		U.active_n = i + 1;
	}


	if (scale > 1)
	{
		B = B / scale; 
		ost = U / scale;
	}
	else ost = U;
	res.active_n = (*this).active_n; // переделать актив
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;
	return ost;
}

LI LI :: operator - ()
{
	LI res = (*this);
	for (int i = 0; i < (*this).active_n; ++i)
		res.digits[i] = -(*this).digits[i];
	return res;
}

LI LI :: operator - (const LI & second)
{
	int z;
	LI res;
	if ((*this) >= second)
	for (int i = 0; i < (*this).active_n; ++i)
	{
		z = (*this).digits[i] - second.digits[i];
		if (z < 0)
		{
			z += p;
			--(*this).digits[i + 1];
		}
		res.digits[i] = z;
		++res.active_n;
	}
	else
	{
		for (int i = 0; i < second.active_n; ++i)
		{
			z = second.digits[i] - (*this).digits[i];
			if (z < 0)
			{
				z += p;
				--(*this).digits[i + 1];
			}
			res.digits[i] = -z;
			++res.active_n;
		}
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;
	return res;
}

bool LI :: operator > (const LI & second)
{
	int N1 = (*this).active_n,
		N2 = second.active_n;
	if (N1 > N2) return true;
	if (N1 < N2) return false;
	for (int i = 1; i <= N1; ++i)
	{
		if ((*this).digits[N1 - i] == second.digits[N1 - i]) continue;
		if ((*this).digits[N1 - i] > second.digits[N1 - i]) return true;
		if ((*this).digits[N1 - i] < second.digits[N1 - i]) return false;
	}
	return false;
}

bool LI :: operator >= (const LI & second)
{
	int N1 = (*this).active_n,
		N2 = second.active_n;
	if (N1 > N2) return true;
	if (N1 < N2) return false;
	for (int i = 1; i <= N1; ++i)
	{
		if ((*this).digits[N1 - i] == second.digits[N1 - i]) continue;
		if ((*this).digits[N1 - i] > second.digits[N1 - i]) return true;
		if ((*this).digits[N1 - i] < second.digits[N1 - i]) return false;
	}
	return true;
}

void LI :: From_String(const string & str)
{
	(*this).Nul();
	int len = 0;
	for (; str[len]; ++len) {}
	int t = p, coef = 0,i, j, c;

	while (t != 1)
	{
		t /= 10;
		++coef;
	}

	for (i = 0; i < len / coef; ++i)
	{
		for (j = 0, c = 0; j < coef; ++j)
			c += (str[len - 1 - (j + i * coef)] - 48) * pow(10, j);
		(*this).digits[i] = c;
	}
	for (j = 0; j < len % coef; ++j)
		(*this).digits[i] += (str[len % coef - 1  - j] - 48) * pow (10, j);

	(*this).active_n = len/coef;
	if (len%coef > 0)
		++(*this).active_n;
	return;
}

void LI :: From_p_system_to_10()
{
	LI dop, a = *this;
	int t = p, coef = 0, i, j, c;
	while (t != 1)
	{
		t /= 10;
		++coef;
	}
	dop.active_n = active_n * coef;
	for (i = 0; i < active_n; ++i)
	for (j = 0; j < coef; ++j)
	{
		c = a.digits[i] % 10;
		a.digits[i] /= 10;
		dop.digits[i * coef + j] = c;
	}
	while (dop.digits[dop.active_n - 1] == 0)
		--dop.active_n;
	(*this) = dop;
}

void LI::Print_String()
{
	int t = p, coef = -1;

	while (t != 0)
	{
		t /= 10;
		++coef;
	}
	cout << digits[active_n - 1];
	for (int i = active_n - 2; i >= 0; --i)
	{
		int c = 0, d = digits[i];
		if (d == 0)
		{
			for (int q = 0; q < coef; ++q)
				cout << 0;
			continue;
		}
		while (d != 0)
		{
			d/=10;
			++c;
		}
		int count = coef - c;
		for (int j = 0; j < count; ++j)
			cout <<0;
		cout << digits[i];
	}
	cout << endl;
}

const LI LI ::  gcd(const LI & second)
{
	LI a = (*this), b = second;
	LI zero;
	zero.From_Int(0);
	if (a == zero) return b;
	if (b == zero) return a;
	while (b != zero)
	{
			a = a % b;
			a.swap(b);
	}
		return a;
}

const bool LI :: operator == (const LI & second)
{
	if ((*this).active_n != second.active_n)
		return false;
	int n = active_n;
	for (int i = 0; i < n; ++i)
		if ((*this).digits[i] != second.digits[i]) return false;
	return true;

}

 
 bool LI :: operator != (LI & second)
{
	return !((*this) == second);
}

void LI :: swap(LI & second)
{
	LI t = (*this);
	(*this) = second; 
	second = t;
	return;
}
/*
const int LI :: Jacobi( LI & n)
{
	cout << "in Jacobi \n";
	LI a = (*this), a1;
	int symb;
	LI zero, one;
	zero.From_Int(0);
	one.From_Int(1);
	symb = 1;
	if (a == zero) return 0;
	if (a == one) return 1;

	LI dop = a, n1;
	int i;
	for (i = 0; dop%2 == 0; ++i)
		dop = dop / 2;
	a1 = dop;

	if (i % 2 == 1) symb = 1;
	else
	{
		if ((n % 8 == 1) || (n % 8 == -1))
			symb = 1;
		if ((n % 8 == 3) || (n % 8 == -3))
			symb = -1;
	}

	if ((n % 4 == 3) && (a1 % 4 == 3))
		symb = -symb;

	n1 = n % a1;
	if (a1 == one)
		return symb;
	else return (n1.Jacobi(a1)) * symb;


}
*/
const int LI::Jacobi(const LI & n)
{
	LI a = (*this), c;
	LI zero, one;
	LI dop = n;
	zero.From_Int(0);
	one.From_Int(1);

	LI nod = a.gcd(n);
	if (nod != one) return 0;
	int symb = 1, i;

	if (!(a >= zero))
	{
		a = -a;
		
		if (dop % 4 == 3)
			symb = -symb;
	}

	do
	{
		for (i = 0; a % 2 == 0; ++i)
			a = a / 2;

		if (i % 2 == 1)
		{
			if (dop % 8 == 3 || dop % 8 == 5)
				symb = -symb;
		}

		if (a % 4 == 3 && dop % 4 == 3)
			symb = -symb;
		c = a; a = dop % a;  dop = c;
		
	}
	while (a != zero);
	return symb;

}


const LI LI :: powLI (const int & n)
{
	LI one;
	one.From_Int(1);
	LI res = one;
	for (int i = 0; i <= n; ++i)
		res = res.Karatsuba (*this);
	return res;
}

const LI LI :: powLI (LI & n)
{
	LI zero, one; 
	zero.From_Int(0);
	one.From_Int(1);
	LI res = one;
	for (LI i = zero; n > i; i = i + one)
		res = res.Karatsuba(*this);
	return res;
}

const void LI::S_S() // Solovay–Strassen primality test
{ 
	int k = 3;
	LI zero, one, two, three;
	zero.From_Int(0);
	one.From_Int(1);
	two.From_Int(2);
	three.From_Int(3);
	
	int N = (*this).active_n;
	LI a;
	srand(time(NULL));
	for (int i = 0; i < k; ++i)
	{ 
		a = three + (*this) / (1 + rand() % k);
		while (three > a) a = a + one;

		LI nod = (*this).gcd(a);
		if (nod > one)
		{
			cout << "The number is composite\n"; return;
		}

		LI st = (*this - one) / 2;
		LI mod, A, JacLI;
		int Jac;
	
		mod = a.pow_mod(st, (*this));
		Jac = a.Jacobi(*this);
		JacLI.From_Int(Jac);
		
		if (JacLI == -one ) JacLI = (*this) - one;

		if (mod != JacLI)
		{
			cout << "The number is composite\n"; return;
		}
	}
	cout << "The number is probably prime\n"; return;
}

const LI LI :: pow_mod(LI & n, const LI & mod)
{
	LI res, zero, one, two, a = (*this);
	zero.From_Int(0);
	one.From_Int(1);
	two.From_Int(2);
	res = one;
	while (n != zero)
	{
		LI d;
		d.From_Int(n%2);
		if (d == zero)
		{
			n = n /	2;
			a = (a * a) % mod;
		}
		else
		{
			n = n - one;
			res = (res * a) % mod;
		}
	}
	return res;
}

const void LI::R_M() //Miller–Rabin primality test
{
	int k = 10;
	
	LI t, s, a, x;
	LI zero, one, two, three;
	zero.From_Int(0);
	one.From_Int(1);
	two.From_Int(2);
	t = zero; s = zero;
	three.From_Int(3);
	srand(time(NULL));
	bool sycl ;
	if ((*this) % 2 == 0)
	{
		cout << "The number is composite\n";	return;
	}
	s = zero;
	
	LI u = (*this) - one;
	while (u % 2 == 0)
	{
		u = u / 2;
		s = s + one;
	}
	t = u;
	
	for (int i = 0; i < k; ++i)
	{
		sycl = true;
		//a.From_Int(3);
		a = u / (1 + rand() % k);
		while (three > a) a = a + one;
		x = a.pow_mod(t, (*this));
		if (x == one || x == (*this) - one) continue;
		for (LI q = zero; !(q >= s - one); q = q + one)
		{
			x = (x * x) % (*this);
			if (x == one)
			{
				cout << "The number is composite\n"; return;
			}
			if (x == (*this) - one)
			{
				sycl = false;
				break;
			}
							
		} if (sycl) { cout << "The number is composite\n"; return; }
		
	}
	cout << "The number is probably prime\n"; return;
}

const void LI::L_L() // Lucas–Lehmer primality test
{
	LI s, zero, one, two, M, pL = (*this), limit;
	zero.From_Int(0);
	one.From_Int(1);
	two.From_Int(2);
	limit.From_Int(3000);
	double time1, time2, time;
	time1 = clock();
	if (pL > limit)
	{
		cout << "Marsen's index is too long, enter enother [3, 3000]\n";
		return;
	}
	s = two + two;
	M = two.powLI(pL);
	M = M - one;
	
	LI sycl = pL - two;
	for (LI i = zero; sycl > i; i = i + one)
		s = (s * s - two) % M;
	if (s == zero) 
	{
		time2 = clock();
		time = time2 - time1;
		cout << "The number  "; M.Print_String(); cout <<" is prime\n";
		cout << "Time = " << time << " msec\n";
		return;
	}
	else
	{
		time2 = clock();
		time = time2 - time1;
		cout << "The number  "; M.Print_String(); cout << " is composite\n"; 
		cout << "Time = " << time << " msec\n";
		return;
	}
}


void LI :: From_String_Double (const string & str)
{
	(*this).From_String(str);
	active_n -= 2;
	digits[active_n] = digits[active_n + 1] = 0;
}

void LI :: From_10_to_binary()
{
LI dop, a = *this, zero;
zero.From_Int(0);
int b;
while (a != zero)
{
	b = a.binmod(2);
	dop.digits[dop.active_n ++] = b;
	a = a.bindiv(2);
}
(*this) = dop;
}


const int LI :: binmod(const int & second)
{
	if (second == 0)
	{
		cout << "Div by 0\n";
		exit(1);
	}
	LI res;
	res.active_n = (*this).active_n;
	int ost = 0, dig;
	for (int i = 0; i < (*this).active_n; ++i)
	{
		dig = ost*10 + (*this).digits[(*this).active_n - 1 - i];
		res.digits[(*this).active_n - 1 - i] = dig / second;
		ost = dig - res.digits[(*this).active_n - 1 - i] * second;
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;

	return ost;
}


const LI LI :: bindiv(const int & second)
{
	LI res;
	LI zero;
	zero.From_Int(0);
	if (second == 0)
	{
		cout << "Div by 0\n";
		exit(1);
	}
	if ((*this).active_n == 1)
	{
		int r = (*this).digits[0] / second;
		res.From_Int(r);
		return res;
	}

	res.active_n = (*this).active_n;
	int ost = 0, dig;
	for (int i = 0; i < (*this).active_n; ++i)
	{
		dig = ost*10 + (*this).digits[(*this).active_n - 1 - i];
		res.digits[(*this).active_n - 1 - i] = dig / second;
		ost = dig - res.digits[(*this).active_n - 1 - i] * second;
	}
	while (res.active_n > 1 && res.digits[res.active_n - 1] == 0)
		--res.active_n;
	return res;
}
/*
const LI LI :: invCook(const string & str)
{
	LI a, z, w;
	int v1, v2, v3;
	w.From_Int(32);
	
	//a.From_String_Double(str);
	//a.From_p_system_to_10();
	//a.From_10_to_binary();
	a.From_String_bin(str);
	if (a.digits[a.active_n - 1] != 1)
	{
		cout << "The wrong number \n"; exit(1);
	}
	v1 = a.digits[a.active_n - 1];
	v2 = a.digits[a.active_n - 2];
	v3 = a.digits[a.active_n - 3];

	z = (w / (v1 * 4 + v2 * 2 + v3)) 
	return a;
}
*/
void LI::From_String_bin(const string & str)
{
	int len  = 0;
	for (; str[len]; ++len) {}
	len -= 2;
	for (int i = 0; i < len; ++i)
		(*this).digits[len - 1 - i] = str[i];
	(*this).active_n = len;
}

typedef complex<double> cd;
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

LI LI :: Strassen( LI & second)
{
	
	int n = 1;
	int N = ident_max_n(*this, second);
	while (n <= N)
		n *= 2;

	vcd a, b;
	a.resize(n);
	b.resize(n);
	for (int i = 0; i < n; ++i)
	{
		a[i] = (*this).digits[i];
		b[i] = second.digits[i];
	}
	vcd temp;
	temp.resize(n);
	fft(a, false);
	fft(b, false);

	for (int i = 0; i < n; ++i)
	{
		temp[i] = a[i] * b[i];
	}
	fft(temp, true);

	LI res;
	for (int i = 0; i < n; ++i)
	{
		res.digits[i] = (int)(temp[i].real() + 0.5);
	}
	
	int ost = 0;
	for (res.active_n = 0; ost != 0 || res.active_n < temp.size(); ++res.active_n)
	{
		int c, z;
		ost = res.digits[res.active_n] / p;
		z = res.digits[res.active_n] % p;
		res.digits[res.active_n] = z + ost;
		
	}

	return res;
}

