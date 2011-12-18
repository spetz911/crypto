#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

import sys
from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce
from random import randrange

def find_eps():
	x = (1<<100) + 521
	i = 0
	y = x*x
	while (x == int(sqrt(y)) ):
		print(x,i)
		x = x*3
		y = x*x
		i+=1
	return x 
		
def NOD( a, b):
	"""Euclidean algorithm"""
	if b == 0:
		return a
	else:
		return NOD(b, a % b)

		
# sqrt max is (1<<50)
def my_sqrt(n):
	if (n < (1<<44)) : return int(sqrt(n))
	l = (1<<22)
	r = n
	while( abs(l-r)>1):
		t = (l+r)//2
		if (t*t>n): r = t
		else: l = t
	return l

def factorization(n):
	sqrt_n = my_sqrt(n)
	z = []
	i = 2
	while(n != 1):
		while (n%i == 0):
			n = n // i
			sqrt_n = my_sqrt(n)
			z.append(i)
		i+=1
		if (i>sqrt_n):break

	if (n!=1): z.append(n)
	return z

def is_square(n):
	t = my_sqrt(n)
	return (t*t == n)

def Ferma(n, itt = 100500):
	"""Fermat's factorization method"""
	a = my_sqrt(n) + 1
	b2 = a*a - n
	while not is_square(b2):
		b2 += 2*a + 1 # equal b2 = a*a - n
		a += 1
		itt -= 1
		if (itt == 0): return None
	return a - my_sqrt(b2)
	# change to recursive

def Lehman(n, itt = 100500666137):  #not safe
	if (n<=8) or (n>100500666137): return None
	n3 = int(n**(1/3))
	
	res = [] # [ n//i;i for i in range(2, n3+2) if (n%i == 0)]
	n1 = n
	i = 2
	while (i<n3+2):
		while (n1%i == 0):
			res.append(i)
			n1//=i
		i+=1
#	pr = product(res)

#	if (product(res)==n): return res
	print (res)
#	print(product(res), n, n/pr)
	for i in res: n//=i
	n3 = int(n**(1/3))
	
	print("n3 =", n3)
	
	B = 0
	for k in range(1, n3+1):
		for d in range(1, int( n**(1/6)/(4*sqrt(k)) ) +2 ):
			A =  int(sqrt(4*k*n))+d
			#print ("A = ", A*A, k, 4*n)
			B = A*A - 4*k*n	
			if is_square( B ): break
	
	B = int(sqrt(B))
	d = NOD(A-B, n)
	if (1 < d)and(d < n): res.append(d)
	return res
	
def Pollard(n, k = 100500):
	"""Pollard's rho algorithm"""
	if( Miller_Rabin(n) ): return n
	x = y = 2
	for i in range(k):
		x = mod(x*x +1, n)
		y = mod(y*y +1, n)
		y = mod(y*y +1, n)
		d = NOD( abs(x-y), n)
		if (1 < d)and(d < n):
			return Pollard(d)
	return None


def Miller_Rabin(n, k = 256):
	"""Primality test"""
	s = 0
	d = n-1
	while d%2 == 0:
		s += 1
		d //= 2
#~	print("s = ",s, "d = ",d)
	for i in range(k):
		a = randrange(2, n-2)
		x = pow(a, d, n)
		if (x==1)or(x==n-1): continue
		for r in range(1, s):
			x = x*x % n
			if x==1   : return False
			if x==n-1 : break
		if (x==n-1) : continue
		return False
	return True

def product(L):
	mul = lambda x, y: x*y
	return reduce(mul, L)
	
def mod (a, b):
	return a % b


def main():
	if (len(sys.argv)==1):
		b = 16417377581384821* 71581384821111231 * 71582384471111233 * 1088093 * 45300907 * 1061323 * 1903619 * 1059259 * 1059323 *997
	elif("h" in sys.argv[1]):
		print("""This programm use Pollard's pho algorithm
for integer factorization.
-help for this text
-primes get some primes number
-primes lim some primes number more than lim
num_1*num_2*...*num_k - factorization of number
none - standart test %""")
		print(sys.argv)
		return
	elif(("prime" in sys.argv[1])and(len(sys.argv) == 3)):
		r = int(sys.argv[2])
		for i in range(r, r+100):
			if(Miller_Rabin(i)): print(i)
		return
	elif("prime" in sys.argv[1]):
		prime_list = [971072423, 971053357, 171086189, 171084391, 171072553, 1059323, 971070613,\
					1971062777, 1971060239,\
					2147644313, 2148793769, 2148568033,\
					1971062777, 1971060239, 1971060239, 1971060239]
		for i in range(7):
			z = randrange(len(prime_list))
			print(prime_list[z])
		return

	elif(len(sys.argv) == 2):
		b = product([ int(x) for x in sys.argv[1].split('*') ])
	else:
		print("error")
		return

	num = b
	
	print("factorization of", num)
	z = []
	primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
	for p in primes:
		while b%p == 0:
			b//=p
			z.append(p)
			print("multiplier", p )

	n = b
	while (n > 1):
		d = Pollard(n , 9000500)
		if (d == None):
			print("Fail, try to Ferma")
			d = Ferma(n, 1000500)
			if (d==None): print("Failx2"); return
			z.append( d)
			z.append( n//d)
			print("multiplier", d )
			print("multiplier", n//d )
			break
		z.append( d)
		print("multiplier", d )
		n//=d
		
	print("Check:", product(z) % 10**10, "==", num % 10**10 , (product(z) == num))

	return

if (__name__ == "__main__"):
	main();
