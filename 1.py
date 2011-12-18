#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

import sys
from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce
from random import randrange

def NOD( a, b):
	"""Euclidean algorithm"""
	if b == 0:
		return a
	else:
		return NOD(b, a % b)

def is_prime(n):
	for i in range(2,n):
		if (n % i== 0): return False
	return True

def primes(n):
	return [2] + [p for p in range(3, n, 2) if is_prime(p) ]

def product(L):
	mul = lambda x, y: x*y
	return reduce(mul, L)
	
def prima(n):
	L = []
	for i in range(2, n+1):
		k = 0
		while (n%i == 0):
			n //= i
			k+=1
		if(k>0):
			L.append( (i,k) )
		if (n == 1): break
	return L

def constr(L):
	return product( [ (p**k) for (p,k) in L ] )

def Lukas_test(L, k = 33):
	n1 = constr(L)
	n = n1+1
	for (p,m) in L:
		find_ai = False
		for i in range(k):
			a_i = randrange(2, n)
			t1 = pow(a_i, n1  , n)
			t2 = pow(a_i, n1//p, n)
			if (t1 == 1)and(t2 != 1):
				find_ai = True
				break
		if (find_ai == False): return False
	return True
	
	
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

def degenerate(q, lim, k = 77):
	""" Get prime number q_i and make q_i+1"""
	if (8*q*q > lim):
		it = range (lim//(2*q)+3, 4*q+2, 1)
		print("last iter, ",lim//(2*q))
	else:
		it = range(4*q+2, 1, -1)
	#~ MOAR bydlokod here :)
	for r in it:
		# q = randrange(2, 4*q+2)
		n = 2*r*q+1
		
		# use function here
		
		print("n =", n%(10**20), "r =", r, Miller_Rabin(n,7))
		for i in range(k):
			a = randrange(3, 2*n, 2)
			t1 = pow(a, n-1, n)
			t2 = pow(a, 2*r, n)
			if (t1==1)and(t2!=1):
				p = 2*q + 1
				if (n == p*p)and(pow(a, p-1, n)==1):
					print("AGHTUNG!")
					return p
				else:
					return n
			#olololo
	return q

def gener( p, lim, st = 1):
	""" st = speed of enlargement """
	st = 2 ** st
	z = 1
	itt = 1
	while (p[-1] < lim):
		#~ print("iter", itt)
		if (itt > 100): return 0
		itt+=1
#		if (product(p) > lim) : p.reverse()
			
		h = (1 << len(p))
		for i in range (h//2+st, h, st):
			l = get_bits(i)
			
			z = [(pp, 1) for (pp,bb) in zip(p,l) if (bb==1)]
			if ( z[0] != (2,1) ):
				z.insert( 0, (2,1) ) #palubas chetnoe
			else:
				z[0] = (2,2)
			
			n1 = constr(z) + 1
			if( Lukas_test(z) ):
				if not(Miller_Rabin(n1) ): print("FAIL"); return 0
				p.append(n1)
				break

	return n1
				
	
def get_bits( n):
	L = []
	while( n>0):
		L.append( n&1 )
		n = n >> 1
	return L 
	

def main():
	
	if (len(sys.argv)==1):
		
		n = 45300907;
		L = prima(n-1)
		print("n =", n)
		print( Lukas_test(L))
	
		primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

		lim = 10050013452347853422436578923456234954278978
		print("lim =", lim)
		n1 = gener( primes , lim)
		print("num =", n1)

		return


	elif("h" in sys.argv[1]):
		print("""Lucas primality test
-help for this text
-test (num) for test number
-gen (lim) for generate number
none - standart test""")
		return
	elif("test" in sys.argv[1]):
		if (len(sys.argv) == 3):
			n =  int(sys.argv[2])
		else:
			n = int(input())
		L = prima(n-1)
		print( Lukas_test(L))		
		return
	elif("gen" in sys.argv[1]):
		if (len(sys.argv) == 3):
			n =  int(sys.argv[2])
		else:
			n = int(input())
		primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
		n1 = gener( primes , n)
		print(n1)
		return
	else:
		print( "error" )
		return

	return

if (__name__ == "__main__"):
	main();
