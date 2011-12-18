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

def ExtendedGCD(a, b):
	'''
		Из книги Т. Корман <<Алгоритмы. Построение и Анализ>>
		стр 966
	'''
	if b == 0:
		return a , 1, 0
	d1, x1, y1 = ExtendedGCD(b, mod(a, b) )
	d, x, y = d1, y1, x1 - (a//b)*y1
	return d, x, y

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

def EulerPhi(num):
	'''
		Функция Эйлера
		\varphi(n), где n — натуральное число, 
		равна количеству натуральных чисел, 
		не больших n и взаимно простых с ним.
	'''
	if (num>13):
		if Miller_Rabin(num, 9000):
			return num-1
	
	res = 1
	i = 2
	while( i*i <= num):
		# пока i^2 <= num
		p = 1
		while(num % i == 0):
			num //= i 		# если не взаимно просты, делим
			p *= i			# произведение делителей i втч и кратных	
		p //= i
		
		if ( p != 0 ):
			# если мы хоть раз делили на текущее i
			# то общее произведение делителей
			# уиножаем на (i - 1)*i^(число раз - 1)
			res *= ( p * (i - 1))
		i += 1
	
	if(num == 1):
		return res
	else:
		return (num - 1) * res
		
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

#======================
def Miller_Rabin(n, k):
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

def mod (a, b):
	return a % b

def get_inv(a, b):
	for i in range(0, min(b, 9000)):
		if ( (a*i)%b == 1): return i
	return 0

def build_abx(alpha, beta, p, phi_n, u,v,x):
	m = (x+1)%3
	if  (m == 0):
		u = mod(2*u, phi_n)
		v = mod(2*v, phi_n)
		x = mod(x*x, p)
	elif(m == 1):
		u = mod(u+1, phi_n)
		x = mod(x*alpha, p)
	elif(m == 2):
		v = mod(v+1, phi_n)
		x = mod(x*beta, p)
	else:
		print("FISK EE")
	return u,v,x
	
def generate1(a, b, p, lim = 900500):
	""" a**x == b (mod p) """
	u = 0
	v = 0
	z = 1
	U = u
	V = v
	Z = z

	phi_n = EulerPhi(p)

	for i in range(lim):
		u,v,z = build_abx(a,b,p,phi_n, u,v,z)
		U,V,Z = build_abx(a,b,p,phi_n, U,V,Z)
		U,V,Z = build_abx(a,b,p,phi_n, U,V,Z)
		#~ print(i+1,u,v,z,U,V,"\t",Z)

		if (z == Z):
			r = (v - V) % p
			print( "\nz[",i,"] =", z,",", Z)
			vx = mod( (V - v), phi_n)
			ux = mod( (u - U), phi_n)
			print("b^vx =", pow(b, vx, p), "; a^ux =", pow(a ,ux, p) )

			x = solve_diofant(ux, vx, p, a, b)
			return x

	print("Over time")

	

def generate2(a, b, p, lim = 900500):
	""" a**x == b (mod p) """
	u = [0]
	v = [0]
	z = [1]
	z1 = [1]
	
	phi_n = EulerPhi(p)
	
	for i in range(lim):
	#	g = (z[-1]*3)//p
		if  (z[-1] <= 1*p/3): # z[-1]
			u.append( mod(u[-1]+1, phi_n) )
			v.append( mod(v[-1],   phi_n) )
		elif(z[-1] <= 2*p/3):
			u.append( mod(2*u[-1], phi_n) )
			v.append( mod(2*v[-1], phi_n) )
		elif(z[-1] <= 3*p/3):
			u.append( mod(u[-1],   phi_n) )
			v.append( mod(v[-1]+1, phi_n) )
		else:
			print("FAIL")
		z.append( (pow(b, u[-1], p) * pow(a, v[-1], p)) % p )

		if(i%2==0)and(z[i] == z[i//2])and(i>2):
			print( "\nz[",i,"] =", z[i//2],",", z[i])
			delta = u[i]-u[i//2]

			ux = mod(v[i] - v[i//2], phi_n)
			vx = mod(u[i//2] - u[i], phi_n)
			print("b^vx =", pow(b, vx, p), "; a^ux =", pow(a, ux, p) )
			
			x = solve_diofant(ux, vx, p, a, b)
			return x

	print("Over time")

	
def solve_diofant(ux, vx, p, a, b):
	
	phi_n = EulerPhi(p)
	
	if( mod(vx , phi_n) == 0 ): return False
	d, nu, mu = ExtendedGCD(vx , phi_n)
	# nu = v^(-1) mod n
	'''
		a^ux == b^vx   mod p
		a^(ux*nu) = b^(v*nu)
		b ^(d - phi_n*nu) = b^(d) = g(x*d) mod p
		x*d =  ux*nu + w*phi_n
	'''
			
	for w in range(d+1):
		x = mod( (ux * nu + w * phi_n)//d, phi_n)
		if( pow(a, x, p) == b):
			print("solved by extend_GCD")
			return x

	# delete this
	l = get_inv(ux, phi_n)
	x = mod(l*vx, phi_n)
	if( pow(a, x, p) == b):
		print("solved by inverse")
		return x

	return False
	
				
def premutive_log(a, b, n, lim = 100500):
	z = []
	for x in range(min(n, lim)):
		if(pow(a,x,n) == b%n):
			z.append(x)
	return z


def main():

	if (len(sys.argv)==1):
		a = 13
		p = 2148568033  # 41*41
		xx = 917771
		b = pow(a, xx, p)
		print("a =", a)
		print("b =", b)
		print("n =", p)


	elif("h" in sys.argv[1]):
		print("""This programm use Pollard's pho algorithm
for solving discrete logarithm problem.

-help for this text
a b n - solve a^x == b (mod n)
none - standart test""")
		print(sys.argv)
		return
	elif(len(sys.argv) == 4):
		a = int(sys.argv[1])
		b = int(sys.argv[2])
		p = int(sys.argv[3])
		
	else:
		print("error")
		return
		
	x1 = generate1(a, b, p)
	x2 = generate2(a, b, p)
	
	x = premutive_log(a, b, p)
	
	print()
	print("x  =", x[:7])
	print("x1 =", x1)
	print("x2 =", x2)
	
	return

if (__name__ == "__main__"):
	main();
