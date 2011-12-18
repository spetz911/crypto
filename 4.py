#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

import sys
from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce
from random import randrange


SHIFT = 25

def get_inv(a, b = (2**16) +1):
	if (a == 0): return 0
	for i in range(0, b):
		if ( (a*i)%b == 1): return i
	return None

def get_sub(a, b = 2**16):
	return (b - a) % b 

def mod(a, b):
	return a%b

def low16(x):
	return ((x) & 0xFFFF)

def low128(x):
	return ((x) & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF)

def println(X):
	for xx in X: print16(xx)

def print64(X):
	for xx in X: print_64(xx, end = " ")
	print()

def print16(X):
	for xx in X: print_16(xx, end = " ")
	print()

def print_128(x, end = "\n"):
	h = hex(x)[2:]
	s = '0'*(32-len(h)) + h
	print(s, end=end)

def print_64(x, end = "\n"):
	h = hex(x)[2:]
	s = '0'*(8-len(h)) + h
	print(s, end = end)

def print_16(x, end = "\n"):
	h = hex(x)[2:]
	s = '0'*(4-len(h)) + h
	print(s, end = end)
	
def mul(a,b):
	return mod( a * b, 2**16 + 1)

def add(a,b):
	return mod( a + b, 2**16)

def ExpandKey(userkey):
	EK = []
	# mixing keys
	for j in range(8):
		EK.append( low16(userkey>>(16*j)) )
	#step 2
	userkey = low128(userkey << SHIFT) | low128(userkey >> (128-SHIFT))
	for j in range(2,8):
		EK.append( low16(userkey>>(16*j)) )
	for j in range(2):
		EK.append( low16(userkey>>(16*j)) )
	#step 3
	userkey = low128(userkey << SHIFT) | low128(userkey >> (128-SHIFT))
	for j in range(6,8):
		EK.append( low16(userkey>>(16*j)) )
	for j in range(6):
		EK.append( low16(userkey>>(16*j)) )
	#step 4
	userkey = low128(userkey << SHIFT) | low128(userkey >> (128-SHIFT))
	for j in range(0,8):
		EK.append( low16(userkey>>(16*j)) )
	#step 5
	userkey = low128(userkey << SHIFT) | low128(userkey >> (128-SHIFT))
	for j in range(4,8):
		EK.append( low16(userkey>>(16*j)) )
	for j in range(4):
		EK.append( low16(userkey>>(16*j)) )
	#step 6
	userkey = low128(userkey << SHIFT) | low128(userkey >> (128-SHIFT))
	for j in range(6,8):
		EK.append( low16(userkey>>(16*j)) )
	for j in range(6):
		EK.append( low16(userkey>>(16*j)) )
	#step 7
	userkey = low128(userkey << SHIFT) | low128(userkey >> (128-SHIFT))
	for j in range(2,8):
		EK.append( low16(userkey>>(16*j)) )

	ZK = [ (EK[6*i:6*(i+1)]) for i in range(8)]
	ZK.append(EK[48:52])

	return ZK

def InverseKey(userkey):
	ZK = ExpandKey(userkey)
	IK = []
	
	for i in range(8, 0, -1):
		IK.append( ZK[i][0:4] + ZK[i-1][4:6] )
	IK.append( ZK[0][0:4] )
	
	for i in range(9):
		IK[i][0] = get_inv( IK[i][0] )
		IK[i][1] = get_sub( IK[i][1] )
		IK[i][2] = get_sub( IK[i][2] )
		IK[i][3] = get_inv( IK[i][3] )
		#swap
		if(0<i<8): (IK[i][1], IK[i][2]) = (IK[i][2], IK[i][1])
	
	return IK
		
		
def IDEA(D, K):
	D = copy(D) #does not matter
	for i in range(8):
		a = mul(D[0], K[i][0])
		b = add(D[1], K[i][1])
		c = add(D[2], K[i][2])
		d = mul(D[3], K[i][3])
		e = a ^ c
		f = b ^ d
		g = mul( add( f, mul(e, K[i][4])) , K[i][5]) 
		h = add( mul(e, K[i][4]), g)
		D[0] = a ^ g
		D[1] = c ^ g
		D[2] = b ^ h
		D[3] = d ^ h
	a = mul(D[0], K[8][0])
	b = add(D[2], K[8][1])
	c = add(D[1], K[8][2])
	d = mul(D[3], K[8][3])
	return [a,b,c,d]

def get_word(word):
	print("word=",word[-16:])
	
	if (word[:2]=='0x'):
		v = int(word, 16)
	else:
		v = int(word, 10)
	Z = []
	for i in range(4):
		Z.append( low16(v>>i*16) )
	if(Z[-1] == 0):
		print("uncorrect word")
		return None
	Z.reverse()
	return Z
	
def get_key(word):
	print("key =", word[-32:])
	if (word[:2]=='0x'):
		v = int(word, 16)
	else:
		v = int(word, 10)
	if(v< (1<<96)):
		print("uncorrect key")
		return None
	return low128(v)

	#low128(v)


		
def main():
	if (len(sys.argv)==1):
		word = [10,11,12,13]
	elif("h" in sys.argv[1]):
		print("""International Data Encryption Algorithm
-help for this text
-encode key128 text64
-decode key128 text64
text64 - test with userr text
none - standart test""")
		print(sys.argv)
		return
	elif("en" in sys.argv[1])and(len(sys.argv) == 4):
		print("encoding")

		k    =  get_key(sys.argv[2])
		word = get_word(sys.argv[3])
		if(k==None)or(word==None): return

		key = ExpandKey(k)
		res = IDEA( word, key)

		print("IDEA:")
		print16(res)
		return
	elif("dec" in sys.argv[1])and(len(sys.argv) == 4):
		print ("decoding")

		k    = get_key(sys.argv[2])
		word = get_word(sys.argv[3])
		
		if(k==None)or(word==None): return
	
		key = InverseKey(k)
		res = IDEA( word, key)

		print("IDEA^-1:")
		print16(res)
		return		
	elif(len(sys.argv) == 2):
		word = get_word(sys.argv[1])
		if(word==None): return
	else:
		print( "error" )
		return

	l = [1,2,3,4,5,6,7,8]
	C = [0x11fb, 0xed2b, 0x0198, 0x6de5]
	k = sum( [ (ll <<(16*i)) for (ll,i) in zip(l, range(8))] )

	key = ExpandKey(k)
	inv_key = InverseKey(k)

	print("key 52:")
	println(key)
	print()
	print("key inverse:")
	println(inv_key)
	print()
	
	z = IDEA( word, key)
	print()

	res = IDEA( z , inv_key)
	
	print("\ntesting...")
	print("IDEA:")
	print16(word)
	
	print("RES:")
	print16(z)
	print("RES-1:")
	print16(res)


	return

if (__name__ == "__main__"):
	main();
