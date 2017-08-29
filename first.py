#Define a factorial function
def factorial(n):
	a = range(1,n+1)
	fact = 1
	for b in a:
		fact = fact * b
		
	return fact
#...............................
number = raw_input("Enter number to calculate factorial: ")
number = int(number)
#m = 10;
print 'factorial of', number, '=', factorial(number)
