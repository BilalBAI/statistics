def Fib1(n): # Output the first n fibonacci series
    fib = [1, 1]
    for k in range(n):
        a = fib[-2] + fib[-1]
        fib.append(a)
    return fib

def Fib2(n): # Output fibonacci series up to n
    result = []
    a, b = 0, 1
    while(b < n):
        result.append(b)
        a, b = b, a + b
    return result

