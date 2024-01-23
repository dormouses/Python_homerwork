def pow (a, n):
  if (n==1): 
    return a
  else:
    aa=pow(a, n/2)
    aa=aa*aa
    if (n%2==0):
      return aa
    else: 
      return aa*a
n=input()
k=input()
p=pow(n, k)
print (p)

