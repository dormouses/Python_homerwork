from __future__ import print_function
def dec(n):
  if (n>10):
    dec(n/10)
  print(n%10, end='')

n=input()
dec(n)

