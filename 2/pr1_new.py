import time
import random
import math

def sort (arr, m, l, r):
  av=arr[(l+r)/2];
  i=l;
  j=r;
  while (i<=j):
    while (arr[i]>av):
      i=i+1;
    while (arr[j]<av):
      j=j-1;
    if (i<=j):
      if (arr[i]<arr[j]):
        arr[i], arr[j]=arr[j], arr[i];
      i=i+1;
      j=j-1;
  if (l < j): sort(arr, m, l, j );
  if (r  > i): sort (arr, m,i ,r);

def count_time (n, m):

  l1=[];
  l2=[];
  l3=[];
  l1=[random.randint(0, 2*n) for x in xrange(m)]
  l2=[random.randint(0, 2*n) for x in xrange(m)]
  start = time.clock()
  sort(l1, m, 0, m-1)
  sort(l2, m, 0, m-1)
  p1=0
  p2=0
  for i in range(n):
    if (p1>=m):
      while (p2<=m):
        l3.append(l2[p2])
        p2=p2+1
      break
    if (p2>=m):
      while (p1<=m):
        l3.append(l1[p1])
        p2=p2+1
      break
      if (l1[p1]>l2[p2]):
        l3.append(l1[p1])
        p1=p1+1
      else:
        l3.append(l2[p2])
        p2=p2+1 
    end = time.clock()
    return end-start


f=open('results_PtorsenkoMA.txt', 'w')
for m, n in [[10, 1], [10, 2], [10, 5], [100, 2], [100, 5], [100, 10], [1000, 5], [1000, 50], [1000, 100],  [1000000, 5], [1000000, 100], [1000000, 1000]]:
  t=count_time(n, m) 
  f.write('N = '+str(n)+' '+'M = '+str(m)+' '+'time = '+str(t)+'\n')

