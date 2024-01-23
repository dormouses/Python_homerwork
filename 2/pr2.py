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

def generator(n, border,  *args):
  l=[]
  for i in xrange (len(args)):
    l.extend(args[i])
  sort (l, len(l), 0, len(l)-1)
  i=0
  last=-100;
  print (l)
  while (l[0]>border) and (i < n) and (i < len(l)):
      if (l[i]!=last):
        last=l[i]
        yield l[i]
      i=i+1

l1=[15, 4, 6, 23, 22]
l2=[11, 93, 18, 123]
l3=[12, 43, 43, 125]

for val in generator(6, 10, l1, l2, l3):
    print(val)
