直接板子了

要把fft序列长度拓展到2的n次方

两个序列长度得一样才能开卷

```c++

#include<complex>
using namespace std;
using ll = long long;
typedef pair<ll, ll> P;
typedef complex<double> cp;

ll f(ll n){					//拓展到2 n次方
    ll bit=1;
    while(1<<bit < n)bit++;
    return 1<<bit;
}
void fft(cp* a, int n, int inv) {
    ll bit=1;
    while(1<<bit < n)bit++;
    for (int i = 0; i < n; i++) {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
        if (i < rev[i])swap(a[i], a[rev[i]]);
    }
    for (int mid = 1; mid  < n; mid *= 2) {
        cp tmp(cos(pi / mid*1.0), inv * sin(pi / mid*1.0));
        for (int i = 0; i < n; i += mid * 2) {
            cp omega(1, 0);
            for (int j = 0; j < mid; j++, omega *= tmp)//只扫左半部分，得到右半部分的答案
            {
                cp x = a[i + j], y = omega * a[i + j + mid];
                a[i + j] = x + y, a[i + j + mid] = x - y;
            }
        }
    }
    if(inv==-1)
        for(int i=0;i<n;i++)a[i].real(a[i].real()/n);
}
```

那么求 a，b卷积就可以(a长度n，b m)

```c++
ll lim=n+m+1;
ll mx=f(lim);
ll c[N];

fft(a,mx,1);
fft(b,mx,1);
for(int i=0;i<mx;i++)a[i]*=b[i];
fft(a,mx,-1);
for(int i=0;i<lim;i++)c[i]=int(a[i].real()+0.5);
```



然后用法啥的：

比如一个数列中两个数的和是否存在

我们直接给a中这些数列的数作为下标打上1，a和自己卷积，可以看出来每个点卷积结果是不是大于0就是和为该值是否存在

差可以a和a的反向序列卷积，也是差不多的意思

又比如问边数 之类的 也可以看成多项式乘法

就能用fft加速