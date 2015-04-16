#include<iostream>
#include<simplecpp>
#include<conio.h>

//functions//

//addition function//
double add()
{
    double a,b,c;
    cout<<"Enter  x ";
    cin>>a;
    cout << endl;
    cout<<"Enter  y ";
    cin>>b;
    cout << endl;

    c=a+b;
    cout<<"Sum : x + y = "<<c;
    return c;

}

//subraction function//
double sub()
{
    double a,b,c;
    cout<<"Enter x ";
    cin>>a;
    cout << endl;
    cout<<"Enter y ";
    cin>>b;
    cout << endl;

    c=a-b;
    cout<<"Difference : x - y = "<<c;
    return c;
}

//multiplication function//
double multi()
{
    double a,b,c;
    cout<<"Enter x ";
    cin>>a;
    cout << endl;
    cout<<"Enter y ";
    cin>>b;
    cout << endl;

    c=a*b;
    cout<<"Product : x * y = "<<c;
    return c;
}

//division function//
double div()
{
    double a,b,c;
    cout<<"Enter x ";
    cin>>a;
    cout << endl;
    cout<<"Enter y ";
    cin>>b;
    cout << endl;

    c=a/b;
    cout<<"division : x / y = "<<c;
    return c;
}

//exponential function//
double exponentialfn(double n)
{
    double sum=0;
    double p=1;
    for (int i=0; i<10000; i++)
    {
        sum+=p;
        p*=n/(i+1);
    }
    return sum;
}

//factorial function//
long int factorial(int n)
{
    if(n==0)
    {
        return 1;
    }
    if (n<0)
    {
        cout<<"Bad input. Factorial of a negative no. is not defined";
        return -1;
    }
    if (n>0)
    {
        int m;
        m=n;
        return n*factorial(n-1);
    }
}

// ln function//
double lnfunction(double x)
{
    if(x==0)
    {
        cout<<"not defined";
    }
    else
    {
        double b,p;
        b=x;
        for(p=0; b>1; p=p+1)
        {
            b=b/10;
        }
        double c,sum=0;
        c=1-b;
        double w;
        w=c;
        for(int i=0; i<10000; i++)
        {
            sum+=w;
            w*=c*(i+1)/(i+2);
        }
        double ans;
        ans=(p-(sum/2.302585093))*2.302585093;
        return ans;
    }
}
//logarithmic function//
double logarithmicfn(double x)
{
    if(x==0)
    {
        cout<<"not defined";
    }
    else
    {
        double b,p;
        b=x;
        for(p=0; b>1; p=p+1)
        {
            b=b/10;
        }
        double c,sum=0;
        c=1-b;
        double w;
        w=c;
        for(int i=0; i<10000; i++)
        {
            sum+=w;
            w*=c*(i+1)/(i+2);
        }
        double ans;
        ans=(p-(sum/2.302585093));
        return ans;
    }
}
//conversion of angle fron degree into radians//
double radianfn(double x)
{

    return x* 3.1415926535897932384626433832795028841971693993751/180;
}

//power function//
double powerfunction(double x,double n)
{
    double a,b;
    if(x>0)
    {
        a=n*lnfunction(x);
        b=exponentialfn(a);
        return b;
    }
    int k=n;
    if(x<0&&(k==n))
    {
        if(k%2==1)
        {
            a=n*lnfunction(-x);
            b=-exponentialfn(a);
            return b;
        }
        if(k%2==0)
        {
            a=n*lnfunction(-x);
            b=exponentialfn(a);
            return b;
        }
    }
    if(x==0)
    {
        if(n>0)
        {
            return 0;
        }

    }
}


//Trignometric functions//
//sine function
double sinefunction(double x) //
{
    double p,y,sum;
    if(x>0)
    {
        for(int j=0; x>180; j++)
        {
            x=x-360;
        }
    }
    else if(x<0)
    {
        for(int k=0; x<-180; k++)
        {
            x=x+360;
        }
    }
    if(x==180||x==-180||x==0)
    {
        x=0;
    }
    y=x* 3.1415926535897932384626433832795028841971693993751/180;
    sum=0;
    p=y;
    for(int i=0; i<10000000; i+=1)
    {
        sum+=p;
        p*=-(y*y)/((2*i+2)*(2*i+3));
    }
    return sum;
}

//cosine function//
double cosinefunction(double x)
{
    float b;
    b=powerfunction(1-(sinefunction(x)*sinefunction(x)),0.5);
    if(x>0)
    {
        for(int j=0; x>180; j++)
        {
            x=x-360;
        }
    }
    else if(x<0)
    {
        for(int k=0; x<-180; k++)
        {
            x=x+360;
        }
    }
    if(x==90||x==-90)
    {
        cout<<"0";
        return 0;
    }
    else if(x>-90&&x<90)
    {
        return b;
    }
    else return -1*b;
}

//tangent function//
double tangentfunction(double x)
{
    double t;
    t=sinefunction(x)/cosinefunction(x);
    return t;
}

//Permutation and combination functions//
int npr(int n,int r)
{
    int a;
    a=factorial(n)/factorial(n-r);
//cout<<a;
    return a;

}

int ncr(int n,int r)
{
    int a;
    a=factorial(n)/(factorial(n-r)*factorial(r));
    return a;
}

//Inverse sine function//
double sineinverse(double x)
{
    double p,sum;
    p=x;
    sum=p;
    if(x<-1||x>1)
    {
        return NULL;
    }
    else
    {
        if(x==-1)return -90;
        else
        {
            if(x==1)return 90;
            else
            {
                for (int i=0; i<10000; i++)
                {
                    p*=x*x*(2*i+1)*(2*i+1)/((2*i+2)*(2*i+3));
                    sum+=p;
                }
            }
        }
    }
    return sum*180/3.141592653589793;
}


//Inverse cosine function//
double cosineinverse(double x)
{
    double b;
    b=90-sineinverse(x);
    return b;
}


//Inverse tangent function//
double taninverse(double x)
{
    double p,sum;
    p=x;
    sum=p;
    if(x==1)
    {
        return 45;
    }
    if(x==-1)
    {
        return -45;
    }
    if (x>1)
    {
        return 90-taninverse(1/x);
    }
    else
    {
        if(x<-1)
        {
            return -90-taninverse(1/x);
        }
        else if(-1<x && x<1)for (int i=0; i<10000; i++)
            {
                p*=-x*x*((2*i)+1)/((2*i)+3);
                sum+=p;
            }
        return sum*180/ 3.1415926535897932384626433832795028841971693993751;
    }
}

//Inverse function//
double inversefn(double x)
{
    double b;
    if(x==0)
    {
        cout<<"bad inputs.Inverse of zero is not defined";
        return -1;
    }
    else
    {
        b=1/x;
        return b;
    }
}

//Modulus function//
double modulusfn(double x,double y)
{
    double b;
    b=powerfunction((x*x)+(y*y),0.5);
    return b;
}

//sineh function//
double sinhfunction(double x)
{
    double b;
    b=(exponentialfn(x)-exponentialfn(-x))/2;
    return b;
}

//cosineh function//
double coshfunction(double x)
{
    double b;
    b=(exponentialfn(x)+exponentialfn(-x))/2;
    return b;
}

double tanhfunction(double x)
{
    double b;
    b=sinhfunction(x)/coshfunction(x);
    return b;
}

//Inverse of sineh function//
double sinhinversefn(double x)
{
    double k;
    k=lnfunction(x+powerfunction(x*x+1,0.5));
    return k;
}

//Inverse of cosineh function//
double coshinversefn(double x)
{
    double k;
    k=lnfunction(x+powerfunction(x*x-1,0.5));
    return k;
}

//Inverse of tangenth function//
double tanhinversefunction(double x)
{
    double k;
    k=0.5*lnfunction(1+x/1-x);
    return k;
}

//Function for solving quadratic equation//
double quadratic_eq_soln(double a,double b,double c)
{
    double x,t,g,h;
    if(a==0&&b==0)
    {
        return NULL;
    }
    if (a==0&&b!=0) return -c/b;
    t=(b*b-4*a*c);
    if(a!=0&&t!=0)
    {
        g=(-b+powerfunction(t,0.5))/(2*a);
        h=(-b-powerfunction(t,0.5))/(2*a);
        cout<<"Solutions are"<<g<<"and"<<h;
    }
    if(a!=0&&t==0)
    {
        return -b/(2*a);
    }
}

//Function involving some triangle properties//
double trianglefn(double x1,double y1,double x2,double y2,double x3,double y3)
{
    cout << "For a triangle ABC";
    cout<<"Enter the coordinates of the point A"<<endl;
    cout << "x1 = ";
    cin>> x1;
    cout << "y1 = ";
    cin>>y1;
    cout<<"Enter the coordinates of the second point B"<<endl;
    cout << "x2 = ";
    cin>> x2;
    cout << "y2 = ";
    cin>>y2;
    cout<<"Enter the coordinates of the third point C"<<endl;
    cout << "x3 = ";
    cin>> x3;
    cout << "y3 = ";
    cin>>y3;
    int x;
    cout<<"Give command"<<endl;
    cout<<"1.sides"<<endl;
    cout<<"2.centroid"<<endl;
    cout<<"3.incentre"<<endl;
    cout<<"4.orthocentre"<<endl;
    cout<<"5.circumcentre"<<endl;
    cout<<"6.area"<<endl;
    cout<<"7.inradius"<<endl;
    cout<<"8.circumradius"<<endl;
    cout<<"9.lengths of medians"<<endl;
    cout<<"10.angles"<<endl;
    cin>>x;
    float a,b,c,g1,g2,i1,i2,s,k,t1,t2,t3,o1,o2,A,B,C,r,R,s1,s2,m1,m2,m3;
    a=powerfunction((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2),0.5);
    b=powerfunction((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3),0.5);
    c=powerfunction((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2),0.5);
    g1=(x1+x2+x3)/3;
    g2=(y1+y2+y3)/3;
    i1=((a*x1)+(b*x2)+(c*x3))/(a+b+c);
    i2=((a*y1)+(b*y2)+(c*y3))/(a+b+c);
    s=(a+b+c)/2;
    k=powerfunction((s-a)*(s-b)*(s-c)*s,0.5);
    R=(a*b*c)/(4*k);
    A=cosineinverse(((b*b)+(c*c)-(a*a))/(2*b*c));
    B=cosineinverse(((a*a)+(c*c)-(b*b))/(2*a*c));
    C=cosineinverse(((b*b)+(a*a)-(c*c))/(2*b*a));
    t1=tangentfunction(A);
    t2=tangentfunction(B);
    t3=tangentfunction(C);
    if(A==90)
    {
        o1=x1;
        o2=y1;
    }
    if(B==90)
    {
        o1=x2;
        o2=y2;
    }
    if(C==90)
    {
        o1=x3;
        o2=y3;
    }
    if (A!=90&&B!=90&C!=90)
    {
        o1=((t1*x1)+(t2*x2)+(t3*x3))/(t1+t2+t3);
        o2=((t1*y1)+(t2*y2)+(t3*y3))/(t1+t2+t3);
    }
    r=k/s;
    s1=((3*g1)-o1)/2;
    s2=((3*g2)-o2)/2;
    m1=powerfunction((b*b)/2+(c*c)/2-(a*a)/4,0.5);
    m2=powerfunction((a*a)/2+(c*c)/2-(b*b)/4,0.5);
    m3=powerfunction((b*b)/2+(a*a)/2-(c*c)/4,0.5);
    switch(x)
    {
    case 1:
        cout<<"BC = "<<a<<" CA = "<<b<<" AB = "<<c;
        break;
    case 2:
        cout<<g1<<" "<<g2;
        break;
    case 3:
        cout<<i1<<" "<<i2;
        break;
    case 4:
        cout<<o1<<" "<<o2;
        break;
    case 5:
        cout<<s1<<" "<<s2;
        break;
    case 6:
        cout<<k;
        break;
    case 7:
        cout<<r;
        break;
    case 8:
        cout<<R;
        break;
    case 9:
        cout<<m1<<" "<<m2<<" "<<m3;
        break;
    case 10:
        cout<<A<<" "<<B<<" "<<C;
        break;
    default:
        cout<<"error";
    }
}

//Function for solving cubic equation//
void cubicrootsequation(float a,float b,float c,float d)
{
    float t,e,f;
    t=(18*a*b*c*d)-(4*b*b*b*d)+(b*b*c*c)-(4*a*c*c*c)-(27*a*a*d*d);
    if(t>0)
    {
        cout<<"Eqn has 3 distinct real roots";
    }
    else if(t=0)cout<<"Eqn has multiple roots and all are real";
    else
    {
        cout<<"1 real root, 2 non real complex conjugate roots";
    }
}

//differetiation function//
void differentialfn()
{
    int x;
    cout<<"Give choice"<<endl;
    cout<<"1.ae^(mx)+b^e(-nx) "<<endl;
    cout<<"2.aSinmx+bCosnx "<<endl;
    cout<<"3.nth degree polynomial "<<endl;
    cout<<"4.alognx "<<endl;
    cout<<"5.aTan(mx)+bSec(nx) "<<endl;
    cout<<"6.c"<<endl;
    cin>>x;
    switch(x)
    {
        double a,m,b,n,c;
    case 1:
        cout<<"Give a,m,b,n";
        cin>>a>>m>>b>>n;
        cout<<a*m<<"e^("<<m<<"x)+"<<-b*n<<"e^(-"<<n<<"x)"<<endl;
        break;
    case 2:
        cout<<"Give a,m,b,n"<<endl;
        cin>>a>>m>>b>>n;
        cout<<a*m<<"Cos"<<m<<"x+"<<-b*n<<"Sin"<<n<<"x"<<endl;
        break;
    case 3:
        int l;
        float k;
        float A[10];
        cout<<"Give the degree of polynomial";
        cin>>l;
        if(l<0||l>10)cout<<"Bad input";
        for(int i=0; i<10; i++)
        {
            A[i]=0;
        }
        cout<<"In polynomial a0+a1x+a2x^2+..... give the values of a0,a1,a2...an"<<endl;
        for(int j=0; j<l+1; j++)
        {
            cin>>A[j];
        }
        for(int j=1; j<l+1; j++)
        {
            if(j<l)cout<<A[j]*j<<"x^"<<j-1<<"+";
            else
            {
                cout<<A[j]*j<<"x^"<<j-1;
            }
        }
        break;
    case 4:
        cout<<"Give a,n";
        cin>>a>>n;
        cout<<a*n<<"/x";
        break;
    case 5:
        cout<<"Give a,m,b,n";
        cin>>a>>m>>b>>n;
        cout<<a*m<<"(Secmx)^2+"<<b*n<<"Secmx.Tannx";
        break;
    case 6:
        cout<<"Give c";
        cin>>c;
        cout<<"0";
        break;
    default :
        cout<<"Error";
        return;
    }
}

//Function for solving simultaneous linear equation//
void lineareqn()
{
    int i, j, k, n;
    float MatA[100][100], MatB[100], X[100];
    float Divisor, Factor, sum;
    cout<<"Enter the number of equations involved"<<endl;
    cin >> n;
//reading matrix A
    cout<<"Give the coefficient matrix of the simultaneous linear equations"<<endl;
    for(i=0; i< n; i++)
    {
        for(j=0; j < n; j++)
        {

            cin >> MatA[i][j];
        }
    }
//reading matrix B
    cout<<"Guve the constant matrix"<<endl;
    for(i=0; i< n; i++)
    {
        cin >> MatB
            [i];
    }
//Gauss elimination
    for (i=0; i< n; i++)
    {
        Divisor = MatA[i][i];
        MatA[i][i] = 1.0;
// divide all values in the row by the divisor
// to recalculate all coefficients in that row
        for (j = i+1; j < n; j++)
        {
            MatA[i][j] = MatA[i][j]/Divisor;
        }
//Also divide the corresponding RHS element
        MatB[i] = MatB[i]/Divisor;
// now replace subsequent rows, by subtracting the
// appropriate portion of the ith equation from it
        if (i+1 < n)
        {
            for (k=i+1; k<n; k++)
            {
                Factor = MatA[k][i];
                MatA[k][i] = 0.0;
                for (j = i+1; j < n; j++)
                {
                    MatA[k][j] = MatA[k][j] - Factor * MatA[i][j];
                }
                MatB[k] = MatB[k] - Factor * MatB[i];
            }
        }
    }
// back substitution starting with last variable
    X[n-1] = MatB[n-1];
    for (i = n-2; i>=0; i--)
    {
// Sum up ith row using values of X already determined
        sum = 0.0;
        for (j = i+1; j < n; j++)
        {
            sum = sum + MatA[i][j] * X[j];
        }
        X[i] = MatB[i] - sum;
    }
//output the results
    for(i=0; i< n; i++)
    {
        for (j = 0; j < n; j++)
        {
            cout << MatA
                 [i][j] << " ";
        }
        cout << " " << MatB
             [i] << endl;
    }
    for (i=0; i<n; i++)
    {
        cout << "X[" << i << "] is: " ;
        cout << X[i] << endl;
    }
    return;
}

//funtion for finding prime numbers below a certain number//
void primeno()
{
    int n,m=0;//m=no. of prime//
    cout<<"Give the value of n below which you want to find primes"<<endl;
    cin>>n;
    int removed[1000];
    if(n<0||n==0||n>1000)
    {
        cout<<"Bad input";
    }
    for (int i=0; i<1000; i++)
    {
        removed[i]=0;
    }
    for (int j=2; j<n; j++)
    {
        for(int k=j+1; k<n; k++)
        {
            if(k%j==0)
            {
                removed[k]=-1;
            }
        }
    }
    for(int i=2; i<n; i++)
    {
        if(removed[i]==0)
        {
            m++;
            cout<<i<<" ";
        }
    }
    cout<<endl<<"no. of prime before n is: "<<m<<endl;
    return;
}

//Function for doing some computation with complex number//
void complexno()
{
    int x;
    float a,b,c,d;
    cout<<"Give the choice:"<<endl;
    cout<<"1.Addition of complex numbers "<<endl;
    cout<<"2.Subtraction of complex numbers "<<endl;
    cout<<"3.Multiplication of complex numbers "<<endl;
    cout<<"4.Division of complex numbers "<<endl;
    cout<<"5.Square root of of complex number"<<endl;
    cin>>x;
    switch(x)
    {
    case 1:
        cout<<"Give the first complex number a+ib.Give a,b"<<endl;
        cin>>a>>b;
        cout<<"Give the second complex number c+id.Give c,d"<<endl;
        cin>>c>>d;
        cout<<"The sum of complex numbers is "<<a+c<<"+i("<<b+d<<")";
        break;
    case 2:
        cout<<"Give the first complex number a+ib.Give a,b"<<endl;
        cin>>a>>b;
        cout<<"Give the second complex number c+id.Give c,d"<<endl;
        cin>>c>>d;
        cout<<"The difference of complex numbers is "<<a-c<<"+i("<<b-d<<")";
        break;
    case 3:
        cout<<"Give the first complex number a+ib.Give a,b"<<endl;
        cin>>a>>b;
        cout<<"Give the second complex number c+id.Give c,d"<<endl;
        cin>>c>>d;
        cout<<"The product of complex numbers is "<<(a*c-b*d)<<"+i("<<(c*b+d*a)<<")";
        break;
    case 4:
        cout<<"Give the first complex number a+ib.Give a,b"<<endl;
        cin>>a>>b;
        cout<<"Give the second complex number c+id.Give c,d"<<endl;
        cin>>c>>d;
        cout<<"The division of complex numbers gives "<<(a*c+b*d)/(c*c+d*d)<<"+i("<<(c*b-d*a)/(c*c+d*d)<<")";
        break;
    case 5:
        cout<<"Give the first complex number a+ib.Give a,b"<<endl;
        cin>>a>>b;
        cout<<"+/-("<<powerfunction((a+powerfunction(a*a+b*b,0.5))/2,0.5)<<"+i"<<(b/modulusfn(b,0))*powerfunction((-a+powerfunction(a*a+b*b,0.5))/2,0.5)<<")";
        break;
    default:
        cout<<"What???";
        break;
    }
    return;

}

//Integral function//
void integralfn()
{
    int x;
    float a,b,m,n,c;
    cout<<"Give command "<<endl;
    cout<<"1.ae^(mx)+be^(-nx) "<<endl;
    cout<<"2.asinmx+bcosnx "<<endl;
    cout<<"3.a/x 4.c 5.a0+a1x+a2(x^2)+..."<<endl;
    cin>>x;
    switch(x)
    {
    case 1:
        cout<<"Give the values of a,m,b,n"<<endl;
        cin>>a>>m>>b>>n;
        cout<<a/m<<"e^("<<m<<"x)+("<<-b/n<<")e^(-"<<n<<")x";
        break;
    case 2:
        cout<<"Give the values of a,m,b,n"<<endl;
        cin>>a>>m>>b>>n;
        cout<<-a/m<<"cos("<<m<<"x)+("<<b/n<<")sin("<<n<<"x)";
        break;
    case 3:
        cout<<"Give a"<<endl;
        cin>>a;
        cout<<a<<"logx";
        break;
    case 4:
        cout<<"Give c";
        cin>>c;
        cout<<c<<"x";
        break;
    case 5:
        int l;
        float k;
        float A[10];
        cout<<"Give the degree of polynomial";
        cin>>l;
        if(l<=0||l>10)
        {
            cout<<"Bad input";
            break;
        }
        for(int i=0; i<10; i++)
        {
            A[i]=0;
        }
        cout<<"In polynomial a0+a1x+a2x^2+..... give the values of a0,a1,a2...an"<<endl;
        for(int j=0; j<l+1; j++)
        {
            cin>>A[j];
        }
        for(int j=0; j<l+1; j++)
        {
            if(j<l)cout<<A[j]/(j+1)<<"x^"<<j+1<<"+";
            else
            {
                cout<<A[j]/(j+1)<<"x^"<<j+1;
            }
        }
        break;
    }
    return;
}

//Function involving some operation related to vectors//
void vectorfn()
{
    int x;
    float v1x,v1y,v1z,v2x,v2y,v2z,vx,vy,vz;
    cout<<"1.Addition "<<endl;
    cout<<"2.subtraction "<<endl;
    cout<<"3.scalar product "<<endl;
    cout<<"4.vector product "<<endl;
    cout<<"5.modulus "<<endl;
    cout<<"6.angles"<<endl;
    cin>>x;
    switch(x)
    {
    case 1:
        cout<<"Give components of first vector"<<endl;
        cout<<"vx = ";
        cin>>v1x;
        cout<<"vy = ";
        cin>>v1y;
        cout<<"vz = ";
        cin>>v1z;
        cout<<"Give components of second vector"<<endl;
        cout<<"vx = ";
        cin>>v2x;
        cout<<"vy = ";
        cin>>v2y;
        cout<<"vz = ";
        cin>>v2z;
        cout<<v1x+v2x<<"i + "<<v1y+v2y<<"j + "<<v1z+v2z<<"k"<<endl;
        break;
    case 2:
        cout<<"Give components of first vector"<<endl;
        cout<<"vx = ";
        cin>>v1x;
        cout<<"vy = ";
        cin>>v1y;
        cout<<"vz = ";
        cin>>v1z;
        cout<<"Give components of second vector"<<endl;
        cout<<"vx = ";
        cin>>v2x;
        cout<<"vy = ";
        cin>>v2y;
        cout<<"vz = ";
        cin>>v2z;
        cout<<v1x-v2x<<"i + "<<v1y-v2y<<"j + "<<v1z-v2z<<"k"<<endl;
        break;
    case 3:
        cout<<"Give components of first vector"<<endl;
        cout<<"vx = ";
        cin>>v1x;
        cout<<"vy = ";
        cin>>v1y;
        cout<<"vz = ";
        cin>>v1z;
        cout<<"Give components of second vector"<<endl;
        cout<<"vx = ";
        cin>>v2x;
        cout<<"vy = ";
        cin>>v2y;
        cout<<"vz = ";
        cin>>v2z;
        cout<<(v1x*v2x+v1y*v2y+v1z*v2z)<<endl;
        break;
    case 4:
        cout<<"Give components of first vector"<<endl;
        cout<<"vx = ";
        cin>>v1x;
        cout<<"vy = ";
        cin>>v1y;
        cout<<"vz = ";
        cin>>v1z;
        cout<<"Give components of second vector"<<endl;
        cout<<"vx = ";
        cin>>v2x;
        cout<<"vy = ";
        cin>>v2y;
        cout<<"vz = ";
        cin>>v2z;
        cout<<v1z*v2y-v2y*v1z<<"i + "<<v1z*v2x-v2z*v1x<<"j + "<<v1x*v2y-v2x*v1y<<"k"<<endl;
        break;
    case 5:
        cout<<"Give components of the vector"<<endl;
        cout<<"vx";
        cin>>vx;
        cout<<"vy";
        cin>>vy;
        cout<<"vz";
        cin>>vz;
        cout<<powerfunction(vx*vx+vy*vy+vz*vz,0.5)<<endl;
        break;
    case 6:
        cout<<"Give components of the vector"<<endl;
        cout<<"vx";
        cin>>vx;
        cout<<"vy";
        cin>>vy;
        cout<<"vz";
        cin>>vz;
        cout<<cosineinverse(vx/powerfunction(vx*vx+vy*vy+vz*vz,0.5))<<" "<<cosineinverse(vy/powerfunction(vx*vx+vy*vy+vz*vz,0.5))<<" "<<cosineinverse(vz/powerfunction(vx*vx+vy*vy+vz*vz,0.5));
        break;
    default:
        cout<<"wrong input";
    }
    return;
}

//deteminant function//
double det(int n, double mat[10][10])

{
    float d;
    int c, subi, i, j, subj;
    double submat[10][10];
    if (n == 2)
    {
        return( (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    }
    else
    {
        for(c = 0; c < n; c++)
        {
            subi = 0;
            for(i = 1; i < n; i++)
            {
                subj = 0;
                for(j = 0; j < n; j++)
                {
                    if (j == c)
                    {
                        continue;
                    }
                    submat[subi][subj] = mat[i][j];
                    subj++;
                }
                subi++;
            }
            d = d + (powerfunction(-1 ,c) * mat[0][c] * det(n - 1 ,submat));
        }
    }
    return d;
}

//Function for some matrix operation//
void matrixfn()
{
    int x;
    int m, n, c, d,first[20][20],second[20][20],sum[20][20],difference[20][20];
    cout<<"Give choice "<<endl;
    cout<<"1.Addition "<<endl;
    cout<<"2.Subtraction "<<endl;
    cout<<"3.Multiplication"<<endl;
    cin>>x;
    switch(x)
    {
    case 1:

        cout << "Enter the number of rows and columns of matrix ";
        cin >> m >> n;
        cout << "Enter the elements of first matrix\n";

        for (  c = 0 ; c < m ; c++ )
        {
            for ( d = 0 ; d < n ; d++ )
                cin >> first[c][d];
            cout<<endl;
        }
        cout<<endl;

        cout << "Enter the elements of second matrix\n";

        for ( c = 0 ; c < m ; c++ )
        {
            for ( d = 0 ; d < n ; d++ )
                cin >> second[c][d];
            cout<<endl;
        }
        cout<<endl;

        for ( c = 0 ; c < m ; c++ )
            for ( d = 0 ; d < n ; d++ )
                sum[c][d] = first[c][d] + second[c][d];

        cout << "Sum of entered matrices:-"<<endl;

        for ( c = 0 ; c < m ; c++ )
        {
            for ( d = 0 ; d < n ; d++ )
                cout << sum[c][d] << "\t";

            cout << endl;
        }
        break;
    case 2:

        cout << "Enter the number of rows and columns of matrix ";
        cin >> m >> n;
        cout << "Enter the elements of first matrix\n";

        for (  c = 0 ; c < m ; c++ )
            for ( d = 0 ; d < n ; d++ )
                cin >> first[c][d];

        cout << "Enter the elements of second matrix\n";

        for ( c = 0 ; c < m ; c++ )
            for ( d = 0 ; d < n ; d++ )
                cin >> second[c][d];

        for ( c = 0 ; c < m ; c++ )
            for ( d = 0 ; d < n ; d++ )
                difference[c][d] = first[c][d] - second[c][d];

        cout << "Difference of entered matrices:-"<<endl;

        for ( c = 0 ; c < m ; c++ )
        {
            for ( d = 0 ; d < n ; d++ )
                cout << difference[c][d] << "\t";

            cout << endl;
        }
        break;
    case 3:
        int a[5][5],b[5][5],c[5][5],m,n,p,q,i,j,k;
        cout<<"Enter rows and columns of first matrix:";
        cin>>m>>n;
        cout<<"Enter rows and columns of second matrix:";
        cin>>p>>q;

        if(n==p)
        {
            cout<<"\nEnter first matrix:\n";
            for(i=0; i<m; ++i)
                for(j=0; j<n; ++j)
                    cin>>a[i][j];

            cout<<"\nEnter second matrix:\n";
            for(i=0; i<p; ++i)
                for(j=0; j<q; ++j)
                    cin>>b[i][j];
            cout<<"\nThe new matrix is:\n";

            for(i=0; i<m; ++i)
            {
                for(j=0; j<q; ++j)
                {
                    c[i][j]=0;
                    for(k=0; k<n; ++k)
                        c[i][j]=c[i][j]+(a[i][k]*b[k][j]);
                    cout<<c[i][j]<<" ";
                }
                cout<<"\n";
            }
        }
        else
            cout<<"\nSorry!!!! Matrix multiplication can't be done";
        getch();
        break;



    }
    return;
}

long int fibo(int n)
{
    int e=0,f=1,c=0;
    cout<<e<<" "<<f<<" ";
    for(int i=1; i<=n-2; i++)
    {
        c=e+f;
        e=f;
        f=c;
        cout<<c<<" ";
    }

}

long int palindrome(int n)
{

    int temp = n, reverse=0;

    while( temp != 0 )
    {
        reverse = reverse * 10;
        reverse = reverse + temp%10;
        temp = temp/10;
    }

    if ( n == reverse )
        cout<<n<<"is a Palindrome Number";
    else
        cout<<n<<"is not a Palindrome Number";

}
long int armstrong(int n)
{
    int s,m,p,r;
    s=0;
    m=n;
    do
    {
        r=n%10;
        n=n/10;
        s=s+r*r*r;
    }
    while (n!=0);
    if(s==m)
        cout<<"This is Armstrong number"<<endl;
    else
        cout<<"This is not Armstrong number"<<endl;

}

void Area()
{
    double r,sh,h,l,b,ar;
    int ch;
    cout<<"\n\n\t#**********************************************#\n";
    cout<<"\t#                                              #\n";
    cout<<"\t#                   input                      #\n";
    cout<<"\t#   1-for area of triangle                     #\n";
    cout<<"\t#   2-for area of square                       #\n";
    cout<<"\t#   3-for area of rectangle                    #\n";
    cout<<"\t#   4-for area of circle                       #\n";
    cout<<"\t#   5-curve surface area of cone               #\n";
    cout<<"\t#   6-curve surface area of hemisphere         #\n";
    cout<<"\t#   7-curve surface area of cylinder           #\n";
    cout<<"\t#   8-surface area of cube                     #\n";
    cout<<"\t#   9-surface area of cuboid                   #\n";
    cout<<"\t#   10-surface area of sphere                  #\n";
    cout<<"\t#   11-total surface area of cone              #\n";
    cout<<"\t#   12-total surface area ofhemisphere         #\n";
    cout<<"\t#   13-total surface area ofcylinder           #\n";
    cout<<"\t#                                              #\n";
    cout<<"\t#**********************************************#\n\n\t";
    cin>>ch;
    switch(ch)
    {
    case 1:
        cout<<"\n\tenter sides of triangle  ";
        cin>>l>>b>>h;
        ar=sqrt(((l+b+h)/2)*((l+b-h)/2)*((l-b+h)/2)*((b+h-l)/2));
        cout<<"\n\trequired area is  "<<ar;
        break;
    case 2:
        cout<<"\n\tenter the side of square ";
        cin>>l;
        ar=l*l;
        cout<<"\n\trequired area is  "<<ar;
        break;
    case 3:
        cout<<"\n\tenter length of rectangle ";
        cin>>l;
        cout<<"\n\tenter the width           ";
        cin>>b;
        ar=l*b;
        cout<<"\n\trequired area is  "<<ar;
        break;
    case 4:
        cout<<"\n\tenter radius of circle\t";
        cin>>r;
        ar=3.14*r*r;
        cout<<"\n\trequired area is  "<<ar;
        break;
    case 5:
        cout<<"\n\tenter the radius of cone   ";
        cin>>r;
        cout<<"\n\tenter the slant height     ";
        cin>>sh;
        ar=3.14*r*sh;
        cout<<"\n\trequired area is           "<<ar;
        break;
    case 6:
        cout<<"\n\tenter the radius of hemisphere   ";
        cin>>r;
        ar=3.14*r*r*2;
        cout<<"\n\trequired area is    "<<ar;
        break;
    case 7:
        cout<<"\n\tenter the radius of cylinder  ";
        cin>>r;
        cout<<"\n\tenter the height of cylinder  ";
        cin>>h;
        ar=2*3.14*r*h;
        cout<<"\n\trequired area is         "<<ar;
        break;
    case 8:
        cout<<"\n\tenter the side of cube     ";
        cin>>l;
        ar=6*l*l;
        cout<<"\n\trequired area is           "<<ar;
        break;
    case 9:
        cout<<"\n\tenter length      ";
        cin>>l;
        cout<<"\n\tenter writh       ";
        cin>>b;
        cout<<"\n\tenter height      ";
        cin>>h;
        ar=2*(l*b+b*h+l*h);
        cout<<"\n\trequired area is  "<<ar;
        break;
    case 10:
        cout<<"\n\tenter the radius   ";
        cin>>r;
        ar=4*3.14*r*r;
        cout<<"\n\trequired area is   "<<ar;
        break;
    case 11:
        cout<<"\n\tenter radius of cone   ";
        cin>>r;
        cout<<"\n\tenter slant height     ";
        cin>>sh;
        ar=3.14*r*(sh+r);
        cout<<"\n\trequired area is       "<<ar;
        break;
    case 12:
        cout<<"\n\tenter the radius of hemispher   ";
        cin>>r;
        ar=3*3.14*r*r;
        cout<<"\n\trequired area is     "<<ar;
        break;
    case 13:
        cout<<"\n\tenter the radius of cylinder     ";
        cin>>r;
        cout<<"\n\tenter the height of cylinder     ";
        cin>>h;
        ar=2*3.14*r*(l+r);
        cout<<"\n\trequired area is                 "<<ar;
        break;
    default :
        cout<<"\n\twrong no\n";
    }
}

void volume()
{
    double l,b,h,r,vol;
    int ch;
    cout<<"\n\n      #*************************************************#\n";
    cout<<"      #                                                 #\n";
    cout<<"      #                       enter                     #\n";
    cout<<"      #      1-for volume of cube                       #\n";
    cout<<"      #      2-for volume of cuboid                     #\n";
    cout<<"      #      3-for cylinder                             #\n";
    cout<<"      #      4-for volume of hemisphere                 #\n";
    cout<<"      #      5-for volume of sphere                     #\n";
    cout<<"      #      6-for volume of cone                       #\n";
    cout<<"      #                                                 #\n";
    cout<<"      #*************************************************#\n\t\t";
    cin>>ch;
    switch(ch)
    {
    case 1:
        cout<<"\n\tenter side of cube   ";
        cin>>l;
        vol=l*l*l;
        cout<<"\n\trequired volume is     "<<vol;
        break;
    case 2:
        cout<<"\n\tenter the length of cuboid    ";
        cin>>l;
        cout<<"\n\tenter the wridth of cuboid    ";
        cin>>b;
        cout<<"\n\tenter the height of cuboid    ";
        cin>>h;
        vol=l*h*b;
        cout<<"\n\trequired volume is            "<<vol;
        break;
    case 3:
        cout<<"\n\tenter the radius of cylinder       \n";
        cin>>r;
        cout<<"\n\tenter the height of the cylinder   \n";
        cin>>h;
        vol=3.14*r*r*h;
        cout<<"\n\trequired volume is    "<<vol;
        break;
    case 4:
        cout<<"\n\tenter the radius of hemisphere       ";
        cin>>r;
        vol=2*3.14*r*r*r/3;
        cout<<"\n\trequired volume is\t"<<vol;
        break;
    case 5:
        cout<<"\n\tenter the radius of sphere     ";
        cin>>r;
        vol=4*3.14*r*r*r/3;
        cout<<"\n\trequired volume is\t"<<vol;
        break;
    case 6:
        cout<<"\n\tenter the radius of the cone   ";
        cin>>r;
        cout<<"\n\tenter the height of the cone   ";
        cin>>h;
        vol=3.14*r*r*h/3;
        cout<<"\n\trequired volume is\t"<<vol;
        break;
    default :
        cout<<"\n\twrong entry\n\t";

    }
}

double simpleinterest()
{
    double pr, ra, ti, answ;
    cout<<"\n\t enter principle     ";
    cin>>pr;
    cout<<"\n\t enter rate          ";
    cin>>ra;
    cout<<"\n\t enter time in year  ";
    cin>>ti;
    answ=pr*ra*ti/100;
    cout<<"\n\tresult is"<<answ;
}

double compoundinterest()
{
    double pr, ra, ti, answ, res;
    cout<<"\n\t enter principle     ";
    cin>>pr;
    cout<<"\n\t enter rate          ";
    cin>>ra;
    cout<<"\n\t enter time in year  ";
    cin>>ti;
    res=1+ra/100;
    answ=pr*pow(res,ti);
    answ=answ-pr;
    cout<<"\n\tresult is"<<answ;
}

main_program


{

    initCanvas("Calculator",600,700);

    wait(1);

//Construction of rectangles on canvas//

    Rectangle button1(50,50,50,50);
    Text t1(50,50,"Add(+)");
    Rectangle button2(150,50,50,50);
    Text t2(150,50,"Sub(-)");
    Rectangle button3(250,50,50,50);
    Text t3(250,50,"Multi(X)");
    Rectangle button4(350,50,50,50);
    Text t4(350,50,"Div(/)");
    Rectangle button5(450,50,50,50);
    Text t5(450,50,"n!");
    Rectangle button6(50,150,50,50);
    Text t6(50,150,"e^x");
    Rectangle button7(150,150,50,50);
    Text t7(150,150,"log");
    Rectangle button8(250,150,50,50);
    Text t8(250,150,"ln");
    Rectangle button9(350,150,50,50);
    Text t9(350,150,"Rad");
    Rectangle button10(450,150,50,50);
    Text t10(450,150,"x^n");

    Rectangle button11(50,250,50,50);
    Text t11(50,250,"sin");
    Rectangle button12(150,250,50,50);
    Text t12(150,250,"cos");
    Rectangle button13(250,250,50,50);
    Text t13(250,250,"tan");
    Rectangle button14(350,225,50,25);
    Text t14(350,225,"nPr");
    Rectangle button36(350,275,50,25);
    Text t36(350,275,"nCr");
    Rectangle button15(450,250,50,50);
    Text t15(450,250,"Matrix");
    Rectangle button16(50,350,50,50);
    Text t16(50,350,"arcsin");
    Rectangle button17(150,350,50,50);
    Text t17(150,350,"arccos");
    Rectangle button18(250,350,50,50);
    Text t18(250,350,"arctan");
    Rectangle button19(350,350,50,50);
    Text t19(350,350,"Inverse");
    Rectangle button20(450,350,50,50);
    Text t20(450,350,"|x|");

    Rectangle button21(50,450,50,50);
    Text t21(50,450,"sinh");
    Rectangle button22(150,450,50,50);
    Text t22(150,450,"cosh");
    Rectangle button23(250,450,50,50);
    Text t23(250,450,"tanh");
    Rectangle button24(350,450,50,50);
    Text t24(350,450,"Quad Eqn");
    Rectangle button25(450,450,50,50);
    Text t25(450,450,"Vectors");
    Rectangle button26(50,550,50,50);
    Text t26(50,550,"arcsinh");
    Rectangle button27(150,550,50,50);
    Text t27(150,550,"arccosh");
    Rectangle button28(250,550,50,50);
    Text t28(250,550,"arctanh");
    Rectangle button29(350,550,50,50);
    Text t29(350,550,"Cubic Eqn");
    Rectangle button30(450,550,50,50);
    Text t30(450,550,"Triangle Prop");

    Rectangle button31(50,650,50,50);
    Text t31(50,650,"Diffn");
    Rectangle button32(150,650,50,50);
    Text t32(150,650,"Integral");
    Rectangle button33(250,650,50,50);
    Text t33(250,650,"Prime");
    Rectangle button34(350,650,50,50);
    Text t34(350,650,"Linear Eqn");
    Rectangle button35(450,650,50,50);
    Text t35(450,650,"Complex No");

    Rectangle button39(550,50,50,50);
    Text t39(550,50,"fibonacci No");
    Rectangle button37(550,150,50,50);
    Text t37(550,150,"palindrome");
    Rectangle button38(550,250,50,50);
    Text t38(550,250,"armstrong");
    Rectangle button40(550,350,50,50);
    Text t40(550,350,"AREA");
    Rectangle button41(550,650,50,50);
    Text t41(550,650,"VOLUME");
    Rectangle button42(550,550,50,50);
    Text t42(550,550,"Simp int");
    Rectangle button43(550,450,50,50);
    Text t43(550,450,"Comp int");

//contruction completed//
    int con=1;

    while(con==1)
    {
//gives the coordinate of the point where the mouse clicks//

        int clickPos = getClick();

        int cx = clickPos/65536;
        int cy = clickPos % 65536;


        //these are 36 conditions//
        //related to coordinate of the point which is clicked//

        if((cx>=25)&&(cx<=75)&&(cy>=25&&(cy<=75)))
        {
            button1.setFill(true);
            button1.setColor(COLOR("red"));
            wait(0.5);
            button1.setColor(COLOR("white"));
            button1.setFill(false);
            button1.setColor(COLOR("black"));
            wait(0.2);

            add();
        }

        if((cx>=125)&&(cx<=175)&&(cy>=25)&&(cy<=75))
        {
            button2.setFill(true);
            button2.setColor(COLOR("red"));
            wait(0.5);
            button2.setColor(COLOR("white"));
            button2.setFill(false);
            button2.setColor(COLOR("black"));
            wait(0.2);

            sub();

        }

        if((cx>=225)&&(cx<=275)&&(cy>=25&&(cy<=75)))
        {
            button3.setFill(true);
            button3.setColor(COLOR("red"));
            wait(0.5);
            button3.setColor(COLOR("white"));
            button3.setFill(false);
            button3.setColor(COLOR("black"));
            wait(0.2);

            multi();

        }


        if((cx>=325)&&(cx<=375)&&(cy>=25&&(cy<=75)))
        {
            button4.setFill(true);
            button4.setColor(COLOR("red"));
            wait(0.5);
            button4.setColor(COLOR("white"));
            button4.setFill(false);
            button4.setColor(COLOR("black"));
            wait(0.2);

            div();

        }
        double n=0;
        if((cx>=425)&&(cx<=475)&&(cy>=25&&(cy<=75)))
        {
            button5.setFill(true);
            button5.setColor(COLOR("red"));
            wait(0.5);
            button5.setColor(COLOR("white"));
            button5.setFill(false);
            button5.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of n"<<endl;
            cin>>n;

            cout<<"n! = "<<factorial(n);

        }

        if((cx>=25)&&(cx<=75)&&(cy>=125&&(cy<=175)))
        {
            button6.setFill(true);
            button6.setColor(COLOR("red"));
            wait(0.5);
            button6.setColor(COLOR("white"));
            button6.setFill(false);
            button6.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"e^x = "<<exponentialfn(n);

        }



        if((cx>=125)&&(cx<=175)&&(cy>=125)&&(cy<=175))
        {
            button7.setFill(true);
            button7.setColor(COLOR("red"));
            wait(0.5);
            button7.setColor(COLOR("white"));
            button7.setFill(false);
            button7.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"log(x)"<<logarithmicfn(n);

        }

        if((cx>=225)&&(cx<=275)&&(cy>=125&&(cy<=175)))
        {
            button8.setFill(true);
            button8.setColor(COLOR("red"));
            wait(0.5);
            button8.setColor(COLOR("white"));
            button8.setFill(false);
            button8.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"ln(x) = "<<lnfunction(n);

        }


        if((cx>=325)&&(cx<=375)&&(cy>=125&&(cy<=175)))
        {
            button9.setFill(true);
            button9.setColor(COLOR("red"));
            wait(0.5);
            button9.setColor(COLOR("white"));
            button9.setFill(false);
            button9.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of angle in degrees"<<endl;
            cin>>n;

            cout<<"The corresponding value of angle in radians is"<<endl<<radianfn(n);



        }
        double x=0;
        if((cx>=425)&&(cx<=475)&&(cy>=125&&(cy<=175)))
        {
            button10.setFill(true);
            button10.setColor(COLOR("red"));
            wait(0.5);
            button10.setColor(COLOR("white"));
            button10.setFill(false);
            button10.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>x;
            cout<<"Enter the value of n"<<endl;
            cin>>n;

            cout<<"x^n = "<<powerfunction(x,n);

        }

        if((cx>=25)&&(cx<=75)&&(cy>=225&&(cy<=275)))
        {
            button11.setFill(true);
            button11.setColor(COLOR("red"));
            wait(0.5);
            button11.setColor(COLOR("white"));
            button11.setFill(false);
            button11.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Sine("<< n << ") = "<<sinefunction(n);

        }



        if((cx>=125)&&(cx<=175)&&(cy>=225)&&(cy<=275))
        {
            button12.setFill(true);
            button12.setColor(COLOR("red"));
            wait(0.5);
            button12.setColor(COLOR("white"));
            button12.setFill(false);
            button12.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Cos("<< n << ") = "<<cosinefunction(n);

        }

        if((cx>=225)&&(cx<=275)&&(cy>=225&&(cy<=275)))
        {
            button13.setFill(true);
            button13.setColor(COLOR("red"));
            wait(0.5);
            button13.setColor(COLOR("white"));
            button13.setFill(false);
            button13.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Tan("<< n << ") = "<<tangentfunction(n);

        }


        if((cx>=325)&&(cx<=375)&&(cy>=212.5&&(cy<=237.5)))
        {
            button14.setFill(true);
            button14.setColor(COLOR("red"));
            wait(0.5);
            button14.setColor(COLOR("white"));
            button14.setFill(false);
            button14.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of n"<<endl;
            cin>>n;
            cout<<"Enter the value of x"<<endl;
            cin>>x;

            cout<<"nPr = "<<npr(n,x);




        }

        if((cx>=425)&&(cx<=475)&&(cy>=225&&(cy<=275)))
        {
            button15.setFill(true);
            button15.setColor(COLOR("red"));
            wait(0.5);
            button15.setColor(COLOR("white"));
            button15.setFill(false);
            button15.setColor(COLOR("black"));
            wait(0.2);

            matrixfn();


        }

        if((cx>=25)&&(cx<=75)&&(cy>=325&&(cy<=375)))
        {
            button16.setFill(true);
            button16.setColor(COLOR("red"));
            wait(0.5);
            button16.setColor(COLOR("white"));
            button16.setFill(false);
            button16.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;
            if(n>1||n<-1)
            {
                cout << "Angle does not exist for this value";
            }
            else
            {

                cout<<"Answer is"<<endl<<sineinverse(n);
            }

        }



        if((cx>=125)&&(cx<=175)&&(cy>=325)&&(cy<=375))
        {
            button17.setFill(true);
            button17.setColor(COLOR("red"));
            wait(0.5);
            button17.setColor(COLOR("white"));
            button17.setFill(false);
            button17.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;
            if(n>1||n<-1)
            {
                cout << "Angle does not exist for this value";
            }
            else
            {
                cout<<"Answer is"<<endl<<cosineinverse(n);
            }
        }

        if((cx>=225)&&(cx<=275)&&(cy>=325&&(cy<=375)))
        {
            button18.setFill(true);
            button18.setColor(COLOR("red"));
            wait(0.5);
            button18.setColor(COLOR("white"));
            button18.setFill(false);
            button18.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Answer is"<<endl<<taninverse(n);

        }


        if((cx>=325)&&(cx<=375)&&(cy>=325&&(cy<=375)))
        {
            button19.setFill(true);
            button19.setColor(COLOR("red"));
            wait(0.5);
            button19.setColor(COLOR("white"));
            button19.setFill(false);
            button19.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Inverse is"<<endl<<inversefn(n);



        }

        if((cx>=425)&&(cx<=475)&&(cy>=325&&(cy<=375)))
        {
            button20.setFill(true);
            button20.setColor(COLOR("red"));
            wait(0.5);
            button20.setColor(COLOR("white"));
            button20.setFill(false);
            button20.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>x;
            cout<<"Enter the value of y"<<endl;
            cin>>n;

            cout<<"Modulus is "<<modulusfn(x,n);

        }


        if((cx>=25)&&(cx<=75)&&(cy>=425&&(cy<=475)))
        {
            button21.setFill(true);
            button21.setColor(COLOR("red"));
            wait(0.5);
            button21.setColor(COLOR("white"));
            button21.setFill(false);
            button21.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"sinh("<<n << ") = "<<sinhfunction(n);

        }



        if((cx>=125)&&(cx<=175)&&(cy>=425)&&(cy<=475))
        {
            button22.setFill(true);
            button22.setColor(COLOR("red"));
            wait(0.5);
            button22.setColor(COLOR("white"));
            button22.setFill(false);
            button22.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"cosh("<<n << ") = "<<coshfunction(n);

        }

        if((cx>=225)&&(cx<=275)&&(cy>=425&&(cy<=475)))
        {
            button23.setFill(true);
            button23.setColor(COLOR("red"));
            wait(0.5);
            button23.setColor(COLOR("white"));
            button23.setFill(false);
            button23.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"tanh("<<n << ") = "<<tanhfunction(n);

        }

        double a=0,b=0,c=0;
        if((cx>=325)&&(cx<=375)&&(cy>=425&&(cy<=475)))
        {
            button24.setFill(true);
            button24.setColor(COLOR("red"));
            wait(0.5);
            button24.setColor(COLOR("white"));
            button24.setFill(false);
            button24.setColor(COLOR("black"));
            wait(0.2);
            cout << "For a quadratic equation ax^2 + bx + c "<< endl;
            cout<<"Enter the value of a"<<endl;
            cin>>a;
            cout<<"Enter the value of b"<<endl;
            cin>>b;
            cout<<"Enter the value of c"<<endl;
            cin>>c;

            quadratic_eq_soln(a,b,c);




        }

        if((cx>=425)&&(cx<=475)&&(cy>=425&&(cy<=475)))
        {
            button25.setFill(true);
            button25.setColor(COLOR("red"));
            wait(0.5);
            button25.setColor(COLOR("white"));
            button25.setFill(false);
            button25.setColor(COLOR("black"));
            wait(0.2);

            vectorfn();

        }

        if((cx>=25)&&(cx<=75)&&(cy>=525&&(cy<=575)))
        {
            button26.setFill(true);
            button26.setColor(COLOR("red"));
            wait(0.5);
            button26.setColor(COLOR("white"));
            button26.setFill(false);
            button26.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Answer is"<<endl<<sinhinversefn(n);

        }



        if((cx>=125)&&(cx<=175)&&(cy>=525)&&(cy<=575))
        {
            button27.setFill(true);
            button27.setColor(COLOR("red"));
            wait(0.5);
            button27.setColor(COLOR("white"));
            button27.setFill(false);
            button27.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Answer is"<<endl<<coshinversefn(n);

        }

        if((cx>=225)&&(cx<=275)&&(cy>=525&&(cy<=575)))
        {
            button28.setFill(true);
            button28.setColor(COLOR("red"));
            wait(0.5);
            button28.setColor(COLOR("white"));
            button28.setFill(false);
            button28.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of x"<<endl;
            cin>>n;

            cout<<"Answer is"<<endl<<tanhinversefunction(n);

        }

        int d=0;
        if((cx>=325)&&(cx<=375)&&(cy>=525&&(cy<=575)))
        {
            button29.setFill(true);
            button29.setColor(COLOR("red"));
            wait(0.5);
            button29.setColor(COLOR("white"));
            button29.setFill(false);
            button29.setColor(COLOR("black"));
            wait(0.2);

            cout<<"For a cubic equation ax^3 + bx^2 + cx + d"<< endl;
            cout<<"Enter the value of a"<<endl;
            cin>>a;
            cout<<"Enter the value of b"<<endl;
            cin>>b;
            cout<<"Enter the value of c"<<endl;
            cin>>c;
            cout<<"Enter the value of d"<<endl;
            cin>>d;


            cubicrootsequation(a,b,c,d);



        }
        float x1,x2,x3,y1,y2,y3;
        if((cx>=425)&&(cx<=475)&&(cy>=525&&(cy<=575)))
        {
            button30.setFill(true);
            button30.setColor(COLOR("red"));
            wait(0.5);
            button30.setColor(COLOR("white"));
            button30.setFill(false);
            button30.setColor(COLOR("black"));
            wait(0.1);



            trianglefn(x1,y1,x2,y2,x3,y3);
        }


        if((cx>=25)&&(cx<=75)&&(cy>=625&&(cy<=675)))
        {
            button31.setFill(true);
            button31.setColor(COLOR("red"));
            wait(0.5);
            button31.setColor(COLOR("white"));
            button31.setFill(false);
            button31.setColor(COLOR("black"));
            wait(0.2);

            differentialfn();
        }

        if((cx>=125)&&(cx<=175)&&(cy>=625)&&(cy<=675))
        {
            button32.setFill(true);
            button32.setColor(COLOR("red"));
            wait(0.5);
            button32.setColor(COLOR("white"));
            button32.setFill(false);
            button32.setColor(COLOR("black"));
            wait(0.2);

            integralfn();

        }

        if((cx>=225)&&(cx<=275)&&(cy>=625&&(cy<=675)))
        {
            button33.setFill(true);
            button33.setColor(COLOR("red"));
            wait(0.5);
            button33.setColor(COLOR("white"));
            button33.setFill(false);
            button33.setColor(COLOR("black"));
            wait(0.2);

            primeno();

        }


        if((cx>=325)&&(cx<=375)&&(cy>=625&&(cy<=675)))
        {
            button34.setFill(true);
            button34.setColor(COLOR("red"));
            wait(0.5);
            button34.setColor(COLOR("white"));
            button34.setFill(false);
            button34.setColor(COLOR("black"));
            wait(0.2);

            lineareqn();

        }

        if((cx>=425)&&(cx<=475)&&(cy>=625&&(cy<=675)))
        {
            button35.setFill(true);
            button35.setColor(COLOR("red"));
            wait(0.5);
            button35.setColor(COLOR("white"));
            button35.setFill(false);
            button35.setColor(COLOR("black"));
            wait(0.2);

            complexno();

        }

        if((cx>=325)&&(cx<=375)&&(cy>=262.5&&(cy<=287.5)))
        {
            button6.setFill(true);
            button6.setColor(COLOR("red"));
            wait(0.5);
            button6.setColor(COLOR("white"));
            button6.setFill(false);
            button6.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of n"<<endl;
            cin>>n;
            cout<<"Enter the value of r"<<endl;
            cin>>x;

            cout<<"nCr = "<<ncr(n,x);

        }

        if((cx>=525)&&(cx<=575)&&(cy>=25)&&(cy<=75))
        {
            button39.setFill(true);
            button39.setColor(COLOR("red"));
            wait(0.5);
            button39.setColor(COLOR("white"));
            button39.setFill(false);
            button39.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of n"<<endl;
            cin>>n;
            cout<<" the number of fibonacci numbers is"<<endl<<fibo(n)<<endl;
        }

        if((cx>=525)&&(cx<=575)&&(cy>=125)&&(cy<=175))
        {
            button37.setFill(true);
            button37.setColor(COLOR("red"));
            wait(0.5);
            button37.setColor(COLOR("white"));
            button37.setFill(false);
            button37.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of n"<<endl;
            cin>>n;
            cout<<palindrome(n);
        }

        if((cx>=525)&&(cx<=575)&&(cy>=225)&&(cy<=275))
        {
            button38.setFill(true);
            button38.setColor(COLOR("red"));
            wait(0.5);
            button38.setColor(COLOR("white"));
            button38.setFill(false);
            button38.setColor(COLOR("black"));
            wait(0.2);

            cout<<"Enter the value of n"<<endl;
            cin>>n;
            cout<<armstrong(n);
        }


        if((cx>=525)&&(cx<=575)&&(cy>=325)&&(cy<=375))
        {
            button40.setFill(true);
            button40.setColor(COLOR("red"));
            wait(0.5);
            button40.setColor(COLOR("white"));
            button40.setFill(false);
            button40.setColor(COLOR("black"));
            wait(0.2);

            Area();
        }


        if((cx>=525)&&(cx<=575)&&(cy>=625)&&(cy<=675))
        {
            button41.setFill(true);
            button41.setColor(COLOR("red"));
            wait(0.5);
            button41.setColor(COLOR("white"));
            button41.setFill(false);
            button41.setColor(COLOR("black"));
            wait(0.2);

            volume();
        }


        if((cx>=525)&&(cx<=575)&&(cy>=525)&&(cy<=575))
        {
            button42.setFill(true);
            button42.setColor(COLOR("red"));
            wait(0.5);
            button42.setColor(COLOR("white"));
            button42.setFill(false);
            button42.setColor(COLOR("black"));
            wait(0.2);

            cout<<simpleinterest();
        }


        if((cx>=525)&&(cx<=575)&&(cy>=425)&&(cy<=475))
        {
            button43.setFill(true);
            button43.setColor(COLOR("red"));
            wait(0.5);
            button43.setColor(COLOR("white"));
            button43.setFill(false);
            button43.setColor(COLOR("black"));
            wait(0.2);

            cout<<compoundinterest();
        }

        cout<<endl;
        cout<<endl;
        cout<<"If you want to continue to work with somemore operations press '1'";
        cout<<endl;
        cout<<"press any other key to exit from the calculator";
        cout<<endl;
        cin>>con;
        cout<<"select an operation in the calculator";
        cout<<endl;
        if(con!=1)
        {
            closeCanvas();
            return 0;

        }

    }

}

