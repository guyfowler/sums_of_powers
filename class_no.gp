/*--------------------------------------------------------------------------------------
****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Computing the discriminants of all quadratic imaginary orders with class number \leq 100

Used in the proof of Proposition 7.3 in
"Some uniform effective results on André--Oort for sums of powers in $\mathbb{C}^n$"
Guy Fowler
30 April 2024
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
****************************************************************************************

Contents
    
    Scripts on the totient function

        is_totient(n)
        totients_le(n)
        largest_totient(n)
        largest_totient_le(n)
    
    Bounds on the discriminant

        fund_class_bds(n)
        upper_bd_discs(n)

    Scripts on discriminants

        LargestFund
        discs(n)
        bdd_discs(n, k)
        discs_cl_le(n)
        write_discs_cl_le(n)
        largest_disc_cl_le(n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
***************************************
--------------------------------------*/

/*-----------------------------
******************************
Scripts on the totient function
*******************************
------------------------------*/

/*-----------------------------------------------------------------
is_totient(n)

Input: integer n>0
Output: 1 if n=totient(k) for some k>0, 0 otherwise 

Note: we use the easy bd totient(n) \geq \sqrt{n/2}
-----------------------------------------------------------------*/

is_totient(n)={
    my(b);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    b=2*n^2;

    for(i=1, b,
    
    if(eulerphi(i)==n, return(1));
    
    );

    return(0);
}

/*-----------------------------------------------------------------
totients_le(n)

Input: integer n>0
Output: vector of all those i \in {1, ..., n} s.t. is_totient(i)==1
-----------------------------------------------------------------*/

totients_le(n)={
    my(u);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    u=[];

    for(i=1, n,

    if(is_totient(i)==1, u=concat(u, [i]));

    );

    return(u);
}

/*-----------------------------------------------------------------
largest_totient(n)

Input: integer n>0
Output: max {k : totient(k) = n} if this set is non-empty, otherwise prints "Not a totient" 
-----------------------------------------------------------------*/

largest_totient(n)={
    my(b, k);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    if(is_totient(n)==0, print("Not a totient"); return());

    b=2*n^2;

    for(i=1, b,
    
    if(eulerphi(i)==n, k=i);
    
    );

    return(k);
}

/*-----------------------------------------------------------------
largest_totient_le(n)

Input: integer n>0
Output: max {k : totient(k) \leq n} 
-----------------------------------------------------------------*/

largest_totient_le(n)={
    my(a);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    for(i=1, n,

    if(is_totient(i)!=1, next());
    
    a=max(a, largest_totient(i));

    );

    return(a);
}

/*-----------------------------
******************************
Bounds on the discriminant
*******************************
------------------------------*/

/*-----------------------------------------------------------------
fund_class_bds(n)

Input: integer n>0
Output: vector with entries [k, f], where k is a positive integer of the form floor(n/totient(l)) for some positive integer l and f is the largest positive integer with floor(n/totient(f)) \leq k
-----------------------------------------------------------------*/

fund_class_bds(n)={
    my(u, l, v, w);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    u=totients_le(n);
    l=length(u);

    v=vector(l, i, floor(n/u[i]));

    w=[];

    for(i=1, l,

    if(i<l&&v[i]==v[i+1], next());

    w=concat(w, [[v[i], largest_totient_le(u[i])]]);

    );

    return(w);
}

/*-----------------------------------------------------------------
upper_bd_discs(n)

Input: integer n \in {1, ..., 100}
Output: upper bound on abs(Delta) for all (fundamental & non-fundamental) discriminants Delta with h(Delta) \leq n
-----------------------------------------------------------------*/

upper_bd_discs(n)={
    my(v, l, b, f, h, D);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);
    if(n>length(LargestFund), print("Input error"); return);

v=fund_class_bds(n);
l=length(v);

b=0;

for(i=1, l,

f=v[i][2];
h=v[i][1];

D=0;

for(j=1, h,

D=max(D, LargestFund[j][2]);

);

b=max(b, D*f^2);

);

f=largest_totient_le(3*n);
a=4*f^2;
b=max(b, a);

if(b>2*10^10, print("Warning: bound too large for PARI qfbclassno"));

return(b);

}

/*-----------------------------
******************************
Scripts on discriminants
*******************************
------------------------------*/

/*---------------------------------------------------------------------------
Vector of length 100 with ith entry [i, max {abs(D) : D is a fundamental discriminant with h(D) = i}]

Source: Table 4 in M. Watkins, "Class numbers of imaginary quadratic fields", Mathematics of Computation, Vol 73 No 246, pp.907--938, 2004.

----------------------------------------------------------------------------*/

LargestFund={
    [[1, 163], [2, 427], [3, 907], [4, 1555],[5, 2683],[6, 3763],[7, 5923],[8, 6307], [9, 10627], [10, 13843], [11, 15667], [12, 17803], [13, 20563], [14, 30067], [15, 34483], [16, 31243], [17, 37123], [18, 48427], [19, 38707], [20, 58507], [21, 61483], [22, 85507], [23, 90787], [24, 111763], [25, 93307], [26, 103027], [27, 103387], [28, 126043], [29, 166147], [30, 134467], [31, 133387], [32, 164803], [33, 222643], [34, 189883], [35, 210907], [36, 217627], [37, 158923], [38, 289963], [39, 253507], [40, 260947], [41, 296587], [42, 280267], [43, 300787], [44, 319867], [45, 308323], [46, 462283], [47, 375523], [48, 335203], [49, 393187], [50, 389467], [51, 546067], [52, 439147], [53, 425107], [54, 532123], [55, 425083], [56, 494323], [57, 615883], [58, 586987], [59, 474307], [60, 662803], [61, 606643], [62, 647707], [63, 991027], [64, 693067], [65, 703123], [66, 958483], [67, 652723], [68, 819163], [69, 888427], [70, 811507], [71, 909547], [72, 947923], [73, 886867], [74, 951043], [75, 916507], [76, 1086187], [77, 1242763], [78, 1004347], [79, 1333963], [80, 1165483], [81, 1030723], [82, 1446547], [83, 1074907], [84, 1225387], [85, 1285747], [86, 1534723], [87, 1261747], [88, 1265587], [89, 1429387], [90, 1548523], [91, 1391083], [92, 1452067], [93, 1475203], [94, 1587763], [95, 1659067], [96, 1684027], [97, 1842523], [98, 2383747], [99, 1480627], [100, 1856563]];
}

/*-----------------------------------------------------------------
discs(n)

Input: integer n>0
Output: the list of all discriminants Delta with abs(Delta) \leq n 
-----------------------------------------------------------------*/

discs(n)={
    my(u, k, r, l, v, m, D);
    u=[];

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    k=floor(n/4);
    r=n%4;

    if(r==3, l=2*k+1, l=2*k);

    v=vector(l, i, 0);

    for(i=1, k,
    
    m=2*i;
    D=-4*i;

   v[m]=D;
   v[m-1]=D+1
    
    );

    if(v[l]==0, v[l]=1-4*(k+1));

    return(v);
}

/*--------------------------------------------------------------------------------
bdd_discs(n, k)

Input: integers n, k>0 (with k < 2*10^10)
Output: vector v of length n such that the ith entry of v is the vector of all discriminants Delta with abs(Delta) \leq k and h(Delta) = i
--------------------------------------------------------------------------------*/

bdd_discs(n, k)={
    my(u, v, l, h);
    u=vector(n, i, []);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);

    if(type(k)!="t_INT", print("Input error"); return);
    if(k<=0, print("Input error"); return);

    v=discs(k);
    l=length(v);

    for(i=1, l,
    
    Delta=v[i];
    h=qfbclassno(Delta);

    if(h>n, next(), u[h]=concat(u[h], [Delta]));
    
    );

    return(u);
}

/*--------------------------------------------------------------------------------
discs_cl_le(n)

Input: integer n \in {1, ..., 100}
Output: vector v of length n such that the ith entry of v is the vector of all discriminants Delta with h(Delta) = i
--------------------------------------------------------------------------------*/

discs_cl_le(n)={
    my(u, v, l, h, b);
    u=vector(n, i, []);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);
    if(n>length(LargestFund), print("Input error"); return);

    b=upper_bd_discs(n);

    v=discs(b);
    l=length(v);

    for(i=1, l,
    
    Delta=v[i];
    h=qfbclassno(Delta);

    if(h>n, next(), u[h]=concat(u[h], [Delta]));
    
    );

    return(u);
}


/*--------------------------------------------------------------------------------
write_discs_cl_le(n)

Input: integer n \in {1, ..., 100}
Output: creates, for each i=1,...,n, a .txt file "discs_with_cl_no_i" containing all the discriminants Delta with h(Delta)=i
--------------------------------------------------------------------------------*/
write_discs_cl_le(n)={
    my(v);

     if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);
    if(n>length(LargestFund), print("Input error"); return);

    v=discs_cl_le(n);

for(i = 1, n,
   filename = concat("discs_with_cl_no_", concat(Str(i), ".txt"));
   
    write(filename, v[i]);   
    );
}


/*--------------------------------------------------------------------------------
largest_disc_cl_le(n)

Input: integer n \in {1, ..., 100}
Output: the maximum value of abs(Delta) for Delta a discriminant with h(Delta) \leq n
--------------------------------------------------------------------------------*/

largest_disc_cl_le(n)={
     my(v, l, b);

    if(type(n)!="t_INT", print("Input error"); return);
    if(n<=0, print("Input error"); return);
    if(n>length(LargestFund), print("Input error"); return);

    v=discs_cl_le(n);
    b=0;

    for(i=1, n,
    
    l=length(v[i]);
    
    for(j=1, l, 
    
    b=max(b, abs(v[i][j]));
    );
    );

    return(b);
}