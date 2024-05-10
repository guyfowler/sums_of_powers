/*---------------------------------------------------------------------------------
***********************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PARI scripts for Section 8.3 in
"Some uniform effective results on André--Oort for sums of powers in $\mathbb{C}^n$"
Guy Fowler
30 April 2024
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
***********************************************************************************

Contents

    General bounds on singular moduli

        j_max(k, r)
        j_min(k, r)

    Bounds for Section 8.3.1

        lower_1_1_1(k)
        upper_1_1_1(k)
        bd_1_1_1()
        lower_1_1_2(k)
        upper_1_1_2(k)
        bd_1_1_2()
        lower_1_2(k)
        upper_1_2(k)
        bd_1_2()

    Bounds for Section 8.3.2

        lower_2_1(k, l)
        upper_2_1(k, l)
        bd_2_1(l)
        lower_2_2(k, l)
        upper_2_2(k, l, a_min)
        bd_2_2(l, a_min)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
******************************************************************
-----------------------------------------------------------------*/



/*-----------------------------------
***********************************
General bounds on singular moduli
**********************************
---------------------------------*/


/*---------------------------------------------------------------------
j_max(k, r)

input: k, r positive numbers
output: upper bound for abs(x), where x is a sing mod of disc Delta and denominator a for some Delta, a s.t. abs(Delta) \leq k and a \geq r
---------------------------------------------------------------------*/

j_max(k, r)={
    my(b);
    
    if(k<=0, print("Input error"); return);

    if(r<=0, print("Input error"); return);
    
    b=(exp(Pi*(k^0.5)/r)+2079);

    return(b);

}

/*--------------------------------------------------------------------
j_min(k, r)

input: k, r positive numbers
output: lower bound for abs(x), where x is a sing mod of disc Delta and denominator a for some Delta, a s.t. abs(Delta) \leq k and a \geq r
----------------------------------------------------------------------*/

j_min(k, r)={
    my(b);
    
    if(k<=0, print("Input error"); return);

    if(r<=0, print("Input error"); return);

    b=max(0,(exp(Pi*(k^0.5)/r)-2079));

return(b);
}



/*****************************************
Bounds for Section 8.3.1
*****************************************/



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Delta_x = 4*Delta and Delta_y = Delta_z = Delta case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*--------------------------------------------------
Subcase where neither of y_1, z_1 is dominant
---------------------------------------------------*/

/*---------------------------------------
lower_1_1_1(k)

input: integer k>0
output: the bound (8.2) when abs(Delta)=k
----------------------------------------*/

lower_1_1_1(k)=j_min(4*k,1)*j_min(k,1)^2;

/*----------------------------------------------------------------
upper_1_1_1(k)

input: integer k>0
output: the bound resulting from (8.3) and (8.4) when abs(Delta)=k
------------------------------------------------------------------*/

upper_1_1_1(k)=((5*j_max(4*k, 1)*j_max(k,1)*j_max(k, 2))+(18*j_max(4*k, 3)*j_max(k,1)^2));

/*------------------------------------------------------
bd_1_1_1(k)

output: integer k such that abs(Delta)<k in this subcase
-------------------------------------------------------*/

bd_1_1_1()={
    my();

    k=1;

    while(lower_1_1_1(k)<upper_1_1_1(k), k=k+1);

    return(k);

}

/*****************************************
Subcase where one of y_1, z_1 is dominant
*****************************************/

/*---------------------------------------
lower_1_1_2(k)

input: integer k>0
output: the bound (8.5) when abs(Delta)=k
----------------------------------------*/

lower_1_1_2(k)=j_min(4*k,1)*j_min(k,2)*j_min(k, 1);

/*----------------------------------------------------------------------------
upper_1_1_2(k)

input: integer k>0
output: the bound resulting from (8.6), (8.7), and (8.8) when abs(Delta)=k
----------------------------------------------------------------------------*/

upper_1_1_2(k)=((18*j_max(4*k, 3)*j_max(k,1)^2)+(2*j_max(4*k, 1)*j_max(k,3)*j_max(k, 1))+(3*j_max(4*k, 1)*j_max(k,2)^2));

/*------------------------------------------------------
bd_1_1_2()

output: integer k such that abs(Delta)<k in this subcase
------------------------------------------------------*/

bd_1_1_2()={
    my();

    k=1;

    while(lower_1_1_2(k)<upper_1_1_2(k), k=k+1);

    return(k);

}

/**************************************************
Delta_x = Delta_y = 4 Delta and Delta_z = Delta case
***************************************************/

/*-----------------------------------------
lower_1_2(k)

input: integer k>0
output: the bound (8.9) when abs(Delta)=k
-----------------------------------------*/

lower_1_2(k)=j_min(4*k,1)^2*j_min(k,2);

/*------------------------------------------------------------------
upper_1_2(k)

input: integer k>0
output: the bound resulting from (8.10) and (8.11) when abs(Delta)=k
--------------------------------------------------------------------*/

upper_1_2(k)=((j_max(4*k, 1)^2*j_max(k,3))+(22*j_max(4*k, 1)*j_max(4*k, 3)*j_max(k,1)));

/*------------------------------------------------------
bd_1_2()

output: integer k such that abs(Delta)<k in this subcase
-------------------------------------------------------*/

bd_1_2()={
    my();

    k=1;

    while(lower_1_2(k)<upper_1_2(k), k=k+1);

    return(k);

}

/*****************************************
Bounds for Section 8.3.2
*****************************************/

/*-------------------------------------------------
Case where at least one of x_1, x_2 is not dominant
--------------------------------------------------*/

/*---------------------------------------------------------
lower_2_1(k, l)

input: l \in {3/2, 2, 3, 4} and integer k > 0
output: the bound (8.14) when abs(Delta)=k for the given l
---------------------------------------------------------*/

lower_2_1(k,l)=j_min(k, 1)*j_min(k*l^2, 1)^2;

/*--------------------------------------------------------------------------
upper_2_1(k, l)

input: l \in {3/2, 2, 3, 4} and integer k > 0
output: the bound resulting from (8.15) and (8.16) when abs(Delta)=k for the given l
---------------------------------------------------------------------------*/

upper_2_1(k, l)=((j_max(k, 2)*j_max(k*l^2,1)^2)+(22*j_max(k,1)*j_max(k*l^2,1)*j_max(k*l^2, 2)));

/*-----------------------------------------------------------------------------
bd_2_1(l)

input: l \in {3/2, 2, 3, 4}
output: integer k such that abs(Delta)<k for the given l in this case
-----------------------------------------------------------------------------*/

bd_2_1(l)={
    my();

    k=1;

    while(lower_2_1(k, l)<upper_2_1(k, l), k=k+1);

    return(k);

}

/*****************************************
Case where x_1, x_2 are both dominant
*****************************************/

/*---------------------------------------------------------
lower_2_2(k, l)

input: l \in {3/2, 2, 3, 4} and integer k > 0
output: the bound (8.18) when abs(Delta)=k for the given l
---------------------------------------------------------*/

lower_2_2(k, l)=800*k^(-4)*j_min(k*l^2, 1)^2;

/*--------------------------------------------------------------------------------
upper_2_2(k, l, a_min)

input: l \in {3/2, 2, 3, 4} and integers k, a_min > 0
output: the bound resulting from (8.19), (8.20), and (8.21) when abs(Delta)=k for the given l and a_min
--------------------------------------------------------------------------------*/

upper_2_2(k, l, a_min)=((8*j_max(k, 1)*j_max(k*l^2, 1)*j_max(k*l^2, a_min))+(2*j_max(k, 2)*j_max(k*l^2, 2)^2)+(12*j_max(k, 1)*j_max(k*l^2, 2)*j_max(k*l^2, a_min)));

/*------------------------------------------------------------------------------
bd_2_2(l, a_min)

input: l \in {3/2, 2, 3, 4} and integer a_min > 0
output: integer k such that abs(Delta)<k for the given l and a_min in this case
------------------------------------------------------------------------------*/

bd_2_2(l, a_min)={
    my(k);

    k=1;

    while(lower_2_2(k, l)<upper_2_2(k, l, a_min), k=k+1);

    return(k);

}