###statistics
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys, string
import unique
import math
import random
import copy

def zscore(associated_in_group,in_static_interval,total,in_flexible_interval):               
    r = float(associated_in_group)       #number of genes with this domain regulated (in one direction)
    _n = float(in_static_interval)       #measured genes in genomic interval - !!from chr_info!!
    N = float(total)                     #total number of genes analyzed
    R = float(in_flexible_interval)      #genes in the hopach interval(not measured) - !!subtract max-min or from hopach_order
    if (R-N) == 0: return 0
    elif r==0 and _n == 0: return 0
    else:
        try:
                #try:
                z = (r - _n*(R/N))/math.sqrt(_n*(R/N)*(1-(R/N))*(1-((_n-1)/(N-1))))
                return z
                #except ZeroDivisionError:
                #print path,'r,_n,R,N: ', r,_n,R,N;kill
        except ValueError:
            print (r - _n*(R/N)), _n*(R/N)*(1-(R/N))*(1-((_n-1)/(N-1)));kill

def choose(n,x):
    """Equation represents the number of ways in which x objects can be selected from a total of n objects without regard to order."""
    #(n x) = n!/(x!(n-x)!)
    f = factorial
    result = f(n)/(f(x)*f(n-x))
    return result

def sum(array):
    sum = 0
    for value in array:
        sum = sum + float(value)
    return sum
    
def avg(array):
    denominator = len(array)
    total = float(sum(array))
    average = total/denominator
    return average

def stdev(array):
    sum_dev = 0
    x_bar = avg(array)
    n = float(len(array))
    for x in array:
        x = float(x)
        sq_deviation = math.pow((x-x_bar),2)
        sum_dev = sum_dev + sq_deviation

    try:
        s_sqr = (1/(n-1))*sum_dev #s squared is the variance
        s = math.sqrt(s_sqr)
    except ZeroDivisionError:
        s = 'null'
    return s

def log_fold_conversion(array):
    try:
        new_array = []
        for log_fold in array:
            try:
                log_fold = float(log_fold)
            except ValueError:
                print log_fold, dog
            if log_fold > 0 or log_fold == 0:
                real_fold = math.pow(2,log_fold)
                new_array.append(real_fold)
            else:
                real_fold = -1/(math.pow(2,log_fold))
                new_array.append(real_fold)
    except TypeError:
        log_fold = float(array)
        if log_fold > 0 or log_fold == 0:
            new_array = math.pow(2,log_fold)
        else:
            new_array = -1/(math.pow(2,log_fold))
    return new_array

def convert_to_log_fold(array):
    list_status = 'yes'
    try:
        if len(array)>1: array = array
    except TypeError: array2 = []; array2.append(array); array = array2; list_status = 'no'
    new_array = []
    for fold in array:
        fold = float(fold)
        if fold < -1: fold = -1/fold
        #elif fold >-1 and fold <1: fold = 1
        log_fold = math.log(fold,2)
        new_array.append(log_fold)
    if list_status == 'no': return new_array[0]
    else: return new_array

def neg_folds_to_fractions(array):
    try:
        new_array = []
        for fold in array:
            try:
                fold = float(fold)
            except ValueError:
                print fold, dog
            if fold > 0:
                fold = fold
                new_array.append(fold)
            else:
                fold = -1/fold
                new_array.append(fold)
    except TypeError:
        fold = float(array)
        if fold > 0:
            new_array = fold
        else:
            new_array = -1/fold
    return new_array

def median(array):
    array.sort()
    len_float = float(len(array))
    len_int = int(len(array))
    if (len_float/2) == (len_int/2):
        try: median_val = avg([array[(len_int/2)-1],array[(len_int/2)]])
        except IndexError: median_val = ''
    else:
        try: median_val = array[len_int/2]
        except IndexError: median_val = ''
    return median_val

def int_check(value):
    val_float = float(value)
    val_int = int(value)
    if val_float == val_int:
        integer_check = 'yes'
    if val_float != val_int:
        integer_check = 'no'
    return integer_check
    
def iqr(array):
    k1 = 75
    k2 = 25
    array.sort()
    n = len(array)
    value1 = float((n*k1)/100)
    value2 = float((n*k2)/100)
    if int_check(value1) == 'no':
        k1_val = int(value1) + 1
    if int_check(value1) == 'yes':
        k1_val = int(value1)
    if int_check(value2) == 'no':
        k2_val = int(value2) + 1
    if int_check(value2) == 'yes':
        k2_val = int(value2)
    median_val = median(array)
    upper75th = array[k1_val]
    lower25th = array[k2_val]
    int_qrt_range = upper75th - lower25th
    return lower25th,median_val,upper75th,int_qrt_range

def ttest(list1,list2,tails,variance): 
    val_list1=[]
    val_list2=[]
    n1 = len(list1)
    n2 = len(list2)
    #make sure values are not strings
    for entry in list1:
        entry = float(entry)
        val_list1.append(entry)
    for entry in list2:
        entry = float(entry)
        val_list2.append(entry)
        
    if variance == 3:
        var1 = math.pow(stdev(val_list1),2)/n1
        var2 = math.pow(stdev(val_list2),2)/n2
        a1 = 1.00/(n1 - 1)
        a2 = 1.00/(n2 - 1)
        
        t = (avg(val_list1) - avg(val_list2))/math.sqrt(var1+var2)
        #calculate the degree's of freedom
        df = math.pow((var1+var2),2)/((math.pow(var1,2)/(n1-1)) + (math.pow(var2,2)/(n2-1)))

    if variance == 2:
        var1 = math.pow(stdev(val_list1),2)*(n1-1)
        var2 = math.pow(stdev(val_list2),2)*(n2-1)
        a1 = 1/n1
        a2 = 1/n2
        sp2 = (var1 + var2)/(n1+n2-2)
        sx = math.sqrt(sp2*(a1-a2))
        t = (avg(val_list1) - avg(val_list2))/sx
        df = (n1 + n2 - 2)
    return t,df,tails


def t_probability(t,df):
    """P(abs(T)<t) is equivalent to the probability between -t and +t.  So the two-sided p value for t is
    1-P(abs(T)<t)."""
    
    t = abs(t)
    df = round(df)
    if df >100:
        df = 100
    pi = 3.141592653589793238    
    if int(df)/2 == float(int(df))/2.0:
        a = 'even'
    else:
        a = 'odd'
    if a == 'even':
        sdf1 = df - 2.0
        x = 2; y = 1; z = 1; w = 1
        while x < sdf1:
            y = y*x; x = x + 2
        sdf2 = df - 3.0
        while z < sdf2:
            w = w*z; z = z + 2.0
    if a == 'odd':
        sdf1 = df - 3.0
        x = 2; y = 1; z = 1; w = 1
        while x < sdf1:
            y = y*x; x = x + 2.0
        sdf2 = df - 2.0
        while z < sdf2:
            w = w*z; z = z + 2.0

    theta = math.atan(t/math.sqrt(df))
    
    if df == 1:
        p = (2.0/pi)*theta
    if df>1 and a =='odd':
        store_var = 0
        while sdf1 > 0:
            var = (((y*(sdf1))/(w*(sdf2)))*math.pow(math.cos(theta),(sdf2)))
            store_var = store_var + var
            sdf1 = sdf1 - 2.0  
            sdf2 = sdf2 - 2.0
            try:
                w = w/sdf2
                y = y/sdf1
            except ZeroDivisionError:
                continue
        p = (2.0/pi)*(theta + math.sin(theta)*(math.cos(theta)+store_var))
        #P(abs(T)<t) = (2/pi) * (theta + sin(theta) * (cos(theta)+ (2/3)*cos(theta)^3 + ... + ((2*4*...*(nu-3))/(1*3*...*(nu-2))) * cos(theta)^(nu-2) ))
    #print w,y #3,8
    if df>1 and a =='even':
        store_var = 0
        while sdf1 > 0:
            var = (((w*(sdf2))/(y*(sdf1)))*math.pow(math.cos(theta),(sdf1)))
            #print 'stats',w,y,sdf1
            store_var = store_var + var
            sdf1 = sdf1 - 2.0  
            sdf2 = sdf2 - 2.0
            try:
                w = w/sdf2
                y = y/sdf1
            except ZeroDivisionError:
                continue

        p = math.sin(theta)*(1.0 + store_var)                  
        #p = math.sin(theta)*(1.0+(1.0/2.0)*math.pow(math.cos(theta),2.0)+((1.0*3.0)/(2.0*4.0))*math.pow(math.cos(theta),4.0) + ((w*(df-3.0))/(y*(df-2.0)))*math.pow(math.cos(theta),(df-2.0)))
        #p= sin(theta)*(1 + 1/2*cos(theta)^2 + ((1*3)/(2*4))*cos(theta)^4 + ... + ((1*3*5*...*(nu-3))/(2*4*6*...*(nu-2))) * cos(theta)^(nu-2) )

        #(1.0/2.0)*math.pow(math.cos(theta),2.0)+   ((1.0*3.0)/(2.0*4.0))*math.pow(math.cos(theta),4.0) + (1.0*3.0*5.0)/(2.0*4.0*6.0)*math.pow(math.cos(theta),(df-2.0))
    p = 1-p
    #print (2.0)/(3.0), ((w*(df-3.0))/(y*(df-2.0)))
    return p

def p_value(z):
    """A formula that is accurate to within 10^(-5) is the following: 
    P(z) = 1 - d(z)*(a1*t + a2*(t^2) + a3*(t^3)), where
    z>=0, 
    P(z) is the standard normal cumulative, 
    d(z) is the standard normal density,
    t = 1/(1+p*z),
    p = 0.33267,
    a1 = 0.4361836,
    a2 = -0.1201676,
    a3 = 0.9372980.
    This is formula 26.2.16 from Abramowitz and Stegun.  If z<0, use P(z) = 1 - P(-z).
    If they need small tail probabilities with low relative error, the 10^(-5) possible error may be too large in some cases.
    For large positive z, try
    1-P(z) = d(z)*(1/(z+1/(z+2/(z+3/(z+4/(z+5/z)))))).
    Check this in R to make sure relative errors are OK for large z.  If not, extend to 6, 7, etc. (it's a continued fractions expansion).
    d(z) = (1/(sqrt(2*pi))) * exp (-(z**2) / 2)"""
    
    p = 0.33267
    a1 = 0.4361836
    a2 = -0.1201676
    a3 = 0.9372980
    t = 1/(1+(p*z))
    pi = 3.141592653589793238

    y = (1/(math.sqrt(2*pi)))* math.exp(-(z**2)/2)

    if z >= 0:
        p_val = 1-(y*((a1*t) + a2*(math.pow(t,2)) + a3*(math.pow(t,3))))

    else:
        z = z*(-1)
        p_val = (y*((a1*t) + a2*(math.pow(t,2)) + a3*(math.pow(t,3))))

    p_val = 2*(1-p_val)
    return p_val
        
def bonferroni_p(z,correction):
    p_val = p_value(z)
    p_val = p_val*correction
    return p_val

def GrandMean(arrays):
    den = 0; num = 0; gn=0
    for array in arrays:
        x = avg(array); n = len(array); den += n; num += n*x; gn += n
        gm = num/den
    return gm,gn

def OneWayANOVA(arrays):
    f,df1,df2 = Ftest(arrays)
    p = fprob(df1,df2,f)
    return p

def Ftest(arrays):
    k = len(arrays); swsq_num=0; swsq_den=(-1)*k; sbsq_num=0; sbsq_den=(k-1); xg,ng = GrandMean(arrays)
    for array in arrays:
        n=len(array); x=avg(array); s=stdev(array); var1=(n-1)*(s**2); var2=n*((x-xg)**2)
        swsq_num += var1; swsq_den += n; sbsq_num += var2
    swsq = swsq_num/swsq_den; sbsq = sbsq_num/sbsq_den
    f = sbsq/swsq
    df1=k-1; df2=ng-k
    return f,df1,df2

if __name__ == '__main__':
    dirfile = unique    

    """where N is the total number of genes measured (all AEIs): 
    R is the total number of genes meeting the criterion (down-regulated AEIs):
    n is the total number of genes in this specific MAPP (all AEIs with miR binding sites): 
    r is the number of genes meeting the criterion in this MAPP (all down-regulated AEIs with miR binding sites): """

    N = 2106
    R = 1432
    n = 272
    r = 205

    z = zscore(r,n,N,R); p = p_value(z)
    print z, p

    N = 2106
    R = 873
    n = 272
    r = 70

    z = zscore(r,n,N,R); p = p_value(z)
    print z, p
