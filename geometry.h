#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>

double dot_product(double x[3], double y[3])
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

void Quadrant0(const double &delta, double *t1min, double *t2min)
{
    double invdet{1 / delta};
    *t1min *= invdet;
    *t2min *= invdet;
    
}

void Quadrant1(const double &b,
               const double &c,
               const double &e,
               double *t1min,
               double *t2min)
{
    *t1min = 1;
    *t2min = -(b + e) / c;
    if (*t2min < 0)
    {
        *t2min = 0;
    }
    else if (*t2min > 1)
    {
        *t2min = 1;
    }
}

void Quadrant2(const double &a,
               const double &b,
               const double &c,
               const double &d,
               const double &e,
               double *t1min,
               double *t2min)
{
    if (a + b + d > 0)
    {
        *t2min = 1;
        *t1min = -(b + d) / a;
        if (*t1min < 0)
        {
            *t1min = 0;
        }
    }
    else
    {
        *t1min = 1;
        if (b + c + e > 0)
        {
            *t2min = -(b + e) / c;
            if (*t2min < 0)
            {
                *t2min = 0;
            }
        }
        else
        {
            *t2min = 1;
        }
    }
}

void Quadrant3(const double &a,
               const double &b,
               const double &d,
               double *t1min,
               double *t2min)
{
    *t2min = 1;
    *t1min = -(b + d) / a;
    if (*t1min < 0)
    {
        *t1min = 0;
    }
    else if (*t1min > 1)
    {
        *t1min = 1;
    }
}

void Quadrant4(const double &a,
               const double &b,
               const double &c,
               const double &d,
               const double &e,
               double *t1min,
               double *t2min)
{
    if (b + d < 0)
    {
        *t2min = 1;
        *t1min = -(b + d) / a;
        if (*t1min > 1)
        {
            *t1min = 1;
        }
    }
    else
    {
        *t1min = 0;
        if (c + e > 0)
        {
            *t2min = -e / c;
            if (*t2min < 0)
            {
                *t2min = 0;
            }
        }
        else
        {
            *t2min = 1;
        }
    }
}

void Quadrant5(const double &c,
               const double &e,
               double *t1min,
               double *t2min)
{
    *t1min = 0;
    *t2min = -e / c;
    if (*t2min < 0)
    {
        *t2min = 0;
    }
    else if (*t2min > 1)
    {
        *t2min = 1;
    }
}

void Quadrant6(const double &a,
               const double &c,
               const double &d,
               const double &e,
               double *t1min,
               double *t2min)
{
    if (d < 0)
    {
        *t2min = 0;
        *t1min = -d / a;
        if (*t1min > 1)
        {
            *t1min = 1;
        }
    }
    else
    {
        *t1min = 0;
        if (e < 0)
        {
            *t2min = -e / c;
            if (*t2min > 1)
            {
                *t2min = 1;
            }
        }
        else
        {
            *t2min = 0;
        }
    }
}

void Quadrant7(const double &a,
               const double &d,
               double *t1min,
               double *t2min)
{
    *t2min = 0;
    *t1min = -d / a;
    if (*t1min < 0)
    {
        *t1min = 0;
    }
    else if (*t1min > 1)
    {
        *t1min = 1;
    }
}

void Quadrant8(const double &a,
               const double &b,
               const double &c,
               const double &d,
               const double &e,
               double *t1min,
               double *t2min)
{
    if (a + d > 0)
    {
        *t2min = 0;
        *t1min = -d / a;
        if (*t1min < 0)
        {
            *t1min = 0;
        }
    }
    else
    {
        *t1min = 1;
        if (b + e < 0)
        {
            *t2min = -(b + e) / c;
            if (*t2min > 1)
            {
                *t2min = 1;
            }
        }
        else
        {
            *t2min = 0;
        }
    }
}

void Parallel(const double &a,
              const double &b,
              const double &d,
              double *t1min,
              double *t2min)
{
    if (b > 0)
    {
        if (d >= 0)
        {
            *t1min = 0;
            *t2min = 0;
        }
        else if (-d <= a)
        {
            *t1min = -d / a;
            *t2min = 0;
        }
        else
        {
            *t1min = 1;
            if (-(a + d) >= b)
            {
                *t2min = 1;
            }
            else
            {
                *t2min = -(a + d) / b;
            }
        }
    }
    else
    {
        if (-d >= a)
        {
            *t1min = 1;
            *t2min = 0;
        }
        else if (d <= 0)
        {
            *t1min = -d / a;
            *t2min = 0;
        }
        else
        {
            *t1min = 0;
            if (d >= -b)
            {
                *t2min = 1;
            }
            else
            {
                *t2min = -d / b;
            }
        }
    }
}

void NonParallel(const double &a,
                 const double &b,
                 const double &c,
                 const double &d,
                 const double &e,
                 const double &delta,
                 double *t1min,
                 double *t2min)
{
    *t1min = (b * e - c * d);
    *t2min = (b * d - a * e);
    if (*t1min >= 0)
    {
        if (*t1min <= delta)
        {
            if (*t2min >= 0)
            {
                if (*t2min <= delta)
                {
                    Quadrant0(delta, t1min, t2min);
                }
                else
                {
                    Quadrant3(a, b, d, t1min, t2min);
                }
            }
            else
            {
                Quadrant7(a, d, t1min, t2min);
            }
        }
        else
        {
            if (*t2min >= 0)
            {
                if (*t2min <= delta)
                {
                    Quadrant1(b, c, e, t1min, t2min);
                }
                else
                {
                    Quadrant2(a, b, c, d, e, t1min, t2min);
                }
            }
            else
            {
                Quadrant8(a, b, c, d, e, t1min, t2min);
            }
        }
    }
    else
    {
        if (*t2min >= 0)
        {
            if (*t2min <= delta)
            {
                Quadrant5(c, e, t1min, t2min);
            }
            else
            {
                Quadrant4(a, b, c, d, e, t1min, t2min);
            }
        }
        else
        {
            Quadrant6(a, c, d, e, t1min, t2min);
        }
    }
}

void distance_between_segments(double *u1, double *u2, double *Delta, double &t1, double &t2, double &dr2)
{
    double a, b, c, d, e, f; // coefficients of distance function
    double delta;            // determinant
    double t1min, t2min; // coordinates of minimum in distance function

    a = dot_product(u1, u1);
    b = -dot_product(u1, u2);
    c = dot_product(u2, u2);
    d = dot_product(u1, Delta);
    e = -dot_product(u2, Delta);
    f = dot_product(Delta, Delta);

    delta = a * c - b * b;

    if (delta * delta < 1e-10)
    {
        Parallel(a, b, d, &t1min, &t2min);
    }
    else
    {
        NonParallel(a, b, c, d, e, delta, &t1min, &t2min);
    }

    t1 = t1min;
    t2 = t2min;
    dr2 = a * t1min * t1min + 2 * b * t1min * t2min + c * t2min * t2min + 2 * d * t1min + 2 * e * t2min + f;
    //cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;
}
