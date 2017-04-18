//************************************************************************
//    sim5distributions.c
//------------------------------------------------------------------------
//    Date   : 10.12.2015
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2015 Michal Bursa
//************************************************************************

#ifndef CUDA

DEVICEFUNC
void distrib_init(sim5distrib* d, double(*pdf)(double), double x_min, double x_max, int N)
//! Creates a distribution based on given probability density function (PDF).
//!
//! Uses given probability density function to initiate the internal data structure with
//! calculated cummulative distribution its inverse.
//!
//! @param d pointer to structure that stores the disribution data
//! @param pdf pointer to probability density function
//! @param x_min left boundary for x
//! @param x_min right boundary for x
//! @param N number of samples to take over the interval [x_min:x_max]
{
    int i;
    double *tmp_x = (double*)calloc(N, sizeof(double));
    double *tmp_w = (double*)calloc(N, sizeof(double));
    double *tmp_pdf = (double*)calloc(N, sizeof(double));
    double *tmp_cdf = (double*)calloc(N, sizeof(double));

    d->x_min = x_min;
    d->x_max = x_max;
    d->norm  = integrate_simpson(pdf, x_min, x_max, 1e-5);

    // the x_min point
    tmp_x[0]   = x_min;
    tmp_pdf[0] = pdf(x_min);
    tmp_cdf[0] = 0.0;

    // set x[] points as nodes of gauss-legendre quadrature formula
    // this is better than linear distrubution as it helps to dump unwanted
    // spline oscillations
    gauleg(x_min, x_max, &tmp_x[1], tmp_w, N-2);

    // middle points
    for (i=1; i<N-1; i++) {
        //we have tmp_x[] already from gauleg(); for linear distibution it is: tmp_x[i]   = x_min + (double)(i)*(x_max-x_min)/(double)(N-1);
        tmp_pdf[i] = pdf(tmp_x[i])/d->norm;
        tmp_cdf[i] = integrate_simpson(pdf, x_min, tmp_x[i], 1e-6)/d->norm;
        if ((i>1) && (tmp_cdf[i]<tmp_cdf[i-1])) fprintf(stderr, "CDF problem at %d/%d (%e  %e)\n", i,N,tmp_cdf[i],tmp_cdf[i-1]);
    }

    // the x_max point
    tmp_x[N-1]   = x_max;
    tmp_pdf[N-1] = pdf(x_max);
    tmp_cdf[N-1] = 1.0;

    sim5_interp_init(&d->pdf, tmp_x, tmp_pdf, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
    sim5_interp_init(&d->cdf, tmp_x, tmp_cdf, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
    sim5_interp_init(&d->icd, tmp_cdf, tmp_x, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);

    free(tmp_x);
    free(tmp_w);
    free(tmp_pdf);
    free(tmp_cdf);
}



DEVICEFUNC
void distrib_done(sim5distrib* d)
//! Frees the internal data for the distribution.
//!
//! @param d pointer to structure that stores the disribution data
{
    sim5_interp_done(&d->pdf);
    sim5_interp_done(&d->cdf);
    sim5_interp_done(&d->icd);
}



DEVICEFUNC INLINE
double distrib_hit(sim5distrib* d)
//! Generates a functional value on interval [x_min:x_max] according to the distribution.
//!
//! With the help of precomputed cummulative distribution function it generates
//! a value form the interval [x_min:x_max] and returns tnis value.
//!
//! @param d pointer to structure that stores the disribution data
//!
//! @result value from the distribution
{
    return sim5_interp_eval(&d->icd, sim5urand());
}

#endif
