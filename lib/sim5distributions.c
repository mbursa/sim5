//************************************************************************
//    sim5distributions.c
//------------------------------------------------------------------------
//    Date   : 10.12.2015
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2015 Michal Bursa
//************************************************************************



DEVICEFUNC
void distrib_init(sim5distrib* d, double(*pdf)(double), double x_min, double x_max)
{
    const int N = 500;
    int i;
    double tmp_x[N];
    double tmp_w[N];
    double tmp_pdf[N];
    double tmp_cdf[N];

    d->x_min = x_min;
    d->x_max = x_max;
    d->norm  = integrate_simpson(pdf, x_min, x_max, 1e-4);

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
        //we have tmp_x[] already; this is for linear distibution: tmp_x[i]   = x_min + (double)(i)*(x_max-x_min)/(double)(N-1);
        tmp_pdf[i] = pdf(tmp_x[i])/d->norm;
        tmp_cdf[i] = integrate_simpson(pdf, x_min, tmp_x[i], 1e-4)/d->norm;
    }

    // the x_max point
    tmp_x[N-1]   = x_max;
    tmp_pdf[N-1] = pdf(x_max);
    tmp_cdf[N-1] = 1.0;

    sim5_interp_init(&d->pdf, tmp_x, tmp_pdf, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
    sim5_interp_init(&d->cdf, tmp_x, tmp_cdf, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
    sim5_interp_init(&d->icd, tmp_cdf, tmp_x, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
}



DEVICEFUNC
void distrib_done(sim5distrib* d)
{
    sim5_interp_done(&d->pdf);
    sim5_interp_done(&d->cdf);
    sim5_interp_done(&d->icd);
}



DEVICEFUNC INLINE 
double distrib_hit(sim5distrib* d)
{
    double u = (double)rand()/(double)(RAND_MAX);
    double r = sim5_interp_eval(&d->icd, u);
    return r;
}





