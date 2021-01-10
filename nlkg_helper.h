#include <vector>

struct NLKG_init {
    interval y[4];
    interval t0;
};

struct NLKG_infty_init {
    interval y[2];
    interval t0;
};

const int ODE_POW = 3;
const interval DIM(3);

float min(float a, float b) {
    if (a < b) return a;
    else return b;
}

const int POS = 1;
const int NEG = -1;
const int UNC = 0;
int get_sign(interval I) {
    if (inf(I) > 0) {
        return POS;
    }
    else if (sup(I) < 0) {
        return NEG;
    }
    else {
        return UNC;
    }
}

template<typename var_type>
var_type my_pow(var_type x, int n) {
    var_type res = x;
    for (int i = 0; i < n-1; i++) {
        res = res * x;
    }
    return res;
} 


// NLKG ODE

// compute energy of a state vector
interval energy(v_blas::iVector& y) {
    interval E = my_pow(y[1], 2) / 2 + my_pow(y[0], ODE_POW+1)/(ODE_POW+1) - my_pow(y[0], 2) / 2;
    return E;
}

template<typename var_type>
var_type f(var_type y) {
    return my_pow(y, ODE_POW) - y;
}

template<typename var_type>
var_type fderiv(var_type y) {
    return ODE_POW * my_pow(y, ODE_POW-1) - 1;
}


/*
defines NLKG ODE vector field
- y'' + (2/t) y' + f(y) = 0, f(y) = y^{alpha} - y
four dimensional phase space vector because we keep track of 
- delta = (d/dy_0) y(t)
- delta'' + (2/t)delta' + f'(y) delta = 0
*/
template<typename var_type>
void NLKG(int n, var_type *yp, const var_type *y, var_type t, void *param) {
    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - f(y[0]);
    yp[2] = y[3];
    yp[3] = - ((DIM-1)/t) * y[3] - fderiv(y[0]);
}

// check if we can guarantee delta > 0 up to time t
// see sec. 3.1
float min_time_valid(interval b) {
    static interval FDERIV_ZERO = sqrt(interval(3)) / interval(3);
    float m1 = inf(sqrt(interval(6) * (b - FDERIV_ZERO) / my_pow(b, 3)));
    float m2 = inf(log(interval(2)) / (interval(3) * my_pow(b, 2)));
    return min(m1, m2);
}

// check maximum errors at time t
// see sec. 3.1
void get_errors_for_time(interval b, interval t, vector<float> &errs) {
    errs[0] = sup(my_pow(b, 5) * my_pow(t, 3) / interval(10));
    errs[1] = sup(my_pow(b, 5) * my_pow(t, 4) / interval(40));
    errs[2] = sup(my_pow(b, 6) * my_pow(t, 6) / interval(84) + 
                     my_pow(b, 4) * my_pow(t, 5) / interval(20) + 
                     my_pow(b, 4) * my_pow(t, 4) / interval(20));
    errs[3] = sup(my_pow(b, 6) * my_pow(t, 5) / interval(14) + 
                     my_pow(b, 4) * my_pow(t, 4) / interval(4) + 
                     my_pow(b, 4) * my_pow(t, 3) / interval(4));
}

// y''(0) = -f(y)/3
// y(t) ~ b - f(y) t^2/6
// y'(t) ~ -f(y)t/3
// delta(t) = 1 - f'(y) t^2/6 
// delta'(t) = -f'(y) t/3
NLKG_init NLKG_init_approx(interval b, float tol) {
    NLKG_init vals;
    interval t0 = interval(min(min_time_valid(b), 0.1));
    float cur_err = 1e6;
    vector<float> errs(4);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time(b, t0, errs);
        float max_err = 0;
        for(int i = 0; i < 4; i++) {
            if (errs[i] > max_err) max_err = errs[i];
        }

        cur_err = max_err;
    }

    // cout << "ERR:" << cur_err << endl;

    vals.t0 = t0;

    // initial calculations based on taylor approx
    vals.y[0] = b - my_pow(t0, 2) * f(b) / interval(6); // y(t0)
    vals.y[1] = -t0 * f(b) / interval(3); // y'(t0)
    vals.y[2] = 1 - my_pow(t0, 2) * fderiv(b) / interval(6); // delta(t0)
    vals.y[3] = -t0 * fderiv(b) / interval(3); // delta'(t0)

    // add in error bars manually
    // TODO: verify in this range, we have estimates in the stated directions
    vals.y[0] = interval(inf(vals.y[0]), sup(vals.y[0] + errs[0]));
    vals.y[1] = interval(inf(vals.y[1]), sup(vals.y[1] + errs[1]));
    vals.y[2] = interval(inf(vals.y[2]), sup(vals.y[2] + errs[2]));
    vals.y[3] = interval(inf(vals.y[3]), sup(vals.y[3] + errs[3]));

    return vals;
}


// NLKG infty ODE

// defines ODE that determines behavior of NLKG at inifnity
// param = s, parameter defining closeness to infinity 
//     (s = 0 is behavior at infinity)
// set y(t) = alpha * v(alpha^k t), where k = (ODE_POW - 1)/2
// then v'' + (2/t)v' + v^{ODE_POW} - alpha^{1-ODE_POW} v = 0
// set v(0) = 1, so alpha = b
// s = b^{1 - ODE_POW}
template<typename var_type>
void NLKG_inf(int n, var_type *yp, const var_type *y, var_type t, void *param) {
    interval s = *((interval *)param);

    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - my_pow(y[0], ODE_POW) + s * y[0];
}


void get_errors_for_time_infty(interval t, vector<float> &errs) {
    errs[0] = sup(my_pow(t, 4) / interval(40));
    errs[1] = sup(my_pow(t, 3) / interval(10));
}


interval energy_infty(v_blas::iVector& y, interval s) {
    interval E = my_pow(y[1], 2) / 2 + my_pow(y[0], 4)/(interval(4)) - s * my_pow(y[0], 2) / 2;
    return E;
}

// w''(0) = -(1-s) / 3
// w(t) ~ 1 - (1-s)t^2/6
// w'(t) ~ -(1-s)t/3
NLKG_infty_init NLKG_infty_init_approx(interval s, float tol) {
    NLKG_infty_init vals;
    interval t0 = 1e-2;
    float cur_err = 1e6;
    vector<float> errs(2);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_infty(t0, errs);
        float max_err = 0;
        for(int i = 0; i < 2; i++) {
            if (errs[i] > max_err) max_err = errs[i];
        }

        cur_err = max_err;
    }

    vals.t0 = t0;

    // initial calculations based on taylor approx
    vals.y[0] = interval(1) - (interval(1)-s) * my_pow(t0, 2) / interval(6); // v(t0)
    vals.y[1] = -t0 * (interval(1) - s) / interval(3); // y'(t0)

    // add in error bars manually
    // TODO: verify in this range, we have estimates in the stated directions
    vals.y[0] = interval(inf(vals.y[0]), sup(vals.y[0] + errs[0]));
    vals.y[1] = interval(inf(vals.y[1]), sup(vals.y[1] + errs[1]));

    return vals;
}
