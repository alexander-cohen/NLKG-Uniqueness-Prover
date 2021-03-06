/**
This project consists of code accompanying the paper 
https://arxiv.org/abs/2101.08356. This is part of joint work done by 
Alex Cohen, Zhenhao Li, and Wilhelm Schlag during the summer and fall of 2020.

The purpose of this code is to prove rigorously that the first several excited 
states of the ODE y'' + (2/t)y' + y^3 - y = 0 are unique. 
*/

#include <ostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <vector>
#include "vnode.h"

using namespace vnodelp;

/*
Approximate bound states:
    b0 =  4.3373877034
    b1 = 14.1035847287
    b2 = 29.1312123369
    b3 = 49.3607112672
    b4 = 74.7723053303
*/

// SEC: CONSTANTS //

// These are used in ODE equations and energy, *not* in rigorous error bounds
const int ODE_POW = 3;
const interval DIM(3);

// tolerance allowed in initializing the algorithm with taylor approx
const double INIT_WIDTH = 1e-8;

// time to start at for equation at infinity
const double TIME_START_INFTY_EQ = 1e-2;

// minimum time to start for finitary equation
const double MIN_TIME_START = 1e-2;

// maximum allowed uncertanties in running of finitary equation
const double MAX_WIDTH_EVENTUALLY_FALLS = 1.0; 
const double MAX_WIDTH_ALGO_RUN = 0.5; 

// time steps for equation at infinity not determined dynamically 
// because we just need to lower bound the number of crossings
const double INFTY_STEP = 1.0; 

// upper bound on step for BOUND_STATE_GOOD, because we need to approach 
// y = 0 slowly
const double MAX_BOUND_STATE_STEP = 1.0; 


// used right now to make buffer for some floating point operations
// correctness is always verified after use
const double EP = 1e-10; 

// SEC: STRUCTS //

// Data structure for storing initial values for NLKG_with_delta 
struct NLKG_init_l4 {
    iVector y;
    interval t0;
    NLKG_init_l4() : y(4) {}
};

// Data structure for storing initial values for NLKG_no_delta and NLKG_inf
struct NLKG_init_l2 {
    iVector y;
    interval t0;
    NLKG_init_l2() : y(2) {}
};

// Structure for handling different intervals with different proof methods
// Used in planning section
enum pf_cases{FALLS,BOUND_GOOD,CROSSES_MANY_INFTY};
struct interval_handler {
    int method;
    int excited_state_num;
    interval the_interval;
    interval_handler(int m, interval I): method(m), the_interval(I) {};
    interval_handler(int m, int en, interval I): 
        method(m), excited_state_num(en), the_interval(I) {};
};

// SEC: INTERFACES // 
// See function implementations for documentation

// Helper functions, all must be rigorous
bool does_intersect(interval I, interval J);
bool split_halfs(interval I, interval &lh, interval &uh);
int get_sign(interval I); 
bool interval_le(interval I, interval J); 
interval min_mag(interval I); 

template<typename var_type>
var_type my_pow(var_type x, int n); 

string pf_text(int method);

// NLKG ODE functions, all must be rigorous
template<typename var_type>
var_type f(var_type y);

template<typename var_type>
var_type fderiv(var_type y);

interval energy(v_blas::iVector& y);

template<typename var_type>
void NLKG_with_delta(int n, var_type *yp, const var_type *y, 
        var_type t, void *param);

template<typename var_type>
void NLKG_no_delta(int n, var_type *yp, const var_type *y, 
        var_type t, void *param);

double min_time_valid_delta(interval b);
void get_errors_for_time_no_delta(interval b, interval t, 
        vector<interval> &errs);
void get_errors_for_time_with_delta(interval b, interval t, 
        vector<interval> &errs);
void set_initial_vals_no_delta(interval b, interval t0, iVector &y);
void set_initial_vals_with_delta(interval b, interval t0, iVector &y);
NLKG_init_l4 NLKG_init_approx_with_delta(interval b, double tol);
NLKG_init_l2 NLKG_init_approx_no_delta(interval b, double tol);

// NLKG infty ODE functions, all must be rigorous
template<typename var_type>
void NLKG_inf(int n, var_type *yp, const var_type *y, var_type t, void *param);
void get_errors_for_time_infty(interval t, vector<interval> &errs);
interval energy_infty(v_blas::iVector& y, interval beta);
NLKG_init_l2 NLKG_infty_init_approx(interval beta, double tol);

// Planning section functions, need not be rigorous
interval infty_step(iVector y, interval t, interval beta);
interval find_nth_excited_state_binsearch(int N, double tol, 
    AD *ad, VNODE *Solver, double lower_bound, double upper_bound);
interval find_nth_excited_state(int N, double tol, AD *ad, VNODE *Solver);
double determine_min_height_for_crossings(int n, AD *ad, VNODE *Solver);
vector<interval_handler> make_n_first_excited_plan(int n);

// Proving section functions, must be rigorous
interval time_max_one_cross(iVector y, interval t);
int crossing_number_smallb(interval b, int maxc, AD *ad, VNODE *Solver);
int crossing_number_infty(interval beta, int max, AD *ad, VNODE *Solver);
bool prove_crosses_many_infty(int n, interval beta, AD *ad, VNODE *Solver);
bool prove_eventually_falls(interval b, AD *ad, VNODE *Solver);
bool prove_eventually_falls_bisection(interval ran_fall, 
        AD *ad, VNODE *Solver);
bool bound_state_good(int n, interval b, AD *ad, VNODE *Solver, bool verbose);

bool verify_drop(interval I);
bool verify_boundstate_good(interval I, int n);
bool verify_crosses_many_infty(interval I, int n);
bool execute_first_n_excited_plan(int n, vector<interval_handler> plan);

// Testing code / code for graphing
void sample_solve(interval b, interval tend);
void NLKG_init_test();
void make_sol_data(interval b, double T, double step);
int make_output_dat_comp_finite_inf();

// Main functions
int run_uniqueness_prover(int argc, char *argv[]);
int make_N3_output_for_graphs();

// SEC: HELPER FUNCTIONS //

/**
(Must be rigorous)
Check if I and J intersect (wrapper function for vnode's `intersect')
*/
bool does_intersect(interval I, interval J) {
    interval _(0);
    return intersect(_, I, J);
}

/**
(Must be rigorous)
split interval I into an upper and lower half who's union contains I
*/
bool split_halfs(interval I, interval &lh, interval &uh) {
    double mp = midpoint(I);
    lh = interval(inf(I), mp); 
    uh = interval(mp, sup(I));
    if (!does_intersect(lh, uh)) return false; // not contiguous

    return true;
}

/**
(Must be rigorous)
Returns the sign of interval I
returns: UNC if I contains 0
         POS if every element of I is > 0
         NEG is every element of I is < 0
*/
enum signs{POS, NEG, UNC};
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

/**
(Must be rigorous)
Checks if x < y for all x in I, y in J
*/
bool interval_le(interval I, interval J) {
    return get_sign(J - I) == POS;
}

/**
(Must be rigorous)
Returns an interval containining a lower bound on the absolute value of I
*/
interval min_mag(interval I) {
    int s = get_sign(I);
    if (s == UNC) {
        return interval(0);
    }
    else if (s == POS) {
        return interval(inf(I));
    }
    else { 
        return -interval(sup(I));
    }
}


/**
(Must be rigorous)
Evaluate a power by repeated multiplication.
Necessary for taking powers of negative numbers in ODE computations
*/
template<typename var_type>
var_type my_pow(var_type x, int n) {
    var_type res = x;
    for (int i = 0; i < n-1; i++) {
        res = res * x;
    }
    return res;
} 


// SEC: NLKG ODE //

/**
(Must be rigorous)
Compute energy of a state vector for NLKG equation
*/
interval energy(v_blas::iVector& y) {
    interval E = my_pow(y[1], 2) / interval(2) + 
                my_pow(y[0], ODE_POW+1)/interval(ODE_POW+1) - 
                my_pow(y[0], 2) / interval(2);
    return E;
}

/**
(Must be rigorous)
Evaluate f(y) = y^3 - y
*/
template<typename var_type>
var_type f(var_type y) {
    return my_pow(y, ODE_POW) - y;
}

/** 
(Must be rigorous)
Evaluate f'(y) = 3y^2 - 1
*/
template<typename var_type>
var_type fderiv(var_type y) {
    return ODE_POW * my_pow(y, ODE_POW-1) - interval(1);
}

/**
(Must be rigorous)
Defines NLKG ODE vector field
y'' + (2/t) y' + f(y) = 0, f(y) = y^3 - y
Four dimensional phase space vector: y, y', delta, delta'
- delta = (d/db) y(t), b = y(0)
- delta'' + (2/t)delta' + f'(y) delta = 0
*/
template<typename var_type>
void NLKG_with_delta(int n, var_type *yp, const var_type *y, 
        var_type t, void *param) {
    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - f(y[0]);
    yp[2] = y[3];
    yp[3] = - ((DIM-1)/t) * y[3] - fderiv(y[0]) * y[2];
}

/**
(Must be rigorous)
Defines NLKG ODE vector field without delta
y'' + (2/t) y' + f(y) = 0, f(y) = y^3 - y
Two dimensional phase space vector: y, y'
*/
template<typename var_type>
void NLKG_no_delta(int n, var_type *yp, const var_type *y, 
        var_type t, void *param) {
    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - f(y[0]);
}

/**
(Must be rigorous)
Helper function for initializing the NLKG ODE past the singularity at t = 0
Finds a time T > 0 such that for all b in given interval, 
    0 < delta(t) < 1 for 0 < t < T,
and we can also verify delta(t), delta'(t) move in the negative direction and 
have appropriate error bounds in [0, T].
See sec 3.3 in the paper
*/
double min_time_valid_delta(interval b) {
    static interval SQRT3 = sqrt(interval(3));
    double m1 = inf( sqrt( interval(6) * (SQRT3 * b - interval(1)) / 
                            (SQRT3 * b * (b*b - 1)) ) );
    double m2 = inf(log(interval(4)) / (SQRT3 * b));
    return min(m1, m2);
}

/**
(Must be rigorous)
Find maximum error bounds up to time T for ODE with no delta
See sec 3.3 in the paper
*/
void get_errors_for_time_no_delta(interval b, interval t, 
        vector<interval> &errs) {
    errs[0] = my_pow(b, 5) * my_pow(t, 4) / interval(40);
    errs[1] = my_pow(b, 5) * my_pow(t, 3) / interval(10);
}

/**
(Must be rigorous)
Find maximum error bounds up to time T for ODE with delta
See sec 3.3 in the paper
*/
void get_errors_for_time_with_delta(interval b, interval t, 
        vector<interval> &errs) {
    get_errors_for_time_no_delta(b, t, errs);
    errs[2] = my_pow(b, 4) * my_pow(t, 4) / interval(8);
    errs[3] = my_pow(b, 4) * my_pow(t, 3) / interval(2);
}

/**
(Must be rigorous)
Set initial values for ODE using second Picard approximation, no delta
See sec 3.3 in the paper
*/
void set_initial_vals_no_delta(interval b, interval t0, iVector &y) {
    y[0] = b - my_pow(t0, 2) * f(b) / interval(6); // y(t0)
    y[1] = -t0 * f(b) / interval(3); // y'(t0)
}

// Set initial values for ODE using second Picard approximation, with delta
// See sec 3.3 in the paper
void set_initial_vals_with_delta(interval b, interval t0, iVector &y) {
    set_initial_vals_no_delta(b, t0, y);
    y[2] = 1 - my_pow(t0, 2) * fderiv(b) / interval(6); // delta(t0)
    y[3] = -t0 * fderiv(b) / interval(3); // delta'(t0)
}


/**
(Must be rigorous)
Find initial time and y interval with explicit error bounds to move past 
the singularity at t = 0.
Given an initial interval B and tolerance tol, 
gives initial time interval T, initial y vector Y, such that y_b(t) in Y
for all b in B, t in T, y in Y. 
Also guarantees that the additional errors incurred from the Picard 
approximation is no more than `tol'

Second Picard approximations:
y(t) ~ b - f(y) t^2/6
y'(t) ~ -f(y)t/3
delta(t) = 1 - f'(y) t^2/6 
delta'(t) = -f'(y) t/3
*/
NLKG_init_l4 NLKG_init_approx_with_delta(interval b, double tol) {
    NLKG_init_l4 vals;
    interval t0 = interval(min(min_time_valid_delta(b), MIN_TIME_START));
    double cur_err = 1e6;
    vector<interval> errs(4);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_with_delta(b, t0, errs);
        double max_err = 0;

        for(int i = 0; i < 4; i++) {
            if (sup(errs[i]) > max_err) max_err = sup(errs[i]); 
        }

        cur_err = max_err;
    }

    vals.t0 = t0;

    // initial calculations based on taylor approx
    set_initial_vals_with_delta(b, t0, vals.y);

    // add in error bars manually
    // Note: we know in this range, approximation is in the stated direction.
    // i.e., for y* the actual y, y[i] < y*[i] < y[i] + errs[i];
    for(int i = 0; i < 4; i++) {
        vals.y[i] = interval(inf(vals.y[i]), sup(vals.y[i] + errs[i]));
    }

    return vals;
}

/**
(Must be rigorous)
Performs the same functionality as 'NLKG_init_approx_with_delta' except 
initializes ODE with no delta. Useful for paarts of the analysis that do not
require keeping track of delta.
*/
NLKG_init_l2 NLKG_init_approx_no_delta(interval b, double tol) {
    NLKG_init_l2 vals;
    interval t0 = MIN_TIME_START;
    double cur_err = 1e6;
    vector<interval> errs(2);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_no_delta(b, t0, errs);
        double max_err = 0;

        for(int i = 0; i < 2; i++) {
            if (sup(errs[i]) > max_err) max_err = sup(errs[i]); 
        }

        cur_err = max_err;
    }

    vals.t0 = t0;

    // initial calculations based on taylor approx
    set_initial_vals_no_delta(b, t0, vals.y);

    // add in error bars manually
    for(int i = 0; i < 2; i++) {
        vals.y[i] = interval(inf(vals.y[i]), sup(vals.y[i] + errs[i]));
    }


    return vals;
}


// SEC: NLKG INFTY ODE //

/**
(Must be rigorous)
Defines ODE that determines behavior of NLKG at inifnity
param = beta, parameter defining closeness to infinity 
        (beta = 0 is behavior at infinity)
        beta = 1/b
set y(t) = alpha * v(alpha^k t), where k = (ODE_POW - 1)/2
then w'' + (2/t)w' + w^{ODE_POW} - alpha^{1-ODE_POW} w = 0
set w(0) = 1, so alpha = b
Set beta = 1/b = 1/alpha
*/
template<typename var_type>
void NLKG_inf(int n, var_type *yp, const var_type *y, var_type t, void *param) {
    interval beta = *((interval *)param);
    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - 
            my_pow(y[0], ODE_POW) + 
            my_pow(beta, ODE_POW-1) * y[0];
}

/**
(Must be rigorous)
Evaluates maximum errors at time t for infty equation, see sec 3.3
*/
void get_errors_for_time_infty(interval t, vector<interval> &errs) {
    errs[0] = my_pow(t, 4) / interval(40);
    errs[1] = my_pow(t, 3) / interval(10);
}

/**
(Must be rigorous)
Evaluates energy of infty equation
*/
interval energy_infty(v_blas::iVector& y, interval beta) {
    interval E = my_pow(y[1], 2) / interval(2) + 
                    my_pow(y[0], 4)/(interval(4)) - 
                    my_pow(beta, ODE_POW-1) * my_pow(y[0], 2) / interval(2);
    return E;
}


/**
(Must be rigorous)
Performs similar functionality to `NLKG_init_approx_with_delta' except for the
equation at infinity. See sec 3.3 in the paper. 

Second order Picard approximants:
w''(0) = -(1-beta^2) / 3
w(t) ~ 1 - (1-beta^2)t^2/6
w'(t) ~ -(1-beta^2)t/3
*/
NLKG_init_l2 NLKG_infty_init_approx(interval beta, double tol) {
    NLKG_init_l2 vals;
    interval t0 = TIME_START_INFTY_EQ;
    double cur_err = 1e6;
    vector<interval> errs(2);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_infty(t0, errs);
        double max_err = 0;
        for(int i = 0; i < 2; i++) {
            if (sup(errs[i]) > max_err) max_err = sup(errs[i]);
        }

        cur_err = max_err;
    }

    vals.t0 = t0;

    // initial calculations based on taylor approx
    vals.y[0] = interval(1) - (interval(1)-my_pow(beta, 2)) * 
                                my_pow(t0, 2) / interval(6); // v(t0)
    vals.y[1] = -t0 * (interval(1) - my_pow(beta, 2)) / interval(3); // y'(t0)

    // add in error bars manually
    vals.y[0] = interval(inf(vals.y[0]), sup(vals.y[0] + errs[0]));
    vals.y[1] = interval(inf(vals.y[1]), sup(vals.y[1] + errs[1]));

    return vals;
}

// SEC: UNIQUENESS PROVER CODE // 

/**
(Must be rigorous)
Determine an amount of time such that after integrating for this long, 
we cross zero at most once. 
Must be rigorous to ensure observed zero crossings match actual zero crossings
*/
interval time_max_one_cross(iVector y, interval t) {
    // if we travel distance to zero + 1 we are at most in the left well so 
    // won't cross zero more than once.  
    static interval MIN_ENERGY = -interval(1/2);
    interval min_d = min_mag(y[0]) + interval(1);

    // overestimate the energy, as this will underestimate the step size
    interval E = interval(sup(energy(y))); 

    interval maxv = sqrt(2 * (E - MIN_ENERGY));
    interval step = min_d / maxv;
    return interval(inf(step));
}

/**
(Need not be rigorous)
Determines step size for equation at infinity
Steps at least INFTY_STEP, possible more if its determined that is safe
Does not need to be rigorous, because we only need to observe a minimum number
of crossings for the infty equation, not verify those are all of the zero 
crossings.
*/
interval infty_step(iVector y, interval t, interval beta) {
    // energy minimized at y^3 - beta^2y = 0, or y = beta;
    interval MIN_ENERGY = -my_pow(beta, 4) / interval(2);

    interval min_d = min_mag(y[0]) + interval(sup(beta));
    interval E = interval(sup(energy_infty(y, beta)));
    interval maxv = sqrt(interval(2) * (E - MIN_ENERGY));
    interval step_guess = min_d / maxv;

    // always take steps of at least length INFTY_STEP (but possibly larger)
    double step_actual = max(inf(step_guess), INFTY_STEP);
    return interval(step_actual);
}

/**
(Must be rigorous)
Computes crossing number with finitary equation
Params:
    - `b' input interval, choices of y(0)
    - `maxc' maximum number of crossings to observe
    - `ad', `Solver', VNODE objects
Returns:
    - maxc if at least maxc zero crossings are observed for all y0 in `b'
    - N if it can be verified that for all y0 in `b', y(t) crosses zero 
        exactly N times before falling into an energy well
    - -1 if integration fails. This will always happen if `b' 
        contains a bound state of order m < maxc
*/
int crossing_number_smallb(interval b, int maxc, AD *ad, VNODE *Solver) {
    // cout << "crossing number smallb: " << b << endl;
    Solver->setFirstEntry();
    NLKG_init_l2 init_val = NLKG_init_approx_no_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    int num_crossings = 0;
    bool pos = true;
    while (sup(energy(y)) > 0 && 
        num_crossings < maxc) {
        if (width(y[0]) > MAX_WIDTH_ALGO_RUN) return -1;

        interval step = time_max_one_cross(y, t);
        interval new_t = t + step;

        Solver->integrate(t, y, new_t);
        
        if(!Solver->successful()) {
            cout << "Crossing number smallb failed: " << 
                t << ", " << new_t << endl;
            for (int i = 0; i < 4; i++) {
                cout << y[i] << endl;
            }
            return -1;
        }

        int sign = get_sign(y[0]);
        if (sign == POS) {
            if (pos == false) num_crossings += 1;
            pos = true;
        }
        else if (sign == NEG) {
            if (pos == true) num_crossings += 1;
            pos = false;
        }
    }

    return num_crossings;
}

/**
(Need not be rigorous)
Uses a binary search to find an interval containing the Nth bound state
N is the number of zero crossings (so N = 0 is ground state)

The output interval must have width < tol
*/
interval find_nth_excited_state_binsearch(int N, double tol, 
    AD *ad, VNODE *Solver, double lower_bound, double upper_bound) {
     while (upper_bound - lower_bound > tol) {
        double mid = (upper_bound + lower_bound) / 2;
        // cout << setprecision(8) << tol << ", " << 
        //     N << "; " << 
        //     lower_bound << ", " << upper_bound << "; " << 
        //     mid << endl;
        int nc_mid = crossing_number_smallb(mid, N+1, ad, Solver);
        // cout << "MID: " << mid << ", " << nc_mid << endl;
        if (nc_mid == -1) {
            cout << "failed to find " << N << "th bound state" << endl;
            return -1;
        }

        if (nc_mid <= N) {
            lower_bound = mid;
        }
        else {
            upper_bound = mid;
        }
    }

    return interval(lower_bound - tol, upper_bound + tol);
}

/**
(Need not be rigorous)
Wrapper function for `find_nth_excited_state_binsearch'
*/
interval find_nth_excited_state(int N, double tol, AD *ad, VNODE *Solver) {
    double upper_bound = determine_min_height_for_crossings(N+1, ad, Solver);
    double lower_bound = 1.0;
    return find_nth_excited_state_binsearch(N, tol, ad, Solver, 
            lower_bound, upper_bound);
}


/**
(Need not be rigorous)
Helper function for `find_nth_excited_state_binsearch'.
Finds an upper bound on the nth bound state
*/
double determine_min_height_for_crossings(int n, AD *ad, VNODE *Solver) {
    double b = 1.0;
    while (crossing_number_smallb(interval(b), n, ad, Solver) < n)  {
        b *= 2;
    }
    return b;
}


/**
(Must be rigorous)
Determines number of crossings using infinitary equation
Guarantees that it will cross *at least* as many times as the return value, 
could potentially cross more times.

Params:
    - `beta' an interval to check for crossing number, corresponds to 
        interval 1/beta of b-values. Should be a small interval,
        will try to integrate.
    - `max' maximum number of crossings to check for
    - `ad', `Solver' VNODE objects

Returns:
    - `max' if for all beta in `beta', w_beta(t) has at least max zero crossings
    - N if for all beta in `beta', w_beta(t) has at least N zero crossings
    - -1 if integration fails
*/
int crossing_number_infty(interval beta, int max, AD *ad, VNODE *Solver) {
    ad->eval(&beta);
    Solver->setFirstEntry();

    NLKG_init_l2 init_val = NLKG_infty_init_approx(beta, INIT_WIDTH);

    interval t = init_val.t0;
    iVector y = init_val.y;

    int num_crossings = 0;
    bool pos = true;
    while (num_crossings < max && sup(energy_infty(y, beta)) > 0) {
        if (width(y[0]) > MAX_WIDTH_ALGO_RUN) return false;

        interval new_t = t + infty_step(y, t, beta);
        Solver->integrate(t, y, new_t);
        
        if(!Solver->successful()) {
            return num_crossings;
        }

        int sign = get_sign(y[0]);
        if (sign == POS) {
            if (pos == false) num_crossings += 1;
            pos = true;
        }
        else if (sign == NEG) {
            if (pos == true) num_crossings += 1;
            pos = false;
        }
    }
    cout << "finished at time: " << t << endl;
    return num_crossings;
}

/** 
(Must be rigorous)
Verify that within the interval beta, there are at least `n' crossings for the 
infty equation. Makes use of `crossing_number_infty'. 

Proceeds by recursively bisecting `beta' and 
trying to run `crossing_number_infty'.
*/
bool prove_crosses_many_infty(int n, interval beta, AD *ad, VNODE *Solver) {
    if (width(beta) < EP * 4) { 
        cout << "Interval too small, cannot verify crosses " << 
                "many times at infinity: " << beta << endl;
        return false; // if interval is too small we can't bisect further
    }
    int nc = crossing_number_infty(beta, n, ad, Solver);
    if (nc == n) {
        cout << "Proved crosses at least " << n << " times: " << beta << endl;
        return true;
    }
    else {
        // bisect interval into two parts
        interval lh;
        interval uh;
        bool succesfull_split = split_halfs(beta, lh, uh);

        if (!succesfull_split) {
            cout << "Could not verify that upper and " << 
                    "lower halves cover interval" << endl;
            return false;
        }

        bool lower_good = prove_crosses_many_infty(n, lh, ad, Solver);
        if (!lower_good) return false;
        bool upper_good = prove_crosses_many_infty(n, uh, ad, Solver);
        if (!upper_good) return false;

        return true;
    }
}

/**
(Must be rigorous)
Prove that starting in the interval height `b', the solution eventually falls 
into one energy well.

Params:
    - `b' initial interval from which to integrate, should be small
    - `ad', `Solver' VNODE objects

Returns:
    - true if for all b in `b', y_b(t) eventually has negative energy
    - false if the above could not be verified
*/
bool prove_eventually_falls(interval b, AD *ad, VNODE *Solver) {
    if (width(b) > MAX_WIDTH_EVENTUALLY_FALLS) {
        return false;
    }

    Solver->setFirstEntry();
    NLKG_init_l2 init_val = NLKG_init_approx_no_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    while (sup(energy(y)) > 0) {
        if (width(y[0]) > MAX_WIDTH_ALGO_RUN) return false;

        interval new_t = t + max(sup(time_max_one_cross(y, t)), 1.0);
        Solver->integrate(t, y, new_t);
        if(!Solver->successful()) {
            return false;
        }
    }

    return true;
}

/**
(Must be rigorous)
Verify that for all b in `ran_fall', y_b(t) eventually has negative energy.
Recursively bisects `ran_fall' and tries to apply `prove_eventually_falls'.
*/
bool prove_eventually_falls_bisection(interval ran_fall, 
        AD *ad, VNODE *Solver) {
    // if interval is too small we can't bisect further
    if (width(ran_fall) < EP * 4) return false;

    // cout << "checking: " << ran_fall << endl;
    // try running directly on this interval first
    // if interval is too big this will fail, so we bisect
    bool res = prove_eventually_falls(ran_fall, ad, Solver);
    if (res) {
        cout << "Proved falls: " << ran_fall << endl;
        return true;
    }

    // bisect interval into two parts
    interval lh;
    interval uh;
    bool succesfull_split = split_halfs(ran_fall, lh, uh);

    if (!succesfull_split) {
        cout << "Could not verify that upper and " << 
                "lower halves cover interval" << endl;
        return false;
    }

    bool lower_good = prove_eventually_falls_bisection(lh, ad, Solver);
    if (!lower_good) return false;

    bool upper_good = prove_eventually_falls_bisection(uh, ad, Solver);  
    if (!upper_good) return false;

    return true;
}

/**
(Must be rigorous)
Verify that within the interval b, there is at most one bound state, and if it 
exists, that bound state must be of order `n'. 

Params:
    - `n' the order of the bound state to prove uniqueness for
    - `b' an interval which should contain the nth bound state, and for which 
        we would like to prove uniqueness
    - `ad', `Solver' VNODE objects
    - `verbose' specifies whether or not to log solution output

Returns:
    - true if it can be verified that within `b' there is at most one bound 
        state, and if it exists, it must be of order `n'
    - false if the above could not be verified
*/
bool bound_state_good(int n, interval b, AD *ad, VNODE *Solver, bool verbose) {
    
    NLKG_init_l4 init_val = NLKG_init_approx_with_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    // cout << y[0] << ", " << t << endl;

    Solver->setHmin(0);
    Solver->setFirstEntry();

    int num_crossings = 0;
    bool pos = true;

    while (true) {
        if (width(y[0]) > MAX_WIDTH_ALGO_RUN) return false;
        // step at most time_max_one_cross, maybe less
        interval step = time_max_one_cross(y, t);
        if (MAX_BOUND_STATE_STEP < inf(step)) {
            step = interval(MAX_BOUND_STATE_STEP);
        }

        interval new_t = t + step;
        Solver->integrate(t,y,new_t);
        interval E = energy(y);

        if (!Solver->successful()) {
            return false;
        }

        // keep track of number of crossings
        int sign = get_sign(y[0]);
        if (sign == POS) {
            if (pos == false) num_crossings += 1;
            pos = true;
        }
        else if (sign == NEG) {
            if (pos == true) num_crossings += 1;
            pos = false;
        }

        // must have exactly n crossings at time of verification
        if (num_crossings < n) continue; 
        if (num_crossings > n) return false; // too many crossings


        // if E < 0, return false
        if (interval_le(E, interval(0))) {
            return false;
        }

        // all variables must have definite signs
        bool correct_signs = true;
        for(int i = 0; i < 4; i++) {
            if(get_sign(y[i]) == UNC) correct_signs = false;
        }
        if (!correct_signs) continue;

        // this is valid because we already checked they have definite signs
        // check that y[1], y[2], y[3] have the opposite sign from y[0]
        for (int i = 1; i < 4; i++) {
            if (get_sign(y[i]) == get_sign(y[0])) correct_signs = false;
        }
        if (!correct_signs) continue;

        // check |y(t)| < 1/2
        if(!interval_le(mag(y[0]), interval(1)/interval(2))) continue;

        // check 0 < E(T) < 1/4
        if(!interval_le(E, interval(1)/interval(4))) continue;

        // energy bound
        interval falls_over_v = (new_t  - interval(2) * log(new_t) + 
            (interval(3)/interval(2))) * E;
        if (!interval_le(falls_over_v, interval(3) / interval(8))) {
            continue;
        }
        if (verbose) {
            cout << "Succesfully proved at most one bound state in range: " << 
                    b << endl;
            cout << "Ending time = " << new_t << ", ending y = " << endl;
            for(int i = 0; i < 4; i++) {
                cout << inf(y[i]) << ", " << sup(y[i]) << endl;
            }
        }

        // if we reach this point, we verified uniqueness of bound state n
        return true;
    }
}

// SEC: PLANNING SECTION // 


string pf_text(int method) {
    if (method == FALLS) {
        return "FALLS";
    }
    else if(method == BOUND_GOOD) {
        return "BOUND STATE GOOD";
    }
    else {
        return "CROSSES MANY INFTY";
    }
}


/**
(Need not be rigorous)
Creates a plan for proving the first `n' excited states are unique

Returns a vector of intervals and proof methods for each of those intervals
*/
vector<interval_handler> make_n_first_excited_plan(int n) {
    static interval ENERGY_ZER0 = sqrt(interval(2));

    vector<interval_handler> pf_plan;

    AD *ad = new FADBAD_AD(2,NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver = new VNODE(ad);

    AD *ad_wd = new FADBAD_AD(4, NLKG_with_delta, NLKG_with_delta);
    VNODE *Solver_wd = new VNODE(ad_wd);

    double b_checked_upto = inf(ENERGY_ZER0);
    double buffer_size = 2.0;

    for(int i = 0; i <= n; i++) {
        cout << "On bound state #" << i << endl;
        cout << "Adaptively determining bound state width" << endl;
        interval nth_excited_approx;
        double bound_state_tol = 0.5;
        double lower = 1.0;
        double upper = determine_min_height_for_crossings(i+1, ad, Solver);
        while (true) {
            // cout << lower << ", " << upper << ", " << bound_state_tol << endl;
            nth_excited_approx = find_nth_excited_state_binsearch(i, 
                bound_state_tol, ad, Solver, lower, upper);
            bool res = bound_state_good(i, nth_excited_approx, 
                        ad_wd, Solver_wd, false);
            if (res) break; // small enough to verify bound state good
            lower = inf(nth_excited_approx);
            upper = sup(nth_excited_approx);
            bound_state_tol /= 2;
        }
        cout << "Decided on bound state width for #" << i << ": " << 
                    width(nth_excited_approx) << endl;
        

        interval drop_before = interval(b_checked_upto, 
            inf(nth_excited_approx) + EP);

        pf_plan.push_back(interval_handler(FALLS, drop_before));
        cout << "created interval to check: FALLS, " << drop_before << endl;

        pf_plan.push_back(interval_handler(BOUND_GOOD, i, nth_excited_approx));
        cout << "created interval to check: BOUND STATE GOOD, " << 
                nth_excited_approx << endl;

        b_checked_upto = sup(nth_excited_approx - EP);
    }

    // check that it falls in a ``buffer interval'' above the last bound state
    interval buffer_int = interval(b_checked_upto - EP, 
                                    b_checked_upto + buffer_size);
    pf_plan.push_back(interval_handler(FALLS, buffer_int));
    cout << "created interval to check: BUFFER FALLS, " << buffer_int << endl;
    b_checked_upto = sup(buffer_int) - EP; 

    // construct interval of s to check
    // need to check b_checked_upto --> infty
    // under inversion, this becomes 0 --> 1/b_checked_upto
    interval s_int = interval(0, sup(1/interval(b_checked_upto)) + 5 * EP);
    pf_plan.push_back(interval_handler(CROSSES_MANY_INFTY, s_int));
    cout << "created interval to check: INFTY CROSSES MANY, " << s_int << endl;

    return pf_plan;
}

// SEC: PROVING SECTION // 

/**
(Must be rigorous)
Verify that in interval I, all solutions eventually fall in an energy well
Wrapper function for `prove_eventually_falls_bisection'
*/
bool verify_drop(interval I) {
    cout << "Will verify drops: " << I << endl;
    AD *ad = new FADBAD_AD(2,NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver = new VNODE(ad);

    return prove_eventually_falls_bisection(I, ad, Solver);
}

/**
(Must be rigorous)
Verify that in the interval I there is at most one bound state, and if it exists
it is of order `n'
Wrapper function for `bound_state_good'
*/
bool verify_boundstate_good(interval I, int n) {
    cout << "Will verify #" << n << " bound state good: " << I << endl;

    // verify that the conditions of Lemmas 6, 7 hold, so at most 
    // one nth bound state, and no (n+1)th bound state
    AD *ad = new FADBAD_AD(4, NLKG_with_delta, NLKG_with_delta);
    VNODE *Solver = new VNODE(ad);
    bool at_most_one = bound_state_good(n, I, ad, Solver, true);

 
    return at_most_one;
}

/**
(Must be rigorous)
Verify that in the interval I we cross at least n+1 times, using infty eq.
Wrapper function for `prove_crosses_many_infty'
*/
bool verify_crosses_many_infty(interval I, int n) {
    cout << "Will verify crosses many infinity: " << I << endl;

    AD *ad = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &I);
    VNODE *Solver = new VNODE(ad);

    return prove_crosses_many_infty(n+1, I, ad, Solver);
}

/**
(Must be rigorous)
Executes plan to prove the first `n' excited states are unique

Params:
    - `n' highest excited state to prove unique. E.g. if n = 1, will prove 
        ground state and first excited state are unique
    - `plan' plan for proving the first `n' excited states are unique, should
        be created by `make_n_first_excited_plan'

Returns:
    - true if it is verified using the plan that the first n excited states 
        are unique
    - false if the above cannot be verified
*/
bool execute_first_n_excited_plan(int n, vector<interval_handler> plan) {
    // verify that b = sqrt(2), zero energy, is covered
    int plan_len = plan.size();

    interval SQRT2 = sqrt(interval(2)); 
    if (!subseteq(SQRT2, plan[0].the_interval))  {
        cout << "b = sqrt(2) not covered." << endl;
        return false;
    }
    cout << "Verified endpoint b = sqrt(2) is covered" << endl;

    // verify that subsequent intervals intersect
    for(int i = 0; i < plan_len - 2; i++) {
        if (!does_intersect(plan[i].the_interval, plan[i+1].the_interval)) {
            cout << "Intervals not contiguous: " << i << endl;
            return false;
        }
    }
    interval infty_interval = plan[plan_len - 1].the_interval;
    interval sub_infty_interval = interval(midpoint(infty_interval), 
        sup(infty_interval) - EP);
    if (!subseteq(sub_infty_interval, infty_interval)) {
        cout << "Failed in making sub infinity interval" << endl;
    }
    interval outer_infty_interval = 1 / sub_infty_interval;
    if (!does_intersect(plan[plan_len - 2].the_interval,outer_infty_interval)) {
        cout << "Intervals not contiguous: end to infinity" << endl;
        cout << plan[plan_len - 2].the_interval << endl;
        cout << outer_infty_interval << endl;
        return false;
    }


    cout << "Verified intervals are contiguous" << endl << endl;

    // Execute proof method for each interval
    for(uint i = 0; i < plan.size(); i++) {
        bool success = false;
        if (plan[i].method == FALLS) {
            success = verify_drop(plan[i].the_interval);
        }
        else if (plan[i].method == BOUND_GOOD) {
            success = verify_boundstate_good(plan[i].the_interval, 
                                            plan[i].excited_state_num);
        }
        else if (plan[i].method == CROSSES_MANY_INFTY) {
            success = verify_crosses_many_infty(plan[i].the_interval, n);
        }

        cout << "In interval: " << plan[i].the_interval << 
                ", tried to prove: " << pf_text(plan[i].method) << 
                ", did: " << (success ? "SUCCEED" : "FAIL") << endl << endl;;
        if (!success) {
            return false;
        }
    }

    cout << "\nSuccesfully proved first " << n << 
            " bound states are unique." << endl;

    return true;
}

// SEC: TESTING & HELPER CODE //

// sample NLKG solver
void sample_solve(interval b, interval tend) {
    // instantiate ODE solver
    AD *ad = new FADBAD_AD(4,NLKG_with_delta,NLKG_with_delta);
    VNODE *Solver = new VNODE(ad);

    Solver->setFirstEntry();
    NLKG_init_l4 init_val = NLKG_init_approx_with_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    // integrate solver up to end time
    Solver->integrate(t,y,tend);

    // if enclosures are too big integration will fail
    if(!Solver->successful()) {
        cout << "VNODE-LP could not reach t = " << tend <<endl;
        return;
    }

    // print succesful enclosure
    cout << "Solution enclosure at t = " << t << endl;
    printVector(y);

    // print ending energy
    interval E = energy(y);
    cout << "Energy: " << E << endl;
}

// function to test NLKG initialization code
void NLKG_init_test() {
    for (double b = 2.0; b < 100.0; b += 2.0) {
        NLKG_init_l4 vals = NLKG_init_approx_with_delta(interval(b, b + 1e-2), 
                                                        INIT_WIDTH);
        cout << b << ", " << vals.t0 << endl;
        for (int i = 0; i < 4; i++) {
            cout << i << ": " << vals.y[i] << endl;
        }
    }
}

// function to make solution data for graphing
void make_sol_data(interval b, double T, double step) {
    AD *ad = new FADBAD_AD(2, NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver = new VNODE(ad);

    Solver->setFirstEntry();
    NLKG_init_l2 init_val = NLKG_init_approx_no_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    int iter = 0;

    cout << "Solution data: " <<
            "b = [" << inf(b) << "," << sup(b) << "], " << 
            "T = " << T << endl;
    cout << "iter,time,inf(y[0]),sup(y[0]),inf(y[1]),sup(y[1])" << endl;

    while (sup(t) < T) {
        // log data
        cout << setprecision(5) << iter << "," << 
                midpoint(t) << "," << 
                inf(y[0]) << "," << sup(y[0]) << "," << 
                inf(y[1]) << "," << sup(y[1]) << endl;


        interval new_t = t + step;
        Solver->integrate(t, y, new_t);

        if(!Solver->successful()) {
            cout << "Failed to integrate, stopping.";
            return;
        }

        iter++;
    }
}

/**
Wrapper function for proving uniqueness 
Takes `N', number of excited states to prove uniqueness for, as a command line
argument.
*/
int run_uniqueness_prover(int argc, char *argv[]) {
    int N = 0;

    if (argc < 2) {
        cerr << "Usage: \n" << argv[0] << " N" << endl <<
                "N is number of bound states to prove uniqueness for" << endl;
        return 0;
    }
    else {
        istringstream ss(argv[1]);
        if (!(ss >> N)) {
            cerr << "Invalid number: " << argv[1] << endl;
        }
        else if (!ss.eof()) {
            cerr << "Trailing characters after number: " << argv[1] << endl;
        }
        else if (N < 0) {
            cerr << "Number of bound states must be nonnegative" << endl;
        }
    }

    clock_t t_start = clock();
    cout << "Making plan to prove first " << N << 
            " bound states are unique." << endl;
    vector<interval_handler> plan = make_n_first_excited_plan(N);
    cout << "\nExecuting plan to prove first " << N << 
            " bound states are unique." << endl;
    execute_first_n_excited_plan(N, plan);

    clock_t t_end = clock();
    cout << "Time to execute = " << 
            (double)(t_end - t_start) / CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

/**
Wrapper function for making output data for graph of first three excited
state uniqueness.
*/
int make_N3_output_for_graphs() {
    double step = 0.001;
    // first bound state
    make_sol_data(interval(4.26611, 4.43298), 1.92140, step); 
    cout << endl;

    // second bound state
    make_sol_data(interval(14.09452, 14.11536), 2.85488, step); 
    cout << endl;

    // third bound state
    make_sol_data(interval(29.09019, 29.17356), 4.97016, step); 
    cout << endl;

    // fourth bound state
    make_sol_data(interval(49.33922, 49.38089), 5.90799, step); 
    cout << endl;

    return 0;
}


/**
Create some data to compare finitary equation to infinity equation
Should have y_b(t) = bw(bt) or beta y_b(beta t) = w(t)
*/
int make_output_dat_comp_finite_inf() {
    interval b(25.0, 25.00001);
    interval beta = 1/b;


    interval step(0.001);
    int num_steps = 10000;

    // integrator for finitary eq
    AD *ad_fin = new FADBAD_AD(2, NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver_fin = new VNODE(ad_fin);

    // integrator for infinitary eq
    AD *ad_inf = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &beta);
    VNODE *Solver_inf = new VNODE(ad_inf);

    NLKG_init_l2 init_val_fin = NLKG_init_approx_no_delta(b, INIT_WIDTH);
    iVector y_fin = init_val_fin.y;
    interval t_fin = init_val_fin.t0;

    NLKG_init_l2 init_val_inf = NLKG_infty_init_approx(beta, INIT_WIDTH);
    iVector y_inf = init_val_inf.y;
    interval t_inf = init_val_inf.t0;

    // cout << t_inf << ", " << y_inf[0] << endl;

    // integrate equations
    cout << "Integrating equations" << endl;
    cout << "step,time,inf(y[0]),sup(y[0]),b*time,inf(b*w(bt)),sup(b*w(bt))" << 
            endl;
    for (int i = 0; i < num_steps; i++) {
        // cout << "Energy: " << energy(y_fin) << endl;
        interval new_t_fin = t_fin + step;
        interval new_t_inf = b * new_t_fin;
        
        Solver_fin->integrate(t_fin, y_fin, new_t_fin);
        Solver_inf->integrate(t_inf, y_inf, new_t_inf);

        if(!Solver_fin->successful() || !Solver_inf->successful()) {
            cout << "Failed to integrate, stopping.";
            return 1;
        }

        if(!does_intersect(y_fin[0], b * y_inf[0])) {
            cout << "DOES NOT INTERSECT" << endl;
            return 1;
        }

        cout << i << "," << 
                midpoint(t_fin) << "," << inf(y_fin[0]) << "," << 
                                          sup(y_fin[0]) << "," << 
                midpoint(t_inf) << "," << inf(b * y_inf[0]) << "," << 
                                          sup(b * y_inf[0]) << endl;
    }

    return 0;
}

int main(int argc, char *argv[]) {
    return run_uniqueness_prover(argc, argv);
}