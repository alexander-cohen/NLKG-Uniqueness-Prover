#include <ostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "vnode.h"

using namespace vnodelp;

double determine_min_height_for_crossings(int n, AD *ad, VNODE *Solver);

/*
excited states from python code:
    b0 =  4.3373877034
    b1 = 14.1035847287
    b2 = 29.1312123369
    b3 = 49.3607112672
    b4 = 74.7723053303
*/


// numerical values of excited states:
interval b0(4.337, 4.338);
// b_1 =  

/*
TODO:
- make graphs comparing boundary ODE at infinity to finitary ODE
- make ``instruction file'' for checking uniqueness (seperate procedure from verification)
- check *every* floating point operation to make sure it is rigorous
- add detailed logging
- don't use `delta' unless we are in that specific case, otherwise only use `y'
- maybe try and section off all the ODE code to a seperate file (like we were doing before)
    then just have three functions that take in initial conditions, 
    AD, Solver, and prepare it
- figure out how to read / write intervals from file

- make a list of which functions are used in rigorous vs non rigorous parts.
*/

// HELPER FUNCTIONS //

#include <vector>


// STRUCTS
struct NLKG_init_l4 {
    iVector y;
    interval t0;
    NLKG_init_l4() : y(4) {}
};

struct NLKG_init_l2 {
    iVector y;
    interval t0;
    NLKG_init_l2() : y(2) {}
};

// CONSTANTS
const int ODE_POW = 3;
const interval DIM(3);

const double MAX_TIME_DET_CROSSINGS = 50; // maximum time to check up to for crossings

const double STEP_TIME = 0.05; // step time in checking for number of crossings
const double EVENTUALLY_FALLS_STEP = 1.0; // we use longer step size for proving falling

// tolerance allowed in initializing the algorithm with taylor approx
const double INIT_WIDTH = 1e-8;

const double TIME_START_INFTY_EQ = 1e-2;
const double MIN_TIME_START = 1e-2;

const double MIN_WIDTH_EVENTUALLY_FALLS = 0.5; 
const double MIN_WIDTH_ALGO_RUN = 2.0; 

const double BOUND_STATE_TOL = 1e-5;

double min(double a, double b) {
    if (a < b) return a;
    else return b;
}

// used right now to make buffer for floating point operations
// ideally should not be used (rounding mode instead)
// however if it is used, right way to do it is verify after the whole procedure 
// has been run that all intervals intersect
const double EP = 1e-9; 


// HELPER FUNCTIONS
interval lower_half(interval I) {
    return interval(inf(I), midpoint(I)+EP);
}

interval upper_half(interval I) {
    return interval(midpoint(I)-EP, sup(I));
}

bool does_intersect(interval I, interval J) {
    interval _(0);
    return intersect(_, I, J);
}

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

bool interval_le(interval I, interval J) {
    return get_sign(J - I) == POS;
}

// returns an interval containining a lower bound on the absolute value of I
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
void NLKG_with_delta(int n, var_type *yp, const var_type *y, var_type t, void *param) {
    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - f(y[0]);
    yp[2] = y[3];
    yp[3] = - ((DIM-1)/t) * y[3] - fderiv(y[0]) * y[2];
}

/*
defines NLKG ODE vector field without delta

*/
template<typename var_type>
void NLKG_no_delta(int n, var_type *yp, const var_type *y, var_type t, void *param) {
    yp[0] = y[1];
    yp[1] = - ((DIM-1)/t) * y[1] - f(y[0]);
}

// check if we can guarantee delta > 0 up to time t
// see sec. 3.1
// TODO: verify that we can also guarantee all values move in the 
// right direction up to this time
double min_time_valid_delta(interval b) {
    static interval SQRT3 = sqrt(interval(3));
    double m1 = inf( sqrt( interval(6) * (SQRT3 * b - interval(1)) / 
                            (SQRT3 * b * (b*b - 1)) ) );
    double m2 = inf(log(interval(4)) / (SQRT3 * b));
    return min(m1, m2);
}

// check maximum errors at time t, no delta
// see sec. 3.1
void get_errors_for_time_no_delta(interval b, interval t, vector<interval> &errs) {
    errs[0] = my_pow(b, 5) * my_pow(t, 4) / interval(40);
    errs[1] = my_pow(b, 5) * my_pow(t, 3) / interval(10);
}

// check maximum errors at time t
// see sec. 3.1
void get_errors_for_time_with_delta(interval b, interval t, vector<interval> &errs) {
    get_errors_for_time_no_delta(b, t, errs);
    errs[2] = my_pow(b, 4) * my_pow(t, 4) / interval(8);
    errs[3] = my_pow(b, 4) * my_pow(t, 3) / interval(2);
}


void set_initial_vals_no_delta(interval b, interval t0, iVector &y) {
    y[0] = b - my_pow(t0, 2) * f(b) / interval(6); // y(t0)
    y[1] = -t0 * f(b) / interval(3); // y'(t0)
}

void set_initial_vals_with_delta(interval b, interval t0, iVector &y) {
    set_initial_vals_no_delta(b, t0, y);
    y[2] = 1 - my_pow(t0, 2) * fderiv(b) / interval(6); // delta(t0)
    y[3] = -t0 * fderiv(b) / interval(3); // delta'(t0)
}


// y''(0) = -f(y)/3
// y(t) ~ b - f(y) t^2/6
// y'(t) ~ -f(y)t/3
// delta(t) = 1 - f'(y) t^2/6 
// delta'(t) = -f'(y) t/3
NLKG_init_l4 NLKG_init_approx_with_delta(interval b, double tol) {
    NLKG_init_l4 vals;
    interval t0 = interval(min(min_time_valid_delta(b), MIN_TIME_START));
    double cur_err = 1e6;
    vector<interval> errs(4);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_with_delta(b, t0, errs);
        double max_err = 0;

        // check if errors are small enough
        // this is floating point arithmetic but that is ok, this doesn't need
        // to be rigorous, just qualitatively need the errors "small"
        // what's important is that the error tolerances computed at the end
        // of this function aare accurate. 
        for(int i = 0; i < 4; i++) {
            if (sup(errs[i]) > max_err) max_err = sup(errs[i]); 
        }

        cur_err = max_err;
    }

    vals.t0 = t0;

    // initial calculations based on taylor approx
    set_initial_vals_with_delta(b, t0, vals.y);

    // add in error bars manually
    // note: we know that in this range, the approximation is in the stated direction
    // i.e., for y* actual y, y[i] < y*[i] < y[i] + errs[i];
    for(int i = 0; i < 4; i++) {
        vals.y[i] = interval(inf(vals.y[i]), sup(vals.y[i] + errs[i]));
    }

    return vals;
}

NLKG_init_l2 NLKG_init_approx_no_delta(interval b, double tol) {
    NLKG_init_l2 vals;
    interval t0 = MIN_TIME_START;
    double cur_err = 1e6;
    vector<interval> errs(2);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_no_delta(b, t0, errs);
        double max_err = 0;

        // check if errors are small enough
        // this is floating point arithmetic but that is ok, this doesn't need
        // to be rigorous, just qualitatively need the errors "small"
        // what's important is that the error tolerances computed at the end
        // of this function aare accurate. 
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


void get_errors_for_time_infty(interval t, vector<interval> &errs) {
    errs[0] = my_pow(t, 4) / interval(40);
    errs[1] = my_pow(t, 3) / interval(10);
}

interval energy_infty(v_blas::iVector& y, interval s) {
    interval E = my_pow(y[1], 2) / 2 + my_pow(y[0], 4)/(interval(4)) - s * my_pow(y[0], 2) / 2;
    return E;
}

// w''(0) = -(1-s) / 3
// w(t) ~ 1 - (1-s)t^2/6
// w'(t) ~ -(1-s)t/3
NLKG_init_l2 NLKG_infty_init_approx(interval s, double tol) {
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
    vals.y[0] = interval(1) - (interval(1)-s) * my_pow(t0, 2) / interval(6); // v(t0)
    vals.y[1] = -t0 * (interval(1) - s) / interval(3); // y'(t0)

    // add in error bars manually
    // TODO: verify in this range, we have estimates in the stated directions
    vals.y[0] = interval(inf(vals.y[0]), sup(vals.y[0] + errs[0]));
    vals.y[1] = interval(inf(vals.y[1]), sup(vals.y[1] + errs[1]));

    return vals;
}

// UNIQUENESS PROVING CODE // 

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
        NLKG_init_l4 vals = NLKG_init_approx_with_delta(interval(b, b + 1e-2), INIT_WIDTH);
        cout << b << ", " << vals.t0 << endl;
        for (int i = 0; i < 4; i++) {
            cout << i << ": " << vals.y[i] << endl;
        }
    }
}

// determine an amount of time such that after integrating for this long
// we cross zero at most once.
// Is rigorous. 
interval time_max_one_cross(iVector y, interval t) {
    // if we travel distance to zero + 1 we are at most in the left well so 
    // won't cross zero more than once.  
    static interval MIN_ENERGY = -interval(1/2);
    interval min_d = min_mag(y[0]) + interval(1);
    interval E = interval(sup(energy(y))); // overestimate the energy, as this will underestimate the step
    interval maxv = sqrt(2 * (E - MIN_ENERGY));
    interval step = min_d / maxv;
    return interval(inf(step));
}

interval time_max_one_cross_infty(iVector y, interval t, interval s) {
    // if we travel distance to zero + 1 we are at most in the left well so 
    // won't cross zero more than once. 

    // energy minimized at y^3 - sy = 0, or y = sqrt(s);
    interval MIN_ENERGY = my_pow(s, 2) / interval(4) - s / interval(2);

    interval min_d = min_mag(y[0]) + interval(sup(sqrt(s)));
    interval E = interval(sup(energy_infty(y, s)));
    interval maxv = sqrt(2 * (E - MIN_ENERGY));
    interval step = min_d / maxv;
    return interval(inf(step));
}

/*
Not meant to be rigorous, just for finding where the excited states actually 
are before we prove that's where they are

Should be rigorous, however.
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
        interval step = time_max_one_cross(y, t);
        interval new_t = t + step;

        Solver->integrate(t, y, new_t);
        
        if(!Solver->successful()) {
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

// will return best estimate of Nth excited state
// performs a binary search
// returns an interval which contains the Nth excited state
// Question: check if this is rigorous
interval find_nth_excited_state(int N, double tol, AD *ad, VNODE *Solver) {
    double upper_bound = determine_min_height_for_crossings(N+1, ad, Solver);
    // cout << "num cross upper: " << 
    //     crossing_number_smallb(upper_bound, ad, Solver) << ", " 
    //     << N << endl;
    double lower_bound = 1.0;
    // cout << "UPPER BOUND: " << upper_bound << endl;
    while (upper_bound - lower_bound > tol) {
        double mid = (upper_bound + lower_bound) / 2;
        int nc_mid = crossing_number_smallb(mid, N+1, ad, Solver);
        // cout << "MID: " << mid << ", " << nc_mid << endl;
        if (nc_mid == -1) {
            cout << "failed to find " << N << "th excited state" << endl;
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


// determine a minimum value of y(0) = b after which we *should* have 
// at least n crossings
// does not prove that there are at least n crossings after this height
double determine_min_height_for_crossings(int n, AD *ad, VNODE *Solver) {
    double b = 1.0;
    while (crossing_number_smallb(interval(b), n, ad, Solver) < n)  {
        b *= 2;
    }
    return b;
}


// determined number of crossings using infinitary equation
// TODO: change step size to be determined dynamically based on energy
// TODO: graph and compare to finitary equations
// guarantees that it will cross *at least* as many times as the return value 
// could potentially cross more times
int crossing_number_infty(interval s, int max, AD *ad, VNODE *Solver) {
    ad->eval(&s);
    Solver->setFirstEntry();

    NLKG_init_l2 init_val = NLKG_infty_init_approx(s, INIT_WIDTH);

    interval t = init_val.t0;
    iVector y = init_val.y;

    int num_crossings = 0;
    bool pos = true;
    while (num_crossings < max && sup(energy_infty(y, s)) > 0) {
        // cout << y[0] << ", " << t << endl;
        if (width(y[0]) > MIN_WIDTH_ALGO_RUN) return false;

        interval new_t = t + time_max_one_cross_infty(y, t, s);
        Solver->integrate(t, y, new_t);
        
        if(!Solver->successful()) {
            // cout << "VNODE-LP could not reach t =  "<< new_t << endl;
            // cout << "Solution enclosure at t = " << t << endl;
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

    return num_crossings;
}

// verify that within the interval s, there are at least `n' crossings 
// applies to equation at infinity
bool prove_crosses_many_infty(int n, interval s, AD *ad, VNODE *Solver) {
    if (width(s) < EP * 4) { 
        cout << "Interval too small, cannot verify crosses many times at infinity: " << 
            s << endl;
        return false; // if interval is too small we can't bisect further
    }
    int nc = crossing_number_infty(s, n, ad, Solver);
    if (nc == n) {
        cout << "Proved crosses at least " << n << " times: " << s << endl;
        return true;
    }
    else {
        bool lower_good = 
            prove_crosses_many_infty(n, lower_half(s), ad, Solver);
        if (!lower_good) return false;
        bool upper_good = 
            prove_crosses_many_infty(n, upper_half(s), ad, Solver);
        if (!upper_good) return false;

        return true;
    }
}

// prove that starting in the interval height b, the solution eventually falls 
// into one well
// treats b as a unified starting interval
// TODO: determine step size dynamically
bool prove_eventually_falls(interval b, AD *ad, VNODE *Solver) {
    if (width(b) > MIN_WIDTH_EVENTUALLY_FALLS) {
        return false;
    }

    Solver->setFirstEntry();
    NLKG_init_l2 init_val = NLKG_init_approx_no_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    while (sup(energy(y)) > 0) {
        if (width(y[0]) > MIN_WIDTH_ALGO_RUN) return false;

        interval new_t = t + max(sup(time_max_one_cross(y, t)), 1.0);
        Solver->integrate(t, y, new_t);
        if(!Solver->successful()) {
            // cout << "VNODE-LP could not reach t = " << new_t << endl;
            return false;
        }
    }

    return true;
}

// bisects the starting interval b into smaller pieces until we can prove 
bool prove_eventually_falls_bisection(interval ran_fall, AD *ad, VNODE *Solver) {
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
    bool lower_good = prove_eventually_falls_bisection(lower_half(ran_fall), 
            ad, Solver);
    bool upper_good = prove_eventually_falls_bisection(upper_half(ran_fall), 
            ad, Solver);  

    return lower_good && upper_good;
}

// verify that within the interval b, there is at most one bound state
// TODO: verify all floating point calculations are rigorous
// TODO: decide ending time dynamically? maybe not...
bool bound_state_good(interval b, AD *ad, VNODE *Solver) {
    
    NLKG_init_l4 init_val = NLKG_init_approx_with_delta(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    // cout << y[0] << ", " << t << endl;

    Solver->setHmin(0);
    Solver->setFirstEntry();

    while (true) {
        interval new_t = t + time_max_one_cross(y, t);
        Solver->integrate(t,y,new_t);
        interval E = energy(y);

        if (!Solver->successful()) {
            cout << "VNODE-LP did not succesfully integrate" << endl;
            return false;
        }

        // if E < 0, return false
        if (interval_le(E, interval(0))) {
            return false;
        }

        // all variables must have definite signs
        for(int i = 0; i < 4; i++) {
            if(get_sign(y[i]) == UNC) continue;
        }

        // this is valid because we already checked they have definite signs
        if(midpoint(y[2]) * midpoint(y[1]) < 0) continue;
        if(midpoint(y[3]) * midpoint(y[1]) < 0) continue;


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

        cout << "Succesfully proved at most one excited state in range: " << b << endl;
        cout << "Ending time = " << new_t << ", ending y = " << endl;
        for(int i = 0; i < 4; i++) {
            cout << "y" << i << " : " << y[i] << endl;
        }

        return true;
    }
}


// make plan for proving the first n excited states are unique
enum pf_cases{FALLS,BOUND_GOOD,CROSSES_MANY_INFTY};
struct interval_handler {
    int method;
    int excited_state_num;
    interval the_interval;
    interval_handler(int m, interval I): method(m), the_interval(I) {};
    interval_handler(int m, int en, interval I): 
        method(m), excited_state_num(en), the_interval(I) {};
};

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

// does not need to be rigorous, just makes a plan
// verification of the plan must be rigorous
vector<interval_handler> make_n_first_excited_plan(int n) {
    static interval ENERGY_ZER0 = sqrt(interval(2));

    vector<interval_handler> pf_plan;

    AD *ad = new FADBAD_AD(2,NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver = new VNODE(ad);

    double b_checked_upto = inf(ENERGY_ZER0);
    double buffer_size = 2.0;

    for(int i = 0; i <= n; i++) {
        cout << "On excited state #" << i << endl;
        interval nth_excited_approx = find_nth_excited_state(i, 
            BOUND_STATE_TOL, ad, Solver);

        interval drop_before = interval(b_checked_upto, 
            inf(nth_excited_approx) + EP);

        pf_plan.push_back(interval_handler(FALLS, drop_before));
        cout << "created interval to check: FALLS, " << drop_before << endl;

        pf_plan.push_back(interval_handler(BOUND_GOOD, i, nth_excited_approx));
        cout << "created interval to check: BOUND STATE GOOD, " << nth_excited_approx << endl;

        b_checked_upto = sup(nth_excited_approx - EP);
    }

    // check that it falls in a ``buffer interval'' above the last excited state
    interval buffer_int = interval(b_checked_upto - EP, b_checked_upto + buffer_size);
    pf_plan.push_back(interval_handler(FALLS, buffer_int));
    cout << "created interval to check: BUFFER FALLS, " << buffer_int << endl;
    b_checked_upto = sup(buffer_int) - EP; 

    // construct interval of s to check
    // need to check b_checked_upto --> infty
    // under inversion, this becomes 0 --> 1/b_checked_upto
    interval s_int = interval(0, sup(1/pow(interval(b_checked_upto), 2)) + 5 * EP);
    pf_plan.push_back(interval_handler(CROSSES_MANY_INFTY, s_int));
    cout << "created interval to check: INFTY CROSSES MANY, " << s_int << endl;

    return pf_plan;
}

// verify that in the interval I, all solutions eventually fall in an energy well
bool verify_drop(interval I) {
    cout << "Will verify drops: " << I << endl;
    AD *ad = new FADBAD_AD(2,NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver = new VNODE(ad);

    return prove_eventually_falls_bisection(I, ad, Solver);
}

// verify that in the interval I there is at most one excited state
bool verify_boundstate_good(interval I, int n) {
    cout << "Will verify #" << n << " bound state good: " << I << endl;

    AD *ad = new FADBAD_AD(2,NLKG_no_delta, NLKG_no_delta);
    VNODE *Solver = new VNODE(ad);
    bool crosses_right_below = crossing_number_smallb(interval(inf(I)), n+2, ad, Solver) == n;
    bool crosses_right_above = crossing_number_smallb(interval(sup(I)), n+2, ad, Solver) == n+1;
    if (crosses_right_below && crosses_right_above) {
        cout << "crosses " << n << " times below #" << n << " excited state, " << 
                "and " << n+1 << " times above." << endl;
    }
    else {
        cout << "Could not verify #" << n << "excited state is in range: " << I << endl;
        return false;
    }

    ad = new FADBAD_AD(4, NLKG_with_delta, NLKG_with_delta);
    Solver = new VNODE(ad);

    bool at_most_one = bound_state_good(I, ad, Solver);

 
    return at_most_one && crosses_right_below && crosses_right_above;
}

// verify that in the interval s we cross at least n+1 times
bool verify_crosses_many_infty(interval I, int n) {
    cout << "Will verify crosses many infinity: " << I << endl;

    AD *ad = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &I);
    VNODE *Solver = new VNODE(ad);

    return prove_crosses_many_infty(n+1, I, ad, Solver);
}

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
    interval outer_infty_interval = sqrt(1 / sub_infty_interval);
    if (!does_intersect(plan[plan_len - 2].the_interval, outer_infty_interval)) {
        cout << "Intervals not contiguous: end to infinity" << endl;
        cout << plan[plan_len - 2].the_interval << endl;
        cmut << outer_infty_interval << endl;
        return false;
    }


    cout << "Verified intervals are contiguous" << endl << endl;

    // verify that intervals do what they say
    for(uint i = 0; i < plan.size(); i++) {
        bool success = false;

        if (plan[i].method == FALLS) {
            success = verify_drop(plan[i].the_interval);
        }
        else if (plan[i].method == BOUND_GOOD) {
            success = verify_boundstate_good(plan[i].the_interval, plan[i].excited_state_num);
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

    cout << "\nSuccesfully proved first " << n << " excited states are unique." << endl;

    return true;
}

/*
verifies that the first n excited states, inclusive, are unique
TODO: make one function to prepare instructions for verifying uniqueness
     (finds times to go to and which intervals to check for what)
and another function to actually verify
the ``verify'' function should output "simple" intervals, i.e., 
as well as a summary of intervals that were bisected, also the input 
intervals that lead to a simple output. 
*/
// bool first_n_excited_states_unique(int n, AD *ad, VNODE *Solver) {
//     double b_checked_upto = 1.0 - EP;
//     double bound_state_tol= 1e-6;
//     double tend = 8;
//     for(int i = 0; i <= n; i++) {
//         interval nth_excited_approx = find_nth_excited_state(i, bound_state_tol, ad, Solver);

//         interval fall_int = interval(b_checked_upto, inf(nth_excited_approx+EP));
//         interval bound_int = nth_excited_approx;
//         double next_b = nth_excited_approx + bound_state_tol - EP;

//         bool falls_good = prove_eventually_falls_bisection(fall_int, ad, Solver);
//         bool unique_near_bound = bound_state_good(bound_int, tend, ad, Solver);

//         if (!falls_good) {
//             cout << "could not prove falls in range: " << fall_int << endl;
//             return false;
//         }
//         else {
//             cout << "proved falls in range: " << fall_int << endl;
//         }
//         if (!unique_near_bound) {
//             cout << "could not prove unique near bound state: " << bound_int << endl;
//             return false;
//         }
//         else {
//             cout << "proved unique near bound state: " << bound_int << endl;
//         }

//         interval fall_int_buff = interval(next_b, next_b + 2);
//         bool falls_good_bugg = prove_eventually_falls_bisection(fall_int_buff, ad, Solver);
//         if (!falls_good_bugg) {
//             cout << "could not prove falls in buffer zone: " << fall_int_buff << endl;
//             return false;
//         }
//         else {
//             cout << "proved falls in buffer zone: " << fall_int_buff << endl;
//         }
//         b_checked_upto = sup(fall_int_buff) - 0.1;
//     }

//     double smax = sup(1/pow(interval(b_checked_upto), 2));
//     interval s_check = interval(0, smax);
//     cout << "s check: " << s_check << endl;
//     bool crosses_many_infty = prove_crosses_many_infty(n+1, s_check);
//     if (!crosses_many_infty) {
//         cout << "could not prove crosses many times at infinity" << endl;
//         return false;
//     }
//     else {
//         cout << "proved crosses many times at infinity" << endl;
//     }

//     return true;
// }

int main() {
    // interval s = interval(0,0.01);
    // AD *ad = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &s);
    // VNODE *Solver = new VNODE(ad);

    // cout << "computing crossing number" << endl;
    // cout << crossing_number_infty(interval(0.1), 5, ad, Solver) << endl;

    // prove_crosses_many_infty(1, s, ad, Solver);
    // interval s = interval(0.0,1e-11);
    // AD *ad = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &s);
    // VNODE *Solver = new VNODE(ad);
    // cout << crossing_number_infty(s, 5, ad, Solver) << endl;

    clock_t t_start = clock();

    int N = 3;
    cout << "Making plan to prove first " << N << " excited states are unique." << endl;
    vector<interval_handler> plan = make_n_first_excited_plan(N);
    cout << "\nExecuting plan to prove first " << N << " excited states are unique." << endl;
    execute_first_n_excited_plan(N, plan);

    clock_t t_end = clock();
    cout << "Time to execute = " << (double)(t_end - t_start) / CLOCKS_PER_SEC << "s" << endl;
    return 0;
}