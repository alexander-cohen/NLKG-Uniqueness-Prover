#include <ostream>
#include <fstream>
#include <math.h>
#include "vnode.h"
using namespace vnodelp;

double determine_min_height_for_crossings(int n);

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
- put NLKG definition and initialization into a seperate file
- make ODE work at zero, because not smooth :(
- make graphs comparing boundary ODE at infinity to finitary ODE
- make ``instruction file'' for checking uniqueness (seperate procedure from verification)
*/

/*
Options for initializing the algorithm:
- Taylor series
- Low order taylor expansion with bounds using Banach contraction
- The above but with automatically computed error bounds
- Approximation with linear equation & an priori bound
*/

// HELPER FUNCTIONS //

#include <vector>

struct NLKG_init {
    iVector y;
    interval t0;
    NLKG_init() : y(4) {}
};

struct NLKG_infty_init {
    iVector y;
    interval t0;
    NLKG_infty_init() : y(2) {}
};

const int ODE_POW = 3;
const interval DIM(3);

const double MAX_TIME_DET_CROSSINGS = 50; 
const double STEP_TIME = 0.1;
const double INIT_WIDTH = 1e-8;

double min(double a, double b) {
    if (a < b) return a;
    else return b;
}

const double MACHINE_EP = 1e-9;
interval lower_half(interval I) {
    return interval(inf(I), midpoint(I)+MACHINE_EP);
}

interval upper_half(interval I) {
    return interval(midpoint(I)-MACHINE_EP, sup(I));
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
    yp[3] = - ((DIM-1)/t) * y[3] - fderiv(y[0]) * y[2];
}

// check if we can guarantee delta > 0 up to time t
// see sec. 3.1
double min_time_valid(interval b) {
    static interval FDERIV_ZERO = sqrt(interval(3)) / interval(3);
    double m1 = inf(sqrt(interval(6) * (b - FDERIV_ZERO) / my_pow(b, 3)));
    double m2 = inf(log(interval(2)) / (interval(3) * my_pow(b, 2)));
    return min(m1, m2);
}

// check maximum errors at time t
// see sec. 3.1
void get_errors_for_time(interval b, interval t, vector<double> &errs) {
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
NLKG_init NLKG_init_approx(interval b, double tol) {
    NLKG_init vals;
    interval t0 = interval(min(min_time_valid(b), 0.1));
    double cur_err = 1e6;
    vector<double> errs(4);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time(b, t0, errs);
        double max_err = 0;
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


void get_errors_for_time_infty(interval t, vector<double> &errs) {
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
NLKG_infty_init NLKG_infty_init_approx(interval s, double tol) {
    NLKG_infty_init vals;
    interval t0 = 1e-2;
    double cur_err = 1e6;
    vector<double> errs(2);
    while (cur_err > tol) {
        t0 = t0 / interval(2);
        get_errors_for_time_infty(t0, errs);
        double max_err = 0;
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

// UNIQUENESS PROVING CODE // 

// sample NLKG solver
void sample_solve(interval b, interval tend) {
    // instantiate ODE solver
    AD *ad = new FADBAD_AD(4,NLKG,NLKG);
    VNODE *Solver = new VNODE(ad);

    Solver->setFirstEntry();
    NLKG_init init_val = NLKG_init_approx(b, INIT_WIDTH);

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

void NLKG_init_test() {
    for (double b = 2.0; b < 100.0; b += 2.0) {
        NLKG_init vals = NLKG_init_approx(interval(b, b + 1e-2), INIT_WIDTH);
        cout << b << ", " << vals.t0 << endl;
        for (int i = 0; i < 4; i++) {
            cout << i << ": " << vals.y[i] << endl;
        }
    }
}


// not meant to be rigorous, just for finding where the excited states actually 
// are before we prove that's where they are
int crossing_number_smallb(interval b, AD *ad, VNODE *Solver) {
    Solver->setFirstEntry();
    NLKG_init init_val = NLKG_init_approx(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    int num_crossings = 0;
    bool pos = true;
    while (sup(t) < MAX_TIME_DET_CROSSINGS && sup(energy(y)) > 0) {
        interval new_t = t + STEP_TIME;
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
// returns a value slightly smaller than the Nth excited state
double find_nth_excited_state(int N, double tol, AD *ad, VNODE *Solver) {
    double upper_bound = determine_min_height_for_crossings(N+1);
    // cout << "num cross upper: " << 
    //     crossing_number_smallb(upper_bound, ad, Solver) << ", " 
    //     << N << endl;
    double lower_bound = 1.0;
    // cout << "UPPER BOUND: " << upper_bound << endl;
    while (upper_bound - lower_bound > tol) {
        double mid = (upper_bound + lower_bound) / 2;
        int nc_mid = crossing_number_smallb(mid, ad, Solver);
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

    return lower_bound;
}


int crossing_number_infty(interval s, int max) {
    NLKG_infty_init init_val = NLKG_infty_init_approx(s, INIT_WIDTH);
    AD *ad = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &s);
    VNODE *Solver = new VNODE(ad);

    interval t = init_val.t0;
    iVector y = init_val.y;

    int num_crossings = 0;
    bool pos = true;
    while (num_crossings < max && sup(energy_infty(y, s)) > 0) {
        interval new_t = t + STEP_TIME;
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

// determine the minimum value of y(0) = b after which we must have 
// at least n crossings
double determine_min_height_for_crossings(int n) {
    interval s(1);
    while (crossing_number_infty(s, n) < n) {
        s = s / interval(2);
    }
    // cout << "min height for " << n << ": " << s << endl;
    return sup(interval(1)/s);
}

// verify that within the interval s, there are at least `n' crossings 
// applies to equation at infinity
bool prove_crosses_many_infty(int n, interval s) {
    if (width(s) < MACHINE_EP * 4) { 
        return false; // if interval is too small we can't bisect further
    }
    int nc = crossing_number_infty(s, n);
    cout << s << nc << endl;
    if (nc == n) {
        return true;
    }
    else {
        bool lower_good = 
            prove_crosses_many_infty(n, lower_half(s));
        bool upper_good = 
            prove_crosses_many_infty(n, upper_half(s));
        return lower_good && upper_good;
    }
}

// prove that starting in the interval height b, the solution eventually falls 
// into one well
// treats b as a unified starting interval
const double EVENTUALLY_FALLS_STEP = 1.0;
bool prove_eventually_falls(interval b, AD *ad, VNODE *Solver) {
    Solver->setFirstEntry();
    NLKG_init init_val = NLKG_init_approx(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    while (sup(energy(y)) > 0) {
        interval new_t = t + EVENTUALLY_FALLS_STEP;
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
    if (width(ran_fall) < MACHINE_EP * 4) return false;

    // cout << "checking: " << ran_fall << endl;
    // try running directly on this interval first
    // if interval is too big this will fail, so we bisect
    bool res = prove_eventually_falls(ran_fall, ad, Solver);
    if (res) {
        // cout << ran_fall << endl;
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
bool bound_state_good(interval b, double tend, AD *ad, VNODE *Solver) {
    
    NLKG_init init_val = NLKG_init_approx(b, INIT_WIDTH);

    iVector y = init_val.y;
    interval t = init_val.t0;

    cout << y[0] << ", " << t << endl;

    Solver->setHmin(0);
    Solver->setFirstEntry();
    Solver->integrate(t,y,tend);

    if (!Solver->successful()) {
        cout << "VNODE-LP did not succesfully integrate" << endl;
        return false;
    }

    cout << b << ", " << tend << ":" << endl;
    for(int i = 0; i < 4; i++) {
        cout << "y" << i << " : " << y[i] << endl;
    }

    // all variables must have definite signs 
    for(int i = 0; i < 4; i++) {
        if(!disjoint(y[i], 0)) return false;
    }

    // this is valid because we already checked they have definite signs
    if(midpoint(y[2]) * midpoint(y[1]) < 0) return false;
    if(midpoint(y[3]) * midpoint(y[1]) < 0) return false;


    double yupb = sup(mag(y[0]));
    if(yupb > 0.5) return false;


    // energy bound
    double E = sup(energy(y));
    if (E > 0.25) {
        return false;
    }
    double T = sup(tend);
    double falls_over_v = T * E - 2 * E * log(T) + 0.82 * E;
    if (falls_over_v >= 0.3) {
        return false;
    }

    return true;
}

// verifies that the first n excited states, inclusive, are unique
bool first_n_excited_states_unique(int n, AD *ad, VNODE *Solver) {
    double b_checked_upto = 1.0 - MACHINE_EP;
    double bound_state_tol= 1e-6;
    double tend = 8;
    for(int i = 0; i <= n; i++) {
        double nth_excited_approx = find_nth_excited_state(i, bound_state_tol, ad, Solver);

        interval fall_int = interval(b_checked_upto, nth_excited_approx);
        interval bound_int = interval(nth_excited_approx - MACHINE_EP, nth_excited_approx + bound_state_tol);
        double next_b = nth_excited_approx + bound_state_tol - MACHINE_EP;

        bool falls_good = prove_eventually_falls_bisection(fall_int, ad, Solver);
        bool unique_near_bound = bound_state_good(bound_int, tend, ad, Solver);

        if (!falls_good) {
            cout << "could not prove falls in range: " << fall_int << endl;
            return false;
        }
        else {
            cout << "proved falls in range: " << fall_int << endl;
        }
        if (!unique_near_bound) {
            cout << "could not prove unique near bound state: " << bound_int << endl;
            return false;
        }
        else {
            cout << "proved unique near bound state: " << bound_int << endl;
        }

        interval fall_int_buff = interval(next_b, next_b + 2);
        bool falls_good_bugg = prove_eventually_falls_bisection(fall_int_buff, ad, Solver);
        if (!falls_good_bugg) {
            cout << "could not prove falls in buffer zone: " << fall_int_buff << endl;
            return false;
        }
        else {
            cout << "proved falls in buffer zone: " << fall_int_buff << endl;
        }
        b_checked_upto = sup(fall_int_buff) - 0.1;
    }

    double smax = sup(1/pow(interval(b_checked_upto), 2));
    // cout << b_checked_upto << ", " << smax << endl;
    interval s_check = interval(0, smax);
    cout << "s check: " << s_check << endl;
    bool crosses_many_infty = prove_crosses_many_infty(n+1, s_check);
    if (!crosses_many_infty) {
        cout << "could not prove crosses many times at infinity" << endl;
        return false;
    }
    else {
        cout << "proved crosses many times at infinity" << endl;
    }

    return true;
}

int main() {
    AD *ad = new FADBAD_AD(4,NLKG,NLKG);
    VNODE *Solver = new VNODE(ad);
    
    first_n_excited_states_unique(2, ad, Solver);
    // prove_crosses_many_infty(1, interval(0.0, 0.1));
    // cout << crossing_number_infty(interval(0.05), 2) << endl;
    // sample_solve(interval(4.337, 4.338), interval(10.0));
    return 0;
}