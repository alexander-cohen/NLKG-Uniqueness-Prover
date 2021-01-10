#include <ostream>
#include <math.h>
#include "vnode.h"
#include "nlkg_helper.h"
using namespace vnodelp;

// numerical values of excited states:
interval b0(4.337, 4.338);
// b_1 =  

/*
TODO:
- put NLKG definition and initialization into a seperate file
- make ODE work at zero, because not smooth :(
*/

/*
Options for initializing the algorithm:
- Taylor series
- Low order taylor expansion with bounds using Banach contraction
- The above but with automatically computed error bounds
- Approximation with linear equation & an priori bound
*/

// sample NLKG solver
void sample_solve() {
    interval t = 0.001, tend = 0.5;
    iVector y(4); // y is a state vector, has 4 entries

    // initial vector y has y = [y(t), y'(t), delta(t), delta'(t)]
    y[0] = interval(4.0); //interval(4.337, 4.338); // interval(4.337387, 4.337389);
    y[1] = 0.0;
    y[2] = 1.0;
    y[3] = 0.0;


    // // alpha is the power in the equation y'' + (2/t)y' + y^{alpha} - y = 0
    // interval alpha = interval(3.0);

    // instantiate ODE solver
    AD *ad = new FADBAD_AD(4,NLKG,NLKG);
    VNODE *Solver = new VNODE(ad);

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
    for (float b = 2.0; b < 100.0; b += 2.0) {
        NLKG_init vals = NLKG_init_approx(interval(b, b + 1e-2), 1e-6);
        cout << b << ", " << vals.t0 << endl;
        for (int i = 0; i < 4; i++) {
            cout << i << ": " << vals.y[i] << endl;
        }
    }
}


int crossing_number(interval s, int max) {
    NLKG_infty_init init_val = NLKG_infty_init_approx(s, 1e-8);
    AD *ad = new FADBAD_AD(2, NLKG_inf, NLKG_inf, &s);
    VNODE *Solver = new VNODE(ad);

    interval t = init_val.t0;
    iVector y(2);
    y[0] = init_val.y[0];
    y[1] = init_val.y[1];

    float tstep = 1.0;

    int num_crossings = 0;
    bool pos = true;
    while (num_crossings < max && sup(energy_infty(y, s)) > 0) {
        interval new_t = t + tstep;
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
float determine_min_height_for_crossings(int n) {
    interval s(1);
    while (crossing_number(s, n) < n) {
        s = s / interval(2);
        cout << s << endl;
    }
    return sup(interval(1)/sqrt(s));
}

bool prove_crosses_many(int n, interval s, float minw) {
    if (width(s) < minw) {
        return false;
    }

    int nc = crossing_number(s, n);
    if (nc == n) {
        return true;
    }
    else {
        bool lower_good = 
            prove_crosses_many(n, interval(inf(s), midpoint(s)), minw);
        bool upper_good = 
            prove_crosses_many(n, interval(midpoint(s), sup(s)), minw);
        return lower_good && upper_good;
    }
}

// prove that starting in the interval height b, the solution eventually falls 
// into one well
bool prove_eventually_falls(interval b, double step, AD *ad, VNODE *Solver) {

}

// prove that starting in the interval height b, crosses heigh at least
// 'ncross' times
bool prove_crosses_many(interval b, double step, AD *ad, VNODE *Solver, 
    int ncross) {

}


int main() {
    // sample_solve();
    // NLKG_init_test();
    // cout << crossing_number(interval(0.001), 10) << endl;
    int n = 2;
    float b = determine_min_height_for_crossings(n);
    float upper_s = sup(1/pow(interval(b), 2));
    bool proven = prove_crosses_many(n, interval(0, upper_s), 1e-8);
    cout << "PROVEN RESULT: " << proven << endl;
    return 0;
}