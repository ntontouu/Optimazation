#ifndef USERPROBLEM_H
#define USERPROBLEM_H

#include <OPTIMUS/problem.h>
#include <cmath>
#include <vector>

using namespace std;

class UserProblem : public Problem
{
public:
    UserProblem();
    virtual ~UserProblem();

    // Η συνάρτηση που υπολογίζει την Salomon
    virtual double funmin(Data &x);

    // Η συνάρτηση που υπολογίζει την κλίση
    virtual Data gradient(Data &x);
};

#endif
