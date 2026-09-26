#include <PROBLEMS/userproblem.h>

UserProblem::UserProblem()
{
    //Ορσιμός του προβλήματος σε 4 διαστάσεις
    setDimension(4);

    //Προσωρινά διανύσματα για τα όρια
    Data left(4);
    Data right(4);

    //Όρια του διανυσμάτος κλειστό -100,100
    for (int i = 0; i < 4; ++i) {
        left[i] = -100.0;
        right[i] = 100.0;
    }

    //Πέρασμα των οριών στην μνήμη του προβλήματος
    setLeftMargin(left);
    setRightMargin(right);
}

UserProblem::~UserProblem()
{
}

//Συνάρτηση Salomon
double UserProblem::funmin(Data &x)
{
    double sum_sq = 0.0;
    // Υπολογισμός του αθροίσματος τετραγώνων
    for (int i = 0; i < getDimension(); ++i) {
        sum_sq += x[i] * x[i];
    }

    // Υπολογισμός της ρίζαςτ του Σ
    double norm = std::sqrt(sum_sq);

    //Τύπος Salomon: 1 - cos(2*pi*norm) + 0.1*norm
    return 1.0 - std::cos(2.0 * M_PI * norm) + (0.1 * norm);
}

//Συνάρτηση Gradient - παράγωγος
Data UserProblem::gradient(Data &x)
{
    //Διάνυσμα που γυρνάει
    Data g(getDimension());

    double sum_sq = 0.0;
    for (int i = 0; i < getDimension(); ++i) {
        sum_sq += x[i] * x[i];
    }
    double norm = std::sqrt(sum_sq);
    if (norm < 1e-9) {
        for (int i = 0; i < getDimension(); ++i) {
            g[i] = 0.0;
        }
        return g;
    }

    // Υπολογισμός μερικής παραγώγου  - chain rule
    //df/du = 2*pi*sin(2*pi*u) + 0.1
    double df_du = (2.0 * M_PI * std::sin(2.0 * M_PI * norm)) + 0.1;

    // Τελικός υπολογισμός για κάθε διάσταση 1,2,3,4
    for (int i = 0; i < getDimension(); ++i) {
        g[i] = df_du * (x[i] / norm);
    }

    return g;
}
