#include <METHODS/usermethod.h>
#include <METHODS/adam.h>

UserMethod::UserMethod()
{
    addParam(Parameter("popsize", "20", "Μέγεθος πληθυσμού φωλιών."));
    addParam(Parameter("max_iters", "1000", "Μέγιστος αριθμός επαναλήψεων."));
    addParam(Parameter("seed", "42", "Τιμή σποράς για την τυχαία γεννήτρια."));
    addParam(Parameter("Pa", "0.25", "Πιθανότητα εγκατάλειψης φωλιάς [0.0, 1.0]."));
    addParam(Parameter("tol", "1e-6", "Ανοχή στην αλλαγή του best fitness."));
    addParam(Parameter("max_stall", "100", "Μέγιστος αριθμός επαναλήψεων χωρίς βελτίωση."));
    addParam(Parameter("top_k", "5", "Αριθμός καλύτερων φωλιών για τοπική αναζήτηση."));
    addParam(Parameter("p_local", "0.5", "Πιθανότητα χρήσης τοπικής αναζήτησης."));
    addParam(Parameter("local_every_R", "1", "Συχνότητα εκτέλεσης τοπικής αναζήτησης."));
}

void    UserMethod::init()
{
    srand((unsigned) time(NULL));
    defaultParametersInit();
    varInit();
    initNests();
}

Data UserMethod::levy_search(Data& current_nest, double beta) {
    // Ανάκτηση ορίων
    Data left_bounds = myProblem->getLeftMargin();
    Data right_bounds = myProblem->getRightMargin();
    double sigma_num = std::tgamma(1.0 + beta) * std::sin(M_PI * beta / 2.0);
    double sigma_den = std::tgamma((1.0 + beta) / 2.0) * beta * std::pow(2.0, (beta - 1.0) / 2.0);
    double sigma = std::pow(sigma_num / sigma_den, 1.0 / beta);

    Data new_nest = current_nest;
    for (int d = 0; d < dim; ++d) {
        double u = (double)(rand() % (int)(sigma+1));
        double v = (double)rand()/(double)RAND_MAX;
        if (std::abs(v) < 1e-9) v = 1e-9;

        double L = u / std::pow(std::abs(v), 1.0 / beta);
        double S = 0.01 * L * (current_nest[d] - var_best_nest[d]);
        new_nest[d] = current_nest[d] + S;
        new_nest[d] = max(left_bounds[d], min(right_bounds[d], new_nest[d]));
    }
    return new_nest;
}

void UserMethod::local_search(int index_to_replace) {
    // Ανάκτηση ορίων
    Data left_bounds = myProblem->getLeftMargin();
    Data right_bounds = myProblem->getRightMargin();

    //Εύρεση των k καλύτερων
    vector<pair<double, int>> fitness_indices(var_num_nests);
    for (int i = 0; i < var_num_nests; ++i){
        fitness_indices[i] = {var_fitness[i], i};
    }
    sort(fitness_indices.begin(), fitness_indices.end());
    int k = min(var_top_k, var_num_nests);
    if (k < 1) k = 1;

    int r1 = rand() % k;
    int idx1 = fitness_indices[r1].second;
    int r2 = rand() % k;
    int idx2 = fitness_indices[r2].second;

    while (idx1 == idx2 && k > 1) {
        idx2 = fitness_indices[(rand() % k)].second;
    }
    //Δημιουργία νέας λύσης (π.χ. με διαφορική εξέλιξη μεταξύ των top-k)
    Data new_nest(dim);
    //βρίσκει τυχαίο αριθμό στο διάστημα [-1,1]
    for (int d = 0; d < dim; ++d) {
        double r = ((double)rand() / (double)RAND_MAX);
        if (rand()%2){
            r = -r;
        }
        double step_size = r;
        new_nest[d] = var_nests[idx1][d] + step_size * (var_nests[idx1][d] - var_nests[idx2][d]);
        new_nest[d] = max(left_bounds[d], min(right_bounds[d], new_nest[d]));
    }
    //Αξιολόγηση
    double new_fit = myProblem->statFunmin(new_nest);
    if (new_fit < var_fitness[index_to_replace]) {
        var_nests[index_to_replace] = new_nest;
        var_fitness[index_to_replace] = new_fit;
        if (new_fit < var_best_fitness) {
            var_best_fitness = new_fit;
            var_best_nest = new_nest;
        }
    }
}

void UserMethod::nest_abandonment() {
    // Ανάκτηση ορίων
    Data left_bounds = myProblem->getLeftMargin();
    Data right_bounds = myProblem->getRightMargin();
    for (int i = 0; i < var_num_nests; ++i) {
        double r = (double)rand() / (double)RAND_MAX;
        if (r < var_Pa) {
            // Φωλιά εγκαταλείφθηκε
            double r2 = (double)rand() / (double)RAND_MAX;
            if (r2 < var_p_local && (var_current_iteration % var_local_every_R == 0)) {
                local_search(i);
            }
            else {
                //τυχαία αντικατάσταση
                for (int d = 0; d < dim; ++d) {
                    double range = right_bounds[d] - left_bounds[d];
                    double rnd = (double)rand() / (double)RAND_MAX;
                    var_nests[i][d] = left_bounds[d] + rnd * range;
                }
                var_fitness[i] = myProblem->statFunmin(var_nests[i]);
                if (var_fitness[i] < var_best_fitness) {
                    var_best_fitness = var_fitness[i];
                    var_best_nest = var_nests[i];
                }
            }
        }
    }
}

void    UserMethod::step()
{
    if (!terminated()){
        double prev_best_fitness = var_best_fitness;
        // Δημιουργία νέας λύσης
        int i = rand() % (int)var_num_nests;
        Data new_nest_pos = levy_search(var_nests[i], 1.5);
        // Αξιολόγηση
        double new_nest_fit = myProblem->statFunmin(new_nest_pos);
        //Σύγκριση
        int j = rand() % (var_num_nests);
        if (new_nest_fit < var_fitness[j]) {
            var_nests[j] = new_nest_pos;
            var_fitness[j] = new_nest_fit;
            if (new_nest_fit < var_best_fitness) {
                var_best_fitness = new_nest_fit;
                var_best_nest = new_nest_pos;
            }
        }
        //Εγκατάλειψη φωλιών
        nest_abandonment();

        //Έλεγχος στασιμότητας
        if (abs(prev_best_fitness - var_best_fitness) < var_tol) {
            var_stall_counter++;
        } else {
            var_stall_counter = 0;
        }
        var_current_iteration++;
    }
}

bool    UserMethod::terminated()
{
    return ((var_current_iteration == var_max_iter) || (var_stall_counter >= var_max_stall));
}

void    UserMethod::done()
{

    //Έλεγχος αν η λύση του CS είναι καλύτερη
    if (var_best_fitness < myProblem->getBesty()) {
        //Ενημέρωση του Problem με τη νέα καλύτερη τιμή
        myProblem->setKnownOptimum(var_best_fitness, var_best_nest);
    }
    cout << "Λύση με Cuckoo Search:" << var_best_fitness << endl;

    //Εφαρμογή Adam optimizer
    Data cuckoo_best_nest = var_best_nest;
    double cuckoo_best_fitness = var_best_fitness;
    Adam *AdamOptimizer = new Adam();
    AdamOptimizer->setProblem(myProblem);
    AdamOptimizer->setPoint(cuckoo_best_nest, cuckoo_best_fitness);
    AdamOptimizer->solve();

    Data local_best_nest;
    double local_best_fitness;
    AdamOptimizer->getPoint(local_best_nest, local_best_fitness);
    //αν προέκυψε καλύτερη λύση
    if (local_best_fitness < myProblem->getBesty()) {
        //ενημέρωση της λύσης
        myProblem->setKnownOptimum(local_best_fitness, local_best_nest);
        var_best_nest = local_best_nest;
        var_best_fitness = local_best_fitness;
        cout << "Λύση με Adam Optimizer:" << var_best_fitness << endl;
    }
    else{
        cout << "Δεν βρέθηκε καλύτερη λύση με Adam Optimizer!" << endl;
    }
    delete AdamOptimizer;
}

void UserMethod::defaultParametersInit(){
    var_num_nests = getParam("popsize").getValue().toInt();
    var_max_iter = getParam("max_iters").getValue().toInt();
    var_Pa = getParam("Pa").getValue().toDouble();
    var_tol = getParam("tol").getValue().toDouble();
    var_max_stall = getParam("max_stall").getValue().toInt();
    var_top_k = getParam("top_k").getValue().toInt();
    var_p_local = getParam("p_local").getValue().toDouble();
    var_local_every_R = getParam("local_every_R").getValue().toInt();
}

void UserMethod::varInit(){
    var_current_iteration = 0;
    var_stall_counter = 0;
    var_best_fitness = std::numeric_limits<double>::max();
    dim = myProblem->getDimension();
    var_nests.resize(var_num_nests, Data(dim));
    var_fitness.resize(var_num_nests);
    var_best_nest.resize(dim);
}

void UserMethod::initNests(){
    Data left_bounds = myProblem->getLeftMargin();
    Data right_bounds = myProblem->getRightMargin();
    for (int i = 0; i < var_num_nests; ++i) {
        for (int d = 0; d < dim; ++d) {
            double range = right_bounds[d] - left_bounds[d];
            double r = (double)rand() / (double)RAND_MAX;
            var_nests[i][d] = left_bounds[d] + r * range;
        }
        // Αξιολόγηση
        var_fitness[i] = myProblem->statFunmin(var_nests[i]);

        if (var_fitness[i] < var_best_fitness) {
            var_best_fitness = var_fitness[i];
            var_best_nest = var_nests[i];
        }
    }
}

UserMethod::~UserMethod()
{

}
