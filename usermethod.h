#ifndef USERMETHOD_H
#define USERMETHOD_H
# include <OPTIMUS/optimizer.h>
#include <climits>
#include<cstdlib>

using namespace std;
/**
 * @brief The UserMethod class is the skeleton class for the used defined optimization methods. The user should extend this class.
 */
class UserMethod :public Optimizer
{
private:
    int random;
    int paramValue;
    int var_num_nests;       // Ο συνολικός αριθμός των πιθανών λύσεων που κρατούνται σε κάθε επανάληψη
    int var_max_iter;        // Ο μέγιστος αριθμός βημάτων που θα εκτελεστούν
    double var_Pa;           // Η πιθανότητα με την οποία μια λύση θα αντικατασταθεί από άλλη
    double var_tol;          // Η ανοχή για τη βέλτιστη τιμή fitness
    int var_max_stall;       // Ο μέγιστος αριθμός επαναλήψεων για μη βελτιωμένη λύση
    int var_top_k;           // Ο αριθμός των καλύτερων λύσεων που χρησιμοποιούνται για τη δημιουργία νέων λύσεων (τοπική αναζήτηση)
    double var_p_local;      // Η πιθανότητα να εφαρμοστεί τοπική αναζήτηση
    int var_local_every_R;   // Ο αριθμός των βημάτων στον οποίο εφαρμόζεται τοπική αναζήτηση
    int dim;
    int var_current_iteration;    // Ο αριθμός της τρέχουσας επανάληψης
    vector<vector<double>> var_nests;
    std::vector<double> var_fitness;
    vector<double> var_best_nest;
    double var_best_fitness;
    int var_stall_counter;


public:
    /**
     * @brief UserMethod
     */
    UserMethod();
    /**
     * @brief init
     */
    virtual void init();
    /**
     * @brief step
     */
    virtual void step();
    /**
     * @brief terminated
     * @return
     */
    virtual bool terminated();
    /**
     * @brief done
     */
    virtual void done();
    /**
     * @brief ~UserMethod
     */
    virtual ~UserMethod();

    /*
     * Υλοποίηση της συνάρτησης δημιουργίας και καθορισμός προκαθορισμένων τιμών για τις παραμέτρους
     */
    void defaultParametersInit();
    /*
     * Αρχικοποίηση
     */
    void varInit();
    /*
     * Αρχικοποίηση των λύσεων (φωλιών)
     */
    void initNests();


    //**************************************
    // ΜΕΘΟΔΟΙ ΓΙΑ CUCKOOSEARCH
    //**************************************

    /*
     * LEVY SEARCH
     */
    Data levy_search(Data& current_nest, double beta);
    /*
     * ΕΓΚΑΤΑΛΕΙΨΗ ΛΥΣΗΣ (ΦΩΛΙΑΣ)
     */
    void nest_abandonment();
    /*
     * ΤΟΠΙΚΗ ΑΝΑΖΗΤΗΣΗ
     */
    void local_search(int index_to_replace);

};

#endif // USERMETHOD_H
