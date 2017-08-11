typedef struct lkmat lkmat;
struct lkmat{
  double* data;
  lkmat* next;
  lkmat* prev;
};

lkmat* dbMarkovPowers(double*, int, int);

double* force(double*, double*, int, double*, int);

void cholInv(double*, int);

void oneStep(double*, double*, int, double*, int, double);
