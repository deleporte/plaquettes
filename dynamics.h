typedef struct lkmat lkmat;
struct lkmat{
  double* data;
  lkmat* next;
  lkmat* prev;
};

lkmat* dbMarkovPowers(double*, int, int);
