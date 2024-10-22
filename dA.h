class dA {

public:
  int N;
  int n_visible;
  int n_hidden;
  double **W;
  double *hbias;
  double *vbias;
  double neglik;
  dA(int, int, int , double**, double*, double*);
  ~dA();
  void get_corrupted_input(int*, int*, double);
  void get_corrupted_input(double*, double*, double);

  void get_hidden_values(int*, double*);
  void get_hidden_values(double*, double*);

  void get_reconstructed_input(double*, double*);
  void train(int*, double, double);
  void train(double*, double, double);

  void reconstruct(int*, double*);
  void reconstruct(double*, double*);

};
