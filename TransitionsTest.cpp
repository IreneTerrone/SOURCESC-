
//static double Flag = -1; 
/*static double param_a;
static double param_u;
static double param_v;
static double param_b;*/
#define MAX_SUM 0.99999999
double param_K = pow(10, 6);
static vector<double> param_a{1.009, 1.008, 1.009, 1.008};
static vector<double> param_b{1,1,1,1};
static vector<double> param_u{3.4*pow(10, -5),3.4*pow(10, -5),3.4*pow(10, -5),3.4*pow(10, -5)};
static vector<double> param_v{0.016,0.016,0.016,0.016};



//da rimettere quando si eseguirà in R ma per ora passo direttamente i valori
/*void read_constant(string fname, double& rate)
{
  ifstream f (fname);
  string line;
  if(f.is_open())
  {
    int i = 1;
    while (getline(f,line))
    {
      switch(i)
      {
      case 1:
        rate = stod(line);
        //cout << fname << ": " << line << "\t" << endl;
        break;
      }
      ++i;
    }
    f.close();
  }
  else
  {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

void init_data_structures()
{
  read_constant("./u_file", param_u);
  read_constant("./v_file", param_v);
  read_constant("./a_file", param_a);
  read_constant("./b_file", param_b);
  read_constant("./K_file", param_K);
  
  Flag = 1; 
}*/


vector<double> getP(double A, double B, double deltaBranch){

  if(A==0 && B ==0){
    return {};
  }

  double lambda = A - B;
  double alpha = (B*exp(lambda*deltaBranch)-B)/(A*exp(lambda*deltaBranch)-B);
  double beta=(A*exp(lambda*deltaBranch)-A)/(A*exp(lambda*deltaBranch)-B);


  vector<double> P(100, 0.0);
  //double P[nPlaces];
  //std::fill(P, P + 100, 0.0);
  double sum_p = alpha;
  P[0] = alpha;
  int p_max = 1;
  for(; p_max<100 && sum_p<MAX_SUM; p_max++){
    P[p_max] = (1-alpha)*(1-beta)*pow(beta, p_max-1);
    sum_p=sum_p+P[p_max];
  }
  auto last_element = P.begin() + p_max -1;
  //double last_element = P[p_max-1];
  vector<double> P_new(P.begin(), last_element);
  double acc = std::accumulate(P.begin(),last_element,0.0);
  P_new.push_back(1-acc);
  cout << 1-acc << "calcolato con l'accumulate " << endl;
  return P_new;

}

double getRate(double total_marking, vector<double> p){

  double rate = 0;

  if(!p.empty()){
    double* p_arr = &p[0];
    int size = p.size();
    unsigned int mult_op[size];
    gsl_ran_multinomial(rng, size, total_marking, p_arr, mult_op);
    double mult = -1;


    for(int i = 0; i<size; i++){
      rate+= mult_op[i] * mult;
      mult++;
    }
  }

  return rate; 

}



double getV(double *Value,
 const struct InfTr* Trans,
 int index_place,
 const int T)
{

  double x = Value[Trans[T].InPlaces[0].Id];
  
  double v = param_v[Trans[T].InPlaces[0].Id] * ( 1 - x / param_K );

  return v;
}


double PMutationFunction(double *Value,
 map <string,int>& NumTrans,
 map <string,int>& NumPlaces,
 const vector<string> & NameTrans,
 const struct InfTr* Trans,
 const int T,
 const double& deltaBranch)
{
  //if( Flag == -1)   init_data_structures();

  double total_marking = Value[Trans[T].InPlaces[0].Id];
  
  double v = param_v[Trans[T].InPlaces[0].Id] * ( 1 - total_marking / param_K );
  
  double rate = 0;

  vector<double> p = getP(v, 0,  deltaBranch);

  rate = getRate(total_marking, p);

  
  return rate;
}


double DMutationFunction(double *Value,
 map <string,int>& NumTrans,
 map <string,int>& NumPlaces,
 const vector<string> & NameTrans,
 const struct InfTr* Trans,
 const int T,
 const double& deltaBranch)
{
  //if( Flag == -1)   init_data_structures();

  double total_marking = Value[Trans[T].InPlaces[0].Id];
  
  double u = param_u[Trans[T].InPlaces[0].Id] * ( 1 - total_marking / param_K );
  
  double rate = 0;

  vector<double> p = getP(u, 0, deltaBranch);


  rate = getRate(total_marking, p);
  return rate;
}


double getU(double *Value,
 const struct InfTr* Trans,
 int index_place,
 const int T)
{
  //if( Flag == -1)   init_data_structures();

  double total_marking = Value[Trans[T].InPlaces[0].Id];
  
  double u = param_u[Trans[T].InPlaces[0].Id] * ( 1 - total_marking / param_K );
  
  return u;
}




double GrowthFunction(double *Value,
 map <string,int>& NumTrans,
 map <string,int>& NumPlaces,
 const vector<string> & NameTrans,
 const struct InfTr* Trans,
 const int T,
 const double& deltaBranch)
{
  //if( Flag == -1)   init_data_structures();


  double total_marking = Value[Trans[T].InPlaces[0].Id];

  
  double u = getU(Value,Trans,Trans[T].InPlaces[0].Id, T);
  double v = getV(Value,Trans,Trans[T].InPlaces[0].Id, T);
  double xx = (1 - total_marking / param_K )/ 2;
  double uv = (1 - u - v)/ 2;

  
  double a = uv + (param_a[Trans[T].InPlaces[0].Id] - param_b[Trans[T].InPlaces[0].Id]) * xx;
  double b = uv - (param_a[Trans[T].InPlaces[0].Id] - param_b[Trans[T].InPlaces[0].Id]) * xx;
  
  vector<double> p = getP(a, b, deltaBranch);

  double rate = getRate(total_marking, p);


  return rate;


}



