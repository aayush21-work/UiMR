#include<fstream>
#include<string>
#include<vector>
#include<omp.h>
#include<iostream>

int rows = 500;
int cols = 500;
int N    = 25;
const int nt = 10;
int sr, sc;

constexpr int MAX_OUTER = 5000;

inline int idx(int i,int j){ return i*sc+j; }

double* LV(std::vector<double>& f){ return f.data(); }
double* RV(std::vector<double>& f){ return f.data()+sr; }
double* TV(std::vector<double>& f){ return f.data()+2*sr; }
double* BV(std::vector<double>& f){ return f.data()+2*sr+sc; }

void Solve_Laplacian(double* __restrict__ A, std::vector<double>& face)
{
  for(int i=0;i<sr;i++){
    A[idx(i,0)]    = LV(face)[i];
    A[idx(i,sc-1)] = RV(face)[i];
  }
  for(int j=0;j<sc;j++){
    A[idx(0,j)]    = TV(face)[j];
    A[idx(sr-1,j)] = BV(face)[j];
  }
  for(int k=0;k<nt;k++)
    for(int i=1;i<sr-1;i++)
      for(int j=1;j<sc-1;j++)
        A[idx(i,j)]=0.25*(A[idx(i+1,j)]+A[idx(i-1,j)]+
                          A[idx(i,j+1)]+A[idx(i,j-1)]);
}

void exchange_boundaries(std::vector<std::vector<double>>& faces,
                         double l,double r,double t,double b)
{
  for(int ti=0;ti<N;ti++)
    for(int tj=0;tj<N-1;tj++){
      int a=ti*N+tj, b_id=ti*N+tj+1;
      double* ra=RV(faces[a]); double* lb=LV(faces[b_id]);
      for(int n=0;n<sr;n++){
        double tmp=(ra[n]+lb[n])*0.5;
        ra[n]=tmp; lb[n]=tmp;
      }
    }
  for(int ti=0;ti<N-1;ti++)
    for(int tj=0;tj<N;tj++){
      int a=ti*N+tj, b_id=(ti+1)*N+tj;
      double* ba=BV(faces[a]); double* tb=TV(faces[b_id]);
      for(int n=0;n<sc;n++){
        double tmp=(ba[n]+tb[n])*0.5;
        ba[n]=tmp; tb[n]=tmp;
      }
    }
  for(int id=0;id<N*N;id++){
    int bx=id%N, by=id/N;
    if(bx==0)   for(int n=0;n<sr;n++) LV(faces[id])[n]=l;
    if(bx==N-1) for(int n=0;n<sr;n++) RV(faces[id])[n]=r;
    if(by==0)   for(int n=0;n<sc;n++) TV(faces[id])[n]=t;
    if(by==N-1) for(int n=0;n<sc;n++) BV(faces[id])[n]=b;
  }
}

void uimr_init(double l,double r,double t,double b,
               std::vector<std::vector<double>>& faces)
{
  sr=rows/N; sc=cols/N;
  int fw=2*sr+2*sc;
  faces.assign(N*N, std::vector<double>(fw,0.0));
  for(int ti=0;ti<N;ti++)
    for(int tj=0;tj<N;tj++){
      int id=ti*N+tj;
      if(tj==0)   for(int n=0;n<sr;n++) LV(faces[id])[n]=l;
      if(tj==N-1) for(int n=0;n<sr;n++) RV(faces[id])[n]=r;
      if(ti==0)   for(int n=0;n<sc;n++) TV(faces[id])[n]=t;
      if(ti==N-1) for(int n=0;n<sc;n++) BV(faces[id])[n]=b;
    }
}

void uimr_solve(std::vector<std::vector<double>>& faces,
                double l,double r,double t,double b)
{
  int Nsub=N*N;
  // precompute d once
  int d, d2;

  std::vector<double*> bufs(Nsub,nullptr);
  for(int id=0;id<Nsub;id++)
    bufs[id] = new double[sr*sc]();

  // precompute after sr/sc are set
  d  = sr/2;
  d2 = sc-1-d;

  for(int k=0;k<MAX_OUTER;k++)
  {
    // static schedule: tiles are equal size so no need for dynamic
    #pragma omp parallel for schedule(static)
    for(int id=0;id<Nsub;id++)
    {
      Solve_Laplacian(bufs[id], faces[id]);
      double* buf = bufs[id];
      double* lv=LV(faces[id]), *rv=RV(faces[id]);
      double* tv=TV(faces[id]), *bv=BV(faces[id]);
      for(int n=0;n<sr;n++){
        lv[n] = buf[idx(n,d)];
        rv[n] = buf[idx(n,d2)];
      }
      for(int n=0;n<sc;n++){
        tv[n] = buf[idx(d,n)];
        bv[n] = buf[idx(sr-1-d,n)];
      }
    }
    exchange_boundaries(faces,l,r,t,b);
  }

  // single merged file with large write buffer
  std::ofstream sol("sol_uimr_omp.dat");
  sol.rdbuf()->pubsetbuf(nullptr, 1<<20);  // 1 MB write buffer
  for(int ti=0;ti<N;ti++)
    for(int i=0;i<sr;i++){
      for(int tj=0;tj<N;tj++){
        int id=ti*N+tj;
        double* buf=bufs[id];
        for(int j=0;j<sc;j++) sol<<buf[idx(i,j)]<<" ";
      }
      sol<<"\n";
    }

  for(int id=0;id<Nsub;id++) delete[] bufs[id];
}

int main()
{
  const double l=100,r=50,t=50,b=80;
  std::vector<std::vector<double>> faces;
  uimr_init(l,r,t,b,faces);

  uimr_solve(faces,l,r,t,b);
  return 0;
}
