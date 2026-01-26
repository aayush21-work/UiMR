#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<algorithm>

int rows = 10;
int cols = 10;
constexpr int nt = 500;
int N = 2;

int idx(int i,int j) {
    return i*cols + j;
}

double* AllocateGrid() {
    return new double[rows*cols]{0.0};
}

void Apply_Boundary(double* A,double l,double r,double t,double b) {
    for(int i=0;i<rows;i++) {
        A[idx(i,0)]       = l;
        A[idx(i,cols-1)] = r;
    }
    for(int j=0;j<cols;j++) {
        A[idx(0,j)]       = t;
        A[idx(rows-1,j)] = b;
    }
}

void Solve_Laplacian(double* A) {
    double* Anew = AllocateGrid();
    const double tol = 1e-6;

    for(int k=0;k<nt;k++) {
        double err = 0.0;

        for(int i=1;i<rows-1;i++) {
            for(int j=1;j<cols-1;j++) {
                Anew[idx(i,j)] = 0.25 * (
                    A[idx(i+1,j)] + A[idx(i-1,j)] +
                    A[idx(i,j+1)] + A[idx(i,j-1)]
                );
                err = std::max(err, std::abs(Anew[idx(i,j)] - A[idx(i,j)]));
            }
        }

        for(int i=1;i<rows-1;i++)
            for(int j=1;j<cols-1;j++)
                A[idx(i,j)] = Anew[idx(i,j)];

        if(err < tol) break;
    }

    delete[] Anew;
}

void uimr_init(int N,double l,double r,double t,double b) {
    rows /= N;
    cols /= N;

    std::ofstream shape("shape.dat");
    shape << rows*cols << " " << 4 << "\n";

    std::ofstream init("init.cache");

    for(int i=1;i<=rows;i++) {
        for(int j=1;j<=cols;j++) {
            init << ((j==1)?l:0) << " "
                 << ((j==cols)?r:0) << " "
                 << ((i==rows)?t:0) << " "
                 << ((i==1)?b:0) << "\n";
        }
    }
}

void uimr_solve() {
    double le,ri,to,bo;
    std::string in_file="init.cache";
    std::string out_file="bounds.dat";
    std::string temp_file="temp.dat";

    const double tol=1e-5;
    int max_iter=1000;
    bool write_out=false;

    for(int k=1;k<=max_iter;k++) {
        std::ifstream ini(in_file);
        std::ofstream out(temp_file);
        bool converge=true;
        int id=0;

        while(ini >> le >> ri >> to >> bo) {
            double le_old=le, ri_old=ri, to_old=to, bo_old=bo;

            double* A = AllocateGrid();
            Apply_Boundary(A,le,ri,to,bo);
            Solve_Laplacian(A);

            double suml=0,sumr=0,sumt=0,sumb=0;

            for(int n=0;n<rows;n++) {
                suml += A[idx(n,1)];
                sumr += A[idx(n,cols-2)];
            }
            for(int n=0;n<cols;n++) {
                sumt += A[idx(1,n)];
                sumb += A[idx(rows-2,n)];
            }

            le = suml/rows;
            ri = sumr/rows;
            to = sumt/cols;
            bo = sumb/cols;

            out << le << " " << ri << " " << to << " " << bo << "\n";

            double diff = std::max({
                std::abs(le-le_old),
                std::abs(ri-ri_old),
                std::abs(to-to_old),
                std::abs(bo-bo_old)
            });

            if(diff > tol) converge=false;

            if(write_out) {
                std::ofstream sol("sol_"+std::to_string(id)+".dat");
                for(int i=0;i<rows;i++) {
                    for(int j=0;j<cols;j++)
                        sol << A[idx(i,j)] << " ";
                    sol << "\n";
                }
            }

            delete[] A;
            id++;
        }

        ini.close();
        out.close();
        std::swap(in_file,out_file);

        if(converge) {
            std::cout << "Converged after " << k << " iterations\n";
            write_out=true;
        }

        if(k==max_iter) {
            std::cout << "Maximum Iteration Reached\n";
            write_out=true;
        }

        if(write_out) break;
    }
}

int main() {
    uimr_init(N,100,50,50,80);
    uimr_solve();
    return 0;
}

