#include<iostream>
#include<fstream>
#include<ostream>
#include<string>
#include<cmath>
#include<algorithm>
#include<vector>

int rows = 100;
int cols = 100;
constexpr int nt = 500;
int N = 10; 

int idx(int i,int j)
{
  return i*cols+j ;
}

int idx(int t, int i, int j)
{
  return t*rows*cols+i*cols+j;
}

double* Boundary_Initilise(int rows,int cols,double l,double r, double t, double b)
{
  double* A= new double [rows*cols];


  //Define Boundary Conditions
  for(int i=0;i<rows;i++) // Set the left and right boundaries
  {
    A[idx(i,0)]=l;
    A[idx(i,rows-1)]=r;
  }

  for(int j=0;j<cols;j++) // Set the top and bottom boundaries
  {
    A[idx(0,j)]=t;
    A[idx(cols-1,j)]=b;
  }
  return A;

}

  //for(int i=0;i<rows;i++)
  //{
    //for(int j=0;j<cols;j++)
    //{
      //out << history[idx(t,i,j)];
    //}
  //}

void Solve_Laplacian(double* A)
{

  //std::ofstream set("shape.dat");
  //set << nt <<" "<< rows << " " << cols ;
  //std::ofstream out("output.dat");
  //double* history = new double[nt*rows*cols];


 
	for(int k=1;k<=nt;k++)
		{
		for (int i = 1; i < rows - 1; i++) {
		    for (int j = 1; j < cols - 1; j++) {
		        A[idx(i,j)] = 0.25 * (A[idx(i+1,j)] + A[idx(i-1,j)] + A[idx(i,j+1)] + A[idx(i,j-1)]);
		    	}
			}
		}

      

      

}
void uimr_init(int N,double l, double r, double t, double b)
{
  rows=rows/N;
  cols=cols/N;
  
  
  
  std::ofstream set("shape.dat");
  set << rows*cols << " " << 4 ;

  std::ofstream init("init.cache");
  for(int i=1;i<=rows;i++)
  {

    for(int j=1;j<=cols;j++)
    {
      if(i==1 && j==1)
      {
        init << l << " " << 0  << " " << 0 << " " << b << std::endl ;
      }

      else if(i==1 && j!=1 && j!=cols)
      {
        
        init << 0 << " " << 0  << " " << 0 << " " << b << std::endl;
      }

      else if(i==1 && j==cols)
      {

        init << 0 << " " << r  << " " << 0 << " " << b << std::endl;
      }
  

      else if(i==rows && j!=cols && j!=1)
      {
        
        init << 0 << " " << 0  << " " << t << " " << 0 << std::endl;
      }

      else if(i==rows && j==1)
      {
        init << l << " " << 0  << " " << t << " " << 0 << std::endl;
      }

      else if(i==rows && j==cols)
      {
        init << 0 << " " << r  << " " << t << " " << 0 << std::endl;
      }

      else if(j==1)
      {
        init << l << " " << 0  << " " << 0 << " " << 0 << std::endl;

      }

      else if(j==cols)
      {
        init << 0 << " " << r  << " " << 0 << " " << 0 << std::endl;

      }
      
      else{
        
        init << 0 << " " << 0  << " " << 0 << " " << 0 << std::endl;
      }
      

  }}}
  
  
  

	

void average_local_boundaries(std::string in_file, std::string rout_file,double l, double r, double t, double b)
{

	  int row, col;
	  std::ifstream shap("shape.dat");
	  shap >> row >> col ;
	  double arr[row][col];

	  std::ifstream file(in_file);
	  for(int i=0;i<row;i++)
	  {
		for(int j=0;j<col;j++)
		{
		  file >> arr[i][j];
		}
	  }

	  double tempx;
	  double tempy;

	  for(int i=0;i<row-1;i++)
	  {
		tempx=(arr[i][1]+arr[i+1][0])/2.0;
		arr[i][1]=tempx;
		arr[i+1][0]=tempx;

	  }


	  for(int j=0;j<col-1;j++)
	  {
		tempy=(arr[j][2]+arr[j][3])/2.0;
		arr[j][2]=tempy;
		arr[j][3]=tempy;
	  }
	  
	  for (int i = 0; i < row; i++) {

		int bx = i % rows;
		int by = i / cols;

		// LEFT boundary
		if (bx == 0)
		    arr[i][0] = 0.5 * (arr[i][0] + l);

		// RIGHT boundary
		if (bx == rows - 1)
		    arr[i][1] = 0.5 * (arr[i][1] + r);

		// TOP boundary
		if (by == cols - 1)
		    arr[i][2] = 0.5 * (arr[i][2] + t);

		// BOTTOM boundary
		if (by == 0)
		    arr[i][3] = 0.5 * (arr[i][3] + b);
		}
		
		
		
	  
	  
	  


	  std::ofstream out_file(rout_file);

	  for(int i=0;i<row;i++)
	  {
		for(int j=0;j<col;j++)
		{
		  out_file << arr[i][j] << " " ;
		}

		out_file << std::endl;
	  }
	  
	  
	  
	  
	  




}





void uimr_solve()
{
  double le,ri,to,bo;  
  std::string in_file  = "init.cache";
  std::string out_file = "bounds.dat";
  std::string temp_file = "temp.dat";
  
  const double tol = 1e-5;
  int max_iter=1000;
  bool write_out=false;
  
  

  for(int k=1;k<=max_iter+1;k++){
  
      std::ifstream ini(in_file);
      std::ofstream out(temp_file);
      bool converge=true;
      int id=0;
      
	  while(ini >> le >> ri >> to >> bo)
	  {
		  double le_old=le;
		  double ri_old=ri;
		  double to_old=to;
		  double bo_old=bo;
		  

	  double* A = Boundary_Initilise(rows, cols,le, ri, to, bo);
	  Solve_Laplacian(A);

	  double suml=0;
	  double sumr=0;

	  for(int n=0;n<rows;n++)
	  {
		suml+=A[idx(n,1)];
		sumr+=A[idx(n,cols-2)];
	  }

	  //le=suml/double(rows);
	  //ri=sumr/double(rows);
	  
	  le = 0.5 * (suml / rows) + 0.5 * A[idx(rows/2, 1)];
	  ri = 0.5 * (sumr / rows) + 0.5 * A[idx(rows/2, cols-2)];

	 double sumt=0;
	 double sumb=0;

	  for(int n=0;n<cols;n++)
	  {
		sumt+=A[idx(1,n)];
		sumb+=A[idx(rows-2,n)];
	  }

	  //to=sumt/double(cols);
	  //bo=sumb/double(cols);
	  
	  to = 0.5 * (sumt / cols) + 0.5 * A[idx(1, cols/2)];
	  bo = 0.5 * (sumb / cols) + 0.5 * A[idx(rows-2, cols/2)];
	  
	  
	  
	  
	  out << le << " " << ri << " " << to << " " << bo  << std::endl ; //writing the local boundaries in the out file 
	  
	  	
	  
	  double diff = std::max({std::abs(le - le_old) , std::abs(ri - ri_old) , std::abs(to - to_old), std::abs(bo - bo_old)});  //checking for convergence
	  
	  if(diff>tol)
	  {
	  	converge=false;
	  }
	  
	 
	  
	  if(write_out)     //checking for write_out, if true write_out
	  {
	  	//std::cout << id;
	  	std::string sol_name = "sol_" + std::to_string(id) + ".dat";
	  	std::ofstream sol(sol_name);

    	for (int i = 0; i < rows; i++) {
        	for (int j = 0; j < cols; j++) {
            	sol << A[idx(i,j)] << " " ;
        	}
        	sol << std::endl;
        	
    	}
    	
    	sol.close();
    	
    	
		}
	  
	  delete [] A;
	  id++;
	  
		}
	  
	  ini.close();
	  out.close();
	  
	  
	  average_local_boundaries(temp_file, out_file,100,50,50,80);
	  std::swap(in_file,out_file);
		
		
	 if(write_out==true){break;} // written out done!!
		
		
		
	  
	  
	  if(converge==true ) {
	  	std::cout << "Converged after " << k << " iterations\n";
	  	write_out=true;
	  	}
	  	
	  if(k==max_iter)
	  { std::cout << "Maximum Iteration Reached" << std::endl;
	    write_out=true;
	  }
	  	
	 
		
		
		
		
			
		
	//ini.close();
    //bounds.close();
    
    //std::cout << k ;
    

	

	}
	
	
	

}





int main()
{

   
  //double* A= Boundary_Initilise(rows/2, cols/2, 100,0,0,50);
  //Solve_Laplacian(A);
  uimr_init(N,100,50,50,80);
  uimr_solve();
 

  //delete [] A;
  return 0;



}
